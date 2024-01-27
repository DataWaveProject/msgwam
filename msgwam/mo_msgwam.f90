!>
!! This module is a collection of subroutines running the MS-GWaM 
!! (Multi Scale Gravity Wave Model) developed at the Goethe Uni.
!! MS-GWaM is a Lagrangian WKB raytracer in phase space based on 
!! the literature below:
!! Achatz et al. (2010), JFM.      -- theory
!! Achatz et al. (2017), QJRMS.    -- theory
!! Muraschko et al. (2015), QJRMS. -- phase space representation
!! Bölöni et al. (2016), JAS.      -- wave breaking
!! Orr et al. (2010), J. Climate   -- background source (Desaubies type)
!! Bölöni et al. (2021), JAS.      -- implementation in ICON
!! Kim et al. (2021),  JAS         -- implementation in ICON 
!!                                    (convective source, intermittency)
!!
!! @author Gergely Bölöni, Goethe Uni Frankfurt
!!
!! @par Revision History
!! Initial version by Gergely Bölöni, Goethe Uni Frankfurt (2016-01-28)
!! updates (listed below) by Young-Ha Kim, Goethe Uni Frankfurt (2019-01-21)
!! -- significant speedup by re-coding
!! -- updated lower boundary condition (launching ray volumes as source) 
!! -- energy based removal of ray volumes
!! -- revised ray volume propagation
!! -- convective sources implemented (Song and Chun, 2005; Kim et al., 2021)

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_msgwam

  USE mo_kind,                  ONLY: wp, vp, sp
  USE mo_mpi,                   ONLY: my_process_is_stdio, p_wait, work_mpi_barrier
  USE mo_exception,             ONLY: message, finish, message_text, &
                                      debug_messages_on, debug_messages_off
  USE mo_physical_constants,    ONLY: grav, rd, cpd, cvd, rd_o_cpd, earth_angular_velocity
  USE mo_grid_config,           ONLY: grid_sphere_radius, is_plane_torus
  USE mo_model_domain,          ONLY: t_patch
  USE mo_sync,                  ONLY: sync_patch_array, sync_patch_array_mult, SYNC_E, &
                                      SYNC_C, SYNC_C1
  USE mo_impl_constants,        ONLY: min_rlcell_int, min_rledge_int, success
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_nonhydro_types,        ONLY: t_nh_diag, t_nh_prog, t_nh_metrics
  USE mo_parallel_config,       ONLY: nproma
  USE mo_run_config,            ONLY: msg_level, ldynamics
  USE mo_dynamics_config,       ONLY: lcoriolis
  USE mo_math_constants,        ONLY: pi, pi2, pi_2, rad2deg, deg2rad
  USE mo_math_gradients,        ONLY: grad_green_gauss_cell
  USE mo_util_vgrid_types,      ONLY: vgrid_buffer
  USE mtime,                    ONLY: datetime, timeDelta, newTimedelta, &
                                      deallocateTimedelta, getTimedeltaFromDatetime, &
                                      getTotalMillisecondsTimedelta
  USE mo_vertical_grid,         ONLY: nrdmax
  USE mo_intp,                  ONLY: cell_avg
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp_lonlat_baryctr,   ONLY: inside_triangle
  USE mo_delaunay_types,        ONLY: t_point
  USE mo_math_divrot,           ONLY: div, div_avg
  USE mo_nh_testcases_nml,      ONLY: nh_test_name
  USE mo_timer
  USE mo_msgwam_config
  USE mo_setup_msgwam_interface
  USE mo_msgwam_util,           ONLY: smooth_vert, smooth_hori
  USE mo_gw_source_config,      ONLY: gws_conv_config
  USE mo_msgwam_diagnostics

  USE, INTRINSIC :: IEEE_ARITHMETIC

  ! Modules for computing the Orr et al. 2010 launch spectrum
  USE data_gwd,    ONLY : nslope, gfluxlaun, &
    &                     ggaussa, ggaussb, ngauss, gcoeff, lozpr

  IMPLICIT NONE

  PRIVATE

  INTEGER  ::  jg
  REAL(wp) ::  coef_sphere, coef_torus

! REAL(wp), PARAMETER :: bvf2_min = 1.E-12_wp    ! minimum for BVF (avoid division by 0)
  REAL(wp), PARAMETER :: bvf2_min = 3.E-8_wp     ! minimum for BVF (avoid division by 0)
                                                 ! set to be larger than f^2 everywhere
  REAL(wp), PARAMETER :: eps_wadens = 1.e-36_wp  ! wadens criterion to remove ray volumes.
! REAL(wp), PARAMETER :: eps_wadens = 0._wp      ! The criterion (even with zero) needs to be
                                                 ! used everywhere wadens is reduced, in order
                                                 ! not to miss being zero due to precision limit.
                                                 ! 1.e-36 must be very small enough to ignore and
                                                 ! still usable even if wp == sp

  PUBLIC  ::  gwdrag_msgwam
  PUBLIC  ::  sync_wave

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE gwdrag_msgwam ( dt_call,                   & ! input
                       &   mtime_datetime,            & ! input
                       &   p_sim_time,                & ! input
                       &   p_patch,p_metrics,         & ! input
                       &   p_int_state,               & ! input
                       &   rho,                       & ! input
                       &   pres,                      & ! input
                       &   pres_sfc,                  & ! input
                       &   temp    ,                  & ! input
                       &   temp_ifc,                  & ! input
                       &   theta,                     & ! input
                       &   theta_ifc,                 & ! input
                       &   u       ,                  & ! input
                       &   v       ,                  & ! input
                       &   p_fld                      ) ! inout

  TYPE(datetime),      POINTER,INTENT(IN)    :: mtime_datetime       ! date/time information
  TYPE(t_patch),        TARGET,INTENT(INOUT) :: p_patch              ! grid/patch info.
  TYPE(t_nh_metrics)          ,INTENT(IN)    :: p_metrics
  TYPE(t_int_state),    TARGET,INTENT(IN)    :: p_int_state          ! interpolation state
  TYPE(t_msgwam),              INTENT(INOUT) :: p_fld                ! the atm phys vars
  REAL(wp),                    INTENT(IN)    :: dt_call              ! time step
  REAL(wp),                    INTENT(IN)    :: p_sim_time           ! elapsed simulation time on this grid level
  REAL(wp),            POINTER,INTENT(IN)    :: rho (:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: pres (:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: pres_sfc (:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: temp(:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: temp_ifc(:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: theta(:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: theta_ifc(:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: u   (:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: v   (:,:,:)
  ! Local array bounds:
  INTEGER                                    :: nlev, nlevp1         ! number of full and half levels
  INTEGER                                    :: rl_start, rl_end
  INTEGER                                    :: i_startblk, i_endblk ! blocks
  INTEGER                                    :: i_startidx, i_endidx ! slices
  INTEGER                                    :: jk,jc,jb             ! block indeces
  INTEGER                                    :: maxindexu(3), maxindexv(3)
  ! Background
  REAL(wp)                                   :: rho_half(nproma,p_patch%nlevp1, p_patch%nblks_c)      ! rho at half levels
  REAL(wp)                                   :: bvf2_full(nproma, p_patch%nlev, p_patch%nblks_c)      ! N**2 at full levels
  REAL(wp)                                   :: bvf2_half(nproma, p_patch%nlevp1, p_patch%nblks_c)    ! N**2 at half levels
  REAL(wp)                                   :: u_fld(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! u at half levels
  REAL(wp)                                   :: v_fld(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! v at half levels
  REAL(wp)                                   :: kvisc(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! kinematic viscosity
  REAL(wp)                                   :: gammash2_full(nproma, p_patch%nlev, p_patch%nblks_c)  ! inverse pinc scale height
                                                                                                      ! pinc:
                                                                                                      ! pseudo-incompressible
  REAL(wp)                                   :: gammash2_half(nproma, p_patch%nlevp1, p_patch%nblks_c)! inverse pinc scale height
                                                                                                      ! pinc:
                                                                                                      ! pseudo-incompressible
  REAL(wp)                                   :: dn2dz(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! vertical grad of bvf2 /2
  REAL(wp)                                   :: dg2dz(nproma, p_patch%nlevp1, p_patch%nblks_c)        ! vertical grad of gammash2 /2
  REAL(wp)                                   :: dudz (nproma, p_patch%nlevp1, p_patch%nblks_c)        ! vertical (radial) grad of u
  REAL(wp)                                   :: dvdz (nproma, p_patch%nlevp1, p_patch%nblks_c)        ! vertical (radial) grad of v
  REAL(wp)                                   :: n2_hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! zonal/meridional grad of N**2
  REAL(wp)                                   :: g2_hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! zonal/meridional grad of gammash2
  REAL(wp)                                   :: u_hgrad (2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! zonal/meridional grad of u
  REAL(wp)                                   :: v_hgrad (2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! zonal/meridional grad of v
  REAL(wp)                                   :: dn2dz_hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! vert & zonal/meridional grad of N**2
  REAL(wp)                                   :: dg2dz_hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! vert & zonal/meridional grad of gammash2
  REAL(wp)                                   :: dudz_hgrad (2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! vert & zonal/meridional grad of u
  REAL(wp)                                   :: dvdz_hgrad (2, nproma, p_patch%nlevp1, p_patch%nblks_c)  ! vert & zonal/meridional grad of v
  REAL(wp)                                   :: fc2(nproma, p_patch%nblks_c)                          ! f_c**2
  REAL(wp)                                   :: fdfdlat(nproma, p_patch%nblks_c)                      ! meridional gradient of f**2
  ! 1D Diagnostics
  REAL(wp)                                   :: kd(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: kd2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: ld(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: ld2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: md(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: md2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: cgz_diag(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: B2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: B2_save(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: A(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: A_save(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: mB2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: mB2_save(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: prec(nproma)                                  ! total precip
! REAL(wp),                      PARAMETER   :: lim_tend = 0.05_wp
  REAL(wp),                      PARAMETER   :: lim_tend_tmp = 1._wp
  ! 
  LOGICAL                                    :: ltimeadd_bg
  LOGICAL                                    :: ltimeadd_cv
  INTEGER                                    :: iter_sm
  INTEGER                                    :: nrays_vacate
  INTEGER                                    :: idiag

#ifndef __msgwam1d
  CALL message('gwdrag_msgwam', 'Transient MS-GWaM called')
#else
  CALL message('gwdrag_msgwam', 'Transient MS-GWaM-1D called')
  ! off : lhsmooth, hard-coded tendency limiter
  ! need to turn on  lmsgwam_pmomflux
  !              off imethod_split/merge (as not implemented)
#endif

  IF (is_plane_torus) THEN
    coef_sphere = 0.     ;  coef_torus = 1._wp
  ELSE
    coef_sphere = 1._wp  ;  coef_torus = 0.
  END IF

  ! Number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  ! Domain ID
  jg     = p_patch%id

  ! Exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c + 1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

  ! Decide whether new waves are emmited (new ray volumes are added)
  ! ltimeadd_bg = .T. --> Desaubies type background source launched
  ! ltimeadd_cv = .T. --> convective source launched
  IF (nrays_add_bg(jg) /= 0) THEN
    IF (ltest_hprop .AND. l1ray) THEN ! horizontal propagation test with single ray
                                                        ! or simulation with idealized source
      ltimeadd_bg = (MOD(p_sim_time,dt_add)<dt_call*0.999_wp .AND. p_sim_time <= 0._wp)
    ELSEIF (ltest_hprop) THEN ! horizontal propagation test
      ltimeadd_bg = (MOD(p_sim_time,dt_add)<dt_call*0.999_wp .AND. p_sim_time <= 43200._wp) ! half day emission
    ELSE
      ltimeadd_bg = (MOD(p_sim_time,dt_add)<dt_call*0.999_wp)
    ENDIF
  ELSE
    ltimeadd_bg = .FALSE.
  ENDIF
  IF ( gws_conv_config%n_source(jg) > 0 ) THEN
    ltimeadd_cv = ( p_sim_time /= 0._wp )
  ELSE
    ltimeadd_cv = .FALSE.
  ENDIF

  ! Idealized gwp test case: no background source, no convective source
  IF (nh_test_name=='gwp') THEN
    ltimeadd_bg = .FALSE.
    ltimeadd_cv = .FALSE.
  ENDIF

! Print whether source is launched
  IF (msg_level >= 12) THEN
    WRITE(message_text,'(a,L8)') 'ltimeadd_bg:', ltimeadd_bg
    IF ( gws_conv_config%n_source(jg) > 0 )  &
      &  WRITE(message_text,'(a,L8)') 'ltimeadd_cv:', ltimeadd_cv
    CALL message('', TRIM(message_text))
  ENDIF

  !=========================== FIELDS/GRADS ===========================
  ! Calculation of all resolved variables and their vertical gradients.
  !====================================================================

  IF (timers_level > 4) CALL timer_start(timer_msgwam_fieldsgrads)

  ! Calculate derived fields from the resolved flow (e.g. Brunt-Väisälä freq and density scale height)
  ! and their vertical gradients (needed for solving the ray equations)
  CALL fieldsgrads(p_patch       = p_patch,                         & ! grid information                    (in)
                   p_int_state   = p_int_state,                     & ! interpolation data of patch         (in)
                   rl_start      = rl_start,                        & !                                     (in)
                   rl_end        = rl_end,                          & !                                     (in)
                   nlev          = nlev,                            & ! no. of full levels                  (in)
                   z             = p_metrics%z_mc(:, :, :),         & ! full level heights                  (in)
                   zhalf         = p_metrics%z_ifc(:, :, :),        & ! half level heights                  (in)
                   dz            = p_metrics%ddqz_z_full(:, :, :),  & ! full layer thickness                (in)
                   dzhalf        = p_metrics%ddqz_z_half(:, :, :),  & ! half layer thickness                (in)
                   wgtfac_c      = p_metrics%wgtfac_c(:, :, :),     & !                                     (in)
                   rho           = rho (:, :, :),                   & ! density at full levels              (in)   
                   temp          = temp(:, :, :),                   & ! temperature at full levels          (in)   
                   u             = u(:, :, :),                      & ! zonal wind (full levels)            (in)
                   v             = v(:, :, :),                      & ! meridional wind (full levels)       (in)
                   fc2           = fc2 (:, :),                      & ! f_c**2                              (out)
                   fdfdlat       = fdfdlat(:, :),                   & ! meridional gradient of f**2         (out)
                   rho_half      = rho_half(:, :, :),               & ! density at half levels              (out)
                   u_fld         = u_fld(:, :, :),                  & ! zonal wind                          (out)
                   v_fld         = v_fld(:, :, :),                  & ! meridional wind                     (out)
                   bvf2_full     = bvf2_full(:, :, :),              & ! Brunt-Väisälä frequency**2          (out)
                   bvf2_half     = bvf2_half(:, :, :),              & ! Brunt-Väisälä frequency**2          (out)
                   gammash2_full = gammash2_full(:, :, :),          & ! inverse pinc scaleheight**2         (out)
                   gammash2_half = gammash2_half(:, :, :),          & ! inverse pinc scaleheight**2         (out)
                   kvisc         = kvisc(:, :, :),                  & ! Kinematic viscosity at half levels  (out)
                   dn2dz         = dn2dz(:, :, :),                  & ! vertical grad of bvf2               (out)
                   dg2dz         = dg2dz(:, :, :),                  & ! vertical grad of gammash2           (out)
                   dudz          = dudz (:, :, :),                  & ! vertical (radial) grad of u         (out)
                   dvdz          = dvdz (:, :, :),                  & ! vertical (radial) grad of v         (out)
                   n2_hgrad      = n2_hgrad(:, :, :, :),            & ! zonal / meridional grad of BVF**2   (out)
                   g2_hgrad      = g2_hgrad(:, :, :, :),            & ! zonal / meridional grad of gamma**2 (out)
                   u_hgrad       = u_hgrad(:, :, :, :),             & ! zonal / meridional grad of u        (out)
                   v_hgrad       = v_hgrad(:, :, :, :),             & ! zonal / meridional grad of v        (out)
                   dn2dz_hgrad   = dn2dz_hgrad(:, :, :, :),         & ! vert & zonal / meridional grad of BVF**2   (out)
                   dg2dz_hgrad   = dg2dz_hgrad(:, :, :, :),         & ! vert & zonal / meridional grad of gamma**2 (out)
                   dudz_hgrad    = dudz_hgrad(:, :, :, :),          & ! vert & zonal / meridional grad of u        (out)
                   dvdz_hgrad    = dvdz_hgrad(:, :, :, :)           ) ! vert & zonal / meridional grad of v        (out)

  IF (timers_level > 4) CALL timer_stop(timer_msgwam_fieldsgrads)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, prec,                      &
!$OMP            kd, kd2, ld, ld2, md, md2, cgz_diag, A, A_save, B2, B2_save, mB2, mB2_save)


  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

    ! Considering N^2 changes between time steps, remove ray volumes where N^2 <= f^2.
    ! This will not happen as long as bvf2_min > f^2 everywhere.
!   DO jc = i_startidx, i_endidx
!     DO jk = 1, nlevp1
!       IF (bvf2_half(jc,jk,jb) > fc2(jc))  CYCLE
!       WHERE ( p_ray(jg)% iexist(jc,:,jb) == jk )
!         p_ray(jg)% iexist(jc,:,jb) = 0
!         p_ray(jg)% specid(jc,:,jb) = 0
!       END WHERE
!     ENDDO
!   ENDDO

    !======================== MERGE/REMOVE RAYS ==========================
    ! Remove ray volumes if there are more than allowed (maxrays_bg(jg) 
    ! or maxrays_cv(jg)) in order to keep computational costs under control.
    ! The removal is done separately for ray volumes coming from the 
    ! background and the convective sources. 
    ! TODO: the removal should be replaced by a merging procedure
    !=====================================================================

    IF (timers_level > 4) CALL timer_start(timer_msgwam_remove_rays)

    !======================== Diagnostics ==================================
    ! WA diagnostic
    !=======================================================================

    CALL idx_rayedge( nlev,i_startidx,i_endidx,1,nrays(jg),& ! (in)
                      p_metrics%z_mc     (:,:,jb),        & ! (in)
                      p_metrics%z_ifc    (:,:,jb),        & ! (in)
                      p_ray(jg)%iexist   (:,:,jb),        & ! (in)
                      p_ray(jg)%specid   (:,:,jb),        & ! (in)
                      p_ray(jg)%jk_active(:,:,jb),        & ! (in)
                      p_ray(jg)%z        (:,:,jb),        & ! (in)
                      p_ray(jg)%dz       (:,:,jb),        & ! (in)
                      p_ray(jg)%jk_full_rtop(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_full_rbot(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_half_rtop(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_half_rbot(:,:,jb)      ) ! (inout)

    CALL project_action(nlev        = nlev,                             & !
                        i_startidx  = i_startidx,                       & !
                        i_endidx    = i_endidx,                         & !
                        jray_start  = 1,                                & !
                        jray_end    = nrays(jg),                        & !
                        z           = p_metrics%z_mc(:,:,jb),           & !
                        zhalf       = p_metrics%z_ifc(:,:,jb),          & !
                        cellarea    = p_patch%cells%area(:,jb),         & !
                        zray        = p_ray(jg)%z(:,:,jb),              & !
                        dzray       = p_ray(jg)%dz(:,:,jb),             & !
                        coslatray   = p_ray(jg)%coslat(:,:,jb),         & !
                        dlatray     = p_ray(jg)%dlat(:,:,jb),           & !
                        dlonray     = p_ray(jg)%dlon(:,:,jb),           & !
                        dkray       = p_ray(jg)%dk(:,:,jb),             & !
                        dlray       = p_ray(jg)%dl(:,:,jb),             & !
                        dmray       = p_ray(jg)%dm(:,:,jb),             & !
                        dens        = p_ray(jg)%wadens(:,:,jb),         & !
                        iexist      = p_ray(jg)%iexist(:,:,jb),         & !
                        specid      = p_ray(jg)%specid(:,:,jb),         & !
                        jk_active   = p_ray(jg)%jk_active(:,:,jb),      & !
                        jkmin_full  = p_ray(jg)%jk_full_rtop(:,:,jb),   & !
                        jkmax_full  = p_ray(jg)%jk_full_rbot(:,:,jb),   & !
                        jkmin_half  = p_ray(jg)%jk_half_rtop(:,:,jb),   & !
                        jkmax_half  = p_ray(jg)%jk_half_rbot(:,:,jb),   & !
                        action      = p_fld%action_mgm_1(:,:,jb),       & !
                        diag        = 0,                                & !
                        jg          = jg                                ) !


    ! Remove ray volumes from convective sources
    IF ( ltimeadd_cv .AND. ( imethod_split <= 0 .AND. imethod_merge <= 0 ) ) THEN
      CALL remove_rays(      nlev         = nlev,                            & ! no. of full levels       (in)
                             i_startidx   = i_startidx,                      & ! first index of the block (in)
                             i_endidx     = i_endidx,                        & ! last index of the block  (in)
                             jray_start   = jray_offset_cv(jg)+1,            &
                             jray_end     = jray_offset_cv(jg)+nrays_cv(jg), &
                             nrays_vacate = nrays_add_cv(jg)*nlaunch_max_cv, &
                             z            = p_metrics%z_mc(:,:,jb),          & ! full level heights       (in)
                             zhalf        = p_metrics%z_ifc(:,:,jb),         & ! half level heights       (in)
                             dlonray      = p_ray(jg)%dlon(:,:,jb),          & ! ray size in lon          (in)
                             dlatray      = p_ray(jg)%dlat(:,:,jb),          & ! ray size in lat          (in)
                             dzray        = p_ray(jg)%dz(:,:,jb),            & ! ray size in z            (in)
                             coslatray    = p_ray(jg)%coslat(:,:,jb),        & ! cosine of ray lat        (in)
                             kray         = p_ray(jg)%k(:,:,jb),             & ! horiz (lon) wavenumber   (in)
                             dkray        = p_ray(jg)%dk(:,:,jb),            & ! horiz (lon) wavenumber   (in)
                             lray         = p_ray(jg)%l(:,:,jb),             & ! horiz (lat) wavenumber   (in)
                             dlray        = p_ray(jg)%dl(:,:,jb),            & ! horiz (lat) wavenumber   (in)
                             mray         = p_ray(jg)%m(:,:,jb),             & ! vertical wavenumber      (in)
                             dmray        = p_ray(jg)%dm(:,:,jb),            & ! vertical wavenumber      (in)
                             dens         = p_ray(jg)%wadens(:,:,jb),        & ! wave action density      (inout)
                             iexist       = p_ray(jg)%iexist(:,:,jb),        & ! existence of ray         (inout)
                             specid       = p_ray(jg)%specid(:,:,jb),        & ! spectral id of rays      (inout)
                             bvf2         = bvf2_half(:,:,jb),               & ! Brunt-Väisala freq**2    (in)
                             gammash2     = gammash2_half(:,:,jb),           & ! inverse pinc scale height(in)
                             fc2          = fc2(:,jb)                        ) ! Coriolis parameter**2    (in)
    END IF

    ! Remove ray volumes from the Desaubies type background source
    IF ( ltimeadd_bg .AND. ( imethod_split <= 0 .AND. imethod_merge <= 0 ) ) THEN
      CALL remove_rays(      nlev         = nlev,                            & ! no. of full levels       (in)
                             i_startidx   = i_startidx,                      & ! first index of the block (in)
                             i_endidx     = i_endidx,                        & ! last index of the block  (in)
                             jray_start   = jray_offset_bg(jg)+1,            &
                             jray_end     = jray_offset_bg(jg)+nrays_bg(jg), &
                             nrays_vacate = nrays_add_bg(jg)*nlaunch_max_bg, &
                             z            = p_metrics%z_mc(:,:,jb),          & ! full level heights       (in)
                             zhalf        = p_metrics%z_ifc(:,:,jb),         & ! half level heights       (in)
                             dlonray      = p_ray(jg)%dlon(:,:,jb),          & ! ray size in lon          (in)
                             dlatray      = p_ray(jg)%dlat(:,:,jb),          & ! ray size in lat          (in)
                             dzray        = p_ray(jg)%dz(:,:,jb),            & ! ray size in z            (in)
                             coslatray    = p_ray(jg)%coslat(:,:,jb),        & ! cosine of ray lat        (in)
                             kray         = p_ray(jg)%k(:,:,jb),             & ! horiz (lon) wavenumber   (in)
                             dkray        = p_ray(jg)%dk(:,:,jb),            & ! horiz (lon) wavenumber   (in)
                             lray         = p_ray(jg)%l(:,:,jb),             & ! horiz (lat) wavenumber   (in)
                             dlray        = p_ray(jg)%dl(:,:,jb),            & ! horiz (lat) wavenumber   (in)
                             mray         = p_ray(jg)%m(:,:,jb),             & ! vertical wavenumber      (in)
                             dmray        = p_ray(jg)%dm(:,:,jb),            & ! vertical wavenumber      (in)
                             dens         = p_ray(jg)%wadens(:,:,jb),        & ! wave action density      (inout)
                             iexist       = p_ray(jg)%iexist(:,:,jb),        & ! existence of ray         (inout)
                             specid       = p_ray(jg)%specid(:,:,jb),        & ! spectral id of rays      (inout)
                             bvf2         = bvf2_half(:,:,jb),               & ! Brunt-Väisala freq**2    (in)
                             gammash2     = gammash2_half(:,:,jb),           & ! inverse pinc scale height(in)
                             fc2          = fc2(:,jb)                        ) ! Coriolis parameter**2    (in)
    END IF

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_remove_rays)


    !========================= INIT GW SOURCES ============================
    ! Launch GW sources (currently 3 kinds):
    ! 1) convective (init_gw_conv)
    ! 2) Desaubies type background (init_gw_orretal)
    ! 3) Idealized gwp (set up in src/testcases/mo_nh_isotherm_rest_atm_gwp.f90)
    ! TODO: mountain waves, jets/fronts
    !======================================================================

    ! Convective source
    IF ( ltimeadd_cv )                                                       &
      &  CALL init_gw_conv( nlev,jb,i_startidx,i_endidx,jray_offset_cv(jg),  &
      &                     nlaunch_max_cv,                                  &
      &                     p_metrics%z_mc(:,:,jb),p_metrics%z_ifc(:,:,jb),  &
      &                     p_patch%cells%center(:,jb)%lon,                  &
      &                     p_patch%cells%center(:,jb)%lat,                  &
      &                     p_patch%cells%area(:,jb),                        &
      &                     p_fld% flag_cgw(:,jb) )

    ! Desaubies type background source
    IF (ltimeadd_bg) THEN

      IF (timers_level > 4) CALL timer_start(timer_msgwam_init_gw_orretal)

      ! Caculate total precipitation (input for the Orr et al., 2010 source)
!     IF ( lozpr )  prec(:) =   prm_diag%rain_gsp_rate(:,jb) &  ! rain_gsp
!       &                     + prm_diag%snow_gsp_rate(:,jb) &  ! snow_gsp
!       &                     + prm_diag%rain_con_rate(:,jb) &  ! rain_con
!       &                     + prm_diag%snow_con_rate(:,jb)    ! snow con
      IF ( lozpr )  prec(:) = 0._wp

      CALL init_gw_orretal(nlev      = nlev,                            & ! no. of full levels       (in)
                          i_startidx = i_startidx,                      & ! first index of the block (in)
                          i_endidx   = i_endidx,                        & ! last index of the block  (in)
                          jray_start = jray_offset_bg(jg)+1,            &
                          jray_end   = jray_offset_bg(jg)+nrays_bg(jg), &
                          mdatetime  = mtime_datetime,                  & ! last index of the block   (in)
                          zhalf      = p_metrics%z_ifc(:,:,jb),         & ! half level heights        (in)
                          dzhalf     = p_metrics%ddqz_z_half(:,:,jb),   & ! half layer thickness      (in)
                          z          = p_metrics%z_mc(:,:,jb),          & ! full level heights        (in)
                          dz         = p_metrics%ddqz_z_full(:,:,jb),   & ! full layer thickness      (in)
                          specid     = p_ray(jg)%specid(:,:,jb),        & ! spectral id of rays       (inout)
                          lonray     = p_ray(jg)%lon(:,:,jb),           & ! ray position lon          (inout)
                          dlonray    = p_ray(jg)%dlon(:,:,jb),          & ! ray size in lon           (inout)
                          latray     = p_ray(jg)%lat(:,:,jb),           & ! ray position lat          (inout)
                          dlatray    = p_ray(jg)%dlat(:,:,jb),          & ! ray size in lat           (inout)
                          zray       = p_ray(jg)%z(:,:,jb),             & ! ray position z            (inout)
                          dzray      = p_ray(jg)%dz(:,:,jb),            & ! ray size in z             (inout)
                          coslatray  = p_ray(jg)%coslat(:,:,jb),        & ! cosine of ray lat         (inout)
                          kray       = p_ray(jg)%k(:,:,jb),             & ! ray horiz (lon) wavenumber(inout)
                          dkray      = p_ray(jg)%dk(:,:,jb),            & ! ray size in k             (inout)
                          lray       = p_ray(jg)%l(:,:,jb),             & ! ray horiz (lat) wavenumber(inout)
                          dlray      = p_ray(jg)%dl(:,:,jb),            & ! ray size in l             (inout)
                          mray       = p_ray(jg)%m(:,:,jb),             & ! ray vertical wavenumber   (inout)
                          dmray      = p_ray(jg)%dm(:,:,jb),            & ! ray size in m             (inout)
                          dens       = p_ray(jg)%wadens(:,:,jb),        & ! ray wave action density   (inout)
                          iexist     = p_ray(jg)%iexist(:,:,jb),        & ! existence of ray          (inout)
                          jk_active  = p_ray(jg)%jk_active(:,:,jb),     & ! jk index for launching    (inout)
                          jr_last    = p_ray(jg)%jr_last(:,:,jb),       & ! last-launched ray index   (inout)
                          jklaunch   = jklaunch_bg(jg),                 & ! launch level index        (in)
                          nlaunch_max= nlaunch_max_bg,                  & ! launch level index        (in)
                          cellarea   = p_patch%cells%area(:,jb),        & ! area of grid cel          (in)
                          lon        = p_patch%cells%center(:,jb)%lon,  & ! latitude of cell center   (in)
                          lat        = p_patch%cells%center(:,jb)%lat,  & ! latitude of cell center   (in)
                          coslat     = p_mgmgrid(jg)%coslat(:,jb),      & ! cosine of latitude        (in)
                          rhol       = rho_half(:,jklaunch_bg(jg),jb),  & ! rho at launch level       (in)
                          bvfl2      = bvf2_half(:,jklaunch_bg(jg),jb), & ! N**2 at launch level      (in)
                          lat_prof_bw= p_lfluxbg(jg)%lat_prof_bw(:,jb), & ! latitudinal factor        (in)
                          lat_prof_bs= p_lfluxbg(jg)%lat_prof_bs(:,jb), & ! latitudinal factor        (in)
                          lat_prof   = p_lfluxbg(jg)%lat_prof(:,jb),    & ! latitudinal factor        (inout)
                          prec       = prec(:)                          ) ! total precipitation       (in)

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_init_gw_orretal)

    ENDIF


    !======================== INDEX OF RAY EDGES ===========================
    ! Calculate the vertical level index (jk) closest to the bottom and top
    ! of the ray volume.
    ! TODO: something similar for horizontal edges if horizontal propagation
    !=======================================================================


    IF (timers_level > 4) CALL timer_start(timer_msgwam_saturation)    ! temporary

    CALL idx_rayedge( nlev,i_startidx,i_endidx,1,nrays(jg),& ! (in)
                      p_metrics%z_mc     (:,:,jb),        & ! (in)
                      p_metrics%z_ifc    (:,:,jb),        & ! (in)
                      p_ray(jg)%iexist   (:,:,jb),        & ! (in)
                      p_ray(jg)%specid   (:,:,jb),        & ! (in)
                      p_ray(jg)%jk_active(:,:,jb),        & ! (in)
                      p_ray(jg)%z        (:,:,jb),        & ! (in)
                      p_ray(jg)%dz       (:,:,jb),        & ! (in)
                      p_ray(jg)%jk_full_rtop(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_full_rbot(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_half_rtop(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_half_rbot(:,:,jb)      ) ! (inout)

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_saturation)     !  temporary

    CALL project_action(nlev        = nlev,                             & !
                        i_startidx  = i_startidx,                       & !
                        i_endidx    = i_endidx,                         & !
                        jray_start  = 1,                                & !
                        jray_end    = nrays(jg),                        & !
                        z           = p_metrics%z_mc(:,:,jb),           & !
                        zhalf       = p_metrics%z_ifc(:,:,jb),          & !
                        cellarea    = p_patch%cells%area(:,jb),         & !
                        zray        = p_ray(jg)%z(:,:,jb),              & !
                        dzray       = p_ray(jg)%dz(:,:,jb),             & !
                        coslatray   = p_ray(jg)%coslat(:,:,jb),         & !
                        dlatray     = p_ray(jg)%dlat(:,:,jb),           & !
                        dlonray     = p_ray(jg)%dlon(:,:,jb),           & !
                        dkray       = p_ray(jg)%dk(:,:,jb),             & !
                        dlray       = p_ray(jg)%dl(:,:,jb),             & !
                        dmray       = p_ray(jg)%dm(:,:,jb),             & !
                        dens        = p_ray(jg)%wadens(:,:,jb),         & !
                        iexist      = p_ray(jg)%iexist(:,:,jb),         & !
                        specid      = p_ray(jg)%specid(:,:,jb),         & !
                        jk_active   = p_ray(jg)%jk_active(:,:,jb),      & !
                        jkmin_full  = p_ray(jg)%jk_full_rtop(:,:,jb),   & !
                        jkmax_full  = p_ray(jg)%jk_full_rbot(:,:,jb),   & !
                        jkmin_half  = p_ray(jg)%jk_half_rtop(:,:,jb),   & !
                        jkmax_half  = p_ray(jg)%jk_half_rbot(:,:,jb),   & !
                        action      = p_fld%action_mgm_2(:,:,jb),       & !
                        diag        = 0,                                & !
                        jg          = jg                                ) !


    !============================ SATURATION ==============================
    ! Wave breaking scheme based on static instability criterion: if the GW 
    ! turns the potential temperture gradient to negative at a certain
    ! height, the wave action of ray volumes are reduced so that static 
    ! stability sets in again. Saturation is (should be!) calculated for 
    ! GWs from all sources (background, convective) together
    !======================================================================

    IF (lsaturation) THEN

      IF (timers_level > 4) CALL timer_start(timer_msgwam_saturation)

      CALL saturation(p_patch    = p_patch,                     & !                          (in)
                      jb         = jb,                          & !                          (in)
                      nlev       = nlev,                        & ! no. of full levels       (in)
                      i_startidx = i_startidx,                  & ! first index of the block (in)
                      i_endidx   = i_endidx,                    & ! last index of the block  (in)
                      jray_start = 1,                           & !                          (in)
                      jray_end   = nrays(jg),                   & !                          (in)
                      z          = p_metrics%z_mc(:,:,jb),      & ! full level heights       (in)
                      zhalf      = p_metrics%z_ifc(:,:,jb),     & ! half level heights       (in)
                      jk_active  = p_ray(jg)%jk_active(:,:,jb), & ! launch level index       (in)
                      zray       = p_ray(jg)%z(:,:,jb),         & ! position                 (in)
                      dzray      = p_ray(jg)%dz(:,:,jb),        & ! size in z-dir            (in)
                      coslatray  = p_ray(jg)%coslat(:,:,jb),    & ! cosine latitude          (in)
                      dlatray    = p_ray(jg)%dlat(:,:,jb),      & ! meridional extent        (in)
                      dlonray    = p_ray(jg)%dlon(:,:,jb),      & ! zonal extent             (in)
                      kray       = p_ray(jg)%k(:,:,jb),         & ! horiz (lon) wavenumber   (in)
                      dkray      = p_ray(jg)%dk(:,:,jb),        & ! size in k-dir            (in)
                      lray       = p_ray(jg)%l(:,:,jb),         & ! horiz (lat) wavenumber   (in)
                      dlray      = p_ray(jg)%dl(:,:,jb),        & ! size in l-dir            (in)
                      mray       = p_ray(jg)%m(:,:,jb),         & ! vertical wavenumber      (in)
                      dmray      = p_ray(jg)%dm(:,:,jb),        & ! size in m-dir            (in)
                      dens       = p_ray(jg)%wadens(:,:,jb),    & ! wave action density      (inout)
                      iexist     = p_ray(jg)%iexist(:,:,jb),    & ! existence of ray         (inout)
                      specid     = p_ray(jg)%specid(:,:,jb),    & ! spectral id of rays      (inout)
                      jkmin      = p_ray(jg)%jk_full_rtop(:,:,jb), & ! jk closest to ray-v top  (in)
                      jkmax      = p_ray(jg)%jk_full_rbot(:,:,jb), & ! jk closest to ray-v bot  (in)
                      bvf2       = bvf2_half(:,:,jb),           & ! Brunt-Väisala freq**2    (in)
                      gammash2   = gammash2_half(:,:,jb),       & ! inverse pinc scale height(in)
                      fc2        = fc2(:,jb),                   & ! Coriolis parameter**2    (in)
                      rho        = rho_half(:,:,jb),            & ! rho at half levels       (in)
                      dt         = dt_call,                     & ! time step for saturation (in)
                      kd         = kd(:,:),                     & ! k diagnostic             (out)
                      kd2        = kd2(:,:),                    & ! k^2 diagnostic           (out)
                      ld         = ld(:,:),                     & ! l diagnostic             (out)
                      ld2        = ld2(:,:),                    & ! l^2 diagnostic           (out)
                      md         = md(:,:),                     & ! m^2 diagnostic           (out)
                      md2        = md2(:,:),                    & ! m^2 diagnostic           (out)
                      cgz_diag   = cgz_diag(:,:),               & ! cgz diagnostic           (out)
                      A          = A(:,:),                      & ! B^2 diagnostic           (out)
                      A_save     = A_save(:,:),                 & ! B^2 diagnostic           (out)
                      B2         = B2(:,:),                     & ! B^2 diagnostic           (out)
                      B2_save    = B2_save(:,:),                & ! B^2 diagnostic           (out)
                      mB2_new    = mB2(:,:),                    & ! m^2B^2 diagnostic        (out)
                      mB2        = mB2_save(:,:)                ) ! m^2B^2 diagnostic        (out)

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_saturation)

    ENDIF


    !============================ WAVE2GRID ===============================
    ! Project quantities of the Lagrangian ray volumes to the Eulerian grid 
    ! and calculate momentum, pseudo-momentum, and temperature fluxes.
    ! The projection is done separately for GWs from the background and 
    ! from the convective sources
    !======================================================================

    IF (timers_level > 4) CALL timer_start(timer_msgwam_wave2grid)

    ! wave2grid
    IF (nrays_add_bg(jg) /= 0) THEN
      CALL wave2grid(nlev            = nlev,                            & ! no. of full levels       (in)
                     i_startidx      = i_startidx,                      & ! first index of the block (in)
                     i_endidx        = i_endidx,                        & ! last index of the block  (in)
                     jray_start      = 1,                               & !                          (in)
                     jray_end        = nrays(jg),                       & !                          (in)
                     z               = p_metrics%z_mc(:,:,jb),          & ! full level heights       (in)
                     zhalf           = p_metrics%z_ifc(:,:,jb),         & ! half level heights       (in)
                     fc              = p_patch%cells%f_c(:,jb),         & ! Coriolis parameter       (in)
                     cellarea        = p_patch%cells%area(:,jb),        & ! area of grid cell        (in)
                     zray            = p_ray(jg)%z(:,:,jb),             & ! position                 (in)
                     dzray           = p_ray(jg)%dz(:,:,jb),            & ! size in z-dir            (in)
                     coslatray       = p_ray(jg)%coslat(:,:,jb),        & ! cosine latitude          (in)
                     dlatray         = p_ray(jg)%dlat(:,:,jb),          & ! meridional extent        (in)
                     dlonray         = p_ray(jg)%dlon(:,:,jb),          & ! zonal extent             (in)
                     kray            = p_ray(jg)%k(:,:,jb),             & ! horiz (lon) wavenumber   (in)
                     dkray           = p_ray(jg)%dk(:,:,jb),            & ! size in k-dir            (in)
                     lray            = p_ray(jg)%l(:,:,jb),             & ! horiz (lat) wavenumber   (in)
                     dlray           = p_ray(jg)%dl(:,:,jb),            & ! size in l-dir            (in)
                     mray            = p_ray(jg)%m(:,:,jb),             & ! vertical wavenumber      (in)
                     dmray           = p_ray(jg)%dm(:,:,jb),            & ! size in m-dir            (in)
                     dens            = p_ray(jg)%wadens(:,:,jb),        & ! wave action density      (in)
                     iexist          = p_ray(jg)%iexist(:,:,jb),        & ! existence of ray         (in)
                     specid          = p_ray(jg)%specid(:,:,jb),        & ! spectral id of rays      (in)
                     jk_active       = p_ray(jg)%jk_active(:,:,jb),     & ! launch level index       (in)
                     jkmin_full      = p_ray(jg)%jk_full_rtop(:,:,jb),  & ! jk closest to ray-v top  (in)
                     jkmax_full      = p_ray(jg)%jk_full_rbot(:,:,jb),  & ! jk closest to ray-v bot  (in)
                     jkmin_half      = p_ray(jg)%jk_half_rtop(:,:,jb),  & ! jk closest to ray-v top  (in)
                     jkmax_half      = p_ray(jg)%jk_half_rbot(:,:,jb),  & ! jk closest to ray-v bot  (in)
                     theta           = theta(:,:,jb),                   & ! potential temperature    (in)
                     bvf2_full       = bvf2_full(:,:,jb),               & ! Brunt-Väisala freq**2    (in)
                     bvf2_half       = bvf2_half(:,:,jb),               & ! Brunt-Väisala freq**2    (in)
                     gammash2_full   = gammash2_full(:,:,jb),           & ! inverse pinc scale height(in)
                     gammash2_half   = gammash2_half(:,:,jb),           & ! inverse pinc scale height(in)
                     fc2             = fc2(:,jb),                       & ! Coriolis parameter**2    (in)
                     lcalc_flux_4dir = lcalc_flux_4dir_bg(jg),          & ! flag to calc. 4-fluxes   (in)
                     u               = u(:, :, jb),                     & ! zonal wind (full levels)      (in)
                     v               = v(:, :, jb),                     & ! meridional wind (full levels) (in)
                     uuflux          = p_fld% uufl_mgm    (:,:,jb),     & ! uu mometum flux          (out)
                     uvflux          = p_fld% uvfl_mgm    (:,:,jb),     & ! uv mometum flux          (out)
                     uwflux          = p_fld% uwfl_mgm    (:,:,jb),     & ! uw mometum flux          (out)
                     vvflux          = p_fld% vvfl_mgm    (:,:,jb),     & ! vv mometum flux          (out)
                     vwflux          = p_fld% vwfl_mgm    (:,:,jb),     & ! vw mometum flux          (out)
                     uupflux         = p_fld% uupfl_mgm   (:,:,jb),     & ! uu pseudo mometum flux   (out)
                     uvpflux         = p_fld% uvpfl_mgm   (:,:,jb),     & ! uv pseudo mometum flux   (out)
                     uwpflux         = p_fld% uwpfl_mgm   (:,:,jb),     & ! uw pseudo mometum flux   (out)
                     vvpflux         = p_fld% vvpfl_mgm   (:,:,jb),     & ! vv pseudo mometum flux   (out)
                     vwpflux         = p_fld% vwpfl_mgm   (:,:,jb),     & ! vw pseudo mometum flux   (out)
                     utflux          = p_fld% utfl_mgm    (:,:,jb),     & ! u theta flux             (out)
                     vtflux          = p_fld% vtfl_mgm    (:,:,jb)      ) ! v theta flux             (out)
    ELSE
      p_fld% uufl_mgm(:,:,jb) = 0._wp   ; p_fld% uvfl_mgm(:,:,jb) = 0._wp  ; p_fld% uwfl_mgm(:,:,jb) = 0._wp
      p_fld% vvfl_mgm(:,:,jb) = 0._wp   ; p_fld% vwfl_mgm(:,:,jb) = 0._wp
      p_fld% uupfl_mgm(:,:,jb) = 0._wp  ; p_fld% uvpfl_mgm(:,:,jb) = 0._wp ; p_fld% uwpfl_mgm(:,:,jb) = 0._wp
      p_fld% vvpfl_mgm(:,:,jb) = 0._wp  ; p_fld% vwpfl_mgm(:,:,jb) = 0._wp
      p_fld% utfl_mgm(:,:,jb) = 0._wp   ; p_fld% vtfl_mgm(:,:,jb) = 0._wp
    END IF

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_wave2grid)

    !============================ DIAGNOSTICS ===============================
    ! Project quantities of the Lagrangian ray volumes to the Eulerian grid 
    ! and calculate everything that the model does not necessarily require.
    ! The projection is done separately for GWs from the background and 
    ! from the convective sources
    !======================================================================

    IF (nrays_add_bg(jg) /= 0 .OR. gws_conv_config%n_source(jg) > 0) THEN
      call project_diagnostics(nlev           = nlev,                             & ! no. of full levels       (in)
                               i_startidx     = i_startidx,                       & ! first index of the block (in)
                               i_endidx       = i_endidx,                         & ! last index of the block  (in)
                               jray_start     = 1,                                & ! first index for ray      (in)
                               jray_end       = nrays(jg),                        & ! last index for ray       (in)
                               z              = p_metrics%z_mc(:,:,jb),           & ! full level heights       (in)
                               zhalf          = p_metrics%z_ifc(:,:,jb),          & ! half level heights       (in)
                               cellarea       = p_patch%cells%area(:,jb),         & ! area of grid cell        (in)
                               jkmin_full     = p_ray(jg)%jk_full_rtop(:,:,jb),   & ! jk closest to ray-v top  (in)
                               jkmax_full     = p_ray(jg)%jk_full_rbot(:,:,jb),   & ! jk closest to ray-v bot  (in)
                               jkmin_half     = p_ray(jg)%jk_half_rtop(:,:,jb),   & ! jk closest to ray-v top  (in)
                               jkmax_half     = p_ray(jg)%jk_half_rbot(:,:,jb),   & ! jk closest to ray-v bot  (in)
                               iexist         = p_ray(jg)%iexist(:,:,jb),         & ! existence of ray         (in)
                               jk_active      = p_ray(jg)%jk_active(:,:,jb),      & ! launch level index       (in)
                               specid         = p_ray(jg)%specid(:,:,jb),         & ! spectral id of rays      (in)
                               zray           = p_ray(jg)%z(:,:,jb),              & ! position                 (in)
                               dzray          = p_ray(jg)%dz(:,:,jb),             & ! size in z-dir            (in)
                               coslatray      = p_ray(jg)%coslat(:,:,jb),         & ! cosine latitude          (in)
                               dlatray        = p_ray(jg)%dlat(:,:,jb),           & ! meridional extent        (in)
                               dlonray        = p_ray(jg)%dlon(:,:,jb),           & ! zonal extent             (in)
                               kray           = p_ray(jg)%k(:,:,jb),              & ! horiz (lon) wavenumber   (in)
                               dkray          = p_ray(jg)%dk(:,:,jb),             & ! size in k-dir            (in)
                               lray           = p_ray(jg)%l(:,:,jb),              & ! horiz (lat) wavenumber   (in)
                               dlray          = p_ray(jg)%dl(:,:,jb),             & ! size in l-dir            (in)
                               mray           = p_ray(jg)%m(:,:,jb),              & ! vertical wavenumber      (in)
                               dmray          = p_ray(jg)%dm(:,:,jb),             & ! size in m-dir            (in)
                               dens           = p_ray(jg)%wadens(:,:,jb),         & ! wave action density      (in)
                               bvf2_full      = bvf2_full(:,:,jb),                & ! Brunt-Väisala freq**2    (in)
                               bvf2_half      = bvf2_half(:,:,jb),                & ! Brunt-Väisala freq**2    (in)
                               gammash2_full  = gammash2_full(:,:,jb),            & ! inverse pinc scale height(in)
                               gammash2_half  = gammash2_half(:,:,jb),            & ! inverse pinc scale height(in)
                               fc             = p_patch%cells%f_c(:,jb),          & ! Coriolis parameter       (in)
                               fc2            = fc2(:,jb),                        & ! Coriolis parameter**2    (in)
                               u              = u(:, :, jb),                      & ! zonal wind (full levels)      (in)
                               v              = v(:, :, jb),                      & ! meridional wind (full levels) (in)
                               theta          = theta(:,:,jb),                    & ! potential temperature    (in)
                               pmflux_e       = p_fld% pmfl_mgm_e   (:,:,jb),      & ! pseudo mometum flux E-ward      (out)
                               pmflux_w       = p_fld% pmfl_mgm_w   (:,:,jb),      & ! pseudo mometum flux W-ward      (out)
                               pmflux_s       = p_fld% pmfl_mgm_s   (:,:,jb),      & ! pseudo mometum flux S-ward      (out)
                               pmflux_n       = p_fld% pmfl_mgm_n   (:,:,jb),      & ! pseudo mometum flux N-ward      (out)
                               apmflux        = p_fld% apmfl_mgm   (:,:,jb),      & ! absolute p-momentum flux (out)
                               mflux_e        = p_fld% mfl_mgm_e   (:,:,jb),      & ! mometum flux E-ward      (out)
                               mflux_w        = p_fld% mfl_mgm_w   (:,:,jb),      & ! mometum flux W-ward      (out)
                               mflux_s        = p_fld% mfl_mgm_s   (:,:,jb),      & ! mometum flux S-ward      (out)
                               mflux_n        = p_fld% mfl_mgm_n   (:,:,jb),      & ! mometum flux N-ward      (out)
                               amflux         = p_fld% amfl_mgm    (:,:,jb),      & ! absolute momentum flux   (out)
                               waflux_u       = p_fld% wafl_mgm_u  (:,:,jb),      & ! wave action flux U-ward  (out)
                               waflux_d       = p_fld% wafl_mgm_d  (:,:,jb),      & ! wave action flux D-ward  (out)
                               waflux_e       = p_fld% wafl_mgm_e  (:,:,jb),      & ! wave action flux E-ward  (out)
                               waflux_w       = p_fld% wafl_mgm_w  (:,:,jb),      & ! wave action flux W-ward  (out)
                               waflux_s       = p_fld% wafl_mgm_s  (:,:,jb),      & ! wave action flux S-ward  (out)
                               waflux_n       = p_fld% wafl_mgm_n  (:,:,jb),      & ! wave action flux N-ward  (out)
                               ptflux_e       = p_fld% ptfl_mgm_e  (:,:,jb),      & ! theta flux E-ward        (out)
                               ptflux_w       = p_fld% ptfl_mgm_w  (:,:,jb),      & ! theta flux W-ward        (out)
                               ptflux_s       = p_fld% ptfl_mgm_s  (:,:,jb),      & ! theta flux S-ward        (out)
                               ptflux_n       = p_fld% ptfl_mgm_n  (:,:,jb),      & ! theta flux N-ward        (out)
                               aptflux        = p_fld% aptfl_mgm   (:,:,jb),      & ! absolute pot temp flux   (out)
                               energy         = p_fld% energy_mgm  (:,:,jb),      & ! GW energy                (out)
                               energy_p       = p_fld% energy_p_mgm(:,:,jb),      & ! GW potential energy      (out)
                               waction        = p_fld% action_mgm_3(:,:,jb),      & ! GW action                (out)
                               active_rays    = p_fld% active_rays_mgm(:,jb)      ) ! number of rays per cell  (out)
                             
    END IF

    !============================== DATOUT ================================
    ! Output profile of subgrid-scale GW quantities for specific lat,lon 
    ! coordinates. Mostly used for idealized tests but can be useful for 
    ! visualizing how ray volumes travel in vertical, etc.
    !======================================================================

    IF ( ldiagprof .AND. jg == 1 ) THEN

      IF (timers_level > 4) CALL timer_start(timer_msgwam_diagprof)

      IF ( ndiag_msgwam > 0 ) THEN
      IF ( ANY(jb_diag(1:ndiag_msgwam) == jb) ) THEN

      CALL debug_messages_on

      DO idiag = 1,ndiag_msgwam

        IF ( jb /= jb_diag(idiag) ) CYCLE

        jc = jc_diag(idiag)

        WRITE(message_text,'(a,i8)') 'No. of existing rays before propagate:',&
                                      COUNT(p_ray(jg)%iexist(jc,:,jb)/=0)
        CALL message('', TRIM(message_text))

        CALL datout(nlev     = nlev,                         &
                    idiag    = idiag,                        &
                    zray     = p_ray(jg)%z(jc,:,jb),         &
                    mray     = p_ray(jg)%m(jc,:,jb),         &
                    dzray    = p_ray(jg)%dz(jc,:,jb),        &
                    dmray    = p_ray(jg)%dm(jc,:,jb),        &
                    dens     = p_ray(jg)%wadens(jc,:,jb),    &
                    specid   = p_ray(jg)%specid(jc,:,jb),    &
                    zz_half  = p_metrics%z_ifc(jc,:,jb),     &
                    zz       = p_metrics%z_mc(jc,:,jb),      &
                    rho      = rho_half(jc,:,jb),            &
                    fld_u    = u(jc,:,jb),                   &
                    fld_v    = v(jc,:,jb),                   &
                    uw_wr    = p_fld% uwfl_mgm(jc,:,jb),     &
                    vw_wr    = p_fld% vwfl_mgm(jc,:,jb),     &
                    bvf2     = bvf2_half(jc,:,jb),           &
                    kd       = kd(jc,:),                     &
                    kd2      = kd2(jc,:),                    &
                    ld       = ld(jc,:),                     &
                    ld2      = ld2(jc,:),                    &
                    md       = md(jc,:),                     &
                    md2      = md2(jc,:),                    &
                    cgz_diag = cgz_diag(jc,:),               &
                    A        = A(jc,:),                      &
                    A_save   = A_save(jc,:),                 &
                    B2       = B2(jc,:),                     &
                    B2_save  = B2_save(jc,:),                &
                    mB2      = mB2(jc,:),                    &
                    mB2_save = mB2_save(jc,:),               &
                    ener     = p_fld% energy_mgm(jc,:,jb)    )

        ! Diagnostics for the wave breaking parametrization (saturation)
        IF (lsaturation) THEN
          DO jk = nlevp1,1,-1
            IF (mB2_save(jc,jk) > (alpha_sat*bvf2_half(jc,jk,jb))**2) THEN
              CALL message('', 'lsaturation=.T. ==> Saturation Parametrized')
              WRITE(message_text,'(a,3E12.4)') 'alpha, N2, height:', &
                                  alpha_sat, bvf2_half(jc,jk,jb), p_metrics%z_ifc(jc,jk,jb)
              CALL message('', TRIM(message_text))
              WRITE(message_text,'(a,3E12.4)') 'N^4*alpha^2, m^2B^2 ori, m^2B^2 reduced:', &
                                  (bvf2_half(jc,jk,jb)*alpha_sat)**2, mB2_save(jc,jk), mB2(jc,jk)
              CALL message('', TRIM(message_text))
              IF (mB2(jc,jk) > mB2_save(jc,jk)) THEN
                WRITE(message_text,'(a,i4)') 'Warning: saturation did not work well for layer:', jk
                CALL message('', TRIM(message_text))
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        IF (msg_level >= 12 ) THEN
          DO jk = nlev,1,-1
            WRITE(message_text,'(a,i6,8E12.4)') 'jk, utflux, vtflux, '// &
                                                    'uuflux, uvflux, uwflux, '// &
                                                    'vvflux, vwflux:', &
                 jk, &
                 p_fld% utfl_mgm(jc,jk,jb), p_fld% vtfl_mgm(jc,jk,jb), &
                 p_fld% uufl_mgm(jc,jk,jb), p_fld% uvfl_mgm(jc,jk,jb), p_fld% uwfl_mgm(jc,jk,jb), &
                 p_fld% vvfl_mgm(jc,jk,jb), p_fld% vwfl_mgm(jc,jk,jb)
            CALL message('', TRIM(message_text))
          ENDDO
        ENDIF

      ENDDO  ! idiag

      CALL debug_messages_off

      ENDIF ! ANY(jb_diag(1:ndiag_msgwam) == jb)
      ENDIF ! ndiag_msgwam > 0

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_diagprof)

    ENDIF ! ldiagprof .AND. jg == 1


    !======================== PROPAGATE GW FIELD ==========================
    ! Ray equations are solved in 6D (x-k space) in order to calculate the 
    ! propagation and deformation of ray volumes in phase-space. 
    !======================================================================

    IF (timers_level > 4) CALL timer_start(timer_msgwam_propagate_wave)

        CALL propagate_wave(    dt_prop    = dt_call,                         & ! dt for prop. routine call  (in)
                                jb         = jb,                              & ! block index (debug)        (in)
                                nlev       = nlev,                            & ! no. of full levels         (in)
                                i_startidx = i_startidx,                      & ! first index of the block   (in)
                                i_endidx   = i_endidx,                        & ! last index of the block    (in)
                                jray_start = 1,                               &
                                jray_end   = nrays(jg),                       &
                                z          = p_metrics%z_mc(:, :, jb),        & ! full level heights         (in)
                                dz         = p_metrics%ddqz_z_full(:, :, jb), & ! full level layer thick.    (in)
                                zhalf      = p_metrics%z_ifc(:, :, jb),       & ! half level heights         (in)
                                dzhalf     = p_metrics%ddqz_z_half(:, :, jb), & ! half level layer thick.    (in)
                                lat        = p_patch%cells%center(:, jb)%lat, & ! latitude of cell center    (in)
                                lon        = p_patch%cells%center(:, jb)%lon, & ! longitude of cell center   (in)
                                fc2        = fc2(:, jb),                      & ! Coriolis parameter**2      (in)
                                fdfdlat    = fdfdlat(:, jb),                  & ! gradiend of fc2            (in)
                                kvisc      = kvisc(:, :, jb),                 & ! Kinematic viscosity        (in)
                                bvf2       = bvf2_half(:, :, jb),             & ! Brunt-Väisala freq**2      (in)
                                gammash2   = gammash2_half(:, :, jb),         & ! inverse pinc scale height  (in)
                                u          = u_fld(:, :, jb),                 & ! zonal wind                 (in)
                                v          = v_fld(:, :, jb),                 & ! meridional wind            (in)
                                dn2dz      = dn2dz(:, :, jb),                 & ! vertical grad of bvf2      (in)
                                dg2dz      = dg2dz(:, :, jb),                 & ! vertical grad of gammash2  (in)
                                dudz       = dudz (:, :, jb),                 & ! vertical grad of u         (in)
                                dvdz       = dvdz (:, :, jb),                 & ! vertical grad of v         (in)
                                n2_hgrad   = n2_hgrad(:, :, :, jb),           & ! horizontal grads of bvf**2 (in)
                                g2_hgrad   = g2_hgrad(:, :, :, jb),           & ! horizontal grads of ga**2  (in)
                                u_hgrad    = u_hgrad(:, :, :, jb),            & ! horizontal grads of u      (in)
                                v_hgrad    = v_hgrad(:, :, :, jb),            & ! horizontal grads of v      (in)
                                dn2dz_hgrad= dn2dz_hgrad(:, :, :, jb),        & ! vert & hori grads of bvf2  (in)
                                dg2dz_hgrad= dg2dz_hgrad(:, :, :, jb),        & ! vert & hori grads of ga**2 (in)
                                dudz_hgrad = dudz_hgrad (:, :, :, jb),        & ! vert & hori grads of u     (in)
                                dvdz_hgrad = dvdz_hgrad (:, :, :, jb),        & ! vert & hori grads of v     (in)
                                lonray     = p_ray(jg)%lon(:, :, jb),         & ! ray position lon           (inout)
                                dlonray    = p_ray(jg)%dlon(:, :, jb),        & ! ray size in lon            (inout)
                                latray     = p_ray(jg)%lat(:, :, jb),         & ! ray position lat           (inout)
                                dlatray    = p_ray(jg)%dlat(:, :, jb),        & ! ray size in lat            (inout)
                                zray       = p_ray(jg)%z(:, :, jb),           & ! position                   (inout)
                                dzray      = p_ray(jg)%dz(:, :, jb),          & ! size in z-dir              (inout)
                                coslatray  = p_ray(jg)%coslat(:, :, jb),      & ! cosine of ray lat          (inout)
                                kray       = p_ray(jg)%k(:, :, jb),           & ! horiz (lon) wavenumber     (in)
                                dkray      = p_ray(jg)%dk(:, :, jb),          & ! size in k-dir              (in)
                                lray       = p_ray(jg)%l(:, :, jb),           & ! horiz (lat) wavenumber     (in)
                                dlray      = p_ray(jg)%dl(:, :, jb),          & ! size in l-dir              (in)
                                mray       = p_ray(jg)%m(:, :, jb),           & ! vertical wavenumber        (inout)
                                dmray      = p_ray(jg)%dm(:, :, jb),          & ! size in m-dir              (inout)
                                dens       = p_ray(jg)%wadens(:, :, jb),      & ! wave action density        (inout)
                                iexist     = p_ray(jg)%iexist(:, :, jb),      & ! existence of ray           (inout)
                                specid     = p_ray(jg)%specid(:, :, jb),      & ! spectral id of rays        (inout)
                                jk_active  = p_ray(jg)%jk_active(:, :, jb),   & ! launch level index         (in)
                                jr_last    = p_ray(jg)%jr_last(:, :, jb)      ) ! last-launched ray index    (inout)

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_propagate_wave)

    !======================== Diagnostics ==================================
    ! WA diagnostic
    !=======================================================================

    CALL idx_rayedge( nlev,i_startidx,i_endidx,1,nrays(jg),& ! (in)
                      p_metrics%z_mc     (:,:,jb),        & ! (in)
                      p_metrics%z_ifc    (:,:,jb),        & ! (in)
                      p_ray(jg)%iexist   (:,:,jb),        & ! (in)
                      p_ray(jg)%specid   (:,:,jb),        & ! (in)
                      p_ray(jg)%jk_active(:,:,jb),        & ! (in)
                      p_ray(jg)%z        (:,:,jb),        & ! (in)
                      p_ray(jg)%dz       (:,:,jb),        & ! (in)
                      p_ray(jg)%jk_full_rtop(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_full_rbot(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_half_rtop(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_half_rbot(:,:,jb)      ) ! (inout)

    CALL project_action(nlev        = nlev,                             & !
                        i_startidx  = i_startidx,                       & !
                        i_endidx    = i_endidx,                         & !
                        jray_start  = 1,                                & !
                        jray_end    = nrays(jg),                        & !
                        z           = p_metrics%z_mc(:,:,jb),           & !
                        zhalf       = p_metrics%z_ifc(:,:,jb),          & !
                        cellarea    = p_patch%cells%area(:,jb),         & !
                        zray        = p_ray(jg)%z(:,:,jb),              & !
                        dzray       = p_ray(jg)%dz(:,:,jb),             & !
                        coslatray   = p_ray(jg)%coslat(:,:,jb),         & !
                        dlatray     = p_ray(jg)%dlat(:,:,jb),           & !
                        dlonray     = p_ray(jg)%dlon(:,:,jb),           & !
                        dkray       = p_ray(jg)%dk(:,:,jb),             & !
                        dlray       = p_ray(jg)%dl(:,:,jb),             & !
                        dmray       = p_ray(jg)%dm(:,:,jb),             & !
                        dens        = p_ray(jg)%wadens(:,:,jb),         & !
                        iexist      = p_ray(jg)%iexist(:,:,jb),         & !
                        specid      = p_ray(jg)%specid(:,:,jb),         & !
                        jk_active   = p_ray(jg)%jk_active(:,:,jb),      & !
                        jkmin_full  = p_ray(jg)%jk_full_rtop(:,:,jb),   & !
                        jkmax_full  = p_ray(jg)%jk_full_rbot(:,:,jb),   & !
                        jkmin_half  = p_ray(jg)%jk_half_rtop(:,:,jb),   & !
                        jkmax_half  = p_ray(jg)%jk_half_rbot(:,:,jb),   & !
                        action      = p_fld%action_mgm_4(:,:,jb),       & !
                        diag        = 0,                                & !
                        jg          = jg                                ) !

  ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL


  !================ AVERAGE GW FIELD IN HORIZONTAL ======================
  ! ...
  ! TODO: 
  !======================================================================

#ifndef __msgwam1d
  IF (lhsmooth) THEN

    IF (timers_level > 4) CALL timer_start(timer_msgwam_smooth_hori)

    ! Background source related fields
    DO iter_sm = 1, ABS(nhsmooth)
      CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_fld% uwfl_mgm, p_fld% vwfl_mgm,  &
        &                           p_fld% uufl_mgm, p_fld% uvfl_mgm, p_fld% vvfl_mgm)
      CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_fld% uwpfl_mgm, p_fld% vwpfl_mgm,  &
        &                          p_fld% uupfl_mgm, p_fld% uvpfl_mgm, p_fld% vvpfl_mgm)
      CALL sync_patch_array_mult(SYNC_C, p_patch, 2, p_fld% utfl_mgm, p_fld% vtfl_mgm)
      !
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%uwfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%vwfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%uufl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%uvfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%vvfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%uwpfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%vwpfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%uupfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%uvpfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%vvpfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%utfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%vtfl_mgm)
    ENDDO

    IF (nhsmooth > 0) THEN
      ! Fluxes that do not affect the tendencies are smoothed only once or not ever,
      ! depending on the sign of nhsmooth.

      CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_fld% apmfl_mgm    , &
        &                                            p_fld% amfl_mgm     , &
        &                                            p_fld% aptfl_mgm    , &
        &                                            p_fld% energy_mgm   , &
        &                                            p_fld% energy_p_mgm )
      CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_fld% action_mgm_1 , &
                                                     p_fld% action_mgm_2 , &
                                                     p_fld% action_mgm_3 , &
                                                     p_fld% action_mgm_4 , &
                                                     p_fld% action_mgm_5 )
      !
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%apmfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%amfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%aptfl_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%energy_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%energy_p_mgm)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%action_mgm_1)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%action_mgm_2)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%action_mgm_3)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%action_mgm_4)
      CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,  p_fld%action_mgm_5)

      IF ( lcalc_flux_4dir_bg(jg) ) THEN
        CALL sync_patch_array_mult(SYNC_C, p_patch, 4, p_fld% mfl_mgm_e  , &
          &                                            p_fld% mfl_mgm_w  , &
          &                                            p_fld% mfl_mgm_n  , &
          &                                            p_fld% mfl_mgm_s  )
        CALL sync_patch_array_mult(SYNC_C, p_patch, 4, p_fld% pmfl_mgm_e  , &
          &                                            p_fld% pmfl_mgm_w  , &
          &                                            p_fld% pmfl_mgm_n  , &
          &                                            p_fld% pmfl_mgm_s  )
        CALL sync_patch_array_mult(SYNC_C, p_patch, 4, p_fld% ptfl_mgm_e , &
          &                                            p_fld% ptfl_mgm_w , &
          &                                            p_fld% ptfl_mgm_n , &
          &                                            p_fld% ptfl_mgm_s )
        !
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%mfl_mgm_e)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%mfl_mgm_w)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%mfl_mgm_n)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%mfl_mgm_s)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%pmfl_mgm_e)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%pmfl_mgm_w)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%pmfl_mgm_n)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlevp1,p_fld%pmfl_mgm_s)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,p_fld%ptfl_mgm_e)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,p_fld%ptfl_mgm_w)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,p_fld%ptfl_mgm_n)
        CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev,p_fld%ptfl_mgm_s)
      ENDIF

    ENDIF  ! nhsmooth > 0

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_smooth_hori)

  ENDIF ! lhsmooth
#endif


  !============================= TENDENCY ===============================
  ! Calculate horizontal wind tendencies based on the vertical divergence
  ! of pseudo-momentum fluxes. The total tendency and the tendency due to 
  ! convective GWs is calculated separately (the latter for diagnostics).
  ! TODO: calculate the 3D divergence of momentum fluxes in case of  
  !       horizontal propagation
  !======================================================================

  IF (timers_level > 4) CALL timer_start(timer_msgwam_tendency)

  CALL tendency(p_patch,p_metrics,p_int_state,rho,temp,theta,p_fld)

  IF (timers_level > 4) CALL timer_stop(timer_msgwam_tendency)


  !====================== SYNCHRONIZE GW FIELD ==========================
  ! ...
  ! TODO: 
  !======================================================================

#ifndef __msgwam1d
  IF (timers_level > 4) CALL timer_start(timer_msgwam_sync_wave)
  CALL sync_wave(p_ray(jg),p_patch)
  IF (timers_level > 4) CALL timer_stop(timer_msgwam_sync_wave)
#endif


  !========================== SPLIT GW FIELD ============================
  ! ...
  ! TODO: 
  !======================================================================

  IF ( imethod_split > 0 .OR. imethod_merge > 0 ) THEN
    IF (timers_level > 4) CALL timer_start(timer_msgwam_split_merge)

    IF (imethod_merge < 10) THEN   ! operate separately for each source

      IF ( gws_conv_config%n_source(jg) > 0 ) THEN
        nrays_vacate = nrays_add_cv(jg)*nlaunch_max_cv
        CALL split_merge_volume( p_ray(jg), p_spl(jg), p_rwork(jg), p_patch,            &
          &                      p_mgmgrid(jg)%c_nb, p_mgmgrid(jg)%ncol_int,            &
          &                      p_metrics%z_mc, p_metrics%z_ifc, nrays_vacate,         &
          &                      fc2, bvf2_half, gammash2_half,                         &
          &                      specid_start = specid_offset_cv(jg)+1,                 &
          &                      specid_end   = specid_offset_cv(jg)+nrays_add_cv(jg),  &
          &                      jray_start = jray_offset_cv(jg)+1,                     &
          &                      jray_end   = jray_offset_cv(jg)+nrays_cv(jg) )
      END IF

      IF (nrays_add_bg(jg) /= 0) THEN
        nrays_vacate = nrays_add_bg(jg)*nlaunch_max_bg
        CALL split_merge_volume( p_ray(jg), p_spl(jg), p_rwork(jg), p_patch,            &
          &                      p_mgmgrid(jg)%c_nb, p_mgmgrid(jg)%ncol_int,            &
          &                      p_metrics%z_mc, p_metrics%z_ifc, nrays_vacate,         &
          &                      fc2, bvf2_half, gammash2_half,                         &
          &                      specid_start = specid_offset_bg(jg)+1,                 &
          &                      specid_end   = specid_offset_bg(jg)+nrays_add_bg(jg),  &
          &                      jray_start = jray_offset_bg(jg)+1,                     &
          &                      jray_end   = jray_offset_bg(jg)+nrays_bg(jg) )
      END IF

    ELSE   ! operate at once

      nrays_vacate =  nrays_add_cv(jg)*nlaunch_max_cv  &
        &           + nrays_add_bg(jg)*nlaunch_max_bg
      CALL split_merge_volume( p_ray(jg), p_spl(jg), p_rwork(jg), p_patch,      &
        &                      p_mgmgrid(jg)%c_nb, p_mgmgrid(jg)%ncol_int,      &
        &                      p_metrics%z_mc, p_metrics%z_ifc, nrays_vacate,   &
        &                      fc2, bvf2_half, gammash2_half,                   &
        &                      specid_start = 1,                                &
        &                      specid_end   = nrays_add(jg),                    &
        &                      jray_start = 1,                                  &
        &                      jray_end   = nrays(jg) )

    END IF

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_split_merge)
  END IF

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

    !======================== Diagnostics ==================================
    ! WA diagnostic
    !=======================================================================

    CALL idx_rayedge( nlev,i_startidx,i_endidx,1,nrays(jg),& ! (in)
                      p_metrics%z_mc     (:,:,jb),        & ! (in)
                      p_metrics%z_ifc    (:,:,jb),        & ! (in)
                      p_ray(jg)%iexist   (:,:,jb),        & ! (in)
                      p_ray(jg)%specid   (:,:,jb),        & ! (in)
                      p_ray(jg)%jk_active(:,:,jb),        & ! (in)
                      p_ray(jg)%z        (:,:,jb),        & ! (in)
                      p_ray(jg)%dz       (:,:,jb),        & ! (in)
                      p_ray(jg)%jk_full_rtop(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_full_rbot(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_half_rtop(:,:,jb),     & ! (inout)
                      p_ray(jg)%jk_half_rbot(:,:,jb)      ) ! (inout)

    CALL project_action(nlev        = nlev,                             & !
                        i_startidx  = i_startidx,                       & !
                        i_endidx    = i_endidx,                         & !
                        jray_start  = 1,                                & !
                        jray_end    = nrays(jg),                        & !
                        z           = p_metrics%z_mc(:,:,jb),           & !
                        zhalf       = p_metrics%z_ifc(:,:,jb),          & !
                        cellarea    = p_patch%cells%area(:,jb),         & !
                        zray        = p_ray(jg)%z(:,:,jb),              & !
                        dzray       = p_ray(jg)%dz(:,:,jb),             & !
                        coslatray   = p_ray(jg)%coslat(:,:,jb),         & !
                        dlatray     = p_ray(jg)%dlat(:,:,jb),           & !
                        dlonray     = p_ray(jg)%dlon(:,:,jb),           & !
                        dkray       = p_ray(jg)%dk(:,:,jb),             & !
                        dlray       = p_ray(jg)%dl(:,:,jb),             & !
                        dmray       = p_ray(jg)%dm(:,:,jb),             & !
                        dens        = p_ray(jg)%wadens(:,:,jb),         & !
                        iexist      = p_ray(jg)%iexist(:,:,jb),         & !
                        specid      = p_ray(jg)%specid(:,:,jb),         & !
                        jk_active   = p_ray(jg)%jk_active(:,:,jb),      & !
                        jkmin_full  = p_ray(jg)%jk_full_rtop(:,:,jb),   & !
                        jkmax_full  = p_ray(jg)%jk_full_rbot(:,:,jb),   & !
                        jkmin_half  = p_ray(jg)%jk_half_rtop(:,:,jb),   & !
                        jkmax_half  = p_ray(jg)%jk_half_rbot(:,:,jb),   & !
                        action      = p_fld%action_mgm_5(:,:,jb),       & !
                        diag        = 0,                                & !
                        jg          = jg                                ) !

  ENDDO  ! jb loop for diagnostic

  !======================== RE-GRID GW FIELD ============================
  ! ...
  ! TODO: 
  !======================================================================

#ifndef __msgwam1d
  IF ( imethod_split <= 0 .AND. imethod_merge <= 0 ) THEN
    IF (timers_level > 4) CALL timer_start(timer_msgwam_split_merge)
    !IF (timers_level > 4) CALL timer_start(timer_msgwam_regrid_wave)
    CALL regrid_wave( p_patch, p_gridinfo4ray(jg), p_metrics%z_ifc, p_metrics%z_mc,  &
      &               nrays(jg), p_ray(jg) )
    IF (timers_level > 4) CALL timer_stop(timer_msgwam_split_merge)
    !IF (timers_level > 4) CALL timer_stop(timer_msgwam_regrid_wave)
  END IF
#endif

  CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev, p_fld%action_mgm_1)
  CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev, p_fld%action_mgm_2)
  CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev, p_fld%action_mgm_3)
  CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev, p_fld%action_mgm_4)
  CALL smooth_hori(p_patch,p_gridinfo4ray(jg),p_int_state,nlev, p_fld%action_mgm_5)


  !======================== TENDENCY LIMITER ============================
  ! Tendency limiter to stabilize high-top runs. This limiter is taken 
  ! from mo_nwp_gw_interface.f90. Normally we would not like to use it as 
  ! MS-GWaM fluxes with a direct wave-meanflow interaction (+ wave 
  ! breaking) should not be out of realistic range. We still keep the 
  ! option...
  !======================================================================

  ! hard-coded tendency checker :  to be removed later
  CALL test_blowup_uvt(p_patch, p_fld%ddt_u_gwd_mgm, p_fld%ddt_v_gwd_mgm, p_fld%ddt_t_gwd_mgm,  &
    &                  lim_tend_tmp, lim_tend_tmp*0.1_wp, '(checkpoint 1)')

  IF (llimittend) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk

      !CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      !  & i_startidx, i_endidx, rl_start, rl_end)
      !
      ! This limiter is taken from mo_nwp_gw_interface.f90. Normally we would not like to use it 
      ! as MS-GWaM fluxes with a direct wave-meanflow interaction (+ wave breaking) should not be out 
      ! of realistic range. We still keep the option...
      !DO jk = 1, nlev
      !  DO jc = i_startidx, i_endidx
      !    p_fld% ddt_u_gwd_mgm(jc,jk,jb) = &
      !      SIGN(MIN(lim_tend,ABS(p_fld% ddt_u_gwd_mgm(jc,jk,jb))),p_fld% ddt_u_gwd_mgm(jc,jk,jb))
      !    p_fld% ddt_v_gwd_mgm(jc,jk,jb) = &
      !      SIGN(MIN(lim_tend,ABS(p_fld% ddt_v_gwd_mgm(jc,jk,jb))),p_fld% ddt_v_gwd_mgm(jc,jk,jb))
      !  ENDDO ! jc
      !ENDDO ! jk

      p_fld%ddt_u_gwd_mgm(:,:,jb) = MAX(-lim_tend_tmp, MIN(lim_tend_tmp,  &
        &                           p_fld%ddt_u_gwd_mgm(:,:,jb)))
      p_fld%ddt_v_gwd_mgm(:,:,jb) = MAX(-lim_tend_tmp, MIN(lim_tend_tmp,  &
        &                           p_fld%ddt_v_gwd_mgm(:,:,jb)))
      p_fld%ddt_u_gwd_pmom_mgm(:,:,jb) = MAX(-lim_tend_tmp, MIN(lim_tend_tmp,  &
        &                                p_fld%ddt_u_gwd_pmom_mgm(:,:,jb)))
      p_fld%ddt_v_gwd_pmom_mgm(:,:,jb) = MAX(-lim_tend_tmp, MIN(lim_tend_tmp,  &
        &                                p_fld%ddt_v_gwd_pmom_mgm(:,:,jb)))

    ENDDO ! jb
!$OMP END DO
!$OMP END PARALLEL

  ENDIF ! llimittend


  !======================== TENDENCY DIAGNOSTIC ============================
  ! Print out maximum of absolute of tendencies
  ! TODO: one would need a syncronization and an mpi_allreduce operation 
  !       here in order to get full domain diagnostics. This is here only a 
  !       preliminary printout to have a first feeling how often the 
  !       tendencies are limited via llimittend = .true.
  !=========================================================================
  IF (msg_level >= 12) THEN
    maxindexu = MAXLOC(ABS(p_fld%ddt_u_gwd_mgm))
    maxindexv = MAXLOC(ABS(p_fld%ddt_v_gwd_mgm))
    WRITE(message_text,'(a,E12.4,a,i4)') 'MAXABS of GW u-tendencies:', & 
      &  MAXVAL(ABS(p_fld% ddt_u_gwd_mgm)), ' at level:', maxindexu(2)
    CALL message('', TRIM(message_text))
    WRITE(message_text,'(a,E12.4,a,i4)') 'MAXABS of GW v-tendencies:', & 
      &  MAXVAL(ABS(p_fld% ddt_v_gwd_mgm)), ' at level:', maxindexv(2)
    CALL message('', TRIM(message_text))
  ENDIF

  ! Record number for subroutine datout (ldiagprof=.T.)
  iout_msgwam = iout_msgwam + 1

#ifndef __msgwam1d
  CALL message('gwdrag_msgwam', 'Transient MS-GWaM finished')
#else
  CALL message('gwdrag_msgwam', 'Transient MS-GWaM-1D finished')
#endif

END SUBROUTINE gwdrag_msgwam
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE fieldsgrads( p_patch, p_int_state, rl_start, rl_end, nlev, z, zhalf, dz, dzhalf, wgtfac_c, &
                        rho, temp, u, v, fc2, fdfdlat, rho_half, u_fld, v_fld, &
                        bvf2_full, bvf2_half, gammash2_full, gammash2_half, kvisc, &
                        dn2dz, dg2dz, dudz, dvdz, &
                        n2_hgrad, g2_hgrad, u_hgrad, v_hgrad, &
                        dn2dz_hgrad, dg2dz_hgrad, dudz_hgrad, dvdz_hgrad )

  TYPE(t_patch), TARGET, INTENT(INOUT) :: p_patch                 ! grid/patch info.
  TYPE(t_int_state),     INTENT(IN)    :: p_int_state             ! interpolation struct
  INTEGER,               INTENT(IN)    :: rl_start, rl_end        ! iteration bounds
  INTEGER,               INTENT(IN)    :: nlev                    ! number of full levels
  REAL(wp),              INTENT(IN)    :: z(:, :, :)              ! height of full levels
  REAL(wp),              INTENT(IN)    :: zhalf(:, :, :)          ! height of half levels
  REAL(wp),              INTENT(IN)    :: dz(:, :, :)             ! thickness of full levels
  REAL(vp),              INTENT(IN)    :: dzhalf(:, :, :)         ! thickness of half levels
  REAL(vp),              INTENT(IN)    :: wgtfac_c(:, :, :)       !
  REAL(wp),              INTENT(IN)    :: rho(:, :, :)            ! density at full levels
  REAL(wp),              INTENT(IN)    :: temp(:, :, :)           ! temperature at full levels
  REAL(wp),              INTENT(IN)    :: u(:, :, :)              ! zonal wind (full levels)
  REAL(wp),              INTENT(IN)    :: v(:, :, :)              ! meridional wind (full levels)
  REAL(wp),              INTENT(  OUT) :: fc2(:, :)               ! fc^2
  REAL(wp),              INTENT(  OUT) :: fdfdlat(:, :)           ! d(fc^2)/dlat
  REAL(wp),              INTENT(  OUT) :: rho_half(:, :, :)       ! density at half levels
  REAL(wp),              INTENT(  OUT) :: u_fld(:, :, :)          ! zonal wind
  REAL(wp),              INTENT(  OUT) :: v_fld(:, :, :)          ! meridional wind
  REAL(wp),              INTENT(  OUT) :: bvf2_full(:, :, :)      ! Brunt-Väisälä freq**2 at full levels
  REAL(wp),              INTENT(  OUT) :: bvf2_half(:, :, :)      ! Brunt-Väisälä freq**2 at half levels
  REAL(wp),              INTENT(  OUT) :: gammash2_full(:, :, :)  ! inverse pinc scale height**2 (i.e. gamma**2) at full levels
  REAL(wp),              INTENT(  OUT) :: gammash2_half(:, :, :)  ! inverse pinc scale height**2 (i.e. gamma**2) at half levels
  REAL(wp),              INTENT(  OUT) :: kvisc(:, :, :)          ! Kinematic viscosity
  REAL(wp),              INTENT(  OUT) :: dn2dz(:, :, :)          ! vertical B-V freq gradients at half levels
  REAL(wp),              INTENT(  OUT) :: dg2dz(:, :, :)          ! vertical gradients of inverse pinc scale height at half levels
  REAL(wp),              INTENT(  OUT) :: dudz (:, :, :)          ! vertical (radial) wind u gradients at half levels
  REAL(wp),              INTENT(  OUT) :: dvdz (:, :, :)          ! vertical (radial) wind v gradients at half levels
  REAL(wp),              INTENT(  OUT) :: n2_hgrad(:, :, :, :)    ! zonal / meridional gradients of BVF**2
  REAL(wp),              INTENT(  OUT) :: g2_hgrad(:, :, :, :)    ! zonal / meridional gradients of gamma**2 at full levels
  REAL(wp),              INTENT(  OUT) :: u_hgrad(:, :, :, :)     ! zonal / meridional wind u gradients at full levels
  REAL(wp),              INTENT(  OUT) :: v_hgrad(:, :, :, :)     ! zonal / meridional wind v gradients at full levels
  REAL(wp),              INTENT(  OUT) :: dn2dz_hgrad(:, :, :, :) ! vert & zonal / meridional gradients of BVF**2
  REAL(wp),              INTENT(  OUT) :: dg2dz_hgrad(:, :, :, :) ! vert & zonal / meridional gradients of gamma**2 at full levels
  REAL(wp),              INTENT(  OUT) :: dudz_hgrad(:, :, :, :)  ! vert & zonal / meridional wind u gradients at full levels
  REAL(wp),              INTENT(  OUT) :: dvdz_hgrad(:, :, :, :)  ! vert & zonal / meridional wind v gradients at full levels

  REAL(wp)                             :: temp_half(nproma, p_patch%nlevp1, p_patch%nblks_c)
  REAL(vp)                             :: hgrad(2, nproma, p_patch%nlevp1, p_patch%nblks_c)
  INTEGER                              :: jk, jc                  ! vertical and block indices
  INTEGER                              :: i_startidx, i_endidx    ! first and last block index
  INTEGER                              :: i_startblk, i_endblk    ! first and last block index
  INTEGER                              :: rl_end_h1               ! iteration bound with 1 halo line
  INTEGER                              :: i_endblk_h1             ! last block index including 1 halo line
  INTEGER                              :: jb                      ! block index
  INTEGER                              :: nlevp1
  REAL(wp)                             :: inv_dz
  REAL(wp)                             :: wgtfac_c_jkm1           ! 1. - wgtfac_c
  REAL(wp)                             :: Hrho                    ! density scale height
  REAL(wp),    PARAMETER               :: grav_o_cpd = grav / cpd !

  IF (msg_level >= 12) CALL message('fieldsgrads', 'MS-GWaM: prepare resolved fields and grads')

  !----------------------------------------------------------------------
  ! Purpose:
  !         Calculate all necessary resolved fields appearing as input to MS-GWaM
  ! Method:
  !         -- smooth wind and temperature (OMP jb loop)
  !         -- calculate Brunt-Väisälä frequency (OMP jb loop)
  !         -- calculate pinc scale height correction term gamma (OMP jb loop)
  !         -- calculate vertcial gradients (OMP jb loop)
  !         -- calculate the horizontal wind gradients (Green-Gauss method)
  !
  !----------------------------------------------------------------------

  ! get indices of first and last block
  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

#ifdef __msgwam1d
  rl_end_h1   = rl_end
  i_endblk_h1 = i_endblk
#else
  rl_end_h1   = rl_end - 1   ! to include one inner line of halo
! rl_end_h1   = rl_end - 2   ! to include one pair (inner + outer) of halo line
  i_endblk_h1 = p_patch%cells%end_block(rl_end_h1)
#endif

  ! set number of vertical levels + 1
  nlevp1 = nlev+1

  IF (is_plane_torus)  fdfdlat(:, :) = 0._wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)

  DO jb = i_startblk, i_endblk

    ! get cell indices
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

    ! Initialize Coriolis frequency
    fc2(:, jb) = p_patch%cells%f_c(:, jb)**2
    IF (.NOT. is_plane_torus)  fdfdlat(:, jb) = (4._wp*earth_angular_velocity)  &
      &            * p_patch%cells%f_c(:, jb) * p_mgmgrid(jg)%coslat(:, jb)

    ! Interpolate the resolved fields to half levels
    DO jk = 2, nlev
      DO jc = i_startidx, i_endidx
        wgtfac_c_jkm1 = 1._wp - wgtfac_c(jc, jk, jb)
        rho_half (jc, jk, jb) = wgtfac_c(jc, jk, jb) * rho (jc, jk ,jb) & 
                                     + wgtfac_c_jkm1 * rho (jc, jk - 1, jb)
        temp_half(jc, jk, jb) = wgtfac_c(jc, jk, jb) * temp(jc, jk, jb) & 
                                     + wgtfac_c_jkm1 * temp(jc, jk - 1, jb)
        u_fld(jc, jk, jb) = wgtfac_c(jc, jk, jb) * u(jc, jk, jb) &
                                 + wgtfac_c_jkm1 * u(jc, jk - 1, jb)
        v_fld(jc, jk, jb) = wgtfac_c(jc, jk, jb) * v(jc, jk, jb) &
                                 + wgtfac_c_jkm1 * v(jc, jk - 1, jb)
      ENDDO  ! jc
    ENDDO  ! jk

    DO jc = i_startidx, i_endidx
      rho_half (jc, 1, jb) = rho (jc, 1, jb)
      rho_half (jc, nlevp1, jb) = rho (jc, nlev, jb)
      temp_half(jc, 1, jb) = temp(jc, 1, jb)
      temp_half(jc, nlevp1, jb) = temp(jc, nlev, jb)
    ENDDO

    DO jc = i_startidx, i_endidx
      ! dA/dz = 0 to prevent extrapolation
      u_fld(jc, 1     , jb) = u_fld(jc, 2   , jb)
      u_fld(jc, nlevp1, jb) = u_fld(jc, nlev, jb)
      v_fld(jc, 1     , jb) = v_fld(jc, 2   , jb)
      v_fld(jc, nlevp1, jb) = v_fld(jc, nlev, jb)
    ENDDO

    ! Smoothing over 2*nsmooth+1 points
    IF (lsmoothb) THEN
      CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .TRUE., temp_half(:, :, jb))
      CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .TRUE., u_fld(:, :, jb))
      CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .TRUE., v_fld(:, :, jb))
      DO jc = i_startidx, i_endidx
        ! dA/dz = 0 to prevent extrapolation
        u_fld(jc, 1     , jb) = u_fld(jc, 2   , jb)
        u_fld(jc, nlevp1, jb) = u_fld(jc, nlev, jb)
        v_fld(jc, 1     , jb) = v_fld(jc, 2   , jb)
        v_fld(jc, nlevp1, jb) = v_fld(jc, nlev, jb)
      ENDDO
    ENDIF

    ! Fields on full levels
    DO jk = 2, nlev - 1

      DO jc = i_startidx, i_endidx

        ! Diag printout
        IF (msg_level >= 15) THEN
          WRITE(message_text, '(a,3i6,E12.4)') 'jc, jk, jb, temp(jc,jk,jb):', &
                jc, jk, jb, temp(jc, jk, jb)
          CALL message('', TRIM(message_text))
        ENDIF

        ! Calculate Brunt-Väisälä freq (full levels)
        ! negative pot temperature gradient is not allowed,
        ! i.e. in that case bvf2 is set to bvf2_min
!       bvf2_full(jc, jk, jb) = MAX(bvf2_min, grav / temp(jc, jk, jb) * (grav_o_cpd + &
!             (temp(jc, jk - 1, jb) - temp(jc, jk + 1, jb)) / (dzhalf(jc, jk, jb) + dzhalf(jc, jk + 1, jb))))
        bvf2_full(jc, jk, jb) = MAX(bvf2_min, grav / temp(jc, jk, jb) * (grav_o_cpd + &
              (temp_half(jc, jk, jb) - temp_half(jc, jk + 1, jb)) / dz(jc, jk, jb)))

        ! Calculate gamma (half levels): inverse pseudo incompressible scaleheight
        IF (lcorrlongwaves) THEN
          ! Density scaleheight assuming "locally" isothermal atmosphere
          Hrho = rd * temp(jc, jk, jb) / grav
          ! Assuming "locally" isothermal atmosphere
          gammash2_full(jc, jk, jb) = (facgamma / Hrho)**2
          ! facgamma is a namelist parameter with the following typical values:
          ! Pseudo-incompressible correction:         facgamma ~= 0.214
          ! Anelastic correction:                     facgamma  = 0.5
          ! Further decrease effects from long waves: facgamma  > 0.5
          ! Diag printout
          IF (msg_level >= 15) THEN
            WRITE(message_text, '(a,3i6,2E12.4)') 'jc, jk, jb, bvf2_full(jc,jk,jb), 1/gammash2_full(jc,jk,jb):', &
                  jc, jk, jb, bvf2_full(jc, jk, jb), 1._wp / gammash2_full(jc, jk, jb)
            CALL message('', TRIM(message_text))
          ENDIF
        ENDIF

      ENDDO ! jc

    ENDDO ! jk

    ! fields on half levels
    DO jk = 2, nlev

      DO jc = i_startidx, i_endidx

        ! Diag printout
        IF (msg_level >= 15) THEN
          WRITE(message_text, '(a,3i6,E12.4)') 'jc, jk, jb, temp_half(jc, jk, jb):', &
                jc, jk, jb, temp_half(jc, jk, jb)
          CALL message('', TRIM(message_text))
        ENDIF

        ! Calculate Brunt-Väisälä freq (half levels)
        ! negative pot temperature gradient is not allowed,
        ! i.e. in that case bvf2 is set to bvf2_min
        bvf2_half(jc, jk, jb) = MAX(bvf2_min, grav / temp_half(jc, jk, jb) * (grav_o_cpd + &
              (temp_half(jc, jk - 1, jb) - temp_half(jc, jk + 1, jb)) / (dz(jc, jk - 1, jb) + dz(jc, jk, jb))))

        ! Calculate gamma (half levels): inverse pseudo incompressible scaleheight
        IF (lcorrlongwaves) THEN
          ! Density scaleheight assuming "locally" isothermal atmosphere
          Hrho = rd * temp_half(jc, jk, jb) / grav
          ! Assuming "locally" isothermal atmosphere
          gammash2_half(jc, jk, jb) = (facgamma / Hrho)**2
          ! facgamma is a namelist parameter with the following typical values:
          ! Pseudo-incompressible correction:         facgamma ~= 0.214
          ! Anelastic correction:                     facgamma  = 0.5
          ! Further decrease effects from long waves: facgamma  > 0.5
        ENDIF

        ! Diag printout
        IF (msg_level >= 15) THEN
          WRITE(message_text, '(a,3i6,2E12.4)') 'jc, jk, jb, bvf2_half(jc, jk, jb), 1 / gammash2_half(jc, jk, jb):', &
                jc, jk, jb, bvf2_half(jc, jk, jb), 1._wp / gammash2_half(jc, jk, jb)
          CALL message('', TRIM(message_text))
        ENDIF

        ! Calculate kinematic viscosity profile
        kvisc(jc, jk, jb) = dyn_visc_sutherland(temp_half(jc, jk, jb)) / rho_half(jc, jk, jb)

      ENDDO ! jc

    ENDDO ! jk

    IF (.NOT. lmvisc)  kvisc(:,:,jb) = 0.

    ! Values at vertical boundaries
    DO jc = i_startidx, i_endidx
      ! dA/dz = 0 to prevent extrapolation
      bvf2_half(jc, 1, jb)          = bvf2_half(jc, 2, jb)
      bvf2_half(jc, nlevp1, jb)     = bvf2_half(jc, nlev, jb)
      gammash2_half(jc, 1, jb)      = gammash2_half(jc, 2, jb)
      gammash2_half(jc, nlevp1, jb) = gammash2_half(jc, nlev, jb)

      bvf2_full(jc, 1, jb)        = bvf2_full(jc, 2, jb)
      bvf2_full(jc, nlev, jb)     = bvf2_full(jc, nlev-1, jb)
      gammash2_full(jc, 1, jb)    = gammash2_full(jc, 2, jb)
      gammash2_full(jc, nlev, jb) = gammash2_full(jc, nlev-1, jb)

      kvisc(jc, 1, jb)    = kvisc(jc, 2, jb)
      kvisc(jc, nlevp1, jb) = kvisc(jc, nlev, jb)
    ENDDO ! jc

  ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL

#ifndef __msgwam1d
  ! Synchronize input fields of horizontal gradient calculations
  CALL sync_patch_array_mult(SYNC_C1, p_patch, 4, u_fld, v_fld, bvf2_half, gammash2_half)
#endif

  ! The vertical gradients calculated below will be used in the equation for Dm/Dt,
  ! after interpolation to the ray-volume center.
  ! It has turned out that the use of vertically interpolated gradient fields
  ! helps to obtain numerical stable results, compared to the use of gradient fields
  ! directly calculated between the top and bottom of a ray volume.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk_h1   ! including one inner line of halo, assuming that sync has
                                    ! already been made above with SYNC_C*, for the 4 variables
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk_h1, &
          & i_startidx, i_endidx, rl_start, rl_end_h1)

    DO jk = 2, nlev
      DO jc = i_startidx, i_endidx

        inv_dz = 1._wp / (dz(jc, jk-1, jb) + dz(jc, jk, jb))  ! two full layers

        ! vertical gradients of N^2, Gamma^2
        dn2dz(jc, jk, jb) = (bvf2_half(jc, jk-1, jb) - bvf2_half(jc, jk+1, jb)) * inv_dz
        dg2dz(jc, jk, jb) = (gammash2_half(jc, jk-1, jb) - gammash2_half(jc, jk+1, jb)) * inv_dz

        ! vertical (radial) gradients of the horizontal wind u, v
        dudz(jc, jk, jb) = (u_fld(jc, jk-1, jb) - u_fld(jc, jk+1, jb)) * inv_dz
        dvdz(jc, jk, jb) = (v_fld(jc, jk-1, jb) - v_fld(jc, jk+1, jb)) * inv_dz

      ENDDO ! jc
    ENDDO ! jk

    ! Values at vertical boundaries
    dn2dz(:, 1, jb) = 0._wp  ;  dn2dz(:, nlevp1, jb) = 0._wp
    dg2dz(:, 1, jb) = 0._wp  ;  dg2dz(:, nlevp1, jb) = 0._wp
    dudz (:, 1, jb) = 0._wp  ;  dudz (:, nlevp1, jb) = 0._wp
    dvdz (:, 1, jb) = 0._wp  ;  dvdz (:, nlevp1, jb) = 0._wp

  ENDDO ! jb
!$OMP END DO
!$OMP END PARALLEL

#ifndef __msgwam1d
  ! get the gradients in lat and lon
  ! the green gauss algorithms are parallelized with OMP internally
  hgrad(:,:,:,:) = 0.
  CALL grad_green_gauss_cell(p_cc = u_fld,           & ! cell centered variable to interpolate
                             ptr_patch = p_patch,    & ! patch info
                             ptr_int = p_int_state,  & ! interpolation info
                             p_grad = hgrad,         & ! horizontal gradients
                             opt_slev = 2,           & ! lowest vertical level
                             opt_elev = nlev,        & ! highest vertical level
                             opt_rlstart = rl_start, &
                             opt_rlend   = rl_end    )
  u_hgrad = REAL(hgrad, KIND=wp)

  hgrad(:,:,:,:) = 0.
  CALL grad_green_gauss_cell(p_cc = v_fld,           & ! cell centered variable to interpolate
                             ptr_patch = p_patch,    & ! patch info
                             ptr_int = p_int_state,  & ! interpolation info
                             p_grad = hgrad,         & ! horizontal gradients
                             opt_slev = 2,           & ! lowest vertical level
                             opt_elev = nlev,        & ! highest vertical level
                             opt_rlstart = rl_start, &
                             opt_rlend   = rl_end    )
  v_hgrad = REAL(hgrad, KIND=wp)

  hgrad(:,:,:,:) = 0.
  CALL grad_green_gauss_cell(p_cc = bvf2_half,       & ! cell centered variable to interpolate
                             ptr_patch = p_patch,    & ! patch info
                             ptr_int = p_int_state,  & ! interpolation info
                             p_grad = hgrad,         & ! horizontal gradients
                             opt_slev = 2,           & ! lowest vertical level
                             opt_elev = nlev,        & ! highest vertical level
                             opt_rlstart = rl_start, &
                             opt_rlend   = rl_end    )
  n2_hgrad = REAL(hgrad, KIND=wp)

  hgrad(:,:,:,:) = 0.
  CALL grad_green_gauss_cell(p_cc = gammash2_half,   & ! cell centered variable to interpolate
                             ptr_patch = p_patch,    & ! patch info
                             ptr_int = p_int_state,  & ! interpolation info
                             p_grad = hgrad,         & ! horizontal gradients
                             opt_slev = 2,           & ! lowest vertical level
                             opt_elev = nlev,        & ! highest vertical level
                             opt_rlstart = rl_start, &
                             opt_rlend   = rl_end    )
  g2_hgrad = REAL(hgrad, KIND=wp)

  hgrad(:,:,:,:) = 0.
  CALL grad_green_gauss_cell(p_cc = dudz,            & ! cell centered variable to interpolate
                             ptr_patch = p_patch,    & ! patch info
                             ptr_int = p_int_state,  & ! interpolation info
                             p_grad = hgrad,         & ! horizontal gradients
                             opt_slev = 2,           & ! lowest vertical level
                             opt_elev = nlev,        & ! highest vertical level
                             opt_rlstart = rl_start, &
                             opt_rlend   = rl_end    )
  dudz_hgrad = REAL(hgrad, KIND=wp)

  hgrad(:,:,:,:) = 0.
  CALL grad_green_gauss_cell(p_cc = dvdz,            & ! cell centered variable to interpolate
                             ptr_patch = p_patch,    & ! patch info
                             ptr_int = p_int_state,  & ! interpolation info
                             p_grad = hgrad,         & ! horizontal gradients
                             opt_slev = 2,           & ! lowest vertical level
                             opt_elev = nlev,        & ! highest vertical level
                             opt_rlstart = rl_start, &
                             opt_rlend   = rl_end    )
  dvdz_hgrad = REAL(hgrad, KIND=wp)

  hgrad(:,:,:,:) = 0.
  CALL grad_green_gauss_cell(p_cc = dn2dz,           & ! cell centered variable to interpolate
                             ptr_patch = p_patch,    & ! patch info
                             ptr_int = p_int_state,  & ! interpolation info
                             p_grad = hgrad,         & ! horizontal gradients
                             opt_slev = 2,           & ! lowest vertical level
                             opt_elev = nlev,        & ! highest vertical level
                             opt_rlstart = rl_start, &
                             opt_rlend   = rl_end    )
  dn2dz_hgrad = REAL(hgrad, KIND=wp)

  hgrad(:,:,:,:) = 0.
  CALL grad_green_gauss_cell(p_cc = dg2dz,           & ! cell centered variable to interpolate
                             ptr_patch = p_patch,    & ! patch info
                             ptr_int = p_int_state,  & ! interpolation info
                             p_grad = hgrad,         & ! horizontal gradients
                             opt_slev = 2,           & ! lowest vertical level
                             opt_elev = nlev,        & ! highest vertical level
                             opt_rlstart = rl_start, &
                             opt_rlend   = rl_end    )
  dg2dz_hgrad = REAL(hgrad, KIND=wp)

  IF (.NOT. is_plane_torus) THEN   ! convert to dA/dlon/cos(lat) or dA/dlat
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          u_hgrad (:,jc,jk,jb) = u_hgrad (:,jc,jk,jb)*grid_sphere_radius
          v_hgrad (:,jc,jk,jb) = v_hgrad (:,jc,jk,jb)*grid_sphere_radius
          n2_hgrad(:,jc,jk,jb) = n2_hgrad(:,jc,jk,jb)*grid_sphere_radius
          g2_hgrad(:,jc,jk,jb) = g2_hgrad(:,jc,jk,jb)*grid_sphere_radius
          dudz_hgrad (:,jc,jk,jb) = dudz_hgrad (:,jc,jk,jb)*grid_sphere_radius
          dvdz_hgrad (:,jc,jk,jb) = dvdz_hgrad (:,jc,jk,jb)*grid_sphere_radius
          dn2dz_hgrad(:,jc,jk,jb) = dn2dz_hgrad(:,jc,jk,jb)*grid_sphere_radius
          dg2dz_hgrad(:,jc,jk,jb) = dg2dz_hgrad(:,jc,jk,jb)*grid_sphere_radius
        ENDDO
      ENDDO
    ENDDO
  END IF

  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      &                i_startidx, i_endidx, rl_start, rl_end)
    DO jc = i_startidx, i_endidx
      ! dA/dz = 0 to prevent extrapolation
      u_hgrad (:,jc,1     ,jb) = u_hgrad (:,jc,2   ,jb)
      u_hgrad (:,jc,nlevp1,jb) = u_hgrad (:,jc,nlev,jb)
      v_hgrad (:,jc,1     ,jb) = v_hgrad (:,jc,2   ,jb)
      v_hgrad (:,jc,nlevp1,jb) = v_hgrad (:,jc,nlev,jb)
      n2_hgrad(:,jc,1     ,jb) = n2_hgrad(:,jc,2   ,jb)
      n2_hgrad(:,jc,nlevp1,jb) = n2_hgrad(:,jc,nlev,jb)
      g2_hgrad(:,jc,1     ,jb) = g2_hgrad(:,jc,2   ,jb)
      g2_hgrad(:,jc,nlevp1,jb) = g2_hgrad(:,jc,nlev,jb)
      dudz_hgrad (:,jc,1     ,jb) = 0.
      dudz_hgrad (:,jc,nlevp1,jb) = 0.
      dvdz_hgrad (:,jc,1     ,jb) = 0.
      dvdz_hgrad (:,jc,nlevp1,jb) = 0.
      dn2dz_hgrad(:,jc,1     ,jb) = 0.
      dn2dz_hgrad(:,jc,nlevp1,jb) = 0.
      dg2dz_hgrad(:,jc,1     ,jb) = 0.
      dg2dz_hgrad(:,jc,nlevp1,jb) = 0.
    ENDDO
  ENDDO
#endif

END SUBROUTINE fieldsgrads
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE init_gw_orretal(nlev,i_startidx,i_endidx,jray_start,jray_end,&
                           mdatetime,jklaunch,nlaunch_max,z,zhalf,dz,dzhalf,&
                           cellarea,lon,lat,coslat,rhol,bvfl2,lat_prof_bw,lat_prof_bs,lat_prof,&
                           lonray,dlonray,latray,dlatray,zray,dzray,coslatray,&
                           kray,dkray,lray,dlray,mray,dmray,dens,&
                           iexist,specid,jk_active,jr_last,prec)
  TYPE(datetime), INTENT(IN)    :: mdatetime
  INTEGER,        INTENT(IN)    :: nlev                 ! no. of model levels
  INTEGER,        INTENT(IN)    :: i_startidx           ! first index of the block
  INTEGER,        INTENT(IN)    :: i_endidx             ! last index of the block
  INTEGER,        INTENT(IN)    :: jray_start
  INTEGER,        INTENT(IN)    :: jray_end
  INTEGER,        INTENT(IN)    :: jklaunch             ! launch level index
  INTEGER,        INTENT(IN)    :: nlaunch_max
  REAL(wp),       INTENT(IN)    :: z(:,:)               ! full level heights
  REAL(wp),       INTENT(IN)    :: zhalf(:,:)           ! half level heights
  REAL(wp),       INTENT(IN)    :: dz(:,:)              ! full level thickness
  REAL(vp),       INTENT(IN)    :: dzhalf(:,:)          ! half level thickness
  REAL(wp),       INTENT(IN)    :: cellarea(:)          ! area of grid cell [m**2]
  REAL(wp),       INTENT(IN)    :: lon(:)               ! longitude of cell-center
  REAL(wp),       INTENT(IN)    :: lat(:)               ! latitude of cell-center
  REAL(wp),       INTENT(IN)    :: coslat(:)            ! cosine of latitude
  REAL(wp),       INTENT(IN)    :: rhol(:)              ! density at launch level
  REAL(wp),       INTENT(IN)    :: bvfl2(:)             ! B-V freq**2 at launch level
  REAL(wp),       INTENT(IN)    :: lat_prof_bw(:)       ! latitudinal factor
  REAL(wp),       INTENT(IN)    :: lat_prof_bs(:)       ! latitudinal factor
  REAL(wp),       INTENT(INOUT) :: lat_prof(:)          ! latitudinal factor
  REAL(wp),       INTENT(IN)    :: prec(:)              ! Total precip
  INTEGER,        INTENT(INOUT) :: iexist(:,:)
  INTEGER,        INTENT(INOUT) :: specid(:,:)          ! spectral id of rays
  INTEGER,        INTENT(INOUT) :: jk_active(:,:)
  INTEGER,        INTENT(INOUT) :: jr_last(:,:)
  REAL(wp),       INTENT(INOUT) :: lonray(:,:), dlonray(:,:)
  REAL(wp),       INTENT(INOUT) :: latray(:,:), dlatray(:,:)
  REAL(wp),       INTENT(INOUT) :: zray(:,:),   dzray(:,:)
  REAL(wp),       INTENT(INOUT) :: coslatray(:,:)
  REAL(wp),       INTENT(INOUT) :: mray(:,:),   dmray(:,:)
  REAL(wp),       INTENT(INOUT) :: kray(:,:),   dkray(:,:)
  REAL(wp),       INTENT(INOUT) :: lray(:,:),   dlray(:,:)
  REAL(wp),       INTENT(INOUT) :: dens(:,:)

  INTEGER  :: jk, jc
  INTEGER  :: jray, jray_offset
  INTEGER  :: jrayprev
  INTEGER  :: vecnew(nproma,nrays_add_bg(jg)) ! vector to retain indices of new rays
  INTEGER  :: khdim    ! dimension of total horizontal wavenumber spectral array

  INTEGER  :: ic, im, ikh, iomega          ! indices for spectral elements
  INTEGER  :: iazi                         ! index of the azimuth angle
  INTEGER  :: ispec                        ! index for spectral location
  INTEGER  :: cdim                         ! dimension of phase speed spectral array
  INTEGER  :: nlaunch                      ! number of launching for each spec element
  INTEGER  :: il                           ! index for launching
  REAL(wp) :: phi(iazidim)                 ! azimuth angle
  REAL(wp) :: dphi                         ! delta azimuth angle
  REAL(wp) :: cin(mdim)                    ! c launch intrinsic phase speed element
  REAL(wp) :: dcin(mdim)                   ! spectrum element size in c (intrinsic)
  REAL(wp) :: m(mdim)                      ! m launch vertical wavenumber element
  REAL(wp) :: dm(mdim)                     ! spectrum element size in m
  REAL(wp) :: kh(mdim,omegadim)            ! kh launch horizontal wavenumber element (kh=sqrt(k^2+l^2))
  REAL(wp) :: dkh(mdim,omegadim)           ! spectrum element size in kh 
  REAL(wp) :: cgz                          ! vertical group velocity at ray position
  REAL(wp) :: fct                          ! factor = rho/bvf
  REAL(wp) :: fnorm(iazidim)               ! normalisation factor (A)
  REAL(wp) :: fluxin(mdim,iazidim)         ! spectral momentum flux at each azimuth
                                           ! as a function of (c_in,phi)
                                           ! this is what is used in the Orr et al. scheme!
  REAL(wp) :: flux1(omegadim,mdim,iazidim) ! spectral momentum flux at each azimuth
                                           ! as a function of (omega,c,phi)
  REAL(wp) :: flux2(omegadim,mdim,iazidim) ! spectral momentum flux at each azimuth
                                           ! as a function of (m,omega,phi)
  REAL(wp) :: flux3(omegadim,mdim,iazidim) ! spectral momentum flux at each azimuth
                                           ! as a function of (k_h,m,phi)
                                           ! remember that omegadim = khdim
  REAL(wp) :: fluxray                      ! spectral momentum flux at each azimuth
                                           ! as a function of (k,l,m)
  REAL(wp) :: pu(nlev,iazidim)             ! total momentum flux
  REAL(wp) :: ccin4, ccin2, ccin3          ! intrisic phase velocity
  REAL(wp) :: bvfl4, bvfl3, bvfl           ! help vars for B-V freq
  REAL(wp) :: mstar, z50s
  REAL(wp) :: latdeg
  REAL(wp) :: grid_coslat
  REAL(wp) :: dlon, dlat
  REAL(wp) :: gauss, fluxlaun
  REAL(wp) :: zray_d_prev
  INTEGER  :: jkmid(0:nlaunch_max)
  REAL(wp) :: zlaunch(nlaunch_max)
  REAL(wp) :: zb_ghost
  REAL(wp) :: tmp3(3)
  REAL(wp) :: tfact
  REAL(wp) :: distance
  INTEGER  :: year_m4k
  REAL(wp), PARAMETER :: tphase_0 = 8520._wp/8760._wp  ! phase for 00 UTC 22 Dec. (for non-leap years)
  REAL(wp), PARAMETER :: ndaysoffset(12) = REAL(  &
    &                      (/0,31,59,90,120,151,181,212,243,273,304,334/), wp )

  IF (msg_level >= 12) CALL message('init_gw_orretal', &
    &                   'MS-GWaM: initialize Desaubies type background GW field')

  !----------------------------------------------------------------------
  ! Purpose: 
  !         Calculate launch spectrum for "background" GWs based on the 
  !         Desaubies spectra. For the time being this source in MS-GWaM 
  !         is meant to represent in a very simple manner 
  !         -- either a simplified all non-orographic GW sources (if 
  !            convective sources turned off) 
  !         -- or GWs emmitted by jets/fronts (if convective sources 
  !            turned on)
  !
  ! Method:
  !         -- PART 1) calculate mom. flux launch spectrum in (c,phi) 
  !                    based on mo_gwd_wms.f90
  !         -- PART 2) translate the launch mom. flux spectra (c,phi) 
  !                    into a launch wave action density spectra (k,l,m)
  !         -- see Orr et al. (2010) and Scinocca et al. (2003)
  !----------------------------------------------------------------------

  ! Note: the launch mom. flux spectrum calculation is based on the assumption
  ! that the POSITIVE FREQUENCY BRANCH IS USED, which also means that a NEGATIVE
  ! VERTICAL WAVENUMBER spectrum will be used. This results in upward propagating 
  ! gravity waves (cgz > 0) injected at launch level.

  jray_offset = jray_start - 1

  ! Calculate time dependent launch flux magnitudes. Currently, only the Gregorian 
  ! calendar with leap years is considered.
  tfact = ( ndaysoffset(mdatetime%date%month) + REAL(mdatetime%date%day-1,wp) )*24._wp  &
    &     + REAL(mdatetime%time%hour,wp) + REAL(mdatetime%time%minute,wp)/60._wp        &
    &     + ( REAL(mdatetime%time%second,wp)                                            &
    &         + 0.001_wp*REAL(mdatetime%time%ms,wp) )/3600._wp
  ! to convert the type to normal integer with a moderate size
  year_m4k = INT( MOD(mdatetime%date%year, INT(4000, KIND=KIND(mdatetime%date%year))) )
  IF ( ( MOD(year_m4k, 4) == 0 .AND. MOD(year_m4k, 100) /= 0 )  &
    &  .OR. MOD(year_m4k, 400) == 0 ) THEN   ! leap year
    IF (mdatetime%date%month > 2)  tfact = tfact + 24._wp
    tfact = tfact/8784._wp   ! hours for 366 days
  ELSE
    tfact = tfact/8760._wp   ! hours for 365 days
  END IF
  tfact = 0.5_wp*(COS((tfact - tphase_0)*(2._wp*pi)) + 1._wp)

  DO jc = i_startidx, i_endidx
    lat_prof(jc) = lat_prof_bs(jc) + (lat_prof_bw(jc)-lat_prof_bs(jc))*tfact
  ENDDO

  ! ----------------------------------------------------------------------------
  ! PART 1: calculate launch spectrum in (c,phi) space (based on mo_gwd_wms.f90)
  ! ----------------------------------------------------------------------------

  ! Set same dimension for m and c spectral elements, we suppose keeping cdim 
  ! (and not simply using mdim everywhere) keeps the code more understandable
  cdim = mdim

  ! mstar
  mstar = 2._wp*pi/2000._wp

  ! Unstreched phase velocity spectrum
  dcin(1) = (cimax-cimin)/REAL(MAX(cdim-1,1),wp)
  dcin(2:) = dcin(1)
  DO ic = 1,cdim
    cin(ic) = cimin + dcin(1)*REAL(ic-1,wp)
  ENDDO

  ! Define of azimuth angles
  dphi = 2._wp*pi/REAL(iazidim,wp)
  DO iazi=1,iazidim
    phi(iazi) = (REAL(iazi,wp)-1._wp)*dphi
  ENDDO

  DO jc = i_startidx, i_endidx

#ifdef __test_awaypole
if(abs(lat(jc))>85.*deg2rad) cycle
#endif

    ! skip the generation if the latitude is ouside the bounds
    IF (abs(lat(jc)) > msgw_source_limit * deg2rad) CYCLE


    ! Horizontal propagation test
    ! Horizontal localization of the source: skip calculations 
    ! if ltest_hprop=.true. and cell is out of the area defined 
    ! by lonrmin, lonrmax, latrmin, latrmax
    IF (ltest_hprop) THEN
      IF (lon(jc)*rad2deg > lonrmax .OR. &
          lon(jc)*rad2deg < lonrmin .OR. &
          lat(jc)*rad2deg > latrmax .OR. &
          lat(jc)*rad2deg < latrmin) THEN
        CYCLE
      ELSE
        ! Write output file for horizontal propagation test
        OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
        WRITE(987,'(a,2x,F14.5,2x,F14.5)') &
              'GW source lon, lat coordinates:', &
               lon(jc)*rad2deg, lat(jc)*rad2deg
        CLOSE(987)
      ENDIF
    ENDIF

    ! Skip calculations for this cell if local 
    ! latitudinal profile is zero
    IF (lat_prof(jc) == 0._wp) CYCLE

    distance = SQRT(cellarea(jc))
    dlat = distance/grid_sphere_radius
    IF (is_plane_torus) THEN
      grid_coslat = 1._wp
      dlon = dlat
    ELSE
      grid_coslat = coslat(jc)
!      dlon = 2._wp*ASIN(SIN(dlat/2._wp)/coslat(jc))  ! argument of asin can be > 1
      dlon = dlat/coslat(jc)
    END IF

    ! Define factor rho/bvf
    bvfl = SQRT(bvfl2(jc))
    fct = rhol(jc)/bvfl

    ! Set mom. flux spectra (c,phi) at launch level: Eq. (25) of Scinocca 2003 (not 
    ! including the 'A' component), and with U-Uo=0. 
    fluxin(:,:) = 0._wp
    IF(nslope==1) THEN
      ! s=1 case
      bvfl4=bvfl2(jc)**2
      DO iazi=1,iazidim
        DO ic=1,cdim
          ccin4=(mstar*cin(ic))**4
          fluxin(ic,iazi)=fct*bvfl4*cin(ic)/(bvfl4+ccin4)
        ENDDO
      ENDDO
    ELSEIF(nslope==-1) THEN
      ! s=-1 case
      DO iazi=1,iazidim
        DO ic=1,cdim
          ccin2=(mstar*cin(ic))**2
          fluxin(ic,iazi)=fct*bvfl2(jc)*cin(ic)/(bvfl2(jc)+ccin2)
        ENDDO
      ENDDO
    ELSEIF(nslope==0) THEN
      ! s=0 case
      bvfl3=bvfl2(jc)*bvfl
      DO iazi=1,iazidim
        DO ic=1,cdim
          ccin3=(mstar*cin(ic))**3
          fluxin(ic,iazi)=fct*bvfl3*cin(ic)/(bvfl3+ccin3)
        ENDDO
      ENDDO
    ENDIF

    ! Integrate (zflux x dx)
    pu(jklaunch,:) = 0._wp
    DO iazi=1,iazidim
      DO ic=1,cdim
        pu(jklaunch,iazi)=pu(jklaunch,iazi)+fluxin(ic,iazi)*dcin(ic)
      ENDDO
    ENDDO
    IF (msg_level >= 12) THEN
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG1 jklaunch, eastward  flux:', jklaunch, pu(jklaunch,1)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG1 jklaunch, northward flux:', jklaunch, pu(jklaunch,2)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG1 jklaunch, westward  flux:', jklaunch, pu(jklaunch,3)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG1 jklaunch, southward flux:', jklaunch, pu(jklaunch,4)
      CALL message('', TRIM(message_text))
    ENDIF

    ! Calculate normalization factor
    fluxlaun=gfluxlaun
    DO iazi=1,iazidim
      fnorm(iazi)=fluxlaun/pu(jklaunch,iazi)
    ENDDO

    ! Introduce latitudional dependence of the launch mom. flux. If lozpr=.true. then 
    ! increase launch flux over tropics fnorm:=A in Scinocca 2003 (independent of height)
    IF (lozpr) THEN
      IF (ngauss==1) THEN
        fluxlaun=gfluxlaun*(1.0_wp+MIN(0.5_wp,gcoeff*prec(jc)))
        DO iazi=1,iazidim
          fnorm(iazi)=fluxlaun/pu(jklaunch,iazi)
        ENDDO
      ELSEIF (ngauss==2) THEN
        latdeg=lat(jc)*rad2deg
        gauss=ggaussb*EXP((-latdeg*latdeg)/(2._wp*ggaussa*ggaussa))
        fluxlaun=(1.0_wp+gauss)*gfluxlaun
        DO iazi=1,iazidim
          fnorm(iazi)=fluxlaun/pu(jklaunch,iazi)
        ENDDO
      ELSEIF (ngauss==4) THEN
        ! Set latitudinal dependence to optimize stratospheric winds
        z50s=-50.0_wp
        latdeg=lat(jc)*rad2deg-z50s
        gauss=ggaussb*EXP((-latdeg*latdeg)/(2._wp*ggaussa*ggaussa))
        fluxlaun=(1.0_wp+gauss)*gfluxlaun
        DO iazi=1,iazidim
          fnorm(iazi)=fluxlaun/pu(jklaunch,iazi)
        ENDDO
      ENDIF
    ENDIF

    ! Normalization factor in case of latitude dependent launch flux
    IF ( lfluxlatsimple .OR. flux_tropics(1) >= 0. ) THEN
      ! Overwrite lozpr options
      fluxlaun=lat_prof(jc)
      DO iazi=1,iazidim
        fnorm(iazi)=fluxlaun/pu(jklaunch,iazi)
      ENDDO
    ENDIF

    ! Renormalize each spectral element
    DO iazi=1,iazidim
      DO ic=1,cdim
        fluxin(ic,iazi)=fnorm(iazi)*fluxin(ic,iazi)
      ENDDO
    ENDDO

    ! Check if normalization leads to an integral of fluxlaunch
    IF (msg_level >= 12) THEN
      pu(jklaunch,:) = 0._wp
      DO iazi=1,iazidim
        DO ic=1,cdim
          pu(jklaunch,iazi)=pu(jklaunch,iazi)+fluxin(ic,iazi)*dcin(ic)
        ENDDO
      ENDDO
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG2 jklaunch, eastward  flux:', jklaunch, pu(jklaunch,1)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG2 jklaunch, northward flux:', jklaunch, pu(jklaunch,2)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG2 jklaunch, westward  flux:', jklaunch, pu(jklaunch,3)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG2 jklaunch, southward flux:', jklaunch, pu(jklaunch,4)
      CALL message('', TRIM(message_text))
    ENDIF

    ! --------------------------------------------------------------------
    ! PART 2: translate the launch mom. flux spectra (c,phi) into a launch 
    !         wave action density spectra (k,l,m). Calculate wave action 
    !         density N(k,l,m) based on fluxes rho*F(c,phi). See steps of 
    !         Scinocca (2003) in a reverse order from Eq.(25) to Eq.(14)
    !         and then define wave action density N = rho*F/kh/cgz/dk/dl/dm
    ! --------------------------------------------------------------------

    ! --------------------------------------------------------------------
    ! Step 1: rho*F(c,phi) --> rho*F(c,omega,phi)
    ! --------------------------------------------------------------------

    ! Calculate rho*F(c,omega,phi)
    flux1(:,:,:) = 0._wp
    DO iazi=1,iazidim
      DO ic=1,cdim
        DO iomega=1,omegadim
          flux1(iomega,ic,iazi) = fluxin(ic,iazi) * sp_omega_l(iomega)
        ENDDO 
      ENDDO 
    ENDDO 

    ! Check if the integral of flux1 is still OK
    IF (msg_level >= 12) THEN
      pu(jklaunch,:) = 0._wp
      DO iazi=1,iazidim
        DO ic=1,cdim
          DO iomega=1,omegadim
            pu(jklaunch,iazi)=pu(jklaunch,iazi)+flux1(iomega,ic,iazi)*domega_l(iomega)*dcin(ic)
          ENDDO
        ENDDO
      ENDDO
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG3 jklaunch, eastward  flux:', jklaunch, pu(jklaunch,1)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG3 jklaunch, northward flux:', jklaunch, pu(jklaunch,2)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG3 jklaunch, westward  flux:', jklaunch, pu(jklaunch,3)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG3 jklaunch, southward flux:', jklaunch, pu(jklaunch,4)
      CALL message('', TRIM(message_text))
    ENDIF

    ! --------------------------------------------------------------------
    ! Step 2: rho*F(c,omega,phi) --> rho*F(m,omega,phi)
    ! --------------------------------------------------------------------

    ! Compute rho*F(m,omega,phi)
    flux2(:,:,:) = 0._wp
    DO iazi=1,iazidim
      DO im=1,mdim
        ! To avoid confusions still use different indices but eventually 
        ! they are the same
        ic = im
        ! Calculate m based on the hydrostatic dispersion relation
        m(im) = -bvfl/cin(ic)
        dm(im) = m(im)**2/bvfl*dcin(ic)
        DO iomega=1,omegadim
          flux2(iomega,im,iazi) = flux1(iomega,ic,iazi) * bvfl/m(im)**2
        ENDDO
      ENDDO
    ENDDO

    ! Check if the integral of flux2 is still OK
    IF (msg_level >= 12) THEN
      pu(jklaunch,:) = 0._wp
      DO iazi=1,iazidim
        DO im=1,mdim
          DO iomega=1,omegadim
            pu(jklaunch,iazi)=pu(jklaunch,iazi)+flux2(iomega,im,iazi)*domega_l(iomega)*dm(im)
          ENDDO
        ENDDO
      ENDDO
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG4 jklaunch, eastward  flux:', jklaunch, pu(jklaunch,1)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG4 jklaunch, northward flux:', jklaunch, pu(jklaunch,2)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG4 jklaunch, westward  flux:', jklaunch, pu(jklaunch,3)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG4 jklaunch, southward flux:', jklaunch, pu(jklaunch,4)
      CALL message('', TRIM(message_text))
    ENDIF

    ! --------------------------------------------------------------------
    ! Step 3: rho*F(m,omega,phi) --> rho*F(m,k_h,phi)
    ! --------------------------------------------------------------------

    ! Define dimension for the total horizontal wavenumber k_h
    khdim = omegadim

    ! Calculate rho*F(m,k_h,phi)
    flux3(:,:,:) = 0._wp
    DO iazi=1,iazidim
      DO im=1,mdim
        DO ikh=1,khdim
          ! To avoid confusions still use different indices but eventually 
          ! they are the same
          iomega = ikh
          ! Calculate kh(m,omega) --> necessary for Step 5
          kh(im,iomega) = ABS(omega_l(iomega)*m(im)/bvfl)
          dkh(im,iomega) = ABS(domega_l(iomega)*m(im)/bvfl)
          flux3(ikh,im,iazi) = flux2(iomega,im,iazi) * (-bvfl/m(im))
        ENDDO
      ENDDO
    ENDDO

    ! Check if the integral of flux3 is still OK
    IF (msg_level >= 12) THEN
      pu(jklaunch,:) = 0._wp
      DO iazi=1,iazidim
        DO im=1,mdim
          DO ikh=1,khdim
            iomega = ikh
            pu(jklaunch,iazi)=pu(jklaunch,iazi)+flux3(ikh,im,iazi)*dkh(im,iomega)*dm(im)
          ENDDO
        ENDDO
      ENDDO
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG5 jklaunch, eastward  flux:', jklaunch, pu(jklaunch,1)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG5 jklaunch, northward flux:', jklaunch, pu(jklaunch,2)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG5 jklaunch, westward  flux:', jklaunch, pu(jklaunch,3)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG5 jklaunch, southward flux:', jklaunch, pu(jklaunch,4)
      CALL message('', TRIM(message_text))
    ENDIF

    ! ----------------------------------------------------------------------
    ! Step 4: rho*F(m,k_h,phi) --> rho*F(k,l,m)
    ! ----------------------------------------------------------------------
    ! Note: This step is not done here explicitly. It is implicitly included 
    ! in the next step (Step 5) by taking into account that 
    ! rho*F(k,l,m) = rho*F(k_h,m,phi)/k_h

    ! --------------------------------------------------------------------
    ! Step 5: rho*F(m,k_h,phi) --> N(k,l,m)
    ! --------------------------------------------------------------------

    zb_ghost = zhalf(jc,jklaunch) - depth_ghost_bg

    ispec = specid_offset_bg(jg)   ! == 0
    jray  = jray_offset
    vecnew(jc,:) = 0
    pu(jklaunch,:) = 0._wp
    DO iazi=1,iazidim ! loop over azimuth angles
      DO ikh=1,khdim  ! loop over horizontal wavenumber
        DO im=1,mdim  ! loop over vertical wavenumber

          ! Index for omega elements (same as for k_h elements)
          iomega = ikh

          ! Counter for spectral element ID
          ispec = ispec + 1

          ! RULE: Launch consecutive ray volumes before the previously launched volume
          !       crosses over z = z_src, the level at which we have the source spectrum
          !       with specified amount of fluxes to be parameterized. For the ray volumes
          !       launched below this level, we do not change (k,l,m) (and thus dz) during
          !       the propagation until their centers arrive at z = z_src in order to have
          !       correct flux amounts and (k,l,m) at this level.

          ! For this, we set a ghost zone from z_src below, with its depth H_g. The timing
          ! of launching a new volume is when the last launched volume has fully crossed
          ! over z = z_src - H_g, provided that H_g > cgz * dt_add to keep the
          ! above rule. (H_g = 3 km might be large enough when dt_add ~ 5 min.)
          ! (See Fig.4 of Bölöni et al. (2021) and corresponding explanations.)

          ! A) For the first launching: one volume is launched such that its top is at
          !    zb_ghost (= z_src - H_g).
          ! B) If the last launched volume has not fully crossed over zb_ghost,
          !    do not launch any.
          ! C) If the last volume has fully crossed over zb_ghost, one or more volumes
          !    are launched consecutively below the former volumes until a launched volume
          !    covers z = zb_ghost. This becomes the new 'last' volume.
          ! D) Rarely, 1) if the last volume crosses over z_src in one dt_add
          !    (i.e., H_g is not large enough for some very fast waves),
          !    or 2) if the last volume has been removed in some reason (saturation),
          !    new volumes are launched consecutively, with the top of the first one at
          !    z - dz/2. (There might be no obvious way in these cases..)

          ! Ray index of this spectral element at the last launch time
          jrayprev = jr_last(jc,ispec)

          !-------------------------------------------------------------------------------------
          ! Step 5.1.: locate the height of ray volume center-points (i.e. zlaunch(:), jkmid(:)) 
          ! for the new ray volumes to be launched depending on the above cases A), B), C), D)
          !-------------------------------------------------------------------------------------

          IF (jrayprev <= 0) THEN   ! rare cases binded

            IF (jrayprev == 0) THEN
              ! A) This spectral element is not present from the last launch time,
              !    i.e. this is the very first launch of this spectral element
              nlaunch = 1
              zlaunch(1) = zb_ghost - 0.5_wp*dzlaunch_def
              jkmid(0) = jklaunch   ! lower limit of jkmid for searching
            ELSE  ! jrayprev == -1
              ! D) This spectral element was launched last time but the
              !    corresponding ray volume has propagated between two launch 
              !    times so far that its bottom is now above the launch level z_src.
              !    In this case several ray volumes are launched with this spectral
              !    ID until the ghost layer is filled.
              nlaunch = nlaunch_max
              DO il = 1, nlaunch
                zlaunch(il) = zhalf(jc,jklaunch) - REAL(il,KIND=wp)*dzlaunch_def
              ENDDO 
              jkmid(0) = jklaunch ! lower limit of jkmid for searching
            END IF

          ! This spectral element was also launched at the last launch time. We
          ! know that from the fact that specid(jc,jrayprev) is negative, which 
          ! is how it will be given a value later in this subroutine. 
          ELSE  ! specid(jc,jrayprev) == -ispec

            ! Calculate height of the bottom of the previously emmitted ray 
            ! volume with this spectral ID
            zray_d_prev = zray(jc,jrayprev) - 0.5_wp*dzray(jc,jrayprev)

            ! B) The ray volume's bottom is still below the ghost layer's bottom
            !    so we do not launch this spectral element this time.
            IF (zray_d_prev < zb_ghost)  CYCLE  ! im (ispec)

            ! C) The ray volume is completely above the the ghost layer's bottom 
            !    so we launch a new one from this spectral element.
            nlaunch = MIN(nlaunch_max,INT((zray_d_prev - zb_ghost)/dzlaunch_def)+1)
            DO il = 1, nlaunch
              zlaunch(il) = zray(jc,jrayprev) - REAL(il,KIND=wp)*dzlaunch_def
            ENDDO 
            jkmid(0) = iexist(jc,jrayprev)   ! lower limit of jkmid for searching

          ENDIF

          ! Calculate jkmid(:) for each zlaunch(:)
          DO il = 1, nlaunch
            jkmid(il) = nlev+1
            DO jk = jkmid(il-1), nlev ! downward search
              IF (z(jc,jk) < zlaunch(il)) THEN
                jkmid(il) = jk  ;  EXIT
              END IF
            ENDDO
          ENDDO

          !-------------------------------------------------------------
          ! Step 5.2.: initialize new ray volumes with their properties
          !-------------------------------------------------------------

          DO il = 1, nlaunch

            ! Ray counter
            jray = jray + 1
 
            ! Do not overwrite existing rays
            DO WHILE (iexist(jc,jray) /= 0)
              jray = jray + 1
              ! Below a check for strange cases indicating problems with the
              ! dimensions for allocating the ray types from the early
              ! development phase. This could be removed at some point.
              IF (jray > jray_end) CALL finish (TRIM('init_gw_orretal'), &
                 'Something is strange in the background GW source: jray>nrays')
            ENDDO

            ! Initialize horizontal position
            lonray(jc,jray) = lon(jc)
            latray(jc,jray) = lat(jc)
            coslatray(jc,jray) = grid_coslat

            ! Initialize horizontal extent. So far simply an extent is taken
            ! of which the square gives the area of the triangle cell. The
            ! extent in lon, lat is calculated from the distance by inverting
            ! the Haversine formula.
            dlonray(jc,jray) = dlon
            dlatray(jc,jray) = dlat

            ! Background GWs are active above their launch level constant in
            ! time: jklaunch (calculated from the namelist parameter plaunch) 
            jk_active(jc,jray) = jklaunch

            ! Retain spectral id of the newly initialized ray. The negative sign 
            ! denotes that the ray volume is "new", i.e. it has been launched at 
            ! the previous launch time. This is turned into its absolute value in 
            ! subroutine propagate_wave when the ray volumes becomes "old" (i.e. 
            ! not launched at the previous launch time but sometime before).
            specid(jc,jray) = -ispec
 
            ! Vector to retain indices of new rays (only diagnostic purpose)
            vecnew(jc,ispec) = jray

            ! Initielize height of the ray volume center-point
            zray(jc,jray) = zlaunch(il)
            ! NOTE: Check whether iexist below is consistent with the zray value
            !       here. If not, afterward routines will fail to track ray
            !       positions.  :  jkmid, jkmin, jkmax, ikhalf

            ! Ray size in z direction
            dzray(jc,jray) = dzlaunch_def

            ! Initialize horizontal wavenumbers and corresponding spectral extents
            ! depending on the direction of the launch momentum flux (so far 4
            ! directions)

            IF (iazidim == 4) THEN
              ! Eastward
              IF (7._wp*pi/4._wp < phi(iazi) .OR. phi(iazi) < 1._wp/4._wp*pi) THEN
                kray(jc,jray)  = kh(im,iomega)
                dkray(jc,jray) = dkh(im,iomega)
                lray(jc,jray)  = 0._wp
                dlray(jc,jray) = kh(im,iomega)*dphi
              ! Northward
              ELSEIF (pi/4._wp < phi(iazi) .AND. phi(iazi) < 3._wp/4._wp*pi) THEN
                kray(jc,jray)  = 0._wp
                dkray(jc,jray) = kh(im,iomega)*dphi
                lray(jc,jray)  = kh(im,iomega)
                dlray(jc,jray) = dkh(im,iomega)
              ! Westward
              ELSEIF (3._wp*pi/4._wp < phi(iazi) .AND. phi(iazi) < 5._wp/4._wp*pi) THEN
                kray(jc,jray)  = -kh(im,iomega)
                dkray(jc,jray) = dkh(im,iomega)
                lray(jc,jray)  = 0._wp
                dlray(jc,jray) = kh(im,iomega)*dphi
              ! Southward
              ELSEIF (5._wp*pi/4._wp < phi(iazi) .AND. phi(iazi) < 7._wp/4._wp*pi) THEN
                kray(jc,jray)  = 0._wp
                dkray(jc,jray) = kh(im,iomega)*dphi
                lray(jc,jray)  = -kh(im,iomega)
                dlray(jc,jray) = dkh(im,iomega)
              ENDIF
            ELSEIF (iazidim == 8) THEN
              ! Eastward
              IF (15._wp/8._wp*pi < phi(iazi) .OR. phi(iazi) < 1._wp/8._wp*pi) THEN
                kray(jc,jray)  = kh(im,iomega)
                dkray(jc,jray) = dkh(im,iomega)
                lray(jc,jray)  = 0._wp
                dlray(jc,jray) = kh(im,iomega)*dphi
              ! North-Eastward
              ELSEIF (3._wp/4._wp*pi < phi(iazi) .OR. phi(iazi) < 1._wp/8._wp*pi) THEN
                kray(jc,jray)  = kh(im,iomega)/SQRT(2._wp)
                dkray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
                lray(jc,jray)  = kh(im,iomega)/SQRT(2._wp)
                dlray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
              ! Northward
              ELSEIF (3._wp/8._wp*pi < phi(iazi) .AND. phi(iazi) < 5._wp/8._wp*pi) THEN
                kray(jc,jray)  = 0._wp
                dkray(jc,jray) = kh(im,iomega)*dphi
                lray(jc,jray)  = kh(im,iomega)
                dlray(jc,jray) = dkh(im,iomega)
              ! North-Westward
              ELSEIF (5._wp/8._wp*pi < phi(iazi) .AND. phi(iazi) < 7._wp/8._wp*pi) THEN
                kray(jc,jray)  = -kh(im,iomega)/SQRT(2._wp)
                dkray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
                lray(jc,jray)  = kh(im,iomega)/SQRT(2._wp)
                dlray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
              ! Westward
              ELSEIF (7._wp*pi/8._wp < phi(iazi) .AND. phi(iazi) < 9._wp/8._wp*pi) THEN
                kray(jc,jray)  = -kh(im,iomega)
                dkray(jc,jray) = dkh(im,iomega)
                lray(jc,jray)  = 0._wp
                dlray(jc,jray) = kh(im,iomega)*dphi
              ! South-Westward
              ELSEIF (9._wp*pi/8._wp < phi(iazi) .AND. phi(iazi) < 11._wp/8._wp*pi) THEN
                kray(jc,jray)  = -kh(im,iomega)/SQRT(2._wp)
                dkray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
                lray(jc,jray)  = -kh(im,iomega)/SQRT(2._wp)
                dlray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
              ! Southward
              ELSEIF (11._wp*pi/8._wp < phi(iazi) .AND. phi(iazi) < 13._wp/8._wp*pi) THEN
                kray(jc,jray)  = 0._wp
                dkray(jc,jray) = kh(im,iomega)*dphi
                lray(jc,jray)  = -kh(im,iomega)
                dlray(jc,jray) = dkh(im,iomega)
              ! South-Eastward
              ELSEIF (13._wp*pi/8._wp < phi(iazi) .AND. phi(iazi) < 15._wp/8._wp*pi) THEN
                kray(jc,jray)  = kh(im,iomega)/SQRT(2._wp)
                dkray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
                lray(jc,jray)  = -kh(im,iomega)/SQRT(2._wp)
                dlray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
              ENDIF
            ELSE
              CALL finish ('init_gw_orretal', 'Invalid iazidim!  Implemented iazidim = 4 or 8')
            ENDIF

            ! Horizontal propagation test
            IF (ltest_hprop .AND. iazidim == 4) THEN
              kray(jc,jray)  = 0._wp
              dkray(jc,jray) = 0._wp
              lray(jc,jray)  = 0._wp
              dlray(jc,jray) = 0._wp
              ! Eastward
              IF (7._wp*pi/4._wp < phi(iazi) .OR. phi(iazi) < 1._wp/4._wp*pi) THEN
                kray(jc,jray)  = kh(im,iomega)
                dkray(jc,jray) = dkh(im,iomega)
                lray(jc,jray)  = 0._wp
                dlray(jc,jray) = kh(im,iomega)*dphi
              ENDIF
              ! Southward
              !IF (5._wp*pi/4._wp < phi(iazi) .AND. phi(iazi) < 7._wp/4._wp*pi) THEN
              !  kray(jc,jray)  = 0._wp
              !  dkray(jc,jray) = kh(im,iomega)*dphi
              !  lray(jc,jray)  = -kh(im,iomega)
              !  dlray(jc,jray) = dkh(im,iomega)
              !ENDIF
              ! One may copy any azimuthal direction here 
              ! if testing other propagation directions
            ELSEIF (ltest_hprop .AND. iazidim == 8) THEN
              kray(jc,jray)  = 0._wp
              dkray(jc,jray) = 0._wp
              lray(jc,jray)  = 0._wp
              dlray(jc,jray) = 0._wp
              ! South-Eastward
              IF (13._wp*pi/8._wp < phi(iazi) .AND. phi(iazi) < 15._wp/8._wp*pi) THEN
                kray(jc,jray)  = kh(im,iomega)/SQRT(2._wp)
                dkray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
                lray(jc,jray)  = -kh(im,iomega)/SQRT(2._wp)
                dlray(jc,jray) = SQRT(dkh(im,iomega)*kh(im,iomega)*dphi)
              ENDIF
              ! One may copy any azimuthal direction here 
              ! if testing other propagation directions
            ENDIF

            ! Calculate vertical wavenumber at ray volume center-point
            mray(jc,jray) = m(im)

            ! Extent in m direction
            dmray(jc,jray) = dm(im)

            ! Vertical group velocity at the ray volume center-point 
            ! TODO: this is based on hydrostatic disp. rel. to be consistent
            !       with Orr et al. (2010) and Scinocca (2003). This might be
            !       changed at some point to the full disp. rel. but then one
            !       should think over Steps 1-3 carefully too (Jacobians).
            cgz = ABS(omega_l(iomega)/m(im))

            ! Initialize ray volumes with their corresponding flux fractions
            ! (This is what was supposedto be done in Step 4 ...)
            fluxray = flux3(ikh,im,iazi)/kh(im,iomega)

            ! Suggestion by Young-Ha: division by dphi for consistency 
            ! with Scinocca (2003), i.e. mo_gwd_wms.f90. 
            ! TODO: This may work for n_direction = 4 (not yet considered 
            !       for the other cases)
            fluxray = fluxray/dphi

            ! Initilize wave action density of ray volumes
            dens(jc,jray) = fluxray/kh(im,iomega)/cgz

            ! Set ray volume to an existing one
            iexist(jc,jray) = jkmid(il)      ! closest half level

            ! Diagnostic outputs
            IF (msg_level >= 15) THEN
              tmp3(:) = ABS(zhalf(jc,iexist(jc,jray)-1:iexist(jc,jray)+1) - zray(jc,jray))
              IF ( tmp3(2) > tmp3(1) .OR. tmp3(2) > tmp3(3) )  &
                &  CALL finish ('init_gw_orretal','re-define iexist')
              WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, kh, L_h, dkh, L_dkh',&
                jc, jray, ikh, iazi, im, kh(im,iomega), 2._wp*pi/kh(im,iomega), dkh(im,iomega), 2._wp*pi/dkh(im,iomega)
              CALL message('', TRIM(message_text))
              IF (kray(jc,jray) /= 0._wp) THEN
                WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, kray(jc,jray), L_x, dkray(jc,jray), L_dx',&
                         jc, jray, ikh, iazi, im, kray(jc,jray), 2._wp*pi/kray(jc,jray), dkray(jc,jray), 2._wp*pi/dkray(jc,jray)
                CALL message('', TRIM(message_text))
              ELSE
                WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, kray(jc,jray), L_x, dkray(jc,jray), L_dx',&
                         jc, jray, ikh, iazi, im, kray(jc,jray), 0._wp, dkray(jc,jray), 2._wp*pi/dkray(jc,jray)
                CALL message('', TRIM(message_text))
              ENDIF
              IF (lray(jc,jray) /= 0._wp) THEN
                WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, lray(jc,jray), L_y, dlray(jc,jray), L_dy',&
                         jc, jray, ikh, iazi, im, lray(jc,jray), 2._wp*pi/lray(jc,jray), dlray(jc,jray), 2._wp*pi/dlray(jc,jray)
                CALL message('', TRIM(message_text))
              ELSE
                WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, lray(jc,jray), L_y, dlray(jc,jray), L_dy',&
                         jc, jray, ikh, iazi, im, lray(jc,jray), 0._wp, dlray(jc,jray), 2._wp*pi/dlray(jc,jray)
                CALL message('', TRIM(message_text))
              ENDIF
              WRITE(message_text,'(a,5i6,5E12.4)') 'jc, jray, ikh, iazi, im, cin, mray(jc,jray), L_z, dmray(jc,jray), L_dz',&
                       jc, jray, ikh, iazi, im, cin(im), mray(jc,jray), 2._wp*pi/mray(jc,jray), dmray(jc,jray), 2._wp*pi/dmray(jc,jray)
              CALL message('', TRIM(message_text))
              WRITE(message_text,'(a,5i6,2E12.4)') 'jc, jray, ikh, iazi, im, zray(jc,jray), dzray(jc,jray)',&
                                                    jc, jray, ikh, iazi, im, zray(jc,jray), dzray(jc,jray)
              CALL message('', TRIM(message_text))
              WRITE(message_text,'(a,5i6,2E12.4)') 'jc, jray, ikh, iazi, im, cgz,         dens(jc,jray)',&
                                                    jc, jray, ikh, iazi, im, cgz, dens(jc,jray)
              CALL message('', TRIM(message_text))
              WRITE(message_text,'(a,5i6,2E12.4)') 'jc, jray, ikh, iazi, im, phi(iazi),         dphi',&
                                                    jc, jray, ikh, iazi, im, phi(iazi), dphi
              CALL message('', TRIM(message_text))
            ENDIF

            ! Diagnostics for DIAG6
            pu(jklaunch,iazi)=pu(jklaunch,iazi)+SQRT(kray(jc,jray)**2+lray(jc,jray)**2) &
                                *cgz*dens(jc,jray)*dkray(jc,jray)*dlray(jc,jray)*dmray(jc,jray)

          ENDDO ! il

          ! Retain ray volume index of this spectral 
          ! element for the next launch time
          jr_last(jc,ispec) = jray

        ENDDO ! im
      ENDDO ! ikh
    ENDDO ! iazi

    IF (msg_level >= 12) THEN
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG6 jklaunch, eastward  flux:', jklaunch, pu(jklaunch,1)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG6 jklaunch, northward flux:', jklaunch, pu(jklaunch,2)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG6 jklaunch, westward  flux:', jklaunch, pu(jklaunch,3)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6,E12.4)') 'DIAG6 jklaunch, southward flux:', jklaunch, pu(jklaunch,4)
      CALL message('', TRIM(message_text))
    ENDIF
    IF (msg_level >= 15) THEN
      WRITE(message_text,'(a,i6)') 'Max no. of rays allowed:', maxrays_bg(jg)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6)') 'No. of newly added rays:', COUNT(vecnew(jc,:)/=0)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a)') 'Init positions of new rays:'
      CALL message('', TRIM(message_text))
      IF (SUM(vecnew(jc,:))/=0) THEN
        DO ispec = 1, nrays_add_bg(jg)
          jray = vecnew(jc,ispec)
          IF (jray /= 0) THEN
            WRITE(message_text,'(a,i8,5E12.4)') 'jray, zray, dzray, mray, dmray, wadens:', &
            & jray, zray(jc,jray), dzray(jc,jray), mray(jc,jray), dmray(jc,jray), dens(jc,jray)
            CALL message('init_gw_orretal', TRIM(message_text))
          ENDIF
        ENDDO
      ENDIF
      WRITE(message_text,'(a,i6)') 'Total no. of rays:', nrays_bg(jg)
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i6)') 'Total no. of existing rays:', &
        &  COUNT(iexist(jc,jray_start:jray_end) /= 0)
      CALL message('', TRIM(message_text))
    ENDIF

  ENDDO ! jc

END SUBROUTINE init_gw_orretal
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE init_gw_conv(nlev,jb,i_startidx,i_endidx,jray_offset,nlaunch_max,  &
  &                     z,zhalf,clon,clat,cellarea,  &
  &                     flag_cgw)

  INTEGER , INTENT(IN)    :: nlev, jb
  INTEGER , INTENT(IN)    :: i_startidx, i_endidx
  INTEGER , INTENT(IN)    :: jray_offset
  INTEGER , INTENT(IN)    :: nlaunch_max
  REAL(wp), INTENT(IN)    :: z    (:,:)
  REAL(wp), INTENT(IN)    :: zhalf(:,:)
  REAL(wp), INTENT(IN)    :: clon(:)
  REAL(wp), INTENT(IN)    :: clat(:)
  REAL(wp), INTENT(IN)    :: cellarea(:)
  INTEGER , INTENT(INOUT) :: flag_cgw(:)

  INTEGER  :: nspec, nlaunch
  INTEGER  :: jc, jk, ispec, il
  INTEGER  :: jray, jrayprev
  INTEGER  :: specid
  INTEGER  :: jkmid(0:nlaunch_max)
  REAL(wp) :: zlaunch(nlaunch_max)
  REAL(wp) :: zb_ghost, zb_ghost_bef
  REAL(wp) :: zray_d_prev, zray_u_prev
  REAL(wp) :: dzlaunch
  REAL(wp) :: dzlaunch_max, half_dzlaunch
  REAL(wp) :: multi_dzl(nlaunch_max)
  REAL(wp) :: dlon, dlat
  REAL(wp) :: grid_coslat
  REAL(wp) :: tmp

  IF (msg_level >= 12)  CALL message('init_gw_conv',  &
    &  'MS-GWaM: initialize GW field from convection')

  !----------------------------------------------------------------------
  ! Purpose: 
  !         Launch ray volumes for the convective GW spectra calculated
  !         in the source scheme.
  !
  ! Method:
  !         -- Once a new convective GW spectrum has been calculated from
  !            a newly diagnosed convection, ray volumes are launched
  !            filling in the ghost layer. Since then, new volumes are
  !            launched continuously, filling in the layer, until the end
  !            of the cumulus-scheme interval.
  !         -- At the end of the cumulus-scheme interval (= at the
  !            beginning of the new calling of the scheme), previous ray
  !            volumes are tailored regardless of the existence of newly
  !            diagnosed convection.
  !         -- Controlling the injection duration and size of ray volumes
  !            integrated over the duration (via the tailoring) is
  !            required to launch an intended amount of the input flux
  !            as calculated from a source scheme.
  !----------------------------------------------------------------------

  nspec = nrays_add_cv(jg)

  dzlaunch      = gws_conv_config%dz_launch
  dzlaunch_max  = gws_conv_config%dz_launch*2.0_wp
  half_dzlaunch = gws_conv_config%dz_launch*0.5_wp
  DO il = 1, nlaunch_max
    multi_dzl(il) = gws_conv_config%dz_launch*REAL(il,KIND=wp)
  ENDDO

  ! Note: the launch mom. flux spectrum calculation is based on the assumption
  ! that the POSITIVE FREQUENCY BRANCH IS USED, which also means that a NEGATIVE
  ! VERTICAL WAVENUMBER spectrum will be used. This results in upward propagating 
  ! gravity waves (cgz > 0) injected at launch level.

  DO jc = i_startidx, i_endidx

    !---------------------------------------------------------------------------
    ! Tailoring the tails of previous ray volumes when necessary
    !---------------------------------------------------------------------------

#ifdef __test_awaypole
if(abs(clat(jc))>85.*deg2rad) cycle
#endif

    ! skip the generation if the latitude is ouside the bounds
    IF (abs(clat(jc)) > msgw_source_limit * deg2rad) CYCLE

    IF ( flag_cgw(jc) == 2 .OR. flag_cgw(jc) == -1 ) THEN
       ! 1) positive (negative): whether rays are to be (not to be) launched in
       !    this convection-call time interval.
       ! 2) flag_cgw == 2 or 20 (-2 or -20) if rays also had been (also had not
       !    been) launched during the previous convection-call time interval.
       ! 3) ABS(flag_cgw) >= 10 if the previous rays were tailored.
       ! 4) flag_cgw will be reset such that ABS(flag_cgw) == 1 or 2 once a new
       !    convection call occurs (even when convection is diagnosed to be not
       !    existing).

      DO jray = jray_offset_cv(jg)+1, jray_offset_cv(jg)+nrays_cv(jg)

        IF (      p_ray(jg)%iexist(jc,jray,jb) == 0  &
          &  .OR. p_ray(jg)%specid(jc,jray,jb) == 0 )  CYCLE

        ! Loop for ghost ray volumes (specid < 0) OR
        ! non-ghost ray volumes (i.e., centers being above the ghost layer)
        ! but whose tails are still in the ghost layer (specid > 0)

        zray_u_prev = p_ray(jg)%z(jc,jray,jb) + 0.5_wp*p_ray(jg)%dz(jc,jray,jb)

        IF (p_ray(jg)%specid(jc,jray,jb) < 0) THEN   ! ghost ray volumes
          IF (zray_u_prev <= zhalf(jc,p_ray(jg)%jk_active(jc,jray,jb))) THEN  ! remove
            p_ray(jg)%iexist(jc,jray,jb) = 0
            p_ray(jg)%specid(jc,jray,jb) = 0
            p_ray(jg)%wadens(jc,jray,jb) = 0._wp
            CYCLE
          END IF
        END IF

        ! Taylor the trailing part of volumes that are partially below the source level
        tmp = zray_u_prev - zhalf(jc,p_ray(jg)%jk_active(jc,jray,jb))
        IF (tmp < p_ray(jg)%dz(jc,jray,jb)) THEN
          p_ray(jg)%dz(jc,jray,jb) = tmp
          p_ray(jg)%z (jc,jray,jb) = zray_u_prev - 0.5_wp*tmp
          DO jk = p_ray(jg)%iexist(jc,jray,jb), 2, -1
            IF (z(jc,jk-1) > p_ray(jg)%z(jc,jray,jb)) THEN
              p_ray(jg)%iexist(jc,jray,jb) = jk  ;  EXIT
            END IF
          ENDDO
        END IF

        ! Make all existing (convective) ray volumes as 'normal' ones
        p_ray(jg)%specid(jc,jray,jb) = 0

      ENDDO

      ! Reset the record of jray
      p_ray(jg)%jr_last(jc,specid_offset_cv(jg)+1:specid_offset_cv(jg)+nspec,jb) = 0

      flag_cgw(jc) = flag_cgw(jc)*10

    END IF

    IF ( flag_cgw(jc) < 0 )  CYCLE

    !---------------------------------------------------------------------------
    ! Launch new rays
    !---------------------------------------------------------------------------

    dlat = SQRT(cellarea(jc))/grid_sphere_radius
    IF (is_plane_torus) THEN
      grid_coslat = 1._wp
      dlon = dlat
    ELSE
      grid_coslat = p_mgmgrid(jg)%coslat(jc,jb)
      dlon = dlat / p_mgmgrid(jg)%coslat(jc,jb)
    END IF

    jray = jray_offset

    zb_ghost = zhalf(jc,p_ray_conv(jg)%jk_source(jc,jb)) - depth_ghost_cv

    DO ispec = 1, nspec

      IF ( p_ray_conv(jg)%wadens(jc,ispec,jb) == 0.0_wp )  CYCLE

      specid = ispec + specid_offset_cv(jg)  ! make distinguishable from other sources

      IF (p_ray(jg)%jr_last(jc,specid,jb) > 0) THEN
         ! p_ray(jg)%specid(jc,jrayprev,jb) == -specid

        jrayprev = p_ray(jg)%jr_last(jc,specid,jb)

        zray_d_prev = p_ray(jg)%z(jc,jrayprev,jb) - 0.5_wp*p_ray(jg)%dz(jc,jrayprev,jb)

        ! B)
        IF (zray_d_prev < zb_ghost)  CYCLE  ! im (ispec)

        ! C)
        nlaunch = MIN(nlaunch_max,INT((zray_d_prev - zb_ghost)/dzlaunch)+1)
        zlaunch(1:nlaunch) = p_ray(jg)%z(jc,jrayprev,jb) - multi_dzl(1:nlaunch)
        jkmid(0) = p_ray(jg)%iexist(jc,jrayprev,jb)   ! lower limit of jkmid for searching

      ELSE   ! p_ray(jg)%jr_last(jc,specid,jb) <= 0)   ! A, D)

        nlaunch = nlaunch_max
        zlaunch(1:nlaunch) = zhalf(jc,p_ray_conv(jg)%jk_source(jc,jb))  &
          &                  + 0.5_wp*dzlaunch - multi_dzl(1:nlaunch)
        jkmid(0) = p_ray_conv(jg)%jk_source(jc,jb)   ! lower limit of jkmid for searching

      END IF

      DO il = 1, nlaunch
        jkmid(il) = nlev+1
        DO jk = jkmid(il-1), nlev   ! downward searching
          IF (z(jc,jk) <= zlaunch(il)) THEN
            jkmid(il) = jk  ;  EXIT
          END IF
        ENDDO
      ENDDO

      DO il = 1, nlaunch

        ! Ray counter
        jray = jray + 1
 
        ! Do not overwrite existing rays
        DO WHILE (p_ray(jg)%iexist(jc,jray,jb) /= 0)
          jray = jray + 1
        ENDDO

        p_ray(jg)%jk_active(jc,jray,jb) = p_ray_conv(jg)%jk_source(jc,jb)

        p_ray(jg)%specid(jc,jray,jb) = -specid
 
        p_ray(jg)%iexist(jc,jray,jb) = jkmid(il)
        p_ray(jg)%z     (jc,jray,jb) = zlaunch(il)
        p_ray(jg)%lon   (jc,jray,jb) = clon(jc)
        p_ray(jg)%lat   (jc,jray,jb) = clat(jc)
        p_ray(jg)%coslat(jc,jray,jb) = grid_coslat
        p_ray(jg)%dz    (jc,jray,jb) = dzlaunch
        p_ray(jg)%dlon  (jc,jray,jb) = dlon
        p_ray(jg)%dlat  (jc,jray,jb) = dlat

        p_ray(jg)%k (jc,jray,jb) = p_ray_conv(jg)%k (jc,ispec,jb)
        p_ray(jg)%l (jc,jray,jb) = p_ray_conv(jg)%l (jc,ispec,jb)
        p_ray(jg)%m (jc,jray,jb) = p_ray_conv(jg)%m (jc,ispec,jb)
        p_ray(jg)%dk(jc,jray,jb) = p_ray_conv(jg)%dk(jc,ispec,jb)
        p_ray(jg)%dl(jc,jray,jb) = p_ray_conv(jg)%dl(jc,ispec,jb)
        p_ray(jg)%dm(jc,jray,jb) = p_ray_conv(jg)%dm(jc,ispec,jb)

        p_ray(jg)%wadens(jc,jray,jb) = p_ray_conv(jg)%wadens(jc,ispec,jb)

      ENDDO  ! il

      ! Ray index is saved to track this spectral element at the next launch time
      p_ray(jg)%jr_last(jc,specid,jb) = jray

    ENDDO  ! ispec

  ENDDO  ! jc

END SUBROUTINE init_gw_conv
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE remove_rays(nlev,i_startidx,i_endidx,jray_start,jray_end,nrays_vacate,&
                       z,zhalf,dlonray,dlatray,dzray,coslatray,&
                       kray,dkray,lray,dlray,mray,dmray,dens,iexist,specid,&
                       bvf2,gammash2,fc2)
  INTEGER,       INTENT(IN)    :: nlev
  INTEGER,       INTENT(IN)    :: i_startidx ! first index of the block
  INTEGER,       INTENT(IN)    :: i_endidx   ! last index of the block
  INTEGER,       INTENT(IN)    :: jray_start
  INTEGER,       INTENT(IN)    :: jray_end
  INTEGER,       INTENT(IN)    :: nrays_vacate
  INTEGER,       INTENT(INOUT) :: iexist(:,:)
  INTEGER,       INTENT(INOUT) :: specid(:,:)
  REAL(wp),      INTENT(IN)    :: z(:,:), zhalf(:,:)
  REAL(wp),      INTENT(IN)    :: dlonray(:,:)
  REAL(wp),      INTENT(IN)    :: dlatray(:,:)
  REAL(wp),      INTENT(IN)    :: dzray(:,:)
  REAL(wp),      INTENT(IN)    :: coslatray(:,:)
  REAL(wp),      INTENT(IN)    :: kray(:,:), dkray(:,:)
  REAL(wp),      INTENT(IN)    :: lray(:,:), dlray(:,:)
  REAL(wp),      INTENT(IN)    :: mray(:,:), dmray(:,:)
  REAL(wp),      INTENT(INOUT) :: dens(:,:)
  REAL(wp),      INTENT(IN)    :: bvf2(:,:)
  REAL(wp),      INTENT(IN)    :: gammash2(:,:)
  REAL(wp),      INTENT(IN)    :: fc2(:)

  INTEGER                      :: jray, jc, jk, jkmid
  INTEGER                      :: jray_offset
  INTEGER                      :: nremove_total
  INTEGER                      :: nrays_active_total
  INTEGER                      :: nremove_layer(nlev+1)
  INTEGER                      :: nrays_active_layer(nlev+1)
  INTEGER                      :: jn
  LOGICAL                      :: l_mask(nlev+1)
  REAL(wp)                     :: rnp1(nlev+1)
  INTEGER                      :: minindex
  REAL(wp)                     :: quantity_ray(jray_start:jray_end)
  REAL(wp)                     :: Kh2, omega2
  REAL(wp)                     :: m2pG2
  LOGICAL                      :: l_warning(nproma)

  INTEGER                      :: min_nrays_z
  INTEGER                      :: nr_aim(nlev+1)
! INTEGER , PARAMETER          :: min_nrays_z_def = 10  ! consistent with subroutine 'remove_vol'
  INTEGER , PARAMETER          :: min_nrays_z_def = 0   ! consistent with the previous paper

  IF (msg_level >= 12) CALL message('remove_rays', 'MS-GWaM: remove rays')

  !----------------------------------------------------------------------
  ! Purpose: 
  !         Remove ray volumes if exceeding a certain amount in order to 
  !         limit CPU costs. Without this removal the only sink of ray 
  !         volumes is the propagation out of the vertical domain, which 
  !         might not compensate the accumulation by GW sources.
  !
  ! Method:
  !         -- remove ray volumes with too small wavelength (<~1m) and 
  !            those placed in B-V frequency regions that are statically 
  !            unstable
  !         -- calculate GW energy fraction carried by each ray volume
  !         -- remove ray volumes with lowest energy in a distributed 
  !            manner over all vertical layers
  !         -- only "old" rays (specid(jc,jray) >= 0) will be removed so 
  !            that no filtering of the GW source is done right away
  !         -- this subroutine is to be called separately for ray volumes 
  !            from each source (convective / background)
  ! TODO:
  !         Change this routine to a merge subroutine where ray volumes 
  !         close to each other in phase space are merged (an expensive 
  !         vesrion of such a scheme is available in the MS-GWaM toymodel
  !         but its application in the ICON implementation is out of 
  !         question).
  !----------------------------------------------------------------------

  jray_offset = jray_start - 1

  l_warning(:) = .FALSE.

  DO jc = i_startidx, i_endidx

    nrays_active_layer(:) = 0

    ! Loop either for background or for convectively emmitted ray volumes
    DO jray = jray_start, jray_end

      ! Ray volume with this spectral ID has just been launched and so we do not
      ! want to remove it
      IF ( specid(jc,jray) < 0 .OR. iexist(jc,jray) == 0 ) CYCLE

      ! Find the closest half level index to the midpoint of the ray
      jkmid = iexist(jc,jray)

      ! Effective horizontal wavenumber squared
      Kh2 = kray(jc,jray)**2 + lray(jc,jray)**2

      ! Vertical wavenumber squared (+ gamma^2 correction)
      m2pG2 = mray(jc,jray)**2 + gammash2(jc,jkmid)

      ! Intrinsic frequency squared with pinc scale height correction
      omega2 = (bvf2(jc,jkmid)*Kh2+fc2(jc)*m2pG2)/(Kh2 + m2pG2)

      ! Calculate energy fraction of ray volume
      quantity_ray(jray) = ABS(                                                          &
        &   (coslatray(jc,jray) * dlonray(jc,jray)) * dlatray(jc,jray) * dzray(jc,jray)  &
        &                       * dkray(jc,jray)    * dlray(jc,jray)   * dmray(jc,jray)  &
        &                       * dens(jc,jray) * SQRT(omega2) )

      ! Count active ray volumes layer by layer
      nrays_active_layer(jkmid) = nrays_active_layer(jkmid) + 1

    ENDDO ! jray

    ! Total no. of active ray volumes
    nrays_active_total = SUM(nrays_active_layer)

    ! No. of ray volumes to be removed
    nremove_total = nrays_vacate - COUNT(iexist(jc,jray_start:jray_end) == 0)

    ! If nothing is to be removed
    IF (nremove_total <= 0)  CYCLE  ! jc

    ! Normal case: the maximum no. of ray volumes to 
    ! be removed is nrays_add_bg(jg) or nrays_add_cv(jg)
    IF (nremove_total < nrays_active_total) THEN

      ! Minimum no. of ray volumes per layer
      ! (This is introduced in order to be consistent with the subroutine 'remove_vol'.
      ! The parameter setting of min_nrays_z_def = 0 corresponds to the previous code
      ! used in Boeloeni et al. and Kim et al.)
      min_nrays_z = MIN( (nrays_active_total - nremove_total)/(2*(nlev+1)), min_nrays_z_def )
        !                ^ half the layer-mean no. of ray volumes

      ! No. of ray volumes to be removed in each layer
      nremove_layer(:) = INT( (REAL(nremove_total)/REAL(nrays_active_total))  &
        &                     *REAL(nrays_active_layer(:)) )
      nremove_layer(:) = MAX(0, MIN(nrays_active_layer(:) - min_nrays_z, nremove_layer(:)))
      nr_aim(:) = nrays_active_layer(:) - nremove_layer(:)
      l_mask(:) = nr_aim(:) > min_nrays_z
      DO jray = 1, nremove_total - SUM(nremove_layer)
        rnp1(:) = REAL(nremove_layer(:)+1, sp)/REAL(MAX(1,nrays_active_layer(:)), sp)
        jk = MINLOC( rnp1, DIM=1, MASK=l_mask )
        nremove_layer(jk) = nremove_layer(jk) + 1
        nr_aim       (jk) = nr_aim       (jk) - 1
        l_mask       (jk) = nr_aim       (jk) > min_nrays_z
      ENDDO

      ! Remove ray volumes with the lowest wave energy
      DO jk = 1,nlev+1
        DO jn = 1,nremove_layer(jk)

          minindex = MINLOC( quantity_ray(:), DIM=1,                            &
            &                MASK=( iexist(jc,jray_start:jray_end) == jk .AND.  &
            &                       specid(jc,jray_start:jray_end) >= 0 ) )     &
            &        + jray_offset

          iexist(jc,minindex) = 0
          specid(jc,minindex) = 0
          dens  (jc,minindex) = 0._wp

        ENDDO ! jn
      ENDDO ! jk

    ELSE  ! nremove_total >= nrays_active_total

      ! It will be very strange if this happens. It may mean that almost 
      ! all of the array were occupied by volumes in the ghost layer. 
      ! A warning is kept for this case.
      ! TODO: if this ever happens, further analysis is needed  
      !       to find out the exact reason and how to avoid it

      l_warning(jc) = .TRUE.
      iexist(jc,jray_start:jray_end) = 0
      specid(jc,jray_start:jray_end) = 0
      dens  (jc,jray_start:jray_end) = 0._wp

    END IF

  ENDDO ! jc

  IF ( ANY(l_warning(:)) )  CALL message('',  &
    &  'remove_rays: nremove_total >= nrays_active_total')

END SUBROUTINE remove_rays
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE idx_rayedge(nlev,i_startidx,i_endidx,jray_start,jray_end,z,zhalf, &
                       iexist,specid,jk_active,zray,dzray, &
                       jk_full_rtop,jk_full_rbot,jk_half_rtop,jk_half_rbot)
  INTEGER,       INTENT(IN)    :: nlev
  INTEGER,       INTENT(IN)    :: i_startidx ! first index of the block
  INTEGER,       INTENT(IN)    :: i_endidx   ! last index of the block
  INTEGER,       INTENT(IN)    :: jray_start
  INTEGER,       INTENT(IN)    :: jray_end
  INTEGER,       INTENT(IN)    :: iexist(:,:)
  INTEGER,       INTENT(IN)    :: specid(:,:)
  INTEGER,       INTENT(IN)    :: jk_active(:,:)
  REAL(wp),      INTENT(IN)    :: z(:,:)
  REAL(wp),      INTENT(IN)    :: zhalf(:,:)
  REAL(wp),      INTENT(IN)    :: zray(:,:), dzray(:,:)
  INTEGER,       INTENT(INOUT) :: jk_full_rtop(:,:)
  INTEGER,       INTENT(INOUT) :: jk_full_rbot(:,:)
  INTEGER,       INTENT(INOUT) :: jk_half_rtop(:,:)
  INTEGER,       INTENT(INOUT) :: jk_half_rbot(:,:)
  ! Counters, indices
  INTEGER                      :: jk, jc
  INTEGER                      :: jkmid
  INTEGER                      :: jray
  REAL(wp)                     :: zray_d, zray_u
 
  !----------------------------------------------------------------------
  ! Purpose: 
  !         Calculate the vertical level indices corresponding to the 
  !         "bottom" and the "top" of a single ray volume.
  !
  ! Method:
  !         -- search for closest full (half) levels below (above)
  !            ray-volume top/bottom
  ! TODO:
  !----------------------------------------------------------------------
  jk_full_rtop(:,:) = 2
  jk_full_rbot(:,:) = nlev     ! nlevp1 - 1

  DO jray = jray_start, jray_end
    DO jc = i_startidx, i_endidx

      IF (iexist(jc,jray) == 0) CYCLE

      ! Ray-volume edges
      zray_d = zray(jc,jray)-dzray(jc,jray)*0.5_wp
      zray_u = zray(jc,jray)+dzray(jc,jray)*0.5_wp

      jkmid = iexist(jc,jray)

      ! Find the index of the highest full level below the top of the ray volume
      ! (== half level index closest to the top of the ray volume).
      ! in case z(jc,2) <= zray_u, jk_full_rtop = 2 as initialized
      DO jk = jkmid-1, 2, -1
        IF (z(jc,jk) > zray_u) THEN
!           jk_full_rtop(jc,jray) = MIN(nlev, jk+1)
            jk_full_rtop(jc,jray) = jk+1   ! possible to become nlevp1 but okay
            EXIT
        ENDIF
      ENDDO

      ! Find the index of the lowest half level above the top of the ray volume.
      IF (zhalf(jc,jk_full_rtop(jc,jray)) > zray_u) THEN
        jk_half_rtop(jc,jray) = jk_full_rtop(jc,jray)
      ELSE
        jk_half_rtop(jc,jray) = MAX(2, jk_full_rtop(jc,jray) - 1)
      END IF

      ! If it is a newly launched ray volume, position its jk_*_rbot to its jk_active
      IF (specid(jc,jray) < 0) THEN
        jk_full_rbot(jc,jray) = jk_active(jc,jray)
        jk_half_rbot(jc,jray) = jk_active(jc,jray)
        CYCLE
      END IF

      ! Find the index of the highest full level below the bottom of the ray volume
      ! (== half level index closest to the bottom of the ray volume).
      ! in case z(jc,nlev-1) >= zray_d, jk_full_rbot = nlev as initialized
      DO jk = jkmid, nlev-1
        IF (z(jc,jk) < zray_d) THEN
            jk_full_rbot(jc,jray) = jk
            EXIT
        ENDIF
      ENDDO

      ! Find the index of the lowest half level above the bottom of the ray volume.
      IF (zhalf(jc,jk_full_rbot(jc,jray)) > zray_d) THEN
        jk_half_rbot(jc,jray) = jk_full_rbot(jc,jray)     ! max: nlev
      ELSE
        jk_half_rbot(jc,jray) = jk_full_rbot(jc,jray) - 1
      END IF

    ENDDO
  ENDDO
 
END SUBROUTINE idx_rayedge
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE saturation(p_patch, jb, nlev,i_startidx,i_endidx,jray_start,jray_end,z,zhalf,jk_active,&
                      zray, dzray, coslatray, dlatray, dlonray, kray, dkray, lray, dlray, mray, dmray, dens, iexist, specid,&
                      jkmin,jkmax,bvf2,gammash2,fc2,rho,dt,&
                      kd,kd2,ld,ld2,md,md2,cgz_diag,A,A_save,B2,B2_save,mB2_new,mB2)
  TYPE(t_patch), TARGET,INTENT(IN)  :: p_patch
  INTEGER,       INTENT(IN)         :: jb, nlev
  INTEGER,       INTENT(IN)         :: i_startidx ! first index of the block
  INTEGER,       INTENT(IN)         :: i_endidx   ! last index of the block
  INTEGER,       INTENT(IN)         :: jray_start
  INTEGER,       INTENT(IN)         :: jray_end
  INTEGER,       INTENT(IN)         :: jk_active(:,:)
  INTEGER,       INTENT(IN)         :: jkmin(:,:)
  INTEGER,       INTENT(IN)         :: jkmax(:,:)
  INTEGER,       INTENT(INOUT)      :: iexist(:,:)
  INTEGER,       INTENT(INOUT)      :: specid(:,:)
  REAL(wp),      INTENT(IN)         :: z(:,:), zhalf(:,:)
  REAL(wp),      INTENT(IN)         :: zray(:,:), dzray(:,:)
  REAL(wp),      INTENT(IN)         :: coslatray(:, :), dlatray(:, :), dlonray(:, :)
  REAL(wp),      INTENT(IN)         :: mray(:,:), dmray(:,:)
  REAL(wp),      INTENT(IN)         :: kray(:,:), dkray(:,:)
  REAL(wp),      INTENT(IN)         :: lray(:,:), dlray(:,:)

  REAL(wp),      INTENT(IN)         :: bvf2(:,:), rho(:,:), gammash2(:,:)
  REAL(wp),      INTENT(IN)         :: fc2(:)          ! Coriolis parameter**2
  REAL(wp),      INTENT(IN)         :: dt              ! timestep in seconds
  REAL(wp),      INTENT(OUT)        :: mB2(:,:), mB2_new(:,:)
  REAL(wp),      INTENT(OUT)        :: kd(:,:), ld(:,:), md(:,:)
  REAL(wp),      INTENT(OUT)        :: kd2(:,:), ld2(:,:), md2(:,:), cgz_diag(:,:)
  REAL(wp),      INTENT(OUT)        :: B2(:,:), B2_save(:,:)
  REAL(wp),      INTENT(OUT)        :: A(:,:), A_save(:,:)
  REAL(wp),      INTENT(INOUT)      :: dens(:,:)
  ! Counters, indices
  INTEGER                    :: nlevp1
  INTEGER                    :: jk, jc
  INTEGER                    :: jk1, jk2
  INTEGER                    :: jray
  REAL(wp)                   :: zray_u
  REAL(wp)                   :: zray_d_eff
  ! Other subquantities
  REAL(wp)                   :: mB2K2(nproma,nlev+1)
  REAL(wp)                   :: diffusion_2dt(nproma,nlev+1)
  REAL(wp)                   :: counter(nproma,nlev+1)
  REAL(wp)                   :: counter_old(nproma,nlev+1)
  REAL(wp)                   :: tmp(nproma,nlev+1)
  REAL(wp)                   :: crit
  REAL(wp)                   :: m2, intg_dens, dkdldm
  REAL(wp)                   :: kappa_2dt, weight
  REAL(wp)                   :: dzi_o_dz, ray_o_cell_area
  REAL(wp)                   :: integral0
  REAL(wp)                   :: integral1
  REAL(wp)                   :: integral2
  REAL(wp)                   :: dens_old(nproma,jray_start:jray_end)
  REAL(wp)                   :: K2(nproma,jray_start:jray_end)
  REAL(wp)                   :: Kh2               ! Local wave vector squares
  REAL(wp)                   :: bvf2_ray          ! Local B-V freq**2
  REAL(wp)                   :: gammash2_ray      ! Local inverse pinc scale height
  REAL(wp)                   :: omega_ray         ! Local intrinsic freq
  REAL(wp)                   :: r_sq, r_sq_o_a

  IF (msg_level >= 12) CALL message('saturation', 'MS-GWaM: wave breaking parametrization')

  !----------------------------------------------------------------------
  ! Purpose: 
  !         Parameterize GW breaking (saturation), i.e. the corresponding 
  !         loss in GW energy and momentum fluxes.         
  !
  ! Method:
  !         The method is based on Lindzen's idea that GWs break when 
  !         they turn the potential temperature gradient to statically 
  !         unstable (d\theta/dz < 0).
  !         -- diagnose static instbility by \int dm m^2B^2 >= N^4
  !         -- calculate turbulent viscosity (in case of instability)
  !         -- calculate reduction in wave action of ray volumes based 
  !            on the turbulent viscosity (in case of instability)
  !         -- check whether static stability is set back after GW 
  !            saturation
  !         -- see Bölöni et al. (2016, 2021)
  !----------------------------------------------------------------------

  nlevp1 = nlev+1

  ! Initialize saturation amplitude (m^2*B^2) with zero
  mB2(:,:)   = 0._wp
  mB2K2(:,:) = 0._wp

  ! For diagnosic printouts
  IF (ldiagprof) THEN
    dens_old(:,:) = dens(:,jray_start:jray_end)
    WHERE (iexist(jray_start:jray_end,:) == 0)
      dens_old(:,:) = 0._wp
    END WHERE
  ENDIF

  IF (is_plane_torus) THEN
    r_sq = 1._wp
  ELSE
    r_sq = grid_sphere_radius**2
  END IF

  DO jc = i_startidx, i_endidx

    r_sq_o_a = r_sq/p_patch%cells%area(jc, jb)

    DO jray = jray_start, jray_end

      ! If ray does not exist nothing is done with it
      IF (iexist(jc,jray) == 0) CYCLE

      ! Define (effective) ray volume edges
      zray_u = zray(jc,jray)+0.5_wp*dzray(jc,jray)
      zray_d_eff = zray(jc,jray)-0.5_wp*dzray(jc,jray)
      IF (specid(jc,jray) < 0) THEN
        IF (zray_u < z(jc,jk_active(jc,jray)))  CYCLE
        zray_d_eff = MAX(z(jc,jk_active(jc,jray)),zray_d_eff)
      ENDIF

      ! Effective horizontal wavenumber squared
      Kh2 = kray(jc,jray)**2 + lray(jc,jray)**2
      ! Vertical wavenumber squared
      m2 = mray(jc,jray)**2
      ! Effective total wavenumber squared
      K2(jc,jray) = Kh2 + m2

      ! Calculate an average N^2 for the ray
      bvf2_ray = 0.5_wp*(bvf2(jc,jkmin(jc,jray)) + bvf2(jc,jkmax(jc,jray)))

      ! Calculate an average pinc scale height for the ray
      gammash2_ray = 0.5_wp*(gammash2(jc,jkmin(jc,jray)) + gammash2(jc,jkmax(jc,jray)))

      ! Intrinsic frequency at ray volume center-point with pinc scale height correction
      omega_ray = branch*SQRT((bvf2_ray*Kh2+fc2(jc)*(m2+gammash2_ray))  &
        &         /(K2(jc,jray)+gammash2_ray))

      ! Calculate integral2
      integral2 = dkray(jc,jray)*dlray(jc,jray)*dmray(jc,jray)*dens(jc,jray)  &
        &         *m2*Kh2/omega_ray

      ! Calculate integral1
      integral1 = integral2/K2(jc,jray)

      ! get the fraction of horizontal area that the ray occupies in the cell for accumulating its effect
      ray_o_cell_area = r_sq_o_a * dlatray(jc, jray) * (dlonray(jc, jray) * coslatray(jc, jray))

      ! Calculate m^2*B^2, K^2*m^2*B^2 with K^2=k^2+l^2+m^2
      DO jk = jkmin(jc,jray), jkmax(jc,jray)   ! jkmin > 1, jkmax < nlevp1

        ! dzi: size of ray in z-direction corresponding to the actual layer
        ! dz : layer depth
        dzi_o_dz = (MIN(zray_u,z(jc,jk-1)) - MAX(zray_d_eff,z(jc,jk)))  &
          &        /(z(jc,jk-1) - z(jc,jk))

        ! Calculate m^2*B^2
        mB2(jc,jk) = mB2(jc,jk) + dzi_o_dz * ray_o_cell_area * integral1

        ! Calculate K^2*m^2*B^2
        mB2K2(jc,jk) = mB2K2(jc,jk) + dzi_o_dz * ray_o_cell_area * integral2

      ENDDO ! jk

    ENDDO ! jray

  ENDDO ! jc


  ! m^2*B^2, K^2*m^2*B^2 multiplied by 2N^4/rho
  DO jk = 2,nlev
    DO jc = i_startidx, i_endidx
      tmp  (jc,jk) = 2._wp*bvf2(jc,jk)**2/rho(jc,jk)
      mB2  (jc,jk) = mB2  (jc,jk)*tmp(jc,jk)
      mB2K2(jc,jk) = mB2K2(jc,jk)*tmp(jc,jk)
    ENDDO
  ENDDO

  ! Calculating the diffusivity coefficient
  diffusion_2dt(:,:) = 0._wp
  DO jk = 2,nlev
    DO jc = i_startidx, i_endidx
      crit = mB2(jc,jk) - (alpha_sat*bvf2(jc,jk))**2
      IF (crit > 0._wp)  diffusion_2dt(jc,jk) = crit/mB2K2(jc,jk)
    ENDDO
  ENDDO
  DO jc = i_startidx, i_endidx
    diffusion_2dt(jc,1) = diffusion_2dt(jc,2)
    diffusion_2dt(jc,nlevp1) = diffusion_2dt(jc,nlev)
  ENDDO

  ! Loop for reducing wave action density 
  ! if m^2*B^2 exceeds saturation threshold
  DO jc = i_startidx, i_endidx

    DO jray = jray_start, jray_end

      ! Skip those ray volumes that were not involved in the above calculation, 
      ! i.e. they do not exist or they did not cross over yet the launch level 
      ! z_src
      IF ( iexist(jc,jray) == 0 ) CYCLE
      zray_u = zray(jc,jray)+0.5_wp*dzray(jc,jray)
      IF ( specid(jc,jray) < 0 .AND. &
           zray_u < z(jc,jk_active(jc,jray)) ) CYCLE

      ! Calculate largest possible diffusion coefficient among the layers 
      ! crossed over by the ray volume
      jk1 = MAX(1     ,jkmin(jc,jray)-2)
      jk2 = MIN(nlevp1,jkmax(jc,jray)+2)
      kappa_2dt = MAXVAL(diffusion_2dt(jc,jk1:jk2))

      ! Calculate weight for reducing wave action density (turbulent viscosity)
      weight = 1._wp - kappa_2dt*K2(jc,jray)

      ! Reduce spectral wave action density proportionally to weight
!      IF (weight > 0._wp) THEN
!        dens(jc,jray) = dens(jc,jray) * weight
!      ELSE ! weight <= 0._wp --> do not change sign of dens, rather
!           !                     suppress ray volume
!       iexist(jc,jray) = 0
!       specid(jc,jray) = 0
!       dens(jc,jray)   = 0._wp
!     ENDIF
      dens(jc,jray) = dens(jc,jray) * weight
      IF (dens(jc,jray) <= eps_wadens) THEN
        iexist(jc,jray) = 0
        specid(jc,jray) = 0
        dens  (jc,jray) = 0.
      ENDIF

    ENDDO ! jray

  ENDDO ! jc


  ! Re-calculate quantities for 
  ! diagnostic purposes
  IF (ldiagprof) THEN

    ! Initialize diag vars 
    counter(:,:) = 1._wp       ! why 1 ?
    counter_old(:,:) = 1._wp
    kd(:,:)      = 0._wp
    kd2(:,:)     = 0._wp
    ld(:,:)      = 0._wp
    ld2(:,:)     = 0._wp
    md(:,:)      = 0._wp
    md2(:,:)     = 0._wp
    cgz_diag(:,:)= 0._wp
    mB2_new(:,:) = 0._wp
    B2_save(:,:) = 0._wp
    B2(:,:)      = 0._wp
    A_save(:,:)  = 0._wp
    A(:,:)       = 0._wp

    DO jc = i_startidx, i_endidx

      r_sq_o_a = r_sq/p_patch%cells%area(jc, jb)

      DO jray = jray_start, jray_end

        ! Skip ray volume if it does not exist
        IF (dens_old(jc,jray) == 0._wp) CYCLE

        ! Define (effective) ray volume edges
        zray_u = zray(jc,jray)+0.5_wp*dzray(jc,jray)
        zray_d_eff = zray(jc,jray)-0.5_wp*dzray(jc,jray)
        IF (specid(jc,jray) < 0) THEN
          IF (zray_u < z(jc,jk_active(jc,jray)))  CYCLE
          zray_d_eff = MAX(z(jc,jk_active(jc,jray)),zray_d_eff)
        ENDIF

        ! Effective horizontal wavenumber squared
        Kh2 = kray(jc,jray)**2 + lray(jc,jray)**2
        ! Vertical wavenumber squared
        m2 = mray(jc,jray)**2

        ! Calculate an average N^2 for the ray
        bvf2_ray = 0.5_wp*(bvf2(jc,jkmin(jc,jray)) + bvf2(jc,jkmax(jc,jray)))

        ! Calculate an average pinc scale height for the ray
        gammash2_ray = 0.5_wp*(gammash2(jc,jkmin(jc,jray)) + gammash2(jc,jkmax(jc,jray)))

        ! Intrinsic frequency at ray position with pinc scale height correction
        omega_ray = branch*SQRT((bvf2_ray*Kh2+fc2(jc)*(m2+gammash2_ray))  &
          &         /(K2(jc,jray)+gammash2_ray))

        dkdldm = dkray(jc,jray)*dlray(jc,jray)*dmray(jc,jray)

        intg_dens = dkdldm*dens_old(jc,jray)

        ! Calculate integral0
        integral0 = intg_dens*Kh2/(omega_ray*K2(jc,jray))

        ! get the fraction of horizontal area that the ray occupies in the cell for accumulating its effect
        ray_o_cell_area = r_sq_o_a * dlatray(jc, jray) * (dlonray(jc, jray) * coslatray(jc, jray))

        DO jk = jkmin(jc,jray), jkmax(jc,jray)   ! jkmin > 1, jkmax < nlevp1

          ! dzi: size of ray in z-direction corresponding to the actual half layer
          ! dz : layer depth
          dzi_o_dz = (MIN(zray_u,z(jc,jk-1)) - MAX(zray_d_eff,z(jc,jk)))  &
            &        /(z(jc,jk-1) - z(jc,jk))

          ! Calculate mean B^2
          B2_save(jc,jk) = B2_save(jc,jk) + dzi_o_dz * ray_o_cell_area * integral0

          ! Wave action density (physical space): A = int N dkdldm
          A_save(jc,jk) = A_save(jc,jk) + dzi_o_dz * ray_o_cell_area * intg_dens

          ! Counter for computing mean
          counter_old(jc,jk) = counter_old(jc,jk) + dzi_o_dz * ray_o_cell_area

        ENDDO ! jk

        ! Skip ray volume if it does not exist
        IF (iexist(jc,jray) == 0) CYCLE  ! jray

        intg_dens = dkdldm*dens(jc,jray)

        ! Calculate integral2
        integral2 = intg_dens*m2*Kh2/omega_ray

        ! Calculate integral1
        integral1 = integral2/K2(jc,jray)

        ! Calculate integral0
        integral0 = integral1/m2

        DO jk = jkmin(jc,jray), jkmax(jc,jray)   ! jkmin > 1, jkmax < nlevp1

          ! dzi: size of ray in z-direction corresponding to the actual half layer
          ! dz : layer depth
          dzi_o_dz = (MIN(zray_u,z(jc,jk-1)) - MAX(zray_d_eff,z(jc,jk)))  &
            &        /(z(jc,jk-1) - z(jc,jk))

          ! Calculate m^2*B^2
          mB2_new(jc,jk) = mB2_new(jc,jk) + dzi_o_dz * ray_o_cell_area * integral1

          ! Calculate mean B^2
          B2(jc,jk) = B2(jc,jk) + dzi_o_dz * ray_o_cell_area * integral0

          ! Calculate mean k
          kd(jc,jk) = kd(jc,jk) + kray(jc,jray)

          ! Calculate mean k^2
          kd2(jc,jk) = kd2(jc,jk) + kray(jc,jray)**2

          ! Calculate mean l
          ld(jc,jk) = ld(jc,jk) + lray(jc,jray)

          ! Calculate mean l^2
          ld2(jc,jk) = ld2(jc,jk) + lray(jc,jray)**2

          ! Calculate mean m
          md(jc,jk) = md(jc,jk) + mray(jc,jray)

          ! Calculate mean m^2
          md2(jc,jk) = md2(jc,jk) + m2

          ! Calculate mean cgz with scale height correction
          cgz_diag(jc,jk) = cgz_diag(jc,jk) - mray(jc,jray)*(omega_ray**2 - &
                            fc2(jc))/(omega_ray*(K2(jc,jray)+gammash2_ray))

          ! Wave action density (physical space): A = int N dkdldm
          A(jc,jk) = A(jc,jk) + dzi_o_dz * ray_o_cell_area * intg_dens

          ! Counter for computing mean
          counter(jc,jk) = counter(jc,jk) + dzi_o_dz * ray_o_cell_area

        ENDDO ! jk

      ENDDO ! jray

    ENDDO ! jc

    ! m^2*B^2, B^2 multiplied by 2N^4/rho
    mB2_new(:,:) = mB2_new(:,:)*tmp(:,:)
    B2     (:,:) = B2     (:,:)*tmp(:,:)
    B2_save(:,:) = B2_save(:,:)*tmp(:,:)

    ! Calculate mean
    kd(:,:)       = kd(:,:)/counter(:,:)
    kd2(:,:)      = kd2(:,:)/counter(:,:)
    ld(:,:)       = ld(:,:)/counter(:,:)
    ld2(:,:)      = ld2(:,:)/counter(:,:)
    md(:,:)       = md(:,:)/counter(:,:)
    md2(:,:)      = md2(:,:)/counter(:,:)
    cgz_diag(:,:) = cgz_diag(:,:)/counter(:,:)
    B2(:,:)       = B2(:,:)/counter(:,:)
    A(:,:)        = A(:,:)/counter(:,:)
    B2_save(:,:)  = B2_save(:,:)/counter_old(:,:)
    A_save(:,:)   = A_save(:,:)/counter_old(:,:)

  ENDIF

END SUBROUTINE saturation
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE wave2grid(nlev, i_startidx, i_endidx, jray_start, jray_end,                        &
                     z, zhalf, cellarea, fc, fc2, lcalc_flux_4dir,                            &
                     theta, u, v,                                                             &
                     bvf2_full, bvf2_half, gammash2_full, gammash2_half,                      &
                     zray, dzray, coslatray, dlatray, dlonray,                                &
                     kray, dkray, lray, dlray, mray, dmray, dens, iexist, specid, jk_active,  &
                     jkmin_full, jkmax_full, jkmin_half, jkmax_half,                          &
                     uupflux, uvpflux, uwpflux, vvpflux, vwpflux,                             &
                     uuflux, uvflux, uwflux, vvflux, vwflux, utflux, vtflux)
  INTEGER,              INTENT(IN)  :: nlev
  INTEGER,              INTENT(IN)  :: i_startidx ! first index of the block
  INTEGER,              INTENT(IN)  :: i_endidx   ! last index of the block
  INTEGER,              INTENT(IN)  :: jray_start
  INTEGER,              INTENT(IN)  :: jray_end
  INTEGER,              INTENT(IN)  :: iexist(:,:)
  INTEGER,              INTENT(IN)  :: specid(:,:)
  INTEGER,              INTENT(IN)  :: jk_active(:,:)
  INTEGER,              INTENT(IN)  :: jkmin_full(:,:)
  INTEGER,              INTENT(IN)  :: jkmax_full(:,:)
  INTEGER,              INTENT(IN)  :: jkmin_half(:,:)
  INTEGER,              INTENT(IN)  :: jkmax_half(:,:)
  REAL(wp),             INTENT(IN)  :: z(:,:), zhalf(:,:)
  REAL(wp),             INTENT(IN)  :: fc(:), cellarea(:)   ! area of grid cell [m**2]
  REAL(wp),             INTENT(IN)  :: zray(:,:), dzray(:,:)
  REAL(wp),             INTENT(IN)  :: coslatray(:, :), dlatray(:, :), dlonray(:, :)
  REAL(wp),             INTENT(IN)  :: mray(:,:), dmray(:,:)
  REAL(wp),             INTENT(IN)  :: kray(:,:), dkray(:,:)
  REAL(wp),             INTENT(IN)  :: lray(:,:), dlray(:,:)
  REAL(wp),             INTENT(IN)  :: dens(:,:)
  REAL(wp),             INTENT(IN)  :: theta(:,:)
  REAL(wp),             INTENT(IN)  :: fc2(:)
  REAL(wp),             INTENT(IN)  :: bvf2_full(:,:)
  REAL(wp),             INTENT(IN)  :: bvf2_half(:,:)
  REAL(wp),             INTENT(IN)  :: gammash2_full(:,:)
  REAL(wp),             INTENT(IN)  :: gammash2_half(:,:)
  LOGICAL,              INTENT(IN)  :: lcalc_flux_4dir
  REAL(wp),             INTENT(IN)  :: u(:, :)
  REAL(wp),             INTENT(IN)  :: v(:, :)
  REAL(wp),             INTENT(OUT) :: uuflux(:,:)
  REAL(wp),             INTENT(OUT) :: uvflux(:,:)
  REAL(wp),             INTENT(OUT) :: uwflux(:,:)
  REAL(wp),             INTENT(OUT) :: vvflux(:,:)
  REAL(wp),             INTENT(OUT) :: vwflux(:,:)
  REAL(wp),             INTENT(OUT) :: utflux(:,:)
  REAL(wp),             INTENT(OUT) :: vtflux(:,:)
  REAL(wp),             INTENT(OUT) :: uupflux(:,:)
  REAL(wp),             INTENT(OUT) :: uvpflux(:,:)
  REAL(wp),             INTENT(OUT) :: uwpflux(:,:)
  REAL(wp),             INTENT(OUT) :: vvpflux(:,:)
  REAL(wp),             INTENT(OUT) :: vwpflux(:,:)

  INTEGER                      :: nlevp1
  INTEGER                      :: jk, jc
  INTEGER                      :: jray
  REAL(wp)                     :: zray_u
  REAL(wp)                     :: zray_d_eff
  REAL(wp)                     :: dzi_o_dz, r_sq_o_a, r_sq
  REAL(wp)                     :: K2, Kh2, Kh                       ! Local wave vector squares
  REAL(wp)                     :: omega, omega2                     ! Local intrinsic freqs
  REAL(wp)                     :: cgz                               ! Local vertical group velocity
  REAL(wp)                     :: cgz_freq_factor_om2               ! Local combined vertical group velocity
  REAL(wp)                     :: a_cgx_ray_jk, a_cgy_ray_jk        ! Local horizontal wave action fluxes (auxiliary)
  REAL(wp)                     :: a_cgz_ray_jk                      ! Local vertical wave action flux (auxiliary)
  REAL(wp)                     :: a_cgx_freq_fac_om2_ray_jk         ! Local combined wave action flux (auxiliary)
  REAL(wp)                     :: a_cgx_freq_fac_fc2_ray_jk         ! Local combined wave action flux (auxiliary)
  REAL(wp)                     :: a_cgy_freq_fac_om2_ray_jk         ! Local combined wave action flux (auxiliary)
  REAL(wp)                     :: a_cgy_freq_fac_fc2_ray_jk         ! Local combined wave action flux (auxiliary)
  REAL(wp)                     :: a_cgz_ray_freq_fac_om2_jk         ! Local combined wave action flux (auxiliary)
  REAL(wp)                     :: m2
  REAL(wp)                     :: K2_p_gam2
  REAL(wp)                     :: a_ray
  REAL(wp)                     :: a_ray_jk
  REAL(wp)                     :: m2_p_gam2
  REAL(wp)                     :: ptf_x_ray
  REAL(wp)                     :: ptf_y_ray
  REAL(wp)                     :: mf_xx_ray
  REAL(wp)                     :: mf_xy_ray   ! identical for pseudo-momentum and momentum fluxes
  REAL(wp)                     :: mf_xz_ray
  REAL(wp)                     :: mf_yy_ray
  REAL(wp)                     :: mf_yz_ray
  REAL(wp)                     :: pmf_xx_ray
  REAL(wp)                     :: pmf_xz_ray
  REAL(wp)                     :: pmf_yy_ray
  REAL(wp)                     :: pmf_yz_ray
  REAL(wp)                     :: fact_grid

  IF (msg_level >= 12) CALL message('wave2grid', 'MS-GWaM: flux and energy calculation')

  !----------------------------------------------------------------------
  ! Purpose: 
  !         Project GW momentum fluxes (and other quantities) carried by 
  !         Lagrangian ray volumes to the Eulerian grid.
  !
  ! Method:
  !         -- loop over ray volumes 
  !         -- calculate fractional flux (energy, etc.) corresponding to 
  !            each vertical layer they cross over (jkmin, jkmax known 
  !            from idx_rayedge)
  !         -- add up these fractional fluxes (energy, etc.) for each 
  !            vertical layer
  !         -- horizontal (vertical) fluxes are calculated on full (half) 
  !            levels
  !         -- see Bölöni et al. (2021) and Voelker et al. (2024)
  !
  ! TODO: 
  !----------------------------------------------------------------------

  nlevp1 = nlev+1

  ! Initialize outputs
  uuflux(:,:) = 0._wp  ; uvflux(:,:) = 0._wp  ; uwflux(:,:) = 0._wp
  vvflux(:,:) = 0._wp  ; vwflux(:,:) = 0._wp
  utflux(:,:) = 0._wp  ; vtflux(:,:) = 0._wp
  uupflux(:,:) = 0._wp ; uvpflux(:,:) = 0._wp ; uwpflux(:,:) = 0._wp
  vvpflux(:,:) = 0._wp ; vwpflux(:,:) = 0._wp

  ! Initialize temporary variables
  ptf_x_ray = 0._wp  ; ptf_y_ray = 0._wp
  mf_xx_ray = 0._wp  ; mf_xy_ray = 0._wp  ; mf_xz_ray = 0._wp
  mf_yy_ray = 0._wp  ; mf_yz_ray = 0._wp
  pmf_xx_ray = 0._wp ; pmf_xz_ray = 0._wp
  pmf_yy_ray = 0._wp ; pmf_yz_ray = 0._wp

  IF (is_plane_torus) THEN
    r_sq = 1._wp
  ELSE
    r_sq = grid_sphere_radius**2
  END IF

  DO jc = i_startidx, i_endidx

    r_sq_o_a = r_sq/cellarea(jc)

    DO jray = jray_start, jray_end

      IF (iexist(jc, jray) == 0) CYCLE

      ! Define (effective) ray volume edges
      zray_u = zray(jc, jray) + 0.5_wp * dzray(jc, jray)
      zray_d_eff = zray(jc, jray) - 0.5_wp * dzray(jc, jray)
      IF (specid(jc, jray) < 0) THEN
        IF (zray_u < z(jc, jk_active(jc, jray)))  CYCLE
        zray_d_eff = MAX(z(jc, jk_active(jc, jray)), zray_d_eff)
      ENDIF

      ! Effective horizontal wavenumber squared
      Kh2 = kray(jc, jray)**2 + lray(jc, jray)**2
!     ! Effective horizontal wavenumber
      Kh = SQRT(Kh2)
      ! Vertical wavenumber squared
      m2 = mray(jc, jray)**2
      ! Effective total wavenumber squared
      K2 = Kh2 + m2

      ! Wave action density (integrated over wavenumber space)
      ! get the fraction of horizontal area that the ray occupies in the cell for accumulating its effect
      a_ray = dkray(jc, jray) * dlray(jc, jray) * dmray(jc, jray) * dens(jc, jray)  &
        &   * dlatray(jc, jray) * (dlonray(jc, jray)*coslatray(jc, jray)) * r_sq_o_a
          !  ^ Here the formular of fractional area should be consistent with other parts
          !    (e.g. source formular, arealonk formular) for conservation and smoothness

      !
      ! Calculate horizontal fluxes on full levels
      !
      DO jk = jkmin_half(jc, jray), jkmax_half(jc, jray) ! jkmin_half > 1, jkmax_half < nlevp1
        ! Note: jkmin_half and jkmax_half denote the half-level indices above the
        ! top and bottom of the ray volume, respectively. An interval between
        ! two half levels, [jk, jk+1], is centered on the full level jk.

        ! dzi: size of ray in z-direction corresponding to the actual full layer
        ! dz : full layer depth
        dzi_o_dz = (MIN(zray_u, zhalf(jc, jk)) - MAX(zray_d_eff, zhalf(jc, jk+1))) / (zhalf(jc, jk) - zhalf(jc, jk+1))

        ! m^2+gamma^2
        m2_p_gam2 = m2 + gammash2_full(jc, jk)

        ! k^2+ml^2+gamma^2
        K2_p_gam2 = K2 + gammash2_full(jc, jk)  ! this can be changed to inv_K2_p_gam2

        ! Intrinsic frequency squared with pinc scale height correction
        omega2 = (bvf2_full(jc, jk) * Kh2 + fc2(jc) * m2_p_gam2) / K2_p_gam2
!       omega2 = MAX(fc2(jc),omega2)   ! where N^2 < f^2 : if it happens, this will
!                                      !   gives zero flux and only small energy

        ! Intrinsic frequency with pinc scale height correction
        omega = branch * SQRT(omega2)

        ! Wave action density (integrated over wavenumber space) relevant to the
        ! vertical level with index jk
        a_ray_jk = dzi_o_dz * a_ray

        ! Calculate potential-temperature flux (rho*u*theta
        ! The grid variables independent of ray variables are multiplied outside the jray loop
        ptf_x_ray = lray(jc,jray) * a_ray_jk * mray(jc, jray) / (omega * K2_p_gam2)
        utflux(jc,jk) = utflux(jc,jk) + ptf_x_ray
        ptf_y_ray = -kray(jc,jray) * a_ray_jk * mray(jc, jray) / (omega * K2_p_gam2)
        vtflux(jc,jk) = vtflux(jc,jk) + ptf_y_ray

        ! wave action fluxes
        ! the fluxes are calculated both individually and as combined
        ! quantities to avoid divergent behavior around omega -> f or N
        a_cgx_ray_jk = a_ray_jk * kray(jc, jray) * (bvf2_full(jc, jk) - omega2) / (omega * K2_p_gam2)
        a_cgy_ray_jk = a_ray_jk * lray(jc, jray) * (bvf2_full(jc, jk) - omega2) / (omega * K2_p_gam2)
        a_cgx_freq_fac_om2_ray_jk = a_ray_jk * kray(jc, jray) * omega * m2_p_gam2 / (Kh2 * K2_p_gam2)
        a_cgx_freq_fac_fc2_ray_jk = a_ray_jk * kray(jc, jray) * fc2(jc) / omega * m2_p_gam2 / (Kh2 * K2_p_gam2)
        a_cgy_freq_fac_om2_ray_jk = a_ray_jk * lray(jc, jray) * omega * m2_p_gam2 / (Kh2 * K2_p_gam2)
        a_cgy_freq_fac_fc2_ray_jk = a_ray_jk * lray(jc, jray) * fc2(jc) / omega * m2_p_gam2 / (Kh2 * K2_p_gam2)

        ! Calculate momentum flux (rho*u*u)
        mf_xx_ray  = a_cgx_freq_fac_om2_ray_jk * kray(jc, jray) + a_cgy_freq_fac_fc2_ray_jk * lray(jc, jray)
        uuflux(jc,jk)  = uuflux(jc,jk)  + mf_xx_ray
        pmf_xx_ray = a_cgx_ray_jk * kray(jc, jray)
        uupflux(jc,jk) = uupflux(jc,jk) + pmf_xx_ray

        ! Calculate momentum flux (rho*u*v)
        mf_xy_ray = a_cgx_ray_jk * lray(jc, jray)
        uvflux(jc,jk) = uvflux(jc,jk) + mf_xy_ray
        uvpflux(jc,jk) = uvpflux(jc,jk) + mf_xy_ray

        ! Calculate momentum flux (rho*v*v)
        mf_yy_ray = a_cgy_freq_fac_om2_ray_jk * lray(jc, jray) + a_cgx_freq_fac_fc2_ray_jk * kray(jc, jray)
        vvflux(jc,jk) = vvflux(jc,jk) + mf_yy_ray
        pmf_yy_ray = a_cgy_ray_jk * lray(jc, jray)
        vvpflux(jc,jk) = vvpflux(jc,jk) + pmf_yy_ray

      ENDDO ! jk

      !
      ! Calculate vertical fluxes on half levels
      !
      DO jk = jkmin_full(jc, jray), jkmax_full(jc, jray)   ! jkmin_full > 1, jkmax_full <= nlev
        ! Note: jkmin_full and jkmax_full denote the full-level indices below the
        ! top and bottom of the ray volume, respectively. An interval between
        ! two full levels, [jk-1, jk], is (approximately) centered on the half level jk.

        ! dzi: size of ray in z-direction corresponding to the actual half layer
        ! dz : layer depth
        dzi_o_dz = (MIN(zray_u, z(jc, jk-1)) - MAX(zray_d_eff, z(jc, jk))) / (z(jc, jk-1) - z(jc, jk))

        ! k^2 + l^2 + m^2 + gamma^2
        K2_p_gam2 = K2 + gammash2_half(jc, jk)  ! this can be changed to inv_K2_p_gam2

        ! Intrinsic frequency squared with pinc scale height correction
        omega2 = (bvf2_half(jc, jk) * Kh2 + fc2(jc) * (m2 + gammash2_half(jc, jk))) / K2_p_gam2
!       omega2 = MAX(fc2(jc),omega2)   ! where N^2 < f^2 : if it happens, this will
!                                      !   gives zero flux and only small energy
        omega = SQRT(omega2)

        ! Vertical group velocity with pinc scale height correction
        cgz = - branch * mray(jc, jray) * (omega2 - fc2(jc)) / (omega * K2_p_gam2)
        cgz_freq_factor_om2 = - branch * mray(jc, jray) * omega / K2_p_gam2

        ! Wave action times cg relevant to the vertical level with index jk
        a_cgz_ray_jk = dzi_o_dz * a_ray * cgz
        a_cgz_ray_freq_fac_om2_jk = dzi_o_dz * a_ray * cgz_freq_factor_om2

        ! Calculate momentum flux (rho*u*w)
        mf_xz_ray = kray(jc, jray) * a_cgz_ray_freq_fac_om2_jk
        uwflux(jc, jk) = uwflux(jc, jk) + mf_xz_ray
        pmf_xz_ray = kray(jc, jray) * a_cgz_ray_jk
        uwpflux(jc, jk) = uwpflux(jc, jk) + pmf_xz_ray

        ! Calculate momentum flux (rho*v*w)
        mf_yz_ray = lray(jc, jray) * a_cgz_ray_freq_fac_om2_jk
        vwflux(jc, jk) = vwflux(jc, jk) + mf_yz_ray
        pmf_yz_ray = lray(jc, jray) * a_cgz_ray_jk
        vwpflux(jc, jk) = vwpflux(jc, jk) + pmf_yz_ray

      ENDDO ! jk

    ENDDO ! jray

  ENDDO ! jc

  ! Finalize heat fluxes
  DO jk = 1, nlev
    DO jc = i_startidx, i_endidx
      fact_grid = fc(jc) * theta(jc, jk) * bvf2_full(jc, jk) / grav
      !
      utflux(jc, jk) = utflux(jc, jk) * fact_grid
      vtflux(jc, jk) = vtflux(jc, jk) * fact_grid
    ENDDO
  ENDDO

  ! Boundary conditions
  DO jc = i_startidx, i_endidx
    ! Boundaries to conserve momentum
    uwflux  (jc, 1) = 0._wp ; uwflux  (jc, nlevp1) = 0._wp
    vwflux  (jc, 1) = 0._wp ; vwflux  (jc, nlevp1) = 0._wp
    uwpflux (jc, 1) = 0._wp ; uwpflux (jc, nlevp1) = 0._wp
    vwpflux (jc, 1) = 0._wp ; vwpflux (jc, nlevp1) = 0._wp
  ENDDO ! jc

  ! Smooth pseudo-momentum fluxes / energy with a Shapiro filter
  IF (lsmooth) THEN
    CALL smooth_vert(nlev,  i_startidx,i_endidx,nsmooth,.FALSE., uuflux)
    CALL smooth_vert(nlev,  i_startidx,i_endidx,nsmooth,.FALSE., uvflux)
    CALL smooth_vert(nlevp1,i_startidx,i_endidx,nsmooth,.FALSE., uwflux)
    CALL smooth_vert(nlev,  i_startidx,i_endidx,nsmooth,.FALSE., vvflux)
    CALL smooth_vert(nlevp1,i_startidx,i_endidx,nsmooth,.FALSE., vwflux)
    CALL smooth_vert(nlev,  i_startidx,i_endidx,nsmooth,.FALSE., uupflux)
    CALL smooth_vert(nlev,  i_startidx,i_endidx,nsmooth,.FALSE., uvpflux)
    CALL smooth_vert(nlevp1,i_startidx,i_endidx,nsmooth,.FALSE., uwpflux)
    CALL smooth_vert(nlev,  i_startidx,i_endidx,nsmooth,.FALSE., vvpflux)
    CALL smooth_vert(nlevp1,i_startidx,i_endidx,nsmooth,.FALSE., vwpflux)
    CALL smooth_vert(nlev,  i_startidx,i_endidx,nsmooth,.FALSE., utflux)
    CALL smooth_vert(nlev,  i_startidx,i_endidx,nsmooth,.FALSE., vtflux)
  END IF

END SUBROUTINE wave2grid
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE tendency(p_patch, p_metrics, p_int_state, rho, temp, theta, p_fld)
  TYPE(t_patch),  TARGET,INTENT(IN)    :: p_patch
  TYPE(t_nh_metrics)    ,INTENT(IN)    :: p_metrics
  TYPE(t_int_state)     ,INTENT(IN)    :: p_int_state
  REAL(wp),      POINTER,INTENT(IN)    :: rho (:,:,:)   ! full lev
  REAL(wp),      POINTER,INTENT(IN)    :: temp (:,:,:)  ! full lev
  REAL(wp),      POINTER,INTENT(IN)    :: theta(:,:,:)  ! full lev
  TYPE(t_msgwam),        INTENT(INOUT) :: p_fld

  ! Local variables
  INTEGER               :: nlev
  INTEGER               :: rl_start, rl_end
  INTEGER               :: i_startblk, i_endblk
  INTEGER               :: i_startidx, i_endidx
  INTEGER               :: jk,jc,jce,jb,ist
  INTEGER               :: nblks_c,nblks_e
  INTEGER,      POINTER :: iidx(:,:,:), iblk(:,:,:)
  REAL(wp)              :: inv_rho, inv_dzrho
  REAL(wp)              :: fc_o_rho_th
  REAL(wp)              :: uufl_c1, uvfl_c1, vvfl_c1
  REAL(wp)              :: uufl_c2, uvfl_c2, vvfl_c2
  REAL(wp)              :: uupfl_c1, uvpfl_c1, vvpfl_c1
  REAL(wp)              :: uupfl_c2, uvpfl_c2, vvpfl_c2
  REAL(wp)              :: utfl_c1, vtfl_c1
  REAL(wp)              :: utfl_c2, vtfl_c2
  REAL(wp)              :: ddt_theta_gwd_mgm
  INTEGER               :: jk0
  REAL(wp)              :: inv_depth, flxgrad
  REAL(wp), ALLOCATABLE :: ufl_hor_e(:,:,:)  ! horizontal flux of u 
                                             ! at cell edge midpoints
  REAL(wp), ALLOCATABLE :: vfl_hor_e(:,:,:)  ! horizontal flux of v 
                                             ! at cell edge midpoints
  REAL(wp), ALLOCATABLE :: upfl_hor_e(:,:,:) ! horizontal pflux of u 
                                             ! at cell edge midpoints
  REAL(wp), ALLOCATABLE :: vpfl_hor_e(:,:,:) ! horizontal pflux of v 
                                             ! at cell edge midpoints
  REAL(wp), ALLOCATABLE :: ptfl_hor_e(:,:,:) ! horizontal flux of pot temp
                                             ! at cell edge midpoints 
  REAL(wp), ALLOCATABLE :: ufldiv_hor(:,:,:) ! horizontal div of u fluxes
  REAL(wp), ALLOCATABLE :: vfldiv_hor(:,:,:) ! horizontal div of v fluxes
  REAL(wp), ALLOCATABLE :: upfldiv_hor(:,:,:) ! horizontal div of u fluxes
  REAL(wp), ALLOCATABLE :: vpfldiv_hor(:,:,:) ! horizontal div of v fluxes
  REAL(wp), ALLOCATABLE :: ptfldiv_hor(:,:,:)! horizontal div of pot temp fluxes

  REAL(wp), PARAMETER :: gam = cpd/cvd

  IF (msg_level >= 12) CALL message('tendency', 'MS-GWaM: wind tendency calculation')

  !----------------------------------------------------------------------
  ! Purpose:
  !         Calculate tendencies based on momentum flux convergence and 
  !         the elastic terms
  !
  ! Method:
  !         tendencies = 1/rho*(vertical centered differences of pseudo-
  !                      momentum fluxes)
  !         staggering:
  !         -- pseudo-momentum fluxes defined on half levels
  !         -- tendencies defined on full levels
  !
  ! TODO:
  ! 
  !
  !----------------------------------------------------------------------

  nlev = p_patch%nlev

  IF ( msgw_lower_bound_opt == 1 ) THEN

    ! Exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    jk0 = nlev+1 - nlev_lbnd(jg) - 1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jk, i_startidx, i_endidx, inv_depth, flxgrad)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk,  &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        inv_depth = 1._wp/(p_metrics%z_mc(jc,jk0,jb) - p_metrics%z_mc(jc,nlev,jb))
        !
#ifndef __msgwam1d
        flxgrad = p_fld%uufl_mgm(jc,jk0,jb)*inv_depth
        DO jk = jk0+1, nlev-1
          p_fld%uufl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = p_fld%uvfl_mgm(jc,jk0,jb)*inv_depth
        DO jk = jk0+1, nlev-1
          p_fld%uvfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = p_fld%vvfl_mgm(jc,jk0,jb)*inv_depth
        DO jk = jk0+1, nlev-1
          p_fld%vvfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = p_fld%utfl_mgm(jc,jk0,jb)*inv_depth
        DO jk = jk0+1, nlev-1
          p_fld%utfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = p_fld%vtfl_mgm(jc,jk0,jb)*inv_depth
        DO jk = jk0+1, nlev-1
          p_fld%vtfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = p_fld%uupfl_mgm(jc,jk0,jb) * inv_depth
        DO jk = jk0+1, nlev-1
          p_fld%uupfl_mgm(jc,jk,jb) = flxgrad * (p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = p_fld%uvpfl_mgm(jc,jk0,jb) * inv_depth
        DO jk = jk0+1, nlev-1
          p_fld%uvpfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = p_fld%vvpfl_mgm(jc,jk0,jb) * inv_depth
        DO jk = jk0+1, nlev-1
          p_fld%vvpfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_mc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        
        p_fld%uupfl_mgm(jc,nlev,jb) = 0._wp
        p_fld%uvpfl_mgm(jc,nlev,jb) = 0._wp
        p_fld%vvpfl_mgm(jc,nlev,jb) = 0._wp
        p_fld%uufl_mgm(jc,nlev,jb) = 0._wp
        p_fld%uvfl_mgm(jc,nlev,jb) = 0._wp
        p_fld%vvfl_mgm(jc,nlev,jb) = 0._wp
        p_fld%utfl_mgm(jc,nlev,jb) = 0._wp
        p_fld%vtfl_mgm(jc,nlev,jb) = 0._wp
        !
#endif
        flxgrad = 0.5_wp*(p_fld%uwfl_mgm(jc,jk0,jb) + p_fld%uwfl_mgm(jc,jk0+1,jb))*inv_depth
        DO jk = jk0+1, nlev
          p_fld%uwfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = 0.5_wp*(p_fld%vwfl_mgm(jc,jk0,jb) + p_fld%vwfl_mgm(jc,jk0+1,jb))*inv_depth
        DO jk = jk0+1, nlev
          p_fld%vwfl_mgm(jc,jk,jb) = flxgrad*(p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = 0.5_wp * (p_fld%uwpfl_mgm(jc,jk0,jb) + p_fld%uwpfl_mgm(jc,jk0+1,jb)) * inv_depth
        DO jk = jk0+1, nlev
          p_fld%uwpfl_mgm(jc,jk,jb) = flxgrad * (p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        flxgrad = 0.5_wp * (p_fld%vwpfl_mgm(jc,jk0,jb) + p_fld%vwpfl_mgm(jc,jk0+1,jb)) * inv_depth
        DO jk = jk0+1, nlev
          p_fld%vwpfl_mgm(jc,jk,jb) = flxgrad * (p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_mc(jc,nlev,jb))
        ENDDO
        p_fld%uwfl_mgm(jc,nlev+1,jb) = 0._wp
        p_fld%vwfl_mgm(jc,nlev+1,jb) = 0._wp
        p_fld%uwpfl_mgm(jc,nlev+1,jb) = 0._wp
        p_fld%vwpfl_mgm(jc,nlev+1,jb) = 0._wp
      ENDDO  ! jc
    ENDDO  ! jb

!$OMP END DO
!$OMP END PARALLEL

  END IF  ! msgw_lower_bound_opt == 1

#ifndef __msgwam1d
  !define pointers
  iidx  => p_patch%edges%cell_idx
  iblk  => p_patch%edges%cell_blk

  ! Local dimensions
  nblks_c = p_patch%nblks_c
  nblks_e = p_patch%nblks_e

  ! Allocate auxiliary fields
  ALLOCATE( ufl_hor_e  (nproma,nlev,nblks_e), &
            vfl_hor_e  (nproma,nlev,nblks_e), &
            upfl_hor_e (nproma,nlev,nblks_e), &
            vpfl_hor_e (nproma,nlev,nblks_e), &
            ptfl_hor_e (nproma,nlev,nblks_e), &
            ufldiv_hor (nproma,nlev,nblks_c), &
            vfldiv_hor (nproma,nlev,nblks_c), &
            upfldiv_hor(nproma,nlev,nblks_c), &
            vpfldiv_hor(nproma,nlev,nblks_c), &
            ptfldiv_hor(nproma,nlev,nblks_c), &
            STAT=ist                        )
  IF (ist /= success) CALL finish ('tendency', 'Allocation of auxiliary fields failed!')

  ! Synchronize fields and fluxes
  CALL sync_patch_array_mult(SYNC_C1, p_patch, 5, p_fld%uufl_mgm, p_fld%uvfl_mgm,  &
    &                        p_fld%vvfl_mgm, p_fld%utfl_mgm, p_fld%vtfl_mgm)
  CALL sync_patch_array_mult(SYNC_C1, p_patch, 3, p_fld%uupfl_mgm, p_fld%uvpfl_mgm, p_fld%vvpfl_mgm)

  ! Initialize horizontal fluxes at cell edges
  ufl_hor_e = 0._wp ; vfl_hor_e = 0._wp ; ptfl_hor_e = 0._wp
  upfl_hor_e = 0._wp ; vpfl_hor_e = 0._wp

  ! Calculate normal horizontal fluxes at edge midpoints
  ! Exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_e + 1
  rl_end   = min_rledge_int

  i_startblk = p_patch%edges%start_block(rl_start)
  i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jce, i_startidx, i_endidx,                     &
!$OMP            uufl_c1, uvfl_c1, uufl_c2, uvfl_c2, vvfl_c1, vvfl_c2,  &
!$OMP            utfl_c1, vtfl_c1, utfl_c2, vtfl_c2)

  DO jb = i_startblk, i_endblk

    CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

#ifdef __LOOP_EXCHANGE
    DO jce = i_startidx, i_endidx
      DO jk = 1, nlev
#else
    DO jk = 1, nlev
      DO jce = i_startidx, i_endidx
#endif

        ! Horizontal flux of u at cell center points divided by rho
        uufl_c1 = p_fld%uufl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
        uvfl_c1 = p_fld%uvfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
        uufl_c2 = p_fld%uufl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
        uvfl_c2 = p_fld%uvfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
        ! Interpolation from cell center to cell edge midpoints
        ufl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                               * ( uufl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                               +   uvfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                             + p_int_state%c_lin_e(jce,2,jb) &
                               * ( uufl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                               +   uvfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

        ! Horizontal flux of v at cell center points divided by rho
        vvfl_c1 = p_fld%vvfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
        vvfl_c2 = p_fld%vvfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
        ! Interpolation from cell center to cell edge midpoints
        vfl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                               * ( uvfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                               +   vvfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                             + p_int_state%c_lin_e(jce,2,jb) &
                               * ( uvfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                               +   vvfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )
        
        ! Horizontal pflux of u at cell center points divided by rho
        uupfl_c1 = p_fld%uupfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
        uvpfl_c1 = p_fld%uvpfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
        uupfl_c2 = p_fld%uupfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
        uvpfl_c2 = p_fld%uvpfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
        ! Interpolation from cell center to cell edge midpoints
        upfl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                                * ( uupfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                                +   uvpfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                              + p_int_state%c_lin_e(jce,2,jb) &
                                * ( uupfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                                +   uvpfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

        ! Horizontal pflux of v at cell center points divided by rho
        vvpfl_c1 = p_fld%vvpfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
        vvpfl_c2 = p_fld%vvpfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
        ! Interpolation from cell center to cell edge midpoints
        vpfl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                                * ( uvpfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                                +   vvpfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                              + p_int_state%c_lin_e(jce,2,jb) &
                                * ( uvpfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                                +   vvpfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

        ! Horizontal flux of pot. temperature at cell center points divided by rho
        utfl_c1 = p_fld%utfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
        vtfl_c1 = p_fld%vtfl_mgm(iidx(jce,jb,1),jk,iblk(jce,jb,1))
        utfl_c2 = p_fld%utfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
        vtfl_c2 = p_fld%vtfl_mgm(iidx(jce,jb,2),jk,iblk(jce,jb,2))
        ! Interpolation from cell center to cell edge midpoints
        ptfl_hor_e(jce,jk,jb) = p_int_state%c_lin_e(jce,1,jb) &
                               * ( utfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v1  &
                               +   vtfl_c1 * p_patch%edges%primal_normal_cell(jce,jb,1)%v2 )&
                              + p_int_state%c_lin_e(jce,2,jb) &
                               * ( utfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v1  &
                               +   vtfl_c2 * p_patch%edges%primal_normal_cell(jce,jb,2)%v2 )

      ENDDO
    ENDDO

  ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL

  ! Synchronize horizontal fluxes at edge mid-points
  CALL sync_patch_array_mult(SYNC_E, p_patch, 3, ufl_hor_e, vfl_hor_e, ptfl_hor_e)
  CALL sync_patch_array_mult(SYNC_E, p_patch, 2, upfl_hor_e, vpfl_hor_e)

  ! Calculate horizontal divergence of u flux
  !CALL div(vec_e     = ufl_hor_e,   & !in
  !         ptr_patch = p_patch,     & !in
  !         ptr_int   = p_int_state, & !in
  !         div_vec_c = ufldiv_hor)    !out
  CALL div_avg(ufl_hor_e,             & !in
               p_patch,               & !in
               p_int_state,           & !in
               p_int_state%c_bln_avg, & !in
               ufldiv_hor)              !out

  ! Calculate horizontal divergence of v flux
  !CALL div(vec_e     = vfl_hor_e,   & !in
  !         ptr_patch = p_patch,     & !in
  !         ptr_int   = p_int_state, & !in
  !         div_vec_c = vfldiv_hor)    !out
  CALL div_avg(vfl_hor_e,             & !in
               p_patch,               & !in
               p_int_state,           & !in
               p_int_state%c_bln_avg, & !in
               vfldiv_hor)              !out

  ! Calculate horizontal divergence of u pflux
  CALL div_avg(upfl_hor_e,            & !in
               p_patch,               & !in
               p_int_state,           & !in
               p_int_state%c_bln_avg, & !in
               upfldiv_hor)             !out

  ! Calculate horizontal divergence of v pflux
  CALL div_avg(vpfl_hor_e,            & !in
               p_patch,               & !in
               p_int_state,           & !in
               p_int_state%c_bln_avg, & !in
               vpfldiv_hor)             !out

  ! Calculate horizontal divergence of pot temperature flux
  !CALL div(vec_e     = ptfl_hor_e,  & !in
  !         ptr_patch = p_patch,     & !in
  !         ptr_int   = p_int_state, & !in
  !         div_vec_c = ptfldiv_hor)   !out
  CALL div_avg(ptfl_hor_e,            & !in
               p_patch,               & !in
               p_int_state,           & !in
               p_int_state%c_bln_avg, & !in
               ptfldiv_hor)             !out
#endif

  ! Exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c + 1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx,  &
!$OMP            inv_rho, inv_dzrho, fc_o_rho_th, ddt_theta_gwd_mgm)

  DO jb = i_startblk, i_endblk
  
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        ! 1/rho, 1/(rho*dz)
        inv_rho   = 1._wp/rho(jc,jk,jb)
        inv_dzrho = 1._wp/(p_metrics%ddqz_z_full(jc,jk,jb)*rho(jc,jk,jb))

        ! f/(rho*theta)
        fc_o_rho_th = p_patch%cells%f_c(jc,jb)/(rho(jc,jk,jb)*theta(jc,jk,jb))

        ! Elastic terms (pot temperature fluxes are at full levels)
        p_fld%etx_mgm(jc,jk,jb) = - fc_o_rho_th * p_fld%vtfl_mgm(jc,jk,jb)
        p_fld%ety_mgm(jc,jk,jb) =   fc_o_rho_th * p_fld%utfl_mgm(jc,jk,jb)

        ! Calculate u tendency at cell centers on full levels
        p_fld%mfcxz_mgm(jc,jk,jb) = -(  p_fld%uwfl_mgm(jc,jk,jb)  &
                                      - p_fld%uwfl_mgm(jc,jk+1,jb) ) * inv_dzrho
#ifdef __msgwam1d
        p_fld%ddt_u_gwd_mgm(jc,jk,jb) = + p_fld%mfcxz_mgm(jc,jk,jb) & ! vert flux conv.
              + p_fld%etx_mgm(jc,jk,jb)                               ! elastic term
#else
        p_fld%ddt_u_gwd_mgm(jc,jk,jb) = + p_fld%mfcxz_mgm(jc,jk,jb) & ! vert flux conv.
              - ufldiv_hor(jc,jk,jb)*inv_rho                        & ! hori flux conv.
              + p_fld%etx_mgm(jc,jk,jb)                               ! elastic term
#endif

        ! Calculate v tendency at cell centers on full levels
        p_fld%mfcyz_mgm(jc,jk,jb) = -(  p_fld%vwfl_mgm(jc,jk,jb)  &
                                      - p_fld%vwfl_mgm(jc,jk+1,jb) ) * inv_dzrho

#ifdef __msgwam1d
        p_fld%ddt_v_gwd_mgm(jc,jk,jb) = p_fld%mfcyz_mgm(jc,jk,jb) & ! vert flux conv.
              + p_fld%ety_mgm(jc,jk,jb)                             ! elastic term
#else
        p_fld%ddt_v_gwd_mgm(jc,jk,jb) = p_fld%mfcyz_mgm(jc,jk,jb) & ! vert flux conv.
              - vfldiv_hor(jc,jk,jb)*inv_rho                      & ! hori flux conv.
              + p_fld%ety_mgm(jc,jk,jb)                             ! elastic term
#endif
        
        ! Calculate u tendency for pmom fluxes at cell centers on full levels
        p_fld%pmfcxz_mgm(jc,jk,jb) = -(  p_fld%uwpfl_mgm(jc,jk,jb)  &
                                       - p_fld%uwpfl_mgm(jc,jk+1,jb) ) * inv_dzrho

#ifdef __msgwam1d
        p_fld%ddt_u_gwd_pmom_mgm(jc,jk,jb) = p_fld%pmfcxz_mgm(jc,jk,jb)   ! vert flux conv.
#else
        p_fld%ddt_u_gwd_pmom_mgm(jc,jk,jb) = p_fld%pmfcxz_mgm(jc,jk,jb) & ! vert flux conv.
              - upfldiv_hor(jc,jk,jb)*inv_rho                             ! hori flux conv.
#endif

        ! Calculate v tendency at cell centers on full levels
        p_fld%pmfcyz_mgm(jc,jk,jb) = -(  p_fld%vwpfl_mgm(jc,jk,jb)  &
                                       - p_fld%vwpfl_mgm(jc,jk+1,jb) ) * inv_dzrho

#ifdef __msgwam1d
        p_fld%ddt_v_gwd_pmom_mgm(jc,jk,jb) = p_fld%pmfcyz_mgm(jc,jk,jb)   ! vert flux conv.
#else
        p_fld%ddt_v_gwd_pmom_mgm(jc,jk,jb) = p_fld%pmfcyz_mgm(jc,jk,jb) & ! vert flux conv.
              - vpfldiv_hor(jc,jk,jb)*inv_rho                             ! hori flux conv.
#endif

#ifndef __msgwam1d
        ! Calculate pot temp tendency at cell centers on full levels
        ddt_theta_gwd_mgm = - ptfldiv_hor(jc,jk,jb)*inv_rho   ! hori flux conv.
!!!!! TEST !!!!! TO BE CLEANED LATER
        p_fld%ddt_pt_gwd_mgm(jc,jk,jb) = ddt_theta_gwd_mgm
!!!!! TEST !!!!! TO BE CLEANED LATER
        ! Calculate temperature tendency from pot temp tendency
        p_fld%ddt_t_gwd_mgm(jc,jk,jb) = ddt_theta_gwd_mgm &
                    *temp(jc,jk,jb)/theta(jc,jk,jb)*gam  ! * pt_prog%exner(jc,jk,jb)*gam
#else
        p_fld%ddt_t_gwd_mgm(jc,jk,jb) = 0._wp
#endif

      ENDDO ! jc
    ENDDO ! jk

    ! Smoothing over 2*nsmooth+1 points
    IF (lsmootht) THEN
      ! Note: usually not necessary to smooth tendencies (defalult is .F.)
      IF (nlevsmootht >= nlev) nlevsmootht = nlev
      CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
        &              p_fld%ddt_u_gwd_mgm(:,1:nlevsmootht+1,jb))
      CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
        &              p_fld%ddt_v_gwd_mgm(:,1:nlevsmootht+1,jb))
      CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
        &              p_fld%ddt_u_gwd_pmom_mgm(:,1:nlevsmootht+1,jb))
      CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
        &              p_fld%ddt_v_gwd_pmom_mgm(:,1:nlevsmootht+1,jb))
      CALL smooth_vert(nlevsmootht,i_startidx,i_endidx,nsmooth,.FALSE.,  &
        &              p_fld%ddt_t_gwd_mgm(:,1:nlevsmootht+1,jb))
    ENDIF

  ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL

#ifndef __msgwam1d
  DEALLOCATE( ufl_hor_e, vfl_hor_e, ptfl_hor_e,    &
              ufldiv_hor, vfldiv_hor, ptfldiv_hor, &
              STAT=ist                             )
  IF (ist /= success) CALL finish ('tendency', 'Deallocation of auxiliary fields failed!')
  DEALLOCATE( upfl_hor_e, vpfl_hor_e,   &
              upfldiv_hor, vpfldiv_hor, &
              STAT=ist                  )
  IF (ist /= success) CALL finish ('tendency', 'Deallocation of auxiliary fields failed!')
#endif

  ! do heuristic checks of tendencies
  CALL heuristic_check_3d(p_patch, p_fld%ddt_u_gwd_mgm, 'ddt_u_gwd_mgm')
  CALL heuristic_check_3d(p_patch, p_fld%ddt_v_gwd_mgm, 'ddt_v_gwd_mgm')
  CALL heuristic_check_3d(p_patch, p_fld%ddt_u_gwd_pmom_mgm, 'ddt_u_gwd_pmom_mgm')
  CALL heuristic_check_3d(p_patch, p_fld%ddt_v_gwd_pmom_mgm, 'ddt_v_gwd_pmom_mgm')

END SUBROUTINE tendency
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE heuristic_check_3d(p_patch, array, array_name)

  TYPE(t_patch), TARGET, INTENT(IN) :: p_patch
  REAL(wp), INTENT(IN)  :: array(:, :, :)
  CHARACTER(len=*), INTENT(IN) :: array_name
  
  INTEGER :: jc, jk, jb             ! loop indices
  INTEGER :: i_startblk, i_endblk   ! block loop boundaries
  INTEGER :: i_startidx, i_endidx   ! cell loop boundaries
  INTEGER :: rl_start, rl_end       !
  INTEGER :: nlev                   ! number of vetical levels

  !----------------------------------------------------------------------
  ! Purpose:
  !         Check 3D array associated to p_patch for NaNs using
  !         IEEE arithmetics
  !
  ! Method:
  !         loop over all members and return message if NaN is found
  !----------------------------------------------------------------------

  nlev = p_patch%nlev
  
  rl_start = grf_bdywidth_c + 1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

  DO jb = i_startblk, i_endblk  ! block loop (jb)
  
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)

    DO jk = 1, nlev  ! vertical level loop (jk)
      DO jc = i_startidx, i_endidx  ! cell loop  (jc)

        IF (IEEE_IS_NAN(array(jc, jk, jb))) THEN
          WRITE(message_text, '(3a,3i5)') &
            & 'Found IEEE NaN in array ', TRIM(array_name), ' at indices (jc, jk, jb): ', jc, jk, jb
          CALL finish('heuristic_check_3d', TRIM(message_text))
        ENDIF ! IEEE NaN check
      ENDDO ! jc
    ENDDO ! jk
  ENDDO ! jb
  
END SUBROUTINE heuristic_check_3d
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE propagate_wave(dt_prop, jb, nlev, i_startidx, i_endidx, jray_start, jray_end, z, dz,   &
                          zhalf, dzhalf, lat, lon, fc2, fdfdlat, kvisc,                           &
                          bvf2, gammash2, u, v, dn2dz, dg2dz, dudz, dvdz,                         &
                          n2_hgrad, g2_hgrad, u_hgrad, v_hgrad,                                   &
                          dn2dz_hgrad, dg2dz_hgrad, dudz_hgrad, dvdz_hgrad,                       &
                          lonray, dlonray, latray, dlatray, zray, dzray, coslatray,               &
                          kray, dkray, lray, dlray, mray, dmray, dens,                            &
                          iexist, specid, jk_active, jr_last                                      )

  INTEGER,          INTENT(IN)    :: jb                           ! block index
  INTEGER,          INTENT(IN)    :: nlev                         ! number of vertical levels
  INTEGER,          INTENT(IN)    :: i_startidx                   ! first index of the block
  INTEGER,          INTENT(IN)    :: i_endidx                     ! last index of the block
  INTEGER,          INTENT(IN)    :: jray_start                   !
  INTEGER,          INTENT(IN)    :: jray_end                     !
  INTEGER,          INTENT(IN)    :: jk_active(:, :)              !
  INTEGER,          INTENT(INOUT) :: iexist(:, :)                 !
  INTEGER,          INTENT(INOUT) :: specid(:, :)                 !
  INTEGER,          INTENT(INOUT) :: jr_last(:, :)                !
  REAL(wp),         INTENT(IN)    :: dt_prop                      ! dt for propagation routine call
  REAL(wp),         INTENT(IN)    :: z(:, :)                      ! vertical grid
  REAL(wp),         INTENT(IN)    :: dz(:, :)                     ! vertical grid spacing
  REAL(wp),         INTENT(IN)    :: zhalf(:, :)                  ! vertical staggered grid (half levels)
  REAL(vp),         INTENT(IN)    :: dzhalf(:, :)                 ! vertical grid spacing of staggered grid (half levels)
  REAL(wp),         INTENT(IN)    :: lat(:)                       ! latitudes of cell centers
  REAL(wp),         INTENT(IN)    :: lon(:)                       ! longitudes of cell centers
  REAL(wp),         INTENT(IN)    :: fc2(:)                       ! Coriolis frequency squared
  REAL(wp),         INTENT(IN)    :: fdfdlat(:)                   ! gradient of fc2
  REAL(wp),         INTENT(IN)    :: kvisc(:, :)                  ! kinematic viscosity
  REAL(wp),         INTENT(IN)    :: bvf2(:, :)                   ! BVF squared
  REAL(wp),         INTENT(IN)    :: gammash2(:, :)               ! vertical gradient of pinc scale height
  REAL(wp),         INTENT(IN)    :: u(:, :)                      ! zonal veliocity (half levels)
  REAL(wp),         INTENT(IN)    :: v(:, :)                      ! meridional velocity (half levels)
  REAL(wp),         INTENT(IN)    :: dn2dz(:, :)                  ! vertical gradient of BVF2
  REAL(wp),         INTENT(IN)    :: dg2dz(:, :)                  ! vertical gradient of gammash2
  REAL(wp),         INTENT(IN)    :: dudz (:, :)                  ! vertical gradient of zonal velocity
  REAL(wp),         INTENT(IN)    :: dvdz (:, :)                  ! vertical gradient of meridional velocity
  REAL(wp),         INTENT(IN)    :: n2_hgrad(:, :, :)            ! horizontal gradients of BVF2
  REAL(wp),         INTENT(IN)    :: g2_hgrad(:, :, :)            ! horizontal gradients of gammash2
  REAL(wp),         INTENT(IN)    :: u_hgrad(:, :, :)             ! horizontal gradients of zonal velocity
  REAL(wp),         INTENT(IN)    :: v_hgrad(:, :, :)             ! horizontal gradients of meridional velocity
  REAL(wp),         INTENT(IN)    :: dn2dz_hgrad(:, :, :)         ! vert & hori gradients of BVF2
  REAL(wp),         INTENT(IN)    :: dg2dz_hgrad(:, :, :)         ! vert & hori gradients of gammash2
  REAL(wp),         INTENT(IN)    :: dudz_hgrad (:, :, :)         ! vert & hori gradients of zonal velocity
  REAL(wp),         INTENT(IN)    :: dvdz_hgrad (:, :, :)         ! vert & hori gradients of meridional velocity
  REAL(wp),         INTENT(INOUT) :: lonray(:, :), dlonray(:, :)  ! ray - zonal position and extent
  REAL(wp),         INTENT(INOUT) :: latray(:, :), dlatray(:, :)  ! ray - meridional position and extent
  REAL(wp),         INTENT(INOUT) :: zray(:, :),   dzray(:, :)    ! ray - vertical position and extent
  REAL(wp),         INTENT(INOUT) :: coslatray(:, :)              ! ray - cosine latitude
  REAL(wp),         INTENT(INOUT) :: kray(:, :),   dkray(:, :)    ! ray - zonal wave number vars
  REAL(wp),         INTENT(INOUT) :: lray(:, :),   dlray(:, :)    ! ray - meridional wave number vars
  REAL(wp),         INTENT(INOUT) :: mray(:, :),   dmray(:, :)    ! ray - vertical wave number vars
  REAL(wp),         INTENT(INOUT) :: dens(:, :)                   ! ray - wave action density
  ! Local array bounds:
  INTEGER                         :: jray                ! ray index
  INTEGER                         :: jk,jc               ! vertical and block indeces
  INTEGER                         :: ikbase              ! indices for interpolation
  INTEGER                         :: nlevp1, itsubstep
  INTEGER                         :: jstage
  INTEGER                         :: jk_lbnd             ! index of the lower boundary below
                                                         ! which the ray volume is considered to
                                                         ! be out of the domain and thus "killed"
  INTEGER                         :: ntsubstep
  ! Physical variabl s
  REAL(wp)                        :: dt_prop_frac                      ! dt of propagation routine call
  REAL(wp)                        :: dt_rk                             ! dt of the RK4 stages
  REAL(wp)                        :: arealonk                          ! ray volume area lon-k
  REAL(wp)                        :: arealatl                          ! ray volume area lat-l
  REAL(wp)                        :: areazm                            ! ray volume area z-m
  REAL(wp)                        :: lonray_st(0:4)                    ! lonray at RK stages
  REAL(wp)                        :: latray_st(0:4)                    ! latray at RK stages
  REAL(wp)                        :: zray_st(0:4)                      ! zray at RK stages
  REAL(wp)                        :: kray_st(0:4)                      ! kray at RK stages
  REAL(wp)                        :: lray_st(0:4)                      ! lray at RK stages
  REAL(wp)                        :: mray_st(0:4)                      ! mray at RK stages
  REAL(wp)                        :: dlonray_st(0:4)                   ! dlonray at RK stages
  REAL(wp)                        :: dlatray_st(0:4)                   ! dlatray at RK stages
  REAL(wp)                        :: dzray_st(0:4)                     ! dzray at RK stages
  REAL(wp)                        :: bvf2_l, bvf2_r                    ! B-V frequency**2 (zonal bounds)
  REAL(wp)                        :: bvf2_f, bvf2_b                    ! B-V frequency**2 (merdional bounds)
  REAL(wp)                        :: bvf2_d, bvf2_u                    ! B-V frequency**2 (vertical bounds)
  REAL(wp)                        :: bvf2_c                            ! B-V frequency**2 (ray locations)
  REAL(wp)                        :: gammash2_l, gammash2_r            ! inverse pinc scale height**2 (zonal bounds)
  REAL(wp)                        :: gammash2_f, gammash2_b            ! inverse pinc scale height**2 (meridional bounds)
  REAL(wp)                        :: gammash2_d, gammash2_u            ! inverse pinc scale height**2 (vertical bounds)
  REAL(wp)                        :: gammash2_c                        ! inverse pinc scale height**2 (ray locations)
  REAL(wp)                        :: dn2dz_c, dg2dz_c                  ! gradients of BVF**2, gamma**2 (ray locations)
  REAL(wp)                        :: dudz_c, dvdz_c                    ! wind gradients (ray locations)
  REAL(wp)                        :: dmdt_ray                          ! dm/dt
  REAL(wp)                        :: dkdt_ray                          ! dk/dt
  REAL(wp)                        :: dldt_ray                          ! dl/dt
  REAL(wp)                        :: K2kvisc0, K2kvisc1                ! molecular viscosity (old, new)
  REAL(wp)                        :: Kh2, m2                           ! Local wave vector squares
  REAL(wp)                        :: m2pG2_l, m2pG2_r                  ! gamma corrected vertical wave number squared (zonal bounds)
  REAL(wp)                        :: m2pG2_f, m2pG2_b                  ! gamma corrected vertical wave number squared (meridional bounds)
  REAL(wp)                        :: m2pG2_d, m2pG2_u                  ! gamma corrected vertical wave number squared (vertical bounds)
  REAL(wp)                        :: m2pG2_c                           ! gamma corrected vertical wave number squared (ray positions)
  REAL(wp)                        :: K2pG2_l, K2pG2_r                  ! gamma corrected total wave number squared (zonal bounds)
  REAL(wp)                        :: K2pG2_f, K2pG2_b                  ! gamma corrected total wave number squared (meridional bounds)
  REAL(wp)                        :: K2pG2_d, K2pG2_u                  ! gamma corrected total wave number squared (vertical bounds)
  REAL(wp)                        :: K2pG2_c                           ! gamma corrected total wave number squared (ray positions)
  REAL(wp)                        :: cglon_l, cglon_r                  ! Local zonal group velocities
  REAL(wp)                        :: cglat_f, cglat_b                  ! Local meridional group velocities
  REAL(wp)                        :: cgz_d, cgz_u, cgz_c               ! Local vertical group velocities
  REAL(wp)                        :: omega2_l, omega2_r                ! local frequency squared (zonal bounds)
  REAL(wp)                        :: omega2_f, omega2_b                ! local frequency squared (merdidional bounds)
  REAL(wp)                        :: omega2_d, omega2_u                ! local frequency squared (vertical bounds)
  REAL(wp)                        :: omega2_c                          ! local frequency squared (ray positions)
! REAL(wp)                        :: u_c, v_c                          ! zonal velocity at the ray volume center
  REAL(wp)                        :: u_l, u_r, u_f, u_b, u_u, u_d      ! zonal velocity at the ray volume faces
  REAL(wp)                        :: v_l, v_r, v_f, v_b, v_u, v_d      ! meridional velocity at the ray volume faces and the center
  REAL(wp)                        :: rdz_d, rdz_u, rdz_c               ! rel. vert. positions of the bottom, top, and center of the ray
  REAL(wp)                        :: dlon_l, dlon_r, dlon_c            ! relative longitude with respect to the cell center
  REAL(wp)                        :: dlat_f, dlat_b, dlat_c            ! relative latitude with respect to the cell center
  REAL(wp)                        :: bvf2_cc_d, n2_hgrad_d(1:2)        ! BVF2 and its gradient in the cell center
  REAL(wp)                        :: bvf2_cc_u, n2_hgrad_u(1:2)        ! BVF2 and its gradient in the	cell center
  REAL(wp)                        :: bvf2_cc_c, n2_hgrad_c(1:2)        ! BVF2 and its gradient in the	cell center
  REAL(wp)                        :: gammash2_cc_d, g2_hgrad_d(1:2)    ! Gamma2 and its gradient in the cell center
  REAL(wp)                        :: gammash2_cc_u, g2_hgrad_u(1:2)    ! Gamma2 and its gradient in the cell center
  REAL(wp)                        :: gammash2_cc_c, g2_hgrad_c(1:2)    ! Gamma2 and its gradient in the cell center
  REAL(wp)                        :: u_cc_d, u_hgrad_d(1:2)            ! horizontal wind and its gradient in the cell center
  REAL(wp)                        :: v_cc_d, v_hgrad_d(1:2)            ! horizontal wind and its gradient in the cell center
  REAL(wp)                        :: u_cc_u, u_hgrad_u(1:2)            ! horizontal wind and its gradient in the cell center
  REAL(wp)                        :: v_cc_u, v_hgrad_u(1:2)            ! horizontal wind and its gradient in the cell center
  REAL(wp)                        :: u_cc_c, u_hgrad_c(1:2)            ! horizontal wind and its gradient in the cell center
  REAL(wp)                        :: v_cc_c, v_hgrad_c(1:2)            ! horizontal wind and its gradient in the cell center
  REAL(wp)                        :: dn2dz_hgrad_c(1:2)                ! vert and hori gradients of BVF2 at the	cell center
  REAL(wp)                        :: dg2dz_hgrad_c(1:2)                ! vert and hori gradients of Gamma2 at the	cell center
  REAL(wp)                        :: dudz_hgrad_c(1:2)                 ! vert and hori gradients of U at the cell center
  REAL(wp)                        :: dvdz_hgrad_c(1:2)                 ! vert and hori gradients of V at the cell center
  REAL(wp)                        :: radray, rcoslat, metric           ! metric auxiliaries (radius, r cos(theta))
  REAL(wp)                        :: tanlatray
  REAL(wp)                        :: factor_c
  REAL(wp)                        :: acc_2K2kvisc
  REAL(wp)                        :: domxmin_ext, domxmax_ext
  REAL(wp)                        :: domymin_ext, domymax_ext
  REAL(wp)                        :: m_abs_max
  LOGICAL                         :: l_outofdomain
  LOGICAL                         :: l_refl_twice
  LOGICAL                         :: l_wavetooshort
  REAL(wp), PARAMETER             :: one_o_6 = 1._wp/6._wp
  INTEGER                         :: ik_st0
  REAL(wp)                        :: lonray_l, lonray_r   ! lon of ray (left and right)
  REAL(wp)                        :: latray_f, latray_b   ! lat of ray (front and back)
  REAL(wp)                        :: zray_d, zray_u       ! z of ray (down and up)
  INTEGER                         :: ik_d, ik_u

  ! TODO: changed definition of ray volume bounds, check consistency and modify through the routine

  IF (msg_level >= 12) CALL message('propagate_wave', 'MS-GWaM: propagate ray volumes')

  !----------------------------------------------------------------------
  ! Purpose:
  !         Propagate ray volumes (their center-points and edges) in z-m
  !         phase space.
  !
  ! Method:
  !         Solve the Ray equations with a RK4 scheme.
  !         -- see Bölöni et al. (2021) (2D)
  !         -- extension for from 2D to 6D
  !
  ! TODO:
  !         Complete this short documentation regarding 
  !         the extention to 6D
  !----------------------------------------------------------------------

  ! Number of vertical levels
  nlevp1 = nlev+1

  ! Index of the lower boundary below which the ray volume is
  ! considered to be out of the domain and thus "killed"
  jk_lbnd = nlevp1 - nlev_lbnd(jg)

  m_abs_max = 2.0_wp*pi/lambda_z_min

  ! Diagnostic printout to follow what happens in the first
  ! element of the first nproma block
  IF (msg_level >= 12) THEN
    DO jray = jray_start, jray_end
      WRITE(message_text,'(a,i6,5E12.4)') 'jray, zray, dzray, mray, dmray, wadens:', &
        &       jray, zray(1, jray), dzray(1, jray), mray(1, jray), dmray(1, jray), dens(1, jray)
      CALL message('propagate_wave', TRIM(message_text))
    ENDDO
    IF (msg_level >= 15) THEN
      DO jk = 1, nlev
        WRITE(message_text,'(a,i4,F14.5)') 'jk, full layer thicknesses:', &
          &       jk, dz(1,jk)
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,i4,F14.5)') 'jk, full level heights:', &
          &       jk, z(1,jk)
        CALL message('', TRIM(message_text))
      ENDDO
      DO jk = 1, nlevp1
        WRITE(message_text,'(a,i4,F14.5)') 'jk, half layer thicknesses:', &
          &       jk, dzhalf(1,jk)
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a,i4,F14.5)') 'jk, half level heights:', &
          &       jk, zhalf(1,jk)
        CALL message('', TRIM(message_text))
      ENDDO
      DO jk = 1, nlevp1
        WRITE(message_text,'(a,i4,F14.5)') 'jk, Brunt-Vaisala^2:', &
          &       jk, bvf2(1,jk)
        CALL message('', TRIM(message_text))
      ENDDO
    ENDIF
  ENDIF

  DO jc = i_startidx, i_endidx

    RAYLOOP: DO jray = jray_start, jray_end

#ifdef __test_awaypole
if(abs(latray(jc,jray))>85.*deg2rad) then
iexist(jc,jray) = 0 ; specid(jc,jray) = 0
endif
#endif

    ! Skip any calculations if ray does not really exist for some reason
    IF (iexist(jc,jray) == 0) CYCLE RAYLOOP

    ! Horizontal/vertical wavenumber squared
    Kh2 = kray(jc,jray)**2 + lray(jc,jray)**2
    m2  = mray(jc,jray)**2

    ! Introduce ray volume area in each slice of the 6D phase space
#ifndef __msgwam1d
    arealonk = (dlonray(jc,jray)*coslatray(jc,jray)) * dkray(jc,jray)
    arealatl = dlatray(jc,jray) * dlray(jc,jray)
#endif
    areazm   = dzray(jc,jray)   * dmray(jc,jray)


    ! The cases below distinguish whether the ray volume center-point is
    ! below/above the launch level (z_src). If below (specid(jc,jray) < 0),
    ! then no real time-stepping is done (see further comments), if above
    ! (specid(jc,jray) >= 0) then the complete ray equations are solved with
    ! an RK4 (nstages == 4) or with a forward Euler (nstages == 1) time step.

    IF (specid(jc,jray) < 0) THEN
      ! This is a newly launched ray volume below z_src. No real time-stepping
      ! is done therefore.

      ! Vertical wavenumber squared (+ gamma^2 correction) at ray volume
      ! center-point
      m2pG2_c = m2 + gammash2(jc, jk_active(jc, jray))

      ! k^2 + l^2 + m^2 + gamma^2 at ray volume center-point
      K2pG2_c = Kh2 + m2pG2_c

      ! Intrinsic frequency with pinc scale height correction at ray volume
      ! center-point
      omega2_c = (bvf2(jc, jk_active(jc, jray)) * Kh2 + fc2(jc) * m2pG2_c) / K2pG2_c

      ! Intrinsic group velocity at ray volume center-point
      cgz_c = ABS(mray(jc, jray) * ABS(omega2_c - fc2(jc)) / (SQRT(omega2_c) * K2pG2_c))
      cgz_c = MAX(1.E-10_wp, cgz_c)

      ! Diagnose whether the ray volume (center-point) would cross the
      ! launch level (z_src) within one time step dt (if yes dt_prop_frac > 0).
      ! Note: here dt_prop_frac is just a diagnostic time, it has nothing to do
      ! with the final time step used in the real propagation scheme
      ! for ray volumes with (specid(jc,jray) >= 0)
      dt_prop_frac = dt_prop - (zhalf(jc, jk_active(jc, jray)) - zray(jc, jray)) / cgz_c

      IF (dt_prop_frac < 0._wp) THEN
        ! If the ray volume does not cross the launch level within dt_prop
        ! then propagate it with dt_prop in a forward Euler step
        zray(jc,jray) = zray(jc,jray) + cgz_c * dt_prop
        ! Fill iexist with the layer index where the ray volume
        ! center-point is located
        DO jk = iexist(jc,jray),jk_active(jc,jray),-1
          IF (z(jc,jk-1) > zray(jc,jray)) THEN
            iexist(jc,jray) = jk  ;  EXIT
          ENDIF
        ENDDO
        CYCLE RAYLOOP
      ELSE
        ! If the ray volume does cross the launch level within dt_prop then
        ! re-initialize its position (center-point) at the launch level
        ! exactly
        zray  (jc,jray) = zhalf(jc,jk_active(jc,jray))
        ! Fill iexist with the layer index where the ray volume
        ! center-point is located
        iexist(jc,jray) = jk_active(jc,jray)
        ! Turn the ray volume to "old" by setting specid(jc,jray) > 0
        specid(jc,jray) = ABS(specid(jc,jray))
        ! Highest half-level index at which zhalf <= zray: launch-layer index
        ! = Highest full-level index at which zfull <= zray
        ikbase = jk_active(jc,jray)
        ! Tag in case the last-launched ray volume is escaping the ghost layer
        IF (jray == jr_last(jc,specid(jc,jray)))  jr_last(jc,specid(jc,jray)) = -1
      END IF

    ELSE ! (specid(jc,jray) >= 0)

      ! The ray volume is already "old", i.e. it has already crossed over
      ! z_src. Simply set the integration length as below and continue to
      ! perform a real time step
      dt_prop_frac = dt_prop

      ! Highest half-level index at which zhalf <= zray
      ikbase = iexist(jc,jray)   ! iexist: closest half level
      IF (zhalf(jc,ikbase) > zray(jc,jray))  ikbase = ikbase + 1
         ! ikbase must be >= 2. But by the above, it would be 1 if zray == zhalf(jc,1).
         !   This case has been excluded by l_outofdomain after integration at the previous call.
         !   (Otherwise, ikbase = MAX(2,ikbase) will be needed.)

    END IF

    ! RK integrations

    ! Interpolate kinematic viscosity to the ray volume center-point
    K2kvisc0 = ( kvisc(jc, ikbase) + (kvisc(jc, ikbase - 1) - kvisc(jc, ikbase))      &
      &              * ((zray(jc, jray) - zhalf(jc, ikbase)) / dz(jc, ikbase - 1)) )  &
      &       * (Kh2 + m2)

    ! Initialize cumulative kinematic viscosity over ntsubstep
    acc_2K2kvisc = 0._wp

    ! Determine the number of substeps
    IF (dt_gw_substep > 0._wp) THEN
      ! The default substep size has been prescribed: dt_gw_substep
      ! A round-off error is considered for, e.g., dt_prop/dt_gw_substep = 2.000000001
      ntsubstep = MAX(1, CEILING(dt_prop_frac / dt_gw_substep - 1.e-3))
    ELSE
      ! Apply the adaptive substep size depending on cgz

      ! Vertical wavenumber squared (+ gamma^2 correction) at ray volume
      ! center-point
      m2pG2_c = m2 + gammash2(jc, iexist(jc,jray))

      ! k^2 + l^2 + m^2 + gamma^2 at ray volume center-point
      K2pG2_c = Kh2 + m2pG2_c

      ! Intrinsic frequency with pinc scale height correction at ray volume
      ! center-point
      omega2_c = (bvf2(jc, iexist(jc,jray)) * Kh2 + fc2(jc) * m2pG2_c) / K2pG2_c

      ! Intrinsic group velocity at ray volume center-point
      cgz_c = ABS(mray(jc, jray) * ABS(omega2_c - fc2(jc)) / (SQRT(omega2_c) * K2pG2_c))

      ! Restrict the substep size such that the propagation distance at one substep
      ! does not exceed D*dz with D <= 1 (currently 0.8), if dz > 500 m.
      ! (The latter condition is considered since dz is very small near the surface ~ 10 m)
      ntsubstep = MAX(1, CEILING( dt_prop_frac * cgz_c / (0.8_wp * MAX(500._vp, dzhalf(jc, iexist(jc,jray)))) ))
    END IF

    ! Substep size for the integration
    dt_rk = dt_prop_frac/REAL(ntsubstep,KIND=wp)

    ! Iteration of RK4 integration substeps
    DO itsubstep = 1, ntsubstep

      ! Save ray volume related values from the previous time step ("0th" RK
      ! stage of the current time step)
      ik_st0 = ikbase
#ifndef __msgwam1d
      lonray_st(0)  = lonray(jc,jray)
      dlonray_st(0) = dlonray(jc,jray)
      kray_st(0)    = kray(jc,jray)
      latray_st(0)  = latray(jc,jray)
      dlatray_st(0) = dlatray(jc,jray)
      lray_st(0)    = lray(jc,jray)
#endif
      zray_st(0)    = zray(jc,jray)
      dzray_st(0)   = dzray(jc,jray)
      mray_st(0)    = mray(jc,jray)

      !                       About the RK4 implementation
      !
      ! The classical presentation of RK4 is that for dx/dt = F(x), we solve
      ! the following sub-integrals x_st(1), ... ,x_st(4) in the RK stages 1,2,3,4:
      !   x_st(1) = dt*F(x[n])
      !   x_st(2) = dt*F(x[n] + 1/2*x_st(1))
      !   x_st(3) = dt*F(x[n] + 1/2*x_st(2))
      !   x_st(4) = dt*F(x[n] + x_st(3))
      ! Then we calculate the time step based on these sub-integrals as:
      !   x[n+1] = x[n] + 1/6*(x_st(1) + 2*x_st(2) + 2*x_st(3) + x_st(4))       (1)
      !
      ! Here, instead the same result is obtained with the following scheme.
      ! We introduce xx to save different terms of eq.(1) above during the RK
      ! stage loop and calculate the following:
      ! stage 1:
      !   x_st(1) = dt*F(x[n])
      !        xx = x[n] + 1/2*x_st(1)
      ! stage 2:
      !   x_st(2) = dt*F(xx) with xx from stage 1
      !        xx = x[n] + 1/2*x_st(2)
      ! stage 3:
      !   x_st(3) = dt*F(xx) with xx from stage 2
      !        xx = x[n] + x_st(3)
      ! stage 4:
      !   x_st(4) = dt*F(xx) with xx from stage 3
      !        xx = x[n] + 1/6*x_st(4)
      !    x[n+1] = xx + 2*x_st(2) + 2*x_st(3) + 1/6*x_st(1)                    (2)
      ! The result of eq.(2) is equivalent to that of eq.(1).
      !
      !
      ! For the ray equations solved here:
      !    x_st(jstage): any of lonray_st(jstage), dlonray_st(jstage), kray(jstage)
      !                         latray_st(jstage), dlatray_st(jstage), lray(jstage)
      !                           zray_st(jstage),   dzray_st(jstage), mray(jstage)
      !    xx          : any of lonray(jc,jray), dlonray(jc,jray), kray(jc,jray)
      !                         latray(jc,jray), dlatray(jc,jray), lray(jc,jray)
      !                           zray(jc,jray),   dzray(jc,jray), mray(jc,jray)

      ! Loop over RK4 stages
      DO jstage = 1, nstages

        ! It is possible that after the 1st stage of RK4 integration, the edges of
        ! a ray volume have flipped over, i.e. negative dzray, dlonray or dlatray.
        ! If this happens even at the end of whole integration stages, we simply
        ! make it positive without changing its length.
        ! (Another option may be to make it to have an infinitesimal length, with
        ! positive sign.)

#ifndef __msgwam1d
        ! Ray volume west (east) edge if dlonray(jc,jray) > 0 (dlonray(jc,jray) < 0)
        lonray_l = lonray(jc, jray) - 0.5_wp * dlonray(jc, jray)
        ! Ray volume east (west) edge if dlonray(jc,jray) > 0 (dlonray(jc,jray) < 0)
        lonray_r = lonray(jc, jray) + 0.5_wp * dlonray(jc, jray)
        ! Ray volume south (north) edge if dlatray(jc,jray) > 0 (dlatray(jc,jray) < 0)
        latray_f = latray(jc, jray) - 0.5_wp * dlatray(jc, jray)
        ! Ray volume north (south) edge if dlatray(jc,jray) > 0 (dlatray(jc,jray) < 0)
        latray_b = latray(jc, jray) + 0.5_wp * dlatray(jc, jray)
#endif
        ! Ray volume bottom (top) if dzray(jc,jray) > 0 (dzray(jc,jray) < 0)
        zray_d = zray(jc, jray) - 0.5_wp * dzray(jc, jray)
        ! Ray volume top (bottom) if dzray(jc,jray) > 0 (dzray(jc,jray) < 0)
        zray_u = zray(jc, jray) + 0.5_wp * dzray(jc, jray)

        ! Find the half level in which the vertical edges of the ray volume are
        ! sitting. This is done separately for not reflected (dzray(jc,jray) >
        ! 0) and reflected (dzray(jc,jray) < 0) cases.
        IF (dzray(jc,jray) > 0._wp) THEN
          ! zray_d: bottom, zray_u: top
          ik_d = nlevp1  ;  ik_u = 2
          DO jk = ikbase,nlev
            IF (zhalf(jc,jk) <= zray_d) THEN
              ik_d = jk ! bottom found
              EXIT
            ENDIF
          ENDDO  ! in case zhalf(jc,nlev) > zray_d, ik_d = nlev+1 as initialized
          DO jk = ikbase-1,2,-1
            IF (zhalf(jc,jk) > zray_u) THEN
              ik_u = jk+1 ! top found
              EXIT
            ENDIF
          ENDDO  ! in case zhalf(jc,2) <= zray_u, ik_u = 2 as initialized
        ELSE
          ! zray_d: top, zray_u: bottom
          ik_u = nlevp1  ;  ik_d = 2
          DO jk = ikbase,nlev
            IF (zhalf(jc,jk) <= zray_u) THEN
              ik_u = jk ! bottom found
              EXIT
            ENDIF
          ENDDO  ! in case z(jc,nlev) > zray_d, ik_u = nlev+1 as initialized
          DO jk = ikbase-1,2,-1
            IF (zhalf(jc,jk) > zray_d) THEN
              ik_d = jk+1 ! top found
              EXIT
            ENDIF
          ENDDO  ! in case zhalf(jc,2) <= zray_u, ik_d = 2 as initialized
        ENDIF

        ! Interpolation coefs for ray volume bottom/top (rdz_d, rdz_u) and center point (rdz_c).
        rdz_d = (zray_d - zhalf(jc, ik_d)) / dz(jc, ik_d - 1)
        rdz_u = (zray_u - zhalf(jc, ik_u)) / dz(jc, ik_u - 1)
        rdz_c = (zray(jc, jray) - zhalf(jc, ikbase)) / dz(jc, ikbase - 1)

#ifndef __msgwam1d
        ! Interpolation coeffs for horizontal positions
        dlon_l = (lonray_l - lon(jc))*coslatray(jc, jray)
        dlon_r = (lonray_r - lon(jc))*coslatray(jc, jray)
        dlon_c = (lonray(jc, jray) - lon(jc))*coslatray(jc, jray)
        dlat_f = latray_f - lat(jc)
        dlat_b = latray_b - lat(jc)
        dlat_c = latray(jc, jray) - lat(jc)

        ! to prevent severe extrapolation
        dlon_l = MAX(-p_mgmgrid(jg)%darc_intpol(jc,jb), MIN(p_mgmgrid(jg)%darc_intpol(jc,jb), dlon_l))
        dlon_r = MAX(-p_mgmgrid(jg)%darc_intpol(jc,jb), MIN(p_mgmgrid(jg)%darc_intpol(jc,jb), dlon_r))
        dlon_c = MAX(-p_mgmgrid(jg)%darc_intpol(jc,jb), MIN(p_mgmgrid(jg)%darc_intpol(jc,jb), dlon_c))
        dlat_f = MAX(-p_mgmgrid(jg)%darc_intpol(jc,jb), MIN(p_mgmgrid(jg)%darc_intpol(jc,jb), dlat_f))
        dlat_b = MAX(-p_mgmgrid(jg)%darc_intpol(jc,jb), MIN(p_mgmgrid(jg)%darc_intpol(jc,jb), dlat_b))
        dlat_c = MAX(-p_mgmgrid(jg)%darc_intpol(jc,jb), MIN(p_mgmgrid(jg)%darc_intpol(jc,jb), dlat_c))
#endif

        ! Brunt-Väisälä frequency at the cell center
        bvf2_cc_d = bvf2(jc, ik_d) + (bvf2(jc, ik_d - 1) - bvf2(jc, ik_d)) * rdz_d
        bvf2_cc_u = bvf2(jc, ik_u) + (bvf2(jc, ik_u - 1) - bvf2(jc, ik_u)) * rdz_u
        bvf2_cc_c = bvf2(jc, ikbase) + (bvf2(jc, ikbase - 1) - bvf2(jc, ikbase)) * rdz_c

#ifndef __msgwam1d
!  ^
        ! Interpolate the horizontal gradients of N**2 in the vertical
        n2_hgrad_d = n2_hgrad(:, jc, ik_d) + (n2_hgrad(:, jc, ik_d - 1) - n2_hgrad(:, jc, ik_d)) * rdz_d
        n2_hgrad_u = n2_hgrad(:, jc, ik_u) + (n2_hgrad(:, jc, ik_u - 1) - n2_hgrad(:, jc, ik_u)) * rdz_u
        n2_hgrad_c = n2_hgrad(:, jc, ikbase) + (n2_hgrad(:, jc, ikbase - 1) - n2_hgrad(:, jc, ikbase)) * rdz_c


        ! bvf2 at the ray volume faces and center
        bvf2_c = MAX(bvf2_min, bvf2_cc_c + n2_hgrad_c(1) * dlon_c + n2_hgrad_c(2) * dlat_c)
        bvf2_l = MAX(bvf2_min, bvf2_cc_c + n2_hgrad_c(1) * dlon_l + n2_hgrad_c(2) * dlat_c)
        bvf2_r = MAX(bvf2_min, bvf2_cc_c + n2_hgrad_c(1) * dlon_r + n2_hgrad_c(2) * dlat_c)
        bvf2_f = MAX(bvf2_min, bvf2_cc_c + n2_hgrad_c(1) * dlon_c + n2_hgrad_c(2) * dlat_f)
        bvf2_b = MAX(bvf2_min, bvf2_cc_c + n2_hgrad_c(1) * dlon_c + n2_hgrad_c(2) * dlat_b)
        bvf2_d = MAX(bvf2_min, bvf2_cc_d + n2_hgrad_d(1) * dlon_c + n2_hgrad_d(2) * dlat_c)
        bvf2_u = MAX(bvf2_min, bvf2_cc_u + n2_hgrad_u(1) * dlon_c + n2_hgrad_u(2) * dlat_c)
#endif

        ! Inverse pinc scale height interpolated to ray volume edges (height) at the cell center-point
        gammash2_cc_d = gammash2(jc, ik_d) + (gammash2(jc, ik_d - 1) - gammash2(jc, ik_d)) * rdz_d
        gammash2_cc_u = gammash2(jc, ik_u) + (gammash2(jc, ik_u - 1) - gammash2(jc, ik_u)) * rdz_u
        gammash2_cc_c = gammash2(jc, ikbase) + (gammash2(jc, ikbase - 1) - gammash2(jc, ikbase)) * rdz_c

#ifndef __msgwam1d
        ! Interpolate the horizontal gradients of gamma in the vertical
        g2_hgrad_d = g2_hgrad(:, jc, ik_d) + (g2_hgrad(:, jc, ik_d - 1) - g2_hgrad(:, jc, ik_d)) * rdz_d
        g2_hgrad_u = g2_hgrad(:, jc, ik_u) + (g2_hgrad(:, jc, ik_u - 1) - g2_hgrad(:, jc, ik_u)) * rdz_u
        g2_hgrad_c = g2_hgrad(:, jc, ikbase) + (g2_hgrad(:, jc, ikbase - 1) - g2_hgrad(:, jc, ikbase)) * rdz_c

        ! gamma**2 at the ray volume faces and center
        gammash2_c = MAX(0._wp, gammash2_cc_c + g2_hgrad_c(1) * dlon_c + g2_hgrad_c(2) * dlat_c)
        gammash2_l = MAX(0._wp, gammash2_cc_c + g2_hgrad_c(1) * dlon_l + g2_hgrad_c(2) * dlat_c)
        gammash2_r = MAX(0._wp, gammash2_cc_c + g2_hgrad_c(1) * dlon_r + g2_hgrad_c(2) * dlat_c)
        gammash2_f = MAX(0._wp, gammash2_cc_c + g2_hgrad_c(1) * dlon_c + g2_hgrad_c(2) * dlat_f)
        gammash2_b = MAX(0._wp, gammash2_cc_c + g2_hgrad_c(1) * dlon_c + g2_hgrad_c(2) * dlat_b)
        gammash2_d = MAX(0._wp, gammash2_cc_d + g2_hgrad_d(1) * dlon_c + g2_hgrad_d(2) * dlat_c)
        gammash2_u = MAX(0._wp, gammash2_cc_u + g2_hgrad_u(1) * dlon_c + g2_hgrad_u(2) * dlat_c)
#endif

#ifndef __msgwam1d
        ! Vertical wavenumber squared (+ gamma^2 correction) at ray volume edges
        ! and center-point
        m2pG2_c = m2 + gammash2_c
        m2pG2_l = m2 + gammash2_l
        m2pG2_r = m2 + gammash2_r
        m2pG2_f = m2 + gammash2_f
        m2pG2_b = m2 + gammash2_b
        m2pG2_d = m2 + gammash2_d
        m2pG2_u = m2 + gammash2_u

        ! k^2 + l^2 + m^2 + gamma^2
        K2pG2_c = Kh2 + m2pG2_c
        K2pG2_l = Kh2 + m2pG2_l
        K2pG2_r = Kh2 + m2pG2_r
        K2pG2_f = Kh2 + m2pG2_f
        K2pG2_b = Kh2 + m2pG2_b
        K2pG2_d = Kh2 + m2pG2_d
        K2pG2_u = Kh2 + m2pG2_u

        ! Intrinsic frequency with pinc scale height correction at ray volume
        ! edges and center-point
        ! fc2 < omega2 < bvf2  since bvf2 >= bvf2_min > fc2
        omega2_c = (bvf2_c * Kh2 + fc2(jc) * m2pG2_c) / K2pG2_c
        omega2_l = (bvf2_l * Kh2 + fc2(jc) * m2pG2_l) / K2pG2_l
        omega2_r = (bvf2_r * Kh2 + fc2(jc) * m2pG2_r) / K2pG2_r
        omega2_f = (bvf2_f * Kh2 + fc2(jc) * m2pG2_f) / K2pG2_f
        omega2_b = (bvf2_b * Kh2 + fc2(jc) * m2pG2_b) / K2pG2_b
        omega2_d = (bvf2_d * Kh2 + fc2(jc) * m2pG2_d) / K2pG2_d
        omega2_u = (bvf2_u * Kh2 + fc2(jc) * m2pG2_u) / K2pG2_u
#else
        ! 1-D mode
        m2pG2_c = m2 + gammash2_cc_c
        K2pG2_c = Kh2 + m2pG2_c
        omega2_c = (bvf2_cc_c * Kh2 + fc2(jc) * m2pG2_c) / K2pG2_c
        m2pG2_d = m2 + gammash2_cc_d
        K2pG2_d = Kh2 + m2pG2_d
        omega2_d = (bvf2_cc_d * Kh2 + fc2(jc) * m2pG2_d) / K2pG2_d
        m2pG2_u = m2 + gammash2_cc_u
        K2pG2_u = Kh2 + m2pG2_u
        omega2_u = (bvf2_cc_u * Kh2 + fc2(jc) * m2pG2_u) / K2pG2_u
#endif

        ! u and v interpolated to ray volume edges (height) at the cell center-point
        u_cc_d = u(jc, ik_d) + (u(jc, ik_d - 1) - u(jc, ik_d)) * rdz_d
        u_cc_u = u(jc, ik_u) + (u(jc, ik_u - 1) - u(jc, ik_u)) * rdz_u
        u_cc_c = u(jc, ikbase) + (u(jc, ikbase - 1) - u(jc, ikbase)) * rdz_c
        v_cc_d = u(jc, ik_d) + (v(jc, ik_d - 1) - v(jc, ik_d)) * rdz_d
        v_cc_u = u(jc, ik_u) + (v(jc, ik_u - 1) - v(jc, ik_u)) * rdz_u
        v_cc_c = u(jc, ikbase) + (v(jc, ikbase - 1) - v(jc, ikbase)) * rdz_c

#ifndef __msgwam1d
        ! Interpolate the horizontal gradients of u and v in the vertical
        u_hgrad_d = u_hgrad(:, jc, ik_d) + (u_hgrad(:, jc, ik_d - 1) - u_hgrad(:, jc, ik_d)) * rdz_d
        u_hgrad_u = u_hgrad(:, jc, ik_u) + (u_hgrad(:, jc, ik_u - 1) - u_hgrad(:, jc, ik_u)) * rdz_u
        u_hgrad_c = u_hgrad(:, jc, ikbase) + (u_hgrad(:, jc, ikbase - 1) - u_hgrad(:, jc, ikbase)) * rdz_c
        v_hgrad_d = v_hgrad(:, jc, ik_d) + (v_hgrad(:, jc, ik_d - 1) - v_hgrad(:, jc, ik_d)) * rdz_d
        v_hgrad_u = v_hgrad(:, jc, ik_u) + (v_hgrad(:, jc, ik_u - 1) - v_hgrad(:, jc, ik_u)) * rdz_u
        v_hgrad_c = v_hgrad(:, jc, ikbase) + (v_hgrad(:, jc, ikbase - 1) - v_hgrad(:, jc, ikbase)) * rdz_c

        ! u at the ray volume faces and center
!       u_c = u_cc_c + u_hgrad_c(1) * dlon_c + u_hgrad_c(2) * dlat_c
        u_l = u_cc_c + u_hgrad_c(1) * dlon_l + u_hgrad_c(2) * dlat_c
        u_r = u_cc_c + u_hgrad_c(1) * dlon_r + u_hgrad_c(2) * dlat_c
        u_f = u_cc_c + u_hgrad_c(1) * dlon_c + u_hgrad_c(2) * dlat_f
        u_b = u_cc_c + u_hgrad_c(1) * dlon_c + u_hgrad_c(2) * dlat_b
        u_d = u_cc_d + u_hgrad_d(1) * dlon_c + u_hgrad_d(2) * dlat_c
        u_u = u_cc_u + u_hgrad_u(1) * dlon_c + u_hgrad_u(2) * dlat_c

        ! v at the ray volume faces and center
!       v_c = v_cc_c + v_hgrad_c(1) * dlon_c + v_hgrad_c(2) * dlat_c
        v_l = v_cc_c + v_hgrad_c(1) * dlon_l + v_hgrad_c(2) * dlat_c
        v_r = v_cc_c + v_hgrad_c(1) * dlon_r + v_hgrad_c(2) * dlat_c
        v_f = v_cc_c + v_hgrad_c(1) * dlon_c + v_hgrad_c(2) * dlat_f
        v_b = v_cc_c + v_hgrad_c(1) * dlon_c + v_hgrad_c(2) * dlat_b
        v_d = v_cc_d + v_hgrad_d(1) * dlon_c + v_hgrad_d(2) * dlat_c
        v_u = v_cc_u + v_hgrad_u(1) * dlon_c + v_hgrad_u(2) * dlat_c

        ! Metric terms related variables
        IF (.NOT. is_plane_torus) THEN
          tanlatray = SIN(latray(jc, jray)) / coslatray(jc, jray)
          ! Calculate r = a + z
          radray  = zray(jc, jray) + grid_sphere_radius
          ! Calculate r * cos(\phi)
          rcoslat = radray * coslatray(jc, jray)
          ! Switch for metric terms
          metric  = 1._wp
          IF (ltest_hprop .AND. (.NOT. ltest_gcircle)) THEN
            ! We need prescribed group velocities
            radray  = 1._wp
            rcoslat = 1._wp
            metric  = 0._wp
          ENDIF
        ELSE ! is_plane_torus
          tanlatray = 0._wp   ! not used (metric term)
          radray  = 1._wp
          rcoslat = 1._wp
          metric  = 0._wp
        ENDIF

        ! Longitudinal group velocity with pinc scale height correction at ray
        ! volume edges
        cglon_l = u_l + kray(jc, jray) * (bvf2_l - omega2_l) / (branch * SQRT(omega2_l) * K2pG2_l)
        cglon_r = u_r + kray(jc, jray) * (bvf2_r - omega2_r) / (branch * SQRT(omega2_r) * K2pG2_r)

        ! Meridional group velocity with pinc scale height correction at ray
        ! volume edges
        cglat_f = v_f + lray(jc, jray) * (bvf2_f - omega2_f) / (branch * SQRT(omega2_f) * K2pG2_f)
        cglat_b = v_b + lray(jc, jray) * (bvf2_b - omega2_b) / (branch * SQRT(omega2_b) * K2pG2_b)
#endif

        ! Vertical group velocity with pinc scale height correction at ray
        ! volume edges
        cgz_d = - mray(jc, jray) * (omega2_d - fc2(jc)) / (branch * SQRT(omega2_d) * K2pG2_d)
        cgz_u = - mray(jc, jray) * (omega2_u - fc2(jc)) / (branch * SQRT(omega2_u) * K2pG2_u)

        factor_c = 0.5_wp*branch/(SQRT(omega2_c) * K2pG2_c)

        ! cgz
        zray_st  (jstage) = 0.5_wp * (cgz_d + cgz_u)
#ifndef __msgwam1d
        ! cgh / r
        lonray_st(jstage) = 0.5_wp * (cglon_l + cglon_r) / radray
        latray_st(jstage) = 0.5_wp * (cglat_f + cglat_b) / radray

        ! Dk/Dt at the ray-volume center
#ifdef __dxy_edge
        dkdt_ray = -( (u_r - u_l) * kray(jc, jray) + (v_r - v_l) * lray(jc, jray) &
                     + factor_c * (Kh2 * (bvf2_r - bvf2_l) - (omega2_c - fc2(jc)) * (gammash2_r - gammash2_l)) ) &
                    / (rcoslat*dlonray(jc, jray)) &
#else
        dkdt_ray = -( u_hgrad_c(1) * kray(jc, jray) + v_hgrad_c(1) * lray(jc, jray) &
                     + factor_c * (Kh2 * n2_hgrad_c(1) - (omega2_c - fc2(jc)) * g2_hgrad_c(1)) ) &
                    / radray &
#endif
              + metric * kray(jc, jray) * ( latray_st(jstage) * tanlatray &
                                            - zray_st(jstage) / radray )

        ! Dl/Dt at the ray-volume center
#ifdef __dxy_edge
        dldt_ray = -( (u_b - u_f) * kray(jc, jray) + (v_b - v_f) * lray(jc, jray) &
                     + factor_c * (Kh2 * (bvf2_b - bvf2_f) - (omega2_c - fc2(jc)) * (gammash2_b - gammash2_f) &
                                    + m2pG2_c * fdfdlat(jc)*dlatray(jc, jray)) ) &
                    / (radray*dlatray(jc, jray)) &
#else
        dldt_ray = -( u_hgrad_c(2) * kray(jc, jray) + v_hgrad_c(2) * lray(jc, jray) &
                     + factor_c * (Kh2 * n2_hgrad_c(2) - (omega2_c - fc2(jc)) * g2_hgrad_c(2) &
                                    + m2pG2_c * fdfdlat(jc)) ) &
                    / radray &
#endif
              - metric * ( lonray_st(jstage) * kray(jc, jray) * tanlatray &
                           + zray_st(jstage) * lray(jc, jray) / radray )

        ! Dm/Dt at the ray-volume center
        !
        ! Vertical gradients of mean-flow variables interpolated
        ! vertically and horizontally to the ray-volume center
        dudz_hgrad_c(:) = dudz_hgrad(:, jc, ikbase)  &
          &            + (dudz_hgrad(:, jc, ikbase - 1) - dudz_hgrad(:, jc, ikbase)) * rdz_c
        dudz_c = dudz(jc, ikbase) + (dudz(jc, ikbase - 1) - dudz(jc, ikbase)) * rdz_c  &
          &      + dudz_hgrad_c(1) * dlon_c + dudz_hgrad_c(2) * dlat_c
        dvdz_hgrad_c(:) = dvdz_hgrad(:, jc, ikbase)  &
          &            + (dvdz_hgrad(:, jc, ikbase - 1) - dvdz_hgrad(:, jc, ikbase)) * rdz_c
        dvdz_c = dvdz(jc, ikbase) + (dvdz(jc, ikbase - 1) - dvdz(jc, ikbase)) * rdz_c  &
          &      + dvdz_hgrad_c(1) * dlon_c + dvdz_hgrad_c(2) * dlat_c
        dn2dz_hgrad_c(:) = dn2dz_hgrad(:, jc, ikbase)  &
          &             + (dn2dz_hgrad(:, jc, ikbase - 1) - dn2dz_hgrad(:, jc, ikbase)) * rdz_c
        dn2dz_c = dn2dz(jc, ikbase) + (dn2dz(jc, ikbase - 1) - dn2dz(jc, ikbase)) * rdz_c  &
          &       + dn2dz_hgrad_c(1) * dlon_c + dn2dz_hgrad_c(2) * dlat_c
        dg2dz_hgrad_c(:) = dg2dz_hgrad(:, jc, ikbase)  &
          &             + (dg2dz_hgrad(:, jc, ikbase - 1) - dg2dz_hgrad(:, jc, ikbase)) * rdz_c
        dg2dz_c = dg2dz(jc, ikbase) + (dg2dz(jc, ikbase - 1) - dg2dz(jc, ikbase)) * rdz_c  &
          &       + dg2dz_hgrad_c(1) * dlon_c + dg2dz_hgrad_c(2) * dlat_c
        !
        dmdt_ray = -( dudz_c * kray(jc, jray) + dvdz_c * lray(jc, jray) &
                     + factor_c * (Kh2 * dn2dz_c - (omega2_c - fc2(jc)) * dg2dz_c) ) &
            + metric * (  kray(jc, jray) * lonray_st(jstage) &
                        + lray(jc, jray) * latray_st(jstage) )
#else
  ! 1-D mode
        dudz_c = dudz(jc, ikbase) + (dudz(jc, ikbase - 1) - dudz(jc, ikbase)) * rdz_c
        dvdz_c = dvdz(jc, ikbase) + (dvdz(jc, ikbase - 1) - dvdz(jc, ikbase)) * rdz_c
        dn2dz_c = dn2dz(jc, ikbase) + (dn2dz(jc, ikbase - 1) - dn2dz(jc, ikbase)) * rdz_c
        dg2dz_c = dg2dz(jc, ikbase) + (dg2dz(jc, ikbase - 1) - dg2dz(jc, ikbase)) * rdz_c
        dmdt_ray = -( dudz_c * kray(jc, jray) + dvdz_c * lray(jc, jray) &
                     + factor_c * (Kh2 * dn2dz_c - (omega2_c - fc2(jc)) * dg2dz_c) )
#endif

#ifndef __msgwam1d
        ! Horizontal propagation tests with prescribed cg
        IF (ltest_hprop .AND. (.NOT. ltest_gcircle)) THEN
          IF (.NOT. is_plane_torus) THEN
            ! Prescribe constant longitudinal group velocity
            !cglon_l = signwn(kray(jc,jray))*30._wp/86400._wp*deg2rad ! 30 degrees/day (in radians/seconds) 
            !cglon_r = signwn(kray(jc,jray))*30._wp/86400._wp*deg2rad ! 30 degrees/day (in radians/seconds)
            cglon_l = signwn(kray(jc,jray))*3._wp/86400._wp*deg2rad ! 30 degrees/day (in radians/seconds) 
            cglon_r = signwn(kray(jc,jray))*3._wp/86400._wp*deg2rad ! 30 degrees/day (in radians/seconds)
            ! Prescribe constant meridional group velocity
            !cglat_f = signwn(lray(jc,jray))*30._wp/86400._wp*deg2rad ! 30 degrees/day (in radians/seconds) 
            !cglat_b = signwn(lray(jc,jray))*30._wp/86400._wp*deg2rad ! 30 degrees/day (in radians/seconds)
            cglat_f = signwn(lray(jc,jray))*3._wp/86400._wp*deg2rad ! 30 degrees/day (in radians/seconds) 
            cglat_b = signwn(lray(jc,jray))*3._wp/86400._wp*deg2rad ! 30 degrees/day (in radians/seconds)
            ! Prescribe zero vertical group velocity
            cgz_d = 0._wp
            cgz_u = 0._wp
            ! Prescribe constant wavenumbers
            dkdt_ray = 0._wp
            dldt_ray = 0._wp
            dmdt_ray = 0._wp
          ELSE ! is_plane_torus
            ! Prescribe constant group velocity in x dir
            cglon_l = signwn(kray(jc,jray))*200._wp
            cglon_r = signwn(kray(jc,jray))*200._wp
            ! Prescribe constant group velocity in y dir
            cglat_f = signwn(lray(jc,jray))*200._wp
            cglat_b = signwn(lray(jc,jray))*200._wp
            ! Prescribe zero vertical group velocity
            cgz_d = 0._wp
            cgz_u = 0._wp
            ! Prescribe constant wavenumbers
            dkdt_ray = 0._wp
            dldt_ray = 0._wp
            dmdt_ray = 0._wp
          ENDIF
          ! cgh / r
          lonray_st(jstage) = 0.5_wp * (cglon_l + cglon_r) / radray
          latray_st(jstage) = 0.5_wp * (cglat_f + cglat_b) / radray
          ! cgz
          zray_st  (jstage) = 0.5_wp * (cgz_d + cgz_u)
        ELSEIF (ltest_hprop .AND. ltest_gcircle) THEN
            ! Prescribe zero vertical group velocity
            cgz_d = 0._wp
            cgz_u = 0._wp
            zray_st(jstage) = 0._wp
        ENDIF
#endif

        ! Calculate x_st(jstage)
        ! jstage == 1     --> x_st(1)      = dt*F(x[n])
        ! jstage == 2,3,4 --> x_st(jstage) = dt*F(xx) with xx from stage jstage-1
        ! Displacement in phase-space (of ray volume center-point) and
        ! vertical deformation
#ifndef __msgwam1d
        lonray_st(jstage)  = dt_rk * lonray_st(jstage) / coslatray(jc, jray)
        dlonray_st(jstage) = dt_rk * (cglon_r - cglon_l) / rcoslat
        kray_st (jstage)   = dt_rk * dkdt_ray
        latray_st(jstage)  = dt_rk * latray_st(jstage)
        dlatray_st(jstage) = dt_rk * (cglat_b - cglat_f) / radray
        lray_st (jstage)   = dt_rk * dldt_ray
#endif
        zray_st (jstage)   = dt_rk * zray_st(jstage)
        dzray_st(jstage)   = dt_rk * (cgz_u - cgz_d)
        mray_st (jstage)   = dt_rk * dmdt_ray

        ! Calculate xx and save it for the next stage
        ! jstage == 1 --> xx = x[n] + (1/2)*x_st(1)
        ! jstage == 2 --> xx = x[n] + (1/2)*x_st(2)
        ! jstage == 3 --> xx = x[n] +       x_st(3)
        ! jstage == 4 --> xx = x[n] + (1/6)*x_st(4)
        ! Position (of ray volume center-point) and extent in phase-space
#ifndef __msgwam1d
        lonray(jc, jray)  = lonray_st (0) + coef_st(jstage) * lonray_st (jstage)
        dlonray(jc, jray) = dlonray_st(0) + coef_st(jstage) * dlonray_st(jstage)
        kray(jc, jray)    = kray_st   (0) + coef_st(jstage) * kray_st   (jstage)
        latray(jc, jray)  = latray_st (0) + coef_st(jstage) * latray_st (jstage)
        dlatray(jc, jray) = dlatray_st(0) + coef_st(jstage) * dlatray_st(jstage)
        lray(jc, jray)    = lray_st   (0) + coef_st(jstage) * lray_st   (jstage)
#endif
        zray(jc, jray)    = zray_st   (0) + coef_st(jstage) * zray_st   (jstage)
        dzray(jc, jray)   = dzray_st  (0) + coef_st(jstage) * dzray_st  (jstage)
        mray(jc, jray)    = mray_st   (0) + coef_st(jstage) * mray_st   (jstage)
        ! coef_st = (/ 1/2 , 1/2 ,  1  , 1/6 /)  for RK4
        !         = (/ 1   ,     ,     ,     /)  for Euler
        ! (defined in setup_msgwam of mo_setup_msgwam_interface.f90
        ! via rk4st_c)

        !
        ! Euler step done at this point with jstage == 1
        !

        IF (jstage == 4) THEN
          ! 4th RK stage, ready to calculate x[n+1]:
          !   x[n+1] = xx + (x_st(1) + 2*x_st(2) + 2*x_st(3))/6
          !   with  xx = x[n] + 1/6*x_st(4)
#ifndef __msgwam1d
          lonray(jc, jray) = lonray(jc, jray) &
                + (lonray_st(1) + 2._wp * (lonray_st(2) + lonray_st(3))) * one_o_6
          dlonray(jc, jray) = dlonray(jc, jray) &
                + (dlonray_st(1) + 2._wp * (dlonray_st(2) + dlonray_st(3))) * one_o_6
          kray(jc, jray) = kray(jc, jray) &
                + (kray_st(1) + 2._wp * (kray_st(2) + kray_st(3))) * one_o_6
          latray(jc, jray) = latray(jc, jray) &
                + (latray_st(1) + 2._wp * (latray_st(2) + latray_st(3))) * one_o_6
          dlatray(jc, jray) = dlatray(jc, jray) &
                + (dlatray_st(1) + 2._wp * (dlatray_st(2) + dlatray_st(3))) * one_o_6
          lray(jc, jray) = lray(jc, jray) &
                + (lray_st(1) + 2._wp * (lray_st(2) + lray_st(3))) * one_o_6
#endif
          zray(jc, jray) = zray(jc, jray) &
                + (zray_st(1) + 2._wp * (zray_st(2) + zray_st(3))) * one_o_6
          dzray(jc, jray) = dzray(jc, jray) &
                + (dzray_st(1) + 2._wp * (dzray_st(2) + dzray_st(3))) * one_o_6
          mray(jc, jray) = mray(jc, jray) &
                + (mray_st(1) + 2._wp * (mray_st(2) + mray_st(3))) * one_o_6
        END IF

#ifndef __msgwam1d
!       coslatray(jc, jray) = MAX(min_coslat, ABS(COS(latray(jc, jray))))
        coslatray(jc, jray) = coef_sphere*MAX(min_coslat, ABS(COS(latray(jc, jray)))) + coef_torus
        !TODO : Check if it is really faster than an explicit IF branching

        ! Horizontal/vertical wavenumber squared
        Kh2 = kray(jc,jray)**2 + lray(jc,jray)**2
#endif
        m2  = mray(jc,jray)**2

        ! Update vertical index of ray volume center-point
        ! for the next stage/substep
        IF (zray(jc,jray) > zray_st(0)) THEN
          ! Ray volume propagated upwards
          ikbase = 2
          DO jk = ik_st0-1, 2, -1
            IF (zhalf(jc,jk) > zray(jc,jray)) THEN
                ikbase = jk+1
                EXIT
            ENDIF
          ENDDO
        ELSE
          ! Ray volume propagated downwards
          ikbase = nlevp1
          DO jk = ik_st0, nlev  ! nlevp1
            IF (zhalf(jc,jk) <= zray(jc,jray)) THEN
                ikbase = jk
                EXIT
            ENDIF
          ENDDO
        END IF
        ! Even if the position of ray volume is out of domain during an
        ! intermediate RK stage, let it just continue to complete the
        ! time step. In this case the mean flow is extrapolated to these
        ! locations using the layer adjacent to the boundary.

      ENDDO ! jstage

      ! Take absolute value of ray extents in physical space so that 
      ! they are > 0 in the next time step. 
#ifndef __msgwam1d
      dlonray(jc,jray) = ABS(dlonray(jc,jray)*coslatray(jc,jray))   ! scaled
      dlatray(jc,jray) = ABS(dlatray(jc,jray))
#endif
      dzray(jc,jray)   = ABS(dzray(jc,jray))

      ! The lower bounds, darcraymin and dzraymin, are introduced for a numerical
      ! reason to avoid zero (or numerically too small) extents which can potentially
      ! induce numerical problems in further calculations. Those bounds should
      ! therefore be very small enough to tolerate as errors.
      ! The upper bounds are applied, for a safe side, considering unphysically
      ! abrupt, large expansion of an extent in ONE time step (e.g, to ~10x grid length).
      ! Although physically this may not occur (or only very rarely), we should
      ! guarantee it in a numerical reason, by explicitly applying some bounds.  
      ! (defined in the namelist)
#ifndef __msgwam1d
      dlonray(jc,jray) = MAX(darcraymin, MIN(darcraymax, dlonray(jc,jray)))
      dlatray(jc,jray) = MAX(darcraymin, MIN(darcraymax, dlatray(jc,jray)))
#endif
      dzray  (jc,jray) = MAX(dzraymin  , MIN(dzraymax  , dzray  (jc,jray)))

      ! Update the ray volume extent in k,l,m
#ifndef __msgwam1d
      dkray(jc,jray) = arealonk / dlonray(jc,jray)   ! scaled dlonray
      dlray(jc,jray) = arealatl / dlatray(jc,jray)
#endif
      dmray(jc,jray) = areazm   / dzray(jc,jray)

#ifndef __msgwam1d
      dlonray(jc,jray) = dlonray(jc,jray) / coslatray(jc,jray)
#endif

      ! Check if ray volume propagates out of vertical domain, and remove if so.
      ! In this respect the lower bound is set at z ~ 1 km because near the
      ! surface, due to the increasing resolution, the resolved flow might vary
      ! on short scales violating the WKB assumption. Some treatment like
      ! smoothing will be needed around the lower boundary if one would use the
      ! full domain to simulate upward reflection at the ground. The upper
      ! boundary is set such that 1) the gridded flux is zero at the uppermost
      ! layer and 2) ray volumes are removed when their entire volumes are
      ! within the sponge layer.
      l_outofdomain = (zray(jc,jray)<zhalf(jc,jk_lbnd) .OR. &
        &              zray(jc,jray)+0.5_wp*dzray(jc,jray)>z(jc,1) .OR. &
        &              zray(jc,jray)-0.5_wp*dzray(jc,jray)>z(jc,nrdmax(jg)-1))

      l_wavetooshort = ABS(mray(jc,jray)) > m_abs_max

      ! Diagnose whether ray volume is reflected twice. Such are removed too.
      l_refl_twice = branch*mray_st(0) > 0._wp .AND. branch*mray(jc,jray) < 0._wp

      ! Removal
      IF ( l_outofdomain .OR. l_refl_twice .OR. l_wavetooshort ) THEN
        iexist(jc,jray) = 0
        specid(jc,jray) = 0
        dens(jc,jray) = 0._wp
        CYCLE RAYLOOP
      END IF

      ! Calculate molecular viscosity of the sub-step
      K2kvisc1 = ( kvisc(jc, ikbase) + (kvisc(jc, ikbase - 1) - kvisc(jc, ikbase))    &
        &            * ((zray(jc, jray) - zhalf(jc, ikbase)) / dz(jc, ikbase - 1)) )  &
        &       * (Kh2 + m2)

      ! Increment the accumulated molecular viscosity
      acc_2K2kvisc = acc_2K2kvisc + K2kvisc0 + K2kvisc1

      ! Save molecular viscosity from the previous sub-step
      K2kvisc0 = K2kvisc1

    ENDDO ! itsubstep

#ifndef __msgwam1d
      IF (.NOT. is_plane_torus) THEN
        ! Treat lon, lat coordinates and latitudinal group velocity sign
        ! in case waves cross poles or cross the zero longitude
        ! Crossing the Greenwich 0 deg longitude
        IF ( ABS(lonray(jc,jray)) > pi ) THEN
          lonray(jc,jray) = lonray(jc,jray) - SIGN(pi2, lonray(jc,jray))
        ENDIF
        ! Pole crossing
        IF ( ABS(latray(jc,jray)) > pi_2 ) THEN
          latray(jc,jray) = - latray(jc,jray) + SIGN(pi, latray(jc,jray))
          lonray(jc,jray) =   lonray(jc,jray) - SIGN(pi, lonray(jc,jray))
          lray(jc,jray)   = - lray(jc,jray)
          kray(jc,jray)   = - kray(jc,jray)
        ENDIF
      ELSE ! is_plane_torus
        ! Treat x coordinates to ensure periodicity
        domxmin_ext = domxmin - 0.5_wp * edge_length_triangle
        domxmax_ext = domxmax + 0.5_wp * edge_length_triangle
        IF (lonray(jc,jray) < domxmin_ext) THEN
          lonray(jc,jray) = domxmax - ABS(domxmin_ext - lonray(jc,jray))
        ELSE IF (lonray(jc,jray) > domxmax_ext) THEN
          lonray(jc,jray) = domxmin + ABS(lonray(jc,jray) - domxmax_ext)
        ENDIF
        ! Treat y coordinates to ensure periodicity
        domymin_ext = domymin - height_triangle
        domymax_ext = domymax + height_triangle
        IF (latray(jc,jray) < domymin_ext) THEN
          latray(jc,jray) = domymax - ABS(domymin_ext - latray(jc,jray))
        ELSE IF (latray(jc,jray) > domymax_ext) THEN
          latray(jc,jray) = domymin + ABS(latray(jc,jray) - domymax_ext)
        ENDIF
      ENDIF ! is_plane_torus?
#endif

#ifndef __msgwam1d
    ! Horizontal propagation test
    IF (ltest_hprop .AND. l1ray .AND. (.NOT. is_plane_torus)) THEN
      OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(987,'(a,i4,2x,i4,2x,F14.5,2x,F14.5)') &
             'Prop wave position: jc, jb, lon, lat:', &
             jc, jb, lonray(jc,jray)*rad2deg, latray(jc,jray)*rad2deg 
      CLOSE(987)
    ELSE IF (ltest_hprop .AND. l1ray .AND. is_plane_torus) THEN
      OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
      WRITE(987,'(a,i4,2x,i4,2x,F14.5,2x,F14.5)') &
             'Prop wave position: jc, jb, x, y:', &
             jc, jb, lonray(jc,jray), latray(jc,jray)
      CLOSE(987)
    ENDIF
#endif

    IF (specid(jc,jray) >= 0) THEN
      ! In case we have an "old" ray volume

      ! Apply molecular viscosity
      dens(jc,jray) = dens(jc,jray) * EXP(-acc_2K2kvisc*dt_rk)
      IF (dens(jc,jray) <= eps_wadens) THEN
        iexist(jc,jray) = 0
        specid(jc,jray) = 0
        dens  (jc,jray) = 0.
        CYCLE RAYLOOP
      END IF

      ! For a ray volume whose center had crossed the source level but bottom had not yet
      IF (specid(jc,jray) > 0) THEN
        ! Once the ray-volume bottom crosses the source level (i.e., the volume is completely
        ! outside the ghost layer), then remove specid (--> 0).
        ! If a vertical reflection occurs before completely escaping the ghost layer,
        ! then delete the ray volume.
        IF ( zray(jc,jray) - 0.5_wp*dzray(jc,jray) > zhalf(jc,jk_active(jc,jray)) ) THEN
          specid(jc,jray) = 0
        ELSE IF (mray(jc,jray) > 0._wp) THEN
          iexist(jc,jray) = 0
          specid(jc,jray) = 0
          dens  (jc,jray) = 0.
          CYCLE RAYLOOP
        END IF
      END IF

      ! Update iexist for next steps
      ! iexist : half level closest to the ray volume center-point
      ! ikbase : half level below the ray volume center-point
      IF (zray(jc,jray) < z(jc,ikbase-1)) THEN
        iexist(jc,jray) = ikbase
      ELSE
        iexist(jc,jray) = ikbase-1
      END IF

#ifdef __test_awaypole
if(abs(latray(jc,jray))>85.*deg2rad) then
iexist(jc,jray) = 0 ; specid(jc,jray) = 0
endif
#endif

    END IF

    ENDDO RAYLOOP

  ENDDO ! jc

END SUBROUTINE propagate_wave
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE split_merge_volume( p_ray, p_spl, p_rwork, p_patch, c_nb, ncol_int, z_mc, z_ifc,  &
  &                            nrays_vacate, fc2, bvf2, gammash2, specid_start, specid_end,  &
  &                            jray_start, jray_end )

  TYPE(t_ray),            INTENT(INOUT) :: p_ray                ! the ray properties
  TYPE(t_spl),            INTENT(INOUT) :: p_spl                ! info for the splitting
  TYPE(t_ray),            INTENT(INOUT) :: p_rwork              ! work array
  TYPE(t_patch),   TARGET,INTENT(INOUT) :: p_patch              ! grid/patch info.
  TYPE(t_mgmgrid_nb),     INTENT(IN   ) :: c_nb(:,:)            ! cell's neighbor info
  INTEGER,                INTENT(IN   ) :: ncol_int             ! no. of interior grid columns
  REAL(wp),       POINTER,INTENT(IN   ) :: z_mc(:,:,:)          ! p_metrics%z_mc
  REAL(wp),       POINTER,INTENT(IN   ) :: z_ifc(:,:,:)         ! p_metrics%z_ifc
  INTEGER,                INTENT(IN   ) :: nrays_vacate
  REAL(wp),               INTENT(IN   ) :: fc2(:,:)
  REAL(wp),               INTENT(IN   ) :: bvf2(:,:,:)
  REAL(wp),               INTENT(IN   ) :: gammash2(:,:,:)
  INTEGER,                INTENT(IN   ) :: specid_start
  INTEGER,                INTENT(IN   ) :: specid_end
  INTEGER,                INTENT(IN   ) :: jray_start
  INTEGER,                INTENT(IN   ) :: jray_end

  INTEGER , PARAMETER :: tp = triprc

  INTEGER  :: nzray                                    ! number of full and half levels
  INTEGER  :: rl_start, rl_end_ext, rl_end_int
  INTEGER  :: i_startblk, i_endblk_ext, i_endblk_int   ! blocks
  INTEGER  :: i_startidx, i_endidx                     ! slices
  INTEGER  :: jc,jk,jb                                 ! slice and block indices
  INTEGER  :: jray
  INTEGER  :: j_ghost
  INTEGER  :: jr_o(p_patch%nlevp1+1)
  INTEGER  :: nr_z(p_patch%nlevp1)
  INTEGER  :: jr_coll(p_patch%nlevp1)
  INTEGER  :: jc_n, jb_n                          ! slice and block indices of neighbors
  INTEGER  :: jseg
  INTEGER  :: jn
  INTEGER  :: nnb
  INTEGER  :: jspec
  INTEGER  :: jiter
  INTEGER  :: nr_reduce_total
  INTEGER  :: nr_active
  INTEGER  :: work_jr_last(nrays_add(jg))
  REAL(tp) :: d(13)      ,  & ! distance between ray volume and cells
    &         ray_coslat ,  & ! cos(lat) of a ray volume (center)
    &         ray_sinlat ,  & ! sin(lat)
    &         ray_x_cart ,  & ! cos(lat)*cos(lon)
    &         ray_y_cart ,  & ! cos(lat)*sin(lon)
    &         ray_coslon ,  & ! cos(lon)
    &         ray_sinlon ,  & ! sin(lon)
    &         lon_we(2)  ,  & ! lon of split ray volumes
    &         lat_sn(2)       ! lat of split ray volumes
  REAL(wp) :: lon_we0(2) ,  &
    &         lat_sn0(2) ,  &
    &         ray_xyz(3)
  REAL(wp) :: coslat_dlon
  INTEGER  :: intriangle
#ifdef __mgm_oog
  INTEGER  :: nnb_n
#endif

  INTEGER , PARAMETER :: imode_skip    = 0
  INTEGER , PARAMETER :: imode_nosplit = -1
  INTEGER , PARAMETER :: imode_x       = 1   ! do not change this
  INTEGER , PARAMETER :: imode_y       = 3   ! do not change this

  ! Number of vertical levels
  nzray = p_patch%nlevp1

  rl_start   = grf_bdywidth_c+1
  rl_end_int = min_rlcell_int       ! Interior grids excluding halo
  rl_end_ext = min_rlcell_int - 2   ! Extended grids including halo

  i_startblk   = p_patch%cells%start_block(rl_start)
  i_endblk_int = p_patch%cells%end_block(rl_end_int)
  i_endblk_ext = p_patch%cells%end_block(rl_end_ext)

  ! Split in the horizontal (do not care out-going rays)
  ! Update location, even for non-split vol.

  ! Check the need for split in the horizontal and
  ! find the new base grid cell(s) of split/unsplit volumes,
  ! based on the distance between the ray volume
  ! and neighboring cells' center points (chord length is used)

  IF ( .NOT. is_plane_torus ) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jray, jseg, jn, i_startidx, i_endidx, nnb, lon_we, lat_sn, d,     &
!$OMP            coslat_dlon, ray_coslat, ray_coslon, ray_sinlon, ray_x_cart, ray_y_cart,  &
#ifdef __mgm_oog
!$OMP            nnb_n,  &
#endif
!$OMP            jc_n, jb_n)

  DO jb = i_startblk, i_endblk_ext
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk_ext, &
      &                 i_startidx, i_endidx, rl_start, rl_end_ext )
    DO jc = i_startidx, i_endidx

      nnb = 1 + c_nb(jc,jb)% n_neighbor

      DO jray = jray_start, jray_end

        IF ( p_ray% iexist(jc,jray,jb) == 0 .OR.  &
          &  p_ray% specid(jc,jray,jb) < 0 ) THEN
          p_spl% mode(jc,jray,jb) = imode_skip
          CYCLE
        END IF

        coslat_dlon = p_ray% coslat(jc,jray,jb)*p_ray% dlon(jc,jray,jb)

        IF ( p_ray% dlat(jc,jray,jb) < p_mgmgrid(jg)%darc_crit(jc,jb) .AND.  &
          &  coslat_dlon             < p_mgmgrid(jg)%darc_crit(jc,jb) .AND.  &
          &  coslat_dlon*p_ray% dlat(jc,jray,jb) < p_mgmgrid(jg)%area_crit(jc,jb) ) THEN
          ! do not split
          !
          ray_coslat = REAL(p_ray% coslat(jc,jray,jb), KIND=tp)
          ray_sinlat = REAL(SIGN(SQRT(1. - p_ray% coslat(jc,jray,jb)**2), &
            &                    p_ray% lat(jc,jray,jb)), KIND=tp)
          ray_x_cart = COS(REAL(p_ray% lon(jc,jray,jb), KIND=tp))*ray_coslat
          ray_y_cart = SIN(REAL(p_ray% lon(jc,jray,jb), KIND=tp))*ray_coslat
          d(:nnb) =  (c_nb(jc,jb)% x_cart(:nnb) - ray_x_cart)**2  &
            &      + (c_nb(jc,jb)% y_cart(:nnb) - ray_y_cart)**2  &
            &      + (c_nb(jc,jb)% sinlat(:nnb) - ray_sinlat)**2
          jn = MINLOC( d(:nnb) , DIM=1 )
          !
          IF (c_nb(jc,jb)%jcol_c(jn) > 0) THEN
            p_spl% jcol(1,jc,jray,jb) = c_nb(jc,jb)%jcol_c(jn)
            p_spl% mode(jc,jray,jb) = imode_nosplit
          ELSE   ! halo
            p_spl% mode(jc,jray,jb) = imode_skip
          END IF

        ELSE IF ( coslat_dlon > p_ray% dlat(jc,jray,jb) ) THEN
          ! split in X
          !
!YHK+ old
!         p_spl% aux(jc,jray,jb) = 0.25_wp*p_ray% dlon(jc,jray,jb)
!         lon_we = (/ REAL(p_ray% lon(jc,jray,jb) - p_spl% aux(jc,jray,jb), KIND=tp),  &
!           &         REAL(p_ray% lon(jc,jray,jb) + p_spl% aux(jc,jray,jb), KIND=tp) /)
!         ray_coslat = REAL(p_ray% coslat(jc,jray,jb), KIND=tp)
!         ray_sinlat = REAL(SIGN(SQRT(1. - p_ray% coslat(jc,jray,jb)**2), &
!           &                    p_ray% lat(jc,jray,jb)), KIND=tp)
!         !
!         DO jseg = 1, 2
!           ray_x_cart = COS(lon_we(jseg))*ray_coslat
!           ray_y_cart = SIN(lon_we(jseg))*ray_coslat
!           d(:nnb) =  (c_nb(jc,jb)% x_cart(:nnb) - ray_x_cart)**2  &
!             &      + (c_nb(jc,jb)% y_cart(:nnb) - ray_y_cart)**2  &
!             &      + (c_nb(jc,jb)% sinlat(:nnb) - ray_sinlat)**2
!           jn = MINLOC( d(:nnb) , DIM=1 )
!           p_spl% jcol(jseg,jc,jray,jb) = c_nb(jc,jb)%jcol_c(jn)
!         ENDDO
!         !
!         IF (p_spl% jcol(1,jc,jray,jb) /= p_spl% jcol(2,jc,jray,jb)) THEN
!           p_spl% mode(jc,jray,jb) = imode_x
!           ! update d_lon
!           p_ray% dlon(jc,jray,jb) = 0.5_wp*p_ray% dlon(jc,jray,jb)
!         ELSE IF (p_spl% jcol(1,jc,jray,jb) > 0) THEN
!           p_spl% mode(jc,jray,jb) = imode_nosplit   ! cancel the plan to split
!         ELSE   ! halo
!           p_spl% mode(jc,jray,jb) = imode_skip
!         END IF
!YHK- old
!YHK+ new
          ! sin(lat) = cos(d_sigma) sin(lat_0) ~= [1 - (1/2)(d_sigma)^2] sin(lat_0)
          !   where d_sigma = (1/4)d_lon cos(lat_0)
          ray_sinlat = REAL(SIGN( (1. - 0.03125_wp * coslat_dlon**2)         &
            &                     *SQRT(1. - p_ray% coslat(jc,jray,jb)**2),  &
            &               p_ray% lat(jc,jray,jb) ), KIND=tp)
          ray_coslat = SQRT(1. - ray_sinlat**2)
          ! lon displacement for split centers:
          !   tan(lon') = tan(d_sigma)/cos(lat_0) ~= (1/4)d_lon
          p_spl% aux(jc,jray,jb) = ATAN(0.25_wp*p_ray% dlon(jc,jray,jb))
          !   approximate of atan(x)
          !     atan(x) ~= x / (1 + 0.273 x^2) for |x| < 1   (within ~1.1% errors)
          !     atan(x) = pi/2 - atan(1/x)     otherwise
!         IF (p_ray% dlon(jc,jray,jb) < 4._wp) THEN
!           p_spl% aux(jc,jray,jb) = p_ray% dlon(jc,jray,jb)  &
!             &                     /(4._wp + 0.06831_wp*p_ray% dlon(jc,jray,jb)**2)
!         ELSE
!           p_spl% aux(jc,jray,jb) = pi_2 - 1._wp/( p_ray% dlon(jc,jray,jb)  &
!             &                         *(0.25_wp + 1.093_wp/p_ray% dlon(jc,jray,jb)**2) )
!         END IF
          !
          lon_we = (/ REAL(p_ray% lon(jc,jray,jb) - p_spl% aux(jc,jray,jb), KIND=tp),  &
            &         REAL(p_ray% lon(jc,jray,jb) + p_spl% aux(jc,jray,jb), KIND=tp) /)
          !
          DO jseg = 1, 2
            ray_x_cart = COS(lon_we(jseg))*ray_coslat
            ray_y_cart = SIN(lon_we(jseg))*ray_coslat
            d(:nnb) =  (c_nb(jc,jb)% x_cart(:nnb) - ray_x_cart)**2  &
              &      + (c_nb(jc,jb)% y_cart(:nnb) - ray_y_cart)**2  &
              &      + (c_nb(jc,jb)% sinlat(:nnb) - ray_sinlat)**2
            jn = MINLOC( d(:nnb) , DIM=1 )
            p_spl% jcol(jseg,jc,jray,jb) = c_nb(jc,jb)%jcol_c(jn)
#ifdef __mgm_oog
IF (jn /= 1) THEN
jc_n = c_nb(jc,jb)% jcb_neighbor(1,jn)  ;  jb_n = c_nb(jc,jb)% jcb_neighbor(2,jn)
nnb_n = 1 + c_nb(jc_n,jb_n)% n_neighbor
d(1) = d(jn)
d(2:nnb_n) =  (c_nb(jc_n,jb_n)% x_cart(2:nnb_n) - ray_x_cart)**2  &
  &         + (c_nb(jc_n,jb_n)% y_cart(2:nnb_n) - ray_y_cart)**2  &
  &         + (c_nb(jc_n,jb_n)% sinlat(2:nnb_n) - ray_sinlat)**2
IF ( MINLOC( d(:nnb_n) , DIM=1 ) /= 1 ) THEN
 print*,'YHK-out-of-gridx [km]', 6371.*SQRT(d(1)), 6371.*SQRT(d(MINLOC( d(:nnb_n), DIM=1 )))
END IF
END IF
#endif
          ENDDO
          !
          IF (p_spl% jcol(1,jc,jray,jb) /= p_spl% jcol(2,jc,jray,jb)) THEN
            p_spl% mode(jc,jray,jb) = imode_x
            ! rotation matrix
            p_spl% rot_c(jc,jray,jb) = p_ray% coslat(jc,jray,jb) / ray_coslat
            p_spl% rot_s(jc,jray,jb) = 0.25_wp*p_ray% dlon(jc,jray,jb)  &   ! tan(lon')
              &                       * p_spl% rot_c(jc,jray,jb) * ray_sinlat
            ! update lat and its dependencies and d_lon
            p_ray% coslat(jc,jray,jb) = REAL(ray_coslat, KIND=wp)
            p_ray% dlon  (jc,jray,jb) = 0.5_wp*coslat_dlon / p_ray% coslat(jc,jray,jb)
            p_ray% lat   (jc,jray,jb) = REAL(ASIN(ray_sinlat), KIND=wp)  ! approximate this
          ELSE IF (p_spl% jcol(1,jc,jray,jb) > 0) THEN
            p_spl% mode(jc,jray,jb) = imode_nosplit   ! cancel the plan to split
          ELSE   ! halo
            p_spl% mode(jc,jray,jb) = imode_skip
          END IF
!YHK- new

        ELSE
          ! split in Y
          !
          lat_sn = (/ REAL(p_ray% lat(jc,jray,jb) - 0.25_wp*p_ray% dlat(jc,jray,jb), KIND=tp),  &
            &         REAL(p_ray% lat(jc,jray,jb) + 0.25_wp*p_ray% dlat(jc,jray,jb), KIND=tp) /)
          ray_coslon = COS(REAL(p_ray% lon(jc,jray,jb), KIND=tp))
          ray_sinlon = SIN(REAL(p_ray% lon(jc,jray,jb), KIND=tp))
          !
          DO jseg = 1, 2
            ray_coslat = COS(lat_sn(jseg))
            ray_sinlat = SIGN(SQRT(1. - ray_coslat**2), lat_sn(jseg))
            ray_x_cart = ray_coslon*ray_coslat
            ray_y_cart = ray_sinlon*ray_coslat
            d(:nnb) =  (c_nb(jc,jb)% x_cart(:nnb) - ray_x_cart)**2  &
              &      + (c_nb(jc,jb)% y_cart(:nnb) - ray_y_cart)**2  &
              &      + (c_nb(jc,jb)% sinlat(:nnb) - ray_sinlat)**2
            jn = MINLOC( d(:nnb) , DIM=1 )
            p_spl% jcol(jseg,jc,jray,jb) = c_nb(jc,jb)%jcol_c(jn)
#ifdef __mgm_oog
IF (jn /= 1) THEN
jc_n = c_nb(jc,jb)% jcb_neighbor(1,jn)  ;  jb_n = c_nb(jc,jb)% jcb_neighbor(2,jn)
nnb_n = 1 + c_nb(jc_n,jb_n)% n_neighbor
d(1) = d(jn)
d(2:nnb_n) =  (c_nb(jc_n,jb_n)% x_cart(2:nnb_n) - ray_x_cart)**2  &
  &         + (c_nb(jc_n,jb_n)% y_cart(2:nnb_n) - ray_y_cart)**2  &
  &         + (c_nb(jc_n,jb_n)% sinlat(2:nnb_n) - ray_sinlat)**2
IF ( MINLOC( d(:nnb_n) , DIM=1 ) /= 1 ) THEN
 print*,'YHK-out-of-gridy [km]', 6371.*SQRT(d(1)), 6371.*SQRT(d(MINLOC( d(:nnb_n), DIM=1 )))
END IF
END IF
#endif
          ENDDO
          !
          IF (p_spl% jcol(1,jc,jray,jb) /= p_spl% jcol(2,jc,jray,jb)) THEN
            p_spl% mode(jc,jray,jb) = imode_y
            ! update d_lat
            p_ray% dlat(jc,jray,jb) = 0.5_wp*p_ray% dlat(jc,jray,jb)
          ELSE IF (p_spl% jcol(1,jc,jray,jb) > 0) THEN
            p_spl% mode(jc,jray,jb) = imode_nosplit   ! cancel the plan to split
          ELSE   ! halo
            p_spl% mode(jc,jray,jb) = imode_skip
          END IF

        END IF  ! 3 cases of the volume splitting

      ENDDO  ! jc
    ENDDO  ! jray
  ENDDO  ! jb

!$OMP END DO
!$OMP END PARALLEL

  ELSE  ! is_plane_torus == .TRUE.   ! TODO: not tested yet

    ! Here the 'inside_triangle_torus' scheme is used. In case a ray-volume segment is
    ! outside of all the considered cells (13 cells: original + neighbor cells),
    ! its grid-cell index cannot be obtained. Currently, such a segment will be removed.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jray, jseg, jn, i_startidx, i_endidx, lon_we0, lat_sn0,  &
!$OMP            coslat_dlon, ray_xyz, intriangle)

  DO jb = i_startblk, i_endblk_ext
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk_ext, &
      &                 i_startidx, i_endidx, rl_start, rl_end_ext )
    DO jc = i_startidx, i_endidx
      DO jray = jray_start, jray_end

        IF ( p_ray% iexist(jc,jray,jb) == 0 .OR.  &
          &  p_ray% specid(jc,jray,jb) < 0 ) THEN
          p_spl% mode(jc,jray,jb) = imode_skip
          CYCLE
        END IF

        coslat_dlon = p_ray% coslat(jc,jray,jb)*p_ray% dlon(jc,jray,jb)

        IF ( p_ray% dlat(jc,jray,jb) < p_mgmgrid(jg)%darc_crit(jc,jb) .AND.  &
          &  coslat_dlon             < p_mgmgrid(jg)%darc_crit(jc,jb) .AND.  &
          &  coslat_dlon*p_ray% dlat(jc,jray,jb) < p_mgmgrid(jg)%area_crit(jc,jb) ) THEN
          ! do not split
          !
          ray_xyz(1) = p_ray% lon(jc,jray,jb)
          ray_xyz(2) = p_ray% lat(jc,jray,jb)

          p_spl% mode(jc,jray,jb) = imode_skip   ! for case it is beyond the neighbors
          DO jn = 1, 1 + c_nb(jc,jb)% n_neighbor
            intriangle = inside_triangle_torus( ray_xyz,   &
              &  c_nb(jc,jb)% v_cart(1,jn), c_nb(jc,jb)% v_cart(2,jn), c_nb(jc,jb)% v_cart(3,jn) )
            IF ( intriangle >= 0 ) THEN
              IF (c_nb(jc,jb)%jcol_c(jn) > 0) THEN
                p_spl% jcol(1,jc,jray,jb) = c_nb(jc,jb)%jcol_c(jn)
                p_spl% mode(jc,jray,jb) = imode_nosplit
              ELSE   ! halo
                p_spl% mode(jc,jray,jb) = imode_skip
              END IF
              EXIT
            END IF
          ENDDO

        ELSE IF ( coslat_dlon > p_ray% dlat(jc,jray,jb) ) THEN
          ! split in X
          !
          p_spl% aux(jc,jray,jb) = 0.25_wp*p_ray% dlon(jc,jray,jb)
          lon_we0 = (/ p_ray% lon(jc,jray,jb) - p_spl% aux(jc,jray,jb),  &
            &          p_ray% lon(jc,jray,jb) + p_spl% aux(jc,jray,jb) /)
          ray_xyz(2) = p_ray% lat(jc,jray,jb)

          p_spl% jcol(1:2,jc,jray,jb) = -999   ! for case it is beyond the neighbors
          DO jseg = 1, 2
            ray_xyz(1) = lon_we0(jseg)
            DO jn = 1, 1 + c_nb(jc,jb)% n_neighbor
              intriangle = inside_triangle_torus( ray_xyz,   &
                &  c_nb(jc,jb)% v_cart(1,jn), c_nb(jc,jb)% v_cart(2,jn), c_nb(jc,jb)% v_cart(3,jn) )
              IF ( intriangle >= 0 ) THEN
                p_spl% jcol(jseg,jc,jray,jb) = c_nb(jc,jb)%jcol_c(jn)
                EXIT
              END IF
            ENDDO
          ENDDO
          IF (p_spl% jcol(1,jc,jray,jb) /= p_spl% jcol(2,jc,jray,jb)) THEN
            p_spl% mode(jc,jray,jb) = imode_x
            ! update d_lon
            p_ray% dlon(jc,jray,jb) = 0.5_wp*p_ray% dlon(jc,jray,jb)
          ELSE IF (p_spl% jcol(1,jc,jray,jb) > 0) THEN
            p_spl% mode(jc,jray,jb) = imode_nosplit   ! cancel the plan to split
          ELSE   ! halo
            p_spl% mode(jc,jray,jb) = imode_skip
          END IF

        ELSE
          ! split in Y
          !
          lat_sn0 = (/ p_ray% lat(jc,jray,jb) - 0.25_wp*p_ray% dlat(jc,jray,jb),  &
            &          p_ray% lat(jc,jray,jb) + 0.25_wp*p_ray% dlat(jc,jray,jb) /)
          ray_xyz(1) = p_ray% lon(jc,jray,jb)

          p_spl% jcol(1:2,jc,jray,jb) = -999   ! for case it is beyond the neighbors
          DO jseg = 1, 2
            ray_xyz(2) = lat_sn(jseg)
            DO jn = 1, 1 + c_nb(jc,jb)% n_neighbor
              intriangle = inside_triangle_torus( ray_xyz,   &
                &  c_nb(jc,jb)% v_cart(1,jn), c_nb(jc,jb)% v_cart(2,jn), c_nb(jc,jb)% v_cart(3,jn) )
              IF ( intriangle >= 0 ) THEN
                p_spl% jcol(jseg,jc,jray,jb) = c_nb(jc,jb)%jcol_c(jn)
                EXIT
              END IF
            ENDDO
          ENDDO
          IF (p_spl% jcol(1,jc,jray,jb) /= p_spl% jcol(2,jc,jray,jb)) THEN
            p_spl% mode(jc,jray,jb) = imode_y
            ! update d_lat
            p_ray% dlat(jc,jray,jb) = 0.5_wp*p_ray% dlat(jc,jray,jb)
          ELSE IF (p_spl% jcol(1,jc,jray,jb) > 0) THEN
            p_spl% mode(jc,jray,jb) = imode_nosplit   ! cancel the plan to split
          ELSE   ! halo
            p_spl% mode(jc,jray,jb) = imode_skip
          END IF

        END IF  ! 3 cases of the volume splitting

      ENDDO  ! jc
    ENDDO  ! jray
  ENDDO  ! jb

!$OMP END DO
!$OMP END PARALLEL

  END IF


  ! Below, tasks are done at the interior grids only, with the access to halo cells
  ! if they are neighbor.

  ! rv at jr_last must be untouched: rather having reserved index

  DO jb = i_startblk, i_endblk_int
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk_int, &
      &                 i_startidx, i_endidx, rl_start, rl_end_int )
    DO jc = i_startidx, i_endidx
      !
      ! Gather ghost ray volumes into the work array from the index (nrays_vacate + 1)
      !   before splitting/merging non-ghost volumes.
      !
      j_ghost = jray_start - 1 + nrays_vacate
      !
      work_jr_last(specid_start:specid_end) = 0
      DO jray = jray_start, jray_end
        IF (p_ray% specid(jc,jray,jb) >= 0)  CYCLE
        !
        j_ghost = j_ghost + 1
        !
        p_rwork% iexist   (jc,j_ghost,jb) = p_ray% iexist   (jc,jray,jb)
        p_rwork% specid   (jc,j_ghost,jb) = p_ray% specid   (jc,jray,jb)
        p_rwork% jk_active(jc,j_ghost,jb) = p_ray% jk_active(jc,jray,jb)
        p_rwork% lon      (jc,j_ghost,jb) = p_ray% lon      (jc,jray,jb)
        p_rwork% lat      (jc,j_ghost,jb) = p_ray% lat      (jc,jray,jb)
        p_rwork% z        (jc,j_ghost,jb) = p_ray% z        (jc,jray,jb)
        p_rwork% dlon     (jc,j_ghost,jb) = p_ray% dlon     (jc,jray,jb)
        p_rwork% dlat     (jc,j_ghost,jb) = p_ray% dlat     (jc,jray,jb)
        p_rwork% dz       (jc,j_ghost,jb) = p_ray% dz       (jc,jray,jb)
        p_rwork% coslat   (jc,j_ghost,jb) = p_ray% coslat   (jc,jray,jb)
        p_rwork% k        (jc,j_ghost,jb) = p_ray% k        (jc,jray,jb)
        p_rwork% l        (jc,j_ghost,jb) = p_ray% l        (jc,jray,jb)
        p_rwork% m        (jc,j_ghost,jb) = p_ray% m        (jc,jray,jb)
        p_rwork% dk       (jc,j_ghost,jb) = p_ray% dk       (jc,jray,jb)
        p_rwork% dl       (jc,j_ghost,jb) = p_ray% dl       (jc,jray,jb)
        p_rwork% dm       (jc,j_ghost,jb) = p_ray% dm       (jc,jray,jb)
        p_rwork% wadens   (jc,j_ghost,jb) = p_ray% wadens   (jc,jray,jb)
        !
        ! Update jr_last
        !
        jspec = ABS(p_ray% specid(jc,jray,jb))
        IF ( p_ray% jr_last(jc,jspec,jb) == jray )  work_jr_last(jspec) = j_ghost
      ENDDO  ! jray
      ! update of jr_last in advance (it does not affect the other variables in this loop)
      p_ray% jr_last(jc,specid_start:specid_end,jb) = work_jr_last(specid_start:specid_end)

      jr_coll(:) = 0

      DO jn = 1, 1+c_nb(jc,jb)% n_neighbor
        jc_n = c_nb(jc,jb)% jcb_neighbor(1,jn)
        jb_n = c_nb(jc,jb)% jcb_neighbor(2,jn)

        DO jray = jray_start, jray_end

          IF ( p_spl% mode(jc_n,jray,jb_n) == imode_skip )  CYCLE
             ! iexist == 0 or specid < 0 (ghost ray volume)
             ! or volume already collected to another column

          IF ( p_spl% mode(jc_n,jray,jb_n) == imode_nosplit ) THEN   ! not split horizontally
            !
            IF (p_spl% jcol(1,jc_n,jray,jb_n) == c_nb(jc,jb)% jcol_c(1)) THEN
              CALL collect_vol( p_ray, p_patch, z_ifc(:,:,jb), z_mc(:,:,jb), z_mc(:,:,jb_n),  &
                &               jray, jc,jb, jc_n,jb_n, jr_coll, 0, p_spl )
              p_spl% mode(jc_n,jray,jb_n) = imode_skip   ! to skip from other grid-columns
            END IF
            !
          ELSE   ! split in either X or Y   # spl_jcol(1) /= spl_jcol(2)
            !
            IF (p_spl% jcol(1,jc_n,jray,jb_n) == c_nb(jc,jb)% jcol_c(1)) THEN
              CALL collect_vol( p_ray, p_patch, z_ifc(:,:,jb), z_mc(:,:,jb), z_mc(:,:,jb_n),  &
                &               jray, jc,jb, jc_n,jb_n, jr_coll,                               &
                &               p_spl% mode(jc_n,jray,jb_n), p_spl )
                                       ! passing imode_x or imode_y
              IF (p_spl% jcol(2,jc_n,jray,jb_n) > 0) THEN
                p_spl% jcol(1,jc_n,jray,jb_n) = 0
              ELSE
                p_spl% mode(jc_n,jray,jb_n) = imode_skip
              END IF
            ELSE IF (p_spl% jcol(2,jc_n,jray,jb_n) == c_nb(jc,jb)% jcol_c(1)) THEN
              CALL collect_vol( p_ray, p_patch, z_ifc(:,:,jb), z_mc(:,:,jb), z_mc(:,:,jb_n),  &
                &               jray, jc,jb, jc_n,jb_n, jr_coll,                               &
                &               p_spl% mode(jc_n,jray,jb_n) + 1, p_spl )
                                       ! passing (imode_x or imode_y) + 1
              IF (p_spl% jcol(1,jc_n,jray,jb_n) > 0) THEN
                p_spl% jcol(2,jc_n,jray,jb_n) = 0
              ELSE
                p_spl% mode(jc_n,jray,jb_n) = imode_skip
              END IF
            END IF
            !
          END IF

        ENDDO  ! jray

      ENDDO  ! jn

      nr_z(:) = jr_coll(:)

      ! In order to constrain the number of ray volumes:
      ! Merge ray volumes located in the grid column (jc,jb), which have been collected in
      ! 'ray_coll'. The number of the collected in each level has been saved in 'nr_z(jk)'.

      IF (imethod_merge > 0) THEN
        DO jiter = 1, 5
          nr_active = SUM(nr_z(1:nzray))
          nr_reduce_total = j_ghost + nr_active - jray_end
          IF (nr_reduce_total <= 0)  EXIT
          !
          ! quite rarely jiter = 2 here ;  jiter > 2 : not yet observed
          CALL merge_vol( jiter, nzray, nr_active, nr_reduce_total,        & ! IN
            &             p_patch%cells%center(jc,jb)%lon,                 & ! IN
            &             p_patch%cells%center(jc,jb)%lat,                 & ! IN
            &             fc2(jc,jb), bvf2(:,:,jb), gammash2(:,:,jb), jc,  & ! IN
            &             nr_z, ray_coll )                                   ! INOUT
          IF (jiter >= 3)  print*, 'YHK- merge done severely:', jiter
        ENDDO
      ELSE
        nr_active = SUM(nr_z(1:nzray))
        nr_reduce_total = j_ghost + nr_active - jray_end
        IF (nr_reduce_total > 0) THEN
          CALL remove_vol( nzray, nr_active, nr_reduce_total,               & ! IN
            &              fc2(jc,jb), bvf2(:,:,jb), gammash2(:,:,jb), jc,  & ! IN
            &              nr_z, ray_coll )                                   ! INOUT
        END IF
      END IF

      ! Put the results into the original ray-volume array 'p_ray'.

      jr_o(1) = j_ghost
      DO jk = 1, nzray
        jr_o(jk+1) = jr_o(jk) + nr_z(jk)
      ENDDO

      IF (jr_o(nzray+1) > jray_end) THEN
        print*, 'YHK- merge (tot): before/after/aim/ghost',  &
          &  nr_active, jr_o(nzray+1) - jr_o(1), nr_active - nr_reduce_total, j_ghost - jray_start + 1
        STOP
      END IF

      DO jk = 1, nzray
        IF (nr_z(jk) == 0)  CYCLE
        p_rwork% iexist   (jc,jr_o(jk)+1:jr_o(jk+1),jb) = jk
        p_rwork% specid   (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% specid
        p_rwork% jk_active(jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% jk_active
        p_rwork% lon      (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% lon
        p_rwork% lat      (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% lat
        p_rwork% z        (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% z
        p_rwork% dlon     (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% dlon
        p_rwork% dlat     (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% dlat
        p_rwork% dz       (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% dz
        p_rwork% coslat   (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% coslat
        p_rwork% k        (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% k
        p_rwork% l        (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% l
        p_rwork% m        (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% m
        p_rwork% dk       (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% dk
        p_rwork% dl       (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% dl
        p_rwork% dm       (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% dm
        p_rwork% wadens   (jc,jr_o(jk)+1:jr_o(jk+1),jb) = ray_coll(1:nr_z(jk),jk)% wadens
      ENDDO  ! jk
      IF (jr_o(nzray+1) < jray_end) THEN
        p_rwork% iexist(jc,jr_o(nzray+1)+1:jray_end,jb) = 0
        p_rwork% specid(jc,jr_o(nzray+1)+1:jray_end,jb) = 0
        p_rwork% wadens(jc,jr_o(nzray+1)+1:jray_end,jb) = 0.
      END IF

    ENDDO  ! jc
  ENDDO  ! jb

  ! Copy to p_ray
  jray = jray_start + nrays_vacate
  IF (nrays_vacate /= 0) THEN
    DO jb = i_startblk, i_endblk_int
      p_ray%iexist(:,jray_start:jray-1,jb) = 0
      p_ray%specid(:,jray_start:jray-1,jb) = 0
      p_ray%wadens(:,jray_start:jray-1,jb) = 0.
    ENDDO
  END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb)

  DO jb = i_startblk, i_endblk_int
    p_ray% iexist   (:,jray:jray_end,jb) = p_rwork% iexist   (:,jray:jray_end,jb)
    p_ray% specid   (:,jray:jray_end,jb) = p_rwork% specid   (:,jray:jray_end,jb)
    p_ray% jk_active(:,jray:jray_end,jb) = p_rwork% jk_active(:,jray:jray_end,jb)
    p_ray% lon      (:,jray:jray_end,jb) = p_rwork% lon      (:,jray:jray_end,jb)
    p_ray% lat      (:,jray:jray_end,jb) = p_rwork% lat      (:,jray:jray_end,jb)
    p_ray% z        (:,jray:jray_end,jb) = p_rwork% z        (:,jray:jray_end,jb)
    p_ray% dlon     (:,jray:jray_end,jb) = p_rwork% dlon     (:,jray:jray_end,jb)
    p_ray% dlat     (:,jray:jray_end,jb) = p_rwork% dlat     (:,jray:jray_end,jb)
    p_ray% dz       (:,jray:jray_end,jb) = p_rwork% dz       (:,jray:jray_end,jb)
    p_ray% coslat   (:,jray:jray_end,jb) = p_rwork% coslat   (:,jray:jray_end,jb)
    p_ray% k        (:,jray:jray_end,jb) = p_rwork% k        (:,jray:jray_end,jb)
    p_ray% l        (:,jray:jray_end,jb) = p_rwork% l        (:,jray:jray_end,jb)
    p_ray% m        (:,jray:jray_end,jb) = p_rwork% m        (:,jray:jray_end,jb)
    p_ray% dk       (:,jray:jray_end,jb) = p_rwork% dk       (:,jray:jray_end,jb)
    p_ray% dl       (:,jray:jray_end,jb) = p_rwork% dl       (:,jray:jray_end,jb)
    p_ray% dm       (:,jray:jray_end,jb) = p_rwork% dm       (:,jray:jray_end,jb)
    p_ray% wadens   (:,jray:jray_end,jb) = p_rwork% wadens   (:,jray:jray_end,jb)
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  IF (.TRUE.) THEN
!  IF (.FALSE.) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jray, jc, i_startidx, i_endidx)
    DO jb = i_startblk, i_endblk_int
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk_int, &
        &                 i_startidx, i_endidx, rl_start, rl_end_int )
      DO jray = jray_start, jray_end
        DO jc = i_startidx, i_endidx
          IF ( p_ray% iexist(jc,jray,jb) == 0 .OR. p_ray% specid(jc,jray,jb) < 0 )  CYCLE
          !
          IF ( p_ray% dz(jc,jray,jb) > 2.5_wp*p_mgmgrid(jg)% dz_crit(p_ray% iexist(jc,jray,jb)) ) THEN
            WRITE(message_text,'(a,2i5)') 'YHK: d_z           /(d_crit) [%]', INT( 100._wp*  &
              &     p_ray% dz(jc,jray,jb)/p_mgmgrid(jg)% dz_crit(p_ray% iexist(jc,jray,jb)) ),  &
              &     p_ray% iexist(jc,jray,jb)
            CALL message('', TRIM(message_text))
          END IF
          IF ( p_ray% dlat(jc,jray,jb) > 2._wp*p_mgmgrid(jg)% darc_crit(jc,jb) ) THEN
            WRITE(message_text,'(a,2i5)') 'YHK: d_lat         /(d_crit) [%]', INT( 100._wp*  &
              &     p_ray% dlat(jc,jray,jb)/p_mgmgrid(jg)% darc_crit(jc,jb) ),  &
              &     p_ray% iexist(jc,jray,jb)
            CALL message('', TRIM(message_text))
          END IF
          IF ( p_ray% coslat(jc,jray,jb) * p_ray% dlon(jc,jray,jb)  &
               > 2._wp*p_mgmgrid(jg)% darc_crit(jc,jb) ) THEN
            WRITE(message_text,'(a,2i5)') 'YHK: cos(lat)*d_lon/(d_crit) [%]', INT( 100._wp*  &
              &     p_ray% coslat(jc,jray,jb) * p_ray% dlon(jc,jray,jb)  &
              &     /p_mgmgrid(jg)% darc_crit(jc,jb) ),  &
              &     p_ray% iexist(jc,jray,jb)
            CALL message('', TRIM(message_text))
          END IF
          !
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END IF

END SUBROUTINE split_merge_volume

SUBROUTINE collect_vol( p_ray, p_patch, z_ifc, z_mc, z_mc_n, jray, jc, jb, jc_n, jb_n, i,  &
  &                     iloc, p_spl )

  TYPE(t_ray),            INTENT(IN   ) :: p_ray
  TYPE(t_patch),   TARGET,INTENT(IN   ) :: p_patch
  REAL(wp),               INTENT(IN   ) :: z_ifc(:,:)
  REAL(wp),               INTENT(IN   ) :: z_mc(:,:)
  REAL(wp),               INTENT(IN   ) :: z_mc_n(:,:)
  INTEGER ,               INTENT(IN   ) :: jray
  INTEGER ,               INTENT(IN   ) :: jc, jb
  INTEGER ,               INTENT(IN   ) :: jc_n, jb_n
  INTEGER ,               INTENT(INOUT) :: i(:)
  INTEGER ,               INTENT(IN   ) :: iloc
  TYPE(t_spl),            INTENT(IN   ) :: p_spl

  INTEGER , POINTER :: jkr
  INTEGER , TARGET  :: jk0
  INTEGER , TARGET  :: jk_ul(2)
  INTEGER           :: jkmp1
  INTEGER           :: jk
  INTEGER           :: ivol
  REAL(wp)          :: z_ul(2)
  REAL(wp)          :: d
  LOGICAL           :: l_1vol_z

  INTEGER           :: jl   ! counter

  REAL(wp), PARAMETER :: two_pi = pi*2._wp
  REAL(wp), PARAMETER :: c(2,0:4) = RESHAPE( (/         &
!                                  1 (lon)   2 (lat)
    &                               0.0_wp ,  0.0_wp ,  & ! 0
    &                              -1.0_wp ,  0.0_wp ,  & ! 1
    &                               1.0_wp ,  0.0_wp ,  & ! 2
    &                               0.0_wp , -0.5_wp ,  & ! 3
    &                               0.0_wp ,  0.5_wp    & ! 4
    &                              /), (/2,5/) )

  ! It is assumed that 'iexist' has already been updated correctly after
  ! the propagation calculation, if the base grid column is unchanged.
  jk0 = p_ray% iexist(jc_n,jray,jb_n)
  jk0m1 = MAX(1, jk0 - 1)

  ! If the vertical grids are different between the grid columns, update jk0.
  ! jk0 must satisfy both of the following conditions:
  !   #  z_mc(jc,jk0) <= zray                   if  jk0 < nlevp1
  !   #                  zray < z_mc(jc,jk0-1)  if  jk0 > 1

  ! To gaurantee that the index does not change, both levels above and below the ray
  ! volume must be consistent between the neighbor cells (identical altitudes at the
  ! lower-level coordinate does not always guarantee the identity at the upper level,
  ! for ICON's vertical coordinate system).

  IF ( z_mc_n(jc_n,jk0) /= z_mc(jc,jk0) .OR. z_mc_n(jc_n,jk0m1) /= z_mc(jc,jk0m1) ) THEN
    !
    ! Some conditions used at the end of the integration (in propagate_wave) need
    ! to be applied again here, as the vertical grids are changed.
    !  - lower bound used in propagate_wave
!   IF (p_ray% z(jc_n,jray,jb_n) < z_ifc(jc, p_patch% nlevp1)) THEN
    IF (p_ray% z(jc_n,jray,jb_n) < z_ifc(jc, p_patch% nlevp1 - nlev_lbnd(jg))) THEN  ! outofdomain
      RETURN
    END IF
    !
!   CALL get_jkhalf_closest_gen( jk0, p_ray% z(jc_n,jray,jb_n), p_patch% nlev, z_mc, jc )
    IF (jk0 < p_patch% nlevp1) THEN
      IF (z_mc(jc,jk0) > p_ray% z(jc_n,jray,jb_n)) THEN
        IF (z_mc(jc, p_patch% nlev) > p_ray% z(jc_n,jray,jb_n)) THEN
          jk0 = p_patch% nlevp1
        ELSE
          DO jk = jk0+1, p_patch% nlev
            IF (z_mc(jc,jk) <= p_ray% z(jc_n,jray,jb_n)) THEN
              jk0 = jk  ;  EXIT
            END IF
          ENDDO
        END IF
      END IF
    END IF
    IF (jk0 > 1) THEN
      IF (z_mc(jc,jk0-1) <= p_ray% z(jc_n,jray,jb_n)) THEN
        IF (z_mc(jc,1) <= p_ray% z(jc_n,jray,jb_n)) THEN   ! may not happen by l_outofdomain
          jk0 = 1
        ELSE
          DO jk = jk0-2, 1, -1
            IF (z_mc(jc,jk) > p_ray% z(jc_n,jray,jb_n)) THEN
              jk0 = jk+1  ;  EXIT
            END IF
          ENDDO
        END IF
      END IF
    END IF
    !
  END IF

  l_1vol_z = p_ray% dz(jc_n,jray,jb_n) < p_mgmgrid(jg)% dz_crit(jk0)

  IF ( .NOT. l_1vol_z ) THEN
    z_ul = (/ p_ray% z(jc_n,jray,jb_n) + 0.25_wp*p_ray% dz(jc_n,jray,jb_n),  &
      &       p_ray% z(jc_n,jray,jb_n) - 0.25_wp*p_ray% dz(jc_n,jray,jb_n) /)

    jk_ul(1) = 2
    DO jl = jk0, 3, -1
      IF (z_mc(jc,jl-1) > z_ul(1)) THEN
        jk_ul(1) = jl
        EXIT
      ENDIF
    ENDDO
    jk_ul(2) = p_patch% nlev   ! nlevp1 - 1
    DO jl = jk0, p_patch% nlev-1
      IF (z_mc(jc,jl) <= z_ul(2)) THEN
        jk_ul(2) = jl
        EXIT
      ENDIF
    ENDDO

    l_1vol_z = jk_ul(1) == jk_ul(2)
  END IF

  IF ( l_1vol_z ) THEN
    ! do not split in the vertical

    jkr => jk0

    i(jkr) = i(jkr) + 1

    ray_coll(i(jkr),jkr)% dz = p_ray% dz(jc_n,jray,jb_n)
    ray_coll(i(jkr),jkr)% z  = p_ray% z (jc_n,jray,jb_n)

  ELSE
    ! split in Z

    ivol = 1

    ! No need to apply the upper boundary conditions again, as they have been applied
    ! to the altitude of top of ray volumes (cf. center of volumes)

!   IF ( upper volume needs to be killed ) THEN
!     ivol = 2   ;   l_1vol_z = .True.
!   ELSE IF ( lower volume needs to be killed ) THEN
!     l_1vol_z = .True.
!   END IF

    jkr => jk_ul(ivol)

    i(jkr) = i(jkr) + 1

    ray_coll(i(jkr),jkr)% dz = p_ray% dz(jc_n,jray,jb_n)*0.5_wp
    ray_coll(i(jkr),jkr)% z  = z_ul(ivol)

  END IF

  ray_coll(i(jkr),jkr)% m      = p_ray% m     (jc_n,jray,jb_n)
  ray_coll(i(jkr),jkr)% dk     = p_ray% dk    (jc_n,jray,jb_n)
  ray_coll(i(jkr),jkr)% dl     = p_ray% dl    (jc_n,jray,jb_n)
  ray_coll(i(jkr),jkr)% dm     = p_ray% dm    (jc_n,jray,jb_n)
  ray_coll(i(jkr),jkr)% wadens = p_ray% wadens(jc_n,jray,jb_n)

  ray_coll(i(jkr),jkr)% jk_active = p_ray% jk_active(jc_n,jray,jb_n)

  ray_coll(i(jkr),jkr)% dlat = p_ray% dlat(jc_n,jray,jb_n)

  IF (c(2,iloc) == 0.) THEN   ! not split or split in lon
    ray_coll(i(jkr),jkr)% lat    = p_ray% lat   (jc_n,jray,jb_n)
    ray_coll(i(jkr),jkr)% coslat = p_ray% coslat(jc_n,jray,jb_n)
    ray_coll(i(jkr),jkr)% dlon   = p_ray% dlon  (jc_n,jray,jb_n)
  ELSE   ! split in lat
    ray_coll(i(jkr),jkr)% lat    =  p_ray% lat (jc_n,jray,jb_n)  &
      &                           + p_ray% dlat(jc_n,jray,jb_n)*c(2,iloc)
    ray_coll(i(jkr),jkr)% coslat =  coef_sphere*MAX(min_coslat, ABS(COS(  &
      &                                                         ray_coll(i(jkr),jkr)% lat)))  &
      &                           + coef_torus
    ray_coll(i(jkr),jkr)% dlon   =  p_ray% coslat(jc_n,jray,jb_n) * p_ray% dlon(jc_n,jray,jb_n)  &
      &                           / ray_coll(i(jkr),jkr)% coslat
  END IF
  IF (c(1,iloc) == 0.) THEN   ! not split or split in lat
    ray_coll(i(jkr),jkr)% lon = p_ray% lon(jc_n,jray,jb_n)
    !
    ray_coll(i(jkr),jkr)% k = p_ray% k(jc_n,jray,jb_n)
    ray_coll(i(jkr),jkr)% l = p_ray% l(jc_n,jray,jb_n)
  ELSE   ! split in lon
    ray_coll(i(jkr),jkr)% lon = p_ray% lon(jc_n,jray,jb_n) + p_spl% aux(jc_n,jray,jb_n)*c(1,iloc)
    ! rotation of (k,l) on spherical coordinates (no change when the torus grid is used)
    ray_coll(i(jkr),jkr)% k =  p_spl% rot_c(jc_n,jray,jb_n) * p_ray% k(jc_n,jray,jb_n)  &
      &            + c(1,iloc)*p_spl% rot_s(jc_n,jray,jb_n) * p_ray% l(jc_n,jray,jb_n)
    ray_coll(i(jkr),jkr)% l =  p_spl% rot_c(jc_n,jray,jb_n) * p_ray% l(jc_n,jray,jb_n)  &
      &            - c(1,iloc)*p_spl% rot_s(jc_n,jray,jb_n) * p_ray% k(jc_n,jray,jb_n)
  END IF

  ! Reset specid based on the new z, dz
  ray_coll(i(jkr),jkr)% specid = 0   ! In most cases, completely outside the ghost layer
  IF (p_ray% specid(jc_n,jray,jb_n) > 0) THEN
    IF ( ray_coll(i(jkr),jkr)% z - 0.5_wp*ray_coll(i(jkr),jkr)% dz  &
      &  < z_ifc(jc, ray_coll(i(jkr),jkr)% jk_active)               &
      &  .AND. jb == jb_n .AND. jc == jc_n )                        &
      &  ray_coll(i(jkr),jkr)% specid = p_ray% specid(jc_n,jray,jb_n)
  END IF

  ! |lat| <= pi/2   (for ray volumes moving or being split across the poles)
  !    CAUTION : The change of lon by pi is INCORRECT for the case of 4-segment splitting.
  !              Simply avoid the 4-seg. splitting across the pole when planning the splitting.
  !              Currently, only 2-seg. splitting is done in this case.

  IF (.NOT. is_plane_torus) THEN

    d = ABS(ray_coll(i(jkr),jkr)% lat) - pi_2
    IF (d > 0.) THEN
      ray_coll(i(jkr),jkr)% lat = SIGN(pi_2 - d, ray_coll(i(jkr),jkr)% lat)
      ray_coll(i(jkr),jkr)% lon = ray_coll(i(jkr),jkr)% lon - SIGN( pi, ray_coll(i(jkr),jkr)% lon )
    END IF

    ! |lon - lon_c| <= pi  [lon_c is the longitude of the (new) base grid-cell (jc,jb) center]

    d = p_patch% cells% center(jc,jb)% lon - ray_coll(i(jkr),jkr)% lon
    IF (ABS(d) > pi)  ray_coll(i(jkr),jkr)% lon = ray_coll(i(jkr),jkr)% lon  &
      &      + SIGN( two_pi*REAL((INT(ABS(d)/pi) + 1)/2, KIND=wp), d )

  ELSE  ! is_plane_torus

    ! TO BE CODED

  END IF

  ! Copy properties to the 2nd volume
  IF ( .NOT. l_1vol_z ) THEN
    i(jk_ul(2)) = i(jk_ul(2)) + 1

    ! Copy all properties
    ray_coll(i(jk_ul(2)),jk_ul(2)) = ray_coll(i(jk_ul(1)),jk_ul(1))

    ! Override z and its index
    ray_coll(i(jk_ul(2)),jk_ul(2))% z = z_ul (2)

    ! Reset specid based on the new z, dz
    IF (p_ray% specid(jc_n,jray,jb_n) > 0) THEN
      IF ( ray_coll(i(jk_ul(2)),jk_ul(2))% z - 0.5_wp*ray_coll(i(jk_ul(2)),jk_ul(2))% dz  &
        &  < z_ifc(jc, ray_coll(i(jk_ul(2)),jk_ul(2))% jk_active)                         &
        &  .AND. jb == jb_n .AND. jc == jc_n )                                            &
        &  ray_coll(i(jk_ul(2)),jk_ul(2))% specid = p_ray% specid(jc_n,jray,jb_n)
    END IF
  END IF

  NULLIFY( jkr )

END SUBROUTINE collect_vol

SUBROUTINE merge_vol( jiter, nz, nr_total_i, nr_reduce_total, clon, clat, fc2,  &
  &                   bvf2, gammash2, jc, nr, rc )

  INTEGER ,               INTENT(IN   ) :: jiter             ! iteration no. for the current call
  INTEGER ,               INTENT(IN   ) :: nz                ! nzray
  INTEGER ,               INTENT(IN   ) :: nr_total_i        ! SUM(nr_z)
  INTEGER ,               INTENT(IN   ) :: nr_reduce_total   ! nr_reduce_total
  REAL(wp),               INTENT(IN   ) :: clon              ! p_patch%cells%center(jc,jb)%lon
  REAL(wp),               INTENT(IN   ) :: clat              ! p_patch%cells%center(jc,jb)%lat
  REAL(wp),               INTENT(IN   ) :: fc2               ! f^2
  REAL(wp),               INTENT(IN   ) :: bvf2(:,:)         ! N^2 at half levels
  REAL(wp),               INTENT(IN   ) :: gammash2(:,:)     ! G^2 at half levels
  INTEGER ,               INTENT(IN   ) :: jc                ! nr_reduce_total
  INTEGER ,               INTENT(INOUT) :: nr(:)             ! nr_z(:)
  TYPE(t_ray_coll),       INTENT(INOUT) :: rc(:,:)           ! ray_coll(:,:)

  INTEGER , PARAMETER :: nmax_m  = 24   ! 24 or 32 (8*x): max/min = 4096 or 65536
  INTEGER , PARAMETER :: nmax_kh = 24   ! 24 or 32 (8*x): max/min = 4096 or 65536
  INTEGER , PARAMETER :: nmax_phi = 8
  INTEGER , PARAMETER :: min_nrays_z_def = 10
  REAL(wp), PARAMETER :: nr_reduce_factor(5) = (/1._wp, 1.25_wp, 1.5_wp, 2._wp, 99999._wp/)

  INTEGER  :: min_nrays_z
  INTEGER  :: nr_reduce(nz)
  INTEGER  :: nr_aim(nz)
  INTEGER  :: nphi
  INTEGER  :: nm
  INTEGER  :: nkh
  INTEGER  :: imerge
  INTEGER  :: jk, jr
  INTEGER  :: jr1
  INTEGER  :: j_m
  INTEGER  :: j_msgn
  INTEGER  :: j_kh
  INTEGER  :: j_phi
  INTEGER  :: jbin
  REAL(wp) :: inv_spmin_m2
  REAL(wp) :: inv_spmin_kh2
  REAL(wp) :: abs_l_o_k
  REAL(wp) :: dv6
  INTEGER  :: jr1_bin       (0:nmax_m-1, 2, 0:nmax_kh-1, 0:nmax_phi)
  LOGICAL  :: l_bin0_filled (0:nmax_m-1, 2, 0:nmax_kh-1, 0:nmax_phi)
  LOGICAL  :: l_binmg_filled(0:nmax_m-1, 2, 0:nmax_kh-1, 0:nmax_phi/2-1, 0:1)
  LOGICAL  :: l_exist(nrays_coll)
  INTEGER  :: nfill(0:1)
  INTEGER  :: jopt
  INTEGER  :: jopt_old
  INTEGER  :: j1, j2, j3, j4
  LOGICAL  :: l_mask(nz)
  REAL(sp) :: rnp1  (nz)
#ifdef __merge_energy
  REAL(wp) :: Kh2, m2pG2
#endif

  REAL(wp), PARAMETER :: spec_max_lz = 128._wp
  REAL(wp), PARAMETER :: spec_max_lh = 2048._wp
  REAL(wp), PARAMETER :: tanpi1_8 = TAN(pi*0.125_wp)
  REAL(wp), PARAMETER :: tanpi3_8 = TAN(pi*0.375_wp)

  inv_spmin_m2  = (spec_max_lz*1.e3_wp/(2._wp*pi))**2
  inv_spmin_kh2 = (spec_max_lh*1.e3_wp/(2._wp*pi))**2

  ! Minimum no. of ray volumes per layer to perform merging (set for efficiency)
  min_nrays_z = MIN( (nr_total_i - nr_reduce_total)/(2*nz), min_nrays_z_def )
    !                ^ (jray_end - j_ghost)/(2*nz) : half the layer-mean no. of ray volumes

  ! No. of ray volumes to be reduced in each layer
  nr_reduce(:) = INT( REAL(nr_reduce_total)/REAL(nr_total_i)*REAL(nr(:))  &
    &                * nr_reduce_factor(jiter) )
  nr_reduce(:) = MAX(0, MIN(nr(:) - min_nrays_z, nr_reduce(:)))
  nr_aim(:) = nr(:) - nr_reduce(:)
  l_mask(:) = nr_aim(:) > min_nrays_z
  DO jr = 1, nr_reduce_total - SUM(nr_reduce)
    rnp1(:) = REAL(nr_reduce(:)+1, sp)/REAL(MAX(1, nr(:)), sp)
    jk = MINLOC( rnp1, DIM=1, MASK=l_mask )
    nr_reduce(jk) = nr_reduce(jk) + 1
    nr_aim   (jk) = nr_aim   (jk) - 1
    l_mask   (jk) = nr_aim   (jk) > min_nrays_z
  ENDDO

  DO jk = 1, nz
    IF (nr_reduce(jk) <= 0)  CYCLE

    l_bin0_filled(:,:,:,:nmax_phi-1) = .FALSE.
    jr1_bin      (:,:,:,:nmax_phi-1) = -1

    DO jr = 1, nr(jk)

      ! Find the index of spectral bin for phi
      abs_l_o_k = ABS(rc(jr,jk)% l) / MAX(1.e-10_wp, ABS(rc(jr,jk)% k))
      IF (abs_l_o_k > tanpi3_8) THEN
        IF (rc(jr,jk)% l > 0.) THEN
          j_phi = 2
        ELSE
          j_phi = 6
        END IF
      ELSE IF (abs_l_o_k < tanpi1_8) THEN
        IF (rc(jr,jk)% k > 0.) THEN
          j_phi = 0
        ELSE
          j_phi = 4
        END IF
      ELSE IF (rc(jr,jk)% k > 0.) THEN
        IF (rc(jr,jk)% l > 0.) THEN
          j_phi = 1
        ELSE
          j_phi = 7
        END IF
      ELSE
        IF (rc(jr,jk)% l > 0.) THEN
          j_phi = 3
        ELSE
          j_phi = 5
        END IF
      END IF

      IF (rc(jr,jk)% m < 0.) THEN
        j_msgn = 1
      ELSE
        j_msgn = 2
      END IF

      ! Find the index of spectral bin for m
      j_m = MIN(nmax_m-1, MAX(0, EXPONENT( rc(jr,jk)% m**2 *inv_spmin_m2 )) )
        ! e.g., for spec_max_lz = 128.
        ! Lz > 128      :   0
        ! [128  , 64  ) :   1 or 2  [bin values increase by a factor of 2^(1/2)]
        ! [64   , 32  ) :   3 or 4

      ! Find the index of spectral bin for Kh
      j_kh = MIN(nmax_kh-1, MAX(0, EXPONENT( (rc(jr,jk)% k**2 + rc(jr,jk)% l**2)  &
        &                                        *inv_spmin_kh2 )) )
        ! e.g., for spec_max_lh = 2048.
        ! Lh > 2048     :   0
        ! [2048 , 1024) :   1 or 2  [bin values increase by a factor of 2^(1/2)]
        ! [1024 , 512 ) :   3 or 4

      IF ( .NOT. l_bin0_filled(j_m, j_msgn, j_kh, j_phi) ) THEN
        l_bin0_filled(j_m, j_msgn, j_kh, j_phi) = .TRUE.
        jr1_bin      (j_m, j_msgn, j_kh, j_phi) = jr
        rc(jr,jk)% l_modif = .FALSE.
      ELSE
        jr1 = jr1_bin(j_m, j_msgn, j_kh, j_phi)
        IF ( .NOT. rc(jr1,jk)% l_modif )  &
#ifdef __merge_energy
          &  CALL merge_task_1st( clon, fc2, bvf2(jc,jk), gammash2(jc,jk), rc(jr1,jk) )
        CALL merge_task_1st( clon, fc2, bvf2(jc,jk), gammash2(jc,jk), rc(jr,jk) )
#else
          &  CALL merge_task_1st( clon, rc(jr1,jk) )
        CALL merge_task_1st( clon, rc(jr,jk) )
#endif
        CALL merge_task_2nd( jr1, jr, rc(:nr(jk),jk) )
      END IF

    ENDDO  ! jr

    l_bin0_filled(:,:,:,nmax_phi) = l_bin0_filled(:,:,:,0)
    jr1_bin      (:,:,:,nmax_phi) = jr1_bin      (:,:,:,0)

    nm   = nmax_m
    nkh  = nmax_kh
    nphi = nmax_phi

    jopt = 0
    nfill(jopt) = COUNT(l_bin0_filled(:,:,:,:nmax_phi-1))

    IF (nfill(jopt) > nr_aim(jk)) THEN
      !
      ! merge the phi bins: 8 --> 4
      l_binmg_filled(:,:,:,:,0) =      l_bin0_filled(:,:,:,0:nphi-2:2)  &
        &                         .OR. l_bin0_filled(:,:,:,1:nphi-1:2)
      l_binmg_filled(:,:,:,:,1) =      l_bin0_filled(:,:,:,1:nphi-1:2)  &
        &                         .OR. l_bin0_filled(:,:,:,2:nphi  :2)
      nphi = nphi/2
      nfill(0) = COUNT(l_binmg_filled(:,:,:,:,0))
      nfill(1) = COUNT(l_binmg_filled(:,:,:,:,1))
      IF (nfill(0) <= nfill(1)) THEN   ! take jopt = 1
        jopt = 0
      ELSE
        jopt = 1
      END IF
      !
      DO j4 = 0, nphi-1
        DO j3 = 0, nkh-1
          DO j2 = 1, 2
            DO j1 = 0, nm-1
              jbin = j4*2 + jopt
              IF ( jr1_bin(j1,j2,j3,jbin) < 0 ) THEN
                jr1_bin(j1,j2,j3,j4) = jr1_bin(j1,j2,j3,jbin+1)
              ELSE IF (jr1_bin(j1,j2,j3,jbin+1) < 0) THEN
                jr1_bin(j1,j2,j3,j4) = jr1_bin(j1,j2,j3,jbin)
              ELSE   ! both are positive
                jr1_bin(j1,j2,j3,j4) = jr1_bin(j1,j2,j3,jbin)
                CALL merge_tasks( clon,                                              &
#ifdef __merge_energy
                      &           fc2, bvf2(jc,jk), gammash2(jc,jk),                 &
#endif
                                  jr1_bin(j1,j2,j3,jbin), jr1_bin(j1,j2,j3,jbin+1),  &
                  &               rc(:nr(jk),jk) )
              END IF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF (nfill(jopt) > nr_aim(jk)) THEN
        DO imerge = 1, 3   ! must be a sufficient number

          ! merge the m bins
          nm = nm/2
          jopt_old = jopt
          jopt = MOD(jopt+1, 2)
          l_binmg_filled(:nm-1,:,:nkh-1,:,jopt)                                    &
            &     =      l_binmg_filled( 0:nm*2-2:2 , : , :nkh-1 , : , jopt_old )  &
            &       .OR. l_binmg_filled( 1:nm*2-1:2 , : , :nkh-1 , : , jopt_old )
          nfill(jopt) = COUNT(l_binmg_filled(:nm-1,:,:nkh-1,:,jopt))
          DO j4 = 0, nphi-1
            DO j3 = 0, nkh-1
              DO j2 = 1, 2
                DO j1 = 0, nm-1
                  jbin = j1*2
                  IF ( jr1_bin(jbin,j2,j3,j4) < 0 ) THEN
                    jr1_bin(j1,j2,j3,j4) = jr1_bin(jbin+1,j2,j3,j4)
                  ELSE IF (jr1_bin(jbin+1,j2,j3,j4) < 0) THEN
                    jr1_bin(j1,j2,j3,j4) = jr1_bin(jbin,j2,j3,j4)
                  ELSE   ! both are positive
                    jr1_bin(j1,j2,j3,j4) = jr1_bin(jbin,j2,j3,j4)
                    CALL merge_tasks( clon,                                              &
#ifdef __merge_energy
                      &               fc2, bvf2(jc,jk), gammash2(jc,jk),                 &
#endif
                                      jr1_bin(jbin,j2,j3,j4), jr1_bin(jbin+1,j2,j3,j4),  &
                      &               rc(:nr(jk),jk) )
                  END IF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF (nfill(jopt) <= nr_aim(jk))  EXIT

          ! merge the kh bins
          nkh = nkh/2
          jopt_old = jopt
          jopt = MOD(jopt+1, 2)
          l_binmg_filled(:nm-1,:,:nkh-1,:,jopt)                                    &
            &     =      l_binmg_filled( :nm-1 , : , 0:nkh*2-2:2 , : , jopt_old )  &
            &       .OR. l_binmg_filled( :nm-1 , : , 1:nkh*2-1:2 , : , jopt_old )
          nfill(jopt) = COUNT(l_binmg_filled(:nm-1,:,:nkh-1,:,jopt))
          DO j4 = 0, nphi-1
            DO j3 = 0, nkh-1
              DO j2 = 1, 2
                DO j1 = 0, nm-1
                  jbin = j3*2
                  IF ( jr1_bin(j1,j2,jbin,j4) < 0 ) THEN
                    jr1_bin(j1,j2,j3,j4) = jr1_bin(j1,j2,jbin+1,j4)
                  ELSE IF (jr1_bin(j1,j2,jbin+1,j4) < 0) THEN
                    jr1_bin(j1,j2,j3,j4) = jr1_bin(j1,j2,jbin,j4)
                  ELSE   ! both are positive
                    jr1_bin(j1,j2,j3,j4) = jr1_bin(j1,j2,jbin,j4)
                    CALL merge_tasks( clon,                                              &
#ifdef __merge_energy
                      &               fc2, bvf2(jc,jk), gammash2(jc,jk),                 &
#endif
                      &               jr1_bin(j1,j2,jbin,j4), jr1_bin(j1,j2,jbin+1,j4),  &
                      &               rc(:nr(jk),jk) )
                  END IF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF (nfill(jopt) <= nr_aim(jk))  EXIT

        ENDDO
      END IF
      !
    END IF

    l_exist(:nr(jk)) = .FALSE.

    DO j4 = 0, nphi-1   ! independent loop: jr's are different for each (j1,j2,j3,j4)
      DO j3 = 0, nkh-1
        DO j2 = 1, 2
          DO j1 = 0, nm-1
            IF ( jr1_bin(j1,j2,j3,j4) < 0 )  CYCLE
            l_exist(jr1_bin(j1,j2,j3,j4)) = .TRUE.
            IF ( .NOT. rc(jr1_bin(j1,j2,j3,j4),jk)% l_modif )  CYCLE
            !
            jr = jr1_bin(j1,j2,j3,j4)
            !
            rc(jr,jk)% dlat = rc(jr,jk)% dlat - rc(jr,jk)% lat_min
            rc(jr,jk)% dz   = rc(jr,jk)% dz   - rc(jr,jk)% z_min
#ifdef __mgm_minmax_center
            rc(jr,jk)% lat  = rc(jr,jk)% lat_min + rc(jr,jk)% dlat*0.5_wp
            rc(jr,jk)% z    = rc(jr,jk)% z_min   + rc(jr,jk)% dz  *0.5_wp
#else
            rc(jr,jk)% lat  = rc(jr,jk)% lat / rc(jr,jk)% wadens
            rc(jr,jk)% z    = rc(jr,jk)% z   / rc(jr,jk)% wadens
!org        rc(jr,jk)% z    = rc(jr,jk)% z_min   + rc(jr,jk)% dz  *0.5_wp
#endif
            !
!           rc(jr,jk)% coslat = ABS(COS(rc(jr,jk)% lat))
            rc(jr,jk)% coslat = coef_sphere*MAX(min_coslat, ABS(COS(rc(jr,jk)% lat))) + coef_torus
            !
            rc(jr,jk)% dlon = (rc(jr,jk)% dlon - rc(jr,jk)% lon_min)/rc(jr,jk)% coslat
#ifdef __mgm_minmax_center
            rc(jr,jk)% lon  = clon + rc(jr,jk)% lon_min/rc(jr,jk)% coslat + rc(jr,jk)% dlon*0.5_wp
#else
            rc(jr,jk)% lon  = clon + rc(jr,jk)% lon / (rc(jr,jk)% wadens * rc(jr,jk)% coslat)
#endif
            !
#ifdef __mgm_minmax_center
            rc(jr,jk)% dk = rc(jr,jk)% dk - rc(jr,jk)% k
            rc(jr,jk)% dl = rc(jr,jk)% dl - rc(jr,jk)% l
            rc(jr,jk)% dm = rc(jr,jk)% dm - rc(jr,jk)% m
            rc(jr,jk)% k  = rc(jr,jk)% k  + rc(jr,jk)% dk*0.5_wp
            rc(jr,jk)% l  = rc(jr,jk)% l  + rc(jr,jk)% dl*0.5_wp
            rc(jr,jk)% m  = rc(jr,jk)% m  + rc(jr,jk)% dm*0.5_wp
            !
#ifdef __merge_energy
            Kh2   = rc(jr,jk)% k**2 + rc(jr,jk)% l**2
            m2pG2 = rc(jr,jk)% m**2 + gammash2(jc,jk)
            rc(jr,jk)% wadens = rc(jr,jk)% wadens  &
              &                 * SQRT(( Kh2 + m2pG2 )/( bvf2(jc,jk)*Kh2 + fc2*m2pG2 ))
#endif
            dv6 = (rc(jr,jk)% coslat * rc(jr,jk)% dlon) * rc(jr,jk)% dlat * rc(jr,jk)% dz  &
              &                      * rc(jr,jk)% dk    * rc(jr,jk)% dl   * rc(jr,jk)% dm
            IF (dv6 /= 0.) THEN
              rc(jr,jk)% wadens = rc(jr,jk)% wadens / dv6
            ELSE   ! rarely
              ! dv could have been truncated to zero if it is so small
              ! e.g., if dk <<< |k|:  k + dk/2 = k - dk/2 = k  --> dk is set to 0
              dv6 = (rc(jr,jk)% coslat * rc(jr,jk)% dlon) * rc(jr,jk)% dlat * rc(jr,jk)% dz
              IF (dv6 == 0.) THEN   ! dlon, dlat, or dz ~= 0 :  ignore
                l_exist(jr) = .FALSE.
              ELSE
                rc(jr,jk)% dk = MAX(1.e-9_wp, rc(jr,jk)% dk)  ! < (2pi/4000km) * 0.1%
                rc(jr,jk)% dl = MAX(1.e-9_wp, rc(jr,jk)% dl)
                rc(jr,jk)% dm = MAX(1.e-9_wp, rc(jr,jk)% dm)
                dv6 = dv6 * rc(jr,jk)% dk * rc(jr,jk)% dl * rc(jr,jk)% dm
                rc(jr,jk)% wadens = rc(jr,jk)% wadens / dv6
              END IF
            END IF
            IF (rc(jr,jk)% wadens <= eps_wadens)  l_exist(jr) = .FALSE.
#else
            ! dk, dl, dm are normalized,
            ! and accordingly wadens includes the factor (dk dl dm) : A
            rc(jr,jk)% dk = 1.
            rc(jr,jk)% dl = 1.
            rc(jr,jk)% dm = 1.
            rc(jr,jk)% k  = rc(jr,jk)% k / rc(jr,jk)% wadens
            rc(jr,jk)% l  = rc(jr,jk)% l / rc(jr,jk)% wadens
            rc(jr,jk)% m  = rc(jr,jk)% m / rc(jr,jk)% wadens
            !
            dv6 = (rc(jr,jk)% coslat * rc(jr,jk)% dlon) * rc(jr,jk)% dlat * rc(jr,jk)% dz
            IF (dv6 /= 0.) THEN
#ifdef __merge_energy
              Kh2   = rc(jr,jk)% k**2 + rc(jr,jk)% l**2
              m2pG2 = rc(jr,jk)% m**2 + gammash2(jc,jk)
              rc(jr,jk)% wadens = rc(jr,jk)% wadens / dv6  &
                &                 * SQRT(( Kh2 + m2pG2 )/( bvf2(jc,jk)*Kh2 + fc2*m2pG2 ))
#else
              rc(jr,jk)% wadens = rc(jr,jk)% wadens / dv6
#endif
              IF (rc(jr,jk)% wadens <= eps_wadens)  l_exist(jr) = .FALSE.
            ELSE   ! rarely
              l_exist(jr) = .FALSE.
            END IF
#endif
!           rc(jr,jk)% l_modif = .FALSE.
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    nfill(jopt) = COUNT(l_exist(:nr(jk)))   ! update

    rc(:nfill(jopt),jk) = PACK(rc(:nr(jk),jk), l_exist(:nr(jk)))

    nr(jk) = nfill(jopt)

  ENDDO  ! jk

END SUBROUTINE merge_vol

#ifdef __merge_energy
SUBROUTINE merge_task_1st( clon, fc2, bvf2, gammash2, rc )
#else
SUBROUTINE merge_task_1st( clon, rc )
#endif

  REAL(wp),               INTENT(IN   ) :: clon
#ifdef __merge_energy
  REAL(wp),               INTENT(IN   ) :: fc2, bvf2, gammash2
#endif
  TYPE(t_ray_coll),       INTENT(INOUT) :: rc

#ifdef __merge_energy
  REAL(wp) :: Kh2, m2pG2

  Kh2   = rc% k**2 + rc% l**2
  m2pG2 = rc% m**2 + gammash2
  rc% wadens = rc% wadens * SQRT(( bvf2*Kh2 + fc2*m2pG2 )/( Kh2 + m2pG2 ))  &
#else
  rc% wadens = rc% wadens  &
#endif
    &  *( (rc% coslat * rc% dlon) * rc% dlat * rc% dz  &
    &                 * rc% dk    * rc% dl   * rc% dm )

  rc% lon = rc% lon - clon
  !
  rc% lon_min = rc% lon - rc% dlon*0.5_wp          ! lon_0 - lon_grid
  rc% lat_min = rc% lat - rc% dlat*0.5_wp          ! lat_0
  rc% z_min   = rc% z   - rc% dz  *0.5_wp          !   z_0
  rc% dlon    = rc% lon_min + rc% dlon             ! lon_1 - lon_grid
  rc% dlat    = rc% lat_min + rc% dlat             ! lat_1
  rc% dz      = rc% z_min   + rc% dz               !   z_1
  !
  rc% lon_min = rc% lon_min * rc% coslat   ! latitude dependency
  rc% dlon    = rc% dlon    * rc% coslat
  !
  rc% lon = rc% lon * rc% wadens * rc% coslat
  rc% lat = rc% lat * rc% wadens
  rc% z   = rc% z   * rc% wadens
  !
#ifdef __mgm_minmax_center
  rc% k  = rc% k - rc% dk*0.5_wp     !   k_0
  rc% l  = rc% l - rc% dl*0.5_wp     !   l_0
  rc% m  = rc% m - rc% dm*0.5_wp     !   m_0
  rc% dk = rc% k + rc% dk            !   k_1
  rc% dl = rc% l + rc% dl            !   l_1
  rc% dm = rc% m + rc% dm            !   m_1
#else
  rc% k  = rc% k * rc% wadens
  rc% l  = rc% l * rc% wadens
  rc% m  = rc% m * rc% wadens
  ! not calculate dk, dl, dm
#endif

  rc% l_modif = .TRUE.

END SUBROUTINE merge_task_1st

SUBROUTINE merge_task_2nd( jr, jr_i, rc )

  INTEGER ,               INTENT(IN   ) :: jr
  INTEGER ,               INTENT(IN   ) :: jr_i
  TYPE(t_ray_coll),       INTENT(INOUT) :: rc(:)

  rc(jr)% wadens = rc(jr)% wadens + rc(jr_i)% wadens
  rc(jr)% lon    = rc(jr)% lon    + rc(jr_i)% lon
  rc(jr)% lat    = rc(jr)% lat    + rc(jr_i)% lat
  rc(jr)% z      = rc(jr)% z      + rc(jr_i)% z
  !
  rc(jr)% lon_min = MIN(rc(jr)% lon_min, rc(jr_i)% lon_min)
  rc(jr)% lat_min = MIN(rc(jr)% lat_min, rc(jr_i)% lat_min)
  rc(jr)% z_min   = MIN(rc(jr)% z_min  , rc(jr_i)% z_min  )
  rc(jr)% dlon = MAX(rc(jr)% dlon, rc(jr_i)% dlon)
  rc(jr)% dlat = MAX(rc(jr)% dlat, rc(jr_i)% dlat)
  rc(jr)% dz   = MAX(rc(jr)% dz  , rc(jr_i)% dz  )
  !
#ifdef __mgm_minmax_center
  rc(jr)% k  = MIN(rc(jr)% k , rc(jr_i)% k )
  rc(jr)% l  = MIN(rc(jr)% l , rc(jr_i)% l )
  rc(jr)% m  = MIN(rc(jr)% m , rc(jr_i)% m )
  rc(jr)% dk = MAX(rc(jr)% dk, rc(jr_i)% dk)
  rc(jr)% dl = MAX(rc(jr)% dl, rc(jr_i)% dl)
  rc(jr)% dm = MAX(rc(jr)% dm, rc(jr_i)% dm)
#else
  rc(jr)% k  = rc(jr)% k + rc(jr_i)% k
  rc(jr)% l  = rc(jr)% l + rc(jr_i)% l
  rc(jr)% m  = rc(jr)% m + rc(jr_i)% m
#endif

!  rc(jr)% specid = specid_merged   ! may not be used
  rc(jr)% specid = 0

!     jk_active : not meaningful now

END SUBROUTINE merge_task_2nd

#ifdef __merge_energy
SUBROUTINE merge_tasks( clon, fc2, bvf2, gammash2, jr, jr_i, rc )
#else
SUBROUTINE merge_tasks( clon, jr, jr_i, rc )
#endif

  REAL(wp),               INTENT(IN   ) :: clon
#ifdef __merge_energy
  REAL(wp),               INTENT(IN   ) :: fc2, bvf2, gammash2
#endif
  INTEGER ,               INTENT(IN   ) :: jr
  INTEGER ,               INTENT(INOUT) :: jr_i
  TYPE(t_ray_coll),       INTENT(INOUT) :: rc(:)

#ifdef __merge_energy
  IF ( .NOT. rc(jr_i)% l_modif )  CALL merge_task_1st( clon, fc2, bvf2, gammash2, rc(jr_i) )
  IF ( .NOT. rc(jr  )% l_modif )  CALL merge_task_1st( clon, fc2, bvf2, gammash2, rc(jr  ) )
#else
  IF ( .NOT. rc(jr_i)% l_modif )  CALL merge_task_1st( clon, rc(jr_i) )
  IF ( .NOT. rc(jr  )% l_modif )  CALL merge_task_1st( clon, rc(jr  ) )
#endif

  CALL merge_task_2nd( jr, jr_i, rc(:) )

  jr_i = -1   ! needed ?

END SUBROUTINE merge_tasks

SUBROUTINE remove_vol( nz, nr_total_i, nr_reduce_total, fc2, bvf2, gammash2, jc, nr, rc )

  INTEGER ,               INTENT(IN   ) :: nz                ! nzray
  INTEGER ,               INTENT(IN   ) :: nr_total_i        ! SUM(nr_z)
  INTEGER ,               INTENT(IN   ) :: nr_reduce_total   ! nr_reduce_total
  REAL(wp),               INTENT(IN   ) :: fc2               ! f^2
  REAL(wp),               INTENT(IN   ) :: bvf2(:,:)         ! N^2 at half levels
  REAL(wp),               INTENT(IN   ) :: gammash2(:,:)     ! G^2 at half levels
  INTEGER ,               INTENT(IN   ) :: jc                ! nr_reduce_total
  INTEGER ,               INTENT(INOUT) :: nr(:)             ! nr_z(:)
  TYPE(t_ray_coll),       INTENT(INOUT) :: rc(:,:)           ! ray_coll(:,:)

  INTEGER , PARAMETER :: min_nrays_z_def = 10   ! consistent with subroutine 'remove_rays'

  INTEGER  :: min_nrays_z
  INTEGER  :: nr_reduce(nz)
  INTEGER  :: nr_aim   (nz)
  LOGICAL  :: l_mask   (nz)
  REAL(sp) :: rnp1     (nz)
  REAL(wp) :: quantity_vol(nrays_coll)
  REAL(wp) :: Kh2, m2pG2, omega2
  INTEGER  :: jr, jk, jn

  ! Minimum no. of ray volumes per layer to perform merging (set for efficiency)
  min_nrays_z = MIN( (nr_total_i - nr_reduce_total)/(2*nz), min_nrays_z_def )
    !                ^ (jray_end - j_ghost)/(2*nz) : half the layer-mean no. of ray volumes

  ! No. of ray volumes to be reduced in each layer
  nr_reduce(:) = INT( REAL(nr_reduce_total)/REAL(nr_total_i)*REAL(nr(:)) )
  nr_reduce(:) = MAX(0, MIN(nr(:) - min_nrays_z, nr_reduce(:)))
  nr_aim(:) = nr(:) - nr_reduce(:)
  l_mask(:) = nr_aim(:) > min_nrays_z
  DO jr = 1, nr_reduce_total - SUM(nr_reduce)
    rnp1(:) = REAL(nr_reduce(:)+1, sp)/REAL(MAX(1, nr(:)), sp)
    jk = MINLOC( rnp1, DIM=1, MASK=l_mask )
    nr_reduce(jk) = nr_reduce(jk) + 1
    nr_aim   (jk) = nr_aim   (jk) - 1
    l_mask   (jk) = nr_aim   (jk) > min_nrays_z
  ENDDO

  ! Remove ray volumes with the lowest wave energy
  DO jk = 1, nz
    IF (nr_reduce(jk) <= 0)  CYCLE

    DO jr = 1, nr(jk)

      ! Horizontal wavenumber squared ; Vertical wavenumber squared + gamma^2
      Kh2   = rc(jr,jk)% k**2 + rc(jr,jk)% l**2
      m2pG2 = rc(jr,jk)% m**2 + gammash2(jc,jk)

      ! Intrinsic frequency squared
      omega2 = ( bvf2(jc,jk)*Kh2 + fc2*m2pG2 )/( Kh2 + m2pG2 )

      ! Energy fraction of ray volume
      quantity_vol(jr) = ABS(                                                        &
        &   (rc(jr,jk)% coslat * rc(jr,jk)% dlon) * rc(jr,jk)% dlat * rc(jr,jk)% dz  &
        &                      * rc(jr,jk)% dk    * rc(jr,jk)% dl   * rc(jr,jk)% dm  &
        &                      * rc(jr,jk)% wadens * SQRT(omega2) )

    ENDDO  ! jr

    DO jn = 1, nr_reduce(jk)
      jr = MINLOC( quantity_vol(:nr(jk)), DIM=1, MASK=( rc(:nr(jk),jk)% wadens /= 0._wp ) )
      rc(jr,jk)% wadens = 0._wp
    ENDDO  ! jn

    rc(:nr_aim(jk),jk) = PACK(rc(:nr(jk),jk), rc(:nr(jk),jk)% wadens /= 0._wp)

  ENDDO  ! jk

  nr(:) = nr_aim(:)   ! nr_aim(:) == nr(:) - nr_reduce(:)

END SUBROUTINE remove_vol
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE sync_wave(p_ray,p_patch)
  TYPE(t_ray),            INTENT(INOUT) :: p_ray                ! the ray properties
  TYPE(t_patch),   TARGET,INTENT(INOUT) :: p_patch              ! grid/patch info.

  IF (msg_level >= 12) CALL message('sync_wave', 'MS-GWaM: synchronize ray volumes')

  !----------------------------------------------------------------------
  ! Purpose:
  !         Synchronize derived type elements of ray volume properties
  !
  ! Method:
  !         -- based on sync_patch_array 
  !            (src/parallel_infrastructure/mo_sync.f90)
  !
  !----------------------------------------------------------------------

  !CALL message('', TRIM(message_text))
  !      WRITE(message_text,'(a,i6)') 'shape p_ray%iexist:', SHAPE(p_ray%iexist)
  !print*, 'SHAPE p_ray%iexist:', SHAPE(p_ray%iexist)

  ! Synchronize integer properties
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%iexist)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%specid)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%jk_active)
  ! For the 6 properties below no synchronization is needed
  ! because they matter either only in the source treatment or
  ! they are anyway recalculated in the new position by idx_rayedge
  ! p_ray_bg%jk_source
  ! p_ray_conv%jk_source
  ! p_ray%jr_last
  ! p_ray%jk_full_rtop
  ! p_ray%jk_full_rbot
  ! p_ray%jk_half_rtop
  ! p_ray%jk_half_rbot

  ! Synchronize real properties
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%lon)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%lat)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%z)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%dlon)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%dlat)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%dz)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%k)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%l)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%m)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%dk)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%dl)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%dm)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%wadens)

  ! Later, directly calculate the cosine here rather than perform MPI-comm. (if it is faster)
  CALL sync_patch_array(SYNC_C, p_patch, p_ray%coslat)

END SUBROUTINE sync_wave
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE regrid_wave( p_patch,        & !in(out)
  &                     p_gridinfo4ray, & !in
  &                     z_ifc,          & !in
  &                     z_mc,           & !in
  &                     nrays,          & !in
  &                     p_ray           ) !inout

  ! In/out variables
  TYPE(t_patch),        INTENT(INOUT) :: p_patch
  TYPE(t_gridinfo4ray), INTENT(IN)    :: p_gridinfo4ray
  REAL(wp),             INTENT(IN)    :: z_ifc(:,:,:)
  REAL(wp),             INTENT(IN)    :: z_mc (:,:,:)
  INTEGER,              INTENT(IN)    :: nrays
  TYPE(t_ray),          INTENT(INOUT) :: p_ray

  ! Local variables
  TYPE(t_point)                       :: v1,v2,v3      !< vertex longitudes/latitudes

  INTEGER  :: rl_start, rl_end
  INTEGER  :: i_startblk, i_endblk
  INTEGER  :: i_startidx, i_endidx
  INTEGER  :: jb, jc, jstencil
  INTEGER  :: jr, jray_loc, jray_nb, jray_loc_old
  INTEGER  :: icn, ibn
  INTEGER  :: jk0, jk0m1

  REAL(wp) :: ray_xyz(3)

  INTEGER  :: intriangle

  !----------------------------------------------------------------------
  ! Purpose:
  !         
  !
  ! Method:
  !         -- 
  !
  !----------------------------------------------------------------------

  ! Prognostic domain
  rl_start   = grf_bdywidth_c + 1
  rl_end     = min_rlcell_int
  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

  IF (.NOT. is_plane_torus) THEN

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb, jc, jray_loc, jray_nb, i_startidx, i_endidx, jstencil, icn, ibn) ICON_OMP_GUIDED_SCHEDULE

    ! Find ray volumes propagated to this cell from neighboring 
    ! cells and initialize them in the local cell
    DO jb = i_startblk, i_endblk ! jb index of the local cell

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx ! jc index of the local cell

        ! Vertex 1 in Cartesian coordinates
        v1%x = p_gridinfo4ray%cellvertices_x(1,jc,jb)
        v1%y = p_gridinfo4ray%cellvertices_y(1,jc,jb)
        v1%z = p_gridinfo4ray%cellvertices_z(1,jc,jb)

        ! Vertex 2 in Cartesian coordinates
        v2%x = p_gridinfo4ray%cellvertices_x(2,jc,jb)
        v2%y = p_gridinfo4ray%cellvertices_y(2,jc,jb)
        v2%z = p_gridinfo4ray%cellvertices_z(2,jc,jb)

        ! Vertex 3 in Cartesian coordinates
        v3%x = p_gridinfo4ray%cellvertices_x(3,jc,jb)
        v3%y = p_gridinfo4ray%cellvertices_y(3,jc,jb)
        v3%z = p_gridinfo4ray%cellvertices_z(3,jc,jb)

        ! Initialize local ray volume index
        jray_loc = 0

        ! Loop over neighboring cell indices
        DO jstencil = 1, p_gridinfo4ray%cellneighbors_nstencil(jc,jb)  ! 11 or 12

          icn = p_gridinfo4ray%cellneighbors_idx(jstencil,jc,jb) ! jc index of neighboring cell
          ibn = p_gridinfo4ray%cellneighbors_blk(jstencil,jc,jb) ! jb index of neighboring cell

          ! Loop over ray volumes in neighboring cell
          DO jray_nb = 1, nrays

            ! Skip calculations if this ray volume in the neighboring cell
            ! does not exist or is below the launch level z_src
            IF (p_ray%iexist(icn,jray_nb,ibn) == 0 .OR. &
                p_ray%specid(icn,jray_nb,ibn) < 0) CYCLE

            ! Ray volume center point in Cartesian coordinates
            ray_xyz(1) = p_ray%coslat(icn,jray_nb,ibn)*COS(p_ray%lon(icn,jray_nb,ibn))
            ray_xyz(2) = p_ray%coslat(icn,jray_nb,ibn)*SIN(p_ray%lon(icn,jray_nb,ibn))
            ray_xyz(3) = SIGN(SQRT(1. - p_ray%coslat(icn,jray_nb,ibn)**2), p_ray%lat(icn,jray_nb,ibn))

            intriangle = inside_triangle(ray_xyz, v1,v2,v3)

            IF (intriangle >= 0) THEN

              jk0 = p_ray%iexist(icn,jray_nb,ibn)
              jk0m1 = MAX(1, jk0 - 1)

              IF ( z_mc(icn,jk0,ibn) /= z_mc(jc,jk0,jb) .OR. z_mc(icn,jk0m1,ibn) /= z_mc(jc,jk0m1,jb) ) THEN
                !
                ! Some conditions used at the end of the integration (in propagate_wave) need
                ! to be applied again here, as the vertical grids are changed.
                !  - lower bound used in propagate_wave
!               IF (p_ray%z(icn,jray_nb,ibn) < z_ifc(jc, p_patch% nlevp1, jb)) THEN
                IF (p_ray%z(icn,jray_nb,ibn) < z_ifc(jc, p_patch% nlevp1 - nlev_lbnd(jg), jb)) THEN
                  ! outofdomain
!                 p_ray%iexist(icn,jray_nb,ibn) = 0      ! CHECK: problematic for OpenMP ?
!                 p_ray%specid(icn,jray_nb,ibn) = 0      !
                  p_ray%wadens(icn,jray_nb,ibn) = 0.     ! - This must be okay
                  CYCLE
                END IF
                !
                ! update jk0
                CALL get_jkhalf_closest_gen( jk0, p_ray% z(icn,jray_nb,ibn),  &
                  &                          p_patch% nlev, z_mc(:,:,jb), jc )
              END IF

              ! Find the next non-existing local ray volume index, which can be used
              ! for hosting the new ray volume moved here from neighboring cell
              jray_loc_old = jray_loc
              DO jr = jray_loc_old+1, nrays
                IF (p_ray%iexist(jc,jr,jb) == 0) THEN
                  jray_loc = jr
                  EXIT
                END IF
              ENDDO
              IF (jray_loc == jray_loc_old)  CALL finish ('regrid_wave', &
                  'Something is strange in horizontal regridding: jray_loc>nrays')

              ! Initialize new ray volume
              ! Synchronized properties: take them from the ray volume that moved 
              ! here from neighboring cell
              p_ray%iexist(jc,jray_loc,jb)    = jk0
              p_ray%specid(jc,jray_loc,jb)    = p_ray%specid(icn,jray_nb,ibn)
              p_ray%jk_active(jc,jray_loc,jb) = p_ray%jk_active(icn,jray_nb,ibn)
              p_ray%lon(jc,jray_loc,jb)       = p_ray%lon(icn,jray_nb,ibn)
              p_ray%lat(jc,jray_loc,jb)       = p_ray%lat(icn,jray_nb,ibn)
              p_ray%z(jc,jray_loc,jb)         = p_ray%z(icn,jray_nb,ibn)
              p_ray%dlon(jc,jray_loc,jb)      = p_ray%dlon(icn,jray_nb,ibn)
              p_ray%dlat(jc,jray_loc,jb)      = p_ray%dlat(icn,jray_nb,ibn)
              p_ray%dz(jc,jray_loc,jb)        = p_ray%dz(icn,jray_nb,ibn)
              p_ray%coslat(jc,jray_loc,jb)    = p_ray%coslat(icn,jray_nb,ibn)
              p_ray%k(jc,jray_loc,jb)         = p_ray%k(icn,jray_nb,ibn)
              p_ray%l(jc,jray_loc,jb)         = p_ray%l(icn,jray_nb,ibn)
              p_ray%m(jc,jray_loc,jb)         = p_ray%m(icn,jray_nb,ibn)
              p_ray%dk(jc,jray_loc,jb)        = p_ray%dk(icn,jray_nb,ibn)
              p_ray%dl(jc,jray_loc,jb)        = p_ray%dl(icn,jray_nb,ibn)
              p_ray%dm(jc,jray_loc,jb)        = p_ray%dm(icn,jray_nb,ibn)
              p_ray%wadens(jc,jray_loc,jb)    = p_ray%wadens(icn,jray_nb,ibn)

              ! Remove from the neighboring cell
!             p_ray%iexist(icn,jray_nb,ibn) = 0      ! CHECK: problematic for OpenMP ?
!             p_ray%specid(icn,jray_nb,ibn) = 0      !
              p_ray%wadens(icn,jray_nb,ibn) = 0.     ! - This must be okay

              ! Horizontal propagation test
              IF (ltest_hprop .AND. l1ray) THEN
                ! Open output file for horizontal propagation test
                OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
                WRITE(987,'(a,3i4,2x,6F14.5)') &
                       'Regrid ray taken from lon, lat:', jstencil, icn, ibn, &
                                          p_gridinfo4ray%cellvertices_lon(1,icn,ibn)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(1,icn,ibn)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lon(2,icn,ibn)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(2,icn,ibn)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lon(3,icn,ibn)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(3,icn,ibn)*rad2deg
                WRITE(987,'(a,2i4,2x,6F14.5)') &
                       'Regrid ray taken   to lon, lat:', jc, jb, &
                                          p_gridinfo4ray%cellvertices_lon(1,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(1,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lon(2,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(2,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lon(3,jc,jb)*rad2deg,&
                                          p_gridinfo4ray%cellvertices_lat(3,jc,jb)*rad2deg
                CLOSE(987)
              ENDIF ! ltest_hprop .AND. l1ray

            ENDIF ! intriangle

          ENDDO  ! jray_nb

        ENDDO  ! jstencil

      ENDDO  ! jc

    ENDDO  ! jb

!!$OMP END DO
!!$OMP END PARALLEL

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb, jc, jray_loc, i_startidx, i_endidx, jstencil, icn, ibn) ICON_OMP_GUIDED_SCHEDULE

    ! Find ray volumes that propagated out of this 
    ! cell and remove them
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx

        ! Vertex 1 in Cartesian coordinates
        v1%x = p_gridinfo4ray%cellvertices_x(1,jc,jb)
        v1%y = p_gridinfo4ray%cellvertices_y(1,jc,jb)
        v1%z = p_gridinfo4ray%cellvertices_z(1,jc,jb)

        ! Vertex 2 in Cartesian coordinates
        v2%x = p_gridinfo4ray%cellvertices_x(2,jc,jb)
        v2%y = p_gridinfo4ray%cellvertices_y(2,jc,jb)
        v2%z = p_gridinfo4ray%cellvertices_z(2,jc,jb)

        ! Vertex 3 in Cartesian coordinates
        v3%x = p_gridinfo4ray%cellvertices_x(3,jc,jb)
        v3%y = p_gridinfo4ray%cellvertices_y(3,jc,jb)
        v3%z = p_gridinfo4ray%cellvertices_z(3,jc,jb)

        DO jray_loc = 1, nrays

          IF (p_ray%iexist(jc,jray_loc,jb) == 0 .OR. &
              p_ray%specid(jc,jray_loc,jb) < 0) CYCLE

          IF (p_ray%wadens(jc,jray_loc,jb) == 0._wp) THEN
            ! Ray volume that has escaped to another prognostic cell
            p_ray%iexist(jc,jray_loc,jb) = 0
            p_ray%specid(jc,jray_loc,jb) = 0
            CYCLE
!         ELSE IF ( this_cell_is_not_a_neighbor_of_halo ) THEN
!           ! Ray volume inside this cell, far from halo: skip the in-triangle test
!           CYCLE
          END IF

          ! Remainder: ray volume either being inside this cell or escaped to a halo cell

          ! Find ray volumes that propagated out of this cell and 
          ! remove them
          ! Ray volume center point in Cartesian coordinates
          ray_xyz(1) = p_ray%coslat(jc,jray_loc,jb)*COS(p_ray%lon(jc,jray_loc,jb))
          ray_xyz(2) = p_ray%coslat(jc,jray_loc,jb)*SIN(p_ray%lon(jc,jray_loc,jb))
          ray_xyz(3) = SIGN(SQRT(1. - p_ray%coslat(jc,jray_loc,jb)**2), p_ray%lat(jc,jray_loc,jb))

          intriangle = inside_triangle(ray_xyz, v1,v2,v3)

          IF (intriangle < 0) THEN  ! not in triangle

            ! Remove ray volume
            p_ray%iexist(jc,jray_loc,jb) = 0
            p_ray%specid(jc,jray_loc,jb) = 0
            p_ray%wadens(jc,jray_loc,jb) = 0._wp

            ! Horizontal propagation test
            IF (ltest_hprop .AND. l1ray) THEN
              ! Open output file for horizontal propagation test
              OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
              WRITE(987,'(a,2i4,2x,6F14.5)') &
                     'Regrid ray removed at lon, lat:', jc, jb, &
                                        p_gridinfo4ray%cellvertices_lon(1,jc,jb)*rad2deg,&
                                        p_gridinfo4ray%cellvertices_lat(1,jc,jb)*rad2deg,&
                                        p_gridinfo4ray%cellvertices_lon(2,jc,jb)*rad2deg,&
                                        p_gridinfo4ray%cellvertices_lat(2,jc,jb)*rad2deg,&
                                        p_gridinfo4ray%cellvertices_lon(3,jc,jb)*rad2deg,&
                                        p_gridinfo4ray%cellvertices_lat(3,jc,jb)*rad2deg
              CLOSE(987)
            ENDIF ! ltest_hprop .AND l1ray

          ENDIF ! not intriangle

        ENDDO  ! jray_loc

      ENDDO  ! jc

    ENDDO  ! jb

!!$OMP END DO
!!$OMP END PARALLEL

  ELSE ! is_plane_torus

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb, jc, jray_loc, jray_nb, i_startidx, i_endidx, jstencil, icn, ibn) ICON_OMP_GUIDED_SCHEDULE

    ! Find ray volumes propagated to this cell from neighboring 
    ! cells and initialize them in the local cell
    DO jb = i_startblk, i_endblk ! jb index of the local cell

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx ! jc index of the local cell

        ! Vertex 1 in Cartesian coordinates
        v1%x = p_gridinfo4ray%cellvertices_x(1,jc,jb)
        v1%y = p_gridinfo4ray%cellvertices_y(1,jc,jb)
        v1%z = p_gridinfo4ray%cellvertices_z(1,jc,jb) ! not used on the plane

        ! Vertex 2 in Cartesian coordinates
        v2%x = p_gridinfo4ray%cellvertices_x(2,jc,jb)
        v2%y = p_gridinfo4ray%cellvertices_y(2,jc,jb)
        v2%z = p_gridinfo4ray%cellvertices_z(2,jc,jb) ! not used on the plane

        ! Vertex 3 in Cartesian coordinates
        v3%x = p_gridinfo4ray%cellvertices_x(3,jc,jb)
        v3%y = p_gridinfo4ray%cellvertices_y(3,jc,jb)
        v3%z = p_gridinfo4ray%cellvertices_z(3,jc,jb) ! not used on the plane

        ! Initialize local ray volume index
        jray_loc = 0

        ! Loop over neighboring cell indices
        DO jstencil = 1, p_gridinfo4ray%cellneighbors_nstencil(jc,jb)  ! 11 or 12

          icn = p_gridinfo4ray%cellneighbors_idx(jstencil,jc,jb) ! jc index of neighboring cell
          ibn = p_gridinfo4ray%cellneighbors_blk(jstencil,jc,jb) ! jb index of neighboring cell

          ! Loop over ray volumes in neighboring cell
          DO jray_nb = 1, nrays

            ! Skip calculations if this ray volume in the neighboring cell
            ! does not exist or is below the launch level z_src
            IF (p_ray%iexist(icn,jray_nb,ibn) == 0 .OR. &
                p_ray%specid(icn,jray_nb,ibn) < 0) CYCLE

            ! Ray volume center point in Cartesian coordinates
            ray_xyz(1) = p_ray%lon(icn,jray_nb,ibn)
            ray_xyz(2) = p_ray%lat(icn,jray_nb,ibn)
            ray_xyz(3) = p_ray%z  (icn,jray_nb,ibn) ! not used on the plane

            intriangle = inside_triangle_torus(ray_xyz, v1,v2,v3)

            IF (intriangle >= 0) THEN

              jk0 = p_ray%iexist(icn,jray_nb,ibn)
              jk0m1 = MAX(1, jk0 - 1)

              IF ( z_mc(icn,jk0,ibn) /= z_mc(jc,jk0,jb) .OR. z_mc(icn,jk0m1,ibn) /= z_mc(jc,jk0m1,jb) ) THEN
                !
                ! Some conditions used at the end of the integration (in propagate_wave) need
                ! to be applied again here, as the vertical grids are changed.
                !  - lower bound used in propagate_wave
!               IF (p_ray%z(icn,jray_nb,ibn) < z_ifc(jc, p_patch% nlevp1, jb)) THEN
                IF (p_ray%z(icn,jray_nb,ibn) < z_ifc(jc, p_patch% nlevp1 - nlev_lbnd(jg), jb)) THEN
                  ! outofdomain
!                 p_ray%iexist(icn,jray_nb,ibn) = 0      ! CHECK: problematic for OpenMP ?
!                 p_ray%specid(icn,jray_nb,ibn) = 0      !
                  p_ray%wadens(icn,jray_nb,ibn) = 0.     ! - This must be okay
                  CYCLE
                END IF
                !
                ! update jk0
                CALL get_jkhalf_closest_gen( jk0, p_ray% z(icn,jray_nb,ibn),  &
                  &                          p_patch% nlev, z_mc(:,:,jb), jc )
              END IF

              ! Find the next non-existing local ray volume index, which can be used
              ! for hosting the new ray volume moved here from neighboring cell
              jray_loc_old = jray_loc
              DO jr = jray_loc_old+1, nrays
                IF (p_ray%iexist(jc,jr,jb) == 0) THEN
                  jray_loc = jr
                  EXIT
                END IF
              ENDDO
              IF (jray_loc == jray_loc_old)  CALL finish ('regrid_wave', &
                  'Something is strange in horizontal regridding: jray_loc>nrays')

              ! Initialize new ray volume
              ! Synchronized properties: take them from the ray volume that moved 
              ! here from neighboring cell
              p_ray%iexist(jc,jray_loc,jb)    = jk0
              p_ray%specid(jc,jray_loc,jb)    = p_ray%specid(icn,jray_nb,ibn)
              p_ray%jk_active(jc,jray_loc,jb) = p_ray%jk_active(icn,jray_nb,ibn)
              p_ray%lon(jc,jray_loc,jb)       = p_ray%lon(icn,jray_nb,ibn)
              p_ray%lat(jc,jray_loc,jb)       = p_ray%lat(icn,jray_nb,ibn)
              p_ray%z(jc,jray_loc,jb)         = p_ray%z(icn,jray_nb,ibn)
              p_ray%dlon(jc,jray_loc,jb)      = p_ray%dlon(icn,jray_nb,ibn)
              p_ray%dlat(jc,jray_loc,jb)      = p_ray%dlat(icn,jray_nb,ibn)
              p_ray%dz(jc,jray_loc,jb)        = p_ray%dz(icn,jray_nb,ibn)
              p_ray%coslat(jc,jray_loc,jb)    = p_ray%coslat(icn,jray_nb,ibn)
              p_ray%k(jc,jray_loc,jb)         = p_ray%k(icn,jray_nb,ibn)
              p_ray%l(jc,jray_loc,jb)         = p_ray%l(icn,jray_nb,ibn)
              p_ray%m(jc,jray_loc,jb)         = p_ray%m(icn,jray_nb,ibn)
              p_ray%dk(jc,jray_loc,jb)        = p_ray%dk(icn,jray_nb,ibn)
              p_ray%dl(jc,jray_loc,jb)        = p_ray%dl(icn,jray_nb,ibn)
              p_ray%dm(jc,jray_loc,jb)        = p_ray%dm(icn,jray_nb,ibn)
              p_ray%wadens(jc,jray_loc,jb)    = p_ray%wadens(icn,jray_nb,ibn)

              ! Remove from the neighboring cell
!             p_ray%iexist(icn,jray_nb,ibn) = 0      ! CHECK: problematic for OpenMP ?
!             p_ray%specid(icn,jray_nb,ibn) = 0      !
              p_ray%wadens(icn,jray_nb,ibn) = 0._wp  ! - This must be okay

              ! Horizontal propagation test
              IF (ltest_hprop .AND. l1ray) THEN
                ! Open output file for horizontal propagation test
                OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
                WRITE(987,'(a,3i4,2x,6F14.5)') &
                       'Regrid ray taken from x, y:', jstencil, icn, ibn, &
                                          p_gridinfo4ray%cellvertices_x(1,icn,ibn),&
                                          p_gridinfo4ray%cellvertices_y(1,icn,ibn),&
                                          p_gridinfo4ray%cellvertices_x(2,icn,ibn),&
                                          p_gridinfo4ray%cellvertices_y(2,icn,ibn),&
                                          p_gridinfo4ray%cellvertices_x(3,icn,ibn),&
                                          p_gridinfo4ray%cellvertices_y(3,icn,ibn)
                WRITE(987,'(a,2i4,2x,6F14.5)') &
                       'Regrid ray taken   to x, y:', jc, jb, &
                                          p_gridinfo4ray%cellvertices_x(1,jc,jb),&
                                          p_gridinfo4ray%cellvertices_y(1,jc,jb),&
                                          p_gridinfo4ray%cellvertices_x(2,jc,jb),&
                                          p_gridinfo4ray%cellvertices_y(2,jc,jb),&
                                          p_gridinfo4ray%cellvertices_x(3,jc,jb),&
                                          p_gridinfo4ray%cellvertices_y(3,jc,jb)
                CLOSE(987)
              ENDIF ! ltest_hprop .AND. l1ray

            ENDIF ! intriangle

          ENDDO  ! jray_nb

        ENDDO  ! jstencil

      ENDDO  ! jc

    ENDDO  ! jb

!!$OMP END DO
!!$OMP END PARALLEL

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb, jc, jray_loc, i_startidx, i_endidx, jstencil, icn, ibn) ICON_OMP_GUIDED_SCHEDULE

    ! Find ray volumes that propagated out of this 
    ! cell and remove them
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx

        ! Vertex 1 in Cartesian coordinates
        v1%x = p_gridinfo4ray%cellvertices_x(1,jc,jb)
        v1%y = p_gridinfo4ray%cellvertices_y(1,jc,jb)
        v1%z = p_gridinfo4ray%cellvertices_z(1,jc,jb) ! not used on the plane

        ! Vertex 2 in Cartesian coordinates
        v2%x = p_gridinfo4ray%cellvertices_x(2,jc,jb)
        v2%y = p_gridinfo4ray%cellvertices_y(2,jc,jb)
        v2%z = p_gridinfo4ray%cellvertices_z(2,jc,jb) ! not used on the plane

        ! Vertex 3 in Cartesian coordinates
        v3%x = p_gridinfo4ray%cellvertices_x(3,jc,jb)
        v3%y = p_gridinfo4ray%cellvertices_y(3,jc,jb)
        v3%z = p_gridinfo4ray%cellvertices_z(3,jc,jb) ! not used on the plane

        DO jray_loc = 1, nrays

          IF (p_ray%iexist(jc,jray_loc,jb) == 0 .OR. &
              p_ray%specid(jc,jray_loc,jb) < 0) CYCLE

          IF (p_ray%wadens(jc,jray_loc,jb) == 0._wp) THEN
            ! Ray volume that has escaped to another prognostic cell
            p_ray%iexist(jc,jray_loc,jb) = 0
            p_ray%specid(jc,jray_loc,jb) = 0
            CYCLE
!         ELSE IF ( this_cell_is_not_a_neighbor_of_halo ) THEN
!           ! Ray volume inside this cell, far from halo: skip the in-triangle test
!           CYCLE
          END IF

          ! Remainder: ray volume either being inside this cell or escaped to a halo cell

          ! Find ray volumes that propagated out of this cell and 
          ! remove them
          ! Ray volume center point in Cartesian coordinates
          ray_xyz(1) = p_ray%lon(jc,jray_loc,jb)
          ray_xyz(2) = p_ray%lat(jc,jray_loc,jb)
          ray_xyz(3) = p_ray%  z(jc,jray_loc,jb) ! not used on the plane

          intriangle = inside_triangle_torus(ray_xyz, v1,v2,v3)

          IF (intriangle < 0) THEN  ! not in triangle

            ! Remove ray volume
            p_ray%iexist(jc,jray_loc,jb) = 0
            p_ray%specid(jc,jray_loc,jb) = 0
            p_ray%wadens(jc,jray_loc,jb) = 0._wp

            ! Horizontal propagation test
            IF (ltest_hprop .AND. l1ray) THEN
              ! Open output file for horizontal propagation test
              OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
              WRITE(987,'(a,2i4,2x,6F14.5)') &
                     'Regrid ray removed at x, y:', jc, jb, &
                                        p_gridinfo4ray%cellvertices_x(1,jc,jb),&
                                        p_gridinfo4ray%cellvertices_y(1,jc,jb),&
                                        p_gridinfo4ray%cellvertices_x(2,jc,jb),&
                                        p_gridinfo4ray%cellvertices_y(2,jc,jb),&
                                        p_gridinfo4ray%cellvertices_x(3,jc,jb),&
                                        p_gridinfo4ray%cellvertices_y(3,jc,jb)
              CLOSE(987)
            ENDIF ! ltest_hprop .AND l1ray

          ENDIF ! not in triangle

        ENDDO  ! jray_loc

      ENDDO  ! jc

    ENDDO  ! jb

!!$OMP END DO
!!$OMP END PARALLEL

  ENDIF ! is_plane_torus?

END SUBROUTINE regrid_wave
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE datout(nlev,idiag,zray,mray,dzray,dmray,dens,specid,zz_half,zz,rho,fld_u,fld_v,&
                  uw_wr,vw_wr,bvf2,kd,kd2,ld,ld2,md,md2,cgz_diag,A,A_save,&
                  B2,B2_save,mB2,mB2_save,ener)

  INTEGER,  INTENT(IN) :: nlev
  INTEGER,  INTENT(IN) :: idiag
  ! Ray volume data
  REAL(wp), INTENT(IN) :: zray(:), mray(:), dens(:)
  REAL(wp), INTENT(IN) :: dzray(:), dmray(:)
  INTEGER , INTENT(IN) :: specid(:)
  ! Gridded data
  REAL(wp), INTENT(IN) :: fld_u(:), fld_v(:)
  REAL(wp), INTENT(IN) :: zz_half(:), zz(:)
  REAL(wp), INTENT(IN) :: uw_wr(:), vw_wr(:)
  REAL(wp), INTENT(IN) :: rho(:), ener(:)
  REAL(wp), INTENT(IN) :: bvf2(:)
  REAL(wp), INTENT(IN) :: kd(:), ld(:), md(:)
  REAL(wp), INTENT(IN) :: kd2(:), ld2(:), md2(:), cgz_diag(:)
  REAL(wp), INTENT(IN) :: B2(:), B2_save(:)
  REAL(wp), INTENT(IN) :: A(:), A_save(:)
  REAL(wp), INTENT(IN) :: mB2(:), mB2_save(:)
  ! Others 
  INTEGER ifield,irc30,ifield2,irc40,jk
  REAL(wp) fld_uinv(nlev+1), fld_vinv(nlev+1)
  REAL(wp) zz_half_inv(nlev+1), zz_inv(nlev)
  REAL(wp) uw_wrinv(nlev+1), vw_wrinv(nlev+1)
  REAL(wp) rho_inv(nlev+1), enerinv(nlev+1)
  REAL(wp) bvf_wrinv(nlev+1)
  REAL(wp) kd_wrinv(nlev+1), ld_wrinv(nlev+1), md_wrinv(nlev+1)
  REAL(wp) kd2_wrinv(nlev+1), ld2_wrinv(nlev+1), md2_wrinv(nlev+1)
  REAL(wp) cgz_diag_wrinv(nlev+1)
  REAL(wp) B2_wrinv(nlev+1), B2_save_wrinv(nlev+1)
  REAL(wp) A_wrinv(nlev+1), A_save_wrinv(nlev+1)
  REAL(wp) mB2_wrinv(nlev+1), mB2_save_wrinv(nlev+1)
  REAL(wp) fldsgl(nrays(jg))
  REAL(wp) fldsgl2(nlev)

  IF (msg_level >= 12) CALL message('datout', 'MS-GWaM: 1D diagnostic output')

  !----------------------------------------------------------------------
  ! Purpose:
  !         Output binary files to track subgrid-scale quantities and the 
  !         movement of ray volume center-points in z-m phase at selected 
  !         lat,lon points.
  !
  ! Method:
  !         -- Invert vertical index
  !         -- output diagnostics (ldiagprof=.T.) to binary files
  !----------------------------------------------------------------------

  ! Invert over vertical layers
  DO jk = 1,nlev+1
    zz_half_inv(jk)    = zz_half(nlev+2-jk)
    rho_inv(jk)        = rho(nlev+2-jk)
    enerinv(jk)        = ener(nlev+2-jk)
    uw_wrinv(jk)       = uw_wr(nlev+2-jk)
    vw_wrinv(jk)       = vw_wr(nlev+2-jk)
    bvf_wrinv(jk)      = SQRT(bvf2(nlev+2-jk))
    kd_wrinv(jk)       = kd(nlev+2-jk)
    kd2_wrinv(jk)      = kd2(nlev+2-jk)
    ld_wrinv(jk)       = ld(nlev+2-jk)
    ld2_wrinv(jk)      = ld2(nlev+2-jk)
    md_wrinv(jk)       = md(nlev+2-jk)
    md2_wrinv(jk)      = md2(nlev+2-jk)
    cgz_diag_wrinv(jk) = cgz_diag(nlev+2-jk)
    B2_wrinv(jk)       = B2(nlev+2-jk)
    B2_save_wrinv(jk)  = B2_save(nlev+2-jk)
    A_wrinv(jk)        = A(nlev+2-jk)
    A_save_wrinv(jk)   = A_save(nlev+2-jk)
    mB2_wrinv(jk)      = mB2(nlev+2-jk)
    mB2_save_wrinv(jk) = mB2_save(nlev+2-jk)
  ENDDO
  DO jk = 1,nlev
    zz_inv(jk)   = zz(nlev+1-jk)
    fld_uinv(jk) = fld_u(nlev+1-jk)
    fld_vinv(jk) = fld_v(nlev+1-jk)
  ENDDO
  fld_uinv(nlev+1) = 0.   ! temporary treatment
  fld_vinv(nlev+1) = 0.

!  ! output ray data
!  OPEN(30,FILE=TRIM(filename_diag1(idiag)),FORM='unformatted',&
!       ACCESS='direct',RECL=2*nrays(jg))
!  irc30=0
  
  ! output data defined on grid
  OPEN(40,FILE=TRIM(filename_diag2(idiag)),FORM='unformatted',&
       ACCESS='direct',RECL=2*nlev)
  irc40=0
  
!  ! Write out ray data
!  ! Position, wavenumber, wave action density, dm
!  irc30=(iout_msgwam-1)*6
!  
!  DO ifield=1,6
!     irc30=irc30+1
!  
!     IF(ifield==1) fldsgl(:) = zray(:)
!     IF(ifield==2) fldsgl(:) = mray(:)
!     IF(ifield==3) fldsgl(:) = dens(:)
!     IF(ifield==4) fldsgl(:) = dzray(:)
!     IF(ifield==5) fldsgl(:) = dmray(:)
!     IF(ifield==6) fldsgl(:) = REAL(specid(:),wp)
!  
!     WRITE(30,REC=irc30)fldsgl
!  ENDDO
  
  ! Write out data on grid
  ! Mean-flow, momentum flux density, energy density 
  irc40=(iout_msgwam-1)*22
  
  DO ifield2=1,22
     irc40=irc40+1
  
     IF(ifield2==1)  fldsgl2(1:nlev)=zz_half_inv(2:nlev+1)
     IF(ifield2==2)  fldsgl2(1:nlev)=zz_inv(1:nlev)
     IF(ifield2==3)  fldsgl2(1:nlev)=rho_inv(2:nlev+1)
     IF(ifield2==4)  fldsgl2(1:nlev)=fld_uinv(2:nlev+1)
     IF(ifield2==5)  fldsgl2(1:nlev)=fld_vinv(2:nlev+1)
     IF(ifield2==6)  fldsgl2(1:nlev)=uw_wrinv(2:nlev+1)
     IF(ifield2==7)  fldsgl2(1:nlev)=vw_wrinv(2:nlev+1)
     IF(ifield2==8)  fldsgl2(1:nlev)=enerinv(2:nlev+1)
     IF(ifield2==9)  fldsgl2(1:nlev)=bvf_wrinv(2:nlev+1)
     IF(ifield2==10) fldsgl2(1:nlev)=kd_wrinv(2:nlev+1)
     IF(ifield2==11) fldsgl2(1:nlev)=kd2_wrinv(2:nlev+1)
     IF(ifield2==12) fldsgl2(1:nlev)=ld_wrinv(2:nlev+1)
     IF(ifield2==13) fldsgl2(1:nlev)=ld2_wrinv(2:nlev+1)
     IF(ifield2==14) fldsgl2(1:nlev)=md_wrinv(2:nlev+1)
     IF(ifield2==15) fldsgl2(1:nlev)=md2_wrinv(2:nlev+1)
     IF(ifield2==16) fldsgl2(1:nlev)=cgz_diag_wrinv(2:nlev+1)
     IF(ifield2==17) fldsgl2(1:nlev)=A_wrinv(2:nlev+1)
     IF(ifield2==18) fldsgl2(1:nlev)=A_save_wrinv(2:nlev+1)
     IF(ifield2==19) fldsgl2(1:nlev)=B2_wrinv(2:nlev+1)
     IF(ifield2==20) fldsgl2(1:nlev)=B2_save_wrinv(2:nlev+1)
     IF(ifield2==21) fldsgl2(1:nlev)=mB2_wrinv(2:nlev+1)
     IF(ifield2==22) fldsgl2(1:nlev)=mB2_save_wrinv(2:nlev+1)
  
     WRITE(40,rec=irc40)fldsgl2
  ENDDO
  
!  CLOSE(30)
  CLOSE(40)

END SUBROUTINE datout
!!
!!-------------------------------------------------------------------------
!!
ELEMENTAL REAL(wp) FUNCTION dyn_visc_sutherland(Ta)
  !
  ! Copied here from mo_2mom_mcrph_util.f90
  !
  ! Calculate dynamic viscosity of air [kg m-1 s-1]
  ! following Sutherland's formula of an ideal
  ! gas with reference temp. T = 291.15 K
  !
  ! There is another alternative in P&K97 on
  ! page 417
  !
  IMPLICIT NONE
  REAL(wp), INTENT(IN) :: Ta   ! ambient temp. [K]
  REAL(wp), PARAMETER :: &
       C = 120._wp     , &     ! Sutherland's constant (for air) [K]
       T0 = 291.15_wp  , &     ! Reference temp. [K]
       eta0 = 1.827E-5_wp      ! Reference dyn. visc. [kg m-1 s-1]
  REAL(wp), PARAMETER :: eta0_a = eta0*(T0 + C)
  REAL(wp) :: r

  r = Ta/T0
  dyn_visc_sutherland = eta0_a/(Ta + C) * r*SQRT(r)

  RETURN
END FUNCTION dyn_visc_sutherland
!!
!!-------------------------------------------------------------------------
!!
ELEMENTAL REAL(wp) FUNCTION signwn(k)
  !
  ! Get sign of the wave number 
  ! or zero if it is zero
  !
  IMPLICIT NONE
  REAL(wp), INTENT(IN) :: k

  IF (k /= 0._wp) THEN
    signwn = SIGN(1._wp, k)
  ELSE
    signwn = 0.
  ENDIF
  
  RETURN
END FUNCTION signwn
!!
!!-------------------------------------------------------------------------
!!
FUNCTION inside_triangle_torus(v, v1,v2,v3)
  INTEGER :: inside_triangle_torus
  REAL(wp),       INTENT(IN)     :: v(3)  ! ray volume center point
  TYPE(t_point),  INTENT(IN)  :: v1,v2,v3 ! grid cell vertices
  ! local variables
  INTEGER :: jvec, jmin, jmax, jthird
  REAL(wp)  :: vecx(3), vecy(3)           ! x, y of grid cell vertices
  REAL(wp)  :: a, b, invdenom

  ! Pick x coordinates of cell vertices
  vecx = (/v1%x, v2%x, v3%x/)

  ! If the cell has vertices corresponding to both sides of the 
  ! domain, take action... some coordinates will have to be replaced
  ! before we actually perform any "in triangle" test.
  IF (MAXVAL(vecx) == domxmax .AND. MINVAL(vecx) == domxmin) THEN
    ! Calculate vertex indices corresponding to domxmin, domxmax 
    ! and that of the third one.
    DO jvec = 1,3
      IF (vecx(jvec) == domxmin) jmin = jvec
      IF (vecx(jvec) == domxmax) jmax = jvec
      IF (vecx(jvec) /= domxmin .AND. &
          vecx(jvec) /= domxmax) jthird = jvec
    ENDDO
    IF (ABS(v(1) - domxmin) < ABS(v(1) - domxmax)) THEN
      ! The ray volume is close to the left side of the domain
      vecx(jmax)   = domxmin - edge_length_triangle * 0.5_wp
      IF (ABS(vecx(jthird) - domxmax) < ABS(vecx(jthird) - domxmin)) THEN
        ! The cell in which we search is on the right side of the domain
        vecx(jthird) = domxmin - edge_length_triangle
      ENDIF
    ELSE ! ABS(v(1) - domxmin) > ABS(v(1) - domxmax)
      ! The ray volume is close to the right side of the domain
      vecx(jmin)   = domxmax + edge_length_triangle * 0.5_wp
      IF (ABS(vecx(jthird) - domxmax) > ABS(vecx(jthird) - domxmin)) THEN
        ! The cell in which we search is on the left side of the domain
        vecx(jthird) = domxmax + edge_length_triangle
      ENDIF
    ENDIF
  ENDIF

  ! Pick y coordinates of cell vertices
  vecy = (/v1%y, v2%y, v3%y/)

  ! If the cell has vertices corresponding to both sides of the 
  ! domain, take action... some coordinates will have to be replaced
  ! before we actually perform any "in triangle" test.
  IF (MAXVAL(vecy) == domymax .AND. MINVAL(vecy) == domymin) THEN
    IF (ABS(v(2) - domymin) < ABS(v(2) - domymax)) THEN
      ! The ray volume is close to the front side of the domain
      DO jvec = 1,3
        IF (vecy(jvec) == domymax) vecy(jvec) = domymin - height_triangle
      ENDDO
    ELSE ! ABS(v(2) - domymin) > ABS(v(2) - domymax)
      ! The ray volume is close to the back side of the domain
      DO jvec = 1,3
        IF (vecy(jvec) == domymin) vecy(jvec) = domymax + height_triangle
      ENDDO
    ENDIF
  ENDIF

  ! Search by barycentric coordinates
  ! Reference: http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html

  invdenom =  1._wp &
           / ((vecy(2) - vecy(3)) * (vecx(1) - vecx(3)) &
           +  (vecx(3) - vecx(2)) * (vecy(1) - vecy(3)))
  a = ((vecy(2) - vecy(3)) * (v(1) - vecx(3)) &
    +  (vecx(3) - vecx(2)) * (v(2) - vecy(3))) * invdenom
  b = ((vecy(3) - vecy(1)) * (v(1) - vecx(3)) &
    +  (vecx(1) - vecx(3)) * (v(2) - vecy(3))) * invdenom

    IF (a >= 0._wp .AND. b >= 0._wp .AND. a + b < 1._wp) THEN
      inside_triangle_torus = 1
    ELSE
      inside_triangle_torus = -1
    ENDIF

END FUNCTION inside_triangle_torus
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_jkhalf_closest_gen(jk0, z_target, nlev, z_mc, jc)

  ! To find the half-level index closest to a target altitude, from its initial guess,
  ! when it is unknown whether the initial guess is higher or lower than the target.
  ! The result must satisfy both of the following conditions:
  !   z_mc(jc,jk0) <= z_target                   (if  jk0 < nlev+1)
  !                   z_target < z_mc(jc,jk0-1)  (if  jk0 > 1     )

  INTEGER , INTENT(INOUT) :: jk0          ! in: initial guess of jk0
  INTEGER , INTENT(IN   ) :: nlev         ! number of full levels
  INTEGER , INTENT(IN   ) :: jc           ! chunk index
  REAL(wp), INTENT(IN   ) :: z_target     ! target altitude
  REAL(wp), INTENT(IN   ) :: z_mc(:,:)    ! full-level altitudes (nproma, nlev)

  INTEGER :: jk

  IF (jk0 <= nlev) THEN
    IF (z_mc(jc,jk0) > z_target) THEN
      IF (z_mc(jc,nlev) > z_target) THEN
        jk0 = nlev+1
      ELSE
        DO jk = jk0+1, nlev
          IF (z_mc(jc,jk) <= z_target) THEN
            jk0 = jk  ;  EXIT
          END IF
        ENDDO
      END IF
    END IF
  END IF
  IF (jk0 > 1) THEN
    IF (z_mc(jc,jk0-1) <= z_target) THEN
      IF (z_mc(jc,1) <= z_target) THEN
        jk0 = 1
      ELSE
        DO jk = jk0-2, 1, -1
          IF (z_mc(jc,jk) > z_target) THEN
            jk0 = jk+1  ;  EXIT
          END IF
        ENDDO
      END IF
    END IF
  END IF

END SUBROUTINE get_jkhalf_closest_gen
!!
!!-------------------------------------------------------------------------
!!
FUNCTION test_jkray( jkhalf, jc, zray, z_mc, jkhmax )
  INTEGER :: test_jkray

  INTEGER , INTENT(IN) :: jkhalf
  INTEGER , INTENT(IN) :: jkhmax
  INTEGER , INTENT(IN) :: jc
  REAL(wp), INTENT(IN) :: zray
  REAL(wp), INTENT(IN) :: z_mc(:,:)

  IF ( jkhalf < 1 .OR. jkhalf > jkhmax ) THEN
    IF (jkhalf < 1) THEN
      test_jkray = -1
    ELSE
      test_jkray = -2
    END IF
    RETURN
  END IF
  test_jkray = 0
  IF (jkhalf > 1) THEN
    IF (zray >= z_mc(jc,jkhalf-1)) THEN
      test_jkray = 1  ;  RETURN
    END IF
  END IF
  IF (jkhalf < jkhmax) THEN
    IF (zray < z_mc(jc,jkhalf)) THEN
      test_jkray = 2  ;  RETURN
    END IF
  END IF

END FUNCTION test_jkray

SUBROUTINE test_blowup_uvt(p_patch, dudt, dvdt, dtdt, thr_u, thr_t, stab_info)

  TYPE(t_patch), TARGET, INTENT(IN) :: p_patch              ! grid/patch info.
  REAL(wp),              INTENT(IN) :: dudt(:,:,:)
  REAL(wp),              INTENT(IN) :: dvdt(:,:,:)
  REAL(wp),              INTENT(IN) :: dtdt(:,:,:)
  REAL(wp),              INTENT(IN) :: thr_u
  REAL(wp),              INTENT(IN) :: thr_t
  CHARACTER (len=*),     INTENT(IN) :: stab_info

  INTEGER  :: rl_start, rl_end, i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER  :: jc, jk, jb

  rl_start   = grf_bdywidth_c + 1
  rl_end     = min_rlcell_int
  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx)

  DO jb = i_startblk, i_endblk
    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)
    DO jk = 1, p_patch%nlev
      DO jc = i_startidx, i_endidx
        IF ( ABS(dudt(jc,jk,jb)) < thr_u .AND. ABS(dvdt(jc,jk,jb)) < thr_u .AND. &
          &  ABS(dtdt(jc,jk,jb)) < thr_t )  CYCLE
        WRITE(message_text,'(a,i4,2f8.2,f8.4,a,2f6.1,2a)') 'MS-GWaM stability check: ',  &
          &  jk, ABS(dudt(jc,jk,jb)), ABS(dvdt(jc,jk,jb)), ABS(dtdt(jc,jk,jb)),  &
          &  ': jk, du/dt, dv/dt, dT/dt |',  &
          &  p_patch%cells%center(jc,jb)%lat*rad2deg, &
          &  p_patch%cells%center(jc,jb)%lon*rad2deg, ': (lat, lon) |', &
          &  stab_info 
        CALL message('', TRIM(message_text))
!       CALL finish('test_blowup_uvt', 'Too large tendencies are detected in MS-GWaM')
      ENDDO
    ENDDO
  ENDDO

!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE test_blowup_uvt
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_msgwam
