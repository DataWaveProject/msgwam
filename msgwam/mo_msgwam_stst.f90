!>
!! This module is a collection of subroutines running a steady state 
!! GW drag parametrization (MS-GWaM-ST). This parametrization is
!! based on:
!! 1) initialization of GW spectral elements (same as in MS-GWaM and 
!!    Orr et al., 2010)
!! 2) equilibrium profiles in a constant background flow
!! 3) a critical layer parametrization
!! 4) a reflection layer parametrization
!! 5) the same wave breaking parametrization as used by MS-GWaM
!! 6) the same tendency calculations as in MS-GWaM
!!
!! @author Young-Ha Kim, Goethe Uni Frankfurt
!!
!! @par Revision History
!! Initial version by Young-Ha Kim, Goethe Uni Frankfurt (2019-06-xx)

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_msgwam_stst

  USE mo_kind,                  ONLY: wp, vp
  USE mo_mpi,                   ONLY: my_process_is_stdio, p_wait, work_mpi_barrier
  USE mo_exception,             ONLY: message, finish, message_text, &
                                      debug_messages_on, debug_messages_off
  USE mo_physical_constants,    ONLY: grav, rd, cpd
  USE mo_grid_config,           ONLY: grid_sphere_radius

  USE mo_model_domain,          ONLY: t_patch
  USE mo_vertical_coord_table,  ONLY: vct_a

  USE mo_impl_constants,        ONLY: min_rlcell_int, success
  USE mo_impl_constants_grf,    ONLY: grf_bdywidth_c
  USE mo_loopindices,           ONLY: get_indices_c

  USE mo_nonhydro_types,        ONLY: t_nh_diag, t_nh_prog, t_nh_metrics
  USE mo_parallel_config,       ONLY: nproma
  USE mo_run_config,            ONLY: msg_level, ldynamics
  USE mo_dynamics_config,       ONLY: lcoriolis
  USE mo_math_constants,        ONLY: pi, rad2deg, deg2rad
  USE mo_util_vgrid_types,      ONLY: vgrid_buffer
  USE mtime,                    ONLY: datetime, timeDelta, newTimedelta, &
                                      deallocateTimedelta, getTimedeltaFromDatetime, &
                                      getTotalMillisecondsTimedelta
  USE mo_nonhydrostatic_config, ONLY: damp_height
  USE mo_vertical_grid,         ONLY: nrdmax
  USE mo_timer
  USE mo_msgwam_config
  USE mo_setup_msgwam_interface
  USE mo_gw_source_config,      ONLY: gws_conv_config

  ! Modules for computing the Orr et al. 2010 launch spectrum
  USE data_gwd,    ONLY : nslope, gfluxlaun, &
    &                     ggaussa, ggaussb, ngauss, gcoeff, lozpr

  IMPLICIT NONE

  PRIVATE

  INTEGER ::  jg

! REAL(wp), PARAMETER :: bvf2_min = 1.E-12_wp    ! minimum for BVF (avoid division by 0)
  REAL(wp), PARAMETER :: bvf2_min = 3.E-8_wp     ! minimum for BVF (avoid division by 0)
                                                 ! set to be larger than f^2 everywhere

  PUBLIC  ::  gwdrag_msgwam_stst

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE gwdrag_msgwam_stst ( dt_call,                   & ! input
                            &   mtime_datetime,            & ! input
                            &   p_sim_time,                & ! input
                            &   p_patch,p_metrics,         & ! input
                            &   rho,                       & ! input
                            &   u       ,                  & ! input
                            &   v       ,                  & ! input
                            &   temp    ,                  & ! input
                            &   temp_ifc,                  & ! input
                            &   p_fld                      ) ! inout

  TYPE(datetime),      POINTER,INTENT(IN)    :: mtime_datetime       ! date/time information
  TYPE(t_patch),        TARGET,INTENT(IN)    :: p_patch              ! grid/patch info.
  TYPE(t_nh_metrics)          ,INTENT(IN)    :: p_metrics
  TYPE(t_msgwam),              INTENT(INOUT) :: p_fld                ! the atm phys vars
  REAL(wp),                    INTENT(IN)    :: dt_call              ! time step
  REAL(wp),                    INTENT(IN)    :: p_sim_time           ! elapsed simulation time on this grid level
  REAL(wp),            POINTER,INTENT(IN)    :: rho (:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: u   (:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: v   (:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: temp(:,:,:)
  REAL(wp),            POINTER,INTENT(IN)    :: temp_ifc(:,:,:)
  ! Local array bounds:
  INTEGER                                    :: nlev, nlevp1         ! number of full and half levels
  INTEGER                                    :: rl_start, rl_end
  INTEGER                                    :: i_startblk, i_endblk ! blocks
  INTEGER                                    :: i_startidx, i_endidx ! slices
  INTEGER                                    :: jk,jc,jb             ! block indeces
  INTEGER                                    :: maxindexu(3), maxindexv(3)
  REAL(wp)                                   :: rho_half(nproma,p_patch%nlevp1) ! rho at half levels
  REAL(wp)                                   :: temp_half(nproma,p_patch%nlevp1) ! temperature at half levels
  REAL(wp)                                   :: u_half(nproma,p_patch%nlevp1) ! u at half levels
  REAL(wp)                                   :: v_half(nproma,p_patch%nlevp1) ! v at half levels
  REAL(wp)                                   :: uwflux_cgw(nproma,p_patch%nlevp1) ! uw momentum flux
  REAL(wp)                                   :: vwflux_cgw(nproma,p_patch%nlevp1) ! vw momentum flux
  ! Background
  REAL(wp)                                   :: bvf2(nproma,p_patch%nlevp1)    ! Brunt-Väisälä freq**2
  REAL(wp)                                   :: kvisc(nproma,p_patch%nlevp1)   ! kinematic viscosity
  REAL(wp)                                   :: gammash2(nproma,p_patch%nlevp1)! inverse pinc scale height
  ! 1D Diagnostics
  REAL(wp)                                   :: kd(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: kd2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: ld(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: ld2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: md(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: md2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: cgz_diag(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: mB2(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: mB2_save(nproma,p_patch%nlevp1)
  REAL(wp)                                   :: prec(nproma)                                  ! total precip
  REAL(wp)                                   :: sponge_cos(p_patch%nlevp1)
! REAL(wp),                      PARAMETER   :: lim_tend = 0.05_wp
  ! 
  INTEGER                                    :: idiag

  CALL message('gwdrag_msgwam_stst', 'Steady-state MS-GWaM-1D called')

  ! Number of vertical levels
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  ! Domain ID
  jg     = p_patch%id

  ! Exclude boundary interpolation zone of nested domains
  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)


  DO jk = 2, nrdmax(jg)
    sponge_cos(jk) = COS( (vct_a(jk+p_patch%nshift_total) - damp_height(jg))  &
      &                   /(vct_a(1) - damp_height(jg))*(0.5_wp*pi) )**2
  ENDDO
  sponge_cos(1) = 0._wp
  sponge_cos(nrdmax(jg)+1:) = 1._wp


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,rho_half,temp_half,u_half,v_half,  &
!$OMP            uwflux_cgw,vwflux_cgw,bvf2,kvisc,gammash2,                      &
!$OMP            kd,kd2,ld,ld2,md,md2,cgz_diag,mB2,mB2_save,prec)

  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
      & i_startidx, i_endidx, rl_start, rl_end)


    !=========================== FIELDS/GRADS ===========================
    ! Calculation of all resolved variables and their vertical gradients.
    !====================================================================

    IF (timers_level > 4) CALL timer_start(timer_msgwam_fieldsgrads)

    ! Calculate derived fields from the resolved flow (e.g. Brunt-Väisälä freq and density scale height)
    ! and their vertical gradients
    CALL fieldsgrads(nlev       = nlev,                          & ! no. of full levels                (in)
                     i_startidx = i_startidx,                    & ! first index of the block          (in)
                     i_endidx   = i_endidx,                      & ! last index of the block           (in)
                     z          = p_metrics%z_mc(:,:,jb),        & ! full level heights                (in)
                     zhalf      = p_metrics%z_ifc(:,:,jb),       & ! half level heights                (in)
                     dz         = p_metrics%ddqz_z_full(:,:,jb), & ! full layer thickness              (in)
                     dzhalf     = p_metrics%ddqz_z_half(:,:,jb), & ! half layer thickness              (in)
                     wgtfac_c   = p_metrics%wgtfac_c(:, :,jb),   & !                                     (in)
                     rho        = rho (:, :,jb),                 & ! density at full levels              (in)
                     temp       = temp(:, :,jb),                 & ! temperature at full levels          (in)
                     u          = u(:, :,jb),                    & ! zonal wind (full levels)            (in)
                     v          = v(:, :,jb),                    & ! meridional wind (full levels)       (in)
                     rho_half   = rho_half(:,:),                 & ! density at half levels         (inout)
                     temp_half  = temp_half(:,:),                & ! temperature at half levels     (inout)
                     u_half     = u_half(:,:),                   & ! zonal wind (half levels)       (inout)
                     v_half     = v_half(:,:),                   & ! merid wind (half levels)       (inout)
                     bvf2       = bvf2(:,:),                     & ! Brunt-Väisälä frequency**2     (inout)
                     kvisc      = kvisc(:,:),                    & ! Kinematic viscosity            (inout)
                     gammash2   = gammash2(:,:)                  ) ! inverse pinc scaleheight**2    (inout)

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_fieldsgrads)


    !========================= INIT GW SOURCES ============================
    ! Launch GW sources (currently 2 kinds):
    ! 1) convective (init_gw_conv)
    ! 2) Desaubies type background (init_gw_orretal)
    ! TODO: mountain waves, jets/fronts
    !======================================================================

    IF ( gws_conv_config%n_source(jg) > 0 ) THEN
      IF ( p_sim_time /= 0._wp ) THEN
        CALL init_gw_conv( jb,i_startidx,i_endidx,jray_offset_cv(jg), p_fld%flag_cgw(:,jb) )
      ELSE
        p_ray(jg)%wadens(:,jray_offset_cv(jg)+1:jray_offset_cv(jg)+nrays_cv(jg),1) = 0._wp
      END IF
    END IF

    IF (nrays_add_bg(jg) /= 0) THEN
      IF (timers_level > 4) CALL timer_start(timer_msgwam_init_gw_orretal)

      ! Calculate total precipitation (input for the Orr et al., 2010 source)
!     IF ( lozpr )  prec(:) =   prm_diag%rain_gsp_rate(:,jb) &  ! rain_gsp
!       &                     + prm_diag%snow_gsp_rate(:,jb) &  ! snow_gsp
!       &                     + prm_diag%rain_con_rate(:,jb) &  ! rain_con
!       &                     + prm_diag%snow_con_rate(:,jb)    ! snow con
      IF ( lozpr )  prec(:) = 0._wp

      CALL init_gw_orretal(nlev        = nlev,                            & ! no. of full levels       (in)
                           i_startidx  = i_startidx,                      & ! first index of the block (in)
                           i_endidx    = i_endidx,                        & ! last index of the block  (in)
                           jray_offset = jray_offset_bg(jg),              &
                           mdatetime   = mtime_datetime,                  & ! last index of the block  (in)
                           kray        = p_ray(jg)%k(:,:,1),              & ! horiz (lon) wavenumber   (inout)
                           lray        = p_ray(jg)%l(:,:,1),              & ! horiz (lat) wavenumber   (inout)
                           mray        = p_ray(jg)%m(:,:,1),              & ! vertical wavenumber      (inout)
                           dens        = p_ray(jg)%wadens(:,:,1),         & ! wave action density      (inout)
                           iexist      = p_ray(jg)%iexist(:,:,1),         & ! existence of ray         (inout)
                           jklaunch    = jklaunch_bg(jg),                 & ! launch level index       (in)
                           lat         = p_patch%cells%center(:,jb)%lat,  & ! latitude at cell center  (in)
                           rhol        = rho_half(:,jklaunch_bg(jg)),     & ! rho at half levels       (in)
                           bvfl2       = bvf2(:,jklaunch_bg(jg)),         & ! Brunt-Väisala freq**2    (in)
                           lat_prof_bw = p_lfluxbg(jg)%lat_prof_bw(:,jb), & ! latitudinal factor     (in)
                           lat_prof_bs = p_lfluxbg(jg)%lat_prof_bs(:,jb), & ! latitudinal factor     (in)
                           lat_prof    = p_lfluxbg(jg)%lat_prof(:,jb),    & ! latitudinal factor     (inout)
                           prec        = prec(:)                          ) ! total precipitation      (in)

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_init_gw_orretal)
    ENDIF


    !============================ SATURATION ==============================
    ! Wave breaking scheme based on static instability criterion: if the GW 
    ! turns the potential temperture gradient to negative at a certain
    ! height, the wave action of ray volumes are reduced so that static 
    ! stability sets in again. Saturation is (should be!) calculated for 
    ! GWs from all sources (background, convective) together
    !======================================================================

    IF (lsaturation) THEN

      IF (timers_level > 4) CALL timer_start(timer_msgwam_saturation)

      ! In MS-GWaM-ST two saturation schemes have been implemented
      ! 1) non-monochromatic scheme: same as in the transient scheme i.e. GWs
      !    from the whole launch spectra are assumed to be present at all
      !    heights at the same time. Therefore the contribution from all these
      !    spectral elements to (mB)^2 is considered when saturation happens
      ! 2) monochromatic scheme: here - similarly to all steady-state schemes -
      !    it is assumed that GWs with different spectral properties are never
      !    collocated and therefore the saturation is performed element by
      !    element of the spectrum

      IF (.NOT. lsat_mono) THEN
        ! Non-monochromatic saturation scheme
        ! In this case saturation is (should be!) calculated for 
        ! GWs from all sources (background, convective) together

        CALL saturation(nlev         = nlev,                             & ! no. of full levels       (in)
                        i_startidx   = i_startidx,                       & ! first index of the block (in)
                        i_endidx     = i_endidx,                         & ! last index of the block  (in)
                        jray_start_cv = jray_offset_cv(jg)+1,            &
                        jray_end_cv   = jray_offset_cv(jg)+nrays_cv(jg), &
                        jray_start_bg = jray_offset_bg(jg)+1,            &
                        jray_end_bg   = jray_offset_bg(jg)+nrays_bg(jg), &
                        dz           = p_metrics%ddqz_z_full(:,:,jb),    & ! full layer thickness     (in)
                        kray         = p_ray(jg)%k(:,:,1),               & ! horiz (lon) wavenumber   (in)
                        lray         = p_ray(jg)%l(:,:,1),               & ! horiz (lat) wavenumber   (in)
                        mray         = p_ray(jg)%m(:,:,1),               & ! vertical wavenumber      (in)
                        dens         = p_ray(jg)%wadens(:,:,1),          & ! wave action density      (in)
                        iexist       = p_ray(jg)%iexist(:,:,1),          & ! existence of ray         (in)
                        u_half       = u_half(:,:),                      & ! zonal wind (half levels) (in)
                        v_half       = v_half(:,:),                      & ! merid wind (half levels) (in)
                        bvf2         = bvf2(:,:),                        & ! Brunt-Väisala freq**2    (in)
                        gammash2     = gammash2(:,:),                    & ! inverse pinc scale height(in)
                        fc           = p_patch%cells%f_c(:,jb),          & ! Coriolis parameter       (in)
                        rho          = rho_half(:,:),                    & ! rho at half levels       (in)
                        kvisc        = kvisc(:,:),                       & ! Kinematic viscosity      (in)
                        flag_cgw     = p_fld% flag_cgw(:,jb),            & ! flag for CGW grid        (in)
                        lcalc_flux_4dir_cv = lcalc_flux_4dir_cv(jg),     & ! flag to calc. 4-fluxes   (in)
                        lcalc_flux_4dir_bg = lcalc_flux_4dir_bg(jg),     & ! flag to calc. 4-fluxes   (in)
                        uwflux_cv     = uwflux_cgw(:,:),                 & ! uw mometum flux          (out)
                        vwflux_cv     = vwflux_cgw(:,:),                 & ! vw mometum flux          (out)
                        uwflux_bg     = p_fld% uwfl_mgm    (:,:,jb),     & ! uw mometum flux          (out)
                        vwflux_bg     = p_fld% vwfl_mgm    (:,:,jb),     & ! vw mometum flux          (out)
                        energy_bg     = p_fld% energy_mgm  (:,:,jb),     & ! wave energy              (out)
                        apmflux_bg    = p_fld% apmfl_mgm   (:,:,jb),     & ! abs momentum flux        (out)
                        energy_p_bg   = p_fld% energy_p_mgm(:,:,jb),     & ! wave potential energy    (out)
                        pmflux_e_bg   = p_fld% mfl_mgm_e  (:,:,jb),      & ! eastward flux            (out)
                        pmflux_w_bg   = p_fld% mfl_mgm_w  (:,:,jb),      & ! westward flux            (out)
                        pmflux_n_bg   = p_fld% mfl_mgm_n  (:,:,jb),      & ! northward flux           (out)
                        pmflux_s_bg   = p_fld% mfl_mgm_s  (:,:,jb),      & ! southward flux           (out)
                        kd           = kd(:,:),                          & ! k diagnostic             (out)
                        kd2          = kd2(:,:),                         & ! k^2 diagnostic           (out)
                        ld           = ld(:,:),                          & ! l diagnostic             (out)
                        ld2          = ld2(:,:),                         & ! l^2 diagnostic           (out)
                        md           = md(:,:),                          & ! m^2 diagnostic           (out)
                        md2          = md2(:,:),                         & ! m^2 diagnostic           (out)
                        cgz_diag     = cgz_diag(:,:),                    & ! cgz diagnostic           (out)
                        mB2          = mB2_save(:,:),                    & ! m^2B^2 diagnostic        (out)
                        mB2_aft      = mB2     (:,:)                     ) ! m^2B^2 diagnostic        (out)

      ELSE ! (lsat_mono)
        ! Monochromatic saturation scheme
        ! In this case saturation is calculated for each GW source 
        ! (background, convective) separately

        IF (nrays_add_bg(jg) /= 0) THEN
          ! Saturation for GWs from background sources
          CALL saturation_mono(nlev  = nlev,                            & ! no. of full levels       (in)
                          i_startidx = i_startidx,                      & ! first index of the block (in)
                          i_endidx   = i_endidx,                        & ! last index of the block  (in)
                          jray_start = 1,                               &
                          jray_end   = nrays(jg),                       &
                          dz         = p_metrics%ddqz_z_full(:,:,jb),   & ! full layer thickness     (in)
                          kray       = p_ray(jg)%k(:,:,1),              & ! horiz (lon) wavenumber   (in)
                          lray       = p_ray(jg)%l(:,:,1),              & ! horiz (lat) wavenumber   (in)
                          mray       = p_ray(jg)%m(:,:,1),              & ! vertical wavenumber      (in)
                          dens       = p_ray(jg)%wadens(:,:,1),         & ! wave action density      (in)
                          iexist     = p_ray(jg)%iexist(:,:,1),         & ! existence of ray         (in)
                          u_half     = u_half(:,:),                     & ! zonal wind (half levels) (in)
                          v_half     = v_half(:,:),                     & ! merid wind (half levels) (in)
                          bvf2       = bvf2(:,:),                       & ! Brunt-Väisala freq**2    (in)
                          gammash2   = gammash2(:,:),                   & ! inverse pinc scale height(in)
                          fc         = p_patch%cells%f_c(:,jb),         & ! Coriolis parameter       (in)
                          rho        = rho_half(:,:),                   & ! rho at half levels       (in)
                          kvisc      = kvisc(:,:),                      & ! Kinematic viscosity      (in)
                          lcalc_flux_4dir = lcalc_flux_4dir_bg(jg),     & ! flag to calc. 4-fluxes   (in)
                          uwflux     = p_fld% uwfl_mgm    (:,:,jb),     & ! uw mometum flux          (out)
                          vwflux     = p_fld% vwfl_mgm    (:,:,jb),     & ! vw mometum flux          (out)
                          energy     = p_fld% energy_mgm  (:,:,jb),     & ! wave energy              (out)
                          apmflux    = p_fld% apmfl_mgm    (:,:,jb),    & ! abs momentum flux        (out)
                          energy_p   = p_fld% energy_p_mgm(:,:,jb),     & ! wave potential energy    (out)
                          pmflux_e   = p_fld% mfl_mgm_e  (:,:,jb),      & ! eastward flux            (out)
                          pmflux_w   = p_fld% mfl_mgm_w  (:,:,jb),      & ! westward flux            (out)
                          pmflux_n   = p_fld% mfl_mgm_n  (:,:,jb),      & ! northward flux           (out)
                          pmflux_s   = p_fld% mfl_mgm_s  (:,:,jb)       ) ! southward flux           (out)
        ELSE
          p_fld% apmfl_mgm(:,:,jb) = 0._wp
          p_fld% uwfl_mgm(:,:,jb) = 0._wp   ;  p_fld% vwfl_mgm(:,:,jb) = 0._wp
          p_fld% mfl_mgm_e(:,:,jb) = 0._wp  ;  p_fld% mfl_mgm_w(:,:,jb) = 0._wp
          p_fld% mfl_mgm_n(:,:,jb) = 0._wp  ;  p_fld% mfl_mgm_s(:,:,jb) = 0._wp
          p_fld% energy_mgm(:,:,jb) = 0._wp ;  p_fld% energy_p_mgm(:,:,jb) = 0._wp
        END IF

      ENDIF  ! lsat_mono

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_saturation)

    ENDIF

    IF ( gws_conv_config%n_source(jg) > 0 ) THEN
      p_fld% uwfl_mgm(:,:,jb) = p_fld% uwfl_mgm(:,:,jb) + uwflux_cgw(:,:)
      p_fld% vwfl_mgm(:,:,jb) = p_fld% vwfl_mgm(:,:,jb) + vwflux_cgw(:,:)
    END IF

    ! Apply a sponge on fluxes 
    ! The other variables are untouched in the sponge layer, thus their output
    ! should be used with caution.
    DO jk = 1, nrdmax(jg)
      p_fld% uwfl_mgm  (:,jk,jb) = p_fld% uwfl_mgm  (:,jk,jb)*sponge_cos(jk)
      p_fld% vwfl_mgm  (:,jk,jb) = p_fld% vwfl_mgm  (:,jk,jb)*sponge_cos(jk)
      p_fld% energy_mgm(:,jk,jb) = p_fld% energy_mgm(:,jk,jb)*sponge_cos(jk)
    ENDDO

    ! Idealized case
    IF ((.NOT.ldynamics) .AND. (.NOT.lcoriolis)) p_fld% vwfl_mgm(:,:,jb) = 0._wp

    !============================ VTENDENCY ===============================
    ! Calculate horizontal wind tendencies based on the vertical divergence
    ! of pseudo-momentum fluxes. The total tendency and the tendency due to 
    ! convective GWs is calculated separately (the latter for diagnostics).
    !======================================================================

    IF (timers_level > 4) CALL timer_start(timer_msgwam_tendency)

    ! Calculate total wind tendency
    CALL tendency(nlev        = nlev,                               & ! no. of full levels       (in)
                  i_startidx  = i_startidx,                         & ! first index of the block (in)
                  i_endidx    = i_endidx,                           & ! last index of the block  (in)
                  dz          = p_metrics%ddqz_z_full(:,:,jb),      & ! full layer thickness     (in)
                  rho         = rho(:,:,jb),                        & ! rho at full levels       (in)
                  uwflux      = p_fld% uwfl_mgm(:,:,jb),            & ! uw mometum flux          (in)
                  vwflux      = p_fld% vwfl_mgm(:,:,jb),            & ! vw mometum flux          (in)
                  utend       = p_fld% ddt_u_gwd_mgm(:,:,jb),       & ! u tendency               (out)
                  vtend       = p_fld% ddt_v_gwd_mgm(:,:,jb)        ) ! v tendency               (out)

    ! Calculate wind tendency from convective GWs
    IF ( lcalc_cgw_tend(jg) ) THEN
      CALL tendency(nlev        = nlev,                             & ! no. of full levels       (in)
                    i_startidx  = i_startidx,                       & ! first index of the block (in)
                    i_endidx    = i_endidx,                         & ! last index of the block  (in)
                    dz          = p_metrics%ddqz_z_full(:,:,jb),    & ! full layer thickness     (in)
                    rho         = rho(:,:,jb),                      & ! rho at full levels       (in)
                    uwflux      = uwflux_cgw(:,:),                  & ! uw mometum flux          (in)
                    vwflux      = vwflux_cgw(:,:),                  & ! vw mometum flux          (in)
                    utend       = p_fld% gwd_conv_u(:,:,jb),        & ! u tendency               (out)
                    vtend       = p_fld% gwd_conv_v(:,:,jb)         ) ! v tendency               (out)
    END IF

    IF (timers_level > 4) CALL timer_stop(timer_msgwam_tendency)


    !======================== DIAG OUTPUT ==========================
    ! Output diagnostics to check whether saturation works correctly
    !===============================================================

    IF ( ldiagprof .AND. ndiag_msgwam > 0 .AND. jg == 1 ) THEN
    IF ( ANY(jb_diag(1:ndiag_msgwam) == jb) ) THEN

      IF (timers_level > 4) CALL timer_start(timer_msgwam_diagprof)

      CALL debug_messages_on

      DO idiag = 1,ndiag_msgwam

        IF ( jb /= jb_diag(idiag) ) CYCLE

          jc = jc_diag(idiag)

          IF (lsaturation) THEN
            DO jk = nlevp1,1,-1
              IF (mB2_save(jc,jk) > (alpha_sat*bvf2(jc,jk))**2) THEN
                CALL message('', 'lsaturation=.T. ==> Saturation Parametrized')
                WRITE(message_text,'(a,3E12.4)') 'alpha, N2, height:', &
                                    alpha_sat, bvf2(jc,jk), p_metrics%z_ifc(jc,jk,jb)
                CALL message('', TRIM(message_text))
                WRITE(message_text,'(a,3E12.4)') 'N^4*alpha^2, m^2B^2 ori, m^2B^2 reduced:', &
                                    (bvf2(jc,jk)*alpha_sat)**2, mB2_save(jc,jk), mB2(jc,jk)
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
              WRITE(message_text,'(a,i6,2E12.4)') 'jk, uwflux, vwflux:', &
                   jk, p_fld% uwfl_mgm(jc,jk,jb), p_fld% vwfl_mgm(jc,jk,jb)
              CALL message('', TRIM(message_text))
              WRITE(message_text,'(a,i6,2E12.4)') 'jk,   utend, vtend:', &
                   jk, p_fld% ddt_u_gwd_mgm(jc,jk,jb), p_fld% ddt_v_gwd_mgm(jc,jk,jb)
              CALL message('', TRIM(message_text))
            ENDDO
          ENDIF

      ENDDO  ! idiag

      CALL debug_messages_off

      IF (timers_level > 4) CALL timer_stop(timer_msgwam_diagprof)

    ENDIF ! ANY(jb_diag(1:ndiag_msgwam) == jb)
    ENDIF ! ldiagprof .AND. ndiag_msgwam > 0 .AND. jg == 1


  ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL

  ! Tendency limiter to stabilize high-top runs (not needed any more ...)
  IF (llimittend) THEN
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

    ENDDO ! jb
  ENDIF ! llimittend

  ! Idealized case
  IF ((.NOT.ldynamics) .AND. (.NOT.lcoriolis)) p_fld% ddt_v_gwd_mgm(:,:,:) = 0._wp

  ! Print out maximum of absolute of tendencies
  ! TODO: one would need a syncronization here in order to get full domain
  !       diagnostics. This is here only a preliminary printout to have a first 
  !       feeling how often the tendencies are limited via llimittend = .true.
  maxindexu = MAXLOC(ABS(p_fld% ddt_u_gwd_mgm))
  maxindexv = MAXLOC(ABS(p_fld% ddt_v_gwd_mgm))
  WRITE(message_text,'(a,E12.4,a,i4)') 'MAXABS of GW u-tendencies:', & 
  MAXVAL(ABS(p_fld% ddt_u_gwd_mgm)), ' at level:', maxindexu(2)
  CALL message('', TRIM(message_text))
  WRITE(message_text,'(a,E12.4,a,i4)') 'MAXABS of GW v-tendencies:', &   
  MAXVAL(ABS(p_fld% ddt_v_gwd_mgm)), ' at level:', maxindexv(2)
  CALL message('', TRIM(message_text))

  CALL message('gwdrag_msgwam_stst', 'Steady-state MS-GWaM-1D called')

END SUBROUTINE gwdrag_msgwam_stst
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE fieldsgrads(nlev,i_startidx,i_endidx,z,zhalf,dz,dzhalf,wgtfac_c,&
                       rho,temp,u,v,rho_half,temp_half,u_half,v_half,bvf2,kvisc,&
                       gammash2)
  INTEGER,     INTENT(IN)    :: nlev       ! number of full levels
  INTEGER,     INTENT(IN)    :: i_startidx ! first index of the block
  INTEGER,     INTENT(IN)    :: i_endidx   ! last index of the block
  REAL(wp),    INTENT(IN)    :: z(:,:)     ! height of full levels
  REAL(wp),    INTENT(IN)    :: zhalf(:,:) ! height of half levels
  REAL(wp),    INTENT(IN)    :: dz(:,:)    ! thickness of full levels
  REAL(vp),    INTENT(IN)    :: dzhalf(:,:)! thickness of half levels
  REAL(vp),    INTENT(IN)    :: wgtfac_c(:,:)           !
  REAL(wp),    INTENT(IN)    :: rho(:,:)                ! density at full levels
  REAL(wp),    INTENT(IN)    :: temp(:,:)               ! temperature at full levels
  REAL(wp),    INTENT(IN)    :: u(:,:)                  ! zonal wind (full levels)
  REAL(wp),    INTENT(IN)    :: v(:,:)                  ! meridional wind (full levels)
  REAL(wp),    INTENT(INOUT) :: rho_half(:,:)   ! density at half levels
  REAL(wp),    INTENT(INOUT) :: temp_half(:,:)  ! temperature at half levels
  REAL(wp),    INTENT(INOUT) :: u_half(:,:)     ! zonal wind
  REAL(wp),    INTENT(INOUT) :: v_half(:,:)     ! meridional wind
  REAL(wp),    INTENT(INOUT) :: bvf2(:,:)  ! Brunt-Väisälä freq**2
  REAL(wp),    INTENT(INOUT) :: kvisc(:,:) ! Kinematic viscosity
  REAL(wp),    INTENT(INOUT) :: gammash2(:,:) ! inverse pinc scale height**2

  INTEGER                    :: jk,jc ! vertical and block indices
  INTEGER                    :: nlevp1
  REAL(wp)                   :: wgtfac_c_jkm1           ! 1. - wgtfac_c
  REAL(wp)                   :: Hrho
  REAL(wp),    PARAMETER     :: grav_o_cpd = grav / cpd

  IF (msg_level >= 12) CALL message('fieldsgrads', 'MS-GWaM-ST: prepare resolved fields and grads')

  !----------------------------------------------------------------------
  ! Purpose:
  !         Calculate all necessary resolved fields appearing as input to 
  !         MS-GWaM
  ! Method:
  !         -- smooth wind and temperature
  !         -- calculate Brunt-Väisälä frequency
  !         -- calculate pinc scale height correction term gamma
  !----------------------------------------------------------------------

  nlevp1 = nlev+1

  bvf2(:,:) = 0._wp
  kvisc(:,:) = 0._wp
  gammash2(:,:) = 0._wp

  DO jk = 2, nlev
    DO jc = i_startidx, i_endidx
      wgtfac_c_jkm1 = 1._wp - wgtfac_c(jc,jk)
      rho_half (jc,jk) = wgtfac_c(jc,jk)*rho (jc,jk  ) &
        &                + wgtfac_c_jkm1*rho (jc,jk-1)
      temp_half(jc,jk) = wgtfac_c(jc,jk)*temp(jc,jk  ) &
        &                + wgtfac_c_jkm1*temp(jc,jk-1)
      u_half(jc,jk) = wgtfac_c(jc,jk)*u(jc,jk  ) &
        &             + wgtfac_c_jkm1*u(jc,jk-1)
      v_half(jc,jk) = wgtfac_c(jc,jk)*v(jc,jk  ) &
        &             + wgtfac_c_jkm1*v(jc,jk-1)
    ENDDO
  ENDDO
  DO jc = i_startidx, i_endidx
    rho_half (jc,1     ) = rho (jc,1   )
    rho_half (jc,nlevp1) = rho (jc,nlev)
    temp_half(jc,1     ) = temp(jc,1   )
    temp_half(jc,nlevp1) = temp(jc,nlev)
    u_half   (jc,1     ) = u   (jc,1   )
    u_half   (jc,nlevp1) = u   (jc,nlev)
    v_half   (jc,1     ) = v   (jc,1   )
    v_half   (jc,nlevp1) = v   (jc,nlev)
  ENDDO

  ! Smoothing over 2*nsmooth+1 points
  IF (lsmoothb) THEN
    CALL smooth_vert(nlevp1,i_startidx,i_endidx,nsmooth, temp_half)
    CALL smooth_vert(nlevp1,i_startidx,i_endidx,nsmooth, u_half   )
    CALL smooth_vert(nlevp1,i_startidx,i_endidx,nsmooth, v_half   )
  ENDIF

  DO jk = 2,nlev

    DO jc = i_startidx, i_endidx

      ! diag printout
      IF (msg_level >= 15) THEN
        WRITE(message_text,'(a,2i6,E12.4)') 'jc, jk, temp_half(jc,jk):', &
                                             jc, jk, temp_half(jc,jk)
        CALL message('', TRIM(message_text))
      ENDIF

      ! Calculate Brunt-Väisälä freq (half levels)
      ! negative pot temperature gradient is not allowed, 
      ! i.e. in that case bvf2 is set to bvf2_min
      bvf2(jc,jk) = MAX( bvf2_min, &
        &                grav/temp_half(jc,jk)*( grav_o_cpd +  &
        &                   (temp_half(jc,jk-1) - temp_half(jc,jk+1)) &
        &                   /(dz(jc,jk-1) + dz(jc,jk)) ) )

      ! calculate gamma (half levels): inverse pseudo incompressible scaleheight
      IF (lcorrlongwaves) THEN
        ! Density scaleheight assuming "locally" isothermal atmosphere 
        Hrho = rd*temp_half(jc,jk)/grav
        ! Assuming "locally" isothermal atmosphere 
        gammash2(jc,jk) = (facgamma/Hrho)**2
        ! facgamma is a namelist parameter with the following typical values:
        ! Pseudo-incompressible correction:         facgamma ~= 0.214
        ! Anelastic correction:                     facgamma  = 0.5
        ! Further decrease effects from long waves: facgamma  > 0.5
      ENDIF

      ! Diag printout
      IF (msg_level >= 15) THEN
        WRITE(message_text,'(a,2i6,2E12.4)') 'jc, jk, bvf2(jc,jk), 1/gammash2(jc,jk):', &
                                              jc, jk, bvf2(jc,jk), 1._wp/gammash2(jc,jk)
        CALL message('', TRIM(message_text))
      ENDIF

      ! Calculate kinematic viscosity profile
      kvisc(jc,jk) = dyn_visc_sutherland(temp_half(jc,jk))/rho_half(jc,jk)

    ENDDO ! jc

  ENDDO ! jk
  bvf2(:,1) = bvf2(:,2)
  bvf2(:,nlevp1) = bvf2(:,nlev)
  gammash2(:,1) = gammash2(:,2)
  gammash2(:,nlevp1) = gammash2(:,nlev)
  kvisc(:,1) = kvisc(:,2)
  kvisc(:,nlevp1) = kvisc(:,nlev)

END SUBROUTINE fieldsgrads
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE init_gw_orretal(nlev,i_startidx,i_endidx,jray_offset,mdatetime,&
                           kray,lray,mray,dens,iexist,&
                           jklaunch,lat,rhol,bvfl2,lat_prof_bw,lat_prof_bs,lat_prof,prec)
  TYPE(datetime),INTENT(IN)    :: mdatetime
  INTEGER,       INTENT(IN)    :: nlev                 ! no. of model levels
  INTEGER,       INTENT(IN)    :: i_startidx           ! first index of the block
  INTEGER,       INTENT(IN)    :: i_endidx             ! last index of the block
  INTEGER,       INTENT(IN)    :: jray_offset
  INTEGER,       INTENT(INOUT) :: iexist(:,:)
  INTEGER,       INTENT(IN)    :: jklaunch             ! launch level index
  REAL(wp),      INTENT(INOUT) :: mray(:,:)
  REAL(wp),      INTENT(INOUT) :: kray(:,:)
  REAL(wp),      INTENT(INOUT) :: lray(:,:)
  REAL(wp),      INTENT(INOUT) :: dens(:,:)
  REAL(wp),      INTENT(IN)    :: lat(:)               ! Latitude
  REAL(wp),      INTENT(IN)    :: rhol(:)              ! density at launch level
  REAL(wp),      INTENT(IN)    :: bvfl2(:)             ! B-V freq**2 at launch level
  REAL(wp),      INTENT(IN)    :: lat_prof_bw(:)       ! latitudinal factor
  REAL(wp),      INTENT(IN)    :: lat_prof_bs(:)       ! latitudinal factor
  REAL(wp),      INTENT(INOUT) :: lat_prof(:)          ! latitudinal factor
  REAL(wp),      INTENT(IN)    :: prec(:)              ! Total precip

  INTEGER  :: jc
  INTEGER  :: jray
  INTEGER  :: khdim    ! dimension of total horizontal wavenumber spectral array

  INTEGER  :: ic, im, ikh, iomega          ! indices for spectral elements
  INTEGER  :: iazi                         ! index of the azimuth angle
  INTEGER  :: ispec                        ! index for spectral location
  INTEGER  :: cdim                         ! dimension of phase speed spectral array
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
  REAL(wp) :: gauss, fluxlaun
  REAL(wp) :: tfact
  INTEGER  :: year_m4k
  REAL(wp), PARAMETER :: tphase_0 = 8520._wp/8760._wp  ! phase for 00 UTC 22 Dec. (for non-leap years)
  REAL(wp), PARAMETER :: ndaysoffset(12) = REAL(  &
    &                      (/0,31,59,90,120,151,181,212,243,273,304,334/), wp )

  IF (msg_level >= 12) CALL message('init_gw_orretal', &
    &                   'MS-GWaM-ST: initialize Desaubies type background GW field')

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

  ! Calculate time dependent launch flux magnitudes. Currently, only the Gregorian 
  ! calendar with leap years is considered.
  tfact = ( ndaysoffset(mdatetime%date%month) + REAL(mdatetime%date%day-1,wp) )*24._wp  &
    &     + REAL(mdatetime%time%hour,wp) + REAL(mdatetime%time%minute,wp)/60._wp        &
    &     + ( REAL(mdatetime%time%second,wp)                                            &
    &         + 0.001_wp*REAL(mdatetime%time%ms,wp) )/3600._wp
  ! to convert the type to normal integer with a moderate size
  year_m4k = INT( MOD(mdatetime%date%year, INT(4000, KIND=KIND(mdatetime%date%year))) )
  IF ( ( MOD(year_m4k,4) == 0 .AND. MOD(year_m4k,100) /= 0 )  &
    &  .OR. MOD(year_m4k,400) == 0 ) THEN   ! leap year
    IF (mdatetime%date%month > 2)  tfact = tfact + 24._wp
    tfact = tfact/8784._wp   ! hours for 366 days
  ELSE
    tfact = tfact/8760._wp   ! hours for 365 days
  END IF
  tfact = 0.5_wp*(COS((tfact - tphase_0)*(2._wp*pi)) + 1._wp)

  DO jc = i_startidx, i_endidx
    lat_prof(jc) = lat_prof_bs(jc) + (lat_prof_bw(jc)-lat_prof_bs(jc))*tfact
  ENDDO

  ! ------------------------------------------------------------------
  ! PART 1: calculate launch spectrum in (c,phi) space (based on mo_gwd_wms.f90)
  ! ------------------------------------------------------------------

  ! Set same dimension for m and c spectral elements, we suppose keeping cdim 
  ! (and not simply using mdim everywhere) keeps the code more understandable
  cdim = mdim

  ! mstar
  mstar = 2._wp*pi/2000._wp

  ! Unstreched phase velocity spectrum
  dcin(1) = (cimax-cimin)/REAL(cdim-1,wp)
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

    IF (lat_prof(jc) == 0._wp) THEN
      dens(jc,jray_offset+1:jray_offset+nrays_add_bg(jg)) = 0.0_wp
      CYCLE
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

    ! Calculate rho*F(m,omega,phi)
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
    jray  = jray_offset
    pu(jklaunch,:) = 0._wp
    DO iazi=1,iazidim ! loop over azimuth angles
      DO ikh=1,khdim  ! loop over horizontal wavenumber
        DO im=1,mdim  ! loop over vertical wavenumber

            ! Index for omega elements (same as for k_h elements)
            iomega = ikh

            ! Ray volume counter
            jray = jray + 1
 
            ! Eastward
            IF (7._wp*pi/4._wp < phi(iazi) .OR. phi(iazi) < 1._wp/4._wp*pi) THEN
              kray(jc,jray)  = kh(im,iomega)
              lray(jc,jray)  = 0._wp
            ! Northward
            ELSEIF (pi/4._wp < phi(iazi) .AND. phi(iazi) < 3._wp/4._wp*pi) THEN
              kray(jc,jray)  = 0._wp
              lray(jc,jray)  = kh(im,iomega)
            ! Westward
            ELSEIF (3._wp*pi/4._wp < phi(iazi) .AND. phi(iazi) < 5._wp/4._wp*pi) THEN
              kray(jc,jray)  = -kh(im,iomega)
              lray(jc,jray)  = 0._wp
            ! Southward
            ELSEIF (5._wp*pi/4._wp < phi(iazi) .AND. phi(iazi) < 7._wp/4._wp*pi) THEN
              kray(jc,jray)  = 0._wp
              lray(jc,jray)  = -kh(im,iomega)
            ENDIF

            ! Calculate m at ray position 
            mray(jc,jray) = m(im)

            ! cgz at ray position based on hydrostatic disp. rel. 
            cgz = -omega_l(iomega)/m(im)

            ! Flux integrated over the spectral segment
            ! rho*F(k_h,m,phi)*dkh*dphi*dm
            !   (but dphi has already been included in flux3)
            fluxray = flux3(ikh,im,iazi)*dkh(im,iomega)*dm(im)

            ! Define wave action density
            ! N(k,l,m) = rho*F(k,l,m)/k_h/cgz
            dens(jc,jray) = fluxray/kh(im,iomega)/cgz

            iexist(jc,jray) = jklaunch   ! closest half level

            IF (msg_level >= 15) THEN
                WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, kh, L_h, dkh, L_dkh',&
                jc, jray, ikh, iazi, im, kh(im,iomega), 2._wp*pi/kh(im,iomega), dkh(im,iomega), 2._wp*pi/dkh(im,iomega)
                CALL message('', TRIM(message_text))
                IF (kray(jc,jray) /= 0._wp) THEN
                  WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, kray(jc,jray), L_x',&
                           jc, jray, ikh, iazi, im, kray(jc,jray), 2._wp*pi/kray(jc,jray)
                  CALL message('', TRIM(message_text))
                ELSE
                  WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, kray(jc,jray), L_x',&
                           jc, jray, ikh, iazi, im, kray(jc,jray), 0._wp
                  CALL message('', TRIM(message_text))
                ENDIF
                IF (lray(jc,jray) /= 0._wp) THEN
                  WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, lray(jc,jray), L_y',&
                           jc, jray, ikh, iazi, im, lray(jc,jray), 2._wp*pi/lray(jc,jray)
                  CALL message('', TRIM(message_text))
                ELSE
                  WRITE(message_text,'(a,5i6,4E12.4)') 'jc, jray, ikh, iazi, im, lray(jc,jray), L_y',&
                           jc, jray, ikh, iazi, im, lray(jc,jray), 0._wp
                  CALL message('', TRIM(message_text))
                ENDIF
                WRITE(message_text,'(a,5i6,5E12.4)') 'jc, jray, ikh, iazi, im, cin, mray(jc,jray), L_z',&
                         jc, jray, ikh, iazi, im, cin(im), mray(jc,jray), 2._wp*pi/mray(jc,jray)
                CALL message('', TRIM(message_text))
                WRITE(message_text,'(a,5i6,2E12.4)') 'jc, jray, ikh, iazi, im',&
                                                      jc, jray, ikh, iazi, im
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
                                *cgz*dens(jc,jray)

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
      WRITE(message_text,'(a)') 'Init positions of new rays:'
      CALL message('', TRIM(message_text))
        DO ispec = 1, nrays_add_bg(jg)
          jray = jray_offset + ispec
          WRITE(message_text,'(a,i8,5E12.4)') 'jray, mray, wadens:', &
            & jray, mray(jc,jray), dens(jc,jray)
          CALL message('init_gw_orretal', TRIM(message_text))
        ENDDO
      WRITE(message_text,'(a,i6)') 'Total no. of rays:', nrays_bg(jg)
      CALL message('', TRIM(message_text))
    ENDIF

  ENDDO ! jc

END SUBROUTINE init_gw_orretal
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE init_gw_conv(jb,i_startidx,i_endidx,jray_offset,flag_cgw)

  INTEGER, INTENT(IN) :: jb
  INTEGER, INTENT(IN) :: i_startidx, i_endidx
  INTEGER, INTENT(IN) :: jray_offset
  INTEGER, INTENT(IN) :: flag_cgw(:)

  INTEGER  :: nspec
  INTEGER  :: jc, ispec
  INTEGER  :: jray

  IF (msg_level >= 12)  CALL message('init_gw_conv',  &
    &  'MS-GWaM-ST: initialize GW field from convection')

  nspec = nrays_add_cv(jg)

  ! Note: the launch mom. flux spectrum computation is based on the assumption
  ! that the POSITIVE FREQUENCY BRANCH IS USED, which also means that a NEGATIVE
  ! VERTICAL WAVENUMBER spectrum will be used. This results in upward propagating 
  ! gravity waves (cgz > 0) injected at launch level.

  DO jc = i_startidx, i_endidx

    IF ( flag_cgw(jc) < 0 ) THEN
      p_ray(jg)%wadens(jc,jray_offset+1:jray_offset+nspec,1) = 0.0_wp
      CYCLE
    END IF

    !---------------------------------------------------------------------------
    ! Launch new rays
    !---------------------------------------------------------------------------

    DO ispec = 1, nspec

      jray = jray_offset + ispec
 
      p_ray(jg)%iexist(jc,jray,1) = p_ray_conv(jg)%jk_source(jc,jb)

      p_ray(jg)%k(jc,jray,1) = p_ray_conv(jg)%k(jc,ispec,jb)
      p_ray(jg)%l(jc,jray,1) = p_ray_conv(jg)%l(jc,ispec,jb)
      p_ray(jg)%m(jc,jray,1) = p_ray_conv(jg)%m(jc,ispec,jb)

      p_ray(jg)%wadens(jc,jray,1) = p_ray_conv(jg)%wadens(jc,ispec,jb)

    ENDDO  ! ispec

  ENDDO  ! jc

END SUBROUTINE init_gw_conv
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE saturation(nlev,i_startidx,i_endidx,jray_start_cv,jray_end_cv,jray_start_bg,jray_end_bg,&
                      dz,kray,lray,mray,dens,iexist,&
                      u_half,v_half,bvf2,gammash2,fc,rho,kvisc,flag_cgw,lcalc_flux_4dir_cv,lcalc_flux_4dir_bg,&
                      uwflux_cv,vwflux_cv,&
                      uwflux_bg,vwflux_bg,energy_bg,apmflux_bg,energy_p_bg,pmflux_e_bg,pmflux_w_bg,pmflux_n_bg,pmflux_s_bg,&
                      kd,kd2,ld,ld2,md,md2,cgz_diag,mB2,mB2_aft)
  INTEGER,       INTENT(IN)    :: nlev
  INTEGER,       INTENT(IN)    :: i_startidx ! first index of the block
  INTEGER,       INTENT(IN)    :: i_endidx   ! last index of the block
  INTEGER,       INTENT(IN)    :: jray_start_cv
  INTEGER,       INTENT(IN)    :: jray_end_cv
  INTEGER,       INTENT(IN)    :: jray_start_bg
  INTEGER,       INTENT(IN)    :: jray_end_bg
  REAL(wp),      INTENT(IN)    :: dz(:,:)
  INTEGER,       INTENT(IN)    :: iexist(:,:)
  REAL(wp),      INTENT(IN)    :: mray(:,:)
  REAL(wp),      INTENT(IN)    :: kray(:,:)
  REAL(wp),      INTENT(IN)    :: lray(:,:)
  REAL(wp),      INTENT(IN)    :: dens(:,:)
  REAL(wp),      INTENT(IN)    :: u_half(:,:), v_half(:,:)
  REAL(wp),      INTENT(IN)    :: bvf2(:,:), rho(:,:), gammash2(:,:)
  REAL(wp),      INTENT(IN)    :: kvisc(:,:)      ! Kinematic viscosity
  REAL(wp),      INTENT(IN)    :: fc(:)           ! Coriolis parameter
  INTEGER,       INTENT(IN)    :: flag_cgw(:)
  LOGICAL,       INTENT(IN)    :: lcalc_flux_4dir_cv
  LOGICAL,       INTENT(IN)    :: lcalc_flux_4dir_bg
  REAL(wp),      INTENT(OUT)   :: uwflux_cv(:,:)
  REAL(wp),      INTENT(OUT)   :: vwflux_cv(:,:)
  REAL(wp),      INTENT(OUT)   :: uwflux_bg(:,:)
  REAL(wp),      INTENT(OUT)   :: vwflux_bg(:,:)
  REAL(wp),      INTENT(OUT)   :: energy_bg(:,:)
  REAL(wp),      INTENT(OUT)   :: apmflux_bg(:,:)
  REAL(wp),      INTENT(OUT)   :: energy_p_bg(:,:)
  REAL(wp),      INTENT(OUT)   :: pmflux_e_bg(:,:), pmflux_w_bg(:,:)
  REAL(wp),      INTENT(OUT)   :: pmflux_n_bg(:,:), pmflux_s_bg(:,:)
  REAL(wp),      INTENT(OUT)   :: mB2(:,:), mB2_aft(:,:)
  REAL(wp),      INTENT(OUT)   :: kd(:,:), ld(:,:), md(:,:)
  REAL(wp),      INTENT(OUT)   :: kd2(:,:), ld2(:,:), md2(:,:), cgz_diag(:,:)
  INTEGER,       PARAMETER   :: max_iter_sat = 1000
  ! Counters, indices
  INTEGER                    :: nlevp1
  INTEGER                    :: jk, jc
  INTEGER                    :: jkmax, jkl_cv, jkl_bg
  INTEGER                    :: jray, jraymin, jraymax
  INTEGER                    :: iter
  ! Other subquantities
  REAL(wp)                   :: mB2_bfr
  REAL(wp)                   :: mB2K2
  REAL(wp)                   :: counter
  REAL(wp)                   :: al2bvf4, tmp
  REAL(wp)                   :: crit
  REAL(wp)                   :: fc2, abs_fc, m2_max
  REAL(wp)                   :: f2_o_bvf2
  REAL(wp)                   :: two_diffusion
  REAL(wp)                   :: omega2
  REAL(wp)                   :: minus_two_kvisc_dz
  REAL(wp)                   :: K2_p_gam2
  REAL(wp)                   :: weight(nrays(jg))
  REAL(wp)                   :: Kh2(nrays(jg)), Kh(nrays(jg))
  REAL(wp)                   :: omega_gr(nrays(jg))
  REAL(wp)                   :: omegaint(nrays(jg))
  REAL(wp)                   :: m2(nrays(jg)), K2(nrays(jg))
  REAL(wp)                   :: cgz(nrays(jg)), inv_cgz(nrays(jg))
  REAL(wp)                   :: ratio(nrays(jg))
  REAL(wp)                   :: integral1(nrays(jg))
  REAL(wp)                   :: integral2(nrays(jg))
  REAL(wp)                   :: wadflxprof (nrays(jg),nlev+1)
  REAL(wp)                   :: wadensprof (nrays(jg),nlev+1)
  REAL(wp)                   :: intgcgz_k  (nrays(jg),nlev+1)
  REAL(wp)                   :: intgcgz_l  (nrays(jg),nlev+1)
  REAL(wp)                   :: intgcgz_h  (nrays(jg),nlev+1)
  REAL(wp)                   :: intgomega  (nrays(jg),nlev+1)
  REAL(wp)                   :: intgomega_p(nrays(jg),nlev+1)

  IF (msg_level >= 12) CALL message('saturation', 'MS-GWaM-ST: wave breaking parametrization')

  !----------------------------------------------------------------------
  ! Purpose: 
  !         Parameterize GW breaking (saturation), i.e. the corresponding 
  !         loss in GW energy and momentum fluxes.         
  !
  ! Method:
  !         The method is based on Lindzen's idea that GWs break when 
  !         they turn the potential temperature gradient to statically 
  !         unstable (d\theta/dz < 0).
  !         -- calculate equilibrium profile
  !         -- calculate effects of critical layer filtering on momflux
  !         -- calculate effects of reflection layers on momflux
  !         -- diagnose static instbility by \int dm m^2B^2 >= N^4
  !         -- calculate turbulent viscosity (in case of instability)
  !         -- calculate effects of wave breaking on momflux
  !         -- ...
  !         -- see Bölöni et al. (2020)
  !----------------------------------------------------------------------

  ratio(:) = 1._wp

  nlevp1 = nlev+1

  m2_max = (2.0_wp*pi/lambda_z_min)**2

  ! Initialization (nproma,nlevp1)

  IF ( gws_conv_config%n_source(jg) > 0 ) THEN
    uwflux_cv(:,:) = 0._wp   ;  vwflux_cv(:,:) = 0._wp
  END IF

  IF (nrays_add_bg(jg) /= 0) THEN
    uwflux_bg(:,:) = 0._wp   ;  vwflux_bg(:,:) = 0._wp  ;  energy_bg(:,:) = 0._wp
    apmflux_bg(:,:) = 0._wp  ;  energy_p_bg(:,:) = 0._wp
    IF ( lcalc_flux_4dir_bg ) THEN
      pmflux_e_bg(:,:) = 0._wp  ;  pmflux_w_bg(:,:) = 0._wp
      pmflux_n_bg(:,:) = 0._wp  ;  pmflux_s_bg(:,:) = 0._wp
    END IF
  END IF

  mB2(:,:) = 0._wp

  IF ( ldiagprof ) THEN
    mB2_aft(:,:) = 0._wp
    kd(:,:)      = 0._wp
    kd2(:,:)     = 0._wp
    ld(:,:)      = 0._wp
    ld2(:,:)     = 0._wp
    md(:,:)      = 0._wp
    md2(:,:)     = 0._wp
    cgz_diag(:,:)= 0._wp
  END IF
 
  ! Outer loop (jc) start

  DO jc = i_startidx, i_endidx

    ! Initialization (nrays(jg),nlevp1)
    wadflxprof (:,:) = 0._wp
    wadensprof (:,:) = 0._wp
    intgcgz_h  (:,:) = 0._wp
    intgomega  (:,:) = 0._wp
    intgomega_p(:,:) = 0._wp

    fc2 = fc(jc)**2
    abs_fc = ABS(fc(jc))

    !=====  SOURCE LEVEL  ======================================================

    ! Calculate the ground-based omega and vertical flux of wave action density
    !   : omega_gr(jray) and wadflxprof(jray,jk_source)

    DO jray = 1, nrays(jg)

      IF (dens(jc,jray) == 0._wp) THEN
        Kh2(jray) = 1._wp   ! non-zero value, as it is used as a denominator
        CYCLE
      END IF

      !K_h^2=k^2+l^2
      Kh2(jray) = kray(jc,jray)**2 + lray(jc,jray)**2
      Kh(jray) = SQRT(Kh2(jray))

      ! at the source level

      jk = iexist(jc,jray)   !  iexist /= 0  IF dens(jc,jray) /= 0._wp

      !K^2=k^2+l^2+m^2
      m2(jray) = mray(jc,jray)**2
      K2(jray) = Kh2(jray) + m2(jray)

      K2_p_gam2 = K2(jray) + gammash2(jc,jk)

      ! Calculate intrinsic frequency at launch level
      omega2 = (bvf2(jc,jk)*Kh2(jray)+fc2*(m2(jray)+gammash2(jc,jk)))/K2_p_gam2
      omegaint(jray) = branch*SQRT(omega2)

      ! Calculate the Doppler shift term at launch level
      omega_gr(jray) = omegaint(jray)  &
        &  + ( kray(jc,jray)*u_half(jc,jk) + lray(jc,jray)*v_half(jc,jk) )

      ! A*cgz
      wadflxprof(jray,jk) = ABS( dens(jc,jray)*( mray(jc,jray)*(omega2 - fc2)  &
        &                        /(omegaint(jray)*K2_p_gam2) ) )

    ENDDO

    !=====  ABOVE THE SOURCE  ==================================================

    ! Calculate the vertical profile of the flux of wave action density and
    !   various variables for diagnostics, given the conserved omega_gr
    !   : wadflxprof(jray,jk), wadensprof(jray,jk), ...

    jkmax = 0

    jkl_bg = 0  ;  jkl_cv = 0
    IF (nrays_add_bg(jg) /= 0) THEN
      jkl_bg = iexist(jc,jray_start_bg)
      jkmax = MAX(jkmax,jkl_bg)
    END IF
    IF (flag_cgw(jc) > 0) THEN
      jkl_cv = iexist(jc,jray_start_cv)
      jkmax = MAX(jkmax,jkl_cv)
    END IF

    jraymin = nrays(jg)+1  ;  jraymax = 0

    VERTLOOP:  DO jk = jkmax, 2, -1

      f2_o_bvf2 = fc2/bvf2(jc,jk)
      al2bvf4 = (alpha_sat*bvf2(jc,jk))**2
      tmp = 2._wp*bvf2(jc,jk)**2/rho(jc,jk)

      mB2_bfr = 0._wp

      ! Initialization (nrays(jg))
      integral1(:) = 0._wp  ;  integral2(:) = 0._wp
      K2(:) = 0._wp
      cgz(:) = 0._wp  ;  inv_cgz(:) = 0._wp

      IF (jk == jkl_cv) THEN
        jraymin = MIN(jraymin,jray_start_cv)
        jraymax = MAX(jraymax,jray_end_cv)
      END IF
      IF (jk == jkl_bg) THEN
        jraymin = MIN(jraymin,jray_start_bg)
        jraymax = MAX(jraymax,jray_end_bg)
      END IF

      ! jraymin <= jraymax always at jk <= jkmax (at least one source exists)

      DO jray = jraymin, jraymax

        ! Skip if
        !   1) this spectral element is empty (no corresponding source), or
        !   2) it is yet below the source level
        IF (wadflxprof(jray,jk) == 0._wp) CYCLE

        ! Calculate the Doppler shift term profile 
        omegaint(jray) = omega_gr(jray)  &
          &  - ( kray(jc,jray)*u_half(jc,jk) + lray(jc,jray)*v_half(jc,jk) )

        omega2 = omegaint(jray)**2
        m2(jray) = (bvf2(jc,jk)-omega2)*Kh2(jray)/(omega2-fc2) - gammash2(jc,jk)

        !=====  FILTERING  =====================================================

        ! 1) critical level if:   omegaint^2-f^2 <= 0 or sign change of omegaint
        !                         or too small Lz
        !    --> remove the spectral element above the critical level
        IF ( omegaint(jray)*branch <= abs_fc .OR. m2(jray) > m2_max ) THEN
          wadensprof(jray,jk) = 0._wp
          wadflxprof(jray,jk) = 0._wp
          CYCLE  ! jray
        ENDIF

        ! 2) reflection level if:   m^2 <= 0
        !    --> remove the complete spectral element as in a steady state
        !        picture a reflected (downward) wave carries exactly the opposite mom
        !        flux as the wave before reflection (upward) so they cancel each other
        IF (m2(jray) <= 0._wp) THEN
          wadflxprof(jray,jk+1:) = wadflxprof(jray,jk+1:) - wadflxprof(jray,jk)
          wadensprof(jray,jk) = 0._wp
          wadflxprof(jray,jk) = 0._wp
          CYCLE  ! jray
        ENDIF

        !K^2_low=k^2+l^2+m^2
        K2(jray) = Kh2(jray) + m2(jray)

        cgz(jray) = SQRT(m2(jray))*(omega2 - fc2)/(omegaint(jray)*(K2(jray)+gammash2(jc,jk)))
        inv_cgz(jray) = 1._wp/cgz(jray)

        wadensprof(jray,jk) = wadflxprof(jray,jk)*inv_cgz(jray)

        integral2(jray) = m2(jray)*Kh2(jray)/omegaint(jray)

        integral1(jray) = integral2(jray)/K2(jray)

!-----------------------------------------------------------
        ! comment out for dt-based formulation
        ratio(jray) = ABS(inv_cgz(jray))   ! dz-based formulation
!-----------------------------------------------------------

      ENDDO  ! jray

      !=====  SATURATION  ======================================================

      IF (ldiagprof)  mB2_bfr =  &
        &  SUM(wadensprof(jraymin:jraymax,jk)*integral1(jraymin:jraymax))*tmp

      L_ITER:  DO iter = 1, MIN(jraymax-jraymin+1, max_iter_sat)

        ! Calculate m^2*B^2
        mB2(jc,jk) = SUM(wadensprof(jraymin:jraymax,jk)  &
          &              *integral1(jraymin:jraymax))*tmp

        crit = mB2(jc,jk) - al2bvf4

        IF (crit < 0._wp)  EXIT  ! L_ITER
!       IF (crit < al2bvf4*1.e-3)  EXIT   ! allow 0.1% errors to reduce iteration

        ! Calculate K^2*m^2*B^2 (if ratio = 1) or K^2*m^2*B^2/cgz
        mB2K2 = SUM(wadensprof(jraymin:jraymax,jk)*integral2(jraymin:jraymax)  &
          &             *ratio(jraymin:jraymax))*tmp

        ! Diffusion : diffusivity * dt (dz)
        two_diffusion = crit/mB2K2    ! 2*(0.5*crit/mB2K2)

        weight(jraymin:jraymax) = 999._wp

        WHERE ( wadflxprof(jraymin:jraymax,jk) /= 0._wp )
          ! Calculate weight for reducing wave action density (turbulent viscosity)
          weight(jraymin:jraymax) = MAX(0._wp, 1._wp - two_diffusion  &
            &                       *K2(jraymin:jraymax)*ratio(jraymin:jraymax))

          wadensprof(jraymin:jraymax,jk) = wadensprof(jraymin:jraymax,jk)  &
            &                                 *weight(jraymin:jraymax)
          wadflxprof(jraymin:jraymax,jk) = wadflxprof(jraymin:jraymax,jk)  &
            &                                 *weight(jraymin:jraymax)
        END WHERE

        IF ( .NOT. ANY( weight(jraymin:jraymax) == 0._wp ) )  EXIT  ! L_ITER

      ENDDO L_ITER

      !=====  CALCFLUX 1 : absolute variables  =================================

      ! The absolute variables are calculated here before 'wadflx' is changed
      !   when the wave reflection occurs somewhere above this level.
      !   : absolute momentum flux and energy

      WHERE ( wadflxprof(jraymin:jraymax,jk) /= 0._wp )
        intgcgz_h(jraymin:jraymax,jk) = Kh(jraymin:jraymax)  &
          &                             *wadflxprof(jraymin:jraymax,jk)
        intgomega(jraymin:jraymax,jk) = omegaint(jraymin:jraymax)  &
          &                             *wadensprof(jraymin:jraymax,jk)
        intgomega_p(jraymin:jraymax,jk) = intgomega(jraymin:jraymax,jk)*0.5_wp  &
          &  /( 1._wp + f2_o_bvf2*(m2(jraymin:jraymax)+gammash2(jc,jk))         &
          &             /Kh2(jraymin:jraymax) )
      END WHERE

      !-------------------------------------------------------------------------

      IF (ldiagprof) THEN
        counter = REAL(COUNT(wadflxprof(jraymin:jraymax,jk) /= 0._wp),wp)
        IF (counter /= 0._wp) THEN
          DO jray = jraymin, jraymax
            IF (wadflxprof(jray,jk) == 0._wp) CYCLE

            mB2_aft(jc,jk) = mB2_aft(jc,jk) + wadensprof(jray,jk)*integral1(jray)
            kd(jc,jk) = kd(jc,jk) + kray(jc,jray)
            kd2(jc,jk) = kd2(jc,jk) + kray(jc,jray)**2
            ld(jc,jk) = ld(jc,jk) + lray(jc,jray)
            ld2(jc,jk) = ld2(jc,jk) + lray(jc,jray)**2
            md(jc,jk) = md(jc,jk) - SQRT(m2(jray))
            md2(jc,jk) = md2(jc,jk) + m2(jray)
            cgz_diag(jc,jk) = cgz_diag(jc,jk) + cgz(jray)
          ENDDO
          mB2_aft(jc,jk) = mB2_aft(jc,jk)*tmp
          kd(jc,jk)       = kd(jc,jk)/counter
          kd2(jc,jk)      = kd2(jc,jk)/counter
          ld(jc,jk)       = ld(jc,jk)/counter
          ld2(jc,jk)      = ld2(jc,jk)/counter
          md(jc,jk)       = md(jc,jk)/counter
          md2(jc,jk)      = md2(jc,jk)/counter
          cgz_diag(jc,jk) = cgz_diag(jc,jk)/counter
        ENDIF
        mB2(jc,jk) = mB2_bfr
      ENDIF

      ! Pass the flux value to the above level.
      IF ( lmvisc ) THEN
        minus_two_kvisc_dz = -(kvisc(jc,jk-1) + kvisc(jc,jk))*dz(jc,jk-1)
        DO jray = jraymin, jraymax
          wadflxprof(jray,jk-1) = wadflxprof(jray,jk)  &
            &                    *EXP(minus_two_kvisc_dz*inv_cgz(jray)*K2(jray))
        ENDDO
      ELSE
        wadflxprof(:,jk-1) = wadflxprof(:,jk)
      END IF

    ENDDO VERTLOOP


    !=====  CALCFLUX 2 : others  ===============================================

    ! The fluxes are calculated after vanishing 'wadflx' for reflected waves.

    IF (flag_cgw(jc) > 0) THEN
      DO jk = 2, nlev
        intgcgz_k(jray_start_cv:jray_end_cv,jk)           &
          &  =        kray(jc,jray_start_cv:jray_end_cv)  &
          &    *wadflxprof(jray_start_cv:jray_end_cv,jk)
        intgcgz_l(jray_start_cv:jray_end_cv,jk)           &
          &  =        lray(jc,jray_start_cv:jray_end_cv)  &
          &    *wadflxprof(jray_start_cv:jray_end_cv,jk)
        uwflux_cv(jc,jk) = SUM(intgcgz_k(jray_start_cv:jray_end_cv,jk))
        vwflux_cv(jc,jk) = SUM(intgcgz_l(jray_start_cv:jray_end_cv,jk))
      ENDDO 

      ! Upper boundary : no transmission to conserve momentum
      uwflux_cv(jc,1) = 0._wp  ;  vwflux_cv(jc,1) = 0._wp
      ! Lower boundary : vertical flux = 0
      uwflux_cv(jc,nlevp1) = 0._wp  ;  vwflux_cv(jc,nlevp1) = 0._wp
    END IF  ! flag_cgw(jc)

    IF (nrays_add_bg(jg) /= 0) THEN
      DO jk = 2, nlev
        intgcgz_k(jray_start_bg:jray_end_bg,jk)           &
          &  =        kray(jc,jray_start_bg:jray_end_bg)  &
          &    *wadflxprof(jray_start_bg:jray_end_bg,jk)
        intgcgz_l(jray_start_bg:jray_end_bg,jk)           &
          &  =        lray(jc,jray_start_bg:jray_end_bg)  &
          &    *wadflxprof(jray_start_bg:jray_end_bg,jk)
        uwflux_bg(jc,jk) = SUM(intgcgz_k(jray_start_bg:jray_end_bg,jk))
        vwflux_bg(jc,jk) = SUM(intgcgz_l(jray_start_bg:jray_end_bg,jk))
        energy_bg(jc,jk) = SUM(intgomega(jray_start_bg:jray_end_bg,jk))
      ENDDO 
      DO jk = 2, nlev
        apmflux_bg (jc,jk) = SUM(intgcgz_h  (jray_start_bg:jray_end_bg,jk))
        energy_p_bg(jc,jk) = SUM(intgomega_p(jray_start_bg:jray_end_bg,jk))
      ENDDO 

      ! Upper boundary : no transmission to conserve momentum
      uwflux_bg(jc,1) = 0._wp  ;  vwflux_bg(jc,1) = 0._wp
      ! Lower boundary : vertical flux = 0
      uwflux_bg(jc,nlevp1) = 0._wp  ;  vwflux_bg(jc,nlevp1) = 0._wp

      ! Calculate energy for 4 dierctions separately
      IF ( lcalc_flux_4dir_bg ) THEN
        DO jray = jray_start_bg, jray_end_bg
          IF (kray(jc,jray) > 0._wp) THEN
            DO jk = 2, nlev
              pmflux_e_bg(jc,jk) = pmflux_e_bg(jc,jk) + intgcgz_k(jray,jk)
            ENDDO
          END IF
          IF (lray(jc,jray) > 0._wp) THEN
            DO jk = 2, nlev
              pmflux_n_bg(jc,jk) = pmflux_n_bg(jc,jk) + intgcgz_l(jray,jk)
            ENDDO
          END IF
        ENDDO
      END IF
    END IF  ! nrays_add_bg(jg)

    ! The other diagnostic outputs which are not directly involved in the tendency
    !   calculation will have zero values at jk = 1 and nlevp1, as initialized.
    !   : absolute flux and fluxes at each direction

  ENDDO ! jc

  ! outer loop (jc) end

  ! Calculate momentum fluxes for 4 dierctions separately
  IF ( lcalc_flux_4dir_bg ) THEN
    DO jk = 2,nlev
      DO jc = i_startidx, i_endidx
        pmflux_w_bg(jc,jk) = uwflux_bg(jc,jk) - pmflux_e_bg(jc,jk)
        pmflux_s_bg(jc,jk) = vwflux_bg(jc,jk) - pmflux_n_bg(jc,jk)
      ENDDO
    ENDDO
  END IF

  ! half levels --> full levels for the variables defined there
  IF (nrays_add_bg(jg) /= 0) THEN
    DO jc = i_startidx, i_endidx
      ! energy :  half levels --> full levels
      energy_bg  (jc,1) = energy_bg  (jc,2)
      energy_p_bg(jc,1) = energy_p_bg(jc,2)
      DO jk = 2, nlev-1
        energy_bg  (jc,jk) = 0.5_wp*(energy_bg  (jc,jk) + energy_bg  (jc,jk+1))
        energy_p_bg(jc,jk) = 0.5_wp*(energy_p_bg(jc,jk) + energy_p_bg(jc,jk+1))
      ENDDO
    ENDDO
  END IF
  IF ( lcalc_flux_4dir_bg ) THEN
    DO jc = i_startidx, i_endidx
      pmflux_e_bg(jc,1) = 0.5_wp*pmflux_e_bg(jc,2)
      pmflux_w_bg(jc,1) = 0.5_wp*pmflux_w_bg(jc,2)
      pmflux_n_bg(jc,1) = 0.5_wp*pmflux_n_bg(jc,2)
      pmflux_s_bg(jc,1) = 0.5_wp*pmflux_s_bg(jc,2)
    ENDDO
    DO jk = 2, nlev-1
      DO jc = i_startidx, i_endidx
        pmflux_e_bg(jc,jk) = 0.5_wp*(pmflux_e_bg(jc,jk) + pmflux_e_bg(jc,jk+1))
        pmflux_w_bg(jc,jk) = 0.5_wp*(pmflux_w_bg(jc,jk) + pmflux_w_bg(jc,jk+1))
        pmflux_n_bg(jc,jk) = 0.5_wp*(pmflux_n_bg(jc,jk) + pmflux_n_bg(jc,jk+1))
        pmflux_s_bg(jc,jk) = 0.5_wp*(pmflux_s_bg(jc,jk) + pmflux_s_bg(jc,jk+1))
      ENDDO
    ENDDO
    DO jc = i_startidx, i_endidx
      pmflux_e_bg(jc,nlev) = 0.5_wp*pmflux_e_bg(jc,nlev)
      pmflux_w_bg(jc,nlev) = 0.5_wp*pmflux_w_bg(jc,nlev)
      pmflux_n_bg(jc,nlev) = 0.5_wp*pmflux_n_bg(jc,nlev)
      pmflux_s_bg(jc,nlev) = 0.5_wp*pmflux_s_bg(jc,nlev)
    ENDDO
  END IF

END SUBROUTINE saturation
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE saturation_mono(nlev,i_startidx,i_endidx,jray_start,jray_end,dz,kray,lray,mray,dens,iexist,&
                      u_half,v_half,bvf2,gammash2,fc,rho,kvisc,lcalc_flux_4dir,&
                      uwflux,vwflux,energy,apmflux,energy_p,pmflux_e,pmflux_w,pmflux_n,pmflux_s,gw_flag)
  INTEGER,       INTENT(IN)    :: nlev
  INTEGER,       INTENT(IN)    :: i_startidx ! first index of the block
  INTEGER,       INTENT(IN)    :: i_endidx   ! last index of the block
  INTEGER,       INTENT(IN)    :: jray_start
  INTEGER,       INTENT(IN)    :: jray_end
  REAL(wp),      INTENT(IN)    :: dz(:,:)
  INTEGER,       INTENT(IN)    :: iexist(:,:)
  REAL(wp),      INTENT(IN)    :: mray(:,:)
  REAL(wp),      INTENT(IN)    :: kray(:,:)
  REAL(wp),      INTENT(IN)    :: lray(:,:)
  REAL(wp),      INTENT(IN)    :: dens(:,:)
  REAL(wp),      INTENT(IN)    :: u_half(:,:), v_half(:,:)
  REAL(wp),      INTENT(IN)    :: bvf2(:,:), rho(:,:), gammash2(:,:)
  REAL(wp),      INTENT(IN)    :: kvisc(:,:)      ! Kinematic viscosity
  REAL(wp),      INTENT(IN)    :: fc(:)           ! Coriolis parameter
  LOGICAL,       INTENT(IN)    :: lcalc_flux_4dir
  REAL(wp),      INTENT(OUT)   :: uwflux(:,:)
  REAL(wp),      INTENT(OUT)   :: vwflux(:,:)
  REAL(wp),      INTENT(OUT)   :: energy(:,:)
  REAL(wp),      INTENT(OUT)   :: apmflux(:,:)
  REAL(wp),      INTENT(OUT)   :: energy_p(:,:)
  REAL(wp),      INTENT(OUT)   :: pmflux_e(:,:), pmflux_w(:,:)
  REAL(wp),      INTENT(OUT)   :: pmflux_n(:,:), pmflux_s(:,:)
  INTEGER, OPTIONAL, INTENT(IN) :: gw_flag(:)
  ! Counters, indices
  INTEGER                    :: nlevp1
  INTEGER                    :: jk, jc
  INTEGER                    :: jray
  ! Other subquantities
  LOGICAL                    :: gw_flag_not(nproma)
  REAL(wp)                   :: f2_o_bvf2(nlev+1)
  REAL(wp)                   :: al2rho05(nlev+1)
  REAL(wp)                   :: al2_05, denssat
  REAL(wp)                   :: fc2, abs_fc, m2_max
  REAL(wp)                   :: omega2
  REAL(wp)                   :: minus_two_kvisc_dz
  REAL(wp)                   :: Kh2, Kh, inv_Kh2
  REAL(wp)                   :: omega_gr
  REAL(wp)                   :: omegaint
  REAL(wp)                   :: m2, K2
  REAL(wp)                   :: inv_cgz
  REAL(wp)                   :: K2_p_gam2
  REAL(wp)                   :: wadflxprof (jray_start:jray_end,nlev+1)
  REAL(wp)                   :: wadensprof (jray_start:jray_end,nlev+1)
  REAL(wp)                   :: intgcgz_k  (jray_start:jray_end,nlev+1)
  REAL(wp)                   :: intgcgz_l  (jray_start:jray_end,nlev+1)
  REAL(wp)                   :: intgcgz_h  (jray_start:jray_end,nlev+1)
  REAL(wp)                   :: intgomega  (jray_start:jray_end,nlev+1)
  REAL(wp)                   :: intgomega_p(jray_start:jray_end,nlev+1)

  IF (msg_level >= 12) CALL message('saturation_mono', 'MS-GWaM-ST: wave breaking parametrization')

  !----------------------------------------------------------------------
  ! Purpose: 
  !         Parameterize GW breaking (saturation), i.e. the corresponding 
  !         loss in GW energy and momentum fluxes.         
  !
  ! Method:
  !         The method is based on Lindzen's idea that GWs break when 
  !         they turn the potential temperature gradient to statically 
  !         unstable (d\theta/dz < 0).
  !         -- calculate equilibrium profile
  !         -- calculate effects of critical layer filtering on momflux
  !         -- calculate effects of reflection layers on momflux
  !         -- diagnose static instbility by \int dm m^2B^2 >= N^4
  !         -- calculate turbulent viscosity (in case of instability)
  !         -- calculate effects of wave breaking on momflux
  !         -- ...
  !         -- see Bölöni et al. (2020)
  !----------------------------------------------------------------------

  nlevp1 = nlev+1

  m2_max = (2.0_wp*pi/lambda_z_min)**2

  al2_05 = 0.5_wp*alpha_sat**2

  ! Initialization (nproma,nlevp1)

  uwflux(:,:) = 0._wp   ;  vwflux(:,:) = 0._wp  ;  energy(:,:) = 0._wp
  apmflux(:,:) = 0._wp  ;  energy_p(:,:) = 0._wp
  IF ( lcalc_flux_4dir ) THEN
    pmflux_e(:,:) = 0._wp  ;  pmflux_w(:,:) = 0._wp
    pmflux_n(:,:) = 0._wp  ;  pmflux_s(:,:) = 0._wp
  END IF

  IF (PRESENT(gw_flag)) THEN
    gw_flag_not(:) = (gw_flag(:) == 0)
  ELSE
    gw_flag_not(:) = .FALSE.
  END IF

  ! Outer loop (jc) start

  DO jc = i_startidx, i_endidx

    IF ( gw_flag_not(jc) )  CYCLE

    ! Initialization (nrays(jg),nlevp1)
    wadflxprof (:,:) = 0._wp
    wadensprof (:,:) = 0._wp
    intgcgz_h  (:,:) = 0._wp
    intgomega  (:,:) = 0._wp
    intgomega_p(:,:) = 0._wp

    fc2 = fc(jc)**2
    abs_fc = ABS(fc(jc))

    f2_o_bvf2(:) = fc2/bvf2(jc,:)
    al2rho05 (:) = al2_05*rho(jc,:)


    SPECLOOP:  DO jray = jray_start, jray_end

      IF (dens(jc,jray) == 0._wp) CYCLE

      !K_h^2=k^2+l^2
      Kh2 = kray(jc,jray)**2 + lray(jc,jray)**2
      Kh = SQRT(Kh2)
      inv_Kh2 = 1._wp/Kh2

      !=====  SOURCE LEVEL  ====================================================

      ! Calculate the ground-based omega and vertical flux of wave action density
      !   : omega_gr and wadflxprof(jray,jk_source)

      ! at the source level

      jk = iexist(jc,jray)   !  iexist /= 0  IF dens(jc,jray) /= 0._wp

      !K^2=k^2+l^2+m^2
      m2 = mray(jc,jray)**2
      K2 = Kh2 + m2

      K2_p_gam2 = K2 + gammash2(jc,jk)

      ! Calculate intrinsic frequency at launch level
      omega2 = (bvf2(jc,jk)*Kh2+fc2*(m2+gammash2(jc,jk)))/K2_p_gam2
      omegaint = branch*SQRT(omega2)

      ! Calculate the Doppler shift term at launch level
      omega_gr = omegaint  &
        &  + ( kray(jc,jray)*u_half(jc,jk) + lray(jc,jray)*v_half(jc,jk) )

      ! A*cgz
      wadflxprof(jray,jk) = ABS( dens(jc,jray)*( mray(jc,jray)*(omega2 - fc2)  &
        &                        /(omegaint*K2_p_gam2) ) )


      !=====  ABOVE THE SOURCE  ================================================

      ! Calculate the vertical profile of the flux of wave action density and
      !   various variables for diagnostics, given the conserved omega_gr
      !   : wadflxprof(jray,jk), wadensprof(jray,jk), ...

      VERTLOOP:  DO jk = iexist(jc,jray), 2, -1

        ! Calculate the Doppler shift term profile 
        omegaint = omega_gr  &
          &  - ( kray(jc,jray)*u_half(jc,jk) + lray(jc,jray)*v_half(jc,jk) )

        omega2 = omegaint**2
        m2 = (bvf2(jc,jk)-omega2)*Kh2/(omega2-fc2) - gammash2(jc,jk)

        !=====  FILTERING  =====================================================

        ! 1) critical level if:   omegaint^2-f^2 <= 0 or sign change of omegaint
        !                         or too small Lz
        !    --> remove the spectral element above the critical level
        IF ( omegaint*branch <= abs_fc .OR. m2 > m2_max ) THEN
          wadensprof(jray,jk) = 0._wp
          wadflxprof(jray,jk) = 0._wp
          EXIT VERTLOOP
        ENDIF

        ! 2) reflection level if:   m^2 <= 0
        !    --> remove the complete spectral element as in a steady state
        !        picture a reflected (downward) wave carries exactly the opposite mom
        !        flux as the wave before reflection (upward) so they cancel each other
        IF (m2 <= 0._wp) THEN
          wadflxprof(jray,jk+1:) = wadflxprof(jray,jk+1:) - wadflxprof(jray,jk)
          wadensprof(jray,jk) = 0._wp
          wadflxprof(jray,jk) = 0._wp
          EXIT VERTLOOP
        ENDIF

        !=====  SATURATION  ====================================================

        !K^2_low=k^2+l^2+m^2
        K2 = Kh2 + m2

        inv_cgz = omegaint*(K2+gammash2(jc,jk))/(SQRT(m2)*(omega2 - fc2))

        wadensprof(jray,jk) = wadflxprof(jray,jk)*inv_cgz

        denssat = al2rho05(jk)*omegaint*(1._wp/m2 + inv_Kh2)

        IF (wadensprof(jray,jk) > denssat) THEN
          wadensprof(jray,jk) = denssat
          wadflxprof(jray,jk) = denssat/inv_cgz
        END IF

        !=====  CALCFLUX 1 : absolute variables  ===============================

        ! The absolute variables are calculated here before 'wadflx' is changed
        !   when the wave reflection occurs somewhere above this level.
        !   : absolute momentum flux and energy

        intgcgz_h(jray,jk) = Kh*wadflxprof(jray,jk)
        intgomega(jray,jk) = omegaint*wadensprof(jray,jk)
        intgomega_p(jray,jk) = intgomega(jray,jk)*0.5_wp  &
          &  /( 1._wp + f2_o_bvf2(jk)*(m2+gammash2(jc,jk))/Kh2 )

        !-----------------------------------------------------------------------

        ! Pass the flux value to the above level.
        IF ( lmvisc ) THEN
          minus_two_kvisc_dz = -(kvisc(jc,jk-1) + kvisc(jc,jk))*dz(jc,jk-1)
          wadflxprof(jray,jk-1) = wadflxprof(jray,jk)  &
            &                    *EXP(minus_two_kvisc_dz*inv_cgz*K2)
        ELSE
          wadflxprof(jray,jk-1) = wadflxprof(jray,jk)
        END IF

      ENDDO VERTLOOP

    ENDDO SPECLOOP


    !=====  CALCFLUX 2 : others  ===============================================

    ! The fluxes are calculated after vanishing 'wadflx' for reflected waves.

    DO jk = 2, nlev
      intgcgz_k(:,jk) = kray(jc,jray_start:jray_end)*wadflxprof(:,jk)
      intgcgz_l(:,jk) = lray(jc,jray_start:jray_end)*wadflxprof(:,jk)
    ENDDO 

    DO jk = 2, nlev
      uwflux(jc,jk) = SUM(intgcgz_k(:,jk))
      vwflux(jc,jk) = SUM(intgcgz_l(:,jk))
      energy(jc,jk) = SUM(intgomega(:,jk))
    ENDDO 
    DO jk = 2, nlev
      apmflux (jc,jk) = SUM(intgcgz_h  (:,jk))
      energy_p(jc,jk) = SUM(intgomega_p(:,jk))
    ENDDO 

    ! Upper boundary : no transmission to conserve momentum
    uwflux(jc,1) = 0._wp  ;  vwflux(jc,1) = 0._wp
    ! Lower boundary : vertical flux = 0
    uwflux(jc,nlevp1) = 0._wp  ;  vwflux(jc,nlevp1) = 0._wp

    ! The other diagnostic outputs which are not directly involved in the tendency
    !   calculation will have zero values at jk = 1 and nlevp1, as initialized.
    !   : absolute flux and fluxes at each direction

    ! Calculate momentum fluxes for 4 dierctions separately
    IF ( lcalc_flux_4dir ) THEN
      DO jray = jray_start, jray_end
        IF (kray(jc,jray) > 0._wp) THEN
          DO jk = 2, nlev
            pmflux_e(jc,jk) = pmflux_e(jc,jk) + intgcgz_k(jray,jk)
          ENDDO
        END IF
        IF (lray(jc,jray) > 0._wp) THEN
          DO jk = 2, nlev
            pmflux_n(jc,jk) = pmflux_n(jc,jk) + intgcgz_l(jray,jk)
          ENDDO
        END IF
      ENDDO
    END IF

  ENDDO ! jc

  ! Outer loop (jc) end

  ! Calculate momentum fluxes for 4 dierctions separately
  IF ( lcalc_flux_4dir ) THEN
    DO jk = 2,nlev
      DO jc = i_startidx, i_endidx
        pmflux_w(jc,jk) = uwflux(jc,jk) - pmflux_e(jc,jk)
        pmflux_s(jc,jk) = vwflux(jc,jk) - pmflux_n(jc,jk)
      ENDDO
    ENDDO
  END IF

  ! half levels --> full levels for the variables defined there
  DO jc = i_startidx, i_endidx
    energy  (jc,1) = energy  (jc,2)
    energy_p(jc,1) = energy_p(jc,2)
    DO jk = 2, nlev-1
      energy  (jc,jk) = 0.5_wp*(energy  (jc,jk) + energy  (jc,jk+1))
      energy_p(jc,jk) = 0.5_wp*(energy_p(jc,jk) + energy_p(jc,jk+1))
    ENDDO
  ENDDO ! jc
  IF ( lcalc_flux_4dir ) THEN
    DO jc = i_startidx, i_endidx
      pmflux_e(jc,1) = 0.5_wp*pmflux_e(jc,2)
      pmflux_w(jc,1) = 0.5_wp*pmflux_w(jc,2)
      pmflux_n(jc,1) = 0.5_wp*pmflux_n(jc,2)
      pmflux_s(jc,1) = 0.5_wp*pmflux_s(jc,2)
    ENDDO
    DO jk = 2, nlev-1
      DO jc = i_startidx, i_endidx
        pmflux_e(jc,jk) = 0.5_wp*(pmflux_e(jc,jk) + pmflux_e(jc,jk+1))
        pmflux_w(jc,jk) = 0.5_wp*(pmflux_w(jc,jk) + pmflux_w(jc,jk+1))
        pmflux_n(jc,jk) = 0.5_wp*(pmflux_n(jc,jk) + pmflux_n(jc,jk+1))
        pmflux_s(jc,jk) = 0.5_wp*(pmflux_s(jc,jk) + pmflux_s(jc,jk+1))
      ENDDO
    ENDDO
    DO jc = i_startidx, i_endidx
      pmflux_e(jc,nlev) = 0.5_wp*pmflux_e(jc,nlev)
      pmflux_w(jc,nlev) = 0.5_wp*pmflux_w(jc,nlev)
      pmflux_n(jc,nlev) = 0.5_wp*pmflux_n(jc,nlev)
      pmflux_s(jc,nlev) = 0.5_wp*pmflux_s(jc,nlev)
    ENDDO
  END IF

END SUBROUTINE saturation_mono
!!!
!!!-------------------------------------------------------------------------
!!!
SUBROUTINE tendency(nlev,i_startidx,i_endidx,dz,rho,uwflux,vwflux,&
                     utend,vtend)
  INTEGER,        INTENT(IN)  :: nlev       ! no. of full levels
  INTEGER,        INTENT(IN)  :: i_startidx ! first index of the block
  INTEGER,        INTENT(IN)  :: i_endidx   ! last index of the block
  REAL(wp),       INTENT(IN)  :: dz(:,:)      ! full layer thickness
  REAL(wp),       INTENT(IN)  :: rho(:,:)     ! density at full levels
  REAL(wp),       INTENT(IN)  :: uwflux(:,:)  ! uw flux
  REAL(wp),       INTENT(IN)  :: vwflux(:,:)  ! vw flux
  REAL(wp),       INTENT(OUT) :: utend(:,:)   ! u tendency
  REAL(wp),       INTENT(OUT) :: vtend(:,:)   ! v tendency

  INTEGER                     :: jc, jk
  REAL(wp)                    :: inv_dzrho

  IF (msg_level >= 12) CALL message('tendency', 'MS-GWaM-ST: wind tendency computation')

  !----------------------------------------------------------------------
  ! Purpose:
  !         Calculate wind tendencies based on pseudo-momentum flux 
  !         convergence (vertical).
  !
  ! Method:
  !         tendencies = 1/rho*(vertical centered differences of pseudo-
  !                      momentum fluxes)
  !         staggering:
  !         -- pseudo-momentum fluxes defined on half levels
  !         -- tendencies defined on full levels
  !----------------------------------------------------------------------

  DO jk = 1,nlev

    DO jc = i_startidx,i_endidx

      ! 1/(rho*dz)
      inv_dzrho = 1._wp/(dz(jc,jk)*rho(jc,jk))

      ! Calculate tendency on full levels
      utend(jc,jk) = - (uwflux(jc,jk) - uwflux(jc,jk+1))*inv_dzrho
      vtend(jc,jk) = - (vwflux(jc,jk) - vwflux(jc,jk+1))*inv_dzrho

    ENDDO ! jc

  ENDDO ! jk

END SUBROUTINE tendency
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE smooth_vert(nlev,i_startidx,i_endidx,npts, var)
  INTEGER,  INTENT(IN)    :: nlev
  INTEGER,  INTENT(IN)    :: i_startidx ! first index of the block
  INTEGER,  INTENT(IN)    :: i_endidx   ! last index of the block
  INTEGER,  INTENT(IN)    :: npts       ! half number of points to smooth
  REAL(wp), INTENT(INOUT) :: var(:,:)

  INTEGER                 :: nstenc
  INTEGER                 :: jc, jk, jks, jk2
  REAL(wp)                :: var_tmp(nproma,nlev)
  REAL(wp), PARAMETER     :: stencil(-3:3,0:3) = RESHAPE( (/                              &
    &  0._wp      ,0._wp      ,0._wp      ,1._wp    ,0._wp      ,0._wp      ,0._wp      , &
    &  0._wp      ,0._wp      ,0.25_wp    ,0.5_wp   ,0.25_wp    ,0._wp      ,0._wp      , &
    &  0._wp      ,-0.0625_wp ,0.25_wp    ,0.625_wp ,0.25_wp    ,-0.0625_wp ,0._wp      , &
    &  0.015625_wp,-0.09375_wp,0.234375_wp,0.6875_wp,0.234375_wp,-0.09375_wp,0.015625_wp  &
    &  /), (/7,4/) )

  !----------------------------------------------------------------------
  ! Purpose: 
  !         Smooth pseudo-momentum fluxes and energy in order to get rid 
  !         of waves with length scales <= 2dz. Such waves in the 
  !         projected fluxes / energy might be present e.g. due to sparse 
  !         ray volume distributions. The preferred smoothing method is a 
  !         0th order Shapiro filter.
  !
  ! Method:
  !         -- see Shapiro (1975)
  !----------------------------------------------------------------------

  ! No. of half points for the smoothing stencil
  nstenc = ABS(npts)
  ! ABS() is kept to avoid problems with older scripts where npts>0 
  ! meant a simpler moving average smooting (removed now).

  ! Keep the boundary conditions
  DO jc = i_startidx, i_endidx
    var_tmp(jc,1) = var(jc,1)  ;  var_tmp(jc,nlev) = var(jc,nlev)
  ENDDO
  ! Near-boundary layers: 3-pt, 5-pt, ...
  DO jk = 2,nstenc
    var_tmp(:,jk) = 0._wp  ;  var_tmp(:,nlev+1-jk) = 0._wp
    DO jks = 1-jk,jk-1
      jk2 = jk + jks
      DO jc = i_startidx, i_endidx
        var_tmp(jc,jk) = var_tmp(jc,jk) + var(jc,jk2)*stencil(jks,jk-1)
        var_tmp(jc,nlev+1-jk) = var_tmp(jc,nlev+1-jk)  &
          &                         + var(jc,nlev+1-jk2)*stencil(jks,jk-1)
      ENDDO ! jc
    ENDDO ! jks
    ! Update before calculating the next inner layer
    DO jc = i_startidx, i_endidx
      var(jc,jk) = var_tmp(jc,jk)
      var(jc,nlev+1-jk) = var_tmp(jc,nlev+1-jk)
    ENDDO ! jc
  ENDDO ! jk
  ! Inside: (nstenc)-point smoothing
  DO jk = 1+nstenc,nlev-nstenc
    DO jc = i_startidx, i_endidx
      var_tmp(jc,jk) = var(jc,jk)
      var(jc,jk) = 0._wp
    ENDDO
  ENDDO
  DO jks = -nstenc,nstenc
    DO jk = 1+nstenc,nlev-nstenc
      jk2 = jk + jks
      DO jc = i_startidx, i_endidx
        var(jc,jk) = var(jc,jk) + var_tmp(jc,jk2)*stencil(jks,nstenc)
      ENDDO ! jc
    ENDDO ! jk
  ENDDO ! jks

END SUBROUTINE smooth_vert
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
  REAL(wp), INTENT(in) :: Ta   ! ambient temp. [K]
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
END MODULE mo_msgwam_stst
