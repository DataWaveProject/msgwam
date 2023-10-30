!>
!! This module is to setup an isothermal atmosphere at rest for 
!! an idealized case to test gravity wave drag parametrization 
!! developed at the Goethe Uni: Lagrangian WKB raytracer in phase space.
!!
!! @author Gergely Bölöni, Goethe Uni Frankfurt 
!!         (based on mo_nh_dcmip_rest_atm.f90)
!!
!! @par Revision History
!! Initial version by Gergely Bölöni (2016-02-25)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_nh_isotherm_rest_atm

   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd, cpd, cvd, grav, rd_o_cpd, p0ref
   USE mo_math_constants,       ONLY: pi
   USE mo_impl_constants,       ONLY: min_rlcell, min_rledge, min_rlvert
   USE mo_parallel_config,      ONLY: nproma, p_test_run
   USE mo_loopindices,          ONLY: get_indices_c, get_indices_e, get_indices_v
   USE mo_model_domain,         ONLY: t_patch
   USE mo_grid_config,          ONLY: grid_sphere_radius
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_ref, t_nh_metrics
   USE mo_intp_data_strc,       ONLY: t_int_state
   USE mo_intp,                 ONLY: cells2edges_scalar
   USE mo_ext_data_types,       ONLY: t_external_data
   USE mo_exception,            ONLY: message, finish, message_text
   USE mo_sync,                 ONLY: sync_patch_array, sync_patch_array_mult, &
     &                                SYNC_C, SYNC_E
   USE mo_nh_init_utils,        ONLY: convert_thdvars  !, init_w 
   USE mo_hydro_adjust,         ONLY: hydro_adjust
   USE mo_mpi,                  ONLY: my_process_is_stdio
   USE mo_run_config,           ONLY: msg_level

   IMPLICIT NONE


   PRIVATE

   PUBLIC :: init_nh_isotherm_rest_atm

!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
  !! isothermal steady state atmosphere at rest
  !!
  !!

  SUBROUTINE init_nh_isotherm_rest_atm(p_patch, p_nh_prog, p_nh_diag, &
    &                                  p_metrics, p_int, l_hydro_adjust )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  p_patch

    TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
      &  p_nh_prog

    TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
      &  p_nh_diag

    TYPE(t_nh_metrics),   INTENT(IN)    :: &  !< NH metrics state
      &  p_metrics

    TYPE(t_int_state),    INTENT(IN)    :: p_int

    LOGICAL,              INTENT(IN)    :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                          ! initial condition

    !local variables
    INTEGER  :: je, jc, jb, jk,   jn, &
                nlen, nblks_e, nblks_c, npromz_c
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend, i_nchdom  
    INTEGER  :: nlev, nlevp1          !< number of full/half levels

    REAL(wp) :: zu(p_patch%nlev), zv(p_patch%nlev)
    REAL(wp) :: hp, htheta, kappa, rho0

!   !DEFINED PARAMETERS for the Schaer-type testcase (DCMIP):
    REAL(wp), PARAMETER :: t0 = 300._wp        ! Temperature (K)

!--------------------------------------------------------------------

    IF (msg_level >= 12) CALL message('init_nh_isothermal_rest_atm', &
      'Initialize an isothermal atmosphere at rest.')

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! number of child domains
    i_nchdom = MAX(1,p_patch%n_childdom)

    !
    ! Init prognostic variables vn, w to zero 
    ! and set up an isothermal atmosphere
    !
!!$OMP PARALLEL PRIVATE(i_startblk,i_endblk,i_rlstart,i_rlend)

    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)

!!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx, zu, zv)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          p_nh_prog%vn(je,jk,jb) = 0._wp

        ENDDO  ! je
      ENDDO  ! jk
    ENDDO ! jb
!!$OMP ENDDO

   i_rlstart = 1
   i_rlend   = min_rlcell

   i_startblk = p_patch%cells%start_blk(i_rlstart,1)
   i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

   ! define density/pressure/pot.temperature scale heights
   kappa  = (cpd-cvd)/cpd
   hp     = rd*t0/grav
   htheta = hp/kappa
   rho0   = p0ref/rd/t0

!!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
   DO jb = i_startblk, i_endblk

     CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                        i_startidx, i_endidx, i_rlstart, i_rlend)

     DO jk = 1, nlev
       DO jc = i_startidx, i_endidx

         ! init temp, pres, rho
         p_nh_diag%temp(jc,jk,jb)  = t0
         p_nh_diag%pres(jc,jk,jb)  = p0ref   * exp(-p_metrics%z_mc(jc,jk,jb)/hp)
         p_nh_prog%rho(jc,jk,jb)   = p_nh_diag%pres(jc,jk,jb) / rd / t0

         ! init theta_v (note that for dry air theta_v=theta)
         p_nh_prog%theta_v(jc,jk,jb) = t0   * (p0ref/p_nh_diag%pres(jc,jk,jb))**rd_o_cpd

         ! initialize vertical velocity
         p_nh_prog%w(jc,jk,jb) = 0._wp

       ENDDO !jc
     ENDDO !jk

     DO jk = 1, nlevp1
       DO jc = i_startidx, i_endidx

         ! init temp, pres, rho
         p_nh_diag%temp_ifc(jc,jk,jb)  = t0
         p_nh_diag%pres_ifc(jc,jk,jb)  = p0ref   * exp(-p_metrics%z_ifc(jc,jk,jb)/hp)
         p_nh_diag%rho_ic(jc,jk,jb)    = p_nh_diag%pres_ifc(jc,jk,jb) / rd / t0

         ! init theta_v (note that for dry air theta_v=theta)
         p_nh_diag%theta_v_ic(jc,jk,jb) = t0   * (p0ref/p_nh_diag%pres_ifc(jc,jk,jb))**rd_o_cpd

       ENDDO !jc
     ENDDO !jk

     ! values at lower boundary
     p_nh_prog%w(i_startidx:i_endidx,nlevp1,jb)    = 0._wp

     ! initialize surface pressure
     p_nh_diag%pres_sfc(i_startidx:i_endidx,jb)    = p0ref

   ENDDO !jb
!!$OMP END DO
!!$OMP END PARALLEL

   ! Initialize Exner pressure
   CALL sync_patch_array_mult(SYNC_C, p_patch, 2, p_nh_diag%temp, p_nh_diag%pres)

   CALL convert_thdvars(p_patch, p_nh_diag%pres, p_nh_diag%temp, &
                      & p_nh_prog%rho, p_nh_prog%exner, p_nh_prog%theta_v  )

   CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_nh_prog%rho, p_nh_diag%rho_ic, &
                                                  p_nh_prog%theta_v, p_nh_diag%theta_v_ic, &
                                                  p_nh_prog%exner)
   CALL sync_patch_array_mult(SYNC_C, p_patch, 3, p_nh_diag%u, p_nh_diag%v, p_nh_prog%w)


   IF (l_hydro_adjust) THEN

     CALL hydro_adjust ( p_patch, p_metrics, p_nh_prog%rho, &
                     & p_nh_prog%exner, p_nh_prog%theta_v   )

     CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_nh_prog%rho, p_nh_diag%rho_ic, &
                                                    p_nh_prog%theta_v, p_nh_diag%theta_v_ic, &
                                                    p_nh_prog%exner)
     CALL sync_patch_array_mult(SYNC_C, p_patch, 3, p_nh_diag%u, p_nh_diag%v, p_nh_prog%w)
   END IF

  END SUBROUTINE init_nh_isotherm_rest_atm
!--------------------------------------------------------------------
  END MODULE mo_nh_isotherm_rest_atm
