!>
!! This module is to setup an isothermal atmosphere at rest for 
!! an idealized case to test gravity wave drag parametrization 
!! developed at the Goethe Uni: Lagrangian WKB raytracer in phase space.
!!
!! @author Gergely Bölöni, Goethe Uni Frankfurt 
!!         (based on mo_nh_isotherm_rest_atm.f90)
!!
!! @par Revision History
!! Initial version by Gergely Bölöni (2021-03-01)
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
MODULE mo_nh_isotherm_rest_atm_gwp

   USE mo_kind,                 ONLY: wp
   USE mo_physical_constants,   ONLY: rd, cpd, cvd, grav, rd_o_cpd
   USE mo_math_constants,       ONLY: pi, rad2deg, deg2rad
   USE mo_impl_constants,       ONLY: min_rlcell, min_rledge, min_rlvert
   USE mo_parallel_config,      ONLY: nproma, p_test_run
   USE mo_loopindices,          ONLY: get_indices_c, get_indices_e, get_indices_v
   USE mo_model_domain,         ONLY: t_patch
   USE mo_grid_config,          ONLY: grid_sphere_radius
   USE mo_grid_geometry_info,   ONLY: planar_torus_geometry, sphere_geometry
   USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_ref, t_nh_metrics
   USE mo_intp_data_strc,       ONLY: t_int_state
   USE mo_intp,                 ONLY: cells2edges_scalar
   USE mo_ext_data_types,       ONLY: t_external_data
   USE mo_exception,            ONLY: message, finish, message_text
   USE mo_run_config,           ONLY: msg_level
   USE mo_sync,                 ONLY: sync_patch_array, sync_patch_array_mult, &
     &                                SYNC_C, SYNC_E
   USE mo_nh_init_utils,        ONLY: convert_thdvars  !, init_w 
   USE mo_hydro_adjust,         ONLY: hydro_adjust
   USE mo_mpi,                  ONLY: my_process_is_stdio, p_barrier, p_comm_work
   USE mo_msgwam,               ONLY: sync_wave
   USE mo_msgwam_config
   USE mo_setup_msgwam_interface

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: init_nh_isotherm_rest_atm_gwp

!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
  !! isothermal steady state atmosphere at rest
  !!
  !!

  SUBROUTINE init_nh_isotherm_rest_atm_gwp(p_patch, p_nh_prog, p_nh_diag, &
    &                                  p_gridinfo4ray, p_ray, p_metrics, &
                                       p_int, l_hydro_adjust )

  TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
    &  p_patch

  TYPE(t_nh_prog),      INTENT(INOUT) :: &  !< prognostic state vector
    &  p_nh_prog

  TYPE(t_nh_diag),      INTENT(INOUT) :: &  !< diagnostic state vector
    &  p_nh_diag

  TYPE(t_gridinfo4ray), INTENT(INOUT) :: &  !< MS-GWaM stencil for neighbor search
    &  p_gridinfo4ray

  TYPE(t_ray),          INTENT(INOUT) :: &  !< MS-GWaM ray volumes
    &  p_ray

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
  INTEGER  :: jkmid, jz, jlon, jlat, jmm, jkk, jll, jray
  INTEGER  :: jvec, jmin, jmax, jthird
  INTEGER  :: jkmin_ini_eff, jkmax_ini_eff

  REAL(wp) :: zu(p_patch%nlev), zv(p_patch%nlev)
  REAL(wp) :: hp, bvf2
  REAL(wp) :: zray, dzray, mray, dmray, lray, dlray, kray, dkray
  REAL(wp) :: lonray, lonray_ref, latray
  REAL(wp) :: dlonray, dlatray, delta_lon, distance, base_small
  REAL(wp) :: omega, Bw, env             ! intrinsic freq, buoy amplitude
  REAL(wp) :: vec_ref(3), vec_thispos(3), dx, dy, dz
  REAL(wp) :: wa(nproma,p_patch%nlev,p_patch%nblks_c) ! wave action density on model levels
  REAL(wp) :: minlat, maxlat
  REAL(wp) :: miny, maxy
  REAL(wp) :: xs, ys
  REAL(wp) :: vecx(3), vecy(3)

  !----------------------------------------------------------------------
  ! Purpose:
  !         Enable idealized test cases of MS-GWaM for i
  !         validation purposes:
  !         1) Set up isothermal atmosphere at rest
  !         2) Initialize idealized gravity wave packet
  ! Method:
  !         -- see Muraschko et al. (2015), Bölöni et al. (2016), 
  !                Wei et al. (2019)
  !----------------------------------------------------------------------

  IF (msg_level >= 12) CALL message('init_nh_isothermal_rest_atm_gwp', &
    'Initialize an isothermal atmosphere at rest and an idealized GW packet.')

  !
  ! Part 1: set up isothermal atmosphere at rest
  !

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

  ! Pressure scale height
  hp = rd*t0/grav

  ! Brunt-Väisälä frequency
  bvf2 = grav**2/cpd/t0

  i_rlstart = 1
  i_rlend   = min_rlcell

  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, i_rlstart, i_rlend)

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        ! init temp, pres, rho
        p_nh_diag%temp(jc,jk,jb) = t0
        p_nh_diag%pres(jc,jk,jb) = p0 * EXP(-p_metrics%z_mc(jc,jk,jb)/hp)
        p_nh_prog%rho(jc,jk,jb)  = p_nh_diag%pres(jc,jk,jb) / rd / t0

        ! init theta_v (note that for dry air theta_v=theta)
        p_nh_prog%theta_v(jc,jk,jb) = t0 * (p0/p_nh_diag%pres(jc,jk,jb))**rd_o_cpd

        ! initialize vertical velocity
        p_nh_prog%w(jc,jk,jb) = 0._wp

      ENDDO !jc
    ENDDO !jk

    DO jk = 1, nlevp1
      DO jc = i_startidx, i_endidx

        ! init temp, pres, rho
        p_nh_diag%temp_ifc(jc,jk,jb) = t0
        p_nh_diag%pres_ifc(jc,jk,jb) = p0 * EXP(-p_metrics%z_ifc(jc,jk,jb)/hp)
        p_nh_diag%rho_ic(jc,jk,jb)   = p_nh_diag%pres_ifc(jc,jk,jb) / rd / t0

        ! init theta_v (note that for dry air theta_v=theta)
        p_nh_diag%theta_v_ic(jc,jk,jb) = t0 * (p0/p_nh_diag%pres_ifc(jc,jk,jb))**rd_o_cpd

      ENDDO !jc
    ENDDO !jk

    ! values at lower boundary
    p_nh_prog%w(i_startidx:i_endidx,nlevp1,jb) = 0._wp

    ! initialize surface pressure
    p_nh_diag%pres_sfc(i_startidx:i_endidx,jb) = p0

  ENDDO !jb
!!$OMP END DO
!!$OMP END PARALLEL

  ! Initialize Exner pressure
  CALL sync_patch_array_mult(SYNC_C, p_patch, 2, p_nh_diag%temp, p_nh_diag%pres)

  CALL convert_thdvars(p_patch, p_nh_diag%pres, p_nh_diag%temp, &
                     & p_nh_prog%rho, p_nh_prog%exner, p_nh_prog%theta_v  )

  ! Synchronize
  CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_nh_prog%rho, p_nh_diag%rho_ic, &
                                                 p_nh_prog%theta_v, p_nh_diag%theta_v_ic, &
                                                 p_nh_prog%exner)
  CALL sync_patch_array_mult(SYNC_C, p_patch, 3, p_nh_diag%u, p_nh_diag%v, p_nh_prog%w)


  ! Hydrostatic adjustment
  IF (l_hydro_adjust) THEN

    CALL hydro_adjust ( p_patch, p_metrics, p_nh_prog%rho, &
                    & p_nh_prog%exner, p_nh_prog%theta_v   )

    CALL sync_patch_array_mult(SYNC_C, p_patch, 5, p_nh_prog%rho, p_nh_diag%rho_ic, &
                                                   p_nh_prog%theta_v, p_nh_diag%theta_v_ic, &
                                                   p_nh_prog%exner)
    CALL sync_patch_array_mult(SYNC_C, p_patch, 3, p_nh_diag%u, p_nh_diag%v, p_nh_prog%w)
  END IF

  !
  ! Part 2: initialize idealized gravity wave packet
  !

  ! Initialize
  wa(:,:,:) = 0._wp

  ! Layer indices for ray volume initialization
  jkmin_ini_eff = MAX(2,jkmin_ini)
  jkmax_ini_eff = MIN(nlev,jkmax_ini)

  ! Define horizontal and vertical wavenumbers
  IF (lambda_lon_ini >= 1.e9) THEN
    k_ini = 0._wp
  ELSE
    k_ini = 2._wp*pi/lambda_lon_ini
  ENDIF
  IF (lambda_lat_ini >= 1.e9) THEN
    l_ini = 0._wp
  ELSE
    l_ini = 2._wp*pi/lambda_lat_ini
  ENDIF
  m_ini = 2._wp*pi/lambda_z_ini

  ! Check if the vertical grid is equidistant and warn if not
  IF (p_metrics%ddqz_z_full(1,nlev,1) /= p_metrics%ddqz_z_full(1,1,1)) THEN
    WRITE(message_text,'(a)') 'Warning: vertical grid not equidistant!  --> no pb, '// &
         'but initial ray geometry does not exactly match the vertical model grid...'
    CALL message  ('init_nh_isotherm_rest_atm_gwp', TRIM(message_text))
  ENDIF

  i_rlstart = 1
  i_rlend   = min_rlcell

  i_startblk = p_patch%cells%start_blk(i_rlstart,1)
  i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

  ! TODO: try new style of index calculations 
  ! It would be:
  ! Prognostic domain
  !rl_start   = grf_bdywidth_c + 1
  !rl_end     = min_rlcell_int
  !i_startblk = p_patch%cells%start_block(rl_start)
  !i_endblk   = p_patch%cells%end_block(rl_end)

  SELECT CASE(p_patch%geometry_info%geometry_type)

  CASE (sphere_geometry)

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev

        DO jc = i_startidx, i_endidx

          ! Intrinsic frequency
          omega = branch*SQRT((bvf2*(k_ini**2+l_ini**2) + &
                               p_patch%cells%f_c(jc,jb)**2*m_ini**2) / &
                              (k_ini**2 + l_ini**2 + m_ini**2))

          ! Buoyancy amplitude
          Bw = a0*bvf2/m_ini

          ! Buoyancy amplitude and wave action density at full levels
          ! CAUTION: consider that the envelope is a squared Gaussian / raised cosine
          ! as it is the envelope of the wave action rather than the bouyancy.

          ! Calculate vertical distance
          dz = p_metrics%z_mc(jc,jk,jb)-z0

          IF (envtype==0) THEN ! no envelop
            env = 1._wp
          ELSEIF (envtype==1) THEN ! 1D Gaussian
            env = EXP(-dz**2/(sigmaz_wp**2))
          ELSEIF (envtype==2) THEN ! 1D Cosine
             IF (ABS(dz)<sigmaz_wp) THEN
               env = (1._wp+COS(pi*dz/sigmaz_wp))**2
             ELSE
               env = 0._wp
             ENDIF
          ELSEIF (envtype==3) THEN ! 2D Gaussian: envelop in latitude, 
                                   ! zonally symmetric
            ! Cartesian vector of this locataion 
            vec_thispos = (/COS(p_patch%cells%center(jc,jb)%lat)*COS(p_patch%cells%center(jc,jb)%lon),&
                            COS(p_patch%cells%center(jc,jb)%lat)*SIN(p_patch%cells%center(jc,jb)%lon),&
                            SIN(p_patch%cells%center(jc,jb)%lat)/)
            ! Cartesian vector of reference position in lat
            vec_ref = (/COS(lat0*deg2rad)*COS(p_patch%cells%center(jc,jb)%lon),&
                        COS(lat0*deg2rad)*SIN(p_patch%cells%center(jc,jb)%lon),&
                        SIN(lat0*deg2rad)/)
            ! Calculate great-circle distance based on chord length in meridional direction
            dy = 2.0_wp*ASIN(MIN(1._wp, 0.5_wp*NORM2(vec_thispos-vec_ref)))*(grid_sphere_radius+z0)
            env = EXP(-dy**2/(sigmay_wp**2) &
                      -dz**2/(sigmaz_wp**2))
          ELSEIF (envtype==4) THEN ! 3D Gaussian
            ! Cartesian vector of this locataion 
            vec_thispos = (/COS(p_patch%cells%center(jc,jb)%lat)*COS(p_patch%cells%center(jc,jb)%lon),&
                            COS(p_patch%cells%center(jc,jb)%lat)*SIN(p_patch%cells%center(jc,jb)%lon),&
                            SIN(p_patch%cells%center(jc,jb)%lat)/)
            ! Cartesian vector of reference position in lat
            vec_ref = (/COS(lat0*deg2rad)*COS(lon0*deg2rad),&
                        COS(lat0*deg2rad)*SIN(lon0*deg2rad),&
                        SIN(lat0*deg2rad)/)
            ! Calculate great-circle distance based on chord length in meridional direction
            dx = 2.0_wp*ASIN(MIN(1._wp, 0.5_wp*NORM2(vec_thispos-vec_ref)))*(grid_sphere_radius+z0)
            dy = dx
            env = EXP(-dx**2/(sigmax_wp**2) &
                      -dy**2/(sigmay_wp**2) &
                      -dz**2/(sigmaz_wp**2))
          ELSE
            CALL finish('init_nh_isothermal_rest_atm_gwp', &
                        'envtype should be 1, 2, 3 or 4!')
          ENDIF
          wa(jc,jk,jb) = env*p_nh_prog%rho(jc,jk,jb)/2._wp * &
                        (1._wp+m_ini**2/(k_ini**2+l_ini**2)) * &
                         omega*Bw**2/bvf2**2

        ENDDO ! jc

      ENDDO ! jk

      DO jc = i_startidx, i_endidx

        ! Horizontal propagation test
        ! Horizontal localization of the source: skip calculations 
        ! if ltest_hprop=.true. and cell is out of the area defined 
        ! by lonmin_bg, lonrmax, latrmin, latrmax
        IF (ltest_hprop) THEN
          IF (p_patch%cells%center(jc,jb)%lon*rad2deg > lonrmax .OR. &
              p_patch%cells%center(jc,jb)%lon*rad2deg < lonrmin .OR. &
              p_patch%cells%center(jc,jb)%lat*rad2deg > latrmax .OR. &
              p_patch%cells%center(jc,jb)%lat*rad2deg < latrmin) THEN
            CYCLE
          ELSE
            ! Write coordinates to file
            OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
            WRITE(987,'(a,2x,F14.5,2x,F14.5)') &
                  'GW source lon, lat coordinates:', &
                   p_patch%cells%center(jc,jb)%lon*rad2deg, &
                   p_patch%cells%center(jc,jb)%lat*rad2deg
            CLOSE(987)
          ENDIF
        ENDIF

        ! Calculate length of the line parallel to the triangle base 
        ! at the height of the center point (1/3h with h being 
        ! the triangle's height). This is the same as calclulating 
        ! the base of a smaller similar equilateral triangle of 
        ! height h_s = 2/3h. Then the "small base" of the smaller 
        ! triangle can be calculated as a_s = 2/sqrt(3)*hs = 4/(3*sqrt(3))*h
        base_small = 4._wp / 3._wp / SQRT(3._wp) * height_triangle

        ! Calculate min, max latitude of cell
        minlat = MINVAL(p_gridinfo4ray%cellvertices_lat(:,jc,jb))
        maxlat = MAXVAL(p_gridinfo4ray%cellvertices_lat(:,jc,jb))

        ! Save distance by inverting for the Haversine formula.
        distance = SQRT(p_patch%cells%area(jc,jb))

        jray = 0

        DO jk = jkmax_ini_eff, jkmin_ini_eff, -1

           ! Loop over initial ray volumes within a layer
           DO jz = 1, nrz

             ! Loop over initial ray volumes within lon direction
             DO jlon = 1, nrlon

               ! Loop over initial ray volumes within lat direction
               DO jlat = 1, nrlat

                  ! Loop over initial ray volumes within dm_ini
                  DO jmm = 1, nrm

                    ! Loop over initial ray volumes within dl_ini
                    DO jll = 1, nrl

                      ! Loop over initial ray volumes within dk_ini
                      DO jkk = 1, nrk

                        ! Vertical extent of ray volume
                        dzray = (p_metrics%z_ifc(jc,jk-1,jb) - &
                                 p_metrics%z_ifc(jc,jk,  jb))/REAL(nrz,wp)

                        ! Altitude of ray volume center point
                        zray = p_metrics%z_ifc(jc,jk,jb) + (jz-0.5_wp) * dzray

                        ! lat position of ray volume center point:
                        latray  = minlat + (jlat-0.5_wp) * &
                                 (maxlat-minlat) / REAL(nrlat,wp)

                        ! Angle difference in lon corresponding to the distance 
                        ! "base_small" based on inverting the Haversine formula
                        delta_lon = 2._wp * &
                                    ASIN(MIN(1._wp, SIN(base_small/grid_sphere_radius/2.) &
                                         / COS(latray)))

                        ! lon position of ray volume center point:
                        lonray  = p_patch%cells%center(jc,jb)%lon &
                                - 0.5_wp * delta_lon              &
                                + (jlon-0.5_wp) * delta_lon / REAL(nrlon,wp)

                        ! For the horizontal extent of the ray volume so far simply 
                        ! an extent is taken of which the square gives the area of 
                        ! the triangle cell. The extent in lon, lat is calculated 
                        ! from the distance by inverting the Haversine formula.
                        ! lat extent of ray volume 
                        dlatray  = distance / grid_sphere_radius / REAL(nrlat,wp)
                        ! lon extent of ray volume 
!                       dlonray  = 2._wp*ASIN(SIN(distance/grid_sphere_radius/2.) &
!                                               / COS(latray)) / REAL(nrlon,wp)
                        dlonray  = distance / grid_sphere_radius / COS(latray) / REAL(nrlon,wp)
                                 ! such that dV = (r0^2 coslat) dlon dlat dz = dA dz (dA: cell area)

                        ! Vertical wavenumber position of ray volume center point
                        mray = m_ini - dm_ini/2._wp + (REAL(jmm,wp)-0.5_wp)*dm_ini/REAL(nrm,wp)

                        ! Vertical wavenumber extent of ray volume
                        dmray = dm_ini/REAL(nrm,wp)

                        ! Meridional wavenumber position of ray volume center point
                        lray = l_ini - dl_ini/2._wp + (REAL(jll,wp)-0.5_wp)*dl_ini/REAL(nrl,wp)

                        ! Meridional wavenumber extent of ray volume
                        dlray = dl_ini/REAL(nrl,wp)

                        ! Latitudinal wavenumber position of ray volume center point
                        kray = k_ini - dk_ini/2._wp + (REAL(jkk,wp)-0.5_wp)*dk_ini/REAL(nrk,wp)

                        ! Latitudinal wavenumber extent of ray volume
                        dkray = dk_ini/REAL(nrk,wp)

                        ! Vertical index of ray volume center point
                        jkmid = jk   ! closest half level
!                        IF (jz > CEILING(nrz*0.5_wp))  jkmid = jk-1
                        IF (jz > nrz/2)  jkmid = jk-1

                        ! Ray volume counter
                        jray = jray + 1

                        ! Initialization of ray volumes
                        p_ray%z     (jc,jray,jb)    = zray
                        p_ray%dz    (jc,jray,jb)    = dzray
                        p_ray%lon   (jc,jray,jb)    = lonray
                        p_ray%dlon  (jc,jray,jb)    = dlonray
                        p_ray%lat   (jc,jray,jb)    = latray
                        p_ray%dlat  (jc,jray,jb)    = dlatray
                        p_ray%coslat(jc,jray,jb)    = COS(latray)
                        p_ray%m     (jc,jray,jb)    = mray
                        p_ray%dm    (jc,jray,jb)    = dmray
                        p_ray%k     (jc,jray,jb)    = kray
                        p_ray%dk    (jc,jray,jb)    = dkray
                        p_ray%l     (jc,jray,jb)    = lray
                        p_ray%dl    (jc,jray,jb)    = dlray
                        p_ray%wadens(jc,jray,jb)    = wa(jc,jk,jb)/dm_ini/dl_ini/dk_ini
                        p_ray%iexist(jc,jray,jb)    = jkmid
                        p_ray%specid(jc,jray,jb)    = 1    ! in idealized gwp case it 
                                                           ! must be kept as >= 0
                        p_ray%jk_active(jc,jray,jb) = nlev ! in idealized gwp case it 
                                                           ! should be kept as high 
                                                           ! (low in altitude) as possible

                      ENDDO ! jkk

                    ENDDO ! jll

                  ENDDO ! jmm

               ENDDO ! jlat

             ENDDO ! jlon

           ENDDO ! jz

        ENDDO ! jk

      ENDDO ! jc

    ENDDO ! jb
!!$OMP END DO
!!$OMP END PARALLEL

  CASE (planar_torus_geometry)

    ! Relocate wave packet center to the closest grid cell center
    !dx = 100000000._wp
    !dy = 100000000._wp
    !DO jb = i_startblk, i_endblk
    !  CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
    !                     i_startidx, i_endidx, i_rlstart, i_rlend)
    !  DO jc = i_startidx, i_endidx
    !    IF (ABS(x0-p_patch%cells%cartesian_center(jc,jb)%x(1)) <= dx) THEN 
    !      dx = ABS(x0-p_patch%cells%cartesian_center(jc,jb)%x(1))
    !      xs = p_patch%cells%cartesian_center(jc,jb)%x(1)
    !    ENDIF
    !    IF (ABS(y0-p_patch%cells%cartesian_center(jc,jb)%x(2)) <= dy) THEN 
    !      dy = ABS(y0-p_patch%cells%cartesian_center(jc,jb)%x(2))
    !      ys = p_patch%cells%cartesian_center(jc,jb)%x(2)
    !    ENDIF
    !  ENDDO ! jc
    !ENDDO ! jb
    !x0 = xs
    !y0 = ys

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jc,jk,jb,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jk = 1, nlev

        DO jc = i_startidx, i_endidx

          ! Cartesian vector of this locataion 
          vec_thispos = (/p_patch%cells%cartesian_center(jc,jb)%x(1),&
                          p_patch%cells%cartesian_center(jc,jb)%x(2),&
                          p_metrics%z_mc(jc,jk,jb)/)
          ! Cartesian vector of reference position
          vec_ref = (/x0,y0,z0/)

          ! Calculate distance in x
          dx = vec_thispos(1)-vec_ref(1)

          ! Calculate distance in y
          dy = vec_thispos(2)-vec_ref(2)

          ! Calculate distance in z
          dz = vec_thispos(3)-vec_ref(3)

          ! Intrinsic frequency
          omega = branch*SQRT((bvf2*(k_ini**2+l_ini**2) + &
                               p_patch%cells%f_c(jc,jb)**2*m_ini**2) / &
                              (k_ini**2 + l_ini**2 + m_ini**2))

          ! Buoyancy amplitude
          Bw = a0*bvf2/m_ini

          ! Buoyancy amplitude and wave action density at full levels
          ! CAUTION: consider that the envelope is a squared Gaussian / raised cosine
          ! as it is the envelope of the wave action rather than the bouyancy.
          IF (envtype==0) THEN ! no envelop
            env = 1._wp
          ELSEIF (envtype==1) THEN ! 1D Gaussian
            env = EXP(-dz**2/(sigmaz_wp**2))
          ELSEIF (envtype==2) THEN ! 1D Cosine
             IF (ABS(dz)<sigmaz_wp) THEN
               env = (1._wp+COS(pi*dz/sigmaz_wp))**2
             ELSE
               env = 0._wp
             ENDIF
          ELSEIF (envtype==3) THEN ! 2D Gaussian: symmetric in y
            env = EXP(-dy**2/(sigmay_wp**2) &
                      -dz**2/(sigmaz_wp**2))
          ELSEIF (envtype==4) THEN ! 3D Gaussian
            env = EXP(-dx**2/(sigmax_wp**2) &
                      -dy**2/(sigmay_wp**2) &
                      -dz**2/(sigmaz_wp**2))
          ELSE
            CALL finish('init_nh_isothermal_rest_atm_gwp', &
                        'envtype should be 1, 2, 3 or 4!')
          ENDIF
          wa(jc,jk,jb) = env*p_nh_prog%rho(jc,jk,jb)/2._wp * &
                        (1._wp+m_ini**2/(k_ini**2+l_ini**2)) * &
                         omega*Bw**2/bvf2**2

        ENDDO ! jc

      ENDDO ! jk

      DO jc = i_startidx, i_endidx

        ! Horizontal propagation test
        ! Horizontal localization of the source: skip calculations 
        ! if ltest_hprop=.true. and cell is out of the area defined 
        ! by xrmin, xrmax, yrmin, yrmax
        IF (ltest_hprop) THEN
          IF (p_patch%cells%cartesian_center(jc,jb)%x(1) > xrmax .OR. &
              p_patch%cells%cartesian_center(jc,jb)%x(1) < xrmin .OR. &
              p_patch%cells%cartesian_center(jc,jb)%x(2) > yrmax .OR. &
              p_patch%cells%cartesian_center(jc,jb)%x(2) < yrmin) THEN
            CYCLE
          ELSE
            ! Write coordinates to file
            OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
            WRITE(987,'(a,2x,F14.5,2x,F14.5)') &
                  'GW source x, y coordinates of cell centers:', &
                   p_patch%cells%cartesian_center(jc,jb)%x(1), &
                   p_patch%cells%cartesian_center(jc,jb)%x(2)
            CLOSE(987)
          ENDIF
        ENDIF

        ! Initialize x coordinate of ray volumes
        lonray_ref  = p_patch%cells%cartesian_center(jc,jb)%x(1)
        ! Correction for lonray_ref in cells at the domain borders in x 
        ! direction is necessary:
        ! For reasons of periodicity, at some of these cells the 
        ! cell center coordinate matches the coordinate of the other 
        ! side of the domain (domxmin, domxmax). In this case one has 
        ! to replace them. 
        vecx = (/p_gridinfo4ray%cellvertices_x(1,jc,jb), &
                 p_gridinfo4ray%cellvertices_x(2,jc,jb), &
                 p_gridinfo4ray%cellvertices_x(3,jc,jb)/)
        IF (MAXVAL(vecx) == domxmax .AND. MINVAL(vecx) == domxmin) THEN
          ! This happens when cells have vertices both with coordinates 
          ! domxmin and domxmax at the same time.
          DO jvec = 1,3
            IF (vecx(jvec) == domxmin) jmin = jvec
            IF (vecx(jvec) == domxmax) jmax = jvec
            IF (vecx(jvec) /= domxmin .AND. &
                vecx(jvec) /= domxmax) jthird = jvec
          ENDDO
          IF (ABS(vecx(jthird) - domxmin) < ABS(vecx(jthird) - domxmax)) THEN
            ! We are on the left side border of the domain
            lonray_ref = domxmin
          ELSE ! (ABS(vecx(jthird) - domxmin) > ABS(vecx(jthird) - domxmax))
            ! We are on the right side border of the domain
            lonray_ref = domxmax
          ENDIF
        ENDIF

        ! Calculate length of the line parallel to the triangle base 
        ! at the height of the center point (1/3h with h being 
        ! the triangle's height). This is the same as calclulating 
        ! the base of a smaller similar equilateral triangle of 
        ! height h_s = 2/3h. Then the "small base" of the smaller 
        ! triangle can be calculated as a_s = 2/sqrt(3)*hs = 4/(3*sqrt(3))*h
        base_small = 4._wp / 3._wp / SQRT(3._wp) * height_triangle

        ! Calculate min, max y of cell
        vecy = (/p_gridinfo4ray%cellvertices_y(1,jc,jb), &
                 p_gridinfo4ray%cellvertices_y(2,jc,jb), &
                 p_gridinfo4ray%cellvertices_y(3,jc,jb)/)
        IF (MAXVAL(vecy) == domymax .AND. MINVAL(vecy) == domymin) THEN
          ! This happens when cells have vertices both with coordinates 
          ! domymin and domymax at the same time.
          IF (ABS(p_patch%cells%cartesian_center(jc,jb)%x(2) - domymin) < &
              ABS(p_patch%cells%cartesian_center(jc,jb)%x(2) - domymax)) THEN
            ! We are on the front side border of the domain
            DO jvec = 1,3
              IF (vecy(jvec) == domymax) vecy(jvec) = domymin - height_triangle
            ENDDO
          ELSE
            ! We are on the back side border of the domain
            DO jvec = 1,3
              IF (vecy(jvec) == domymin) vecy(jvec) = domymax + height_triangle
            ENDDO
          ENDIF
        ENDIF
        miny = MINVAL(vecy)
        maxy = MAXVAL(vecy)

        ! Calculate length of the line parallel to the base 
        ! of which the square gives the area of the triangle cell.
        dlonray  = SQRT(p_patch%cells%area(jc,jb)) / REAL(nrlon,wp)
        dlatray  = SQRT(p_patch%cells%area(jc,jb)) / REAL(nrlat,wp)

        jray = 0

        DO jk = jkmax_ini_eff, jkmin_ini_eff, -1

           ! Loop over initial ray volumes within a layer
           DO jz = 1, nrz

             ! Loop over initial ray volumes in x direction
             DO jlon = 1, nrlon

               ! Loop over initial ray volumes in y direction
               DO jlat = 1, nrlat

                  ! Loop over initial ray volumes within dm_ini
                  DO jmm = 1, nrm

                    ! Loop over initial ray volumes within dl_ini
                    DO jll = 1, nrl

                      ! Loop over initial ray volumes within dk_ini
                      DO jkk = 1, nrk

                        ! Vertical extent of ray volume
                        dzray = (p_metrics%z_ifc(jc,jk-1,jb) - &
                                 p_metrics%z_ifc(jc,jk,  jb))/REAL(nrz,wp)

                        ! Altitude of ray volume center point
                        zray = p_metrics%z_ifc(jc,jk,jb) + (REAL(jz,wp)-0.5_wp) * dzray

                        ! y position of ray volume center point:
                        ! (height_triangle = maxy-miny)
                        latray  = miny + (jlat-0.5_wp) * height_triangle &
                                                       / REAL(nrlat,wp)

                        ! x position of ray volume center point:
                        lonray  = lonray_ref - 0.5_wp * base_small &
                                + (jlon-0.5_wp) * base_small / REAL(nrlon,wp)

                        ! Vertical wavenumber position of ray volume center point
                        mray = m_ini - dm_ini/2._wp + (REAL(jmm,wp)-0.5_wp)*dm_ini/REAL(nrm,wp)

                        ! Vertical wavenumber extent of ray volume
                        dmray = dm_ini/REAL(nrm,wp)

                        ! Meridional wavenumber position of ray volume center point
                        lray = l_ini - dl_ini/2._wp + (REAL(jll,wp)-0.5_wp)*dl_ini/REAL(nrl,wp)

                        ! Meridional wavenumber extent of ray volume
                        dlray = dl_ini/REAL(nrl,wp)

                        ! Latitudinal wavenumber position of ray volume center point
                        kray = k_ini - dk_ini/2._wp + (REAL(jkk,wp)-0.5_wp)*dk_ini/REAL(nrk,wp)

                        ! Latitudinal wavenumber extent of ray volume
                        dkray = dk_ini/REAL(nrk,wp)

                        ! Vertical index of ray volume center point
                        jkmid = jk  ! closest half level
!                        IF (jz > CEILING(nrz*0.5_wp))  jkmid = jk-1
                        IF (jz > nrz/2)  jkmid = jk-1

                        ! Ray volume counter
                        jray = jray + 1

                        ! Initialization of ray volumes
                        p_ray%z     (jc,jray,jb)    = zray
                        p_ray%dz    (jc,jray,jb)    = dzray
                        p_ray%lon   (jc,jray,jb)    = lonray
                        p_ray%dlon  (jc,jray,jb)    = dlonray
                        p_ray%lat   (jc,jray,jb)    = latray
                        p_ray%dlat  (jc,jray,jb)    = dlatray
                        p_ray%m     (jc,jray,jb)    = mray
                        p_ray%dm    (jc,jray,jb)    = dmray
                        p_ray%k     (jc,jray,jb)    = kray
                        p_ray%dk    (jc,jray,jb)    = dkray
                        p_ray%l     (jc,jray,jb)    = lray
                        p_ray%dl    (jc,jray,jb)    = dlray
                        p_ray%wadens(jc,jray,jb)    = wa(jc,jk,jb)/dm_ini/dl_ini/dk_ini
                        p_ray%iexist(jc,jray,jb)    = jkmid
                        p_ray%specid(jc,jray,jb)    = 1    ! in idealized gwp case it 
                                                           ! must be kept as >= 0
                        p_ray%jk_active(jc,jray,jb) = nlev ! in idealized gwp case it 
                                                           ! should be kept as high 
                                                           ! (low in altitude) as possible

                      ENDDO ! jkk

                    ENDDO ! jll

                  ENDDO ! jmm

               ENDDO ! jlat

             ENDDO ! jlon

           ENDDO ! jz

        ENDDO ! jk

      ENDDO ! jc

    ENDDO ! jb
!!$OMP END DO
!!$OMP END PARALLEL

  CASE DEFAULT

    CALL finish('init_nh_isothermal_rest_atm_gwp', &
         'Geometry not recognized! Should be spherical or torus!')

  END SELECT
  
  ! Synchronize GW field
  CALL sync_wave(p_ray,p_patch)

  END SUBROUTINE init_nh_isotherm_rest_atm_gwp
!--------------------------------------------------------------------
  END MODULE mo_nh_isotherm_rest_atm_gwp
