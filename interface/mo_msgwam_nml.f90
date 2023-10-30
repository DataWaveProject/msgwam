!>
!! Namelist for MS-GWaM
!!
!! @author Gergely Bölöni, Goethe Uni Frankfurt
!!
!! @par Revision History
!! Initial revision by Gergely Bölöni, Goethe Uni Frankfurt (2016-08-05)
!! 
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_msgwam_nml

  USE mo_msgwam_config,       ONLY: config_lmsgwam           => lmsgwam,            &
                                  & config_imethod_split     => imethod_split,      &
                                  & config_imethod_merge     => imethod_merge,      &
                                  & config_branch            => branch,             &
                                  & config_envtype           => envtype,            &
                                  & config_t0                => t0,                 &
                                  & config_p0                => p0,                 &
                                  & config_a0                => a0,                 &
                                  & config_sigmax_wp         => sigmax_wp,          &
                                  & config_sigmay_wp         => sigmay_wp,          &
                                  & config_sigmaz_wp         => sigmaz_wp,          &
                                  & config_lon0              => lon0,               &
                                  & config_lat0              => lat0,               &
                                  & config_x0                => x0,                 &
                                  & config_y0                => y0,                 &
                                  & config_z0                => z0,                 &
                                  & config_nrz               => nrz,                &
                                  & config_nrlon             => nrlon,              &
                                  & config_nrlat             => nrlat,              &
                                  & config_nrk               => nrk,                &
                                  & config_nrl               => nrl,                &
                                  & config_nrm               => nrm,                &
                                  & config_dk_ini            => dk_ini,             &
                                  & config_dl_ini            => dl_ini,             &
                                  & config_dm_ini            => dm_ini,             &
                                  & config_lambda_lon_ini    => lambda_lon_ini,     &
                                  & config_lambda_lat_ini    => lambda_lat_ini,     & 
                                  & config_lambda_z_ini      => lambda_z_ini,       &
                                  & config_k_ini             => k_ini,              &
                                  & config_l_ini             => l_ini,              &
                                  & config_m_ini             => m_ini,              &
                                  & config_omegamin          => omegamin,           &
                                  & config_omegamax          => omegamax,           &
                                  & config_cimin             => cimin,              &
                                  & config_cimax             => cimax,              &
                                  & config_flux_bw           => flux_bw,            &
                                  & config_flux_bs           => flux_bs,            &
                                  & config_lambda_z_min      => lambda_z_min,       &
                                  & config_darcraymin        => darcraymin,         &
                                  & config_darcraymax        => darcraymax,         &
                                  & config_dzraymin          => dzraymin,           &
                                  & config_dzraymax          => dzraymax,           &
                                  & config_lhsmooth          => lhsmooth,           &
                                  & config_lsmooth           => lsmooth,            &
                                  & config_lsmootht          => lsmootht,           &
                                  & config_lsmoothb          => lsmoothb,           &
                                  & config_nhsmooth          => nhsmooth,           &
                                  & config_nsmooth           => nsmooth,            &
                                  & config_nlevsmootht       => nlevsmootht,        &
                                  & config_lsteady           => lsteady,            &
                                  & config_lsat_mono         => lsat_mono,          &
                                  & config_ldiagprof         => ldiagprof,          &
                                  & config_ltest_hprop       => ltest_hprop,        &
                                  & config_ltest_gcircle     => ltest_gcircle,      &
                                  & config_l1ray             => l1ray,              &
                                  & config_lsaturation       => lsaturation,        &
                                  & config_llimittend        => llimittend,         &
                                  & config_lmsgwam_pmomflux  => lmsgwam_pmomflux,   &
                                  & config_lmsgwam_noforce   => lmsgwam_noforce,    &
                                  & config_lmsgwam_offline   => lmsgwam_offline,    &
                                  & config_lmsgwam_fricheat  => lmsgwam_fricheat,   &
                                  & config_lmvisc            => lmvisc,             &
                                  & config_lcorrlongwaves    => lcorrlongwaves,     &
                                  & config_lfluxlatsimple    => lfluxlatsimple,     &
                                  & config_lonrmin           => lonrmin,            &
                                  & config_lonrmax           => lonrmax,            &
                                  & config_latrmin           => latrmin,            &
                                  & config_latrmax           => latrmax,            &
                                  & config_xrmin             => xrmin,              &
                                  & config_xrmax             => xrmax,              &
                                  & config_yrmin             => yrmin,              &
                                  & config_yrmax             => yrmax,              &
                                  & config_zrmin             => zrmin,              &
                                  & config_zrmax             => zrmax,              & 
                                  & config_mylon             => mylon,              &
                                  & config_mylat             => mylat,              &
                                  & config_alpha_sat         => alpha_sat,          &
                                  & config_facgamma          => facgamma,           &
                                  & config_flux_tropics      => flux_tropics,       &
                                  & config_plaunch           => plaunch,            &
                                  & config_dt_add            => dt_add,             &
                                  & config_dt_gw_substep     => dt_gw_substep,      &
                                  & config_dh_factor_intpol  => dh_factor_intpol,   &
                                  & config_dz_crit_min       => dz_crit_min,        &
                                  & config_dz_factor_split   => dz_factor_split,    &
                                  & config_dh_factor_split   => dh_factor_split,    &
                                  & config_da_factor_split   => da_factor_split,    &
                                  & config_nstages           => nstages,            &
                                  & config_limfactor         => limfactor,          &
                                  & config_iazidim           => iazidim,            &
                                  & config_omegadim          => omegadim,           &
                                  & config_mdim              => mdim,               &
                                  & config_lrestart_from_nothing => lrestart_from_nothing, &
                                  & config_ltest_restart     => ltest_restart,       &
                                  & config_dir_restartfiles  => dir_restartfiles,    &
                                  & config_lvertinterpolation => lvertinterpolation, &
                                  & config_msgw_lower_bound   => msgw_lower_bound,   &
                                  & config_msgw_lower_bound_opt => msgw_lower_bound_opt, &
                                  & config_msgw_source_limit  => msgw_source_limit
  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: nnml
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_config,       ONLY: isRestart
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,   &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_impl_constants,      ONLY: max_dom, inwp
  USE mo_run_config,          ONLY: iforcing
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_exception,           ONLY: finish, message 
  USE mo_nh_testcases_nml,    ONLY: nh_test_name

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_msgwam_namelist

  !--------------------
  ! Namelist variables   
  !--------------------

  LOGICAL :: lmsgwam(max_dom)     ! logical switch for MS-GWaM (instead of default GWD scheme)
  LOGICAL :: lhsmooth             ! logical switch for horizontal smoothing GW fluxes
  LOGICAL :: lsmooth              ! logical switch for smoothing GW fluxes
  LOGICAL :: lsmootht             ! logical switch for smoothing GW tendencies
  LOGICAL :: lsmoothb             ! logical switch for smoothing background
  LOGICAL :: ldiagprof            ! logical switch for diagnostic profile output
  LOGICAL :: lsteady              ! logical switch for steady state version
  LOGICAL :: lsat_mono            ! logical switch for monochromatic saturation method
  LOGICAL :: ltest_hprop          ! logical switch for horizontal propagation test
  LOGICAL :: ltest_gcircle        ! logical switch for horizontal propagation test (great circle)
  LOGICAL :: l1ray                ! logical switch for horizontal propagation test (single ray)
  LOGICAL :: lsaturation          ! logical switch for wave breaking
  LOGICAL :: llimittend           ! logical switch for limiting tendencies
  LOGICAL :: lmsgwam_pmomflux     ! logical switch for using pseudomomentum fluxes
  LOGICAL :: lmsgwam_noforce      ! logical switch for no forcing
  LOGICAL :: lmsgwam_offline      ! logical switch for offline simulations
  LOGICAL :: lmsgwam_fricheat     ! logical switch for frictional heating by GW tendencies
  LOGICAL :: lmvisc               ! logical switch for molecular viscosity
  LOGICAL :: lcorrlongwaves       ! logical switch for psinc correction in dispersion relation
  LOGICAL :: lfluxlatsimple       ! logical switch for a tanh latitude profile of launch fluxes
  REAL(wp):: t0                   ! background temperature
  REAL(wp):: p0                   ! background surface pressure
  REAL(wp):: a0                   ! wave amplitude with respect to static instability
  REAL(wp):: sigmax_wp            ! half width of the wavepacket in x
  REAL(wp):: sigmay_wp            ! half width of the wavepacket in y
  REAL(wp):: sigmaz_wp            ! half width of the wavepacket in z
  REAL(wp):: lon0                 ! center of the wave packet in longitude (deg)
  REAL(wp):: lat0                 ! center of the wave packet in latitude  (deg)
  REAL(wp):: x0                   ! center of the wave packet in x (m)
  REAL(wp):: y0                   ! center of the wave packet in y (m)
  REAL(wp):: z0                   ! center of the wave packet in z (m)
  REAL(wp):: dk_ini               ! area around rays in wavenumber space (1/m)
  REAL(wp):: dl_ini               ! area around rays in wavenumber space (1/m)
  REAL(wp):: dm_ini               ! area around rays in wavenumber space (1/m)
  REAL(wp):: lambda_lon_ini       ! initial Llon of the wave packet (m)
  REAL(wp):: lambda_lat_ini       ! initial Llat of the wave packet (m)
  REAL(wp):: lambda_z_ini         ! initial Lz of the wave packet   (m)
  REAL(wp):: k_ini                ! initial k of the wave packet
  REAL(wp):: l_ini                ! initial l of the wave packet
  REAL(wp):: m_ini                ! initial m of the wave packet
  REAL(wp):: omegamin             ! minimum intrinsic frequency
  REAL(wp):: omegamax             ! maximum intrinsic frequency
  REAL(wp):: cimin                ! minimum horizontal phase speed in source
  REAL(wp):: cimax                ! maximum horizontal phase speed in source
  REAL(wp):: flux_bw(2)           ! min/max launch flux at boreal winter
  REAL(wp):: flux_bs(2)           ! min/max launch flux at boreal summer
  REAL(wp):: lambda_z_min         ! minimum vertical wavelength allowed (m)
  REAL(wp):: darcraymin           ! minimum ray size in the horizontal (rad)
  REAL(wp):: darcraymax           ! maximum ray size in the horizontal (rad)
  REAL(wp):: dzraymin             ! minimum ray size in z (m)
  REAL(wp):: dzraymax             ! maximum ray size in z (m)
  REAL(wp):: lonrmin              ! min longitude for bg source
  REAL(wp):: lonrmax              ! max longitude for bg source
  REAL(wp):: latrmin              ! min latitude for bg source
  REAL(wp):: latrmax              ! max latitude for bg source
  REAL(wp):: xrmin                ! min x for source
  REAL(wp):: xrmax                ! max x for source
  REAL(wp):: yrmin                ! min y for source
  REAL(wp):: yrmax                ! max y for source
  REAL(wp):: zrmin                ! min z for source
  REAL(wp):: zrmax                ! max z for source
  REAL(wp):: mylon                ! longitude coordinate for diag profile output
  REAL(wp):: mylat                ! latitude coordinate for diag profile output
  REAL(wp):: alpha_sat            ! tuning parameter for wave breaking
  REAL(wp):: facgamma             ! tuning parameter for gamma in dispersion relation
  REAL(wp):: dt_add               ! timestep for adding new rays
  REAL(wp):: dt_gw_substep        ! sub-timestep for gw propagation
  REAL(wp):: dh_factor_intpol     ! factor to the cell size dh (interpolation criterion)
  REAL(wp):: dz_crit_min          ! dz criterion minimum [m]   (splitting criterion)
  REAL(wp):: dz_factor_split      ! factor to the cell size dz (splitting criterion)
  REAL(wp):: dh_factor_split      ! factor to the cell size dh (splitting criterion)
  REAL(wp):: da_factor_split      ! factor to the cell area dA (splitting criterion)
  REAL(wp):: msgw_lower_bound     ! lower bound for discarding ray volumes propagating downwards
  REAL(wp):: msgw_source_limit    ! latitude limit for the gravity wave sources
  INTEGER :: imethod_split        ! option for split of ray volumes
  INTEGER :: imethod_merge        ! option for merge of ray volumes
  INTEGER :: envtype              ! envelope type 1:Gaussian, 2:Cosine
  INTEGER :: branch               ! frequency branch (disperion relation)
  INTEGER :: nrz                  ! no. of rays initialized within one vertical layer
  INTEGER :: nrlon                ! no. of rays initialized within one triangle in lon dir
  INTEGER :: nrlat                ! no. of rays initialized within one triangle in lat dir
  INTEGER :: nrk                  ! no. of rays initialized within dk_ini
  INTEGER :: nrl                  ! no. of rays initialized within dl_ini
  INTEGER :: nrm                  ! no. of rays initialized within dm_ini
  INTEGER :: nhsmooth             ! number of smoothing of fluxes in the horizontal
  INTEGER :: nsmooth              ! half number of points for vertical smoothing
  INTEGER :: nlevsmootht          ! level index above which tendency smoothing is applied
  INTEGER :: limfactor            ! limit factor for maxrays
  INTEGER :: iazidim              ! dimension of azimuthal angles
  INTEGER :: omegadim             ! dimension of intrinsic frequency spectral array
  INTEGER :: mdim                 ! dimension of vertical wavenumber spectral array
                                  ! in the Orr et al (2010) type source
  INTEGER :: nstages              ! 1: forward Euler, 4: Runge-Kutta 4th order
  INTEGER :: msgw_lower_bound_opt ! option for tendency distribution below the lower bound
  REAL(wp):: flux_tropics(4)      ! flux value at the equator, inflection latitude
                                  ! and width
  REAL(wp):: plaunch              ! launch-level pressure
  LOGICAL :: lrestart_from_nothing ! flag to restart without previous ray volumes
  LOGICAL :: ltest_restart        ! flag to test restart-file reading
  CHARACTER(LEN=1024) ::  &
    &  dir_restartfiles           ! directory to search/write restart files for ray volumes
  LOGICAL :: lvertinterpolation   ! flag fpr vertical interpolation of gradients for the propagation of rays


  NAMELIST /msgwam_nml/ lmsgwam, imethod_split, imethod_merge, t0, p0, a0, sigmax_wp, sigmay_wp, sigmaz_wp, lon0, lat0, x0, y0, z0, &
                        & dk_ini, dl_ini, dm_ini, lambda_lon_ini, lambda_lat_ini, lambda_z_ini, k_ini, &
                        & l_ini, m_ini, omegamin, omegamax, cimin, cimax, lambda_z_min, &
                        & darcraymin, darcraymax, dzraymin, dzraymax, &
                        & envtype, branch, nrz, nrlon, nrlat, nrk, nrl, nrm, &
                        & lhsmooth, lsmooth, lsmootht, lsmoothb, nhsmooth, nsmooth, nlevsmootht, &
                        & lsteady, lsat_mono, ldiagprof, lmsgwam_pmomflux, lvertinterpolation, &
                        & lsaturation, llimittend, lmsgwam_offline, lmsgwam_fricheat, &
                        & lmvisc, lcorrlongwaves, mylon, mylat, lonrmin, lonrmax, latrmin, latrmax, &
                        & xrmin, xrmax, yrmin, yrmax, zrmin, zrmax, msgw_lower_bound, msgw_lower_bound_opt, msgw_source_limit, &
                        & alpha_sat, facgamma, dt_add, dt_gw_substep, flux_tropics, plaunch, &
                        & dh_factor_intpol, dz_crit_min, dz_factor_split, dh_factor_split, da_factor_split, &
                        & nstages, limfactor, iazidim, omegadim, mdim, lfluxlatsimple, flux_bw, flux_bs, &
                        & lrestart_from_nothing, ltest_restart, dir_restartfiles, ltest_hprop, ltest_gcircle, l1ray, lmsgwam_noforce

CONTAINS
  !>
  !!
  SUBROUTINE read_msgwam_namelist( filename )

    ! IN/OUT variables
    CHARACTER(LEN=*), INTENT(IN) :: filename

    ! Local variables
    INTEGER :: ist, funit
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_msgwam_nml: read_msgwam_namelist'


    !----------------------------------------------------------------
    ! Default values
    !----------------------------------------------------------------

    lmsgwam(:) = .FALSE.      ! logical switch for MS-GWaM (instead of default GWD scheme)

    imethod_split = 1         ! option for split of ray volumes

    imethod_merge = 1         ! option for merge of ray volumes

    lambda_lon_ini = 1000._wp ! initial lon-dir wavelength of the wave packet (m)

    lambda_lat_ini = 1000._wp ! initial lat-dir wavelength of the wave packet (m)

    lambda_z_ini = 1000._wp   ! initial z-dir wavelength of the wave packet (m)

    lon0 =  0._wp             ! center of the wave packet in longitude (deg)

    lat0 = 45._wp             ! center of the wave packet in latitude (deg)

    x0 = 0._wp                ! center of the wave packet in x (m)

    y0 = 0._wp                ! center of the wave packet in y (m)

    z0 = 10000._wp            ! center of the wave packet in z (m)

    nrz = 1                   ! no. of rays initialized within one vertical layer

    nrlon = 1                 ! no. of rays initialized within one triangle in lon dir

    nrlat = 1                 ! no. of rays initialized within one triangle in lat dir

    dk_ini = 0.0001_wp        ! area around rays in wavenumber space (1/m)

    dl_ini = 0.0001_wp        ! area around rays in wavenumber space (1/m)

    dm_ini = 0.0001_wp        ! area around rays in wavenumber space (1/m)

    nrk = 1                   ! no. of rays initialized within dk_ini or dk

    nrl = 1                   ! no. of rays initialized within dl_ini or dl

    nrm = 1                   ! no. of rays initialized within dm_ini or dm

    envtype = 1               ! envelope type 1:Gaussian, 2:Cosine

    t0 = 240._wp              ! background temperature

    p0 = 100000.0_wp          ! background temperature (same as p0ref)

    a0 = 0.1_wp               ! amplitude of the wave packet wrt static instab

    sigmax_wp = 2000._wp      ! width of the wave packet in x (m)

    sigmay_wp = 2000._wp      ! width of the wave packet in y (m)

    sigmaz_wp = 2000._wp      ! width of the wave packet in z (m)

    branch = 1                ! frequency branch (dispersion relation)

    lhsmooth = .FALSE.        ! logical switch for horizontal smoothing GW fluxes

    lsmooth = .FALSE.         ! logical switch for smoothing GW fluxes

    lsmootht = .FALSE.        ! logical switch for smoothing GW tendencies

    lsmoothb = .FALSE.        ! logical switch for smoothing background

    nhsmooth =  2             ! number of smoothing of fluxes in the horizontal

    nsmooth =  1              ! half number of points for vertical smoothing

    nlevsmootht = 120         ! level index above which tendency smoothing is applied

    lsteady = .FALSE.         ! logical switch for steady state version

    lsat_mono = .FALSE.       ! logical switch for monochromatic saturation method

    ldiagprof = .TRUE.        ! logical switch for diagnostic profile output

    ltest_hprop = .FALSE.     ! logical switch for horizontal propagation test

    ltest_gcircle = .FALSE.   ! logical switch for horizontal propagation test (great circle)

    l1ray = .FALSE.           ! logical switch for horizontal propagation test (single ray)

    lsaturation = .FALSE.     ! logical switch for wave breaking

    llimittend = .FALSE.      ! logical switch for limiting tendencies

    lmsgwam_pmomflux = .FALSE.! logical switch for using pseudomomentum fluxes

    lmsgwam_noforce = .FALSE. ! logical switch for no forcing

    lmsgwam_offline = .FALSE. ! logical switch for offline simulations

    lmsgwam_fricheat = .FALSE.! logical switch for frictional heating by GW tendencies

    lmvisc = .FALSE.          ! logical switch for molecular viscosity

    lcorrlongwaves = .FALSE.  ! logical switch for psinc correction in dispersion relation

    lfluxlatsimple = .FALSE.  ! logical switch for a tanh latitude profile of launch fluxes

    lonrmin = -180._wp        ! min longitude for bg source

    lonrmax =  180._wp        ! max longitude for bg source

    latrmin = -90._wp         ! min latitude for bg source

    latrmax =  90._wp         ! max latitude for bg source

    xrmin = -1.e15_wp         ! min x for source

    xrmax =  1.e15_wp         ! max x for source

    yrmin = -1.e15_wp         ! min y for source

    yrmax =  1.e15_wp         ! max y for source

    zrmin = 2000._wp          ! min z for source

    zrmax = 75000._wp         ! max z for source

    mylon = 19._wp            ! longitude coordinate for diag profile output

    mylat = 47._wp            ! latitude coordinate for diag profile output

    alpha_sat = 1._wp         ! tuning parameter for wave breaking

    facgamma = 1._wp          ! tuning parameter for gamma in dispersion relation
                              ! psinc:     facgamma = 0.5_wp-2._wp/7._wp 
                              ! anelastic: facgamma = 0.5_wp
                              ! tuning:    facgamma > 0.5_wp e.g. 1._wp

    dt_add = 60._wp           ! timestep for adding new rays (s)

    dt_gw_substep = 30._wp    ! sub-timestep for gw propagation (s)

    nstages = 4               ! RK4 time stepping scheme

    limfactor = 1             ! limit factor for maxrays

    mdim = 20                 ! dimension of vertical wavenumber spectral array
                              ! in the Orr et al (2010) type source

    iazidim = 4               ! dimension of azimuthal angles

    omegadim = 1              ! dimension of intrinsic frequency spectral array

    omegamin = 1.0E-3_wp      ! minimum intrinsic frequency
    omegamax = 1.1E-3_wp      ! maximum intrinsic frequency

    cimin = 0.25              ! minimum horizontal phase speed in source (as in Orr et al)
    cimax = 100.              ! maximum horizontal phase speed in source (as in Orr et al)

    flux_bw = (/ 1.5E-3_wp, 2.5E-3_wp /) ! min/max launch flux at boreal winter
    flux_bs = (/ 1.5E-3_wp, 2.5E-3_wp /) ! min/max launch flux at boreal summer

    lambda_z_min = 1.         ! minimum horizontal wavelength allowed (m)

    darcraymin =  1.e-6       ! minimum ray size in the horizontal (rad)  (~6 m for Earth)
    darcraymax =  1.          ! maximum ray size in the horizontal (rad)  (1 rad ~ 57 deg)

    dzraymin =     1.         ! minimum ray size in z (m)
    dzraymax = 15000.         ! maximum ray size in z (m)

    flux_tropics = (/ -999._wp, 35._wp, 35._wp, 15._wp /)
                              ! flux value at the equator,
                              !   inflection latitudes (winter and summer hemispheres),
                              !   and width (-999: turn off this setting)
    plaunch = 45000.          ! launch-level pressure (Pa)

    ! Below for dh_factor_* [horizontal length-scale factors to SQRT(cell area)],
    ! some reference values are :
    ! 1.5 ~= 1x l_edge  ;  1.3 ~= 1x l_height
    ! in case of equilateral triangle cells

    dh_factor_intpol = 1.0_wp  ! factor to the cell size dh (interpolation criterion)

    dz_crit_min     = 1000._wp ! dz criterion minimum (splitting criterion)
    dz_factor_split = 2.5_wp   ! factor to the cell size dz (splitting criterion)
    dh_factor_split = 2.5_wp   ! factor to the cell size dh (splitting criterion)
    da_factor_split = 2.5_wp   ! factor to the cell area dA (splitting criterion)

    lrestart_from_nothing = .FALSE. ! flag to restart without previous ray volumes
    ltest_restart = .FALSE.   ! flag to test restart-file reading
    dir_restartfiles = ''

    lvertinterpolation = .False.

    msgw_lower_bound     = 1000._wp  ! set the lower bound to consider ray volumes
    msgw_lower_bound_opt = 0         ! the tendency is calculated as if fluxes reduce to zero
                                     ! 0: at the bound
                                     ! 1: linearly from the bound to the ground

    msgw_source_limit = 85         ! (absolute) latitude limit for the gravity wave sources

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above 
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (isRestart()) THEN
      funit = open_and_restore_namelist('msgwam_nml')
      READ(funit,NML=msgwam_nml)
      CALL close_tmpfile(funit)
    END IF    

    !---------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !---------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml('msgwam_nml',STATUS=ist)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, msgwam_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, msgwam_nml)   ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, msgwam_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Consistency check
    !-----------------------------------------------------

    IF ( ANY(lmsgwam) ) THEN
      ! inform that MS-GWaM has been switched on
      CALL message(TRIM(routine),'MS-GWaM switched on.')

      IF ( lambda_z_min == 0._wp )  &
        &  CALL finish(TRIM(routine), 'lambda_z_min zero!')

      IF (nh_test_name /= 'gwp' .AND. branch < 0) &
          CALL finish(TRIM(routine),  &
          &  'negative frequency branch not allowed with the Orr et al.'  &
          &  //' 2010 source scheme')

      IF (nh_test_name == 'gwp' .AND. lcorrlongwaves) THEN
          lcorrlongwaves = .FALSE.
          CALL message(TRIM(routine),'Test case gwp --> overwrite lcorrlongwaves to .FALSE.')
      ENDIF

      IF ( .NOT. ( nstages == 1 .OR. nstages == 4 ) )  &
        &  CALL finish(TRIM(routine),  &
        &     'Namelist variable nstages must be 4 (RK4) or 1 (forward Euler)!')

      IF (lsteady .AND. iazidim == 8) &
        & CALL finish(TRIM(routine), &
        &      'For lsteady=.true. iazidim > 4 is not implemented!')

      IF (iforcing /= 3) &
        & CALL finish(TRIM(routine), &
        &      '3D MS-GWaM implemented only for iforcing == 3!')


      IF (lhsmooth .AND. nhsmooth == 0)  lhsmooth = .FALSE.
      IF (lsmooth .AND. nsmooth == 0)  lsmooth = .FALSE.
      IF (lsmootht .AND. nsmooth == 0)  lsmootht = .FALSE.
      IF (lsmoothb .AND. nsmooth == 0)  lsmoothb = .FALSE.

      ! Horizontal propagation test of great circle propagation 
      ! (e.g. Hasha et al., 2008)
      IF (ltest_gcircle) ltest_hprop = .TRUE.

      ! Horizontal propagation test with a single ray
      IF (l1ray) THEN
        ltest_hprop = .TRUE.
        mdim      = 1
        omegadim  = 1
        lonrmin = 155.0
        lonrmax = 157.0
        latrmin = 45.0
        latrmax = 47.0
      ENDIF

      ! Horizontal propagation tests --> no forcing
      IF (ltest_hprop) THEN
        CALL message(TRIM(routine),'ltest_hprop=.T. '// &
        '--> lmsgwam_noforce=.T., lmsgwam_pmomflux=.T.')
        lmsgwam_noforce = .TRUE.
        lmsgwam_pmomflux = .TRUE.
      ENDIF

    ENDIF

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=msgwam_nml)
      CALL store_and_close_namelist(funit, 'msgwam_nml')
    ENDIF

    !-----------------------------------------------------
    ! Fill the configuration state
    !-----------------------------------------------------
    config_lmsgwam(:)        = lmsgwam(:)
    config_imethod_split     = imethod_split
    config_imethod_merge     = imethod_merge
    config_branch            = branch
    config_envtype           = envtype
    config_t0                = t0 
    config_p0                = p0 
    config_a0                = a0
    config_sigmax_wp         = sigmax_wp
    config_sigmay_wp         = sigmay_wp
    config_sigmaz_wp         = sigmaz_wp
    config_lon0              = lon0
    config_lat0              = lat0
    config_x0                = x0
    config_y0                = y0
    config_z0                = z0
    config_nrz               = nrz
    config_nrlon             = nrlon
    config_nrlat             = nrlat
    config_nrk               = nrk
    config_nrl               = nrl
    config_nrm               = nrm
    config_dk_ini            = dk_ini
    config_dl_ini            = dl_ini
    config_dm_ini            = dm_ini
    config_lambda_lon_ini    = lambda_lon_ini
    config_lambda_lat_ini    = lambda_lat_ini
    config_lambda_z_ini      = lambda_z_ini
    config_k_ini             = k_ini
    config_l_ini             = l_ini
    config_m_ini             = m_ini
    config_omegamin          = omegamin
    config_omegamax          = omegamax
    config_cimin             = cimin
    config_cimax             = cimax
    config_flux_bw(:)        = flux_bw(:)
    config_flux_bs(:)        = flux_bs(:)
    config_lambda_z_min      = lambda_z_min
    config_darcraymin        = darcraymin
    config_darcraymax        = darcraymax
    config_dzraymin          = dzraymin
    config_dzraymax          = dzraymax
    config_lhsmooth          = lhsmooth
    config_lsmooth           = lsmooth
    config_lsmootht          = lsmootht
    config_lsmoothb          = lsmoothb
    config_nhsmooth          = nhsmooth
    config_nsmooth           = nsmooth
    config_nlevsmootht       = nlevsmootht
    config_lsteady           = lsteady
    config_lsat_mono         = lsat_mono
    config_ldiagprof         = ldiagprof
    config_ltest_hprop       = ltest_hprop
    config_ltest_gcircle     = ltest_gcircle
    config_l1ray             = l1ray
    config_lsaturation       = lsaturation
    config_llimittend        = llimittend
    config_lmsgwam_pmomflux  = lmsgwam_pmomflux
    config_lmsgwam_noforce   = lmsgwam_noforce
    config_lmsgwam_offline   = lmsgwam_offline
    config_lmsgwam_fricheat  = lmsgwam_fricheat
    config_lmvisc            = lmvisc
    config_lcorrlongwaves    = lcorrlongwaves
    config_lfluxlatsimple    = lfluxlatsimple
    config_lonrmin           = lonrmin
    config_lonrmax           = lonrmax
    config_latrmin           = latrmin
    config_latrmax           = latrmax
    config_xrmin             = xrmin
    config_xrmax             = xrmax
    config_yrmin             = yrmin
    config_yrmax             = yrmax
    config_zrmin             = zrmin
    config_zrmax             = zrmax
    config_mylon             = mylon
    config_mylat             = mylat
    config_alpha_sat         = alpha_sat
    config_facgamma          = facgamma
    config_dt_add            = dt_add
    config_dt_gw_substep     = dt_gw_substep
    config_nstages           = nstages
    config_limfactor         = limfactor
    config_iazidim           = iazidim
    config_omegadim          = omegadim
    config_mdim              = mdim
    config_flux_tropics(:)   = flux_tropics(:)
    config_dh_factor_intpol  = dh_factor_intpol
    config_dz_crit_min       = dz_crit_min
    config_dz_factor_split   = dz_factor_split
    config_dh_factor_split   = dh_factor_split
    config_da_factor_split   = da_factor_split
    config_plaunch           = plaunch
    config_lrestart_from_nothing = ( lrestart_from_nothing .OR. lsteady )
    config_ltest_restart     = (/ ltest_restart, .FALSE. /)
    config_dir_restartfiles  = dir_restartfiles
    config_lvertinterpolation = lvertinterpolation
    config_msgw_lower_bound   = msgw_lower_bound
    config_msgw_lower_bound_opt = msgw_lower_bound_opt
    config_msgw_source_limit  = msgw_source_limit

  END SUBROUTINE read_msgwam_namelist

END MODULE mo_msgwam_nml
