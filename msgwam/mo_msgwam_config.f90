!>
!! configuration setup for MS-GWaM
!!
!! configuration setup for deep atmosphere
!!
!! @author Gergely Bölöni, Goethe Uni Frankfurt
!!
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
MODULE mo_msgwam_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: max_dom

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: lmsgwam
  PUBLIC :: imethod_split
  PUBLIC :: imethod_merge
  PUBLIC :: branch
  PUBLIC :: envtype
  PUBLIC :: t0
  PUBLIC :: p0
  PUBLIC :: a0
  PUBLIC :: sigmax_wp
  PUBLIC :: sigmay_wp
  PUBLIC :: sigmaz_wp
  PUBLIC :: lon0
  PUBLIC :: lat0
  PUBLIC :: x0
  PUBLIC :: y0
  PUBLIC :: z0
  PUBLIC :: nrz
  PUBLIC :: nrlon
  PUBLIC :: nrlat
  PUBLIC :: nrk
  PUBLIC :: nrl
  PUBLIC :: nrm
  PUBLIC :: dk_ini
  PUBLIC :: dl_ini
  PUBLIC :: dm_ini
  PUBLIC :: lambda_lon_ini
  PUBLIC :: lambda_lat_ini
  PUBLIC :: lambda_z_ini
  PUBLIC :: k_ini
  PUBLIC :: l_ini
  PUBLIC :: m_ini
  PUBLIC :: omegamin
  PUBLIC :: omegamax
  PUBLIC :: cimin
  PUBLIC :: cimax
  PUBLIC :: flux_bw
  PUBLIC :: flux_bs
  PUBLIC :: lambda_z_min
  PUBLIC :: darcraymin
  PUBLIC :: darcraymax
  PUBLIC :: dzraymin
  PUBLIC :: dzraymax
  PUBLIC :: lhsmooth
  PUBLIC :: lsmooth
  PUBLIC :: lsmootht
  PUBLIC :: lsmoothb
  PUBLIC :: nhsmooth
  PUBLIC :: nsmooth
  PUBLIC :: nlevsmootht
  PUBLIC :: lsteady
  PUBLIC :: lsat_mono
  PUBLIC :: ldiagprof
  PUBLIC :: ltest_hprop
  PUBLIC :: ltest_gcircle
  PUBLIC :: l1ray
  PUBLIC :: lsaturation
  PUBLIC :: llimittend
  PUBLIC :: lmsgwam_pmomflux
  PUBLIC :: lmsgwam_noforce
  PUBLIC :: lmsgwam_offline
  PUBLIC :: lmsgwam_fricheat
  PUBLIC :: lmvisc
  PUBLIC :: lcorrlongwaves
  PUBLIC :: lfluxlatsimple
  PUBLIC :: lonrmin
  PUBLIC :: lonrmax
  PUBLIC :: latrmin
  PUBLIC :: latrmax
  PUBLIC :: xrmin
  PUBLIC :: xrmax
  PUBLIC :: yrmin
  PUBLIC :: yrmax
  PUBLIC :: zrmin
  PUBLIC :: zrmax
  PUBLIC :: mylon
  PUBLIC :: mylat
  PUBLIC :: alpha_sat
  PUBLIC :: facgamma
  PUBLIC :: dt_add
  PUBLIC :: dt_gw_substep
  PUBLIC :: dh_factor_intpol
  PUBLIC :: dz_crit_min
  PUBLIC :: dz_factor_split
  PUBLIC :: dh_factor_split
  PUBLIC :: da_factor_split
  PUBLIC :: nstages
  PUBLIC :: limfactor
  PUBLIC :: iazidim
  PUBLIC :: omegadim
  PUBLIC :: mdim
  PUBLIC :: flux_tropics
  PUBLIC :: plaunch
  PUBLIC :: lvertinterpolation
  PUBLIC :: msgw_lower_bound
  PUBLIC :: msgw_lower_bound_opt
  PUBLIC :: msgw_source_limit
  PUBLIC :: nlev_lbnd

  PUBLIC :: nrays
  PUBLIC :: nrays_add
  PUBLIC :: maxrays
  PUBLIC :: nrays_bg
  PUBLIC :: nrays_add_bg
  PUBLIC :: maxrays_bg
  PUBLIC :: jray_offset_bg
  PUBLIC :: specid_offset_bg
  PUBLIC :: jklaunch_bg
  PUBLIC :: nrays_cv
  PUBLIC :: nrays_add_cv
  PUBLIC :: maxrays_cv
  PUBLIC :: jray_offset_cv
  PUBLIC :: specid_offset_cv
  PUBLIC :: nrays_coll

  PUBLIC :: lcalc_cgw_tend
  PUBLIC :: lcalc_flux_4dir_bg
  PUBLIC :: lcalc_flux_4dir_cv

  PUBLIC :: lrestart_from_nothing
  PUBLIC :: ltest_restart
  PUBLIC :: dir_restartfiles

  ! Namelist variables
  LOGICAL :: lmsgwam(max_dom)     ! logical switch for MS-GWaM (instead of default GWD scheme)
  LOGICAL :: lhsmooth             ! logical switch for horizontal smoothing GW fluxes
  LOGICAL :: lsmooth              ! logical switch for smoothing GW fluxes
  LOGICAL :: lsmootht             ! logical switch for smoothing GW tendencies
  LOGICAL :: lsmoothb             ! logical switch for smoothing background
  LOGICAL :: lsteady              ! logical switch for steady state version
  LOGICAL :: lsat_mono            ! logical switch for monochromatic saturation method
  LOGICAL :: ldiagprof            ! logical switch for diagnostic profile output
  LOGICAL :: ltest_hprop          ! logical switch for horizontal propagation test
  LOGICAL :: ltest_gcircle        ! logical switch for horizontal propagation test (great circle)
  LOGICAL :: l1ray                ! logical switch for horizontal propagation test (single ray)
  LOGICAL :: lsaturation          ! logical switch for wave breaking
  LOGICAL :: llimittend           ! logical switch for limiting tendencies
  LOGICAL :: lmsgwam_pmomflux     ! logical switch for using pseudomomentum fluxes
  LOGICAL :: lmsgwam_noforce      ! logical switch for no forcing
  LOGICAL :: lmsgwam_offline      ! logical switch for offline simulation
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
  REAL(wp):: lat0                 ! center of the wave packet in latitude (deg)
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
  REAL(wp):: omegamin             ! minimum horizontal phase speed in source
  REAL(wp):: omegamax             ! maximum horizontal phase speed in source
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
                                  !   and width
  REAL(wp):: plaunch              ! launch-level pressure (Pa)
  LOGICAL :: lvertinterpolation   ! flag fpr vertical interpolation of gradients for the propagation of rays
  REAL(wp):: msgw_lower_bound     ! lower bound for discarding ray volumes propagating downwards
  REAL(wp):: msgw_source_limit    ! (absolute) latitude limit for the gravity wave sources
  INTEGER :: nlev_lbnd(max_dom)   ! index for lower bound to discard in ray volumes

  ! Constants determined based on the above namelist values
  INTEGER :: nrays    (max_dom)         ! no. of ray volumes
  INTEGER :: nrays_add(max_dom)         ! no. of spectral bins of sources (total)
  INTEGER :: maxrays  (max_dom)         ! max no. of ray volumes
  INTEGER :: nrays_bg     (max_dom)     ! no. of ray volumes (background source)
  INTEGER :: nrays_add_bg  (max_dom)    ! no. of spectral bins of background source
  INTEGER :: maxrays_bg    (max_dom)    ! max no. of ray volumes (background source)
  INTEGER :: jray_offset_bg(max_dom)    ! offset to distinguish ray volume indices from different sources
                                        ! (background source)
  INTEGER :: specid_offset_bg(max_dom)  ! offset in spectral ID to distinguish different sources
  INTEGER :: jklaunch_bg   (max_dom)    ! launch-level index (background source)
  INTEGER :: nrays_cv      (max_dom)    ! no. of ray volumes (convective source)
  INTEGER :: nrays_add_cv  (max_dom)    ! no. of spectral bins of convective source
  INTEGER :: maxrays_cv    (max_dom)    ! max no. of ray volumes (convective source)
  INTEGER :: jray_offset_cv(max_dom)    ! offset to distinguish ray volume indices from different sources
                                        ! (convective source)
  INTEGER :: specid_offset_cv(max_dom)  ! offset in spectral ID to distinguish different sources
  INTEGER :: nrays_coll                 !
  LOGICAL :: lcalc_cgw_tend   (max_dom)  ! whether or not output convective GW tendencies
                                         ! (decision based on the output namelist)
  LOGICAL :: lcalc_flux_4dir_bg(max_dom) ! whether or not calculate 4 directional flux diagnostics
                                         ! (background source, decision based on the output namelist)
  LOGICAL :: lcalc_flux_4dir_cv(max_dom) ! whether or not calculate 4 directional flux diagnostics
                                         ! (convective source, decision based on the output namelist)

  LOGICAL :: lrestart_from_nothing      ! flag to restart without previous ray volumes
  LOGICAL :: ltest_restart(2)           ! flag to test restart-file reading
  CHARACTER(LEN=1024) :: dir_restartfiles

END MODULE mo_msgwam_config
