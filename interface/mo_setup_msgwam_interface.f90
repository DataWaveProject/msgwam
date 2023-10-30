!>
!! This module is an additional interface between nwp_gw_interface to the
!! gravity wave drag parametrization developed at the Goethe Uni:
!! Lagrangian WKB raytracer in phase space.
!!
!! @author Gergely Bölöni, Goethe Uni Frankfurt
!!
!! @par Revision History
!! Initial version by Gergely Bölöni, Goethe Uni Frankfurt (2016-03-02)
!! updates by Young-Ha Kim, Goethe Uni Frankfurt (2019-01-21)

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_setup_msgwam_interface

  USE mo_model_domain,         ONLY: t_patch
  USE mo_kind,                 ONLY: wp, vp, dp, sp
  USE mo_mpi,                  ONLY: get_my_mpi_work_id, get_my_mpi_work_comm_size, &
  &                                  my_process_is_work, p_barrier, p_comm_work
  USE mo_exception,            ONLY: message, finish, message_text, &
  &                                  debug_messages_on, debug_messages_off
  USE mo_vertical_coord_table, ONLY: vct_a
  USE mo_impl_constants,       ONLY: success, max_char_length, HINTP_TYPE_LONLAT_NNB, &
  &                                  min_rlcell_int, min_rledge_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,          ONLY: get_indices_c, get_indices_e
  USE mo_parallel_config,      ONLY: nproma
  USE mo_var_list,             ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,    ONLY: vlr_add, vlr_del
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_zaxis_type,           ONLY: ZA_SURFACE, ZA_REFERENCE, ZA_REFERENCE_HALF
  USE mo_var_groups,           ONLY: groups
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_cdi,                  ONLY: TSTEP_MIN, TSTEP_MAX, TSTEP_INSTANT, TSTEP_CONSTANT, &
  &                                  TSTEP_AVG, TSTEP_ACCUM, DATATYPE_PACK16, &
  &                                  DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d, t_ptr_i2d3d
  USE mo_delaunay_types,       ONLY: t_point
  USE mo_util_string,          ONLY: int2string
  USE mo_msgwam_config
  USE mo_name_list_output_config, ONLY: is_variable_in_output

  USE mo_dynamics_config,      ONLY: lcoriolis
  USE mo_math_constants,       ONLY: pi, pi2, pi_2, rad2deg, deg2rad
  USE mo_physical_constants,   ONLY: grav, rd, p0sl_bg, dtdz_standardatm
  USE mo_grid_config,          ONLY: grid_sphere_radius
  USE mo_grid_geometry_info,   ONLY: planar_torus_geometry, sphere_geometry

  USE mo_gw_source_config,     ONLY: gws_conv_config, max_nscale_cgw

  USE mtime,                   ONLY: datetime
  USE mo_nh_testcases_nml,     ONLY: nh_test_name

  IMPLICIT NONE

  PRIVATE

  INTEGER , PARAMETER :: triprc = sp   ! precision for trigonometric ftn calculation

  ! Define type for ray volumes
  TYPE t_ray

    INTEGER, POINTER  &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS     &
#endif
      ::         &
      iexist(:,:,:),      & ! existence / vertical level index of ray
                            ! iexist <= 0 --> ray does not exist
                            ! if iexist > 0 this integer is used for
                            ! storing the jk index of the ray volume
                            ! center point
      specid(:,:,:),      & ! spectral ID of ray in the launch spectrum
      jk_active(:,:,:),   & ! jk index above which the ray volume after
                            ! launch becomes active: prognostically moves
                            ! in phase-space and contributes to the mom-flux
                            ! calculation. This is always the same for
                            ! the background source (jklaunch) but it changes
                            ! every launch time for the convective source
                            ! where the launch height is based on the cloud top
      jr_last(:,:,:),     & ! last-launched ray index for a certain specid
      jk_source(:,:),     & ! jk index of the launch level (for the background
                            ! source same as jk_active)
      jk_full_rtop(:,:,:),& ! full level index closest to the top of the ray volume
      jk_full_rbot(:,:,:),& ! full level index closest to the bottom of the ray volume
      jk_half_rtop(:,:,:),& ! half level index closest to the top of the ray volume
      jk_half_rbot(:,:,:)   ! half level index closest to the bottom of the ray volume

    REAL(wp), POINTER  &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS     &
#endif
      ::         &
      lon   (:,:,:)         ,& ! ray position (longitude)          (rad)
      lat   (:,:,:)         ,& ! ray position (latitude)           (rad)
      z     (:,:,:)         ,& ! ray position (height)             (m)
      dlon  (:,:,:)         ,& ! ray size (longitude)              (rad)
      dlat  (:,:,:)         ,& ! ray size (latitude)               (rad)
      dz    (:,:,:)         ,& ! ray size (height)                 (m)
      coslat(:,:,:)         ,& ! cosine of ray latitude            (rad)
      k     (:,:,:)         ,& ! wavenumber along lon              (m^-1)
      l     (:,:,:)         ,& ! wavenumber along lat              (m^-1)
      m     (:,:,:)         ,& ! wavenumber along z                (m^-1)
      dk    (:,:,:)         ,& ! wavenumber spectrum along lon     (m^-1)
      dl    (:,:,:)         ,& ! wavenumber spectrum along lat     (m^-1)
      dm    (:,:,:)         ,& ! wavenumber spectrum along z       (m^-1)
      wadens(:,:,:)            ! wave action density (phase space) (m^3s^-1)
  END TYPE t_ray

  TYPE t_ray_coll

    INTEGER  ::   &
      specid     ,&
      jk_active

    REAL(wp) ::   &
      lon        ,&
      lat        ,&
      z          ,&
      dlon       ,&
      dlat       ,&
      dz         ,&
      coslat     ,&
      k          ,&
      l          ,&
      m          ,&
      dk         ,&
      dl         ,&
      dm         ,&
      wadens     ,&
      lon_min    ,&
      lat_min    ,&
      z_min

    LOGICAL  ::   &
      l_modif

  END TYPE t_ray_coll

  TYPE t_spl

    INTEGER , POINTER  &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS     &
#endif
      ::         &
      mode(:,:,:)       ! temporary array to be used in split routine

    INTEGER , ALLOCATABLE  &
      ::                   &
      jcol(:,:,:,:)     ! temporary array to be used in split routine

    REAL(wp), POINTER  &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS     &
#endif
      ::         &
      rot_c(:,:,:)  ,& ! temporary array to be used in split routine
      rot_s(:,:,:)  ,& ! temporary array to be used in split routine
      aux(:,:,:)       ! temporary array to be used in split routine

  END TYPE t_spl

  ! Define type for latitude dependent launch flux for the background source
  TYPE t_lfluxbg
    REAL(wp), POINTER  &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS     &
#endif
      ::         &
      lat_prof_bw(:,:)      ,& ! flux factor depending on latitude
      lat_prof_bs(:,:)      ,& ! flux factor depending on latitude
      lat_prof(:,:)            ! flux factor depending on latitude
  END TYPE t_lfluxbg

  TYPE t_msgwam

    TYPE(t_ptr_2d3d), ALLOCATABLE :: pars_cgw_ptr (:)
    TYPE(t_ptr_2d3d), ALLOCATABLE :: ctmfl_cgw_ptr(:)

    REAL(wp), POINTER             &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      &     , CONTIGUOUS          &
#endif
      &  ::                       &
      &  uwfl             (:,:,:) ,  & !! momentum flux from the ORR scheme             (Pa)
      &  vwfl             (:,:,:) ,  & !! momentum flux from the ORR scheme             (Pa)
      &  apmfl_mgm        (:,:,:) ,  & !! abs pseudomomentum flux from MS-GWaM          (Pa)
      &  amfl_mgm         (:,:,:) ,  & !! abs momentum flux from MS-GWaM                (Pa)
      &  uufl_mgm         (:,:,:) ,  & !! momentum flux from MS-GWaM                    (Pa)
      &  uvfl_mgm         (:,:,:) ,  & !! momentum flux from MS-GWaM                    (Pa)
      &  uwfl_mgm         (:,:,:) ,  & !! momentum flux from MS-GWaM                    (Pa)
      &  vufl_mgm         (:,:,:) ,  & !! momentum flux from MS-GWaM                    (Pa)
      &  vvfl_mgm         (:,:,:) ,  & !! momentum flux from MS-GWaM                    (Pa)
      &  vwfl_mgm         (:,:,:) ,  & !! momentum flux from MS-GWaM                    (Pa)
      &  mfl_mgm_e        (:,:,:) ,  & !! eastward momentum flux from MS-GWaM           (Pa)
      &  mfl_mgm_w        (:,:,:) ,  & !! westward momentum flux from MS-GWaM           (Pa)
      &  mfl_mgm_n        (:,:,:) ,  & !! northward momentum flux from MS-GWaM          (Pa)
      &  mfl_mgm_s        (:,:,:) ,  & !! southward momentum flux from MS-GWaM          (Pa)
      &  wafl_mgm_u       (:,:,:) ,  & !! wave action flux, (radially) upward           (Pa m)
      &  wafl_mgm_d       (:,:,:) ,  & !! wave action flux, (radially) downward         (Pa m)
      &  wafl_mgm_e       (:,:,:) ,  & !! wave action flux, eastward                    (Pa m)
      &  wafl_mgm_w       (:,:,:) ,  & !! wave action flux, westward                    (Pa m)
      &  wafl_mgm_n       (:,:,:) ,  & !! wave action flux, northward                   (Pa m)
      &  wafl_mgm_s       (:,:,:) ,  & !! wave action flux, southward                   (Pa m)
      &  aptfl_mgm        (:,:,:) ,  & !! abs potential temperature flux from MS-GWaM
      &  utfl_mgm         (:,:,:) ,  & !! potential temperature flux from MS-GWaM
      &  vtfl_mgm         (:,:,:) ,  & !! potential temperature flux from MS-GWaM
      &  ptfl_mgm_e       (:,:,:) ,  & !! eastward pot temperature flux from MS-GWaM
      &  ptfl_mgm_w       (:,:,:) ,  & !! westward pot temperature flux from MS-GWaM
      &  ptfl_mgm_n       (:,:,:) ,  & !! northward pot temperature flux from MS-GWaM
      &  ptfl_mgm_s       (:,:,:) ,  & !! southward pot temperature flux from MS-GWaM
      &  energy_mgm       (:,:,:) ,  & !! energy from MS-GWaM                           (kg/m/s2)
      &  energy_p_mgm     (:,:,:) ,  & !! potential energy from MS-GWaM                 (kg/m/s2)
      &  action_mgm       (:,:,:) ,  & !! action from MS-GWaM                           (kg/m/s)
!!!!! TEST !!!!! TO BE CLEANED LATER
      &  mfcxz_mgm     (:,:,:) ,  & ! test
      &  mfcyz_mgm     (:,:,:) ,  & ! test
      &  etx_mgm     (:,:,:) ,  & ! test
      &  ety_mgm     (:,:,:) ,  & ! test
!!!!! TEST !!!!! TO BE CLEANED LATER
      &  gwd_conv_u       (:,:,:) ,  & !! ZonalW-tendency from gravity wave drag (convective sources)
      &  gwd_conv_v       (:,:,:) ,  & !! MeridW-tendency from gravity wave drag (convective sources)
      &  gwd_conv_t       (:,:,:) ,  & !! PotTemp-tendency from gravity wave drag (convective sources)
      &  entr_cgw         (:,:,:) ,  & !! entropy profile for conv. GW parameterization
      &  heat_cgw         (:,:,:) ,  & !! heating profile for conv. GW parameterization
      &  tupd_cgw         (:,:,:) ,  & !! temperature of convective updraft parcel
      &  pars_cgw         (:,:,:) ,  & !! used parameters (output diagnosis)
      &  ctmfl_cgw      (:,:,:) ,  & !! cloud-top GW momentum flux, E,W,N,S-ward & absolute
      &  ctmfl_specc_cgw(:,:,:) ,  & !! momentum flux spectra, absolute
      &  ctmfl_spec0_cgw(:,:,:) ,  & !! momentum flux spectra, absolute
      &  ctmfl_spec1_cgw(:,:,:) ,  & !! momentum flux spectra, absolute
      &  apmfl_cgw      (:,:,:) ,  & !! pseudomomentum flux, absolute
      &  amfl_cgw       (:,:,:) ,  & !! momentum flux, absolute
      &  aptfl_cgw      (:,:,:) ,  & !! potential temperature flux, absolute
      &  mfl_cgw_e        (:,:,:) ,  & !! momentum flux, eastward
      &  mfl_cgw_w        (:,:,:) ,  & !! momentum flux, westward
      &  mfl_cgw_n        (:,:,:) ,  & !! momentum flux, northward
      &  mfl_cgw_s        (:,:,:) ,  & !! momentum flux, southward
      &  wafl_cgw_u       (:,:,:) ,  & !! wave action flux, (radially) upward           (Pa m)
      &  wafl_cgw_d       (:,:,:) ,  & !! wave action flux, (radially) downward         (Pa m)
      &  wafl_cgw_e       (:,:,:) ,  & !! wave action flux, eastward                    (Pa m)
      &  wafl_cgw_w       (:,:,:) ,  & !! wave action flux, westward                    (Pa m)
      &  wafl_cgw_n       (:,:,:) ,  & !! wave action flux, northward                   (Pa m)
      &  wafl_cgw_s       (:,:,:) ,  & !! wave action flux, southward                   (Pa m)
      &  ptfl_cgw_e       (:,:,:) ,  & !! theta flux, eastward
      &  ptfl_cgw_w       (:,:,:) ,  & !! theta flux, westward
      &  ptfl_cgw_n       (:,:,:) ,  & !! theta flux, northward
      &  ptfl_cgw_s       (:,:,:) ,  & !! theta flux, southward
      &  energy_cgw       (:,:,:) ,  & !! energy
      &  energy_p_cgw     (:,:,:) ,  & !! potential energy
      &  action_cgw       (:,:,:) ,  & !! action
      &  test_1_cgw       (:,:,:) ,  & !!
      &  test_2_cgw       (:,:,:) ,  & !!
      &  test_cgw         (:,:,:) ,  & !! temporary
      &  ll_dz_cgw        (:,:)   ,  &
      &  ml_dz_cgw        (:,:)

    INTEGER, POINTER              &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      &    , CONTIGUOUS           &
#endif
      &  ::                       &
      &  flag_cgw      (:,:)   ,  &
      &  mbas_con_cgw  (:,:)   ,  & !< cloud base level index
      &  mtop_con_cgw  (:,:)   ,  & !< cloud top  level index
      &  ktype_cgw     (:,:)   ,  & !< Type of convection
      &  ll_k_cgw      (:,:)   ,  &
      &  ml_k_cgw      (:,:)

    REAL(wp), POINTER            &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      &     , CONTIGUOUS          &
#endif
      &  ::                       &
      &  ddt_u_gwd_mgm (:,:,:) ,  & !! ZonalW-tendency from gravity wave drag
      &  ddt_v_gwd_mgm (:,:,:) ,  & !! MeridW-tendency from gravity wave drag
      &  ddt_t_gwd_mgm (:,:,:),   & !! Temp-tendency from gravity wave drag
!!!!! TEST !!!!! TO BE CLEANED LATER
      &  ddt_pt_gwd_mgm (:,:,:)      !! PotTemp-tendency from gravity wave drag
!!!!! TEST !!!!! TO BE CLEANED LATER

  END TYPE t_msgwam

  ! Define type for grid info for rays (based on Sebastian Borchert's idea)
  TYPE t_gridinfo4ray
    REAL(wp), ALLOCATABLE             &
      &  ::                           &
      &  cellvertices_lon(:,:,:),     &  ! (jvertex=1:3, jc=1:nproma, jb=1:nblks_c) longitudes
                                         ! of a given cell's 3 vertices
      &  cellvertices_lat(:,:,:),     &  ! (jvertex=1:3, jc=1:nproma, jb=1:nblks_c) latitudes
                                         ! of a given cell's 3 vertices
      &  cellvertices_x(:,:,:),       &  ! (jvertex=1:3, jc=1:nproma, jb=1:nblks_c) x (Cartesian) 
                                         ! of a given cell's 3 vertices
      &  cellvertices_y(:,:,:),       &  ! (jvertex=1:3, jc=1:nproma, jb=1:nblks_c) y (Cartesian)
                                         ! of a given cell's 3 vertices
      &  cellvertices_z(:,:,:)           ! (jvertex=1:3, jc=1:nproma, jb=1:nblks_c) y (Cartesian)
                                         ! of a given cell's 3 vertices

    ! For nstencil_nom = 3, the following information
    ! would be available via:
    ! * p_patch%cells%neighbor_idx(jc=1:nproma, jb=1:nblks_c, jneighbor=1:3)
    ! * p_patch%cells%neighbor_blk(jc=1:nproma, jb=1:nblks_c, jneighbor=1:3)
    INTEGER,  ALLOCATABLE             &
      & ::                            &
      & cellneighbors_nstencil(:,:),  &  ! (jc=1:nproma, jb=1:nblks_c) actual number of cells
                                         ! surrounding a given cell (might be 11
                                         ! instead of nstencil_nom = 12,
                                         ! if one of the given cell's vertices
                                         ! is one of the 12 pentagon points
                                         ! of the global icosahedral grid)
      & cellneighbors_idx(:,:,:),     &  ! (jneighbor=1:nstencil_nom, jc=1:nproma, jb=1:nblks_c)
                                         ! index jc of neighbor cell jneighbor
                                         ! of given cell (jc, jb)
      & cellneighbors_blk(:,:,:)         ! (jneighbor=1:nstencil_nom, jc=1:nproma, jb=1:nblks_c)
                                         ! index jb of neighbor cell jneighbor
                                         ! of given cell (jc, jb)
  END TYPE t_gridinfo4ray

  TYPE t_mgmgrid_nb    ! info of each cell and its neighbor cells

    INTEGER                       &
      &  ::                       &
      &  n_neighbor            ,  & !! no. of neighbor cells (12 or 11)
      &  jcb_neighbor(2,13)    ,  & !! (jc, jb)  of each cell and its neighbors
      &  jcol_c(13)                 !! 1-d index of ...

    REAL(triprc)                  &
      &  ::                       &
!cl   &  coslat(13)            ,  & !! cos(lat) of centers of each cell and its neighbors
      &  sinlat(13)            ,  & !! sin(lat) of ...
      &  x_cart(13)            ,  & !! x (Cartesian) of ...
      &  y_cart(13)                 !! y (Cartesian) of ...
                                    !! z = sinlat

    TYPE(t_point) :: v_cart(3,13)   !! (x,y,z) of vertices of each cell and its neighbors

  END TYPE t_mgmgrid_nb

  TYPE t_mgmgrid
    TYPE(t_mgmgrid_nb), ALLOCATABLE ::  &
      &  c_nb(:,:)                       !! cell's neighbor info (nproma,nblks_c)
    REAL(wp), ALLOCATABLE ::      &
      &  coslat(:,:)           ,  &
      &  darc_intpol(:,:)      ,  &
      &  darc_crit(:,:)        ,  &
      &  area_crit(:,:)        ,  &
      &  dz_crit(:)
    INTEGER  ::                   &
      &  ncol_int              ,  &      !! no. of interior grid columns
      &  ncol_halo                       !! no. of halo columns
  END TYPE t_mgmgrid

  ! Declare ray volumes
  TYPE(t_ray),          ALLOCATABLE :: p_ray(:)          ! total
  TYPE(t_var_list_ptr), ALLOCATABLE :: p_ray_list(:)     ! total
  TYPE(t_spl),          ALLOCATABLE :: p_spl(:)          ! total
  TYPE(t_var_list_ptr), ALLOCATABLE :: p_spl_list(:)     ! total
  TYPE(t_ray),          ALLOCATABLE :: p_rwork(:)        ! total
  TYPE(t_var_list_ptr), ALLOCATABLE :: p_rwork_list(:)   ! total
  TYPE(t_ray),          ALLOCATABLE :: p_ray_conv(:)     ! from convective sources
  TYPE(t_var_list_ptr), ALLOCATABLE :: p_ray_conv_list(:)! from convective sources
  TYPE(t_ray),          ALLOCATABLE :: p_ray_bg(:)       ! from bg sources
  TYPE(t_var_list_ptr), ALLOCATABLE :: p_ray_bg_list(:)  ! from bg sources
  TYPE(t_msgwam),       ALLOCATABLE :: p_msgwam(:)       ! field var (integrated over rays)
  TYPE(t_var_list_ptr), ALLOCATABLE :: p_msgwam_list(:)  ! field var (integrated over rays)
  ! Declare array to collect and manipulate ray volumes
  TYPE(t_ray_coll),     ALLOCATABLE :: ray_coll(:,:)
  ! Declare grid info for rays
  TYPE(t_gridinfo4ray), ALLOCATABLE :: p_gridinfo4ray(:) ! (1:n_dom)
  TYPE(t_mgmgrid),      ALLOCATABLE :: p_mgmgrid(:)    ! (1:n_dom)
  ! Declare launch flux
  TYPE(t_lfluxbg),      ALLOCATABLE :: p_lfluxbg(:)
  TYPE(t_var_list_ptr), ALLOCATABLE :: p_lfluxbg_list(:)
  ! Other declarations
  INTEGER                           :: jkmin_ini, jkmax_ini ! vertical level indices
                                                            ! corresponding to
                                                            ! zrmin, zrmax
  INTEGER                           :: nlaunch_max_cv       ! no. of ray volumes
                                                            ! launched at a time
                                                            ! (convective)
  INTEGER                           :: nlaunch_max_bg       ! no. of ray volumes
                                                            ! launched at a time
                                                            ! (background)
  INTEGER                           :: iout_msgwam, ndiag_msgwam, &
                                       jc_diag(9999), jb_diag(9999) ! variables for
                                                                    ! diagnostic outputs
  INTEGER,           PARAMETER      :: specid_merged    = 9999   ! offset in spectral ID to distinguish
  REAL(wp),          PARAMETER      :: dzlaunch_def = 1000._wp   ! initial vertical ray volume size
  REAL(wp),          PARAMETER      :: depth_ghost_cv = 3000._wp ! ghost layer depth (convective source)
  REAL(wp),          PARAMETER      :: depth_ghost_bg = 3000._wp ! ghost layer depth (background source)
  REAL(wp),          PARAMETER      :: min_coslat = 1.e-10  ! small enough (1-mm error)
                                                            ! and still treatable for cos^2
  REAL(wp),        ALLOCATABLE      :: omega_l(:)     ! intrinsic freq. at launch level
  REAL(wp),        ALLOCATABLE      :: domega_l(:)    ! spectral width for omegal
  REAL(wp),        ALLOCATABLE      :: sp_omega_l(:)  ! omega_l spectrum normalized
  REAL(wp)                          :: coef_st(4)     ! coeff to be used in the time scheme (RK4 or Euler)
  REAL(wp)                          :: domxmin        ! domain min on the torus in x
  REAL(wp)                          :: domxmax        ! domain max on the torus in x
  REAL(wp)                          :: domymin        ! domain min on the torus in y
  REAL(wp)                          :: domymax        ! domain max on the torus in y
  REAL(wp)                          :: edge_length_triangle ! edge length of triangles
  REAL(wp)                          :: height_triangle      ! height of triangles
  CHARACTER(len=128)                :: filename_diag1(9999), filename_diag2(9999) ! for diagnostic outputs
  CHARACTER(len=22)                 :: filename_ltest_hprop ! for horizontal propagation test with single ray

  INTEGER, PARAMETER :: nstencil_nom = 12  ! nominal number of cells surrounding a given cell
                                           ! defining its neighborhood (e.g. 3,
                                           ! 9 or 12)

  ! constants for restart I/O :  variable IDs (being unique, among 1-99)
  INTEGER, PARAMETER ::  iv_gidx =  1, iv_glon =  2, iv_glat =  3
  INTEGER, PARAMETER ::  iv_lon  =  4, iv_lat  =  5, iv_z    =  6
  INTEGER, PARAMETER ::  iv_dlon =  7, iv_dlat =  8, iv_dz   =  9
  INTEGER, PARAMETER ::  iv_k    = 10, iv_l    = 11, iv_m    = 12
  INTEGER, PARAMETER ::  iv_dk   = 13, iv_dl   = 14, iv_dm   = 15
  INTEGER, PARAMETER ::  iv_wad = 16
  INTEGER, PARAMETER ::  iv_iex = 17
  INTEGER, PARAMETER ::  iv_sid = 18
  INTEGER, PARAMETER ::  iv_kac = 19
  INTEGER, PARAMETER ::  iv_jrl = 91

  PUBLIC  :: setup_msgwam_interface
  PUBLIC  :: setup_msgwam
  PUBLIC  :: cleanup_msgwam
  PUBLIC  :: msgwam_read_restartfiles
  PUBLIC  :: msgwam_write_restartfiles
  PUBLIC  :: p_ray, t_ray
  PUBLIC  :: p_ray_list
  PUBLIC  :: p_spl, t_spl
  PUBLIC  :: p_spl_list
  PUBLIC  :: p_rwork
  PUBLIC  :: p_rwork_list
  PUBLIC  :: p_ray_conv, p_ray_conv_list
  PUBLIC  :: p_ray_bg, p_ray_bg_list
  PUBLIC  :: t_gridinfo4ray, p_gridinfo4ray
  PUBLIC  :: ray_coll, t_ray_coll
  PUBLIC  :: p_mgmgrid, t_mgmgrid, t_mgmgrid_nb, triprc
  PUBLIC  :: p_lfluxbg, p_lfluxbg_list
  PUBLIC  :: p_msgwam, t_msgwam
  PUBLIC  :: p_msgwam_list
  PUBLIC  :: jkmin_ini, jkmax_ini
  PUBLIC  :: dzlaunch_def
  PUBLIC  :: specid_merged
  PUBLIC  :: depth_ghost_cv
  PUBLIC  :: nlaunch_max_cv
  PUBLIC  :: depth_ghost_bg
  PUBLIC  :: nlaunch_max_bg
  PUBLIC  :: min_coslat
  PUBLIC  :: omega_l
  PUBLIC  :: domega_l
  PUBLIC  :: sp_omega_l
  PUBLIC  :: coef_st
  PUBLIC  :: domxmin
  PUBLIC  :: domxmax
  PUBLIC  :: domymin
  PUBLIC  :: domymax
  PUBLIC  :: edge_length_triangle
  PUBLIC  :: height_triangle
  PUBLIC  :: iout_msgwam
  PUBLIC  :: ndiag_msgwam, jc_diag, jb_diag
  PUBLIC  :: filename_diag1, filename_diag2
  PUBLIC  :: filename_ltest_hprop
  PUBLIC  :: nstencil_nom

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE setup_msgwam_interface(n_dom,p_patch,defcase)

  INTEGER,                INTENT(IN)   :: n_dom      ! number of model domains
  CHARACTER(LEN=*),       INTENT(IN)   :: defcase    ! construction or destruction?
  TYPE(t_patch),   TARGET,INTENT(IN)   :: p_patch(:) ! grid/patch info.

! Construct / destruct MS-GWaM
! Later use this ifdef
!#ifdef msgwam
    IF (TRIM(defcase) == 'construct') THEN
      CALL setup_msgwam(n_dom,p_patch(1:))
    END IF

    IF (TRIM(defcase) == 'destruct') THEN
      CALL cleanup_msgwam(n_dom)
    END IF
!#endif

END SUBROUTINE setup_msgwam_interface
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE setup_msgwam(n_dom,p_patch)
  INTEGER,                INTENT(IN) :: n_dom      ! number of model domains
  TYPE(t_patch),   TARGET,INTENT(IN) :: p_patch(:) ! grid/patch info.
  ! Local array bounds:
  INTEGER                        :: nlev, nlevp1, nblks_c, maxraysalloc
  INTEGER                        :: jk, jks, jg, ist
  INTEGER                        :: nlevp1_coll
  INTEGER                        :: rl_start, rl_end
  INTEGER                        :: i_startblk, i_endblk ! blocks
  INTEGER                        :: i_startidx, i_endidx ! slices
  INTEGER                        :: jc,jb
  INTEGER                        :: vidx(3), vblk(3)
  INTEGER                        :: cidx(3), cblk(3)
  INTEGER                        :: cidx_v(6,3), cblk_v(6,3)
  INTEGER                        :: jv, jec, jtri
  INTEGER                        :: istencil, jstencil
  INTEGER                        :: icn, ibn
  INTEGER                        :: nbcells_v
  INTEGER                        :: jcol_int, jcol_halo, jn, n_neighbor
  INTEGER                        :: jc_n, jb_n, jc21(21), jb21(21)
  INTEGER                        :: icheck
  INTEGER                        :: iomega
  INTEGER                        :: nsrc_cv
  REAL(wp)                       :: lat0_winter
  REAL(wp)                       :: lat0_summer
  REAL(wp)                       :: alat_deg, phi_lat
  REAL(wp)                       :: domega_work
  REAL(wp)                       :: mylon_rad, mylat_rad
  REAL(wp)                       :: distance, myvec(3), searchvec(3), radius
  REAL(wp)                       :: ttropo, ptropo, temp, zfull, pref
  REAL(wp)                       :: alpha
  REAL(wp), PARAMETER            :: pexp = 5._wp/3._wp   ! exponent p in Scinocca (2003)
                                                         ! 1 <= p <= 2
  REAL(wp), PARAMETER            :: one_o_6 = 1._wp/6._wp
  REAL(wp), PARAMETER            :: rk4st_c(4) = (/0.5_wp,0.5_wp,1.0_wp,one_o_6/)
  ! Reference atmosphere parameters
  REAL(wp), PARAMETER            :: htropo = 11000._wp   ! [m]    tropopause height
  REAL(wp), PARAMETER            :: t00    = 288.15_wp   ! [m]    temperature at sea level
  REAL(wp), PARAMETER            :: tp          = 0._wp  ! turning point for tanh
  REAL(wp), PARAMETER            :: sm          = 11._wp ! smoothness for tanh
  REAL(wp), PARAMETER            :: shiftalpha  = 0.5_wp ! center of alpha function
  LOGICAL                        :: lsimloc
  CHARACTER(len=28)              :: str_latlon
  CHARACTER(len=max_char_length) :: listname

  IF ( (nh_test_name /= 'gwp' .OR. (.NOT. ltest_hprop)) .AND. .NOT. lcoriolis )  &
    & CALL finish('setup_msgwam',  &
    & 'lcoriolis should be TRUE if the Orr et al. 2010 source scheme is used.')

  ! Allocate fields for tendency and diagnostic variables,
  ! regardless of lmsgwam
  ALLOCATE (p_msgwam(n_dom), p_msgwam_list(n_dom), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam', &
    &                            'MS-GWaM: allocation of MS-GWaM field variables failed')
  DO jg = 1, n_dom
    nblks_c = p_patch(jg)%nblks_c
    nlev    = p_patch(jg)%nlev
    WRITE(listname,'(a,i2.2)') 'p_msgwam_list_in_D',jg
    CALL new_field_list(jg, nlev, nblks_c, TRIM(listname), p_msgwam_list(jg), p_msgwam(jg))
  ENDDO

  ! As this setup is always called we make sure that no
  ! calculations are done if MS-GWaM will not be used on
  ! any subdomains
  IF (ALL( .NOT. lmsgwam(:) ))  RETURN

  ! Initialize counter for diagnostic output (ldiagprof=.T.)
  iout_msgwam = 1

  ! Allocate WKB ray type
  ALLOCATE (p_ray(n_dom), p_ray_list(n_dom), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam',&
      &          'MS-GWaM: allocation of WKB ray type failed')

  ALLOCATE (p_spl(n_dom), p_spl_list(n_dom), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam',&
      &          'MS-GWaM: allocation of WKB spl type failed')

  ALLOCATE (p_rwork(n_dom), p_rwork_list(n_dom), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam',&
      &          'MS-GWaM: allocation of WKB ray work type failed')

  IF (ANY( gws_conv_config%n_source(:) > 0 )) THEN
    ! Allocate WKB ray type (conv)
    ALLOCATE (p_ray_conv(n_dom), p_ray_conv_list(n_dom), STAT=ist)
    IF (ist/=success) CALL finish ('setup_msgwam',&
        &          'MS-GWaM: allocation of WKB ray type failed (conv)')
  END IF

  ! Allocate WKB ray type (bg)
  ALLOCATE (p_ray_bg(n_dom), p_ray_bg_list(n_dom), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam',&
      &          'MS-GWaM: allocation of WKB ray type failed (bg)')

  ! Allocate grid info for rays
  ALLOCATE (p_gridinfo4ray(n_dom), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam',&
      &          'MS-GWaM: allocation of gridinfo4ray failed')

  ALLOCATE (p_mgmgrid(n_dom), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam',&
      &          'MS-GWaM: allocation of mgmgrid failed')

  ! Allocate grid-dependent variable
  ALLOCATE (p_lfluxbg(n_dom), p_lfluxbg_list(n_dom), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam',&
      &          'MS-GWaM: allocation of WKB ray grid failed')

  ALLOCATE (omega_l(omegadim), domega_l(omegadim), sp_omega_l(omegadim), STAT=ist)
  IF (ist/=success) CALL finish ('setup_msgwam',&
      &          'MS-GWaM: allocation of omega_l failed')

  nlaunch_max_cv = CEILING(depth_ghost_cv/gws_conv_config%dz_launch - 0.5_wp)
  nlaunch_max_bg = CEILING(depth_ghost_bg/dzlaunch_def - 0.5_wp)

  ! Set coef according to the ray propagator time integration scheme (RK4 or Euler)
  IF (nstages == 4) THEN
    coef_st(1:4) = rk4st_c(:)
  ELSE IF (nstages == 1) THEN
    coef_st(1) = 1._wp
  END IF

  ! Define spectral elements for omega. Due to the hydrostatic assumption in
  ! Scinocca (2003) one should keep omega such as bvf >> abs(omega) >> fc.
  ! Discretization: equidistant
  domega_work = (omegamax-omegamin)/REAL(omegadim,wp)
  DO iomega = 1, omegadim
    omega_l(iomega) = omegamin + (REAL(iomega,wp) - 0.5_wp)*domega_work
  ENDDO
  domega_l(:) = domega_work

  DO iomega = 1, omegadim
    sp_omega_l(iomega) = omega_l(iomega)**(1._wp-pexp)
  ENDDO
  domega_work = SUM(sp_omega_l(:)*domega_l(:))
  sp_omega_l(:) = sp_omega_l(:)/domega_work

  ! Total number of ray volumes : zero for lmsgwam(jg) == .F.
  nrays(:) = 0

  DO jg = 1, n_dom

    ! Cycle domain if MS-GWaM will not be used for the current one
    IF ( .NOT. lmsgwam(jg) )  CYCLE

    ! Determine size of arrays
    nblks_c = p_patch(jg)%nblks_c

    ! Number of vertical levels
    nlev   = p_patch(jg)%nlev
    nlevp1 = p_patch(jg)%nlevp1

    ! Determine the launch level for constant-source (based on Orr et al. 2010) waves: jklaunch_bg
    ! (imported from atm_phy_nwp/mo_nwp_phy_init.f90 and atm_phy_schemes/data_gwd.f90)
    jklaunch_bg(jg) = 0
    ttropo = t00 + dtdz_standardatm*htropo
    ptropo = p0sl_bg*(ttropo/t00)**(-grav/(rd*dtdz_standardatm))
    DO jk = nlev, 2, -1
      jks = jk + p_patch(jg)%nshift_total
      zfull = 0.5_wp*(vct_a(jks) + vct_a(jks+1))
      IF (zfull < htropo) THEN
        temp = t00 + dtdz_standardatm*zfull
        pref = p0sl_bg*(temp/t00)**(-grav/(rd*dtdz_standardatm))
      ELSE
        pref = ptropo*EXP(-grav*(zfull-htropo)/(rd*ttropo))
      ENDIF
      IF (pref > plaunch)  jklaunch_bg(jg) = jk
    ENDDO
    IF (jklaunch_bg(jg) == 0)  CALL finish('setup_msgwam',  &
      &                                   'MS-GWaM: invalid plaunch [Pa]')

    ! find index of lower bound for consideration of ray volumes
    ! downward propagating rays below will be discarded
    ! note that nlev_lbnd has a minimum = 1
    nlev_lbnd(jg) = 1
    DO jk = nlev, 1, -1
      jks = jk + p_patch(jg)%nshift_total
      IF (vct_a(jks) <= msgw_lower_bound) THEN
        nlev_lbnd(jg) = nlevp1 - jk
      ELSE
        EXIT
      END IF
    ENDDO

    ! Number of rays to be added if ltimeadd=.true.
    IF (nh_test_name == 'gwp') THEN

      ! Check if zrmin and zrmax is within the vertical domain
      IF (zrmin < vct_a(nlev+p_patch(jg)%nshift_total) .OR. &
          zrmax > vct_a(1+p_patch(jg)%nshift_total)) THEN
        CALL finish('setup_msgwam', &
             'zrmin or zrmax lies out of the vertical domain!')
      ENDIF

      ! Find layer indices in between rays are initialized
      DO jk = nlev, 1, -1
        jks = jk + p_patch(jg)%nshift_total
        IF (vct_a(jks) <= zrmin) jkmax_ini = jk
        IF (vct_a(jks) <= zrmax) jkmin_ini = jk
      ENDDO

      nrays_add_bg(jg) = (jkmax_ini-jkmin_ini+1)*nrz*nrlon*nrlat*nrm*nrl*nrk

    ELSE ! for the time being source based on Orr et al. 2010

      IF ( lfluxlatsimple .OR. flux_tropics(1) >= 0. ) THEN
        nrays_add_bg(jg) = iazidim*omegadim*mdim
      ELSE
        ! i.e. you can switch this source off by setting lfluxlatsimple = .F. and
        ! flux_tropics(1) < 0. (default for both) in the namelist
        nrays_add_bg(jg) = 0
      ENDIF

    ENDIF

    IF ( gws_conv_config%n_source(jg) > 0 ) THEN
      nsrc_cv = gws_conv_config%n_source(jg)
      nrays_add_cv(jg) = SUM(gws_conv_config%nm(1:nsrc_cv)*gws_conv_config%nphi(1:nsrc_cv))
    ELSE
      ! i.e. you can switch this source off by setting n_source = 0
      ! (default) in the namelist
      nrays_add_cv(jg) = 0
    END IF

    nrays_add(jg) = nrays_add_bg(jg) + nrays_add_cv(jg)  ! + nrays_add_??(jg)

    specid_offset_bg(jg) = 0
    specid_offset_cv(jg) = nrays_add_bg(jg)
!   specid_offset_??(jg) = specid_offset_cv(jg)

    ! Set maximum number of ray volumes per column. The key parameter is
    ! limfactor / limfactor_conv to be set in the namelist.
    IF (.NOT. lsteady) THEN
      ! Transient case
      maxrays_bg(jg) = nrays_add_bg(jg) * limfactor
      maxrays_cv(jg) = nrays_add_cv(jg) * gws_conv_config%limfactor_conv
      ! Up to two volumes can be launched at once (for each spectral element)
      nrays_bg(jg) = maxrays_bg(jg) + nrays_add_bg(jg)*nlaunch_max_bg
      nrays_cv(jg) = maxrays_cv(jg) + nrays_add_cv(jg)*nlaunch_max_cv
    ELSE
      ! Steady-state case
      maxrays_bg(jg) = nrays_add_bg(jg)
      maxrays_cv(jg) = nrays_add_cv(jg)
      nrays_bg  (jg) = nrays_add_bg(jg)
      nrays_cv  (jg) = nrays_add_cv(jg)
    ENDIF
    maxrays(jg) = maxrays_bg(jg) + maxrays_cv(jg)
    ! Note: maxrays and maxrays_?? are not used any more (except in output
    ! messages.
    nrays  (jg) = nrays_bg  (jg) + nrays_cv  (jg)

    IF (nrays(jg) == 0)  CALL finish ('setup_msgwam',  &
      &                               'MS-GWaM: no spectral bin for GW sources')

    ! Offset values to distinguish between bg and conv sources
    jray_offset_cv(jg) = 0
    jray_offset_bg(jg) = nrays_cv(jg)

    ! Allocation size for rays
    maxraysalloc = nrays(jg)

    IF ( gws_conv_config%n_source(jg) > 0 ) THEN
      WRITE(message_text,'(a,i8)') 'Initial no. of rays (conv):', nrays_add_cv(jg)
      CALL message('setup_msgwam', TRIM(message_text))
    END IF

    WRITE(message_text,'(a,i8)') 'Initial no. of rays (background):', nrays_add_bg(jg)
    CALL message('setup_msgwam', TRIM(message_text))

    WRITE(message_text,'(a,i8)') 'Maximum no. of rays allowed:', maxrays(jg)
    CALL message('setup_msgwam', TRIM(message_text))

    WRITE(message_text,'(a,i8)') 'Ray array allocated with size:', maxraysalloc
    CALL message('setup_msgwam', TRIM(message_text))

    WRITE(listname,'(a,i2.2)') 'p_ray_list_in_D',jg
    IF (.NOT. lsteady) THEN
      CALL new_ray_list(jg, maxraysalloc, nrays_add(jg), nblks_c, TRIM(listname),  &
        &  p_ray_list(jg), p_ray(jg), 'ray')
    ELSE
      CALL new_ray_list(jg, maxraysalloc, 0, 1, TRIM(listname), p_ray_list(jg), p_ray(jg), 'spec')
    END IF

    IF (.NOT. lsteady) THEN
      WRITE(listname,'(a,i2.2)') 'p_spl_list_in_D',jg
      CALL new_spl_list(jg, maxraysalloc, nblks_c, TRIM(listname),  &
        &  p_spl_list(jg), p_spl(jg), 'spl')
    END IF

    IF (.NOT. lsteady) THEN
      WRITE(listname,'(a,i2.2)') 'p_rwork_list_in_D',jg
      CALL new_ray_list(jg, maxraysalloc, 0, nblks_c, TRIM(listname),  &
        &  p_rwork_list(jg), p_rwork(jg), 'rwork')
    END IF

    IF ( gws_conv_config%n_source(jg) > 0 ) THEN
      WRITE(listname,'(a,i2.2)') 'p_ray_conv_list_in_D',jg
      CALL new_ray_list(jg, nrays_add_cv(jg), 0, nblks_c, TRIM(listname),  &
        &  p_ray_conv_list(jg), p_ray_conv(jg), 'src_conv')
    END IF

    IF ( nrays_add_bg(jg) /= 0 ) THEN
      WRITE(listname,'(a,i2.2)') 'p_ray_bg_list_in_D',jg
      CALL new_ray_list(jg, nrays_add_bg(jg), 0, nblks_c, TRIM(listname),  &
        &  p_ray_bg_list(jg), p_ray_bg(jg), 'src_bg')
    END IF

    WRITE(listname,'(a,i2.2)') 'p_lfluxbg_list_in_D',jg
    CALL new_lfluxbg_list(jg, nblks_c, TRIM(listname), p_lfluxbg_list(jg), p_lfluxbg(jg))

    ! Initialize ray volumes. None of them "exists" initially
    p_ray(jg)%iexist(:,:,:) = 0
    p_ray(jg)%lon(:,:,:)    = 0.0_wp
    p_ray(jg)%lat(:,:,:)    = 0.0_wp
    p_ray(jg)%z(:,:,:)      = 0.0_wp
    p_ray(jg)%dlon(:,:,:)   = 0.0_wp
    p_ray(jg)%dlat(:,:,:)   = 0.0_wp
    p_ray(jg)%dz(:,:,:)     = 0.0_wp
    p_ray(jg)%coslat(:,:,:) = 1.0_wp
    p_ray(jg)%k(:,:,:)      = 0.0_wp
    p_ray(jg)%dk(:,:,:)     = 0.0_wp
    p_ray(jg)%l(:,:,:)      = 0.0_wp
    p_ray(jg)%dl(:,:,:)     = 0.0_wp
    p_ray(jg)%m(:,:,:)      = 0.0_wp
    p_ray(jg)%dm(:,:,:)     = 0.0_wp
    p_ray(jg)%wadens(:,:,:) = 0.0_wp
    IF (.NOT. lsteady) THEN
      p_ray(jg)%specid(:,:,:)       = 0
      p_ray(jg)%jk_active(:,:,:)    = 0
      p_ray(jg)%jr_last(:,:,:)      = 0
      p_ray(jg)%jk_full_rtop(:,:,:) = 0
      p_ray(jg)%jk_full_rbot(:,:,:) = 0
      p_ray(jg)%jk_half_rtop(:,:,:) = 0
      p_ray(jg)%jk_half_rbot(:,:,:) = 0
      !
      p_spl(jg)%rot_c(:,:,:) = 1.
      p_spl(jg)%rot_s(:,:,:) = 0.
      p_spl(jg)%aux  (:,:,:) = 0.
      !
      p_rwork(jg)%iexist(:,:,:)    = 0
      p_rwork(jg)%specid(:,:,:)    = 0
      p_rwork(jg)%jk_active(:,:,:) = 0
      p_rwork(jg)%lon(:,:,:)    = 0.
      p_rwork(jg)%lat(:,:,:)    = 0.
      p_rwork(jg)%z(:,:,:)      = 0.
      p_rwork(jg)%dlon(:,:,:)   = 0.
      p_rwork(jg)%dlat(:,:,:)   = 0.
      p_rwork(jg)%dz(:,:,:)     = 0.
      p_rwork(jg)%coslat(:,:,:) = 1.
      p_rwork(jg)%k(:,:,:)      = 0.
      p_rwork(jg)%dk(:,:,:)     = 0.
      p_rwork(jg)%l(:,:,:)      = 0.
      p_rwork(jg)%dl(:,:,:)     = 0.
      p_rwork(jg)%m(:,:,:)      = 0.
      p_rwork(jg)%dm(:,:,:)     = 0.
      p_rwork(jg)%wadens(:,:,:) = 0.
    END IF

    ! Initialize conv rays
    IF ( gws_conv_config%n_source(jg) > 0 ) THEN
      p_ray_conv(jg)%jk_source(:,:) = 0
      p_ray_conv(jg)%wadens(:,:,:) = 0.0_wp
    END IF

    ! Initialize latitude dependent flux profile
    p_lfluxbg(jg)%lat_prof(:,:) = 1.0_wp
    p_lfluxbg(jg)%lat_prof_bw(:,:) = 1.0_wp
    p_lfluxbg(jg)%lat_prof_bs(:,:) = 1.0_wp

    ! initialize tendencies
    p_msgwam(jg)%ddt_u_gwd_mgm(:,:,:) = 0.0_wp
    p_msgwam(jg)%ddt_v_gwd_mgm(:,:,:) = 0.0_wp
    p_msgwam(jg)%ddt_t_gwd_mgm(:,:,:) = 0.0_wp
    p_msgwam(jg)%ddt_pt_gwd_mgm(:,:,:) = 0.0_wp
    !
    p_msgwam(jg)%mfcxz_mgm(:,:,:) = 0._wp
    p_msgwam(jg)%mfcyz_mgm(:,:,:) = 0._wp
    p_msgwam(jg)%etx_mgm  (:,:,:) = 0._wp
    p_msgwam(jg)%ety_mgm  (:,:,:) = 0._wp

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch(jg)%cells%start_block(rl_start)
    i_endblk   = p_patch(jg)%cells%end_block(rl_end)

    ! Background launch flux profile with smooth transition between bw and bs
    ! profiles
    IF (lfluxlatsimple) THEN

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          alat_deg = p_patch(jg)%cells%center(jc,jb)%lat*rad2deg

          ! Turning point in latitude degrees
          alpha = shiftalpha+shiftalpha*TANH((alat_deg-tp)/sm)
          ! Boreal winter
          p_lfluxbg(jg)%lat_prof_bw(jc,jb) = (1._wp-alpha)*flux_bw(1) + alpha*flux_bw(2)
          ! Boreal summer
          p_lfluxbg(jg)%lat_prof_bs(jc,jb) = (1._wp-alpha)*flux_bs(1) + alpha*flux_bs(2)

          ! linearly reduce the wave amplitudes from 100% to 0 in the 5 degrees below msgw_source_limit
          IF (msgw_source_limit <= 90.) THEN
            p_lfluxbg(jg)%lat_prof_bw(jc, jb) = p_lfluxbg(jg)%lat_prof_bw(jc, jb) &
                  * min(1._wp, max(0._wp, 1. - (abs(alat_deg) - msgw_source_limit + 5.) / 5.))
            p_lfluxbg(jg)%lat_prof_bs(jc, jb) = p_lfluxbg(jg)%lat_prof_bs(jc, jb) &
                  * min(1._wp, max(0._wp, 1. - (abs(alat_deg) - msgw_source_limit + 5.) / 5.))
          END IF

        ENDDO ! jc

      ENDDO ! jb

    ELSE IF (flux_tropics(1) >= 0.) THEN

      lat0_winter = flux_tropics(2)
      lat0_summer = flux_tropics(3)

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          alat_deg = p_patch(jg)%cells%center(jc,jb)%lat*rad2deg

          ! Boreal winter profile
          IF ( alat_deg >= lat0_winter + flux_tropics(4) ) THEN
            p_lfluxbg(jg)%lat_prof_bw(jc,jb) = flux_bw(2)
          ELSE IF ( alat_deg <= -(lat0_summer + flux_tropics(4)) ) THEN
            p_lfluxbg(jg)%lat_prof_bw(jc,jb) = flux_bw(1)
          ELSE IF ( alat_deg > lat0_winter - flux_tropics(4) ) THEN
            phi_lat = MAX(-0.5_wp,MIN(0.5_wp,  &
              &       (ABS(alat_deg) - lat0_winter)/(2._wp*flux_tropics(4))))
            p_lfluxbg(jg)%lat_prof_bw(jc,jb) = 0.5_wp  &
              &  *(  (flux_bw(2) + flux_tropics(1))     &
              &    + (flux_bw(2) - flux_tropics(1))*SIN(phi_lat*pi) )
          ELSE IF ( alat_deg < -(lat0_summer - flux_tropics(4)) ) THEN
            phi_lat = MAX(-0.5_wp,MIN(0.5_wp,  &
              &       (ABS(alat_deg) - lat0_summer)/(2._wp*flux_tropics(4)) ))
            p_lfluxbg(jg)%lat_prof_bw(jc,jb) = 0.5_wp  &
              &  *(  (flux_bw(1) + flux_tropics(1))     &
              &    + (flux_bw(1) - flux_tropics(1))*SIN(phi_lat*pi) )
          ELSE
            p_lfluxbg(jg)%lat_prof_bw(jc,jb) = flux_tropics(1)
          END IF

          ! Boreal summer profile
          IF ( alat_deg >= lat0_summer + flux_tropics(4) ) THEN
            p_lfluxbg(jg)%lat_prof_bs(jc,jb) = flux_bs(2)
          ELSE IF ( alat_deg <= -(lat0_winter + flux_tropics(4)) ) THEN
            p_lfluxbg(jg)%lat_prof_bs(jc,jb) = flux_bs(1)
          ELSE IF ( alat_deg > lat0_summer - flux_tropics(4) ) THEN
            phi_lat = MAX(-0.5_wp,MIN(0.5_wp,  &
              &       (ABS(alat_deg) - lat0_summer)/(2._wp*flux_tropics(4))))
            p_lfluxbg(jg)%lat_prof_bs(jc,jb) = 0.5_wp  &
              &  *(  (flux_bs(2) + flux_tropics(1))     &
              &    + (flux_bs(2) - flux_tropics(1))*SIN(phi_lat*pi) )
          ELSE IF ( alat_deg < -(lat0_winter - flux_tropics(4)) ) THEN
            phi_lat = MAX(-0.5_wp,MIN(0.5_wp,  &
              &       (ABS(alat_deg) - lat0_winter)/(2._wp*flux_tropics(4)) ))
            p_lfluxbg(jg)%lat_prof_bs(jc,jb) = 0.5_wp  &
              &  *(  (flux_bs(1) + flux_tropics(1))     &
              &    + (flux_bs(1) - flux_tropics(1))*SIN(phi_lat*pi) )
          ELSE
            p_lfluxbg(jg)%lat_prof_bs(jc,jb) = flux_tropics(1)
          END IF

          ! linearly reduce the wave amplitudes from 100% to 0 in the 5 degrees below msgw_source_limit
          IF (msgw_source_limit <= 90.) THEN
            p_lfluxbg(jg)%lat_prof_bw(jc, jb) = p_lfluxbg(jg)%lat_prof_bw(jc, jb) &
                  * min(1._wp, max(0._wp, 1. - (abs(alat_deg) - msgw_source_limit + 5.) / 5.))
            p_lfluxbg(jg)%lat_prof_bs(jc, jb) = p_lfluxbg(jg)%lat_prof_bs(jc, jb) &
                  * min(1._wp, max(0._wp, 1. - (abs(alat_deg) - msgw_source_limit + 5.) / 5.))
          END IF

        ENDDO ! jc

      ENDDO ! jb

!   ELSE  ! another setting can be made here.

    ENDIF ! lfluxlatsimple, flux_tropics(1)

  ENDDO ! jg

  ! Allocate array to collect and manipulate ray volumes
! nrays_coll = MAXVAL(CEILING( REAL(nrays(1:n_dom)*39)  &
!   &                         /REAL(MAX(1,p_patch(1:n_dom)%nlevp1)) ))
  nrays_coll = MAXVAL(nrays(1:n_dom))*39
  nlevp1_coll = MAXVAL(p_patch(1:n_dom)%nlevp1, MASK=(nrays(1:n_dom) > 0))
  ALLOCATE( ray_coll(nrays_coll,nlevp1_coll), STAT=ist )
  IF (ist /= success) CALL finish('setup_msgwam', 'MS-GWaM: allocation of ray_coll failed')

  ndiag_msgwam = 0

  ! Diagnostic output for selected single grid columns
  IF ( ldiagprof .AND. lmsgwam(1) ) THEN

    jg = 1   ! currently, not yet coded for n_dom > 1 for the diag printing

    ! Are we on the Sphere or on the Torus?
    SELECT CASE(p_patch(jg)%geometry_info%geometry_type)

    CASE(sphere_geometry)

      jc_diag(:) = -999  ;  jb_diag(:) = -999

      ! Exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch(jg)%cells%start_block(rl_start)
      i_endblk   = p_patch(jg)%cells%end_block(rl_end)

      CALL debug_messages_on

      ! Convert to radians
      mylon_rad = mylon * deg2rad
      mylat_rad = mylat * deg2rad

      myvec = (/COS(mylat_rad)*COS(mylon_rad),&
                COS(mylat_rad)*SIN(mylon_rad),&
                SIN(mylat_rad)/)

      DIAGLOOP1: DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          searchvec = (/COS(p_patch(jg)%cells%center(jc,jb)%lat)*COS(p_patch(jg)%cells%center(jc,jb)%lon),&
                        COS(p_patch(jg)%cells%center(jc,jb)%lat)*SIN(p_patch(jg)%cells%center(jc,jb)%lon),&
                        SIN(p_patch(jg)%cells%center(jc,jb)%lat)/)

          ! Calculate great-circle distance based on chord length
          distance = 2.0_wp*ASIN(0.5_wp*NORM2(searchvec-myvec))*grid_sphere_radius

          !radius = 180.e3_wp  ! 180 km
          radius = 300.e3_wp  ! 300 km

          ! Ensure that the point of interest is at the same hemisphere
          lsimloc = (mylat_rad*p_patch(jg)%cells%center(jc,jb)%lat > 0._wp .AND. &
                     mylon_rad*p_patch(jg)%cells%center(jc,jb)%lon > 0._wp)

          IF (distance <= radius .AND. lsimloc) THEN

            ndiag_msgwam = ndiag_msgwam + 1

            jc_diag(ndiag_msgwam) = jc  ;  jb_diag(ndiag_msgwam) = jb

            ! Include lat/lon in filename to avoid overwriting when multiple PEs found such locations..
            WRITE(str_latlon,'(f0.1,a,f0.1)') p_patch(jg)%cells%center(jc,jb)%lat * rad2deg, &
              &  '_', p_patch(jg)%cells%center(jc,jb)%lon * rad2deg
            WRITE(filename_diag1(ndiag_msgwam),'(a)') 'raytracer_'//TRIM(str_latlon)//'.dat'
            WRITE(filename_diag2(ndiag_msgwam),'(a)') 'fields_'   //TRIM(str_latlon)//'.dat'

            ! Print mylat, mylon, my_jc, my_jb
            WRITE(message_text,'(a,2F10.3)') 'ldiagprof=.T. ==> print MS-GWaM diagnostic profiles'// &
              &   ' to fields.dat at lat, lon:', mylat_rad*rad2deg, mylon_rad*rad2deg
            CALL message('', TRIM(message_text))
            WRITE(message_text,'(a,2F10.3)') 'closest grid point found at lat, lon:', &
              &   p_patch(jg)%cells%center(jc,jb)%lat*rad2deg, &
              &   p_patch(jg)%cells%center(jc,jb)%lon*rad2deg
            CALL message('', TRIM(message_text))

            ! Print coriolis parameter
            WRITE(message_text,'(a,2i4,E12.4)') 'jb, jc, Coriolis parameter:', &
              &   jb, jc, p_patch(jg)%cells%f_c(jc,jb)
            CALL message('', TRIM(message_text))

            IF (ndiag_msgwam == 9999) THEN
              CALL message('', 'CAUTION: ndiag_msgwam reaches its maximum: stop searching')
              EXIT DIAGLOOP1
            ENDIF

          ENDIF

        ENDDO ! jc

      ENDDO DIAGLOOP1 ! jb

      CALL debug_messages_off

    CASE(planar_torus_geometry)

      jc_diag(:) = -999  ;  jb_diag(:) = -999

      ! Exclude boundary interpolation zone of nested domains
      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int

      i_startblk = p_patch(jg)%cells%start_block(rl_start)
      i_endblk   = p_patch(jg)%cells%end_block(rl_end)

      CALL debug_messages_on

      ! Convert to radians
      mylon_rad = mylon
      mylat_rad = mylat

      myvec = (/mylon,mylat,0._wp/)

      DIAGLOOP2: DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          & i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          searchvec = (/p_patch(jg)%cells%cartesian_center(jc,jb)%x(1),&
                        p_patch(jg)%cells%cartesian_center(jc,jb)%x(2),&
                        p_patch(jg)%cells%cartesian_center(jc,jb)%x(3) /)

          ! Calculate great-circle distance based on chord length
          distance = SQRT((myvec(1)-searchvec(1))**2 &
                   + (myvec(2)-searchvec(2))**2)

          !radius = 180.e3_wp  ! 180 km
          radius = 300.e3_wp  ! 300 km

          IF (distance <= radius) THEN

            ndiag_msgwam = ndiag_msgwam + 1

            jc_diag(ndiag_msgwam) = jc  ;  jb_diag(ndiag_msgwam) = jb

            ! Include lat/lon in filename to avoid overwriting when multiple PEs found such locations..
            WRITE(str_latlon,'(F12.2,a,F12.2)') p_patch(jg)%cells%cartesian_center(jc,jb)%x(2)/1000._wp, & ! km
              &                            '_', p_patch(jg)%cells%cartesian_center(jc,jb)%x(1)/1000._wp    ! km
            WRITE(filename_diag1(ndiag_msgwam),'(a)') 'raytracer_'//TRIM(str_latlon)//'.dat'
            WRITE(filename_diag2(ndiag_msgwam),'(a)') 'fields_'   //TRIM(str_latlon)//'.dat'

            ! Print mylat, mylon, my_jc, my_jb
            WRITE(message_text,'(a,2F7.3)') 'ldiagprof=.T. ==> print MS-GWaM diagnostic profiles'// &
              &   ' to fields.dat at x, y:', mylon/1000._wp, mylat/1000._wp  ! km
            CALL message('', TRIM(message_text))
            WRITE(message_text,'(a,2F12.2)') 'closest grid point found at x, y:', &
              &   p_patch(jg)%cells%cartesian_center(jc,jb)%x(1)/1000._wp, & !  km
              &   p_patch(jg)%cells%cartesian_center(jc,jb)%x(2)/1000._wp    !  km
            CALL message('', TRIM(message_text))

            ! Print coriolis parameter
            WRITE(message_text,'(a,2i4,E12.4)') 'jb, jc, Coriolis parameter:', &
              &   jb, jc, p_patch(jg)%cells%f_c(jc,jb)
            CALL message('', TRIM(message_text))

            IF (ndiag_msgwam == 9999) THEN
              CALL message('', 'CAUTION: ndiag_msgwam reaches its maximum: stop searching')
              EXIT DIAGLOOP2
            ENDIF

          ENDIF

        ENDDO ! jc

      ENDDO DIAGLOOP2 ! jb

      CALL debug_messages_off

    CASE DEFAULT

      CALL finish('setup_msgwam', &
           'Geometry not recognized! Should be spherical or torus!')

    END SELECT

  ENDIF

  ! Decide whether gw fluxes are put out for
  ! all 4 directions (E, W, N, S) separately
  lcalc_cgw_tend   (:) = .FALSE.
  lcalc_flux_4dir_cv(:) = .FALSE.
  lcalc_flux_4dir_bg(:) = .FALSE.

  DO jg = 1, n_dom

    IF (.NOT. lmsgwam(jg))  CYCLE

    IF ( gws_conv_config%n_source(jg) > 0 ) THEN

      lcalc_cgw_tend(jg) =  &
        &  is_variable_in_output(var_name='gwd_conv_u') .OR.  &
        &  is_variable_in_output(var_name='gwd_conv_v') .OR.  &
        &  is_variable_in_output(var_name='gwd_conv_t')

      lcalc_flux_4dir_cv(jg) =  &
        &  is_variable_in_output(var_name='mfl_cgw_e')  .OR.  &
        &  is_variable_in_output(var_name='mfl_cgw_w')  .OR.  &
        &  is_variable_in_output(var_name='mfl_cgw_n')  .OR.  &
        &  is_variable_in_output(var_name='mfl_cgw_s')  .OR.  &
        &  is_variable_in_output(var_name='wafl_cgw_u')  .OR.  &
        &  is_variable_in_output(var_name='wafl_cgw_d')  .OR.  &
        &  is_variable_in_output(var_name='wafl_cgw_e')  .OR.  &
        &  is_variable_in_output(var_name='wafl_cgw_w')  .OR.  &
        &  is_variable_in_output(var_name='wafl_cgw_n')  .OR.  &
        &  is_variable_in_output(var_name='wafl_cgw_s')  .OR.  &
        &  is_variable_in_output(var_name='ptfl_cgw_e') .OR.  &
        &  is_variable_in_output(var_name='ptfl_cgw_w') .OR.  &
        &  is_variable_in_output(var_name='ptfl_cgw_n') .OR.  &
        &  is_variable_in_output(var_name='ptfl_cgw_s')
    END IF

    IF ( nrays_add_bg(jg) /= 0 ) THEN
      lcalc_flux_4dir_bg(jg) =  &
        &  is_variable_in_output(var_name='mfl_mgm_e')  .OR.  &
        &  is_variable_in_output(var_name='mfl_mgm_w')  .OR.  &
        &  is_variable_in_output(var_name='mfl_mgm_n')  .OR.  &
        &  is_variable_in_output(var_name='mfl_mgm_s')  .OR.  &
        &  is_variable_in_output(var_name='wafl_mgm_u')  .OR.  &
        &  is_variable_in_output(var_name='wafl_mgm_d')  .OR.  &
        &  is_variable_in_output(var_name='wafl_mgm_e')  .OR.  &
        &  is_variable_in_output(var_name='wafl_mgm_w')  .OR.  &
        &  is_variable_in_output(var_name='wafl_mgm_n')  .OR.  &
        &  is_variable_in_output(var_name='wafl_mgm_s')  .OR.  &
        &  is_variable_in_output(var_name='ptfl_mgm_e') .OR.  &
        &  is_variable_in_output(var_name='ptfl_mgm_w') .OR.  &
        &  is_variable_in_output(var_name='ptfl_mgm_n') .OR.  &
        &  is_variable_in_output(var_name='ptfl_mgm_s')
    END IF

    ! Determine size of arrays
    nblks_c = p_patch(jg)%nblks_c

    ! Number of vertical levels
    nlevp1 = p_patch(jg)%nlevp1

    ! Setup of grid info for rays
    ALLOCATE( p_mgmgrid(jg)%c_nb       (nproma,nblks_c),  &
      &       p_mgmgrid(jg)%coslat     (nproma,nblks_c),  &
      &       p_mgmgrid(jg)%darc_intpol(nproma,nblks_c),  &
      &       p_mgmgrid(jg)%darc_crit  (nproma,nblks_c),  &
      &       p_mgmgrid(jg)%area_crit  (nproma,nblks_c),  &
      &       p_mgmgrid(jg)%dz_crit(nlevp1), STAT=ist )
    IF (ist /= success) CALL finish('setup_msgwam', &
      &          'MS-GWaM: allocation of p_mgmgrid elements failed')

    ! Initialize grid info for ray
    ALLOCATE( p_gridinfo4ray(jg)%cellvertices_lon      ( 3, nproma, nblks_c ), &
      &       p_gridinfo4ray(jg)%cellvertices_lat      ( 3, nproma, nblks_c ), &
      &       p_gridinfo4ray(jg)%cellvertices_x        ( 3, nproma, nblks_c ), &
      &       p_gridinfo4ray(jg)%cellvertices_y        ( 3, nproma, nblks_c ), &
      &       p_gridinfo4ray(jg)%cellvertices_z        ( 3, nproma, nblks_c ), &
      &       p_gridinfo4ray(jg)%cellneighbors_nstencil( nproma, nblks_c ),    &
      &       p_gridinfo4ray(jg)%cellneighbors_idx     ( nstencil_nom, nproma, nblks_c ), &
      &       p_gridinfo4ray(jg)%cellneighbors_blk     ( nstencil_nom, nproma, nblks_c ), &
      &       STAT=ist)
    IF (ist/=success) CALL finish ('setup_msgwam',&
        &          'MS-GWaM: allocation of gridinfo4ray elements failed')

!YHK+
    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int - 2
    i_startblk = p_patch(jg)%cells%start_block(rl_start)
    i_endblk   = p_patch(jg)%cells%end_block(rl_end)

    ! p_mgmgrid(jg)%c_nb(jc,jb)%XXXX_c each have a size of 13.
    ! In XXXX_c(1), the cell's own info is saved, while
    ! in XXXX_c(2:13), its neighbor cells' info are saved.
    ! For some cells, the number of neighbors is 11 instead of 12 (around the vertices
    ! of the original icosahedron) or even less (around the outer halo lines).
    ! In these cases, the remaining portion of XXXX_c(:) is filled with a special value
    ! (see below).

    jcol_int = 0   ;   jcol_halo = 0

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx

        IF ( p_patch(jg)%cells%decomp_info%decomp_domain(jc,jb) == 0 ) THEN  ! interior cells
          jcol_int = jcol_int + 1
          p_mgmgrid(jg)%c_nb(jc,jb)% jcol_c(1) = jcol_int
        ELSE IF ( p_patch(jg)%cells%decomp_info%decomp_domain(jc,jb) > 0 ) THEN  ! halos
          jcol_halo = jcol_halo + 1
          p_mgmgrid(jg)%c_nb(jc,jb)% jcol_c(1) = -jcol_halo
        ELSE   ! This must not happen
          p_mgmgrid(jg)%c_nb(jc,jb)% jcol_c(1) = 0
        END IF
        p_mgmgrid(jg)%c_nb(jc,jb)% jcb_neighbor(:,1) = (/jc,jb/)

        p_mgmgrid(jg)%c_nb(jc,jb)% x_cart(1)  &
          &  = REAL(p_patch(jg)%cells%cartesian_center(jc,jb)% x(1), KIND=triprc)
        p_mgmgrid(jg)%c_nb(jc,jb)% y_cart(1)  &
          &  = REAL(p_patch(jg)%cells%cartesian_center(jc,jb)% x(2), KIND=triprc)
        DO jv = 1, 3
          vidx(jv) = p_patch(jg)%cells%vertex_idx(jc,jb,jv)
          vblk(jv) = p_patch(jg)%cells%vertex_blk(jc,jb,jv)
          p_mgmgrid(jg)%c_nb(jc,jb)% v_cart(jv,1)% x  &
            &  = p_patch(jg)%verts%cartesian(vidx(jv),vblk(jv))%x(1)
          p_mgmgrid(jg)%c_nb(jc,jb)% v_cart(jv,1)% y  &
            &  = p_patch(jg)%verts%cartesian(vidx(jv),vblk(jv))%x(2)
          p_mgmgrid(jg)%c_nb(jc,jb)% v_cart(jv,1)% z  &
            &  = p_patch(jg)%verts%cartesian(vidx(jv),vblk(jv))%x(3)
        ENDDO

      ENDDO  ! jc
    ENDDO  ! jb

    p_mgmgrid(jg)% ncol_int  = jcol_int
    p_mgmgrid(jg)% ncol_halo = jcol_halo

    SELECT CASE ( p_patch(jg)%geometry_info%geometry_type )

    CASE ( sphere_geometry )
      !
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          p_mgmgrid(jg)%coslat(jc,jb)  &
            &  = MAX(min_coslat, ABS(COS(p_patch(jg)%cells%center(jc,jb)%lat)))
!cl       p_mgmgrid(jg)%c_nb(jc,jb)% coslat(1)  &
!cl         &  = REAL(p_mgmgrid(jg)%coslat(jc,jb), KIND=triprc)
          p_mgmgrid(jg)%c_nb(jc,jb)% sinlat(1)  &
            &  = REAL(SIGN(SQRT(1. - p_mgmgrid(jg)%coslat(jc,jb)**2),  &
            &              p_patch(jg)%cells%center(jc,jb)%lat), KIND=triprc)
        ENDDO
      ENDDO

    CASE ( planar_torus_geometry )
      !
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          p_mgmgrid(jg)%coslat(jc,jb) = 1.
!cl       p_mgmgrid(jg)%c_nb(jc,jb)% coslat(1) = 1.
          p_mgmgrid(jg)%c_nb(jc,jb)% sinlat(1) = 0.
        ENDDO
      ENDDO
      darcraymin = 0.     ! reset to de-activate
      darcraymax = 1.e10

    CASE DEFAULT
      !
      CALL finish('setup_msgwam', 'Geometry should be spherical or torus')

    END SELECT

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        p_mgmgrid(jg)%c_nb(jc,jb)% jcol_c(2:)         = -(p_mgmgrid(jg)%ncol_halo+1)
        p_mgmgrid(jg)%c_nb(jc,jb)% jcb_neighbor(:,2:) = -9999
!cl     p_mgmgrid(jg)%c_nb(jc,jb)% coslat(2:) = -9999.0_triprc
        p_mgmgrid(jg)%c_nb(jc,jb)% sinlat(2:) = -9999.0_triprc
        p_mgmgrid(jg)%c_nb(jc,jb)% x_cart(2:) = -9999.0_triprc
        p_mgmgrid(jg)%c_nb(jc,jb)% y_cart(2:) = -9999.0_triprc
        p_mgmgrid(jg)%c_nb(jc,jb)% v_cart(:,2:)% x = -9999.0_wp
        p_mgmgrid(jg)%c_nb(jc,jb)% v_cart(:,2:)% y = -9999.0_wp
        p_mgmgrid(jg)%c_nb(jc,jb)% v_cart(:,2:)% z = -9999.0_wp
      ENDDO
    ENDDO

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx

        ! Find neighbor cells  -->  jc21, jb21 (1:n_neighbor)
        jc21(1:3) = p_patch(jg)%cells%neighbor_idx(jc,jb,:)
        jb21(1:3) = p_patch(jg)%cells%neighbor_blk(jc,jb,:)
        DO jv = 1, 3
          jn = (jv-1)*6 + 4
          vidx(jv) = p_patch(jg)%cells%vertex_idx(jc,jb,jv)
          vblk(jv) = p_patch(jg)%cells%vertex_blk(jc,jb,jv)
          jc21(jn:jn+5) = p_patch(jg)%verts%cell_idx(vidx(jv),vblk(jv),:)
          jb21(jn:jn+5) = p_patch(jg)%verts%cell_blk(vidx(jv),vblk(jv),:)
        ENDDO
!       DO jn = 4, 21
!         IF (     ( jc21(jn) == jc      .AND. jb21(jn) == jb      )    &
!           & .OR. ( jc21(jn) == jc21(1) .AND. jb21(jn) == jb21(1) )    &
!           & .OR. ( jc21(jn) == jc21(2) .AND. jb21(jn) == jb21(2) )    &
!           & .OR. ( jc21(jn) == jc21(3) .AND. jb21(jn) == jb21(3) ) )  &
!           &  jc21(jn) = 0
!       ENDDO
        ! The above loop does not work where the number of neighbors is less than 12
        ! (i.e., in 60 interior cells around the 12 vertices of the original icosahedron:
        ! 5 x 12, or in halo lines). It is because around such a cell, the vertex-based
        ! cell searching (using p_patch%verts%cell_XXX) currently gives 6 cells somehow,
        ! with duplicated index pair(s). Thus, all the duplications should be removed to
        ! get a correct number of neighbors, as the following.
        DO jn = 2, 21
          IF ( ( jc21(jn) == jc .AND. jb21(jn) == jb )                                &
            &  .OR. ANY( jc21(jn) == jc21(1:jn-1) .AND. jb21(jn) == jb21(1:jn-1) ) )  &
            &  jc21(jn) = 0
        ENDDO

        ! For halo cells, further check is required, as it has been observed that
        ! the vertex-based cell searching (using p_patch%verts%cell_XXX) gives some
        ! index pairs which are not really for the neighbors !

        IF (p_patch(jg)%cells%decomp_info%decomp_domain(jc,jb) > 0) THEN
          ! Halo cells

          DO jn = 1, 21
            IF (jc21(jn) == 0)  CYCLE
            icheck = 0
            DO jv = 1, 3
              vidx(jv) = p_patch(jg)%cells%vertex_idx(jc21(jn),jb21(jn),jv)
              vblk(jv) = p_patch(jg)%cells%vertex_blk(jc21(jn),jb21(jn),jv)
              IF (ANY(      vidx(jv) == p_patch(jg)%cells%vertex_idx(jc,jb,:)     &
                &     .AND. vblk(jv) == p_patch(jg)%cells%vertex_blk(jc,jb,:) ))  &
                &  icheck = icheck + 1
            ENDDO
            IF ( icheck /= 1 .AND. icheck /= 2 )  jc21(jn) = 0
          ENDDO

        ELSE IF (p_patch(jg)%cells%decomp_info%decomp_domain(jc,jb) < 0) THEN
          CALL finish('setup_msgwam', 'wrong index range')
!       ELSE   ! p_patch(jg)%cells%decomp_info%decomp_domain(jc,jb) == 0
          ! Interior cells: Nothing to do additionally
        END IF

        ! Count the number of neighbors and shift the array elements removing missings
        n_neighbor = COUNT(jc21 > 0)
        jb21(:n_neighbor) = PACK(jb21, jc21 > 0)  ! do for jb21 first, before jc21 is updated
        jc21(:n_neighbor) = PACK(jc21, jc21 > 0)

        IF ( p_patch(jg)%cells%decomp_info%decomp_domain(jc,jb) == 0  &  ! interior cells
          &  .AND. n_neighbor /= 11 .AND. n_neighbor /= 12 )          &
          &  CALL finish('setup_msgwam', 'Error in searching neighbor cells (i)')
        IF ( p_patch(jg)%cells%decomp_info%decomp_domain(jc,jb) > 0   &  ! halos
          &  .AND. ( n_neighbor > 12 .OR. n_neighbor < 3 ) )          &
          &  CALL finish('setup_msgwam', 'Error in searching neighbor cells of halos')

        ! Save info of neighbor cells
        p_mgmgrid(jg)%c_nb(jc,jb)% n_neighbor = n_neighbor
        DO jn = 1, n_neighbor
          jc_n = jc21(jn)
          jb_n = jb21(jn)
          p_mgmgrid(jg)%       c_nb(jc  ,jb  )% jcol_c(1+jn)  &
            &  = p_mgmgrid(jg)%c_nb(jc_n,jb_n)% jcol_c(1)
          p_mgmgrid(jg)%       c_nb(jc  ,jb  )% jcb_neighbor(:,1+jn)  &
            &  = p_mgmgrid(jg)%c_nb(jc_n,jb_n)% jcb_neighbor(:,1)
!cl       p_mgmgrid(jg)%       c_nb(jc  ,jb  )% coslat(1+jn)  &
!cl         &  = p_mgmgrid(jg)%c_nb(jc_n,jb_n)% coslat(1)
          p_mgmgrid(jg)%       c_nb(jc  ,jb  )% sinlat(1+jn)  &
            &  = p_mgmgrid(jg)%c_nb(jc_n,jb_n)% sinlat(1)
          p_mgmgrid(jg)%       c_nb(jc  ,jb  )% x_cart(1+jn)  &
            &  = p_mgmgrid(jg)%c_nb(jc_n,jb_n)% x_cart(1)
          p_mgmgrid(jg)%       c_nb(jc  ,jb  )% y_cart(1+jn)  &
            &  = p_mgmgrid(jg)%c_nb(jc_n,jb_n)% y_cart(1)
          p_mgmgrid(jg)%       c_nb(jc  ,jb  )% v_cart(:,1+jn)  &
            &  = p_mgmgrid(jg)%c_nb(jc_n,jb_n)% v_cart(:,1)
        ENDDO

      ENDDO ! jc
    ENDDO ! jb

    ! TODO: for test
    IF (imethod_split > 0) THEN
      DO jk = 2, nlevp1-1
        jks = jk + p_patch(jg)%nshift_total
        p_mgmgrid(jg)%dz_crit(jk) = MAX( dz_crit_min,  &
          &                           dz_factor_split*0.5_wp*(vct_a(jks-1) - vct_a(jks+1)) )
      ENDDO
      p_mgmgrid(jg)%dz_crit(nlevp1) = p_mgmgrid(jg)%dz_crit(nlevp1-1)
      p_mgmgrid(jg)%dz_crit(1)      = p_mgmgrid(jg)%dz_crit(2)
    ELSE
      p_mgmgrid(jg)%dz_crit(:) = 1.e20_wp
    END IF

    IF (imethod_split > 0) THEN
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          p_mgmgrid(jg)%darc_crit(jc,jb) = dh_factor_split  &
            &      *SQRT(p_patch(jg)%cells%area(jc,jb))
          p_mgmgrid(jg)%area_crit(jc,jb) = da_factor_split  &
            &      *p_patch(jg)%cells%area(jc,jb)
        ENDDO
        IF ( p_patch(jg)%geometry_info%geometry_type /= planar_torus_geometry ) THEN
          DO jc = i_startidx, i_endidx
            p_mgmgrid(jg)%darc_crit(jc,jb) = p_mgmgrid(jg)%darc_crit(jc,jb)  &
              &                             /grid_sphere_radius
            p_mgmgrid(jg)%area_crit(jc,jb) = p_mgmgrid(jg)%area_crit(jc,jb)  &
              &                             /grid_sphere_radius**2
          ENDDO
        END IF
      ENDDO
    ELSE
      p_mgmgrid(jg)%darc_crit(:,:) = 1.e20_wp
      p_mgmgrid(jg)%area_crit(:,:) = 1.e20_wp
    END IF
!YHK-

    ! Prognostic domain
    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int
    i_startblk = p_patch(jg)%cells%start_block(rl_start)
    i_endblk   = p_patch(jg)%cells%end_block(rl_end)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        p_mgmgrid(jg)%darc_intpol(jc,jb) = dh_factor_intpol  &
          &                               *SQRT(p_patch(jg)%cells%area(jc,jb))
      ENDDO
      IF ( p_patch(jg)%geometry_info%geometry_type /= planar_torus_geometry ) THEN
        DO jc = i_startidx, i_endidx
          p_mgmgrid(jg)%darc_intpol(jc,jb) = p_mgmgrid(jg)%darc_intpol(jc,jb)  &
            &                               /grid_sphere_radius
        ENDDO
      END IF
    ENDDO

    ! Initialization
    p_gridinfo4ray(jg)%cellneighbors_idx(:,:,:)    = 0
    p_gridinfo4ray(jg)%cellneighbors_blk(:,:,:)    = 0
    p_gridinfo4ray(jg)%cellneighbors_nstencil(:,:) = 0
    p_gridinfo4ray(jg)%cellvertices_lon(:,:,:)     = -9999._wp
    p_gridinfo4ray(jg)%cellvertices_lat(:,:,:)     = -9999._wp
    p_gridinfo4ray(jg)%cellvertices_x(:,:,:)       = -9999._wp
    p_gridinfo4ray(jg)%cellvertices_y(:,:,:)       = -9999._wp
    p_gridinfo4ray(jg)%cellvertices_z(:,:,:)       = -9999._wp

    IF (ltest_hprop) THEN
      ! Name output file for horizontal propagation test
      filename_ltest_hprop = 'output_ltest_hprop.txt'
      ! Open output file for horizontal propagation test
      OPEN(987,FILE=filename_ltest_hprop,STATUS='UNKNOWN',POSITION='APPEND')
    ENDIF

    ! Are we on the Sphere or on the Torus?
    SELECT CASE(p_patch(jg)%geometry_info%geometry_type)

    CASE(sphere_geometry)
    ! Sphere:
    ! The function inside_triangle will need "global" cartesian
    ! coordinates with the Earth's center as origin. Here we
    ! trasnsform spherical coordinates to "global" cartesian
    ! coordinates.

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb, jc, jv, jec, jtri, i_startidx, i_endidx, istencil, vidx, vblk, cidx_v, cblk_v, cidx, cblk) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          istencil = 0

          ! Get line and block indices of cell vertices
          vidx(1:3) = p_patch(jg)%cells%vertex_idx(jc,jb,1:3)
          vblk(1:3) = p_patch(jg)%cells%vertex_blk(jc,jb,1:3)

          ! Store longitudes and latitudes of cell vertices
          DO jv = 1,3
            p_gridinfo4ray(jg)%cellvertices_lon(jv,jc,jb) = p_patch(jg)%verts%vertex(vidx(jv),vblk(jv))%lon
            p_gridinfo4ray(jg)%cellvertices_lat(jv,jc,jb) = p_patch(jg)%verts%vertex(vidx(jv),vblk(jv))%lat

            ! Transform spherical to Cartesian coordinates and store them
            p_gridinfo4ray(jg)%cellvertices_x(jv,jc,jb) = COS(p_gridinfo4ray(jg)%cellvertices_lat(jv,jc,jb))&
                                                         *COS(p_gridinfo4ray(jg)%cellvertices_lon(jv,jc,jb))
            p_gridinfo4ray(jg)%cellvertices_y(jv,jc,jb) = COS(p_gridinfo4ray(jg)%cellvertices_lat(jv,jc,jb))&
                                                         *SIN(p_gridinfo4ray(jg)%cellvertices_lon(jv,jc,jb))
            p_gridinfo4ray(jg)%cellvertices_z(jv,jc,jb) = SIN(p_gridinfo4ray(jg)%cellvertices_lat(jv,jc,jb))
          ENDDO

          ! For each vertex: get all the cells which share this vertex
          DO jv = 1,3
            nbcells_v    = p_patch(jg)%verts%num_edges(vidx(jv),vblk(jv))
            cidx_v(:,jv) = -9999
            cblk_v(:,jv) = -9999
            cidx_v(1:nbcells_v,jv) = p_patch(jg)%verts%cell_idx(vidx(jv),vblk(jv),1:nbcells_v)
            cblk_v(1:nbcells_v,jv) = p_patch(jg)%verts%cell_blk(vidx(jv),vblk(jv),1:nbcells_v)
          ENDDO

          ! 1st add the 3 direct neighbors to the stencil
          DO jec = 1, 3

            istencil = istencil + 1

            ! Get line and block indices of direct neighbors
            cidx(jec) = p_patch(jg)%cells%neighbor_idx(jc,jb,jec)
            cblk(jec) = p_patch(jg)%cells%neighbor_blk(jc,jb,jec)

            ! Store indices of direct neighbor cells
            p_gridinfo4ray(jg)%cellneighbors_idx(istencil,jc,jb) = cidx(jec)
            p_gridinfo4ray(jg)%cellneighbors_blk(istencil,jc,jb) = cblk(jec)
          ENDDO

          ! 2nd loop over the vertices and add all the cells
          ! that are no direct neighbors and not the actual cell
          ! we are sitting in
          DO jv = 1,3     ! loop over vertices
            DO jtri = 1,6   ! loop over cells around each vertex

              IF (.NOT.( (cidx_v(jtri,jv) == cidx(1) .AND. cblk_v(jtri,jv) == cblk(1))  &
                &  .OR.  (cidx_v(jtri,jv) == cidx(2) .AND. cblk_v(jtri,jv) == cblk(2))  &
                &  .OR.  (cidx_v(jtri,jv) == cidx(3) .AND. cblk_v(jtri,jv) == cblk(3))  &
                &  .OR.  (cidx_v(jtri,jv) == jc      .AND. cblk_v(jtri,jv) == jb) &
                &  .OR.  (cidx_v(jtri,jv) == 0       .AND. cblk_v(jtri,jv) == 0 ) &
                &  .OR.  (cidx_v(jtri,jv) == -9999   .AND. cblk_v(jtri,jv) == -9999 )) ) THEN

                istencil = istencil + 1

                ! Store indices of non-direct neighbor cells
                p_gridinfo4ray(jg)%cellneighbors_idx(istencil,jc,jb) = cidx_v(jtri,jv)
                p_gridinfo4ray(jg)%cellneighbors_blk(istencil,jc,jb) = cblk_v(jtri,jv)
              ENDIF

            ENDDO  !jtri
          ENDDO  !jv

          ! Store number of neighboring cells (either 12 or 11)
          p_gridinfo4ray(jg)%cellneighbors_nstencil(jc,jb) = istencil

        ENDDO  !jc
      ENDDO  !jb
!!$OMP END DO
!!$OMP END PARALLEL

      !IF (ltest_hprop .AND. l1ray) THEN
      IF (ltest_hprop) THEN

        ! Store lat, lon coordinates of halo
        ! cells for diagnostic purpose only

        ! Halo cells
        rl_start = min_rlcell_int - 1
        rl_end   = min_rlcell_int - 2
        i_startblk = p_patch(jg)%cells%start_block(rl_start)
        i_endblk   = p_patch(jg)%cells%end_block(rl_end)

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb, jc, jv, jec, jtri, i_startidx, i_endidx, istencil, vidx, vblk, cidx_v, cblk_v, cidx, cblk) ICON_OMP_GUIDED_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            ! Get line and block indices of cell vertices
            vidx(1:3) = p_patch(jg)%cells%vertex_idx(jc,jb,1:3)
            vblk(1:3) = p_patch(jg)%cells%vertex_blk(jc,jb,1:3)

            ! Store longitudes and latitudes of cell vertices
            DO jv = 1,3
              p_gridinfo4ray(jg)%cellvertices_lon(jv,jc,jb) = p_patch(jg)%verts%vertex(vidx(jv),vblk(jv))%lon
              p_gridinfo4ray(jg)%cellvertices_lat(jv,jc,jb) = p_patch(jg)%verts%vertex(vidx(jv),vblk(jv))%lat
            ENDDO ! jv

          ENDDO  !jc
        ENDDO  !jb
!!$OMP END DO
!!$OMP END PARALLEL

        ! Diagnose if we managed to get a comprehensive stencil:
        ! print out lat, lon coordinates of all cells and those
        ! of their neighbors

        ! Wait until all PEs reach this point
        CALL p_barrier(p_comm_work) ! Not sure it is really necessary
                                    ! but this is a rarely used diagnostic
                                    ! anyway, so no worries if this slows
                                    ! anything down

        ! Prognostic domain
        rl_start = grf_bdywidth_c + 1
        rl_end   = min_rlcell_int
        i_startblk = p_patch(jg)%cells%start_block(rl_start)
        i_endblk   = p_patch(jg)%cells%end_block(rl_end)

        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            WRITE(987,'(a,3i4,2x,6F14.5)') &
                   'gridinfo4ray CELL: nstencil, jc, jb, lon, lat:', &
                                      p_gridinfo4ray(jg)%cellneighbors_nstencil(jc,jb),&
                                                                                jc,jb ,&
                                      p_gridinfo4ray(jg)%cellvertices_lon(1,jc,jb)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lat(1,jc,jb)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lon(2,jc,jb)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lat(2,jc,jb)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lon(3,jc,jb)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lat(3,jc,jb)*rad2deg

            DO jstencil = 1, p_gridinfo4ray(jg)%cellneighbors_nstencil(jc,jb)

              icn = p_gridinfo4ray(jg)%cellneighbors_idx(jstencil,jc,jb) ! jc index of neighboring cell
              ibn = p_gridinfo4ray(jg)%cellneighbors_blk(jstencil,jc,jb) ! jb index of neighboring cell

              WRITE(987,'(a,3i4,2x,6F14.5)') &
                   'gridinfo4ray NEIGHBORS: jstencil, icn, ibn, lon, lat:', &
                                                                           jstencil, &
                                                                            icn,ibn, &
                                      p_gridinfo4ray(jg)%cellvertices_lon(1,icn,ibn)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lat(1,icn,ibn)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lon(2,icn,ibn)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lat(2,icn,ibn)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lon(3,icn,ibn)*rad2deg,&
                                      p_gridinfo4ray(jg)%cellvertices_lat(3,icn,ibn)*rad2deg

            ENDDO ! jstencil

          ENDDO ! jc

        ENDDO ! jb

      ENDIF ! (ltest_hprop .AND. l1ray)

    CASE(planar_torus_geometry)
    ! Torus:
    ! The function inside_triangle_torus will need
    ! horizontal cartesian coordinates of the tangential
    ! plane

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb, jc, jv, jec, jtri, i_startidx, i_endidx, istencil, vidx, vblk, cidx_v, cblk_v, cidx, cblk) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx

          istencil = 0

          ! Get line and block indices of cell vertices
          vidx(1:3) = p_patch(jg)%cells%vertex_idx(jc,jb,1:3)
          vblk(1:3) = p_patch(jg)%cells%vertex_blk(jc,jb,1:3)

          ! Store longitudes and latitudes of cell vertices
          DO jv = 1,3
            p_gridinfo4ray(jg)%cellvertices_lon(jv,jc,jb) = p_patch(jg)%verts%vertex(vidx(jv),vblk(jv))%lon
            p_gridinfo4ray(jg)%cellvertices_lat(jv,jc,jb) = p_patch(jg)%verts%vertex(vidx(jv),vblk(jv))%lat

            ! Transform spherical to Cartesian coordinates and store them
            p_gridinfo4ray(jg)%cellvertices_x(jv,jc,jb) = p_patch(jg)%verts%cartesian(vidx(jv),vblk(jv))%x(1)
            p_gridinfo4ray(jg)%cellvertices_y(jv,jc,jb) = p_patch(jg)%verts%cartesian(vidx(jv),vblk(jv))%x(2)
            p_gridinfo4ray(jg)%cellvertices_z(jv,jc,jb) = p_patch(jg)%verts%cartesian(vidx(jv),vblk(jv))%x(3)
          ENDDO

          ! For each vertex: get all the cells which share this vertex
          DO jv = 1,3
            nbcells_v    = p_patch(jg)%verts%num_edges(vidx(jv),vblk(jv))
            cidx_v(:,jv) = -9999
            cblk_v(:,jv) = -9999
            cidx_v(1:nbcells_v,jv) = p_patch(jg)%verts%cell_idx(vidx(jv),vblk(jv),1:nbcells_v)
            cblk_v(1:nbcells_v,jv) = p_patch(jg)%verts%cell_blk(vidx(jv),vblk(jv),1:nbcells_v)
          ENDDO

          ! 1st add the 3 direct neighbors to the stencil
          DO jec = 1, 3

            istencil = istencil + 1

            ! Get line and block indices of direct neighbors
            cidx(jec) = p_patch(jg)%cells%neighbor_idx(jc,jb,jec)
            cblk(jec) = p_patch(jg)%cells%neighbor_blk(jc,jb,jec)

            ! Store indices of direct neighbor cells
            p_gridinfo4ray(jg)%cellneighbors_idx(istencil,jc,jb) = cidx(jec)
            p_gridinfo4ray(jg)%cellneighbors_blk(istencil,jc,jb) = cblk(jec)
          ENDDO

          ! 2nd loop over the vertices and add all the cells
          ! that are no direct neighbors and not the actual cell
          ! we are sitting in
          DO jv = 1,3     ! loop over vertices
            DO jtri = 1,6   ! loop over cells around each vertex

              IF (.NOT.( (cidx_v(jtri,jv) == cidx(1) .AND. cblk_v(jtri,jv) == cblk(1))  &
                &  .OR.  (cidx_v(jtri,jv) == cidx(2) .AND. cblk_v(jtri,jv) == cblk(2))  &
                &  .OR.  (cidx_v(jtri,jv) == cidx(3) .AND. cblk_v(jtri,jv) == cblk(3))  &
                &  .OR.  (cidx_v(jtri,jv) == jc      .AND. cblk_v(jtri,jv) == jb) &
                &  .OR.  (cidx_v(jtri,jv) == 0       .AND. cblk_v(jtri,jv) == 0 ) &
                &  .OR.  (cidx_v(jtri,jv) == -9999   .AND. cblk_v(jtri,jv) == -9999 )) ) THEN

                istencil = istencil + 1

                ! Store indices of non-direct neighbor cells
                p_gridinfo4ray(jg)%cellneighbors_idx(istencil,jc,jb) = cidx_v(jtri,jv)
                p_gridinfo4ray(jg)%cellneighbors_blk(istencil,jc,jb) = cblk_v(jtri,jv)
              ENDIF

            ENDDO  !jtri
          ENDDO  !jv

          ! Store number of neighboring cells (either 12 or 11)
          p_gridinfo4ray(jg)%cellneighbors_nstencil(jc,jb) = istencil

        ENDDO  !jc
      ENDDO  !jb
!!$OMP END DO
!!$OMP END PARALLEL

      !IF (ltest_hprop .AND. l1ray) THEN
      IF (ltest_hprop) THEN

        ! Store x,y coordinates of halo
        ! cells for diagnostic purpose only

        ! Halo cells
        rl_start = min_rlcell_int - 1
        rl_end   = min_rlcell_int - 2
        i_startblk = p_patch(jg)%cells%start_block(rl_start)
        i_endblk   = p_patch(jg)%cells%end_block(rl_end)

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb, jc, jv, jec, jtri, i_startidx, i_endidx, istencil, vidx, vblk, cidx_v, cblk_v, cidx, cblk) ICON_OMP_GUIDED_SCHEDULE
        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            ! Get line and block indices of cell vertices
            vidx(1:3) = p_patch(jg)%cells%vertex_idx(jc,jb,1:3)
            vblk(1:3) = p_patch(jg)%cells%vertex_blk(jc,jb,1:3)

            ! Store x, y coords of cell vertices
            DO jv = 1,3
              p_gridinfo4ray(jg)%cellvertices_x(jv,jc,jb) = p_patch(jg)%verts%cartesian(vidx(jv),vblk(jv))%x(1)
              p_gridinfo4ray(jg)%cellvertices_y(jv,jc,jb) = p_patch(jg)%verts%cartesian(vidx(jv),vblk(jv))%x(2)
            ENDDO ! jv

          ENDDO  !jc
        ENDDO  !jb
!!$OMP END DO
!!$OMP END PARALLEL

        ! Diagnose if we managed to get a comprehensive stencil:
        ! print out x, y coordinates of all cells and those
        ! of their neighbors

        ! Wait until all PEs reach this point
        CALL p_barrier(p_comm_work) ! Not sure it is really necessary
                                    ! but this is a rarely used diagnostic
                                    ! anyway, so no worries if this slows
                                    ! anything down

        ! Prognostic domain
        rl_start = grf_bdywidth_c + 1
        rl_end   = min_rlcell_int
        i_startblk = p_patch(jg)%cells%start_block(rl_start)
        i_endblk   = p_patch(jg)%cells%end_block(rl_end)

        DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

          DO jc = i_startidx, i_endidx

            WRITE(987,'(a,3i4,2x,6F14.5)') &
                   'gridinfo4ray CELL: nstencil, jc, jb, x, y:', &
                                      p_gridinfo4ray(jg)%cellneighbors_nstencil(jc,jb),&
                                                                                jc,jb ,&
                                      p_gridinfo4ray(jg)%cellvertices_x(1,jc,jb),&
                                      p_gridinfo4ray(jg)%cellvertices_y(1,jc,jb),&
                                      p_gridinfo4ray(jg)%cellvertices_x(2,jc,jb),&
                                      p_gridinfo4ray(jg)%cellvertices_y(2,jc,jb),&
                                      p_gridinfo4ray(jg)%cellvertices_x(3,jc,jb),&
                                      p_gridinfo4ray(jg)%cellvertices_y(3,jc,jb)

            DO jstencil = 1, p_gridinfo4ray(jg)%cellneighbors_nstencil(jc,jb)

              icn = p_gridinfo4ray(jg)%cellneighbors_idx(jstencil,jc,jb) ! jc index of neighboring cell
              ibn = p_gridinfo4ray(jg)%cellneighbors_blk(jstencil,jc,jb) ! jb index of neighboring cell

              WRITE(987,'(a,3i4,2x,6F14.5)') &
                   'gridinfo4ray NEIGHBORS: jstencil, icn, ibn, x, y:', &
                                                                         jstencil, &
                                                                          icn,ibn, &
                                      p_gridinfo4ray(jg)%cellvertices_x(1,icn,ibn),&
                                      p_gridinfo4ray(jg)%cellvertices_y(1,icn,ibn),&
                                      p_gridinfo4ray(jg)%cellvertices_x(2,icn,ibn),&
                                      p_gridinfo4ray(jg)%cellvertices_y(2,icn,ibn),&
                                      p_gridinfo4ray(jg)%cellvertices_x(3,icn,ibn),&
                                      p_gridinfo4ray(jg)%cellvertices_y(3,icn,ibn)

            ENDDO ! jstencil

          ENDDO ! jc

        ENDDO ! jb

      ENDIF ! (ltest_hprop .AND. l1ray)

      ! Store edge length and height of triangle
      rl_start = grf_bdywidth_e + 1
      rl_end   = min_rledge_int
      i_startblk = p_patch(jg)%edges%start_block(rl_start)
      i_endblk   = p_patch(jg)%edges%end_block  (rl_end)

      CALL get_indices_e(p_patch(jg), i_startblk, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)

      IF (p_patch(jg)%edges%primal_edge_length(i_startidx,i_startblk) /= &
          p_patch(jg)%edges%primal_edge_length(  i_endidx,i_startblk))   &
        CALL message('setup_msgwam', 'Warning: triangle edges vary in in horizontal!')
      edge_length_triangle = p_patch(jg)%edges%primal_edge_length(i_startidx,i_startblk)
      height_triangle      = edge_length_triangle * 0.5_wp * SQRT(3._wp)

      ! Store x, y coordinates of domain edges
      domxmin = p_patch(jg)%geometry_info%center%x(1) &
              - p_patch(jg)%geometry_info%domain_length*0.5_wp &
              - edge_length_triangle*0.5_wp
      domxmax = p_patch(jg)%geometry_info%center%x(1) &
              + p_patch(jg)%geometry_info%domain_length*0.5_wp &
              - edge_length_triangle
      domymin = p_patch(jg)%geometry_info%center%x(2) &
              - p_patch(jg)%geometry_info%domain_height*0.5_wp
      domymax = p_patch(jg)%geometry_info%center%x(2) &
              + p_patch(jg)%geometry_info%domain_height*0.5_wp &
              - height_triangle
      ! The domain is not completely symmetric to the "center"
      ! the domain stored in geometry_info%center%x. The "rule"
      ! above to get the min, max values of vertex x,y coordinates
      ! is somewhat empirical.

    CASE DEFAULT

      CALL finish('setup_msgwam', &
           'Geometry not recognized! Should be spherical or torus!')

    END SELECT

    IF (ltest_hprop) CLOSE(987)

  ENDDO ! jg

  CALL message ('setup_msgwam:', 'MS-GWaM: allocation of WKB ray type completed')

END SUBROUTINE setup_msgwam
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE cleanup_msgwam(n_dom)
  INTEGER,      INTENT(IN)   :: n_dom      ! number of model domains
  INTEGER                    :: jg, ist

  ! Deallocate WKB ray type
  DO jg = 1, n_dom

    CALL vlr_del( p_msgwam_list(jg) )

    IF (.NOT. lmsgwam(jg))  CYCLE

    IF ( ALLOCATED(p_spl(jg)%jcol) )  DEALLOCATE( p_spl(jg)%jcol )

    CALL vlr_del( p_ray_list(jg) )

    CALL vlr_del( p_spl_list(jg) )

    CALL vlr_del( p_rwork_list(jg) )

    IF ( gws_conv_config%n_source(jg) > 0 )  &
      &  CALL vlr_del( p_ray_conv_list(jg) )

    CALL vlr_del( p_ray_bg_list(jg) )

    CALL vlr_del( p_lfluxbg_list(jg) )

  ENDDO ! jg

  DEALLOCATE(p_msgwam, p_msgwam_list, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of MS-GWaM field variables')
  
  IF (ALL( .NOT. lmsgwam(:) ))  RETURN

  DEALLOCATE(ray_coll, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of ray_coll failed')

  DEALLOCATE(p_ray, p_ray_list, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of WKB ray type failed')

  DEALLOCATE(p_spl, p_spl_list, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of WKB spl type failed')

  DEALLOCATE(p_rwork, p_rwork_list, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of WKB ray work type failed')

  IF ( ALLOCATED(p_ray_conv) ) THEN
    DEALLOCATE(p_ray_conv, p_ray_conv_list, STAT=ist)
    IF (ist/=success) CALL finish ('cleanup_msgwam',&
         & 'MS-GWaM: deallocation of WKB ray type failed (conv)')
  END IF

  DEALLOCATE(p_ray_bg, p_ray_bg_list, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of WKB ray type failed (bg)')

  DEALLOCATE(p_mgmgrid, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of mgmgrid failed')

  DEALLOCATE(p_lfluxbg, p_lfluxbg_list, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of WKB ray grid failed')

  DEALLOCATE(omega_l, domega_l, sp_omega_l, STAT=ist)
  IF (ist/=success) CALL finish ('cleanup_msgwam',&
       & 'MS-GWaM: deallocation of omega_l failed')

  CALL message ('cleanup_msgwam', 'MS-GWaM: deallocation of WKB ray type completed')

END SUBROUTINE cleanup_msgwam
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE new_ray_list( k_jg, nray_nspec, nspec, kblks, listname, p_list, p_var, name_head )

  TYPE(t_var_list_ptr), INTENT(INOUT) :: p_list
  TYPE(t_ray),          INTENT(INOUT) :: p_var
  INTEGER,              INTENT(IN)    :: k_jg, kblks
  INTEGER,              INTENT(IN)    :: nray_nspec  ! nrays or nspec
  INTEGER,              INTENT(IN)    :: nspec       ! for jr_last when name_head == 'ray'
  CHARACTER(len=*),     INTENT(IN)    :: listname
  CHARACTER(len=*),     INTENT(IN)    :: name_head

  ! Local variables

  TYPE(t_cf_var)    ::    cf_desc
  TYPE(t_grib2_var) :: grib2_desc

  INTEGER           :: shape3d_ray(3)
  INTEGER           :: shape3d_spec(3)
  INTEGER           :: shape2d_grid(2)
  INTEGER           :: ibits
  INTEGER           :: datatype_flt

  IF ( lnetcdf_flt64_output ) THEN
    datatype_flt = DATATYPE_FLT64
  ELSE
    datatype_flt = DATATYPE_FLT32
  ENDIF

  ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

  shape3d_ray = (/nproma, nray_nspec, kblks/)
  shape3d_spec = (/nproma, nspec, kblks/)
  shape2d_grid = (/nproma, kblks/)

  CALL vlr_add( p_list, TRIM(listname), patch_id=k_jg, lrestart=.FALSE. )

  !--------------------------------
  ! Add ray properties as variables
  !--------------------------------

  ! &      lon(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%lon', 'deg', 'lon of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%lon', p_var%lon, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      lat(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%lat', 'deg', 'lat of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%lat', p_var%lat, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      z(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%z', 'm', 'height of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%z', p_var%z, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      dlon(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%dlon', 'deg', 'dlon of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%dlon', p_var%dlon, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      dlat(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%dlat', 'deg', 'dlat of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%dlat', p_var%dlat, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      dz(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%dz', 'm', 'vertical size of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%dz', p_var%dz, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      coslat(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%coslat', '', 'cosine lat of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%coslat', p_var%coslat, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      k(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%k', 'm^-1', 'latitude wavenumber of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%k', p_var%k, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      l(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%l', 'm^-1', 'meridional wavenumber of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%l', p_var%l, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      m(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%m', 'm^-1', 'vertical wavenumber of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%m', p_var%m, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      dk(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%dk', 'm^-1', 'dk spectral size of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%dk', p_var%dk, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      dl(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%dl', 'm^-1', 'dl spectral size of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%dl', p_var%dl, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      dm(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%dm', 'm^-1', 'dm spectral size of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%dm', p_var%dm, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  ! &      wadens(nproma,nray_nspec,nblks)
  cf_desc    = t_cf_var(name_head//'%wadens', 'm^3s^-1', 'wave action density of GW ray', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, name_head//'%wadens', p_var%wadens, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
  !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
              & ldims=shape3d_ray, &
              & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  IF ( name_head == 'ray' .OR. name_head == 'rwork' .OR. name_head == 'spec' ) THEN
    ! &      iexist(nproma,nray_nspec,nblks)
    cf_desc    = t_cf_var(name_head//'%iexist', '0/1', 'existence of GW ray', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%iexist', p_var%iexist, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
  ENDIF

  IF ( name_head == 'ray' .OR. name_head == 'rwork' ) THEN

    ! &      specid(nproma,nray_nspec,nblks)
    cf_desc    = t_cf_var(name_head//'%specid', 'N/A', 'spectral id of GW ray', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%specid', p_var%specid, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! &      jk_active(nproma,nray_nspec,nblks)
    cf_desc    = t_cf_var(name_head//'%jk_active', 'N/A', 'jk_active', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%jk_active', p_var%jk_active, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  END IF

  IF ( name_head == 'ray' ) THEN

    ! &      jr_last(nproma,nspec,nblks)
    cf_desc    = t_cf_var(name_head//'%jr_last', 'N/A', 'last-launched ray index', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%jr_last', p_var%jr_last, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape3d_spec, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_spec, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! &      jk_full_rtop(nproma,nray_nspec,nblks)
    cf_desc    = t_cf_var(name_head//'%jk_full_rtop', 'N/A', 'full lev index of ray-volume top', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%jk_full_rtop', p_var%jk_full_rtop, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! &      jk_full_rbot(nproma,nray_nspec,nblks)
    cf_desc    = t_cf_var(name_head//'%jk_full_rbot', 'N/A', 'full lev index of ray-volume bottom', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%jk_full_rbot', p_var%jk_full_rbot, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! &      jk_half_rtop(nproma,nray_nspec,nblks)
    cf_desc    = t_cf_var(name_head//'%jk_half_rtop', 'N/A', 'half lev index of ray-volume top', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%jk_half_rtop', p_var%jk_half_rtop, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! &      jk_half_rbot(nproma,nray_nspec,nblks)
    cf_desc    = t_cf_var(name_head//'%jk_half_rbot', 'N/A', 'half lev index of ray-volume bottom', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%jk_half_rbot', p_var%jk_half_rbot, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

  END IF

  IF ( name_head == 'src_conv' ) THEN

    ! &      jk_source(nproma,nblks)
    cf_desc    = t_cf_var(name_head//'%jk_source', 'N/A', 'source level index', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%jk_source', p_var%jk_source, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
    !           & ldims=shape2d_grid, in_group=groups("ray_prog_vars"),     &
                & ldims=shape2d_grid, &
                & lcontainer=.FALSE., lrestart=.FALSE., loutput=.FALSE.)

  END IF

END SUBROUTINE new_ray_list
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE new_spl_list( k_jg, nray, kblks, listname, p_list, p_var, name_head )

  TYPE(t_var_list_ptr), INTENT(INOUT) :: p_list
  TYPE(t_spl),          INTENT(INOUT) :: p_var
  INTEGER,              INTENT(IN)    :: k_jg, kblks
  INTEGER,              INTENT(IN)    :: nray
  CHARACTER(len=*),     INTENT(IN)    :: listname
  CHARACTER(len=*),     INTENT(IN)    :: name_head

  ! Local variables

  TYPE(t_cf_var)    ::    cf_desc
  TYPE(t_grib2_var) :: grib2_desc

  INTEGER           :: shape3d_ray(3)
  INTEGER           :: ibits
  INTEGER           :: datatype_flt

  IF ( lnetcdf_flt64_output ) THEN
    datatype_flt = DATATYPE_FLT64
  ELSE
    datatype_flt = DATATYPE_FLT32
  ENDIF

  ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

  shape3d_ray = (/nproma, nray, kblks/)

  CALL vlr_add( p_list, TRIM(listname), patch_id=k_jg, lrestart=.FALSE. )

  !--------------------------------
  ! Add ray properties as variables
  !--------------------------------

  IF ( name_head == 'spl' ) THEN

    ! &      mode(nproma,nray,nblks)
    cf_desc    = t_cf_var(name_head//'%mode', 'N/A', 'horizontal split mode', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%mode', p_var%mode,          &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
!               & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! &      rot_c(nproma,nray,nblks)
    cf_desc    = t_cf_var(name_head//'%rot_c', 'N/A', 'horizontal split rot_c', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%rot_c', p_var%rot_c,          &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
!               & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! &      rot_s(nproma,nray,nblks)
    cf_desc    = t_cf_var(name_head//'%rot_s', 'N/A', 'horizontal split rot_s', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%rot_s', p_var%rot_s,          &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
!               & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! &      aux(nproma,nray,nblks)
    cf_desc    = t_cf_var(name_head//'%aux', 'N/A', 'horizontal split aux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, name_head//'%aux', p_var%aux,          &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
!               & ldims=shape3d_ray, in_group=groups("ray_prog_vars"),     &
                & ldims=shape3d_ray, &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ALLOCATE( p_var%jcol(2,nproma,1:nray,kblks) )

  END IF

END SUBROUTINE new_spl_list
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE new_lfluxbg_list( k_jg, kblks, listname, p_list, p_var )

  TYPE(t_var_list_ptr), INTENT(INOUT) :: p_list
  TYPE(t_lfluxbg),      INTENT(INOUT) :: p_var
  INTEGER,              INTENT(IN)    :: k_jg, kblks !< dimension sizes
  CHARACTER(len=*),     INTENT(IN)    :: listname

  ! Local variables

  TYPE(t_cf_var)    ::    cf_desc
  TYPE(t_grib2_var) :: grib2_desc

  INTEGER           :: shape2d_grid(2)
  INTEGER           :: ibits
  INTEGER           :: datatype_flt

  IF ( lnetcdf_flt64_output ) THEN
    datatype_flt = DATATYPE_FLT64
  ELSE
    datatype_flt = DATATYPE_FLT32
  ENDIF

  ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

  shape2d_grid = (/nproma, kblks/)

  CALL vlr_add( p_list, TRIM(listname), patch_id=k_jg, lrestart=.FALSE. )

  !---------------------------------------
  ! Add launch flux profiles as variables
  !---------------------------------------

  ! &      lat_prof_bw(nproma,nblks)
  cf_desc    = t_cf_var('lflux%lat_prof_bw', 'lat_prof_bw', &
    &                   'latitudinal profile of GW flux factor (boreal winter)', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'lflux%lat_prof_bw', p_var%lat_prof_bw, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
              & ldims=shape2d_grid,                                        &
              & lrestart=.FALSE., loutput=.FALSE.)

  ! &      lat_prof_bs(nproma,nblks)
  cf_desc    = t_cf_var('lflux%lat_prof_bs', 'lat_prof_bs', &
    &                   'latitudinal profile of GW flux factor (boreal summer)', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'lflux%lat_prof_bs', p_var%lat_prof_bs, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
              & ldims=shape2d_grid,                                        &
              & lrestart=.FALSE., loutput=.FALSE.)

  ! &      lat_prof(nproma,nblks)
  cf_desc    = t_cf_var('lflux%lat_prof', 'lat_prof', &
    &                   'latitudinal profile of GW flux factor', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'lflux%lat_prof', p_var%lat_prof, &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
              & ldims=shape2d_grid,                                        &
              & lrestart=.FALSE., loutput=.FALSE.)

END SUBROUTINE new_lfluxbg_list
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE new_field_list( k_jg, klev, kblks, listname, p_list, p_var )

  TYPE(t_var_list_ptr), INTENT(INOUT) :: p_list
  TYPE(t_msgwam),       INTENT(INOUT) :: p_var
  INTEGER,              INTENT(IN)    :: k_jg, klev, kblks !< dimension sizes
  CHARACTER(len=*),     INTENT(IN)    :: listname

  ! Local variables

  TYPE(t_cf_var)    ::    cf_desc
  TYPE(t_grib2_var) :: grib2_desc

  INTEGER           :: shape3d(3), shape3dkp1(3)
  INTEGER           :: shape2d(2)
  INTEGER           :: ibits
  INTEGER           :: datatype_flt

  CHARACTER(LEN=1)  :: dir4p1(5)
  CHARACTER(LEN=32) :: str
  INTEGER           :: nsrc
  INTEGER           :: j, js

  IF ( lnetcdf_flt64_output ) THEN
    datatype_flt = DATATYPE_FLT64
  ELSE
    datatype_flt = DATATYPE_FLT32
  ENDIF

  ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

  shape3d      = (/nproma, klev  , kblks/)
  shape3dkp1   = (/nproma, klev+1, kblks/)
  shape2d      = (/nproma, kblks/)

  CALL vlr_add( p_list, TRIM(listname), patch_id=k_jg, lrestart=.TRUE. )

  !-------------------------------
  ! add flux profiles as variables
  !-------------------------------

  ! Pseudomomentum fluxes
  !           apmfl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('apmfl_mgm', 'Pa', 'absolute pseudomomentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'apmfl_mgm'  , p_var%apmfl_mgm,         &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
    & ldims=shape3dkp1, lrestart=.FALSE., loutput=.TRUE.)

  ! Momentum fluxes
  !           amfl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('amfl_mgm', 'Pa', 'absolute momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'amfl_mgm'  , p_var%amfl_mgm,           &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           uufl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('uufl_mgm', 'Pa', 'uu momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'uufl_mgm', p_var%uufl_mgm,             &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           uvfl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('uvfl_mgm', 'Pa', 'uv momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'uvfl_mgm', p_var%uvfl_mgm,             &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           uwfl_mgm(nproma,nlevp1,nblks_c)
  cf_desc    = t_cf_var('uwfl_mgm', 'Pa', 'uw momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'uwfl_mgm', p_var%uwfl_mgm,             &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
    & ldims=shape3dkp1, lrestart=.FALSE., loutput=.TRUE.)

  !           vufl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('vufl_mgm', 'Pa', 'vu momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'vufl_mgm', p_var%vufl_mgm,             &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           vvfl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('vvfl_mgm', 'Pa', 'vv momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'vvfl_mgm', p_var%vvfl_mgm,             &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           vwfl_mgm(nproma,nlevp1,nblks_c)
  cf_desc    = t_cf_var('vwfl_mgm', 'Pa', 'vw momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'vwfl_mgm', p_var%vwfl_mgm,             &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
    & ldims=shape3dkp1, lrestart=.FALSE., loutput=.TRUE.)

  ! Momentum fluxes (directional)
  !           mfl_mgm_e(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('mfl_mgm_E', 'Pa', 'eastward momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'mfl_mgm_e', p_var%mfl_mgm_e,           &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           mfl_mgm_w(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('mfl_mgm_W', 'Pa', 'westward momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'mfl_mgm_w', p_var%mfl_mgm_w,           &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           mfl_mgm_n(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('mfl_mgm_N', 'Pa', 'northward momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'mfl_mgm_n', p_var%mfl_mgm_n,           &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           mfl_mgm_s(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('mfl_mgm_S', 'Pa', 'southward momentum fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'mfl_mgm_s', p_var%mfl_mgm_s,           &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! wave action fluxes (directional)
  ! wafl_mgm_u(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('wafl_mgm_U', 'Pa m', 'upward wave action fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'wafl_mgm_u', p_var%wafl_mgm_u,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
        & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! wafl_mgm_d(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('wafl_mgm_D', 'Pa m', 'dowmward wave action fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'wafl_mgm_d', p_var%wafl_mgm_d,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
        & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! wafl_mgm_e(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('wafl_mgm_E', 'Pa m', 'eastward wave action fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'wafl_mgm_e', p_var%wafl_mgm_e,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
        & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! wafl_mgm_w(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('wafl_mgm_W', 'Pa m', 'westward wave action fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'wafl_mgm_w', p_var%wafl_mgm_w,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
        & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! wafl_mgm_n(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('wafl_mgm_N', 'Pa m', 'northward wave action fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'wafl_mgm_n', p_var%wafl_mgm_n,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
        & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! wafl_mgm_s(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('wafl_mgm_S', 'Pa m', 'southward wave action fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'wafl_mgm_s', p_var%wafl_mgm_s,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
        & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! Potential temperature fluxes
  !           aptfl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('aptfl_mgm', 'Pa K s/m', 'absolute pot temp fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'aptfl_mgm'  , p_var%aptfl_mgm,         &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           utfl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('utfl_mgm', 'Pa K s/m', 'zonal pot temp fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'utfl_mgm', p_var%utfl_mgm,             &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           vtfl_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('vtfl_mgm', 'Pa K s/m', 'meridional pot temp fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'vtfl_mgm', p_var%vtfl_mgm,             &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! Potential temperature fluxes (directional)
  !           ptfl_mgm_e(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ptfl_mgm_E', 'Pa K s/m', 'eastward pot temp fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ptfl_mgm_e', p_var%ptfl_mgm_e,         &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           ptfl_mgm_w(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ptfl_mgm_W', 'Pa K s/m', 'westward pot temp fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ptfl_mgm_w', p_var%ptfl_mgm_w,         &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           ptfl_mgm_n(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ptfl_mgm_N', 'Pa K s/m', 'northward pot temp fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ptfl_mgm_n', p_var%ptfl_mgm_n,         &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  !           ptfl_mgm_s(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ptfl_mgm_S', 'Pa K s/m', 'southward pot temp fluxes', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ptfl_mgm_s', p_var%ptfl_mgm_s,         &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! Energy
  !           energy_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('energy_mgm', 'kg/m/s**2', 'gw energy', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'energy_mgm', p_var%energy_mgm,         &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! Potential energy
  !           energy_p_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('energy_p_mgm', 'kg/m/s**2', 'gw potential energy', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'energy_p_mgm', p_var%energy_p_mgm,     &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

  ! Potential energy
  !           action_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('action_mgm', 'kg/m/s', 'gw action', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'action_mgm', p_var%action_mgm,     &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

!!!!! TEST !!!!! TO BE CLEANED LATER
  ! mom flux convergence
  !           mfcxz_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('mfcxz_mgm', 'ms-1', 'mfcxz', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'mfcxz_mgm', p_var%mfcxz_mgm,           &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)
  !           mfcyz_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('mfcyz_mgm', 'ms-1', 'mfcyz', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'mfcyz_mgm', p_var%mfcyz_mgm,           &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)
  ! elastic terms
  !           etx_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('etx_mgm', 'ms-1', 'etx', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'etx_mgm', p_var%etx_mgm,               &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)
  !           ety_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ety_mgm', 'ms-1', 'ety', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ety_mgm', p_var%ety_mgm,               &
    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)
!!!!! TEST !!!!! TO BE CLEANED LATER

  ! Convective source related
  !       gwd_conv_u(nproma,nlev,nblks)
  cf_desc    = t_cf_var('gwd_conv_u', 'm s-2', 'GWD tendency of zonal wind (conv)', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'gwd_conv_u', p_var%gwd_conv_u,                    &
              & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
              & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       gwd_conv_v(nproma,nlev,nblks)
  cf_desc    = t_cf_var('gwd_conv_v', 'm s-2', 'GWD tendency of meridional wind (conv)', datatype_flt)
  grib2_desc = grib2_var(192, 128, 221, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'gwd_conv_v', p_var%gwd_conv_v,                    &
              & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
              & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       gwd_conv_t(nproma,nlev,nblks)
  cf_desc    = t_cf_var('gwd_conv_t', 'm s-2', 'GWD tendency of pot temp (conv)', datatype_flt)
  grib2_desc = grib2_var(192, 128, 221, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'gwd_conv_t', p_var%gwd_conv_t,                    &
              & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
              & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       entr_cgw(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('entr_cgw'  , 'W K-1 kg-1' , 'entropy due to convection', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'entr_cgw'  , p_var%entr_cgw  ,                    &
    &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    &  ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       heat_cgw(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('heat_cgw'  , 'W kg-1' , 'heating due to convection', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'heat_cgw'  , p_var%heat_cgw  ,                    &
    &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    &  ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       tupd_cgw(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('tupd_cgw'  , 'K'      , 'temperature of updraft parcel', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'tupd_cgw'  , p_var%tupd_cgw  ,                    &
    &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    &  ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       mbas_con_cgw(nproma,nblks_c)
  cf_desc    = t_cf_var('mbas_con_cgw', '', 'cloud base level index', datatype_flt)
  grib2_desc = grib2_var(0, 6, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'mbas_con_cgw', p_var%mbas_con_cgw,                &
    &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
    &           ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )

  !       mtop_con_cgw(nproma,nblks_c)
  cf_desc    = t_cf_var('mtop_con_cgw', '', 'cloud top level index', datatype_flt)
  grib2_desc = grib2_var(0, 6, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'mtop_con_cgw', p_var%mtop_con_cgw,                &
    &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
    &           ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )

  !       ktype_cgw(nproma,nblks_c)
  cf_desc    = t_cf_var('ktype_cgw', '', 'type of convection', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ktype_cgw', p_var%ktype_cgw,                      &
    &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
    &           ldims=shape2d, lrestart=.FALSE., loutput=.TRUE. )

  !       test_1_cgw(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('test_1_cgw', 'x' , 'cgw test variable', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'test_1_cgw', p_var%test_1_cgw,                    &
    &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    &  ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       test_2_cgw(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('test_2_cgw', 'x' , 'cgw test variable', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'test_2_cgw', p_var%test_2_cgw,                    &
    &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    &  ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       test_cgw(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('test_cgw'  , 'x' , 'cgw test variable', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'test_cgw'  , p_var%test_cgw  ,                    &
    &  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc,     &
    &  ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !       flag_cgw(nproma,nblks_c)
  cf_desc    = t_cf_var('flag_cgw', '', 'cgw flag', datatype_flt)
  grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'flag_cgw', p_var%flag_cgw,                        &
    &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
    &           ldims=shape2d, lrestart=.TRUE., loutput=.TRUE., initval=-1)

  IF ( gws_conv_config%n_source(k_jg) > 0 ) THEN

    nsrc = gws_conv_config%n_source(k_jg)

    !       pars_cgw(nproma,nblks_c,13+nsrc*2)
    cf_desc    = t_cf_var('pars_cgw', '', 'cgw parameters' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'pars_cgw', p_var%pars_cgw,                        &
      &  GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      &  ldims=(/nproma,kblks,13+nsrc*2/), lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE. )
    ALLOCATE(p_var%pars_cgw_ptr(13+nsrc*2))
    DO j = 1, 13 + nsrc*2
      WRITE(str,'(i2)') j
      CALL add_ref( p_list, 'pars_cgw', 'pars_cgw_'//TRIM(ADJUSTL(str)),                 &
        &  p_var%pars_cgw_ptr(j)%p_2d, GRID_UNSTRUCTURED_CELL, ZA_SURFACE,               &
        &  t_cf_var('pars_cgw_'//TRIM(ADJUSTL(str)), '', 'cgw parameter', datatype_flt), &
        &  grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),                &
        &  ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                              &
        &  opt_var_ref_pos=3, ref_idx=j )
    ENDDO

    !       ctmfl_cgw(nproma,nblks_c,5+nsrc*5)
    cf_desc    = t_cf_var('ctmfl_cgw', 'Pa',                             &
      &                   'cloud-top cgw momentum fluxes, EWNS+A', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ctmfl_cgw', p_var%ctmfl_cgw,                  &
      &  GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
      &  ldims=(/nproma,kblks,5+nsrc*5/), lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE. )
    ALLOCATE(p_var%ctmfl_cgw_ptr(5+nsrc*5))
    dir4p1 = (/'E','W','N','S','a'/)
    DO js = 0, nsrc
      DO j = 1, 5
        IF (js == 0) THEN
          WRITE(str,'(a)') dir4p1(j)
        ELSE
          WRITE(str,'(a,i1)') dir4p1(j), js
        END IF
        CALL add_ref( p_list, 'ctmfl_cgw', 'ctmfl_cgw_'//TRIM(ADJUSTL(str)),           &
          &  p_var%ctmfl_cgw_ptr(js*5+j)%p_2d, GRID_UNSTRUCTURED_CELL, ZA_SURFACE,     &
          &  t_cf_var('ctmfl_cgw_'//TRIM(ADJUSTL(str)), 'Pa',                          &
          &     'cloud-top cgw momentum fluxes, '//TRIM(ADJUSTL(str)), datatype_flt),  &
          &  grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),            &
          &  ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                          &
          &  opt_var_ref_pos=3, ref_idx=js*5+j )
      ENDDO
    ENDDO

    !       apmfl_cgw(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('apmfl_cgw', 'Pa', 'cgw absolute pseudomomentum fluxes' , datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'apmfl_cgw', p_var%apmfl_cgw,                           &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
      &           ldims=shape3dkp1, lrestart=.FALSE., loutput=.TRUE. )

    !       amfl_cgw(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('amfl_cgw', 'Pa', 'cgw absolute momentum fluxes' , datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'amfl_cgw', p_var%amfl_cgw,                             &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       mfl_cgw_e(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('cgw_mfl_E', 'Pa', 'cgw eastward momentum fluxes' , datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'mfl_cgw_e', p_var%mfl_cgw_e,                           &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       mfl_cgw_w(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('cgw_mfl_W', 'Pa', 'cgw westward momentum fluxes' , datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'mfl_cgw_w', p_var%mfl_cgw_w,                           &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       mfl_cgw_n(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('cgw_mfl_N', 'Pa', 'cgw northward momentum fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'mfl_cgw_n', p_var%mfl_cgw_n,                           &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       mfl_cgw_s(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('cgw_mfl_S', 'Pa', 'cgw southward momentum fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'mfl_cgw_s', p_var%mfl_cgw_s,                           &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    ! wafl_cgw_u(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('wafl_cgw_U', 'Pa m', 'upward wave action fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'wafl_cgw_u', p_var%wafl_cgw_u,           &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
          & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

    ! wafl_cgw_d(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('wafl_cgw_D', 'Pa m', 'dowmward wave action fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'wafl_cgw_d', p_var%wafl_cgw_d,           &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
          & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

    ! wafl_cgw_e(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('wafl_cgw_E', 'Pa m', 'eastward wave action fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'wafl_cgw_e', p_var%wafl_cgw_e,           &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
          & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

    ! wafl_cgw_w(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('wafl_cgw_W', 'Pa m', 'westward wave action fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'wafl_cgw_w', p_var%wafl_cgw_w,           &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
          & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

    ! wafl_cgw_n(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('wafl_cgw_N', 'Pa m', 'northward wave action fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'wafl_cgw_n', p_var%wafl_cgw_n,           &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
          & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

    ! wafl_cgw_s(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('wafl_cgw_S', 'Pa m', 'southward wave action fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'wafl_cgw_s', p_var%wafl_cgw_s,           &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,     &
          & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE.)

    !       aptfl_cgw(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('aptfl_cgw', 'Pa K s/m', 'cgw absolute pot temp fluxes' , datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'aptfl_cgw', p_var%aptfl_cgw,                           &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       ptfl_cgw_e(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('cgw_ptfl_E', 'Pa K s/m', 'cgw eastward pot temp fluxes' , datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ptfl_cgw_e', p_var%ptfl_cgw_e,                         &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       ptfl_cgw_w(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('cgw_ptfl_W', 'Pa K s/m', 'cgw westward pot temp fluxes' , datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ptfl_cgw_w', p_var%ptfl_cgw_w,                         &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       ptfl_cgw_n(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('cgw_ptfl_N', 'Pa K s/m', 'cgw northward pot temp fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ptfl_cgw_n', p_var%ptfl_cgw_n,                         &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       ptfl_cgw_s(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('cgw_ptfl_S', 'Pa K s/m', 'cgw southward pot temp fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ptfl_cgw_s', p_var%ptfl_cgw_s,                         &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       energy_cgw(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('energy_cgw', 'Pa', 'cgw energy', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'energy_cgw', p_var%energy_cgw,                         &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       energy_p_cgw(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('energy_p_cgw', 'Pa', 'cgw potential energy', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'energy_p_cgw', p_var%energy_p_cgw,                     &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       action_cgw(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('action_cgw', 'kg/m/s', 'cgw action', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'action_cgw', p_var%action_cgw,                         &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE     , cf_desc, grib2_desc, &
      &           ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

    !       ll_k_cgw (nproma,nblks_c)
    cf_desc    = t_cf_var('ll_k_cgw' , '', 'cgw ll k'          , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'll_k_cgw' , p_var%ll_k_cgw,                      &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,  &
      &           ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    !       ml_k_cgw (nproma,nblks_c)
    cf_desc    = t_cf_var('ml_k_cgw' , '', 'cgw ml k'          , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ml_k_cgw' , p_var%ml_k_cgw,                      &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,  &
      &           ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    !       ll_dz_cgw(nproma,nblks_c)
    cf_desc    = t_cf_var('ll_dz_cgw', '', 'cgw ll dz'         , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'll_dz_cgw', p_var%ll_dz_cgw,                     &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,  &
      &           ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    !       ml_dz_cgw(nproma,nblks_c)
    cf_desc    = t_cf_var('ml_dz_cgw', '', 'cgw ml dz'         , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ml_dz_cgw', p_var%ml_dz_cgw,                     &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,  &
      &           ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    !       ctmfl_specc_cgw(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('ctmfl_specc_cgw', 'Pa', 'cgw absolute momentum fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ctmfl_specc_cgw', p_var%ctmfl_specc_cgw,               &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
      &           ldims=shape3dkp1, lrestart=.FALSE., loutput=.TRUE. )

    !       ctmfl_spec0_cgw(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('ctmfl_spec0_cgw', 'Pa', 'cgw absolute momentum fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ctmfl_spec0_cgw', p_var%ctmfl_spec0_cgw,               &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
      &           ldims=shape3dkp1, lrestart=.FALSE., loutput=.TRUE. )

    !       ctmfl_spec1_cgw(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('ctmfl_spec1_cgw', 'Pa', 'cgw absolute momentum fluxes', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_list, 'ctmfl_spec1_cgw', p_var%ctmfl_spec1_cgw,               &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
      &           ldims=shape3dkp1, lrestart=.FALSE., loutput=.TRUE. )

  ENDIF ! gws_conv_config%n_source(k_jg)

  ! Tendencies from MS-GWaM
  !        ddt_u_gwd_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ddt_u_gwd_mgm', 'm s-2', &
       &                            'GWD tendency of zonal wind', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ddt_u_gwd_mgm', p_var%ddt_u_gwd_mgm,              &
              & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
              & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !        ddt_v_gwd_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ddt_v_gwd_mgm', 'm s-2', &
       &                            'GWD tendency of meridional wind', datatype_flt)
  grib2_desc = grib2_var(192, 128, 221, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ddt_v_gwd_mgm', p_var%ddt_v_gwd_mgm,              &
              & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
              & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !        ddt_t_gwd_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ddt_t_gwd_mgm', 'K s-1', &
       &                            'GWD tendency of temp', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ddt_t_gwd_mgm', p_var%ddt_t_gwd_mgm,              &
              & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
              & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

  !        ddt_pt_gwd_mgm(nproma,nlev,nblks_c)
  cf_desc    = t_cf_var('ddt_pt_gwd_mgm', 'K s-1', &
       &                            'GWD tendency of pot temp', datatype_flt)
  grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
  CALL add_var( p_list, 'ddt_pt_gwd_mgm', p_var%ddt_pt_gwd_mgm,            &
              & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
              & ldims=shape3d, lrestart=.FALSE., loutput=.TRUE. )

END SUBROUTINE new_field_list
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE msgwam_write_restartfiles( mtime_current, p_patch )

  TYPE(datetime), POINTER            :: mtime_current   !< current datetime (mtime)
  TYPE(t_patch) , TARGET, INTENT(IN) :: p_patch

  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=1024) :: r_filename

  INTEGER , ALLOCATABLE :: gidx(:)
  REAL(wp), ALLOCATABLE :: glat(:), glon(:)

  INTEGER , ALLOCATABLE :: nidx(:), jcol1(:)

  INTEGER  :: rl_start, rl_end
  INTEGER  :: i_startblk, i_endblk  ! blocks
  INTEGER  :: i_startidx, i_endidx  ! slices
  INTEGER  :: jg
  INTEGER  :: my_work_id, my_work_size
  INTEGER  :: ncol

  ! netCDF-related variables
  INTEGER  :: ncid, ncvarid(99), ncstart(2), nccount(2), ncimap(2),  &
    &         ncdimid(7), ncrecid, ncdimid_v(2)
  INTEGER, PARAMETER ::  ncstri(2) = (/1,1/)

  INTEGER  :: jc,jb

  jg = p_patch%id

  IF ( .NOT. my_process_is_work() )  RETURN
  IF ( lsteady )  RETURN

  CALL message('', ' ')
  CALL message('msgwam_write_restartfiles',  &
    &          'MS-GWaM: production of restart files for rays')

  my_work_id   = get_my_mpi_work_id()
  my_work_size = get_my_mpi_work_comm_size()

  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

  ! file name
  WRITE(r_filename,'(a,i0,5(a,i2.2),a,i2.2)')                             &
    &    'msgwam_ray_restart_', mtime_current%date%year,                  &
    &    '-', mtime_current%date%month,  '-', mtime_current%date%day,     &
    &    '_', mtime_current%time%hour,   ':', mtime_current%time%minute,  &
    &    ':', mtime_current%time%second, '_D', jg
  IF (TRIM(dir_restartfiles) /= '')  &
    &  r_filename = TRIM(dir_restartfiles)//'/'//TRIM(r_filename)
  IF (my_work_size == 1) THEN
    r_filename = TRIM(r_filename)//'.nc'
  ELSE
    WRITE(r_filename,'(a,i6.6,a)') TRIM(r_filename)//'_P', my_work_id, '.nc'
  END IF
  IF ( ltest_restart(2) )  r_filename = TRIM(r_filename)//'_test'

  ! number of grid columns  &  1-D column index for (i_startidx,jb)

  ALLOCATE( nidx(i_startblk:i_endblk), jcol1(i_startblk:i_endblk) )

  DO jb = i_startblk, i_endblk
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
      &                 i_startidx, i_endidx, rl_start, rl_end )
    nidx(jb) = MAX(0, i_endidx - i_startidx + 1)
  ENDDO
  ncol = SUM(nidx(:))

  jcol1(i_startblk) = 1
  DO jb = i_startblk+1, i_endblk
    jcol1(jb) = jcol1(jb-1) + nidx(jb-1)
  ENDDO

  ! The below will not happen, at least in the current ICON.
  IF ( ANY(nidx(i_startblk:i_endblk) == 0) )  &
    &  CALL finish('msgwam_write_restartfiles', 'nidx == 0 somewhere')

  !---  Create the file  &  Define dimensions/variables  -----------------------

  CALL nf( nf_create(TRIM(r_filename), NF_64BIT_OFFSET, ncid) )
! CALL nf( nf_create(TRIM(r_filename), IOR(NF_NOCLOBBER, NF_64BIT_OFFSET), ncid) )

  CALL nf( nf_def_dim(ncid, 'ray', nrays(jg), ncdimid(1)) )

  CALL nf( nf_def_dim(ncid, 'spec', nrays_add(jg), ncdimid(2)) )

  CALL nf( nf_def_dim(ncid, 'col', NF_UNLIMITED, ncrecid) )

  ! [0]  grid index and lat/lon

  CALL nf( nf_def_var(ncid, 'gidx', NF_INT   , 1, (/ncrecid/), ncvarid(iv_gidx)) )
  CALL nf( nf_def_var(ncid, 'glat', NF_DOUBLE, 1, (/ncrecid/), ncvarid(iv_glat)) )
  CALL nf( nf_def_var(ncid, 'glon', NF_DOUBLE, 1, (/ncrecid/), ncvarid(iv_glon)) )

  ! [1]  ray variables

  ncdimid_v = (/ ncdimid(1), ncrecid /)

  CALL nf( nf_def_var(ncid, 'lon'      , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_lon  )) )
  CALL nf( nf_def_var(ncid, 'lat'      , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_lat  )) )
  CALL nf( nf_def_var(ncid, 'z'        , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_z    )) )
  CALL nf( nf_def_var(ncid, 'dlon'     , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_dlon )) )
  CALL nf( nf_def_var(ncid, 'dlat'     , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_dlat )) )
  CALL nf( nf_def_var(ncid, 'dz'       , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_dz   )) )
  CALL nf( nf_def_var(ncid, 'k'        , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_k    )) )
  CALL nf( nf_def_var(ncid, 'l'        , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_l    )) )
  CALL nf( nf_def_var(ncid, 'm'        , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_m    )) )
  CALL nf( nf_def_var(ncid, 'dk'       , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_dk   )) )
  CALL nf( nf_def_var(ncid, 'dl'       , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_dl   )) )
  CALL nf( nf_def_var(ncid, 'dm'       , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_dm   )) )
  CALL nf( nf_def_var(ncid, 'wadens'   , NF_DOUBLE, 2, ncdimid_v, ncvarid(iv_wad  )) )

    ! (The variable 'area' is not saved but will be newly calculated.)

  CALL nf( nf_def_var(ncid, 'iexist'   , NF_INT , 2, ncdimid_v, ncvarid(iv_iex)) )
  CALL nf( nf_def_var(ncid, 'specid'   , NF_INT , 2, ncdimid_v, ncvarid(iv_sid)) )
  CALL nf( nf_def_var(ncid, 'jk_active', NF_INT , 2, ncdimid_v, ncvarid(iv_kac)) )

  ! [2]  (spectral) source variables

  ncdimid_v = (/ ncdimid(2), ncrecid /)
  CALL nf( nf_def_var(ncid, 'spec_jr_last', NF_INT, 2, ncdimid_v, ncvarid(iv_jrl)) )

  !---  Finalize the Def mode  -------------------------------------------------

  CALL nf( nf_enddef(ncid) )

  !---  Write variables  -------------------------------------------------------

  ! [0]  index and lat/lon

  ALLOCATE( gidx(ncol), glat(ncol), glon(ncol) )

  DO jb = i_startblk, i_endblk
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
      &                 i_startidx, i_endidx, rl_start, rl_end )
    DO jc = i_startidx, i_endidx
      gidx(jcol1(jb)+jc-i_startidx) = (jb-1)*nproma + jc
      glat(jcol1(jb)+jc-i_startidx) = p_patch%cells%center(jc,jb)%lat
      glon(jcol1(jb)+jc-i_startidx) = p_patch%cells%center(jc,jb)%lon
    ENDDO
  ENDDO

  CALL nf( nf_put_vara_int   (ncid, ncvarid(iv_gidx), (/1/), (/ncol/), gidx) )
  CALL nf( nf_put_vara_double(ncid, ncvarid(iv_glat), (/1/), (/ncol/), glat) )
  CALL nf( nf_put_vara_double(ncid, ncvarid(iv_glon), (/1/), (/ncol/), glon) )

  DEALLOCATE( gidx, glat, glon )

  ! [1]  ray variables

  ncstart(1) = 1   ;   nccount(1) = nrays(jg)
  ncimap(2) = 1

  DO jb = i_startblk, i_endblk
!   IF (nidx(jb) == 0)  CYCLE
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
      &                 i_startidx, i_endidx, rl_start, rl_end )

    ncstart(2) = jcol1(jb)
    nccount(2) = nidx (jb)
    ncimap (1) = nidx (jb)

    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_lon ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%lon       (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_lat ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%lat       (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_z   ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%z         (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_dlon), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dlon      (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_dlat), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dlat      (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_dz  ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dz        (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_k   ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%k         (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_l   ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%l         (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_m   ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%m         (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_dk  ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dk        (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_dl  ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dl        (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_dm  ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dm        (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_double(ncid, ncvarid(iv_wad ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%wadens    (i_startidx:i_endidx,:,jb)) )

    CALL nf( nf_put_varm_int(ncid, ncvarid(iv_iex), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%iexist    (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_int(ncid, ncvarid(iv_sid), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%specid    (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_put_varm_int(ncid, ncvarid(iv_kac), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%jk_active (i_startidx:i_endidx,:,jb)) )

  ENDDO

  ! [2]  (spectral) source variables

  ncstart(1) = 1   ;   nccount(1) = nrays_add(jg)
  DO jb = i_startblk, i_endblk
!   IF (nidx(jb) == 0)  CYCLE
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
      &                 i_startidx, i_endidx, rl_start, rl_end )
    ncstart(2) = jcol1(jb)
    nccount(2) = nidx (jb)
    ncimap (1) = nidx (jb)
    CALL nf( nf_put_varm_int(ncid, ncvarid(iv_jrl), ncstart, nccount, ncstri, ncimap,  &
      &                      p_ray(jg)%jr_last(i_startidx:i_endidx,:,jb)) )
  ENDDO

  !---  Close the file and finalize  -------------------------------------------

  CALL nf( nf_close(ncid) )

  DEALLOCATE( nidx, jcol1 )

  ltest_restart(2) = .FALSE.   ! flag to test restart-file reading

  CALL message('msgwam_write_restartfiles', 'MS-GWaM: finished writing')
  CALL message('', ' ')

END SUBROUTINE msgwam_write_restartfiles
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE read_restartfile_specified( ncid, ncvarid, jcol_1st, p_patch,  &
        &                              jb1, jc1, jb2, jc2, l_first_read )

  TYPE(t_patch) , TARGET, INTENT(IN   ) :: p_patch
  INTEGER ,               INTENT(IN   ) :: ncid
  INTEGER ,               INTENT(IN   ) :: jcol_1st
  INTEGER ,               INTENT(IN   ) :: jb1, jc1, jb2, jc2
  LOGICAL ,               INTENT(IN   ) :: l_first_read
  INTEGER ,               INTENT(INOUT) :: ncvarid(:)

  INCLUDE 'netcdf.inc'

  INTEGER  :: nidx(jb1:jb2), jcol1(jb1:jb2)

  INTEGER  :: rl_start, rl_end
  INTEGER  :: i_startblk, i_endblk  ! blocks
  INTEGER  :: i_startidx, i_endidx  ! slices
  INTEGER  :: jg

  ! netCDF-related variables
  INTEGER  :: ncstart(2), nccount(2), ncimap(2)
  INTEGER, PARAMETER ::  ncstri(2) = (/1,1/)

  INTEGER  :: jb

  jg = p_patch%id

  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int - 2

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

  ! number of grid columns  &  1-D column index for (i_startidx,jb)

  DO jb = jb1, jb2
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
      &                 i_startidx, i_endidx, rl_start, rl_end )
    IF (jb == jb1)  i_startidx = jc1
    IF (jb == jb2)  i_endidx   = jc2
    nidx(jb) = MAX(0, i_endidx - i_startidx + 1)
  ENDDO

  jcol1(jb1) = jcol_1st
  DO jb = jb1+1, jb2
    jcol1(jb) = jcol1(jb-1) + nidx(jb-1)
  ENDDO

  ! The below will not happen, at least in the current ICON.
  IF ( ANY(nidx(jb1:jb2) == 0) )  &
    &  CALL finish('msgwam_read_restartfiles', 'nidx == 0 somewhere')

  !---  Read Var IDs  ----------------------------------------------------------

  IF ( l_first_read ) THEN

    ! [1]  ray variables

    CALL nf( nf_inq_varid(ncid, 'lon'      , ncvarid(iv_lon )) )
    CALL nf( nf_inq_varid(ncid, 'lat'      , ncvarid(iv_lat )) )
    CALL nf( nf_inq_varid(ncid, 'z'        , ncvarid(iv_z   )) )
    CALL nf( nf_inq_varid(ncid, 'dlon'     , ncvarid(iv_dlon)) )
    CALL nf( nf_inq_varid(ncid, 'dlat'     , ncvarid(iv_dlat)) )
    CALL nf( nf_inq_varid(ncid, 'dz'       , ncvarid(iv_dz  )) )
    CALL nf( nf_inq_varid(ncid, 'k'        , ncvarid(iv_k   )) )
    CALL nf( nf_inq_varid(ncid, 'l'        , ncvarid(iv_l   )) )
    CALL nf( nf_inq_varid(ncid, 'm'        , ncvarid(iv_m   )) )
    CALL nf( nf_inq_varid(ncid, 'dk'       , ncvarid(iv_dk  )) )
    CALL nf( nf_inq_varid(ncid, 'dl'       , ncvarid(iv_dl  )) )
    CALL nf( nf_inq_varid(ncid, 'dm'       , ncvarid(iv_dm  )) )
    CALL nf( nf_inq_varid(ncid, 'wadens'   , ncvarid(iv_wad )) )

      ! (The variable 'area' is not saved but will be newly calculated.)

    CALL nf( nf_inq_varid(ncid, 'iexist'   , ncvarid(iv_iex)) )
    CALL nf( nf_inq_varid(ncid, 'specid'   , ncvarid(iv_sid)) )
    CALL nf( nf_inq_varid(ncid, 'jk_active', ncvarid(iv_kac)) )

    ! [2]  (spectral) source variables

    CALL nf( nf_inq_varid(ncid, 'spec_jr_last', ncvarid(iv_jrl)) )

  END IF

  !---  Read variables  --------------------------------------------------------

  ! [1]  ray variables

  ncstart(1) = 1   ;   nccount(1) = nrays(jg)
  ncimap(2) = 1

  DO jb = jb1, jb2
!   IF (nidx(jb) == 0)  CYCLE
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
      &                 i_startidx, i_endidx, rl_start, rl_end )

    IF (jb == jb1)  i_startidx = jc1
    IF (jb == jb2)  i_endidx   = jc2
    ncstart(2) = jcol1(jb)
    nccount(2) = nidx (jb)
    ncimap (1) = nidx (jb)

    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_lon ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%lon       (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_lat ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%lat       (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_z   ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%z         (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_dlon), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dlon      (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_dlat), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dlat      (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_dz  ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dz        (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_k   ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%k         (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_l   ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%l         (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_m   ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%m         (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_dk  ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dk        (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_dl  ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dl        (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_dm  ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%dm        (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_double(ncid, ncvarid(iv_wad ), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%wadens    (i_startidx:i_endidx,:,jb)) )

    CALL nf( nf_get_varm_int(ncid, ncvarid(iv_iex), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%iexist    (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_int(ncid, ncvarid(iv_sid), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%specid    (i_startidx:i_endidx,:,jb)) )
    CALL nf( nf_get_varm_int(ncid, ncvarid(iv_kac), ncstart, nccount, ncstri, ncimap,  &
      &                  p_ray(jg)%jk_active (i_startidx:i_endidx,:,jb)) )

    p_ray(jg)%coslat(i_startidx:i_endidx,:,jb) = MAX(min_coslat,  &
      &  ABS(COS(p_ray(jg)%lat(i_startidx:i_endidx,:,jb))))

  ENDDO

  IF ( p_patch%geometry_info%geometry_type == planar_torus_geometry )  &
    &  p_ray(jg)%coslat(:,:,:) = 1._wp   ! override

  ! [2]  (spectral) source variables

  ncstart(1) = 1   ;   nccount(1) = nrays_add(jg)
  DO jb = jb1, jb2
!   IF (nidx(jb) == 0)  CYCLE
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
      &                 i_startidx, i_endidx, rl_start, rl_end )
    IF (jb == jb1)  i_startidx = jc1
    IF (jb == jb2)  i_endidx   = jc2
    ncstart(2) = jcol1(jb)
    nccount(2) = nidx (jb)
    ncimap (1) = nidx (jb)
    CALL nf( nf_get_varm_int(ncid, ncvarid(iv_jrl), ncstart, nccount, ncstri, ncimap,  &
      &                      p_ray(jg)%jr_last(i_startidx:i_endidx,:,jb)) )
  ENDDO

END SUBROUTINE read_restartfile_specified
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE msgwam_read_restartfiles( mtime_current, p_patch )

  TYPE(datetime), POINTER            :: mtime_current   !< current datetime (mtime)
  TYPE(t_patch) , TARGET, INTENT(IN) :: p_patch

  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=1024)   :: r_filename, r_filename_sngl_conc
  CHARACTER(LEN=1024)   :: r_filename_head

  LOGICAL , ALLOCATABLE :: l_target(:,:)

! INTEGER , ALLOCATABLE :: gidx(:)
  REAL(wp), ALLOCATABLE :: glat(:), glon(:)

  INTEGER , ALLOCATABLE :: jfile_search(:)

  INTEGER  :: rl_start, rl_end
  INTEGER  :: i_startblk, i_endblk  ! blocks
  INTEGER  :: i_startidx, i_endidx  ! slices
  INTEGER  :: jg
  INTEGER  :: my_work_id, my_work_size
  INTEGER  :: ncol_in               ! number of grid columns in an input file
  INTEGER  :: jb_start
  INTEGER  :: nfile_total
  LOGICAL  :: ex, l_sngl_conc, l_dimcheck, l_read_data, l_first_read

  ! netCDF-related variables
  INTEGER  :: ncid, ncvarid(99), ncdummy, ncdimid, ncdimlen

  INTEGER  :: jc,jb, jf, jcol_in
  INTEGER  :: jcol_1st, jb1, jc1, jb2, jc2

  REAL(wp), PARAMETER ::  ll_tol = 1.e-5_wp   ! [rad]

  jg = p_patch%id

  IF ( .NOT. my_process_is_work() )  RETURN
  IF ( lrestart_from_nothing .OR. lsteady )  RETURN

  CALL message('', ' ')
  CALL message('msgwam_read_restartfiles',  &
    &          'MS-GWaM: ray-restart file reading')

  my_work_id   = get_my_mpi_work_id()
  my_work_size = get_my_mpi_work_comm_size()

  rl_start = grf_bdywidth_c+1
  rl_end   = min_rlcell_int - 2

  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

  ! filename head of the set of restart files
  WRITE(r_filename_head,'(a,i0,5(a,i2.2),a,i2.2)')                        &
    &    'msgwam_ray_restart_', mtime_current%date%year,                  &
    &    '-', mtime_current%date%month,  '-', mtime_current%date%day,     &
    &    '_', mtime_current%time%hour,   ':', mtime_current%time%minute,  &
    &    ':', mtime_current%time%second, '_D', jg
  IF (TRIM(dir_restartfiles) /= '')  &
    &  r_filename_head = TRIM(dir_restartfiles)//'/'//TRIM(r_filename_head)

  ! get the number of files in the set
  nfile_total = 0  ;  l_sngl_conc = .FALSE.
  DO
    WRITE(r_filename,'(a,i6.6,a)') TRIM(r_filename_head)//'_P', nfile_total, '.nc'
    INQUIRE(FILE=TRIM(r_filename), EXIST=ex)
    IF ( .NOT. ex )  EXIT
    nfile_total = nfile_total + 1
  ENDDO
  ! in case none of them exist, check if the single, concatenated file exists
  IF (nfile_total == 0) then
    WRITE(r_filename_sngl_conc,'(a)') TRIM(r_filename_head)//'.nc'
    INQUIRE(FILE=TRIM(r_filename_sngl_conc), EXIST=l_sngl_conc)
    IF ( .NOT. l_sngl_conc )  CALL finish('msgwam_read_restartfiles',  &
      &                                    'MS-GWaM: restart file not found')
    nfile_total = 1
  END IF

  ! make the order of files to search for the very first grid (jc1,jb1)
  IF ( .NOT. l_sngl_conc ) THEN
    ALLOCATE( jfile_search(nfile_total) )

    ! first guess of the target file
    jfile_search(1) = my_work_id*nfile_total/my_work_size

    ! searching order : jfile_search(1) + 0, +1, -1, +2, -2, ...
    jfile_search(2) = jfile_search(1) + 1
    IF (jfile_search(2) == nfile_total)  jfile_search(2) = 0
    DO jf = 3, nfile_total, 2
      jfile_search(jf) = jfile_search(jf-2) - 1
      IF (jfile_search(jf) == -1)  jfile_search(jf) = nfile_total-1
    ENDDO
    DO jf = 4, nfile_total, 2
      jfile_search(jf) = jfile_search(jf-2) + 1
      IF (jfile_search(jf) == nfile_total)  jfile_search(jf) = 0
    ENDDO
  END IF

  !---  Open the file  &  Check idx/lat/lon  -----------------------------------

  ! tar for grids to get data
  ALLOCATE( l_target(nproma,i_startblk:i_endblk) )
  l_target(:,:) = .FALSE.
  DO jb = i_startblk, i_endblk
    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
      &                 i_startidx, i_endidx, rl_start, rl_end )
    l_target(i_startidx:i_endidx,jb) = .TRUE.
  ENDDO

  ! the 1st grid
  jb_start = i_startblk

  LOOP_FILES:  DO jf = 1, nfile_total

    ! filename
    IF ( l_sngl_conc ) THEN
      r_filename = TRIM(r_filename_sngl_conc)
    ELSE
      WRITE(r_filename,'(a,i6.6,a)') TRIM(r_filename_head)//'_P', jfile_search(jf), '.nc'
    END IF

    ! open the file
    CALL nf( nf_open(TRIM(r_filename), NF_NOWRITE, ncid) )

    ! dimension length (number of ray volumes) check -- do only once
    IF (jf == 1) THEN
      CALL nf( nf_inq_dimid(ncid, 'ray', ncdimid) )
      CALL nf( nf_inq_dimlen(ncid, ncdimid, ncdimlen) )
      l_dimcheck = ( ncdimlen == nrays(jg) )
      CALL nf( nf_inq_dimid(ncid, 'spec', ncdimid) )
      CALL nf( nf_inq_dimlen(ncid, ncdimid, ncdimlen) )
      l_dimcheck = l_dimcheck .AND. ( ncdimlen == nrays_add(jg) )
      IF ( .NOT. l_dimcheck )  CALL finish('msgwam_read_restartfiles',  &
        &                                  'MS-GWaM: inconsistent dimension')
    END IF

    ! number of grid columns in this file
    CALL nf( nf_inq_dimid(ncid, 'col', ncdimid) )
    CALL nf( nf_inq_dimlen(ncid, ncdimid, ncol_in) )

    ! grid info
!   ALLOCATE( gidx(ncol_in) )
!   CALL nf( nf_inq_varid(ncid, 'gidx', ncdummy) )
!   CALL nf( nf_get_var_int(ncid, ncdummy, gidx) )
    ALLOCATE( glat(ncol_in), glon(ncol_in) )
    CALL nf( nf_inq_varid(ncid, 'glat', ncdummy) )
    CALL nf( nf_get_var_double(ncid, ncdummy, glat) )
    CALL nf( nf_inq_varid(ncid, 'glon', ncdummy) )
    CALL nf( nf_get_var_double(ncid, ncdummy, glon) )

    l_first_read = .TRUE.   ! per file

    ! start from the 1st written grid in the input file
    jcol_in = 1
    DO WHILE ( jcol_in <= ncol_in )

      l_read_data = .FALSE.

      J_BLK:  DO jb = jb_start, i_endblk   ! jb_start may change
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
          &                 i_startidx, i_endidx, rl_start, rl_end )
        DO jc = i_startidx, i_endidx
          IF ( l_target(jc,jb) ) THEN   ! this and the following IF conditions are split
                                        !   to minimize the number of calculations
            IF (   ABS(p_patch%cells%center(jc,jb)%lat - glat(jcol_in)) < ll_tol .AND.  &
                &  ABS(p_patch%cells%center(jc,jb)%lon - glon(jcol_in)) < ll_tol ) THEN
              ! grid matched
              IF ( .NOT. l_read_data ) THEN   ! the 1st grid to read, for this iteration
                l_read_data = .TRUE.
                jcol_1st = jcol_in
                jb1 = jb  ;  jc1 = jc
              END IF
              jb2 = jb  ;  jc2 = jc
              l_target(jc,jb) = .FALSE.   ! no more a target, from the next iteration
              jcol_in = jcol_in + 1
              IF ( jcol_in == ncol_in + 1 )  EXIT J_BLK
            ELSE   ! not matched with the 'jcol_in' grid in the input file
              IF ( l_read_data )  EXIT J_BLK  ! the successive grid matching is stopped here
            END IF
          ELSE   ! not a target grid
            IF ( l_read_data )  EXIT J_BLK  ! the successive grid matching is stopped here
          END IF
        ENDDO
      ENDDO J_BLK

      IF ( l_read_data ) THEN

        ! sequential data reading for the given range of successive grids
        CALL read_restartfile_specified( ncid, ncvarid, jcol_1st, p_patch,  &
          &     jb1, jc1, jb2, jc2, l_first_read )
        l_first_read = .FALSE.

        ! shorten the index range for grid searching in order to reduce iteration
        DO jb = jb_start, i_endblk
          IF ( ANY(l_target(:,jb)) )  EXIT
          IF (jb == i_endblk) THEN   ! DONE :  ALL the grids have got data
            DEALLOCATE( glat, glon )
!           DEALLOCATE( gidx )
            CALL nf( nf_close(ncid) )
            EXIT LOOP_FILES
          END IF
          jb_start = jb + 1
        ENDDO

      ELSE   ! NO grids are matched with jcol_in

        jcol_in = jcol_in + 1

      END IF

    ENDDO  ! jcol_in <= ncol_in

    DEALLOCATE( glat, glon )
!   DEALLOCATE( gidx )
    CALL nf( nf_close(ncid) )

  ENDDO LOOP_FILES

  !---  Close the file and finalize  -------------------------------------------

  DEALLOCATE( l_target )
  IF ( .NOT. l_sngl_conc )  DEALLOCATE( jfile_search )

  IF ( ltest_restart(1) )  ltest_restart(2) = .TRUE.   ! flag to test restart-file reading

  CALL message('msgwam_read_restartfiles', 'MS-GWaM: finished reading')
  CALL message('', ' ')

END SUBROUTINE msgwam_read_restartfiles
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE nf( i_status )
  INTEGER, INTENT(in) :: i_status
  INCLUDE 'netcdf.inc'
  IF (i_status /= NF_NOERR)  CALL finish(  &
    &   'msgwam_read/write_restartfiles: netCDF error', nf_strerror(i_status) )
END SUBROUTINE nf
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_setup_msgwam_interface
