!>
!! Namelist for the gravity wave source schemes in MS-GWaM
!!
!! @author Young-Ha Kim, Goethe Uni Frankfurt
!!
!! @par Revision History
!! Initial revision by Young-Ha Kim, Goethe Uni Frankfurt (2018-06-11)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_gw_source_nml

  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_exception,           ONLY: finish, message, message_text

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom, inwp

  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_gw_source_config,    ONLY: gws_conv_config, max_nscale_cgw
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  USE mo_run_config,          ONLY: iforcing
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: read_gw_source_namelist

  !-----------------------------------!
  ! gw_source_nml namelist variables  !
  !-----------------------------------!

  INTEGER  ::  cgw_n_source(max_dom)
  INTEGER  ::  cgw_conv_type
  REAL(wp) ::  cgw_dc       (2)
  REAL(wp) ::  cgw_c_dc_vari(2)
  REAL(wp) ::  cgw_lzmax_scale
  REAL(wp) ::  cgw_r_maxmin_m
  INTEGER  ::  cgw_nm           (max_nscale_cgw)
  INTEGER  ::  cgw_nphi         (max_nscale_cgw)
  REAL(wp) ::  cgw_scale_h      (max_nscale_cgw)
  REAL(wp) ::  cgw_scale_t      (max_nscale_cgw)
  REAL(wp) ::  cgw_heat_areafrac(max_nscale_cgw)
  REAL(wp) ::  cgw_heat_min
  REAL(wp) ::  cgw_crit_zctop
  REAL(wp) ::  cgw_dz_launch
  INTEGER  ::  cgw_limfactor
! REAL(wp) ::  cgw_flx_modul(7)

  NAMELIST /gw_source_nml/ cgw_n_source, cgw_conv_type, cgw_dc, cgw_c_dc_vari,  &
    &      cgw_lzmax_scale, cgw_r_maxmin_m,                                     &
    &      cgw_nm, cgw_nphi, cgw_scale_h, cgw_scale_t, cgw_heat_areafrac,       &
    &      cgw_heat_min, cgw_crit_zctop, cgw_dz_launch, cgw_limfactor  !,          &
!   &      cgw_flx_modul

CONTAINS


  SUBROUTINE read_gw_source_namelist( filename )

    CHARACTER(len=*), INTENT(IN) ::  filename

    CHARACTER(len=*), PARAMETER ::  routine  &
      &  = 'mo_gw_source_nml: read_gw_source_namelist'

    INTEGER ::  istat, funit
    INTEGER ::  jg 
    INTEGER ::  iunit
    INTEGER ::  js

    !------------------------------------------------------------------
    ! 1. Set default values
    !------------------------------------------------------------------

    cgw_n_source(:)        = 0
    cgw_conv_type          = 1
    cgw_dc       (1:2)     = (/  2.0_wp ,  4.0_wp /)
    cgw_c_dc_vari(1:2)     = (/ 15.0_wp , 40.0_wp /)
    cgw_lzmax_scale        = 1.0_wp
    cgw_r_maxmin_m         = 128.0_wp
    cgw_nm  (1:3)          = (/ 7*4 + 1   , 7*4 + 1          , 7*4 + 1           /)
    cgw_nphi(1:3)          = (/ 4         , 4                , 4                 /)
    cgw_scale_h(1:3)       = (/ 5.e3_wp   , 100.e3_wp        , 1000.e3_wp        /)
    cgw_scale_t(1:3)       = (/ 1200.0_wp , 6.0_wp*3600.0_wp , 48.0_wp*3600.0_wp /)
    cgw_heat_areafrac(1:3) = (/ 0.01      , 0.01             , 0.01              /)
    cgw_heat_min           = 1.e-6_wp/86400.0_wp
    cgw_crit_zctop         = 3.e3_wp
    cgw_dz_launch          = 1.e3_wp
    cgw_limfactor          = 0
!   cgw_flx_modul(1:7)     = (/ 1.0_wp , 15.0_wp , 10.0_wp ,  &
!     &                                  15.e3_wp , 2.e3_wp , 35.e3_wp , 3.e3_wp /)

    !-------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !-------------------------------------------------------------------

    IF ( use_restart_namelists() ) THEN
      funit = open_and_restore_namelist('gw_source_nml')
      READ(funit, NML=gw_source_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------
    ! 3. Read user's (new) specifications (done by all MPI processes)
    !-------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml('gw_source_nml', status=istat)
    IF ( my_process_is_stdio() ) THEN
      iunit = temp_defaults()
      WRITE(iunit, gw_source_nml)        ! write defaults to temporary text file
    END IF
    SELECT CASE ( istat )
    CASE ( positioned )
      READ (nnml, gw_source_nml)         ! overwrite default settings
      IF ( my_process_is_stdio() ) THEN
        iunit = temp_settings()
        WRITE(iunit, gw_source_nml)      ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !-------------------------------------------------------------------
    ! 4. Sanity Check
    !-------------------------------------------------------------------

    ! convective GW calculation should be turned off when convection is off
    IF (iforcing == inwp) THEN
      DO jg = 1, max_dom
        IF ( atm_phy_nwp_config(jg)%inwp_convection == 0 .AND. cgw_n_source(jg) /= 0 ) THEN
          cgw_n_source(jg) = 0
          CALL message(TRIM(routine), 'Turning off convective GW scheme !')
        END IF
      ENDDO
    END IF

    IF (ANY( cgw_n_source(:) > 0 )) THEN
      IF ( cgw_conv_type <= 0 .OR. cgw_conv_type > 4 )  &
        &  CALL finish(TRIM(routine), 'invalid cgw_conv_type (1-4)')
      IF (ANY( MOD(cgw_nphi(:),2) /= 0 ))  &
        &  CALL finish(TRIM(routine), 'nphi should be a multiple of 2')
    END IF

    !-------------------------------------------------------------------
    ! 5. Store the namelist for restart
    !-------------------------------------------------------------------

    IF ( my_process_is_stdio() ) THEN
      funit = open_tmpfile()
      WRITE(funit, NML=gw_source_nml)                    
      CALL store_and_close_namelist(funit, 'gw_source_nml') 
    ENDIF

    !-------------------------------------------------------------------
    ! 6. Write the namelist to an ASCII file
    !-------------------------------------------------------------------

    IF ( my_process_is_stdio() )  WRITE(nnml_output, nml=gw_source_nml)

    !-------------------------------------------------------------------
    ! 7. Fill the configuration state
    !-------------------------------------------------------------------

    DO jg = 1, max_dom
      gws_conv_config% n_source(jg) = cgw_n_source(jg)
    ENDDO
    gws_conv_config% conv_type        = cgw_conv_type
    gws_conv_config% dc           (:) = cgw_dc           (:)
    gws_conv_config% c_dc_vari    (:) = cgw_c_dc_vari    (:)
    gws_conv_config% lzmax_scale      = cgw_lzmax_scale
    gws_conv_config% r_maxmin_m       = cgw_r_maxmin_m
    gws_conv_config% nm           (:) = cgw_nm           (:)
    gws_conv_config% nphi         (:) = cgw_nphi         (:)
    gws_conv_config% scale_h      (:) = cgw_scale_h      (:)
    gws_conv_config% scale_t      (:) = cgw_scale_t      (:)
    gws_conv_config% heat_areafrac(:) = cgw_heat_areafrac(:)
    gws_conv_config% heat_min         = cgw_heat_min
    gws_conv_config% crit_zctop       = cgw_crit_zctop
    gws_conv_config% dz_launch        = cgw_dz_launch
    gws_conv_config% limfactor_conv   = cgw_limfactor
!   gws_conv_config% flx_modul    (:) = cgw_flx_modul    (:)

    IF (ANY( gws_conv_config%n_source(:) > 0 )) THEN
      WRITE(message_text,'(a)') 'parameter values :'
      CALL message(TRIM(routine),message_text)
      WRITE(message_text,'(a)') '|  nm  | nphi | heat_areafrac |'
      CALL message(' ',message_text)
      DO js = 1, MAXVAL(cgw_n_source(:))
        WRITE(message_text,'(a,i3,a,i3,a,f13.6,a)')  &
          &  '| ', cgw_nm(js), '  | ', cgw_nphi(js), '  | ', cgw_heat_areafrac(js), ' |'
        CALL message(' ',message_text)
      ENDDO
    END IF

  END SUBROUTINE read_gw_source_namelist

END MODULE mo_gw_source_nml
