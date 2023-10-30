!>
!! configuration setup for gravity wave source schemes
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author Young-Ha Kim, Goethe Uni Frankfurt
!!
!!
!! @par Revision History
!! Initial writing by Young-Ha Kim, Goethe Uni Frankfurt (2018-06-11)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_gw_source_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: max_dom

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  t_gws_conv_config, gws_conv_config, max_nscale_cgw

  !--------------------------------------------------------------------------
  ! Configuration setup for convective gravity wave drag scheme
  !--------------------------------------------------------------------------

  INTEGER, PARAMETER ::  max_nscale_cgw = 3

  TYPE t_gws_conv_config

    INTEGER  ::  n_source(max_dom)

    INTEGER  ::  conv_type
    REAL(wp) ::  dc       (2)
    REAL(wp) ::  c_dc_vari(2)
    REAL(wp) ::  lzmax_scale
    REAL(wp) ::  r_maxmin_m
    INTEGER  ::  nm           (max_nscale_cgw)
    INTEGER  ::  nphi         (max_nscale_cgw)
    REAL(wp) ::  scale_h      (max_nscale_cgw)
    REAL(wp) ::  scale_t      (max_nscale_cgw)
    REAL(wp) ::  heat_areafrac(max_nscale_cgw)
    REAL(wp) ::  heat_min
    REAL(wp) ::  crit_zctop
    REAL(wp) ::  dz_launch
    INTEGER  ::  limfactor_conv
!   REAL(wp) ::  flx_modul(7)
    INTEGER  ::  diag_code

  END TYPE t_gws_conv_config

  TYPE(t_gws_conv_config) ::  gws_conv_config

END MODULE mo_gw_source_config
