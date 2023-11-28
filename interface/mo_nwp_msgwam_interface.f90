!>
!! This module is an additional interface between nwp_gw_interface to the
!! gravity wave drag parametrization developed at the Goethe Uni:
!! Lagrangian WKB raytracer in phase space.
!!
!! @author Gergely Bölöni, Goethe Uni Frankfurt
!!
!! @par Revision History
!! Initial version by Gergely Bölöni, Goethe Uni Frankfurt (2016-03-02)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_nwp_msgwam_interface

  USE mo_kind,                   ONLY: wp
  USE mo_model_domain,           ONLY: t_patch
  USE mo_nonhydro_types,         ONLY: t_nh_diag, t_nh_prog, t_nh_metrics
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_nwp_phy_types,          ONLY: t_nwp_phy_tend
  USE mtime,                     ONLY: datetime
  USE mo_impl_constants,         ONLY: min_rlcell_int
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_exception,              ONLY: message
!#ifdef msgwam
  USE mo_msgwam_config,          ONLY: lsteady, nrays, lmsgwam_offline, &
                                       lmsgwam_noforce, ltest_restart, lmsgwam_pmomflux
  USE mo_msgwam,                 ONLY: gwdrag_msgwam
  USE mo_msgwam_stst,            ONLY: gwdrag_msgwam_stst
  USE mo_setup_msgwam_interface, ONLY: p_msgwam, msgwam_write_restartfiles
!#endif
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nwp_msgwam_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE nwp_msgwam_interface( dt_call,                   & !>input
      &                          mtime_datetime,            & !>input
      &                          p_sim_time,                & !>input
      &                          p_patch, p_metrics,        & !>input
      &                          p_int_state,               & !>input
      &                          p_diag,p_prog,             & !>input
      &                          prm_nwp_tend               ) !>inout
  TYPE(datetime),      POINTER,INTENT(IN)   :: mtime_datetime       !< date/time information
  TYPE(t_patch),        TARGET,INTENT(INOUT):: p_patch              ! grid/patch info.
  TYPE(t_nh_metrics)          ,INTENT(IN)   :: p_metrics
  TYPE(t_int_state),    TARGET,INTENT(IN)   :: p_int_state
  TYPE(t_nh_prog),      TARGET,INTENT(IN)   :: p_prog               ! the dyn prog vars
  TYPE(t_nh_diag),      TARGET,INTENT(IN)   :: p_diag               ! the dyn diag vars
  TYPE(t_nwp_phy_tend), TARGET,INTENT(INOUT):: prm_nwp_tend         ! atm tend vars
  REAL(wp),                    INTENT(IN)   :: dt_call              ! time step
  REAL(wp),                    INTENT(IN)   :: p_sim_time           ! elapsed simulation time

  INTEGER  :: jg
  INTEGER  :: rl_start, rl_end
  INTEGER  :: i_startblk, i_endblk ! blocks
  INTEGER  :: i_startidx, i_endidx ! slices
  INTEGER  :: jb                   ! block indeces


! Later use this ifdef
!#ifdef msgwam

  IF (nrays(p_patch%id) == 0)  RETURN

  ! Domain ID
  jg = p_patch%id

  ! Decide whether we run MS-GWaM or its steady state version
  IF (.NOT. lsteady) THEN

    IF ( ltest_restart(2) )  &   ! flag to test restart-file reading
      &  CALL msgwam_write_restartfiles( mtime_datetime, p_patch )

    CALL gwdrag_msgwam ( dt_call                    , & !>input
      &                  mtime_datetime, p_sim_time , & !>input
      &                  p_patch, p_metrics         , & !>input
      &                  p_int_state                , & !>input
      &                  p_prog% rho                , & !>input
      &                  p_diag% pres               , & !>input
      &                  p_diag% pres_sfc           , & !>input
      &                  p_diag% temp               , & !>input
      &                  p_diag% temp_ifc           , & !>input
      &                  p_prog% theta_v            , & !>input
      &                  p_diag% theta_v_ic         , & !>input
      &                  p_diag% u                  , & !>input
      &                  p_diag% v                  , & !>input
      &                  p_msgwam(jg)               )   !>inout

  ELSE

    CALL gwdrag_msgwam_stst ( dt_call                    , & !>input
      &                       mtime_datetime, p_sim_time , & !>input
      &                       p_patch, p_metrics         , & !>input
      &                       p_prog% rho                , & !>input
      &                       p_diag% u                  , & !>input
      &                       p_diag% v                  , & !>input
      &                       p_diag% temp               , & !>input
      &                       p_diag% temp_ifc           , & !>input
      &                       p_msgwam(jg)               )   !>inout

  ENDIF

  ! Transfer the resultant GWD to the NWP physics tendency array prm_nwp_tend
  IF (lmsgwam_offline .OR. lmsgwam_noforce) THEN
    CALL message('nwp_msgwam_interface', 'MS-GWaM in offline / noforce mode!')
  ELSE
    ! Exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    IF (lmsgwam_pmomflux) THEN
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
          &                 i_startidx, i_endidx, rl_start, rl_end )

        prm_nwp_tend%        ddt_u_gwd         (i_startidx:i_endidx,:,jb)  &
          &  = p_msgwam(jg)% ddt_u_gwd_pmom_mgm(i_startidx:i_endidx,:,jb)
        prm_nwp_tend%        ddt_v_gwd         (i_startidx:i_endidx,:,jb)  &
          &  = p_msgwam(jg)% ddt_v_gwd_pmom_mgm(i_startidx:i_endidx,:,jb)
        prm_nwp_tend%        ddt_temp_drag(i_startidx:i_endidx,:,jb) = 0
      ENDDO
    ELSE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
          &                 i_startidx, i_endidx, rl_start, rl_end )

        prm_nwp_tend%        ddt_u_gwd    (i_startidx:i_endidx,:,jb)  &
          &  = p_msgwam(jg)% ddt_u_gwd_mgm(i_startidx:i_endidx,:,jb)
        prm_nwp_tend%        ddt_v_gwd    (i_startidx:i_endidx,:,jb)  &
          &  = p_msgwam(jg)% ddt_v_gwd_mgm(i_startidx:i_endidx,:,jb)
        prm_nwp_tend%        ddt_temp_drag(i_startidx:i_endidx,:,jb)  &
          &  = p_msgwam(jg)% ddt_t_gwd_mgm(i_startidx:i_endidx,:,jb)
      ENDDO
    ENDIF
  END IF

!#endif

END SUBROUTINE nwp_msgwam_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_nwp_msgwam_interface
