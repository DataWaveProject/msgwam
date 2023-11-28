!>
!! utilities module for for MS-GWaM
!!
!! @author Georg Sebastian Voelker, Goethe Uni Frankfurt
!!
!!
!! @par Revision History
!! Initial revision by Georg Sebastian Voelker, Goethe Uni Frankfurt (2023-11-28)
!!
!! @par Copyright and License
!!
!! This code is subject to the Creative Commons Licnese, version 4.0 with Attribution
!! extension (CC-BY 4.0). It is distributed as is without any warranty. 
!!
MODULE mo_msgwam_util

    USE mo_kind,                    ONLY: wp, vp, sp
    USE mo_model_domain,            ONLY: t_patch
    USE mo_intp_data_strc,          ONLY: t_int_state
    USE mo_setup_msgwam_interface,  ONLY: t_gridinfo4ray
    USE mo_parallel_config,         ONLY: nproma
    USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c
    USE mo_impl_constants,          ONLY: min_rlcell_int
    USE mo_loopindices,             ONLY: get_indices_c

    PUBLIC  ::  smooth_vert
    PUBLIC  ::  smooth_hori

CONTAINS

SUBROUTINE smooth_vert(nlev,i_startidx,i_endidx,npts,lchange_bdy, var)
  INTEGER,  INTENT(IN)    :: nlev
  INTEGER,  INTENT(IN)    :: i_startidx ! first index of the block
  INTEGER,  INTENT(IN)    :: i_endidx   ! last index of the block
  INTEGER,  INTENT(IN)    :: npts       ! half number of points to smooth
  LOGICAL,  INTENT(IN)    :: lchange_bdy  ! flag for changing boundary values
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
SUBROUTINE smooth_hori( p_patch,        & !inout
  &                     p_gridinfo4ray, & !in
  &                     p_int_state,    & !in
  &                     nlev,           & !in
  &                     var             ) !inout

  ! In/out variables
  TYPE(t_patch),        INTENT(INOUT) :: p_patch
  TYPE(t_gridinfo4ray), INTENT(IN)    :: p_gridinfo4ray
  TYPE(t_int_state),    INTENT(IN)    :: p_int_state
  INTEGER,              INTENT(IN)    :: nlev
  REAL(wp),             INTENT(INOUT) :: var(:,:,:)

  INTEGER                 :: rl_start, rl_end
  INTEGER                 :: i_startblk, i_endblk
  INTEGER                 :: i_startidx, i_endidx
  INTEGER                 :: jb, jc, jk, jv, jc_v
  INTEGER                 :: vidx(3), vblk(3)
  INTEGER                 :: cidx_v(6,3), cblk_v(6,3)
  INTEGER                 :: nbcells_v
  REAL(wp)                :: var_tmp(nproma,nlev,p_patch%nblks_c)
  REAL(wp)                :: avg_v(nlev,3)

  !----------------------------------------------------------------------
  ! Purpose:
  !         
  !
  ! Method:
  !         -- 
  !
  !----------------------------------------------------------------------

  ! Synchronize var to smooth
!  CALL sync_patch_array(SYNC_C,p_patch,var)  ! placed outside to apply for multiple variables

  IF (.TRUE.) THEN
  ! All neighbor smoothing with weights suggested by Young-Ha:
  ! 1) calculate non-weighted average of cells corresponding to 
  !    each vertex of the actual cell
  ! 2) calculate non-weighted average of these 3 average values
  !
  ! Equivalent to using weightings of
  !   1/6 for the target (center) cell,
  !   1/9 for 3 direct neighbors, and
  !   1/18 for 9 non-direct neighbors.

  ! Prognostic domain
  rl_start = grf_bdywidth_c + 1
  rl_end   = min_rlcell_int
  i_startblk = p_patch%cells%start_block(rl_start)
  i_endblk   = p_patch%cells%end_block(rl_end)

  ! TODO : How about using a loop for vertices (jv) directly.
  !        It will do less calculation (~1/6) than now.

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jv, jc_v, jk, i_startidx, i_endidx, nbcells_v,  &
!$OMP            avg_v, vidx, vblk, cidx_v, cblk_v) ICON_OMP_GUIDED_SCHEDULE

  ! Loop over cells
  DO jb = i_startblk, i_endblk
    var_tmp(:,:,jb) = 0._wp

    CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

    DO jc = i_startidx, i_endidx

      avg_v(:,:) = 0._wp

      ! Get line and block indices of cell vertices
      vidx(1:3) = p_patch%cells%vertex_idx(jc,jb,1:3)
      vblk(1:3) = p_patch%cells%vertex_blk(jc,jb,1:3)

      ! For each vertex: get all the cells which share this vertex
      DO jv = 1, 3
        nbcells_v = p_patch%verts%num_edges(vidx(jv),vblk(jv))
        cidx_v(1:nbcells_v,jv) = p_patch%verts%cell_idx(vidx(jv),vblk(jv),1:nbcells_v)
        cblk_v(1:nbcells_v,jv) = p_patch%verts%cell_blk(vidx(jv),vblk(jv),1:nbcells_v)

        DO jc_v = 1, nbcells_v
          DO jk = 1, nlev
            avg_v(jk,jv) = avg_v(jk,jv) + var(cidx_v(jc_v,jv),jk,cblk_v(jc_v,jv))
          ENDDO ! jk
        ENDDO ! jc_v
        avg_v(:,jv) = avg_v(:,jv)/REAL(nbcells_v, KIND=wp)

      ENDDO ! jv

      DO jk = 1, nlev
        var_tmp(jc,jk,jb) = SUM(avg_v(jk,:))/3._wp
      ENDDO ! jk

    ENDDO ! jc

  ENDDO ! jb

!$OMP END DO
!$OMP END PARALLEL

  var(:,:,i_startblk:i_endblk) = var_tmp(:,:,i_startblk:i_endblk)

  ENDIF

END SUBROUTINE smooth_hori

END MODULE