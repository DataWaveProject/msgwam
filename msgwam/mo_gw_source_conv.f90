!>
!! gw_source_conv: parameters and initialization for gravity wave source schemes
!!                 in MS-GWaM
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

MODULE mo_gw_source_conv

  USE mo_kind                   , ONLY: wp
  USE mo_exception              , ONLY: finish, print_value
  USE mo_parallel_config        , ONLY: nproma
  USE mo_model_domain           , ONLY: t_patch
  USE mo_impl_constants         , ONLY: min_rlcell_int
  USE mo_impl_constants_grf     , ONLY: grf_bdywidth_c
  USE mo_loopindices            , ONLY: get_indices_c
  USE mo_run_config             , ONLY: msg_level
  USE mo_timer
  USE mo_name_list_output_config, ONLY: is_variable_in_output

  USE mo_math_constants         , ONLY: pi
  USE mo_physical_constants     , ONLY: grav, cpd, rd_o_cpd
  USE mo_gw_source_config       , ONLY: cfg => gws_conv_config,  &
    &                                   max_nscale_cgw

  USE mo_msgwam_config          , ONLY: nrays, lsteady
  USE mo_setup_msgwam_interface , ONLY: t_msgwam, p_msgwam, p_ray_conv

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  gws_conv_init, gw_source_conv

  INTEGER , PARAMETER ::  ncmax = 50
  REAL(wp), PARAMETER ::  zag_ll = 3000.0_wp
  REAL(wp), PARAMETER ::  zag_ml = 6000.0_wp
  REAL(wp), PARAMETER ::  ri_bg_limit = 1.0_wp/3.0_wp   ! > 1/4

  INTEGER , PARAMETER ::  dp = KIND(1.0d0)   ! needed regardless of wp

  REAL(wp), PARAMETER ::  rdz_ll = 1.0_wp/zag_ll
  REAL(wp), PARAMETER ::  rdz_ml = 1.0_wp/(zag_ml - zag_ll)

  REAL(wp), ALLOCATABLE ::  entrmax(:,:),  &
    &                       zm     (:,:),  &
    &                       zd     (:,:),  &
    &                       cqx    (:,:),  &
    &                       cqy    (:,:),  &
    &                       nbv_ct (:,:),  &
    &                       n2_q   (:,:),  &
    &                       rho0   (:,:),  &
    &                       u_ct   (:,:),  &
    &                       v_ct   (:,:),  &
    &                       u_cb   (:,:),  &
    &                       v_cb   (:,:)
  REAL(wp), ALLOCATABLE ::  zctop  (:,:)
  INTEGER , ALLOCATABLE ::  kcta   (:,:)
  INTEGER , ALLOCATABLE ::  flag_last(:,:)

  REAL(wp) ::  coef_gridmeanflux(max_nscale_cgw)
  REAL(wp) ::  entr_crit
  REAL(dp) ::  eipi_r, eipi_i
  REAL(wp) ::  c_intr_ct_s(ncmax)
  REAL(wp) ::  a_mwn_o_n_s(ncmax)
  REAL(wp) ::  dm_o_n_s   (ncmax)

  INTEGER  ::  j_initialized = -999


CONTAINS

  !>
  !! <Short description of the subroutine for listings and indices>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
  SUBROUTINE gws_conv_init ( z_ifc, p_patch,  &
    &                        ll_k, ml_k, ll_dz, ml_dz, test )

    REAL(wp)     , POINTER, INTENT(in   ) ::  z_ifc(:,:,:)
    TYPE(t_patch), TARGET , INTENT(in   ) ::  p_patch
    INTEGER               , INTENT(  out) ::  ll_k (:,:)
    INTEGER               , INTENT(  out) ::  ml_k (:,:)
    REAL(wp)              , INTENT(  out) ::  ll_dz(:,:)
    REAL(wp)              , INTENT(  out) ::  ml_dz(:,:)
    REAL(wp)              , INTENT(inout), OPTIONAL ::  test(:,:,:)

    REAL(wp), PARAMETER ::  c0_low = 1.0_wp

    REAL(dp), PARAMETER ::  pi_tr = 3.14159265358979323846_dp   ! truncated

    REAL(wp) ::  dc1, dc2, c1, c2
    REAL(wp) ::  c_intr_ct_s_low, dc, dc2mdc1_half, pi_o_c2mc1

    INTEGER  ::  nlev, nlevp1
    INTEGER  ::  rl_start, rl_end
    INTEGER  ::  i_startblk, i_endblk ! blocks
    INTEGER  ::  i_startidx, i_endidx ! slices

    INTEGER  ::  jc, jk, jb, ic

    ! define the constants

    eipi_r = COS(pi_tr)
    eipi_i = SIN(pi_tr)
    IF (eipi_i <= 0.0_dp)  CALL finish('mo_gw_source_conv:gws_conv_init',  &
      &                                'pi_tr should be slightly smaller than pi')

    cfg% crit_zctop = MAX(zag_ll, cfg% crit_zctop)

    ! target phase-speed components
    dc1 = cfg% dc(1)  ;  c1 = cfg% c_dc_vari(1)
    dc2 = cfg% dc(2)  ;  c2 = cfg% c_dc_vari(2)
    IF ( dc2 == dc1 .OR. dc2 <= 0.0_wp ) THEN  ! use of a constant dc
      dc2 = dc1
      c2 = c1 + 100.0_wp   ! meaningless, but letting it differ from c1
    END IF
    c_intr_ct_s_low = c0_low
    dc2mdc1_half = (dc2 - dc1)*0.5_wp
    pi_o_c2mc1 = pi/(c2 - c1)
    DO ic = 1, ncmax
      dc = dc1
      IF (c_intr_ct_s_low >= c2) THEN
        dc = dc2
      ELSE IF (c_intr_ct_s_low > c1) THEN
        ! dc = dc1 + (dc2 - dc1)/2 * {1 - cos[(c_low - c1)/(c2 - c1)*pi]}
        dc = dc1 + dc2mdc1_half*(1.0_wp - COS((c_intr_ct_s_low - c1)*pi_o_c2mc1))
      END IF
      c_intr_ct_s(ic) = c_intr_ct_s_low + 0.5_wp*dc
      a_mwn_o_n_s(ic) = 1.0_wp/c_intr_ct_s(ic)
      dm_o_n_s(ic) = dc*(a_mwn_o_n_s(ic)*a_mwn_o_n_s(ic))
      c_intr_ct_s_low = c_intr_ct_s_low + dc
    ENDDO

    coef_gridmeanflux(:) =  8.0_wp*cfg% heat_areafrac(:)       &
      &                    /(  pi*SQRT(pi)*cfg% scale_h(:)**2  &
      &                       *cfg% scale_t(:) )

    entr_crit = MAX(0.0_wp, cfg% heat_min/300.0_wp)

    ! define the domain-dependent (jg) values

    ll_k (:,:) = 0
    ml_k (:,:) = 0
    ll_dz(:,:) = 0.0_wp
    ml_dz(:,:) = 0.0_wp

    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1
!   jg     = p_patch%id       ! domain ID

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
        &                 i_startidx, i_endidx, rl_start, rl_end )

      DO jc = i_startidx, i_endidx
        DO jk = nlev-1, 1, -1
          IF ( z_ifc(jc,jk,jb) > zag_ll + z_ifc(jc,nlevp1,jb) ) THEN
            ll_k (jc,jb) = jk + 1
            ll_dz(jc,jb) = zag_ll + z_ifc(jc,nlevp1,jb) - z_ifc(jc,jk+1,jb)
            EXIT
          END IF
        ENDDO
      ENDDO
      DO jc = i_startidx, i_endidx
        DO jk = ll_k(jc,jb)-1, 1, -1
          IF ( z_ifc(jc,jk,jb) > zag_ml + z_ifc(jc,nlevp1,jb) ) THEN
            ml_k (jc,jb) = jk + 1
            ml_dz(jc,jb) = zag_ml + z_ifc(jc,nlevp1,jb) - z_ifc(jc,jk+1,jb)
            EXIT
          END IF
        ENDDO
      ENDDO

test(:,1,jb) = float(ll_k(:,jb))
test(:,2,jb) = ll_dz(:,jb)
    ENDDO

    ! define output code

    cfg% diag_code = 0
    IF (is_variable_in_output(var_name='ctmfl_specc_cgw'))  &
      &  cfg% diag_code = cfg% diag_code + 10
    IF (is_variable_in_output(var_name='ctmfl_spec0_cgw'))  &
      &  cfg% diag_code = cfg% diag_code + 100
    IF (is_variable_in_output(var_name='ctmfl_spec1_cgw'))  &
      &  cfg% diag_code = cfg% diag_code + 200

    j_initialized = p_patch%id

  END SUBROUTINE gws_conv_init

  !>
  !! <Short description of the subroutine for listings and indices>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
  SUBROUTINE alloc_local_var(nblks_c)

    INTEGER, INTENT(in) ::  nblks_c

    ALLOCATE( entrmax  (nproma,nblks_c),  &
      &       zm       (nproma,nblks_c),  &
      &       zd       (nproma,nblks_c),  &
      &       cqx      (nproma,nblks_c),  &
      &       cqy      (nproma,nblks_c),  &
      &       nbv_ct   (nproma,nblks_c),  &
      &       n2_q     (nproma,nblks_c),  &
      &       rho0     (nproma,nblks_c),  &
      &       u_ct     (nproma,nblks_c),  &
      &       v_ct     (nproma,nblks_c),  &
      &       u_cb     (nproma,nblks_c),  &
      &       v_cb     (nproma,nblks_c),  &
      &       zctop    (nproma,nblks_c),  &
      &       kcta     (nproma,nblks_c),  &
      &       flag_last(nproma,nblks_c) )

  END SUBROUTINE alloc_local_var

  SUBROUTINE dealloc_local_var()

    DEALLOCATE( entrmax, zm, zd, cqx, cqy, nbv_ct, n2_q, rho0,  &
      &         u_ct, v_ct, u_cb, v_cb, zctop, kcta, flag_last )

  END SUBROUTINE dealloc_local_var

  !>
  !! <Short description of the subroutine for listings and indices>
  !!
  !! <Describe the purpose of the subroutine and its algorithm(s).>
  !! <Include any applicable external references inline as module::procedure,>
  !! <external_procedure(), or by using @see.>
  !! <Don't forget references to literature.>
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
  SUBROUTINE gw_source_conv ( tcall_gws_conv_jg, js,        &
    &                         p_patch, z_ifc, z_mc, dz,     &
    &                         rho, u, v, temp_ifc, pres_ifc )

    TYPE(t_patch)  , TARGET, INTENT(in)    ::  p_patch            !< grid/patch info.
    REAL(wp)               , INTENT(in)    ::  tcall_gws_conv_jg  !< time interval for
                                                                  !< convective GWs
    INTEGER                , INTENT(in)    ::  js                 !< index for source
    REAL(wp)      , POINTER, INTENT(in)    ::  z_ifc(:,:,:)
    REAL(wp)      , POINTER, INTENT(in)    ::  z_mc (:,:,:)
    REAL(wp)      , POINTER, INTENT(in)    ::  dz   (:,:,:)   ! p_metrics%ddqz_z_full
    REAL(wp)      , POINTER, INTENT(in)    ::  rho  (:,:,:)
    REAL(wp)      , POINTER, INTENT(in)    ::  u    (:,:,:)
    REAL(wp)      , POINTER, INTENT(in)    ::  v    (:,:,:)
    REAL(wp)      , POINTER, INTENT(in)    ::  temp_ifc(:,:,:)
    REAL(wp)      , POINTER, INTENT(in)    ::  pres_ifc(:,:,:)

    REAL(wp) ::  mfs_ct(cfg%nm(js),cfg%nphi(js),nproma,p_patch%nblks_c),  &
      &          kh_ct (cfg%nm(js),cfg%nphi(js),nproma,p_patch%nblks_c),  &
      &          dkh_ct(cfg%nm(js),cfg%nphi(js),nproma,p_patch%nblks_c)
    REAL(wp) ::  a_mwn(cfg%nm(js),nproma,p_patch%nblks_c),  &
      &          dm   (cfg%nm(js),nproma,p_patch%nblks_c)

    REAL(wp) ::  cosphi(cfg%nphi(js),nproma,p_patch%nblks_c),  &
      &          sinphi(cfg%nphi(js),nproma,p_patch%nblks_c)

    INTEGER  ::  nlev, nlevp1            !< number of full and half levels
    INTEGER  ::  rl_start, rl_end
    INTEGER  ::  i_startblk, i_endblk    !< blocks
    INTEGER  ::  i_startidx, i_endidx    !< slices
    INTEGER  ::  i_nchdom                !< domain index

    REAL(wp) ::  entr  (p_patch%nlev)
    REAL(wp) ::  entrdz(p_patch%nlev)
    REAL(wp) ::  n2_ct, n2_ct_limit
    REAL(wp) ::  phi0, phi, dphi
    REAL(wp) ::  dmdphi(cfg%nm(js)), sqr_dkdl(cfg%nm(js)), jac(cfg%nm(js))

    INTEGER  ::  nphi05
    INTEGER  ::  kctop0, kcbas0, kentrmax, kll, kml
    INTEGER  ::  kctop2, kcbas2, kct, kcb, kcm
    REAL(wp) ::  tmp, tmps, sum05_entrdz

    INTEGER  ::  jray_offset
    INTEGER  ::  jc, jk, jb, jg, iphi, iray1, iray9, jss, jd

    INTEGER  ::  diag_code_10, diag_code_100

    !----------------------------------------------------------------------
    ! Purpose:
    !         Calculate forcing of GWs by convection for MS-GWaM
    ! Method:
    !         ???
    !----------------------------------------------------------------------

    CALL timer_start(timer_gws_test1)

    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1
    jg     = p_patch%id       ! domain ID

    ! memory allocation and initialization
    IF (js == 1) THEN

      CALL alloc_local_var(p_patch%nblks_c)

      entrmax(:,:) = 0.0_wp   ! mandatory

      zm     (:,:) = 0.0_wp
      zd     (:,:) = 0.0_wp
      cqx    (:,:) = 0.0_wp
      cqy    (:,:) = 0.0_wp
      nbv_ct (:,:) = 0.0_wp
      n2_q   (:,:) = 1.0_wp   ! denominator
      rho0   (:,:) = 0.0_wp
      u_ct   (:,:) = 0.0_wp
      v_ct   (:,:) = 0.0_wp
      u_cb   (:,:) = 0.0_wp
      v_cb   (:,:) = 0.0_wp
      kcta   (:,:) = 0

      zctop  (:,:) = 0.0_wp

      flag_last(:,:) = p_msgwam(jg)% flag_cgw(:,:)

      p_msgwam(jg)% flag_cgw(:,:) = -1

    END IF

    nphi05 = cfg% nphi(js)/2
    dphi   = pi/REAL(cfg% nphi(js)/2, kind=wp)

    ! setup
    IF (j_initialized /= jg)                            &
      &  CALL gws_conv_init( z_ifc, p_patch           , &
      &                      p_msgwam(jg)% ll_k_cgw   , &
      &                      p_msgwam(jg)% ml_k_cgw   , &
      &                      p_msgwam(jg)% ll_dz_cgw  , &
      &                      p_msgwam(jg)% ml_dz_cgw  , &
      &                      p_msgwam(jg)% test_2_cgw )

    cosphi(1,:,:) = 1.0_wp
    sinphi(1,:,:) = 0.0_wp
    DO iphi = 2, nphi05
      phi = dphi*REAL(iphi-1,kind=wp)
      cosphi(iphi,:,:) = COS(phi)
      sinphi(iphi,:,:) = SIN(phi)
    ENDDO

    p_msgwam(jg)% ctmfl_cgw(:,:,js*5+1:js*5+5) = 0.0_wp

    a_mwn(:,:,:) = -1.0_wp
    dm   (:,:,:) = 0.0_wp


    !:::::   LOOP BEGINS   :::::


    i_nchdom = MAX(1,p_patch%n_childdom)

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c + 1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1       )
    i_endblk   = p_patch%cells%end_blk  (rl_end  ,i_nchdom)


    IF (js == 1) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jk, i_startidx, i_endidx, entr, entrdz, n2_ct, n2_ct_limit,  &
!$OMP            kctop0, kcbas0, kentrmax, kctop2, kcbas2, kct, kcb, kcm,             &
!$OMP            tmp, tmps, sum05_entrdz)

    L_JB:  DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
        &                 i_startidx, i_endidx, rl_start, rl_end )

    L_JC:  DO jc = i_startidx, i_endidx

      kcbas0 = p_msgwam(jg)% mbas_con_cgw(jc,jb)-1
      kctop0 = p_msgwam(jg)% mtop_con_cgw(jc,jb)

      entr(:) = p_msgwam(jg)% entr_cgw(jc,:,jb)

      kentrmax = MAXLOC(entr(kctop0:kcbas0),1) + kctop0 - 1
      kentrmax = MAX(1, kentrmax)  ! for the below IF statement when ktype == 0

      IF ( p_msgwam(jg)% ktype_cgw(jc,jb) /= 0 .AND.  &
        &  entr(kentrmax) > entr_crit ) THEN

        !-----------------------------------------------------------------------
        ! Determine the center height, half depth, maximum of the entropy
        ! profile, and cloud-top vertical index.
        !-----------------------------------------------------------------------

        entrdz(:) = 0.0_wp

        entrdz(kctop0:kcbas0) = entr(kctop0:kcbas0)*dz(jc,kctop0:kcbas0,jb)
        DO jk = kctop0, kcbas0
          entrdz(jk) = MAX( 0.0_wp, entrdz(jk) )
        ENDDO

        sum05_entrdz = SUM(entrdz(kctop0:kcbas0))*0.5_wp

        kctop2 = kctop0   ;   kcbas2 = kcbas0
        tmp = 0.0_wp      ;   tmps = 0.0_wp
        DO jk = kctop0, kcbas0
          IF (entrdz(jk) /= 0.0_wp) THEN   !  > 0
            tmp  = tmp  + entrdz(jk)*z_mc(jc,jk,jb)
            tmps = tmps + entrdz(jk)
          ELSE IF (tmp == 0.0_wp) THEN       !  no positive values yet
            kctop2 = jk + 1
          ELSE                               !  a positive-valued profile found
            kcbas2 = jk - 1
            ! The profile should contain more than 50% of the entropy integrated
            !   over the whole layer or include the layer of entropy maximum
            !   (kentrmax). Repeat, if not. The latter condition, along with
            !   entr_crit >= 0, guarantees tmps /= 0 when the loop ends.
            IF ( tmps > sum05_entrdz .or.  &
              &  (kctop2 - kentrmax)*(kcbas2 - kentrmax) <= 0 )  EXIT
            kctop2 = jk + 1   ;   kcbas2 = kcbas0
            tmp = 0.0_wp      ;   tmps = 0.0_wp
          END IF
        ENDDO

        zm(jc,jb) = tmp/tmps
        zd(jc,jb) = SQRT( 5.0_wp*SUM(  entrdz(kctop2:kcbas2)        &
          &                           *( z_mc(jc,kctop2:kcbas2,jb)  &
          &                              - zm(jc,jb) )**2 )/tmps )
        ! zd >~ one-layer depth, at least
        zd(jc,jb) = MAX( dz(jc,kcbas2,jb), zd(jc,jb) )

        zctop(jc,jb) = zm(jc,jb) + zd(jc,jb)

        IF ( zctop(jc,jb) > z_ifc(jc,nlevp1,jb) + cfg% crit_zctop ) THEN
          entrmax(jc,jb) = MAXVAL(entr(kctop2:kcbas2))
        ELSE
          entrmax(jc,jb) = 0.0_wp
!         p_msgwam(jg)% flag_cgw(jc,jb) = -1  ! already initialized
          CYCLE
        END IF

        kcb = 2  ;  kcm = 2  ;  kct = 2

        tmp = zm(jc,jb) - zd(jc,jb)
        DO jk = nlev-1,2,-1
          IF (tmp < z_mc(jc,jk,jb)) THEN
            kcb = jk+1  ;  EXIT   ! closest interface level
          END IF
        ENDDO

!       DO jk = nlev-1,2,-1
        DO jk = kcb-1,2,-1
          IF (zm(jc,jb) < z_mc(jc,jk,jb)) THEN
            kcm = jk+1  ;  EXIT   ! closest interface level
          END IF
        ENDDO

!       DO jk = nlev-1,2,-1
        DO jk = kcm-1,2,-1
          IF (zctop(jc,jb) < z_mc(jc,jk,jb)) THEN
            kct = jk+1  ;  EXIT   ! closest interface level
          END IF
        ENDDO
        kcta(jc,jb) = kct

        zm(jc,jb) = zm(jc,jb) - z_ifc(jc,nlevp1,jb)

        ! force z_cbas >= 100 m to simplify the numerical calculation of spectra
        !   (cf. z_cbas = z_s = 0, which will be used as denominators).
        !   The slight change (100 m, in the PBL) will not affect the physics of
        !   our concern but rather help to avoid assigning unrealistically large
        !   values of the low-level shear.
        zd(jc,jb) = MIN( zm(jc,jb) - 100.0_wp, zd(jc,jb) )

        !-----------------------------------------------------------------------
        ! Obtain the parameters which construct the basic-state profile of SC05.
        !-----------------------------------------------------------------------

!       rho0(jc,jb) = 0.5_wp*(rho(jc,kct-1,jb) + rho(jc,kct,jb))
        rho0(jc,jb) = 0.5_wp*(rho(jc,kcm-1,jb) + rho(jc,kcm,jb))
        u_ct(jc,jb) = 0.5_wp*(u  (jc,kct-1,jb) + u  (jc,kct,jb))
        v_ct(jc,jb) = 0.5_wp*(v  (jc,kct-1,jb) + v  (jc,kct,jb))
        u_cb(jc,jb) = 0.5_wp*(u  (jc,kcb-1,jb) + u  (jc,kcb,jb))
        v_cb(jc,jb) = 0.5_wp*(v  (jc,kcb-1,jb) + v  (jc,kcb,jb))

        n2_ct = grav*LOG(    temp_ifc(jc,kct-1,jb)/temp_ifc(jc,kct+1,jb)                &
          &              *(  pres_ifc(jc,kct+1,jb)/pres_ifc(jc,kct-1,jb) )**rd_o_cpd )  &
          &         /(       z_ifc   (jc,kct-1,jb) - z_ifc (jc,kct+1,jb) )

        ! assuring kcb /= nlevp1 (guaranteed if dz < 100 m at jk = nlev)
        !   and kcb /= kct (guaranteed, given zd >~ dz)
        n2_q(jc,jb) = grav*LOG(    temp_ifc(jc,kct,jb)/temp_ifc(jc,kcb,jb)                &
          &                    *(  pres_ifc(jc,kcb,jb)/pres_ifc(jc,kct,jb) )**rd_o_cpd )  &
          &               /(       z_ifc   (jc,kct,jb) - z_ifc (jc,kcb,jb) )
        ! n2_q will be used as a denominator later
        IF (n2_q(jc,jb) == 0.0_wp)  n2_q(jc,jb) = -1.e-20_wp

        n2_ct_limit =  ri_bg_limit                               &
          &           *(   (u(jc,kct-1,jb) - u(jc,kct,jb))**2    &
          &              + (v(jc,kct-1,jb) - v(jc,kct,jb))**2 )  &
          &           /( z_mc(jc,kct-1,jb) - z_mc(jc,kct,jb) )**2

        IF (n2_ct <= n2_ct_limit) THEN
          entrmax(jc,jb) = 0.0_wp
!         p_msgwam(jg)% flag_cgw(jc,jb) = -1  ! already initialized
          CYCLE
        END IF

        nbv_ct(jc,jb) = SQRT(n2_ct)

        ! flag for output diagnosis
        p_msgwam(jg)% flag_cgw(jc,jb) = 1   ! entrmax(jc,jb) /= 0.0_wp

      END IF

    ENDDO  L_JC
    ENDDO  L_JB

!$OMP END DO
!$OMP END PARALLEL

    END IF   ! js == 1


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, iphi, i_startidx, i_endidx, phi0, phi, kll, kml)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
        &                 i_startidx, i_endidx, rl_start, rl_end )
      DO jc = i_startidx, i_endidx

        IF (p_msgwam(jg)% flag_cgw(jc,jb) < 0)  CYCLE

        !-----------------------------------------------------------------------
        ! Estimate moving velocity of convective cell (Choi and Chun, 2011, JAS)
        !-----------------------------------------------------------------------

        kll = p_msgwam(jg)% ll_k_cgw(jc,jb)

        ! low-level wind
        cqx(jc,jb) = (  SUM(  u(jc,kll:nlev,jb)*dz(jc,kll:nlev,jb) )     &
          &           + u(jc,kll-1,jb)*p_msgwam(jg)% ll_dz_cgw(jc,jb) )  &
          &          *rdz_ll
        cqy(jc,jb) = (  SUM(  v(jc,kll:nlev,jb)*dz(jc,kll:nlev,jb) )     &
          &           + v(jc,kll-1,jb)*p_msgwam(jg)% ll_dz_cgw(jc,jb) )  &
          &          *rdz_ll

        IF (cfg% scale_h(js) > 30.e3_wp) THEN
          ! mid-level wind PLUS ( mid-level wind MINUS low-level wind ) -- Corfidi (2003)
          kml = p_msgwam(jg)% ml_k_cgw(jc,jb)
          cqx(jc,jb) = (  SUM(  u(jc,kml:kll-2,jb)*dz(jc,kml:kll-2,jb) )                       &
            &           + u(jc,kml-1,jb)*p_msgwam(jg)% ml_dz_cgw(jc,jb)                        &
            &           + u(jc,kll-1,jb)*(dz(jc,kll-1,jb) - p_msgwam(jg)% ll_dz_cgw(jc,jb)) )  &
            &          *rdz_ml*2._wp - cqx(jc,jb)
          cqy(jc,jb) = (  SUM(  v(jc,kml:kll-2,jb)*dz(jc,kml:kll-2,jb) )                       &
            &           + v(jc,kml-1,jb)*p_msgwam(jg)% ml_dz_cgw(jc,jb)                        &
            &           + v(jc,kll-1,jb)*(dz(jc,kll-1,jb) - p_msgwam(jg)% ll_dz_cgw(jc,jb)) )  &
            &          *rdz_ml*2._wp - cqy(jc,jb)
        END IF

        IF ( cqx(jc,jb) /= u_ct(jc,jb) .or. cqy(jc,jb) /= v_ct(jc,jb) ) THEN
          phi0 = ATAN2( cqy(jc,jb) - v_ct(jc,jb), cqx(jc,jb) - u_ct(jc,jb) )
          DO iphi = 1, nphi05
            phi = phi0 + dphi*REAL(iphi-1,kind=wp)
            cosphi(iphi,jc,jb) = COS(phi)
            sinphi(iphi,jc,jb) = SIN(phi)
            IF (cosphi(iphi,jc,jb) < 0.0_wp) THEN
              ! -pi/2 <= phi <= pi/2 (i.e., cosphi(1:nphi05) >= 0)
              cosphi(iphi,jc,jb) = -cosphi(iphi,jc,jb)
              sinphi(iphi,jc,jb) = -sinphi(iphi,jc,jb)
            END IF
          ENDDO
        END IF

        a_mwn(:,jc,jb) = a_mwn_o_n_s(1:cfg%nm(js))*nbv_ct(jc,jb)
        dm   (:,jc,jb) = dm_o_n_s   (1:cfg%nm(js))*nbv_ct(jc,jb)

      ENDDO
    ENDDO

!$OMP END DO
!$OMP END PARALLEL

    cosphi(nphi05+1:cfg%nphi(js),:,:) = -cosphi(1:nphi05,:,:)
    sinphi(nphi05+1:cfg%nphi(js),:,:) = -sinphi(1:nphi05,:,:)


    !:::::   LOOP ENDS   :::::
    CALL timer_stop(timer_gws_test1)


    CALL timer_start(timer_gws_test2)

    CALL gw_source_conv_calc ( nproma, p_patch%nblks_c, i_startblk, i_endblk,  &
      &                        cfg% nm(js), cfg% nphi(js),                     &
      &                        cfg% scale_h(js), cfg% scale_t(js),             &
      &                        cfg% heat_areafrac(js),                         &
      &                        coef_gridmeanflux(js),                          &
      &                        cosphi, sinphi, a_mwn,                          &
      &                        entrmax, zm, zd, cqx, cqy,                      &
      &                        nbv_ct, n2_q, rho0,                             &
      &                        u_ct, v_ct, u_cb, v_cb,                         &
      &                        u(:,nlev,:), v(:,nlev,:),                       &
      &                        mfs_ct, kh_ct, dkh_ct )
    CALL timer_stop(timer_gws_test2)

    IF (msg_level >= 15) THEN
      CALL print_value('CGW flux spectral max. : ',MAXVAL(mfs_ct))
    END IF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, iphi, i_startidx, i_endidx, dmdphi)

    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
        &                 i_startidx, i_endidx, rl_start, rl_end )

    DO jc = i_startidx, i_endidx

      dmdphi(:) = dm(:,jc,jb)*dphi

      DO iphi = 1, nphi05
        ! -pi/2 <= phi <= pi/2 (i.e., cosphi(1:nphi05) >= 0)
        p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+1) = p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+1)  &
          &  + SUM(mfs_ct(:,iphi,jc,jb)*dmdphi(:))*cosphi(iphi,jc,jb)
        p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+2) = p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+2)  &
          &  + SUM(mfs_ct(:,nphi05+iphi,jc,jb)*dmdphi(:))*cosphi(iphi,jc,jb)
        IF (sinphi(iphi,jc,jb) > 0.0_wp) THEN
          p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+3) = p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+3)  &
            &  + SUM(mfs_ct(:,iphi,jc,jb)*dmdphi(:))*sinphi(iphi,jc,jb)
          p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+4) = p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+4)  &
            &  + SUM(mfs_ct(:,nphi05+iphi,jc,jb)*dmdphi(:))*sinphi(iphi,jc,jb)
        ELSE
          p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+3) = p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+3)  &
            &  + SUM(mfs_ct(:,nphi05+iphi,jc,jb)*dmdphi(:))*(-sinphi(iphi,jc,jb))
          p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+4) = p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+4)  &
            &  + SUM(mfs_ct(:,iphi,jc,jb)*dmdphi(:))*(-sinphi(iphi,jc,jb))
        END IF
      ENDDO

      p_msgwam(jg)% ctmfl_cgw(jc,jb,js*5+5) = SUM( SUM(mfs_ct(:,:,jc,jb), dim=2)*dmdphi(:) )

    ENDDO
    ENDDO

!$OMP END DO
!$OMP END PARALLEL

    IF (js == cfg% n_source(jg)) THEN

      p_msgwam(jg)% ctmfl_cgw(:,:,1:5) = 0._wp

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, jd, jss, i_startidx, i_endidx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
          &                 i_startidx, i_endidx, rl_start, rl_end )
        DO jd = 1, 5
          DO jss = 1, cfg% n_source(jg)
            DO jc = i_startidx, i_endidx
              p_msgwam(jg)% ctmfl_cgw(jc,jb,jd) =  p_msgwam(jg)% ctmfl_cgw(jc,jb,jd)  &
                &                                + p_msgwam(jg)% ctmfl_cgw(jc,jb,jss*5+jd)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      ! output diagnosis (not used during the integration)
      p_msgwam(jg)% pars_cgw(:,:,1 ) = entrmax(:,:)
      p_msgwam(jg)% pars_cgw(:,:,2 ) = zm     (:,:)
      p_msgwam(jg)% pars_cgw(:,:,3 ) = zd     (:,:)
      p_msgwam(jg)% pars_cgw(:,:,4 ) = REAL(kcta(:,:),kind=wp)
      p_msgwam(jg)% pars_cgw(:,:,5 ) = nbv_ct(:,:)**2
      p_msgwam(jg)% pars_cgw(:,:,6 ) = n2_q   (:,:)
      p_msgwam(jg)% pars_cgw(:,:,7 ) = rho0   (:,:)
      p_msgwam(jg)% pars_cgw(:,:,8 ) = u_ct   (:,:)
      p_msgwam(jg)% pars_cgw(:,:,9 ) = v_ct   (:,:)
      p_msgwam(jg)% pars_cgw(:,:,10) = u_cb   (:,:)
      p_msgwam(jg)% pars_cgw(:,:,11) = v_cb   (:,:)
      p_msgwam(jg)% pars_cgw(:,:,12) = u(:,nlev,:)
      p_msgwam(jg)% pars_cgw(:,:,13) = v(:,nlev,:)

    END IF

    p_msgwam(jg)% pars_cgw(:,:,12+js*2) = cqx(:,:)
    p_msgwam(jg)% pars_cgw(:,:,13+js*2) = cqy(:,:)


    CALL timer_start(timer_gws_test3)

    !   ctmfl_cgw_spec
    diag_code_10  = MOD(cfg% diag_code, 100)/10
    diag_code_100 = MOD(cfg% diag_code, 1000)/100

    IF ( diag_code_10 /= 0 .OR. diag_code_100 /= 0 ) THEN

      IF (js == 1) THEN
        p_msgwam(jg)% ctmfl_specc_cgw(:,:,:) = 0._wp
        p_msgwam(jg)% ctmfl_spec0_cgw(:,:,:) = 0._wp
        p_msgwam(jg)% ctmfl_spec1_cgw(:,:,:) = 0._wp
      END IF

      CALL diag_spectrum ( p_patch, p_msgwam(jg),                           &
        &                  rl_start, rl_end, i_startblk, i_endblk, nlevp1,  &
        &                  diag_code_10, diag_code_100,                     &
        &                  cfg% nm(js), cfg% nphi(js),                      &
        &                  cfg% scale_h(js), cfg% scale_t(js),              &
        &                  cosphi, sinphi,                                  &
        &                  mfs_ct, kh_ct, dm,                               &
        &                  cqx, cqy, u_ct, v_ct )
    END IF

    CALL timer_stop(timer_gws_test3)

    IF (nrays(jg) == 0) THEN   ! lmsgwam == .FALSE.
      IF (js == cfg% n_source(jg))  CALL dealloc_local_var()
      RETURN
    END IF


    IF (js == 1) THEN
      jray_offset = 0
    ELSE
      jray_offset = SUM( cfg% nm(1:js-1) * cfg% nphi(1:js-1) )
    END IF

    p_ray_conv(jg)%wadens(:, jray_offset+1:jray_offset+cfg%nm(js)*cfg%nphi(js), :) = 0.0_wp

    ! p_ray_conv(jg)%jr_leader should not be reset yet since it will be used to
    !   tailor the tails of the previous rays

    IF (.NOT. lsteady) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, iphi, i_startidx, i_endidx, iray1, iray9, jac, sqr_dkdl)

      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
          &                 i_startidx, i_endidx, rl_start, rl_end )

      DO jc = i_startidx, i_endidx

!       IF (p_msgwam(jg)% flag_cgw(jc,jb) < 0)  CYCLE   ! entrmax(jc,jb) == 0.0_wp
        IF (entrmax(jc,jb) == 0.0_wp)  CYCLE

        DO iphi = 1, cfg% nphi(js)

          iray1 = jray_offset + (iphi-1)*cfg% nm(js) + 1
          iray9 = jray_offset + iphi*cfg% nm(js)

          ! F(phi,m) --> F(k,l,m)
          jac(:) = kh_ct(:,iphi,jc,jb)*dkh_ct(:,iphi,jc,jb)

          sqr_dkdl(:) = SQRT(jac(:)*dphi)

          p_ray_conv(jg)%k (jc,iray1:iray9,jb) = kh_ct(:,iphi,jc,jb)  &
            &                                    *cosphi(iphi,jc,jb)
          p_ray_conv(jg)%l (jc,iray1:iray9,jb) = kh_ct(:,iphi,jc,jb)  &
            &                                    *sinphi(iphi,jc,jb)
          p_ray_conv(jg)%dk(jc,iray1:iray9,jb) = sqr_dkdl(:)
          p_ray_conv(jg)%dl(jc,iray1:iray9,jb) = sqr_dkdl(:)
          p_ray_conv(jg)%m (jc,iray1:iray9,jb) = -a_mwn(:,jc,jb)
          p_ray_conv(jg)%dm(jc,iray1:iray9,jb) = dm    (:,jc,jb)

          ! A = F/(kh*cgz) ~ F*(m/kh)**2/N
          !   where cgz ~ N*kh/(m**2) (hydrostatic approx.)
          p_ray_conv(jg)%wadens(jc,iray1:iray9,jb) = mfs_ct(:,iphi,jc,jb)  &
            &      /jac(:)*(a_mwn(:,jc,jb)/kh_ct(:,iphi,jc,jb))**2         &
            &             /nbv_ct(jc,jb)

        ENDDO

      ENDDO
      ENDDO

!$OMP END DO
!$OMP END PARALLEL

      IF (js == cfg% n_source(jg)) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx)
        DO jb = i_startblk, i_endblk
          CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
            &                 i_startidx, i_endidx, rl_start, rl_end )
          DO jc = i_startidx, i_endidx
            IF ( p_msgwam(jg)% flag_cgw(jc,jb) == SIGN(1,flag_last(jc,jb)) )  &
              &  p_msgwam(jg)% flag_cgw(jc,jb) = p_msgwam(jg)% flag_cgw(jc,jb)*2
          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END IF

    ELSE  ! lsteady

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, iphi, i_startidx, i_endidx, iray1, iray9)

      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
          &                 i_startidx, i_endidx, rl_start, rl_end )

      DO jc = i_startidx, i_endidx

!       IF (p_msgwam(jg)% flag_cgw(jc,jb) < 0)  CYCLE   ! entrmax(jc,jb) == 0.0_wp
        IF (entrmax(jc,jb) == 0.0_wp)  CYCLE

        DO iphi = 1, cfg% nphi(js)

          iray1 = jray_offset + (iphi-1)*cfg% nm(js) + 1
          iray9 = jray_offset + iphi*cfg% nm(js)

          p_ray_conv(jg)%k(jc,iray1:iray9,jb) = kh_ct(:,iphi,jc,jb)  &
            &                                   *cosphi(iphi,jc,jb)
          p_ray_conv(jg)%l(jc,iray1:iray9,jb) = kh_ct(:,iphi,jc,jb)  &
            &                                   *sinphi(iphi,jc,jb)
          p_ray_conv(jg)%m(jc,iray1:iray9,jb) = -a_mwn(:,jc,jb)

          ! A = F/(kh*cgz) ~ F*(m/kh)**2/N
          !   where cgz ~ N*kh/(m**2) (hydrostatic approx.)
          p_ray_conv(jg)%wadens(jc,iray1:iray9,jb) = mfs_ct(:,iphi,jc,jb)  &
            &             *(a_mwn(:,jc,jb)/kh_ct(:,iphi,jc,jb))**2         &
            &             /nbv_ct(jc,jb)                                   &
            &             *(dphi*dm(:,jc,jb))

        ENDDO

      ENDDO
      ENDDO

!$OMP END DO
!$OMP END PARALLEL

    END IF  ! lsteady


    IF (js == 1) THEN

      p_ray_conv(jg)%jk_source(:,:) = 1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx)
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
          &                 i_startidx, i_endidx, rl_start, rl_end )
        DO jc = i_startidx, i_endidx

!         IF (p_msgwam(jg)% flag_cgw(jc,jb) < 0)  CYCLE   ! entrmax(jc,jb) == 0.0_wp
          IF (entrmax(jc,jb) == 0.0_wp)  CYCLE

          ! jk_source : closest interface-level index closest to the launch level
          p_ray_conv(jg)%jk_source(jc,jb) = kcta(jc,jb)

!YH       ! cloud center or cloud top ?
!         p_ray_conv(jg)%jk_source(jc,jb) =

        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    END IF

    IF (js == cfg% n_source(jg))  CALL dealloc_local_var()

  END SUBROUTINE gw_source_conv

  SUBROUTINE gw_source_conv_calc ( ncol, nblk, i_startblk, i_endblk,         &
    &                              nm, nphi,                                 &
    &                              hscale, tscale, areafrac, coef_gridmean,  &
    &                              cosphi, sinphi, a_mwn,                    &
    &                              entrmax, zm, zd, cqx, cqy,                &
    &                              nbv_ct, n2_q, rho0,                       &
    &                              u_ct, v_ct, u_cb, v_cb, u_sfc, v_sfc,     &
    &                              mfs, kh, dkh )

    INTEGER , INTENT(in)  ::  ncol, nblk, i_startblk, i_endblk
    INTEGER , INTENT(in)  ::  nm, nphi
    REAL(wp), INTENT(in)  ::  hscale, tscale, areafrac, coef_gridmean
    REAL(wp), INTENT(in)  ::  cosphi(:,:,:),  &   ! (nphi,ncol,nblk)
      &                       sinphi(:,:,:)
    REAL(wp), INTENT(in)  ::  a_mwn(:,:,:)        ! (nm,ncol,nblk)
    REAL(wp), INTENT(in)  ::  entrmax(:,:),  &    ! (ncol,nblk)
      &                       zm     (:,:),  &
      &                       zd     (:,:),  &
      &                       cqx    (:,:),  &
      &                       cqy    (:,:),  &
      &                       nbv_ct (:,:),  &
      &                       n2_q   (:,:),  &
      &                       rho0   (:,:),  &
      &                       u_ct   (:,:),  &
      &                       v_ct   (:,:),  &
      &                       u_cb   (:,:),  &
      &                       v_cb   (:,:),  &
      &                       u_sfc  (:,:),  &
      &                       v_sfc  (:,:)
    REAL(wp), INTENT(out) ::  mfs(:,:,:,:),  &    ! (nm,nphi,ncol,nblk)
      &                       kh (:,:,:,:),  &
      &                       dkh(:,:,:,:)

    ! data arrays
    REAL(wp) ::  x_sq(nm,2), c_intr(nm,2), th_ftn(nm,2)
    REAL(wp) ::  c_intr_twist, tmp_m(nm)
    REAL(wp) ::  q0sqc(ncol,nblk), zcta(ncol,nblk), zcba(ncol,nblk)

    ! variables used to calculate |X|^2
    REAL(wp)    ::  ub_sfc, ub_cb, zcsa, shear, shear_sq
    REAL(dp)    ::  n2dn1, inv_zd, ri_cs, mu, ztm_zs, ztqs, ztqdzb, c2drip2dzd2
    REAL(dp)    ::  mzalph, zb_mzalph, zs_mzalph, lm1, c2zs, c2zb, ztss,  &
      &             inv_zdlm1, ztu, y2r, xpn(2)
    COMPLEX(dp) ::  bfac, twoimu, c05p_aimu, c05m_aimu
    COMPLEX(dp) ::  denom, x0, x2, x3, x4, x2c, x3c, x4c, xa, xac, term_xa,  &
      &             y1, y3, y1c, term_y1, coef2, coef3, ctmp

    ! work arrays
    REAL(wp) ::  cqh, c0, coef0, tmp, tmp2, coef_kh
    INTEGER  ::  nphi05
    INTEGER  ::  im, iphi, jc, jb, iposneg   ! loop counters

    !----------------------------------------------------------------------
    ! Purpose:
    !         Calculate forcing of GWs by convection for MS-GWaM
    ! Method:
    !         -- see Song and Chun (2005)
    !----------------------------------------------------------------------

    ! initialization
    mfs(:,:,:,:) = 0.0_wp
    kh (:,:,:,:) = 1.0e20_wp    ! can be used as a denominator
    dkh(:,:,:,:) = 1.0e-20_wp   ! can be used as a denominator

    ! constants
    nphi05  = nphi/2
    c0      = hscale/tscale
    coef0   = 0.125_wp*(hscale*tscale)*(grav/cpd)/areafrac
    coef_kh = 2.0_wp/hscale           ! use of the kh peak
!   coef_kh = SQRT(pi*2.0_wp)/hscale  ! use of the kh mean

    c_intr(:,1) = c_intr_ct_s(1:nm)
    c_intr(:,2) = -c_intr_ct_s(1:nm)

    ! loop for jc,jb
    q0sqc(:,:) = coef_gridmean*rho0(:,:)*(coef0*entrmax(:,:)/n2_q(:,:))**2
    zcta (:,:) = zm(:,:) + zd(:,:)
    zcba (:,:) = zm(:,:) - zd(:,:)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, iphi, im, iposneg,                                                &
!$OMP            x_sq, th_ftn, c_intr_twist, tmp_m, ub_sfc, ub_cb, zcsa, shear, shear_sq,  &
!$OMP            n2dn1, inv_zd, ri_cs, mu, ztm_zs, ztqs, ztqdzb, c2drip2dzd2, mzalph,      &
!$OMP            zb_mzalph, zs_mzalph, lm1, c2zs, c2zb, ztss, inv_zdlm1, ztu, y2r, xpn,    &
!$OMP            bfac, twoimu, c05p_aimu, c05m_aimu,                                       &
!$OMP            denom, x0, x2, x3, x4, x2c, x3c, x4c, xa, xac, term_xa,                   &
!$OMP            y1, y3, y1c, term_y1, coef2, coef3, ctmp, cqh, tmp, tmp2)

    L_JB:  DO jb = i_startblk, i_endblk
    L_JC:  DO jc = 1, ncol

      IF (entrmax(jc,jb) == 0.0_wp)  CYCLE  ! mfs(:,:,jc,jb) = 0 (initialized)


    L_PHI:  DO iphi = 1, nphi05

      !-------------------------------------------------------------------------
      ! Obtain the basic-state profile for the analytic solution of SC05.
      !-------------------------------------------------------------------------

      ub_cb  =   (u_cb (jc,jb) - u_ct(jc,jb))*cosphi(iphi,jc,jb)  &
        &      + (v_cb (jc,jb) - v_ct(jc,jb))*sinphi(iphi,jc,jb)
      ub_sfc =   (u_sfc(jc,jb) - u_ct(jc,jb))*cosphi(iphi,jc,jb)  &
        &      + (v_sfc(jc,jb) - v_ct(jc,jb))*sinphi(iphi,jc,jb)

      zcsa = zcta(jc,jb)
      tmp = ub_cb - ub_sfc
      IF ( ABS(tmp) > 0.1_wp ) THEN
        tmp2 = -ub_sfc/tmp
        IF (tmp2 > 0.0_wp) THEN
          zcsa = tmp2*zcba(jc,jb)
          zcsa = MIN( zcta(jc,jb), MAX( zcsa, zcba(jc,jb) ) )
        END IF
      END IF
      shear = -ub_sfc/zcsa
      shear_sq = shear*shear

      ub_cb = ub_sfc + shear*zcba(jc,jb)   ! adjustment

      !-------------------------------------------------------------------------
      ! Calculate |X|^2, theta function, and momentum flux
      !-------------------------------------------------------------------------

      x_sq  (:,:) = 0.0_wp
      th_ftn(:,:) = 0.0_wp

      IF ( n2_q(jc,jb) > 1000.0_wp*shear_sq ) THEN

        ! |X|^2 for the uniform wind

        n2dn1  = REAL(nbv_ct(jc,jb)/SQRT(n2_q(jc,jb)),kind=dp)
        inv_zd = 1.0_dp/REAL(zd(jc,jb),kind=dp)
        ztm_zs = REAL(zd(jc,jb) + zd(jc,jb),kind=dp)

        DO im = 1, nm
          lm1        = REAL(a_mwn(im,jc,jb),kind=dp)/n2dn1
          inv_zdlm1  = inv_zd/lm1
          y1         = CMPLX(inv_zdlm1*inv_zdlm1,-inv_zdlm1,kind=dp)
          x2         = EXP(CMPLX(0.0_dp,lm1*ztm_zs,kind=dp))
          x3         = EXP(CMPLX(0.0_dp,lm1*REAL(zcta(jc,jb),kind=dp),kind=dp))
          x_sq(im,1) = REAL( 4.0_dp*REAL(x3*(y1 - CONJG(x2*y1)))**2  &
            &                /( REAL(x3)**2 + (n2dn1*AIMAG(x3))**2 ), kind=wp )
        ENDDO

!im = maxloc(x_sq(:,1),1)
!if (x_sq(im,1) > 30. .and. nbv_ct(jc,jb) > 0.5e-2_wp) then
!print*,"YH,X2 ",x_sq(im,1),c_intr(im,1)
!print*,"YH, zd,zb",zd(jc,jb),zcba(jc,jb)
!print*,"YH, nt,n2q",nbv_ct(jc,jb),n2_q(jc,jb)
!print*,"YH, c_min",c_intr(nm,1)
!end if

        ! Since (intrinsic) cqh = 0, spectra become symmetric.

        th_ftn(:,1) = 1.0_wp/(1.0_wp + (c_intr(:,1)/c0)**2)

        mfs(:,iphi,jc,jb) = x_sq(:,1)*th_ftn(:,1)*c_intr(:,1)*q0sqc(jc,jb)

        kh (:,iphi,jc,jb) = coef_kh*SQRT(th_ftn(:,1))
        dkh(:,iphi,jc,jb) = kh(:,iphi,jc,jb)

        mfs(:,nphi05+iphi,jc,jb) = mfs(:,iphi,jc,jb)
        kh (:,nphi05+iphi,jc,jb) = kh (:,iphi,jc,jb)
        dkh(:,nphi05+iphi,jc,jb) = dkh(:,iphi,jc,jb)

      ELSE IF ( n2_q(jc,jb) > ri_bg_limit*shear_sq ) THEN

        ! |X|^2 for the sheared wind

        n2dn1       = REAL(nbv_ct(jc,jb)/SQRT(n2_q(jc,jb)),kind=dp)
        inv_zd      = 1.0_dp/REAL(zd(jc,jb),kind=dp)
        bfac        = CMPLX(eipi_r,SIGN(eipi_i,REAL(-shear,kind=dp)),kind=dp)
        ri_cs       = REAL(n2_q(jc,jb)/shear_sq,kind=dp)
        mu          = SQRT(ri_cs - 0.25_dp)
        twoimu      = CMPLX(0.0_dp,mu+mu,kind=dp)
        c05p_aimu   = CMPLX(0.5_dp,mu ,kind=dp)
        c05m_aimu   = CMPLX(0.5_dp,-mu,kind=dp)
        ztm_zs      = REAL(zcta(jc,jb) - zcsa,kind=dp)
        ztqs        = 1.0_dp - (REAL(zcsa - zm(jc,jb),kind=dp)*inv_zd)**2
        ztqdzb      = inv_zd+inv_zd
        c2drip2dzd2 = ztqdzb*inv_zd/(ri_cs + 2.0_dp)

        xpn = (/n2dn1,-n2dn1/)

        DO iposneg = 1, 2
        DO im = 1, nm

          ! --- Singular point ---
          ! c_intr values are changed (by a negligible amount) if it is too
          !   close to ub_cb or ub_sfc in order to avoid singularity
          !   (zs_mzalph = 0 or mzalph = 0, respectively).
          c_intr_twist = c_intr(im,iposneg)
          IF ( ABS(c_intr(im,iposneg) - ub_sfc) < 0.999e-3_wp ) THEN
            c_intr_twist = ub_sfc + SIGN(0.999e-3_wp,-shear)
          ELSE IF ( ABS(c_intr(im,iposneg) - ub_cb) < 0.999e-3_wp ) THEN
            c_intr_twist = ub_cb + SIGN(0.999e-3_wp,-shear)
          END IF

          mzalph    = REAL((ub_sfc - c_intr_twist)/shear,kind=dp)
          zb_mzalph = REAL(zcba(jc,jb),kind=dp) + mzalph
          zs_mzalph = REAL(zcsa,kind=dp) + mzalph
          lm1       = REAL(a_mwn(im,jc,jb),kind=dp)/n2dn1
          ctmp      = CMPLX(0.0_dp,lm1*zs_mzalph,kind=dp)
          coef2     = ctmp - c05m_aimu
          c2zs      = c2drip2dzd2*zs_mzalph
          c2zb      = c2drip2dzd2*zb_mzalph
          ztss      = c2zs*zs_mzalph
          inv_zdlm1 = inv_zd/lm1
          ztu       = inv_zdlm1*inv_zdlm1*2.0_dp
          term_y1   = (   c05m_aimu*(c2zb*zb_mzalph)  &
            &           - (c2zb+c2zb + ztqdzb)*zb_mzalph )/twoimu
          IF (mzalph > 0.0_dp) THEN
            x0 = mzalph**twoimu
          ELSE
            x0 = ((-mzalph)*bfac)**twoimu
          END IF
          IF (zb_mzalph > 0.0_dp) THEN
            y1  = zb_mzalph**(-c05p_aimu)*term_y1
            y1c = CONJG(y1)
          ELSE
            y1  = ((-zb_mzalph)*bfac)**(-c05p_aimu)*term_y1
            y1c = ((-zb_mzalph)*bfac)**(-c05m_aimu)*CONJG(term_y1)
          END IF
          IF (zs_mzalph > 0.0_dp) THEN
            x2  = zs_mzalph**c05p_aimu
            x3  = CONJG(x2)*(ctmp - c05p_aimu)
            x2  = x2*coef2
            x2c = CONJG(x2)
            x3c = CONJG(x3)
          ELSE
            coef3 = ctmp - c05p_aimu
            x2  = ((-zs_mzalph)*bfac)**c05p_aimu
            x3  = ((-zs_mzalph)*bfac)**c05m_aimu
            x2c = x3*CONJG(coef2)
            x3c = x2*CONJG(coef3)
            x2  = x2*coef2
            x3  = x3*coef3
          END IF
          y2r      = (ztqs + ztss) - zs_mzalph*(c2zs+c2zs)
          y3       = CMPLX(ztu,inv_zdlm1+inv_zdlm1,kind=dp)
          term_xa  = (y2r - c05m_aimu*(ztss - ztu))/coef2
          xa       = x2 *term_xa        + twoimu*zs_mzalph*y1c
          xac      = x2c*CONJG(term_xa) - twoimu*zs_mzalph*y1
          x4       = 2.0_dp*EXP(CMPLX(0.0_dp,lm1*ztm_zs,kind=dp))
          x4c      = CONJG(x4)*(x3c - x2c*x0)
          x4       = x4*(x2 - x3 *x0)
          denom    = (1.0_dp + xpn(iposneg))*x4 + (1.0_dp - xpn(iposneg))*x4c
          x_sq(im,iposneg) = REAL( ABS((x4*CONJG(y3) + x4c*y3 + 4.0_dp*(xa - xac*x0))  &
            &                      /denom)**2, kind=wp )

        ENDDO  ! im
        ENDDO  ! iposneg

!do iposneg = 1, 2
!im = maxloc(x_sq(:,iposneg),1)
!if (x_sq(im,iposneg) > 30. .and. nbv_ct(jc,jb) > 0.5e-2_wp) then
!print*,"YH,sh,X2 ",x_sq(im,iposneg),c_intr(im,iposneg)
!print*,"YH,sh zd,zb",zd(jc,jb),zcba(jc,jb)
!print*,"YH,sh nt,n2q",nbv_ct(jc,jb),n2_q(jc,jb)
!print*,"YH,sh c_min",c_intr(nm,iposneg)
!end if
!enddo

        cqh =   (cqx(jc,jb) - u_ct(jc,jb))*cosphi(iphi,jc,jb)  &
          &   + (cqy(jc,jb) - v_ct(jc,jb))*sinphi(iphi,jc,jb)

        th_ftn(:,:) = 1.0_wp/(1.0_wp + ((c_intr(:,:) - cqh)/c0)**2)

        tmp_m(:) = c_intr(:,1)*q0sqc(jc,jb)
        ! c_intr(:,1) = ABS(c_intr(:,1)) = ABS(c_intr(:,2))

        mfs(:,iphi       ,jc,jb) = x_sq(:,1)*th_ftn(:,1)*tmp_m(:)
        mfs(:,nphi05+iphi,jc,jb) = x_sq(:,2)*th_ftn(:,2)*tmp_m(:)

        kh (:,iphi       ,jc,jb) = coef_kh*SQRT(th_ftn(:,1))
        kh (:,nphi05+iphi,jc,jb) = coef_kh*SQRT(th_ftn(:,2))
        dkh(:,iphi       ,jc,jb) = kh(:,iphi       ,jc,jb)
        dkh(:,nphi05+iphi,jc,jb) = kh(:,nphi05+iphi,jc,jb)

      END IF  ! n2_q(jc,jb) ; 1000.0_wp*shear_sq ; ri_bg_limit*shear_sq

      ! If the case does not conform to any one above
      !   (i.e., n2_q(jc,jb) <= ri_bg_limit*shear_sq), no available solution
      !   exists: mfs(:,(/iphi,nphi05+iphi/),jc,jb) = 0 (already initialized).

    ENDDO  L_PHI

    ENDDO  L_JC
    ENDDO  L_JB

!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE gw_source_conv_calc

  SUBROUTINE diag_spectrum ( p_patch, field,                                  &
    &                        rl_start, rl_end, i_startblk, i_endblk, nlevp1,  &
    &                        diag_code_10, diag_code_100,                     &
    &                        nm, nphi,                                        &
    &                        hscale, tscale,                                  &
    &                        cosphi, sinphi,                                  &
    &                        mfs_ct, kh_ct, dm,                               &
    &                        cqx, cqy, u_ct, v_ct )

    TYPE(t_patch)  , TARGET, INTENT(in)    ::  p_patch            !< grid/patch info.
    TYPE(t_msgwam)         , INTENT(inout) ::  field

    REAL(wp), INTENT(in) ::  mfs_ct(:,:,:,:),  &
      &                      kh_ct (:,:,:,:)
    REAL(wp), INTENT(in) ::  dm(:,:,:)

    REAL(wp), INTENT(in) ::  cqx (:,:),  &
      &                      cqy (:,:),  &
      &                      u_ct(:,:),  &
      &                      v_ct(:,:)

    REAL(wp), INTENT(in) ::  cosphi(:,:,:),  &
      &                      sinphi(:,:,:)

    INTEGER , INTENT(in) ::  rl_start, rl_end
    INTEGER , INTENT(in) ::  i_startblk, i_endblk    !< blocks
    INTEGER , INTENT(in) ::  nlevp1
    INTEGER , INTENT(in) ::  diag_code_10, diag_code_100
    INTEGER , INTENT(in) ::  nm, nphi
    REAL(wp), INTENT(in) ::  hscale, tscale

    REAL(wp), PARAMETER ::  dc_bin = 2.0_wp
    INTEGER , PARAMETER ::  nc_bin = 40
    INTEGER , PARAMETER ::  nc_bin_neg = 20
    INTEGER , PARAMETER ::  nkh_bin = 40
    INTEGER , PARAMETER ::  noi_bin = 40
    REAL(wp), PARAMETER ::  kh_edge0 = 1.e-7_wp
    REAL(wp), PARAMETER ::  oi_edge0 = 1.e-6_wp

    INTEGER  ::  i_startidx(i_startblk:i_endblk)
    INTEGER  ::  i_endidx  (i_startblk:i_endblk)

    REAL(wp) ::  dphi
    REAL(wp) ::  sp_dphi_o_dc
    REAL(wp) ::  sp_a_inv_cosphi
    REAL(wp) ::  sp_uk
    REAL(wp) ::  sp_hscale2o8, sp_tscale2o8
    REAL(wp) ::  sp_cqh
    REAL(wp) ::  sp_kh_bin_edge (0:nkh_bin)
    REAL(wp) ::  sp_kh2_bin_edge(0:nkh_bin)
    REAL(wp) ::  sp_intg_kh     (0:nkh_bin)
    REAL(wp) ::  sp_dphi_o_dkh(nkh_bin)
    REAL(wp) ::  sp_oi_bin_edge (0:noi_bin)
    REAL(wp) ::  sp_oi2_bin_edge(0:noi_bin)
    REAL(wp) ::  sp_intg_oi     (0:noi_bin)
    REAL(wp) ::  sp_dphi_o_doi(noi_bin)
    REAL(wp) ::  sp_term1k(nm), sp_term1o(nm), sp_termc(nm)
    REAL(wp) ::  sp_kh_bin_edge2(0:noi_bin,nm)
    REAL(wp) ::  tmp
    INTEGER  ::  jk_add

    INTEGER  ::  jc, jk, jb, iphi, im

    INTEGER , PARAMETER ::  ic_bin_neg_offset = nc_bin_neg*2 + nc_bin
    REAL(wp), PARAMETER ::  inv_dc_bin = 1._wp/dc_bin
    REAL(wp), PARAMETER ::  kh_edge0_rad = 2._wp*pi*kh_edge0
    REAL(wp), PARAMETER ::  oi_edge0_rad = 2._wp*pi*oi_edge0

    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,  &
        &                 i_startidx(jb), i_endidx(jb), rl_start, rl_end )
    ENDDO

    dphi = pi/REAL(nphi/2, kind=wp)

    IF (diag_code_10 == 1) THEN

      sp_dphi_o_dc = dphi*inv_dc_bin

      DO jb = i_startblk, i_endblk
        DO jc = i_startidx(jb), i_endidx(jb)
        DO iphi = 1, nphi
          sp_a_inv_cosphi = 1._wp/MAX(1.e-30_wp, ABS(cosphi(iphi,jc,jb)))
          sp_uk = u_ct(jc,jb)*cosphi(iphi,jc,jb) + v_ct(jc,jb)*sinphi(iphi,jc,jb)
          IF (cosphi(iphi,jc,jb) > 0._wp) THEN
            jk_add = nc_bin_neg
          ELSE
            jk_add = ic_bin_neg_offset
          END IF
          DO im = 1, nm
            tmp = (c_intr_ct_s(im) + sp_uk)*sp_a_inv_cosphi*inv_dc_bin
            IF ( ABS(tmp) > 16000. )  CYCLE  ! as ceiling does not work when cosphi ~ 0
            jk = CEILING(tmp)
            IF ( jk > -nc_bin_neg .AND. jk <= nc_bin )          &
              &  field%     ctmfl_specc_cgw(jc,jk+jk_add,jb)  &
              &    = field% ctmfl_specc_cgw(jc,jk+jk_add,jb)  &
              &      + mfs_ct(im,iphi,jc,jb)*dm(im,jc,jb)*sp_dphi_o_dc
          ENDDO
        ENDDO
        ENDDO
      ENDDO

    END IF

!   ctmfl_cgw_spec
    IF (diag_code_100 /= 0) THEN
      DO jk = 0, nkh_bin
        sp_kh_bin_edge(jk) = kh_edge0_rad*10._wp**(REAL(jk,kind=wp)*0.1_wp)
      ENDDO
      DO jk = 0, noi_bin
        sp_oi_bin_edge(jk) = oi_edge0_rad*10._wp**(REAL(jk,kind=wp)*0.1_wp)
      ENDDO
      sp_dphi_o_dkh(:) = dphi/(   sp_kh_bin_edge(1:nkh_bin  )  &
        &                       - sp_kh_bin_edge(0:nkh_bin-1) )
      sp_dphi_o_doi(:) = dphi/(   sp_oi_bin_edge(1:noi_bin  )  &
        &                       - sp_oi_bin_edge(0:noi_bin-1) )
    END IF
    IF ( diag_code_100 == 1 .OR. diag_code_100 == 3 ) THEN

      sp_hscale2o8 = 0.125*hscale**2
      sp_tscale2o8 = 0.125*tscale**2
      sp_kh2_bin_edge(:) = sp_kh_bin_edge(:)**2
      sp_oi2_bin_edge(:) = sp_oi_bin_edge(:)**2

      DO jb = i_startblk, i_endblk
        DO jc = i_startidx(jb), i_endidx(jb)
        DO iphi = 1, nphi
          sp_cqh =   (cqx(jc,jb) - u_ct(jc,jb))*cosphi(iphi,jc,jb)  &
            &      + (cqy(jc,jb) - v_ct(jc,jb))*sinphi(iphi,jc,jb)
          sp_term1k(:) = sp_hscale2o8 + sp_tscale2o8*(c_intr_ct_s(1:nm) - sp_cqh)**2
          sp_term1o(:) = sp_term1k(:)/c_intr_ct_s(1:nm)**2
          sp_termc(:) = mfs_ct(:,iphi,jc,jb)*dm(:,jc,jb)
          DO jk = 0, nkh_bin
            sp_intg_kh(jk) = SUM(sp_termc(:)*EXP(-sp_kh2_bin_edge(jk)*sp_term1k(:)))
          ENDDO
          DO jk = 1, MIN(nkh_bin, nlevp1)
            field% ctmfl_spec0_cgw(jc,jk,jb) = field% ctmfl_spec0_cgw(jc,jk,jb)  &
              &    + (sp_intg_kh(jk-1) - sp_intg_kh(jk))*sp_dphi_o_dkh(jk)
          ENDDO
          DO jk = 0, noi_bin
            sp_intg_oi(jk) = SUM(sp_termc(:)*EXP(-sp_oi2_bin_edge(jk)*sp_term1o(:)))
          ENDDO
          DO jk = 1, MIN(noi_bin, nlevp1-nkh_bin)
            field%        ctmfl_spec0_cgw(jc,nkh_bin+jk,jb)  &
              &  = field% ctmfl_spec0_cgw(jc,nkh_bin+jk,jb)  &
              &    + (sp_intg_oi(jk-1) - sp_intg_oi(jk))*sp_dphi_o_doi(jk)
          ENDDO
        ENDDO
        ENDDO
      ENDDO

    END IF
    IF ( diag_code_100 == 2 .OR. diag_code_100 == 3 ) THEN

      DO im = 1, nm
        sp_kh_bin_edge2(:,im) = sp_oi_bin_edge(:)/c_intr_ct_s(im)
      ENDDO

      DO jb = i_startblk, i_endblk
        DO jc = i_startidx(jb), i_endidx(jb)
        DO iphi = 1, nphi
          sp_termc(:) = mfs_ct(:,iphi,jc,jb)*dm(:,jc,jb)
          DO im = 1, nm
          DO jk = 1, nkh_bin
            IF ( kh_ct(im,iphi,jc,jb) >= sp_kh_bin_edge(jk-1) .AND.  &
              &  kh_ct(im,iphi,jc,jb) <  sp_kh_bin_edge(jk  ) ) THEN
              field%        ctmfl_spec1_cgw(jc,jk,jb)  &
                &  = field% ctmfl_spec1_cgw(jc,jk,jb)  &
                &    + sp_termc(im)*sp_dphi_o_dkh(jk)
              EXIT
            END IF
          ENDDO
          ENDDO
          DO im = 1, nm
          DO jk = 1, noi_bin
            IF ( kh_ct(im,iphi,jc,jb) >= sp_kh_bin_edge2(jk-1,im) .AND.  &
              &  kh_ct(im,iphi,jc,jb) <  sp_kh_bin_edge2(jk  ,im) ) THEN
              field%        ctmfl_spec1_cgw(jc,nkh_bin+jk,jb)  &
                &  = field% ctmfl_spec1_cgw(jc,nkh_bin+jk,jb)  &
                &    + sp_termc(im)*sp_dphi_o_doi(jk)
              EXIT
            END IF
          ENDDO
          ENDDO
        ENDDO
        ENDDO
      ENDDO

    END IF

  END SUBROUTINE diag_spectrum

END MODULE mo_gw_source_conv
