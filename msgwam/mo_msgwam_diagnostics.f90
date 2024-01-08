!>
!! disgnostics module for for MS-GWaM
!!
!! @author Georg Sebastian Voelker, Goethe Uni Frankfurt
!!
!!
!! @par Revision History
!! Initial revision by Georg Sebastian Voelker, Goethe Uni Frankfurt (2016-08-05)
!!
!! @par Copyright and License
!!
!! This code is subject to the Creative Commons Licnese, version 4.0 with Attribution
!! extension (CC-BY 4.0). It is distributed as is without any warranty. 
!!
MODULE mo_msgwam_diagnostics

    ! use msgwam / ICON libraries 
    USE mo_kind,                  ONLY: wp
    USE mo_parallel_config,       ONLY: nproma
    USE mo_grid_config,           ONLY: grid_sphere_radius, is_plane_torus
    USE mo_exception,             ONLY: message, message_text, debug_messages_off, debug_messages_on
    USE mo_run_config,            ONLY: msg_level
    USE mo_physical_constants,    ONLY: grav
    USE mo_timer
    USE mo_msgwam_config
    USE mo_msgwam_util,           ONLY: smooth_vert
    USE mo_setup_msgwam_interface
  
    IMPLICIT NONE
  
    ! declaration of private variables
    PRIVATE
  
    ! declaration of public variables
    PUBLIC :: project_action
    PUBLIC :: project_diagnostics
  
    ! concretization of variables
    INTEGER :: msgwam__diagnostics(10)


CONTAINS
!!
!!-------------------------------------------------------------------------
!!
! SUBROUTINE setup_msgwam_diagnostics()
!     
! END SUBROUTINE setup_msgwam_diagnostics
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE project_diagnostics(nlev, i_startidx, i_endidx, jray_start, jray_end,    &   ! loop indices
    z, zhalf, cellarea,                                                             &   ! grid variables
    jkmin_full, jkmax_full, jkmin_half, jkmax_half, iexist, jk_active, specid,      &   ! ray indices and flags
    zray, dzray, coslatray, dlatray, dlonray,                                       &   ! ray position and extent
    kray, dkray, lray, dlray, mray, dmray, dens,                                    &   ! phase space position and action desity
    bvf2_full, bvf2_half, gammash2_full, gammash2_half, fc, fc2, u, v, theta,       &   ! gridded input
    pmflux_e, pmflux_w, pmflux_n, pmflux_s, apmflux,                                &   ! pseudo-momentum fluxes
    mflux_e, mflux_w, mflux_n, mflux_s, amflux,                                     &   ! momentum-fluxes
    waflux_u, waflux_d, waflux_e, waflux_w, waflux_s, waflux_n,                     &   ! wave action fluxes
    ptflux_e, ptflux_w, ptflux_s, ptflux_n, aptflux,                                &   ! potential tempertature fluxes
    energy, energy_p, waction, active_rays                                          )   ! energies and action

    INTEGER,              INTENT(IN)  :: nlev
    INTEGER,              INTENT(IN)  :: i_startidx ! first index of the block
    INTEGER,              INTENT(IN)  :: i_endidx   ! last index of the block
    INTEGER,              INTENT(IN)  :: jray_start
    INTEGER,              INTENT(IN)  :: jray_end
    INTEGER,              INTENT(IN)  :: iexist(:,:)
    INTEGER,              INTENT(IN)  :: specid(:,:)
    INTEGER,              INTENT(IN)  :: jk_active(:,:)
    INTEGER,              INTENT(IN)  :: jkmin_full(:,:)
    INTEGER,              INTENT(IN)  :: jkmax_full(:,:)
    INTEGER,              INTENT(IN)  :: jkmin_half(:,:)
    INTEGER,              INTENT(IN)  :: jkmax_half(:,:)
    REAL(wp),             INTENT(IN)  :: z(:,:), zhalf(:,:)
    REAL(wp),             INTENT(IN)  :: fc(:), cellarea(:)   ! area of grid cell [m**2]
    REAL(wp),             INTENT(IN)  :: zray(:,:), dzray(:,:)
    REAL(wp),             INTENT(IN)  :: coslatray(:, :), dlatray(:, :), dlonray(:, :)
    REAL(wp),             INTENT(IN)  :: mray(:,:), dmray(:,:)
    REAL(wp),             INTENT(IN)  :: kray(:,:), dkray(:,:)
    REAL(wp),             INTENT(IN)  :: lray(:,:), dlray(:,:)
    REAL(wp),             INTENT(IN)  :: dens(:,:)
    REAL(wp),             INTENT(IN)  :: theta(:,:)
    REAL(wp),             INTENT(IN)  :: fc2(:)
    REAL(wp),             INTENT(IN)  :: bvf2_full(:,:)
    REAL(wp),             INTENT(IN)  :: bvf2_half(:,:)
    REAL(wp),             INTENT(IN)  :: gammash2_full(:,:)
    REAL(wp),             INTENT(IN)  :: gammash2_half(:,:)
    REAL(wp),             INTENT(IN)  :: u(:, :)
    REAL(wp),             INTENT(IN)  :: v(:, :)
    REAL(wp),             INTENT(OUT) :: pmflux_e(:,:)   
    REAL(wp),             INTENT(OUT) :: pmflux_w(:,:) 
    REAL(wp),             INTENT(OUT) :: pmflux_s(:,:) 
    REAL(wp),             INTENT(OUT) :: pmflux_n(:,:) 
    REAL(wp),             INTENT(OUT) :: apmflux(:,:)   
    REAL(wp),             INTENT(OUT) :: mflux_e(:,:)   
    REAL(wp),             INTENT(OUT) :: mflux_w(:,:) 
    REAL(wp),             INTENT(OUT) :: mflux_s(:,:) 
    REAL(wp),             INTENT(OUT) :: mflux_n(:,:) 
    REAL(wp),             INTENT(OUT) :: amflux (:,:) 
    REAL(wp),             INTENT(OUT) :: waflux_u(:,:)
    REAL(wp),             INTENT(OUT) :: waflux_d(:,:)
    REAL(wp),             INTENT(OUT) :: waflux_e(:,:)
    REAL(wp),             INTENT(OUT) :: waflux_w(:,:)
    REAL(wp),             INTENT(OUT) :: waflux_s(:,:)
    REAL(wp),             INTENT(OUT) :: waflux_n(:,:)
    REAL(wp),             INTENT(OUT) :: ptflux_e(:,:)
    REAL(wp),             INTENT(OUT) :: ptflux_w(:,:)
    REAL(wp),             INTENT(OUT) :: ptflux_s(:,:)
    REAL(wp),             INTENT(OUT) :: ptflux_n(:,:)
    REAL(wp),             INTENT(OUT) :: aptflux (:,:)
    REAL(wp),             INTENT(OUT) :: energy(:,:)  
    REAL(wp),             INTENT(OUT) :: energy_p(:,:)
    REAL(wp),             INTENT(OUT) :: waction(:,:) 
    REAL(wp),             INTENT(OUT) :: active_rays(:)

    INTEGER                      :: nlevp1
    INTEGER                      :: jk, jc
    INTEGER                      :: jray
    REAL(wp)                     :: zray_u, zray_d_eff
    REAL(wp)                     :: dzi_o_dz, r_sq_o_a, r_sq
    REAL(wp)                     :: K2, Kh2, Kh, m2, K2_p_gam2  ! Local wave vector squares
    REAL(wp)                     :: omega, omega2               ! Local intrinsic freqs
    REAL(wp)                     :: cglon, cglat, cgr, cgr_mod  ! Local group velocities
    REAL(wp)                     :: a_ray, a_ray_jk             ! local wave actions
    REAL(wp)                     :: m2_p_gam2(nlev+1)
    REAL(wp)                     :: e_tot_ray(nlev)
    REAL(wp)                     :: f2_o_bvf2(nlev)
    REAL(wp)                     :: fact_grid(nproma,nlev)

    IF (msg_level >= 12) CALL message('project_diagnostics', 'MS-GWaM: diagnostics on full levels')

    !----------------------------------------------------------------------
    ! Purpose:
    !         Project GW fluxes and other quantities carried by
    !         Lagrangian ray volumes to the Eulerian grid.
    !
    ! Method:
    !         -- loop over ray volumes
    !         -- calculate fractional flux (energy, etc.) corresponding to
    !            each vertical layer they cross over (jkmin, jkmax known
    !            from idx_rayedge)
    !         -- add up these fractional fluxes (energy, etc.) for each
    !            vertical layer
    !         -- all horizontal flux diagnostics and scalar diagnostics are
    !            calculated on full levels; vertical fluxes on half levels
    !         -- see Bölöni et al. (2021); Voelker et al. (2023)
    !
    ! TODO:
    !----------------------------------------------------------------------

    nlevp1 = nlev+1

    ! Initialize outputs
    apmflux(:, :) = 0._wp ; amflux(:, :) = 0._wp ; aptflux(:, :) = 0._wp ;
    energy(:, :) = 0._wp  ; energy_p(:, :) = 0._wp ; waction(:, :) = 0._wp
    active_rays(:) = 0._wp

    ! Initialize fluxes in E, W, N, S directions
    mflux_e(:,:) = 0._wp  ; mflux_w(:,:) = 0._wp
    mflux_n(:,:) = 0._wp  ; mflux_s(:,:) = 0._wp
    pmflux_e(:,:) = 0._wp ; pmflux_w(:,:) = 0._wp
    pmflux_n(:,:) = 0._wp ; pmflux_s(:,:) = 0._wp
    waflux_d(:,:) = 0._wp ; waflux_u(:,:) = 0._wp
    waflux_e(:,:) = 0._wp ; waflux_w(:,:) = 0._wp
    waflux_n(:,:) = 0._wp ; waflux_s(:,:) = 0._wp
    ptflux_e(:,:) = 0._wp ; ptflux_w(:,:) = 0._wp
    ptflux_n(:,:) = 0._wp ; ptflux_s(:,:) = 0._wp

    ! Initialize temporary variables
    e_tot_ray(:) = 0._wp ; f2_o_bvf2(:) = 0._wp

    IF (is_plane_torus) THEN
        r_sq = 1._wp
    ELSE
        r_sq = grid_sphere_radius**2
    END IF

    DO jc = i_startidx, i_endidx

        f2_o_bvf2(:) = fc2(jc) / bvf2_full(jc,:)    ! f^2/N^2
        r_sq_o_a = r_sq / cellarea(jc)              ! fractional area

        DO jray = jray_start, jray_end

            IF (iexist(jc,jray) == 0) CYCLE

            active_rays(jc) = active_rays(jc) + 1._wp

            ! Define (effective) ray volume edges
            zray_u = zray(jc, jray) + 0.5_wp*dzray(jc, jray)
            zray_d_eff = zray(jc, jray) - 0.5_wp*dzray(jc, jray)
            IF (specid(jc, jray) < 0) THEN
                IF (zray_u < z(jc, jk_active(jc, jray)))  CYCLE
                zray_d_eff = MAX(z(jc, jk_active(jc, jray)), zray_d_eff)
            ENDIF

            Kh2 = kray(jc, jray)**2 + lray(jc, jray)**2 ! Effective horizontal wavenumber squared
            Kh = SQRT(Kh2)                              ! Effective horizontal wavenumber
            m2 = mray(jc, jray)**2                      ! Vertical wavenumber squared
            K2 = Kh2 + m2                               ! Effective total wavenumber squared

            ! Wave action density (integrated over wavenumber space)
            ! get the fraction of horizontal area that the ray occupies in the cell for accumulating its effect
            ! Here the formular of fractional area should be consistent with other parts
            ! (e.g. source formular, arealonk formular) for conservation and smoothness

            a_ray = dens(jc,jray) * dkray(jc,jray) * dlray(jc,jray) * dmray(jc,jray)  &
                &    * dlatray(jc,jray) * (dlonray(jc,jray) * coslatray(jc,jray)) * r_sq_o_a

            ! calculate horizontal diagnostic fluxes and scalar fields on full levels 
            DO jk = jkmin_half(jc,jray), jkmax_half(jc,jray)

                ! calculate locally corrected wave number squares
                m2_p_gam2(jk) = m2 + gammash2_full(jc, jk)
                K2_p_gam2 = K2 + gammash2_full(jc, jk)

                ! calculate the local intrinsic frequency
                omega2 = (bvf2_full(jc, jk) * Kh2 + fc2(jc) * m2_p_gam2(jk)) / K2_p_gam2
                omega = SQRT(omega2)

                ! calculate the corresponding local intrinsic group velocities
                cglon = branch * kray(jc, jray) / (SQRT(omega2) * K2_p_gam2) * (bvf2_full(jc, jk) - omega2) + u(jc, jk)
                cglat = branch * lray(jc, jray) / (SQRT(omega2) * K2_p_gam2) * (bvf2_full(jc, jk) - omega2) + v(jc, jk)
                cgr = - branch * mray(jc, jray) / (SQRT(omega2) * K2_p_gam2) * (omega2 - fc2(jc))

                ! get partial height of ray volume in layer
                dzi_o_dz = (MIN(zray_u, z(jc, jk-1)) - MAX(zray_d_eff, z(jc, jk)))  / (z(jc, jk-1) - z(jc, jk))

                ! calculate partial wave action associated to the vertical layer
                a_ray_jk = ABS(dzi_o_dz * a_ray)

                ! calculate horizontal wave action fluxes
                waflux_w(jc, jk) = waflux_w(jc, jk) + MIN(0._wp, cglon) * a_ray_jk      ! westward wave action flux
                waflux_e(jc, jk) = waflux_e(jc, jk) + MAX(0._wp, cglon) * a_ray_jk      ! eastward wave action flux
                waflux_s(jc, jk) = waflux_s(jc, jk) + MIN(0._wp, cglat) * a_ray_jk      ! southward wave action flux
                waflux_n(jc, jk) = waflux_n(jc, jk) + MAX(0._wp, cglat) * a_ray_jk      ! northward wave action flux

                ! calculate potential temperature fluxes
                ! grid related parameters are multiplied later
                ptflux_w(jc, jk) = ptflux_w(jc, jk) + MIN(0._wp, -lray(jc, jray) * cgr) * a_ray_jk / (omega2 - fc2(jc))   ! westward pt flux
                ptflux_e(jc, jk) = ptflux_e(jc, jk) + MAX(0._wp, -lray(jc, jray) * cgr) * a_ray_jk / (omega2 - fc2(jc))   ! eastward pt flux
                ptflux_s(jc, jk) = ptflux_s(jc, jk) + MIN(0._wp,  kray(jc, jray) * cgr) * a_ray_jk / (omega2 - fc2(jc))   ! southward pt flux
                ptflux_n(jc, jk) = ptflux_n(jc, jk) + MAX(0._wp,  kray(jc, jray) * cgr) * a_ray_jk / (omega2 - fc2(jc))   ! northward pt flux
                
                aptflux(jc, jk) = aptflux(jc, jk) + SQRT(           &
                    &     (ptflux_e(jc, jk) + ptflux_w(jc, jk))**2  &
                    &   + (ptflux_n(jc, jk) + ptflux_s(jc, jk))**2)

                ! calculate the energies and wave actions
                e_tot_ray(jk) = omega * a_ray_jk
                energy(jc,jk) = energy(jc,jk) + e_tot_ray(jk)
                waction(jc,jk) = waction(jc,jk) + a_ray_jk
                energy_p(jc,jk) = energy_p(jc,jk) + e_tot_ray(jk) / (1._wp + f2_o_bvf2(jk) * m2_p_gam2(jk) / Kh2)
            
            ENDDO ! jk

            ! calculate vertical fluxes on half levels 
            DO jk = jkmin_full(jc,jray), jkmax_full(jc,jray)

                ! calculate locally corrected wave number squares
                m2_p_gam2(jk) = m2 + gammash2_half(jc, jk)
                K2_p_gam2 = K2 + gammash2_half(jc, jk)

                ! calculate the local intrinsic frequency
                omega2 = (bvf2_half(jc, jk) * Kh2 + fc2(jc) * m2_p_gam2(jk)) / K2_p_gam2
                omega = SQRT(omega2)

                ! calculate the corresponding local intrinsic group velocities
                cgr = - branch * mray(jc, jray) / (SQRT(omega2) * K2_p_gam2) * (omega2 - fc2(jc))
                cgr_mod = cgr * omega2 / (omega2 + fc2(jc))

                ! get partial height of ray volume in layer
                dzi_o_dz = (MIN(zray_u, z(jc, jk-1)) - MAX(zray_d_eff, z(jc, jk)))  / (z(jc, jk-1) - z(jc, jk))

                ! calculate partial wave action associated to the vertical layer
                a_ray_jk = ABS(dzi_o_dz * a_ray)

                ! calculate vertical fluxes of directed horizontal pseodo-momentum
                pmflux_w(jc, jk) = pmflux_w(jc, jk) + MIN(0._wp, cgr * kray(jc, jray)) * a_ray_jk       ! vertical flux of westward pseudo-momentum
                pmflux_e(jc, jk) = pmflux_e(jc, jk) + MAX(0._wp, cgr * kray(jc, jray)) * a_ray_jk       ! vertical flux of eastward pseudo-momentum
                pmflux_s(jc, jk) = pmflux_s(jc, jk) + MIN(0._wp, cgr * lray(jc, jray)) * a_ray_jk       ! vertical flux of southward pseudo-momentum
                pmflux_n(jc, jk) = pmflux_n(jc, jk) + MAX(0._wp, cgr * lray(jc, jray)) * a_ray_jk       ! vertical flux of northward pseudo-momentum

                ! calculate absolute pseudo-momentum flux
                apmflux(jc, jk) = apmflux(jc, jk) + SQRT(           &
                    &     (pmflux_e(jc, jk) + pmflux_w(jc, jk))**2  &
                    &   + (pmflux_n(jc, jk) + pmflux_s(jc, jk))**2)

                ! calculate vertical fluxes of directed horizontal momentum
                mflux_w(jc, jk) = mflux_w(jc, jk) + MIN(0._wp, cgr_mod * kray(jc, jray)) * a_ray_jk     ! vertical flux of westward momentum
                mflux_e(jc, jk) = mflux_e(jc, jk) + MAX(0._wp, cgr_mod * kray(jc, jray)) * a_ray_jk     ! vertical flux of eastward momentum
                mflux_s(jc, jk) = mflux_s(jc, jk) + MIN(0._wp, cgr_mod * lray(jc, jray)) * a_ray_jk     ! vertical flux of southward momentum
                mflux_n(jc, jk) = mflux_n(jc, jk) + MAX(0._wp, cgr_mod * lray(jc, jray)) * a_ray_jk     ! vertical flux of northward momentum

                ! calculate absolute momentum flux
                amflux(jc, jk) = amflux(jc, jk) + SQRT(           &
                    &     (mflux_e(jc, jk) + mflux_w(jc, jk))**2  &
                    &   + (mflux_n(jc, jk) + mflux_s(jc, jk))**2)

                ! calculate wave action fluxes
                waflux_d(jc, jk) = waflux_d(jc, jk) + MIN(0._wp, cgr) * a_ray_jk                        ! downward wave action flux
                waflux_u(jc, jk) = waflux_u(jc, jk) + MAX(0._wp, cgr) * a_ray_jk                        ! upward wave action flux
            
            ENDDO ! jk

        ENDDO ! jray

    ENDDO ! jc

    ! Finalize heat fluxes
    DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
            ! calculate gridded prefactor
            fact_grid(jc, jk) = fc(jc) * theta(jc, jk) * bvf2_full(jc, jk) / grav
            
            ! multiply the factor into the fluxes
            ptflux_w(jc, jk) = ptflux_w(jc, jk) * fact_grid(jc, jk)
            ptflux_e(jc, jk) = ptflux_e(jc, jk) * fact_grid(jc, jk) 
            ptflux_s(jc, jk) = ptflux_s(jc, jk) * fact_grid(jc, jk)
            ptflux_n(jc, jk) = ptflux_n(jc, jk) * fact_grid(jc, jk)
        ENDDO
    ENDDO

    ! Finalize potential energy calculation (use energy equipartition)
    DO jk = 2, nlev-1
        DO jc = i_startidx, i_endidx
            energy_p(jc,jk) = energy_p(jc,jk) * 0.5_wp
        ENDDO
    ENDDO

    ! Smooth fluxes / action / energies with a Shapiro filter
    IF (lsmooth) THEN
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., energy)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., energy_p)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., waction)
        
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., ptflux_e)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., ptflux_w)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., ptflux_n)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., ptflux_s)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., aptflux)
        
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., waflux_e)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., waflux_w)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., waflux_n)
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., waflux_s)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., waflux_u)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., waflux_d)
        
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., pmflux_e)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., pmflux_w)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., pmflux_n)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., pmflux_s)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., apmflux)
        
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., mflux_e)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., mflux_w)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., mflux_n)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., mflux_s)
        CALL smooth_vert(nlevp1, i_startidx, i_endidx, nsmooth, .FALSE., amflux)
    END IF


    END SUBROUTINE project_diagnostics

    SUBROUTINE project_action(nlev, i_startidx, i_endidx, jray_start, jray_end, z, zhalf, &
        cellarea, zray, dzray, coslatray, dlatray, dlonray, dkray, dlray, dmray, dens, iexist, &
        specid, jk_active, jkmin_full, jkmax_full, jkmin_half, jkmax_half, action, diag, jg)
        
    INTEGER,              INTENT(IN)  :: nlev
    INTEGER,              INTENT(IN)  :: i_startidx ! first index of the block
    INTEGER,              INTENT(IN)  :: i_endidx   ! last index of the block
    INTEGER,              INTENT(IN)  :: jray_start
    INTEGER,              INTENT(IN)  :: jray_end
    INTEGER,              INTENT(IN)  :: iexist(:,:)
    INTEGER,              INTENT(IN)  :: specid(:,:)
    INTEGER,              INTENT(IN)  :: jk_active(:,:)
    INTEGER,              INTENT(IN)  :: jkmin_full(:,:)
    INTEGER,              INTENT(IN)  :: jkmax_full(:,:)
    INTEGER,              INTENT(IN)  :: jkmin_half(:,:)
    INTEGER,              INTENT(IN)  :: jkmax_half(:,:)
    INTEGER,              INTENT(IN)  :: diag, jg
    REAL(wp),             INTENT(IN)  :: z(:,:), zhalf(:,:)
    REAL(wp),             INTENT(IN)  :: cellarea(:)   ! area of grid cell [m**2]
    REAL(wp),             INTENT(IN)  :: zray(:,:), dzray(:,:)
    REAL(wp),             INTENT(IN)  :: coslatray(:, :), dlatray(:, :), dlonray(:, :)
    REAL(wp),             INTENT(IN)  :: dmray(:,:), dkray(:,:), dlray(:,:)
    REAL(wp),             INTENT(IN)  :: dens(:,:)
    
    
    REAL(wp),             INTENT(OUT) :: action(:,:)
    
    INTEGER                      :: jk, jc
    INTEGER                      :: jray
    REAL(wp)                     :: zray_u
    REAL(wp)                     :: zray_d_eff
    REAL(wp)                     :: dzi_o_dz, r_sq_o_a, r_sq
    REAL(wp)                     :: a_ray
    REAL(wp)                     :: a_ray_jk
    
    IF (msg_level >= 12) CALL message('project_action', 'MS-GWaM: wave action diagnostic')
    
    !----------------------------------------------------------------------
    ! Purpose: 
    !         Project GW wave action carried by Lagrangian ray volumes to
    !         the Eulerian grid.
    !
    ! Method:
    !         -- loop over ray volumes 
    !         -- calculate fractional action corresponding to 
    !            each vertical layer they cross over (jkmin, jkmax known 
    !            from idx_rayedge)
    !         -- add up these fractional action for each vertical layer
    !----------------------------------------------------------------------
    
    ! Initialize outputs
    action(:,:) = 0._wp
    
    IF (is_plane_torus) THEN
        r_sq = 1._wp
    ELSE
        r_sq = grid_sphere_radius**2
    END IF
    
    CALL debug_messages_on
    
    DO jc = i_startidx, i_endidx
    
        r_sq_o_a = r_sq / cellarea(jc)
    
        DO jray = jray_start, jray_end
    
        IF (iexist(jc, jray) == 0) CYCLE
        IF ((jkmin_half(jc, jray) < 1) .OR. (jkmax_half(jc, jray) > nlev)) cycle
    
        ! Define (effective) ray volume edges
        zray_u = zray(jc, jray) + 0.5_wp*dzray(jc, jray)
        zray_d_eff = zray(jc, jray) - 0.5_wp*dzray(jc, jray)
        IF (specid(jc,jray) < 0) THEN
            IF (zray_u < z(jc,jk_active(jc,jray)))  CYCLE
            zray_d_eff = MAX(z(jc,jk_active(jc,jray)),zray_d_eff)
        ENDIF
    
        ! Wave action density (integrated over wavenumber space)
        ! get the fraction of horizontal area that the ray occupies in the cell for accumulating its effect
        a_ray = dkray(jc,jray)*dlray(jc,jray)*dmray(jc,jray)*dens(jc,jray)  &
            &    *dlatray(jc,jray)*(dlonray(jc,jray)*coslatray(jc,jray))*r_sq_o_a
    
        if (diag > 0 .AND. jg == 1) then
            if (jc == 5) then
            WRITE(message_text,'(a,i4,2E12.4)') 'ray diag: jray, dens(jc=5, jray), a_ray', &
                    & jray, dens(jc, jray), a_ray
            CALL message('', TRIM(message_text))
            endif
        endif
        
        ! Calculate wave action on full levels
        DO jk = jkmin_half(jc,jray), jkmax_half(jc,jray) ! jkmin_half > 1, jkmax_half < nlevp1
            ! Note: jkmin_half and jkmax_half denote the half level index above the
            ! ray volume top/bottom
    
            ! dzi: size of ray in z-direction corresponding to the actual full layer
            ! dz : full layer depth
            dzi_o_dz = (MIN(zray_u, zhalf(jc, jk)) - MAX(zray_d_eff, zhalf(jc, jk+1)))  &
                / (zhalf(jc, jk) - zhalf(jc, jk+1))
                
            ! Wave action density (integrated over wavenumber space) relevant to the
            ! vertical level with index jk
            a_ray_jk = dzi_o_dz * a_ray
            action(jc,jk) = action(jc,jk) + a_ray_jk
        ENDDO ! jk
    
        ENDDO ! jray
    
    ENDDO ! jc
    
    CALL debug_messages_off
    
    ! smooth wave action
    IF (lsmooth) THEN
        CALL smooth_vert(nlev, i_startidx, i_endidx, nsmooth, .FALSE., action)
    END IF
  
  END SUBROUTINE project_action
  
END MODULE mo_msgwam_diagnostics
  