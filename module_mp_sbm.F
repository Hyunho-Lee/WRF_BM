module module_mp_sbm
!-------------------------------------------------------------------------!
!                                                                         !
!  Author : Alexander P. Khain and many coauthors                         !
!    Hebrew University of Jerusalem, Israel                               !
!                                                                         !
!  Coupling to WRF : Hyunho Lee                                           !
!    Seoul National University, South Korea                               !
!                                                                         !
!-------------------------------------------------------------------------!

integer, private, save :: inan
integer, private, save :: iremap       ! 1 for KO, 2 for new
integer, private, save :: ibreakup     ! 1 for implicit, 2 for explicit
integer, private, save :: iturb        ! 1 on, 2 off
integer, private, save :: icollision   ! 1 for NQS, 2 for Poisson
real,    private, save :: n0ccn1       ! N0 on the land [# m-3]
real,    private, save :: kccn1        ! k on the land
real,    private, save :: n0ccn2       ! N0 on the ocean [# m-3]
real,    private, save :: kccn2        ! k on the ocean
real,    private, save :: urb_ccnf     ! Moltiplier for N0 on the urban

real, parameter, private :: tol = 5.0e-3

integer, parameter, private :: nkr        = 43
integer, parameter, private :: kr_droplet = 14 ! cloud droplet: less than 40 micro m
integer, parameter, private :: kr_big_ice = 29 ! type conversion
integer, parameter, private :: kr_snow    = 23 ! type conversion
integer, parameter, private :: kr_graupel = 34 ! type conversion
integer, parameter, private :: kr_snow1   = 31 ! snow breakup
integer, parameter, private :: kr_snow2   = 34 ! snow breakup
real,    parameter, private :: lwf_upper_limit = 0.5
real,    parameter, private :: rf_upper_limit = 0.95
real,    parameter, private :: rho_upper_limit = 200.

integer, parameter, private :: nhydro   = 7
integer, parameter, private :: WATER    = 1
integer, parameter, private :: COLUMN   = 2
integer, parameter, private :: PLATE    = 3
integer, parameter, private :: DENDRITE = 4
integer, parameter, private :: SNOW     = 5
integer, parameter, private :: GRAUPEL  = 6
integer, parameter, private :: HAIL     = 7

integer, parameter, private :: scal = 1
real,    parameter, private :: dlnr = log(2.0)/3.0/scal

real, dimension(nkr),        save, private :: mass_hydro
real, dimension(nhydro,nkr), save, private :: vt_hydro
real, dimension(nhydro,nkr), save, private :: rho_hydro
real, dimension(nhydro,nkr), save, private :: r_hydro
real, dimension(nhydro,nkr), save, private :: cap_hydro
real, dimension(nhydro,nkr), save, private :: bulkrho
real, dimension(nhydro,nkr), save, private :: bulkvt
real, dimension(nhydro,nkr), save, private :: vt_ref
real, dimension(nkr),        save, private :: rho_ice_snow
real, dimension(nhydro,nhydro,  nkr,nkr), save, private :: ckxx
real, dimension(nhydro,nhydro,3,nkr,nkr), save, private :: ywxx

real, parameter, private :: rccn_min = 0.003e-4 ! cm
real, parameter, private :: smax     = 0.20     ! 120% of RH
real, parameter, private :: z0_ccn   = 2000.    ! 2.0 km
real, parameter, private :: dz_ccn   = 2000.    ! 2.0 km

real, parameter, private :: t00      = 273.15
real, parameter, private :: p00      = 101325.
real, parameter, private :: p1000hpa = 100000.

real, parameter, private :: Rd    = 287.
real, parameter, private :: Rv    = 461.51
real, parameter, private :: eps   = Rd/Rv
real, parameter, private :: cp    = 1004.7
real, parameter, private :: Ki    = 2.25        ! heat conductivity of ice
real, parameter, private :: Ls    = 2.8345e6    ! latent heat for sublimation
real, parameter, private :: Lm    = 0.3337e6    ! latent heat for melting
real, parameter, private :: Lv    = 2.5008e6    ! latent heat for evaporation
real, parameter, private :: difv   = 2.11e-5    ! diffusivity of vapor
real, parameter, private :: thca   = 4.38e-3    ! thermal condunctivity of air
real, parameter, private :: visa   = 1.827e-5   ! dynamic viscosity of air

real, parameter, private :: b_w = Lv/Rv
real, parameter, private :: a_w = b_w/t00
real, parameter, private :: b_i = Ls/Rv
real, parameter, private :: a_i = b_i/t00

real, parameter, private :: pi    = 3.141593
real, parameter, private :: grav  = 9.81

real, parameter, private :: rho_water = 1000.
real, parameter, private :: rho_ice   =  920.
real, parameter, private :: rho_rim   =  900.
real, parameter, private :: rho_crit  =  910.

real,    parameter, private :: tt_coal = t00 - 40.0
integer, dimension(nkr,nkr), save, private :: ima
real,    dimension(nkr,nkr), save, private :: cour
integer, parameter, private :: num_coll_type_lt0 = 17
integer, parameter, private :: num_coll_type_ge0 = 5
integer, dimension(num_coll_type_lt0+num_coll_type_ge0,4), private :: coll_table

integer, parameter, private :: kr_icempl = 9
real,    parameter, private :: t_icempl0 = t00 - 5.0
real,    parameter, private :: t_icempl1 = t00 - 8.0
real,    parameter, private :: t_icempl2 = t00 - 3.0

integer, parameter, private :: krbreak = 28
real, dimension(krbreak,krbreak,krbreak), save, private :: pkij
real, dimension(krbreak,krbreak),         save, private :: qkj

real, parameter, private :: ventpl_max = 5.0

real, parameter, private :: ro_solute = 2.16
real, parameter, private :: akoe = 3.3e-05
real, parameter, private :: bkoe = 2.0*4.3/(22.9+35.5)*(4.0/3.0)*pi*ro_solute

real, parameter, private :: a1_mey = -0.639
real, parameter, private :: b1_mey = 0.1296
real, parameter, private :: c1_mey = 1.0e+3
real, parameter, private :: a2_mey = -2.8
real, parameter, private :: b2_mey = 0.262
real, parameter, private :: c2_mey = 0.0
real, parameter, private :: supi_max = 30.0
real, parameter, private :: dsupi_max = 1.0/b1_mey

real, parameter, private :: tt_nucl_drop_min = -45.0
real, parameter, private :: tt_nucl_ice_min = -45.0

real, parameter, private :: tau_ccn = 3600.0 ! relaxation time is i h

real, parameter, private :: dt_cond_ref = 1.0
real, save,      private :: dt_cond
integer, save,   private :: n_cond

real, parameter :: a_frz  = 1.0e-1
real, parameter :: b_frz  = 0.66
real, parameter :: n0_frz = 1.0e+7 ! parameters an0 and gamma obtained from (Vali,1975) n0 = 10**7/m^3
real, parameter :: gam    = 4.4
real, parameter :: mass_frz_bnd = 2.0e-9

integer, parameter, private :: krminl_gl=1
integer, parameter, private :: krming_gl=1
integer, parameter, private :: krmaxl_gl=24
integer, parameter, private :: krmaxg_gl=33
real, save, private :: turbgl_grid_400(krming_gl:krmaxg_gl,krminl_gl:krmaxl_gl)
real, save, private :: turbgl_grid_800(krming_gl:krmaxg_gl,krminl_gl:krmaxl_gl)
real, save, private :: turbll_grid(13,13,11)

contains
!--------------------------------------------------------------------------
subroutine sbm (u_phy, v_phy, w_phy, th_phy, qv, pi_phy, p_phy, &
                dx, dy, dz8w,                                   &
                tke_opt, tke3d, n2_3d,                          &
                xland, ivgtyp, isurban,                         &
                scalar, n_scalar,                               &
                qc, qr, qi, qs, qg, qh,                         &
                nc, nr, ni, ns, ng, nh,  na,                    &
                rainnc,    rainncv,    snownc, snowncv,         &
                graupelnc, graupelncv, hailnc, hailncv,         &
                dTfm3d, dTc3d, dTn3d, dTce3d, dTds3d,           &
                eps_phy, rey_phy,                               &
                qr_au, qr_ac,                                   &
                its, ite, jts, jte, kts, kte,                   &
                ims, ime, jms, jme, kms, kme,                   &
                ids, ide, jds, jde, kds, kde,                   &
                p_start, dt, itimestep                          )
!--------------------------------------------------------------------------
use ieee_arithmetic
implicit none

integer, intent(in) ::              &
      its, ite, jts, jte, kts, kte, &
      ims, ime, jms, jme, kms, kme, &
      ids, ide, jds, jde, kds, kde, &
      n_scalar, isurban, tke_opt
integer, intent(in) :: p_start, itimestep
real, intent(in) :: dx, dy, dt
real, dimension(ims:ime,kms:kme,jms:jme), intent(in   ) :: &
      u_phy                                                &
    , v_phy                                                &
    , w_phy                                                &
    , dz8w                                                 &
    , pi_phy                                               &
    , p_phy                                                &
    , tke3d                                                &
    , n2_3d
real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: &
      th_phy                                               &
    , qv
real, dimension(ims:ime,kms:kme,jms:jme), intent(out  ) :: &
      qc                                                   &
    , qr                                                   &
    , qi                                                   &
    , qs                                                   &
    , qg                                                   &
    , qh                                                   &
    , nc                                                   &
    , nr                                                   &
    , ni                                                   &
    , ns                                                   &
    , ng                                                   &
    , nh                                                   &
    , na                                                   &
    , dTc3d                                                &
    , dTn3d                                                &
    , dTfm3d                                               &
    , dTce3d                                               &
    , dTds3d                                               &
    , qr_au                                                &
    , qr_ac                                                &
    , eps_phy                                              &
    , rey_phy
real, dimension(ims:ime,kms:kme,jms:jme,n_scalar), intent(inout) :: &
      scalar
real, dimension(ims:ime,jms:jme), intent(in) :: &
      xland
real, dimension(ims:ime,jms:jme), intent(inout) :: &
      rainnc,  snownc,  graupelnc,  hailnc
real, dimension(ims:ime,jms:jme), intent(out  ) :: &
      rainncv, snowncv, graupelncv, hailncv
integer, dimension(ims:ime,jms:jme), intent(in) :: &
      ivgtyp

! Local variables
real, dimension(ims:ime,kms:kme,jms:jme) :: tt, supi, dtemp, dsupi, hgt
integer :: i, j, k, m, n, ii, kk
real,    dimension(ite-its+1) :: rain1d, snow1d, graupel1d, hail1d, xland1d
integer, dimension(ite-its+1) :: ivgtyp1d
real,    dimension(ite-its+1,kte-kts+1) ::         &
      th2d, qv2d, pp2d, pi2d, hgt2d, dz2d          &
    , qc2d, qr2d, qi2d, qs2d, qg2d, qh2d           &
    , nc2d, nr2d, ni2d, ns2d, ng2d, nh2d, na2d     &
    , dTfm2d, dTc2d, dTn2d, dTce2d, dTds2d         &
    , tke2d, n2_2d, eps2d, rey2d, qr_au2d, qr_ac2d &
    , dtemp2d, dsupi2d
real, dimension(ite-its+1,kte-kts+1,       nkr) :: fccn
real, dimension(ite-its+1,kte-kts+1,nhydro,nkr) :: lw
real, dimension(ite-its+1,kte-kts+1,       nkr) :: rm
real, dimension(ite-its+1,kte-kts+1,nhydro,nkr) :: ffx

! tke option
if(itimestep.eq.1) then
  if(tke_opt.eq.1) then
    call wrf_message("tke is inputed normally.")
  else if(tke_opt.eq.-1) then
    call wrf_message("pbl_physics should predict tke. input tke is zero.")
  else if(tke_opt.eq.-2) then
    call wrf_message("km_opt should be 1.5th TKE option. input tke is zero.")
  else
    call wrf_message("Something wrong in processing.")
  end if
end if

do j=jts, jte
do i=its, ite

  hgt(i,kts,j) = dz8w(i,kts,j)/2.0
  do k=kts+1, kte
    hgt(i,k,j) = hgt(i,k-1,j) + (dz8w(i,k-1,j)+dz8w(i,k,j))/2.0
  end do

end do
end do

! calculate total derivatives
dtemp = 0.
dsupi = 0.

do j=jms, jme
do k=kms, kme
do i=ims, ime
  tt  (i,k,j) = th_phy(i,k,j)*pi_phy(i,k,j)
  qv  (i,k,j) = max(qv(i,k,j),0.)
  supi(i,k,j) = (p_phy(i,k,j)*qv(i,k,j)/(eps+qv(i,k,j)))/esi(tt(i,k,j))-1.0
end do
end do
end do

do j=max(jts,jds+1), min(jte,jde-2)
do k=max(kts,kds+1), min(kte,kde-2)
do i=max(its,ids+1), min(ite,ide-2)
  dtemp(i,k,j) = (tt(i+1,k,j)*u_phy(i+1,k,j)           &
                 -tt(i-1,k,j)*u_phy(i-1,k,j))/(2.0*dx) &
                +(tt(i,k,j+1)*v_phy(i,k,j+1)           &
                 -tt(i,k,j-1)*v_phy(i,k,j-1))/(2.0*dy) &
                +(tt(i,k+1,j)*w_phy(i,k+1,j)           &
                 -tt(i,k-1,j)*w_phy(i,k-1,j))/(hgt(i,k+1,j)-hgt(i,k-1,j))

  dtemp(i,k,j) = dtemp(i,k,j)*dt

  dsupi(i,k,j) = (supi(i+1,k,j)*u_phy(i+1,k,j)     &
                 -supi(i-1,k,j)*u_phy(i-1,k,j))/(2.0*dx) &
                +(supi(i,k,j+1)*v_phy(i,k,j+1)     &
                 -supi(i,k,j-1)*v_phy(i,k,j-1))/(2.0*dy) &
                +(supi(i,k+1,j)*w_phy(i,k+1,j)     &
                 -supi(i,k-1,j)*w_phy(i,k-1,j))/(hgt(i,k+1,j)-hgt(i,k-1,j))

  dsupi(i,k,j) = dsupi(i,k,j)*dt
end do
end do
end do

! From 3-D to 2-D
do j=jts, jte

  do k=kts, kte
  do i=its, ite

    ii = i-its+1
    kk = k-kts+1

    th2d(ii,kk) = th_phy(i,k,j)
    qv2d(ii,kk) = max(qv(i,k,j),0.)
    pp2d(ii,kk) = p_phy(i,k,j)
    pi2d(ii,kk) = pi_phy(i,k,j)
    hgt2d(ii,kk) = hgt(i,k,j)
    dz2d(ii,kk) = dz8w(i,k,j)

    tke2d(ii,kk) = tke3d(i,k,j)
    n2_2d(ii,kk) = n2_3d(i,k,j)

    xland1d (ii) = xland (i,j)
    ivgtyp1d(ii) = ivgtyp(i,j)

    dtemp2d(ii,kk) = dtemp(i,k,j)
    dsupi2d(ii,kk) = dsupi(i,k,j)

! ff for hydrometeors (kg kg-1)
    do m=1, nhydro
    do n=1, nkr
      ffx(ii,kk,m,n) = scalar(i,k,j,(p_start-1)+(m-1)*nkr+n)
    end do
    end do

! ff for aerosol (# kg-1)
    do n=1, nkr
      fccn(ii,kk,n) = scalar(i,k,j,(p_start-1)+nhydro*nkr+n)
    end do

! liquid water for water
    do n=1, nkr
      lw(ii,kk,WATER,n) = ffx(ii,kk,WATER,n)
    end do

! liquid water fraction for ice crystals
    do m=COLUMN, DENDRITE
    do n=1, nkr
      lw(ii,kk,m,n) = 0.0
    end do
    end do

! liquid water fraction for snow, graupel, and hail
    do m=SNOW, HAIL
    do n=1, nkr
      lw(ii,kk,m,n) = scalar(i,k,j,(p_start-1)+(nhydro+(m-SNOW)+1)*nkr+n)
    end do
    end do

! rimed fraction for snow
    do n=1, nkr
      rm(ii,kk,n) = scalar(i,k,j,(p_start-1)+(nhydro+4)*nkr+n)
    end do

    if(inan.eq.1) then
      if(ieee_is_nan(th2d(ii,kk)).or. &
         ieee_is_nan(qv2d(ii,kk)).or. &
         ieee_is_nan(supi(i,k,j)).or. &
         ieee_is_nan(dtemp2d(ii,kk)).or. &
         ieee_is_nan(dsupi2d(ii,kk)).or. &
         ieee_is_nan(dz2d(ii,kk)).or. &
         ieee_is_nan(sum(ffx(ii,kk,:,:))).or. &
         ieee_is_nan(sum(fccn(ii,kk,:))).or. &
         ieee_is_nan(sum(lw(ii,kk,:,:))).or. &
         ieee_is_nan(sum(rm(ii,kk,:)))) then

        write(6,*) i, k, j
        write(6,*) th2d(ii,kk), qv2d(ii,kk), dz2d(ii,kk)
        write(6,*) supi(i,k,j), dtemp2d(ii,kk), dsupi2d(ii,kk)
        write(6,*) sum(ffx(ii,kk,:,:)), sum(fccn(ii,kk,:)), sum(lw(ii,kk,:,:)), sum(rm(ii,kk,:))
        call wrf_error_fatal("nan to mp.")
      end if

      if(iturb.and.(ieee_is_nan(tke2d(ii,kk)).or.ieee_is_nan(n2_2d(ii,kk)))) then
        write(6,*) i, k, j
        call wrf_error_fatal("nan turbulent to mp.")
      end if
    end if

  end do
  end do

!--------------------------------------------------------------------------
  call sbm2d(ite-its+1, kte-kts+1, dx, dy, dt,           &
             th2d, qv2d, pp2d, pi2d, hgt2d, dz2d,        &
             xland1d, ivgtyp1d, isurban,                 &
             dtemp2d, dsupi2d,                           &
             ffx, fccn, rm, lw,                          &
             qc2d, qr2d, qi2d, qs2d, qg2d, qh2d,         &
             nc2d, nr2d, ni2d, ns2d, ng2d, nh2d, na2d,   &
             rain1d, snow1d, graupel1d, hail1d,          &
             dTfm2d, dTc2d, dTn2d, dTce2d, dTds2d,       &
             qr_au2d, qr_ac2d,                           &
             tke2d, n2_2d, eps2d, rey2d,                 &
             its, j, itimestep                           )
!--------------------------------------------------------------------------

! From 2-D to 3-D
  do k=kts, kte
  do i=its, ite

    ii = i-its+1
    kk = k-kts+1

    if(inan.eq.1) then
      if(ieee_is_nan(th2d(ii,kk)).or. &
         ieee_is_nan(qv2d(ii,kk)).or. &
         ieee_is_nan(qc2d(ii,kk)).or. &
         ieee_is_nan(qr2d(ii,kk)).or. &
         ieee_is_nan(qi2d(ii,kk)).or. &
         ieee_is_nan(qg2d(ii,kk)).or. &
         ieee_is_nan(qh2d(ii,kk)).or. &
         ieee_is_nan(nc2d(ii,kk)).or. &
         ieee_is_nan(nr2d(ii,kk)).or. &
         ieee_is_nan(ni2d(ii,kk)).or. &
         ieee_is_nan(ns2d(ii,kk)).or. &
         ieee_is_nan(ng2d(ii,kk)).or. &
         ieee_is_nan(nh2d(ii,kk)).or. &
         ieee_is_nan(na2d(ii,kk)).or. &
         ieee_is_nan(dTc2d(ii,kk)).or. &
         ieee_is_nan(dTn2d(ii,kk)).or. &
         ieee_is_nan(dTfm2d(ii,kk)).or. &
         ieee_is_nan(dTce2d(ii,kk)).or. &
         ieee_is_nan(dTds2d(ii,kk)).or. &
         ieee_is_nan(sum(ffx(ii,kk,:,:))).or. &
         ieee_is_nan(sum(fccn(ii,kk,:))).or. &
         ieee_is_nan(sum(lw(ii,kk,:,:))).or. &
         ieee_is_nan(sum(rm(ii,kk,:))).or. &
         ieee_is_nan(rain1d(ii)).or. &
         ieee_is_nan(snow1d(ii)).or. &
         ieee_is_nan(graupel1d(ii)).or. &
         ieee_is_nan(hail1d(ii))) then

         write(6,*) i, k, j
         write(6,*) th2d(ii,kk), qv2d(ii,kk), qc2d(ii,kk), qr2d(ii,kk)
         write(6,*) qi2d(ii,kk), qs2d(ii,kk), qg2d(ii,kk), qh2d(ii,kk)
         write(6,*) nc2d(ii,kk), nr2d(ii,kk), ni2d(ii,kk), ns2d(ii,kk)
         write(6,*) ng2d(ii,kk), nh2d(ii,kk), na2d(ii,kk)
         write(6,*) dTn2d(ii,kk),  dTc2d(ii,kk)
         write(6,*) dTfm2d(ii,kk), dTce2d(ii,kk), dTds2d(ii,kk)
         write(6,*) sum(ffx(ii,kk,:,:)), sum(fccn(ii,kk,:)), sum(lw(ii,kk,:,:)), sum(rm(ii,kk,:))
         write(6,*) rain1d(ii), snow1d(ii), graupel1d(ii), hail1d(ii)
         call wrf_error_fatal("nan from mp.")
      end if
    end if

    th_phy(i,k,j) = th2d(ii,kk)

    qv(i,k,j) = qv2d(ii,kk)
    qc(i,k,j) = qc2d(ii,kk)
    qr(i,k,j) = qr2d(ii,kk)
    qi(i,k,j) = qi2d(ii,kk)
    qs(i,k,j) = qs2d(ii,kk)
    qg(i,k,j) = qg2d(ii,kk)
    qh(i,k,j) = qh2d(ii,kk)

    nc(i,k,j) = nc2d(ii,kk)
    nr(i,k,j) = nr2d(ii,kk)
    ni(i,k,j) = ni2d(ii,kk)
    ns(i,k,j) = ns2d(ii,kk)
    ng(i,k,j) = ng2d(ii,kk)
    nh(i,k,j) = nh2d(ii,kk)
    na(i,k,j) = na2d(ii,kk)

    dTc3d(i,k,j)  = dTc2d(ii,kk)
    dTn3d(i,k,j)  = dTn2d(ii,kk)
    dTfm3d(i,k,j) = dTfm2d(ii,kk)
    dTce3d(i,k,j) = dTce2d(ii,kk)
    dTds3d(i,k,j) = dTds2d(ii,kk)

    eps_phy(i,k,j) = eps2d(ii,kk)
    rey_phy(i,k,j) = rey2d(ii,kk)

    qr_au(i,k,j) = qr_au2d(ii,kk)
    qr_ac(i,k,j) = qr_ac2d(ii,kk)

    do m=1, nhydro
    do n=1, nkr
      scalar(i,k,j,(p_start-1)+(m-1)*nkr+n) = ffx(ii,kk,m,n)
    end do
    end do

    do n=1, nkr
      scalar(i,k,j,(p_start-1)+nhydro*nkr+n) = fccn(ii,kk,n)
    end do

    do m=SNOW, HAIL
    do n=1, nkr
      scalar(i,k,j,(p_start-1)+(nhydro+(m-SNOW)+1)*nkr+n) = lw(ii,kk,m,n)
    end do
    end do

    do n=1, nkr
      scalar(i,k,j,(p_start-1)+(nhydro+4)*nkr+n) = rm(ii,kk,n)
    end do

    if(k==kts) then
      rainncv(i,j)    = rain1d(ii)
      snowncv(i,j)    = snow1d(ii)
      graupelncv(i,j) = graupel1d(ii)
      hailncv(i,j)    = hail1d(ii)

      rainnc(i,j)    = rainnc(i,j)    + rainncv(i,j)
      snownc(i,j)    = snownc(i,j)    + snowncv(i,j)
      graupelnc(i,j) = graupelnc(i,j) + graupelncv(i,j)
      hailnc(i,j)    = hailnc(i,j)    + hailncv(i,j)
    end if

  end do
  end do

end do

end subroutine sbm

!--------------------------------------------------------------------------
subroutine sbm2d(ni, nk, dx, dy, dt,                          &
                 th, qq, pp, pii, alt2d, dz2d,                &
                 xland1d, ivgtyp1d, isurban,                  &
                 dtemp2d, dsupi2d,                            &
                 ffx, fccnx, rmx, lwx,                        &
                 qcloud, qrain, qice, qsnow, qgraupel, qhail, &
                 ncloud, nrain, nice, nsnow, ngraupel, nhail, &
                 naerosol,                                    &
                 rain_x, snow_x, graupel_x, hail_x,           &
                 dT_fm, dT_c, dT_n, dT_ce, dT_ds,             &
                 qr_au2d, qr_ac2d,                            &
                 tke, n2, eps2d, rey2d,                       &
                 i_id, j_id, itimestep                        )
use ieee_arithmetic
implicit none

integer, intent(in) :: ni, nk, i_id, j_id, isurban, itimestep
real,    intent(in) :: dx, dy, dt

real,    dimension(ni    ),      intent(out  ) ::         &
   rain_x, snow_x, graupel_x, hail_x
real,    dimension(ni    ),      intent(in   ) :: xland1d
integer, dimension(ni    ),      intent(in   ) :: ivgtyp1d
real,    dimension(ni, nk), intent(in   ) :: pp, pii, tke, n2, alt2d, dz2d, dtemp2d, dsupi2d
real,    dimension(ni, nk), intent(inout) :: th, qq
real,    dimension(ni, nk), intent(out  ) ::              &
   qcloud, qrain, qice, qsnow, qgraupel, qhail,           &
   ncloud, nrain, nice, nsnow, ngraupel, nhail, naerosol, &
   dT_fm, dT_c, dT_n, dT_ce, dT_ds,                       &
   eps2d, rey2d, qr_au2d, qr_ac2d
real, dimension(ni, nk, nhydro, nkr), intent(inout) :: ffx
real, dimension(ni, nk,         nkr), intent(inout) :: fccnx
real, dimension(ni, nk, nhydro, nkr), intent(inout) :: lwx
real, dimension(ni, nk,         nkr), intent(inout) :: rmx

! local variables
integer :: i, j, k, m, n, nn1, nn2, nn1s, nn1e, nn2s, nn2e, &
           cs, ce, t1, t2, t31, t32, t3, iter, ivgtype
real, dimension(ni,nk) :: tt
real, dimension(nhydro,nkr) :: ff, lwf, rf, g, ff_temp, m_new, vt_sfc
real, dimension(nkr) :: fccn, ff_hydro1, ff_hydro2, ff_hydro3, m_temp
real, dimension(ni) :: precp
real, dimension(nk) :: rhoaz, altz, dzz, varz, vtz
real :: tt1, pp1, qq1, rhoa1, tke1, n2_1, eps1, rey1,           &
        dT_ce1, dT_ds1, T_old,                                  &
        ffi, ffs_i, ffl_i, lwfl_i, massi, vti, lwfi, rfi, rhoi, &
        fccni, rccni, ffwi, ffi1i, ffi2i, ffi3i, m_newi,        &
        etaa, Nsc, Npr, Nre, m_shed, lw_temp1, lw_temp2,        &
        vtw, m_ref, ak, dz, dtemp, dsupi,                       &
        drop_num_init, drop_mass_init, drop_mass_finl, deldrop, &
        rain_mass_init, rain_mass_finl,                         &
        island, alt

! initialization of output array
rain_x = 0.; snow_x = 0.; graupel_x = 0.; hail_x = 0.
dT_n  = 0.; dT_c  = 0.
dT_fm = 0.; dT_ce = 0.; dT_ds = 0.
qr_au2d = 0.; qr_ac2d = 0.

do i=1, ni
do k=1, nk

  tt(i,k) = th(i,k)*pii(i,k)
  rhoa1 = pp(i,k)/tt(i,k)/Rd

  ! unit conversion ---------------------------------------------------------
  ! ffx : kg kg-1, ff : # m-3
  ! fccnx : # kg-1, fccn : # m-3

  ! caution: at the initialization, the unit of fccn is # m-3, not # kg-1.
  do m=1, nhydro
    ff(m,:) = max(ffx(i,k,m,:)/mass_hydro*rhoa1,0.0)
  end do
  if(itimestep.gt.1) then 
    fccn = max(fccnx(i,k,:)*rhoa1,0.0)
  else
    fccn = max(fccnx(i,k,:),0.0)
  end if

  do m=1, nhydro
  do n=1, nkr
    if(m.eq.WATER) then
      lwf(m,n) = 1.0
    else if(m.ge.COLUMN.and.m.le.DENDRITE) then
      lwf(m,n) = 0.0
    else
      if(ffx(i,k,m,n).gt.0.0) then
        lwf(m,n) = min(max(lwx(i,k,m,n)/ffx(i,k,m,n),0.0),1.0)
      else
        lwf(m,n) = 0.0
      end if
    end if
  end do
  end do

  do m=1, nhydro
  do n=1, nkr
    if(m.eq.SNOW) then
      if(ffx(i,k,m,n).gt.0.0.and.lwf(m,n).lt.1.0) then
        rf(m,n) = min(max(rmx(i,k,n)/(ffx(i,k,m,n)*(1.0-lwf(m,n))),0.0),1.0)
      else
        rf(m,n) = 0.0
      end if
    else
      rf(m,n) = 0.0
    end if
  end do
  end do

  if(inan.eq.1) then
    if(ieee_is_nan(tt(i,k)).or. &
       ieee_is_nan(pp(i,k)).or. &
       ieee_is_nan(qq(i,k))) then
      write(6,*) i_id+i-1, j_id, k
      write(6,*) tt(i,k), pp(i,k), qq(i,k)
      call wrf_error_fatal("nan in tt before integration")
    end if
  end if
   
  call check_validity(ff, lwf, rf, "start")

  call vt_pressure(pp(i,k))

  call hydrometeor_density(lwf, rf)

  ! freezing and melting -----------------------------------------------------
  T_old = tt(i,k)

  if(tt(i,k).lt.t00) then

    pp1    = pp(i,k)
    qq1    = qq(i,k)
    rhoa1  = pp1/T_old/Rd
    dtemp  = dtemp2d(i,k)

    ! drop freezing
    do n=1, nkr

      tt1    = T_old
      massi  = mass_hydro(n)
      ffi    = ff(WATER,n)
      ffs_i  = ff(COLUMN,n)
      ffl_i  = ff(HAIL,n)
      lwfl_i = lwf(HAIL,n)

      if(ffi.gt.0.0) then
        call freeze_drop(tt1, rhoa1, massi, ffi, ffs_i, ffl_i, lwfl_i, dt, dtemp)

        ff(WATER,n) = ffi
        ff(COLUMN,n) = ffs_i
        ff(HAIL,n) = ffl_i
        lwf(HAIL,n) = lwfl_i

        dT_fm(i,k) = dT_fm(i,k) + (tt1 - T_old)
      end if

    end do

    call check_validity(ff, lwf, rf, "freeze_drop")

    ! reduce lwf
    do m=COLUMN, HAIL
    do n=1, nkr

      tt1    = T_old
      massi  = mass_hydro(n)
      vtw    = vt_hydro(WATER,n)
      vti    = vt_hydro(m,n)
      rhoi   = rho_hydro(m,n)
      ffi    = ff(m,n)
      lwfi   = lwf(m,n)
 
      if(ffi.gt.0.0.and.lwfi.gt.0.0) then
        call freeze_lwf(dt, tt1, pp1, qq1, rhoa1, massi, vtw, vti, rhoi, ffi, lwfi)

        vt_hydro(m,n) = vti
        lwf(m,n) = lwfi
        dT_fm(i,k) = dT_fm(i,k) + (tt1 - T_old)
      end if

    end do
    end do

    call check_validity(ff, lwf, rf, "freeze_lwf")

    tt(i,k) = tt(i,k) + dT_fm(i,k)

  ! melting -----------------------------------------------------------------
  else

    do n=1, nkr
      lw_temp1 = ff(SNOW,n)*lwf(SNOW,n)

      do m=COLUMN, DENDRITE
        ff(SNOW,n) = ff(SNOW,n)+ff(m,n)
        ff(m,n) = 0.
      end do

      if(ff(SNOW,n).gt.0.0) then
        lwf(SNOW,n) = lw_temp1/ff(SNOW,n)
      else
        lwf(SNOW,n) = 0.0
      end if
    end do

    pp1   = pp(i,k)
    qq1   = qq(i,k)
    rhoa1 = pp1/T_old/Rd

    do m=SNOW, HAIL
    do n=1, nkr

      tt1   = T_old
      massi = mass_hydro(n)
      vtw   = vt_hydro(WATER,n)
      vti   = vt_ref(m,n)
      if(m.eq.SNOW) then
        rhoi = rho_ice_snow(n)
      else
        rhoi = rho_hydro(m,n)
      end if
      ffi   = ff(m,n)
      lwfi  = lwf(m,n)

      if(ffi.gt.0.0.and.lwfi.lt.lwf_upper_limit) then
        call melting(dt, tt1, pp1, qq1, rhoa1, massi, vtw, vti, rhoi, ffi, lwfi, m, Nre, m_shed)

        vt_hydro(m,n) = vti
        lwf(m,n) = lwfi

        dT_fm(i,k) = dT_fm(i,k) + (tt1 - T_old)

        if(inan.eq.1) then
          if(ieee_is_nan(lwfi).or..not.ieee_is_finite(lwfi)) then
            write(6,*) "Before shedding", m, n
            write(6,*) Nre, m_shed, ffi, lwfi, massi
            call wrf_error_fatal("nan")
          end if
        end if

!        if(m_shed.gt.0.0) then
!          lwfi      = lwf(m,n)
!          ff_hydro1 = ff(WATER,:)
!          ff_hydro2 = ff(m,:)
!          ff_hydro3 = lwf(m,:)
!          call shedding(Nre, m_shed, ffi, lwfi, massi, ff_hydro1, ff_hydro2, ff_hydro3, n)

!          do nn1=1, nkr
!            if(ieee_is_nan(ff_hydro3(nn1))) then
!              write(6,*) "After shedding", m, n, nn1
!              write(6,*) Nre, m_shed, ffi, lwfi, massi
!              write(6,*) ff(m,:)
!              write(6,*) lwf(m,:)
!              write(6,*) ff_hydro2
!              write(6,*) ff_hydro3
!              call wrf_error_fatal("nan")
!            end if
!          end do

!          ff(WATER,:) = ff_hydro1
!          ff(m,:)     = ff_hydro2
!          lwf(m,:)    = ff_hydro3
!        end if

      end if

      if(lwfi.ge.lwf_upper_limit) then
        ff(WATER,n) = ff(WATER,n) + ff(m,n)
        ff(m,n) = 0.0
        lwf(m,n) = 0.0
        rf(m,n) = 0.0
      end if

    end do
    end do

    call check_validity(ff, lwf, rf, "melting")

    tt(i,k) = tt(i,k) + dT_fm(i,k)

  end if

  ! collision ---------------------------------------------------------------
  tt1   = tt(i,k)
  pp1   = pp(i,k)
  qq1   = qq(i,k)
  rhoa1 = pp1/tt1/Rd
  tke1  = tke(i,k)
  n2_1  = n2(i,k)
  dz    = dz2d(i,k)

  call kernel_pressure(pp1)

  call kernel_mod(tt1, pp1, qq1, tke1, n2_1, dx, dy, dz, rhoa1, eps1, rey1)

  eps2d(i,k) = eps1
  rey2d(i,k) = rey1

  do m=1, nhydro
    g(m,:) = ff(m,:)*mass_hydro/dlnr ! unit: kg m-3 / dlnr
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! water-water: sepatated collision of c-c, c-r, r-r
  t1 = WATER; t2 = WATER; t3 = WATER

! cloud-cloud
  rain_mass_init = sum(g(WATER,kr_droplet+1:))*dlnr
  do nn1 =     1, kr_droplet
  do nn2 = nn1+1, kr_droplet
    if(g(t1,nn1).gt.0.0.and.g(t2,nn2).gt.0.0.and.ckxx(t1,t2,nn1,nn2).gt.0.0) then
      if(icollision.eq.1) then
        call collision(g, lwf, rf, t1, t2, t3, nn1, nn2, tt1)
      else if(icollision.eq.2) then
        call collision_poisson(g, lwf, rf, t1, t2, t3, nn1, nn2, tt1)
      else
        call wrf_error_fatal("something wrong in icollision")
      end if
    end if
  end do
  end do
  rain_mass_finl = sum(g(WATER,kr_droplet+1:))*dlnr
  qr_au2d(i,k) = (rain_mass_finl - rain_mass_init)/dt

! cloud-rain
  rain_mass_init = rain_mass_finl
  do nn1 =            1, kr_droplet
  do nn2 = kr_droplet+1, nkr
    if(g(t1,nn1).gt.0.0.and.g(t2,nn2).gt.0.0.and.ckxx(t1,t2,nn1,nn2).gt.0.0) then
      if(icollision.eq.1) then
        call collision(g, lwf, rf, t1, t2, t3, nn1, nn2, tt1)
      else if(icollision.eq.2) then
        call collision_poisson(g, lwf, rf, t1, t2, t3, nn1, nn2, tt1)
      else
        call wrf_error_fatal("something wrong in icollision")
      end if
    end if
  end do
  end do
  rain_mass_finl = sum(g(WATER,kr_droplet+1:))*dlnr
  qr_ac2d(i,k) = (rain_mass_finl - rain_mass_init)/dt

! rain-rain
  do nn1 = kr_droplet+1, nkr
  do nn2 =        nn1+1, nkr
    if(g(t1,nn1).gt.0.0.and.g(t2,nn2).gt.0.0.and.ckxx(t1,t2,nn1,nn2).gt.0.0) then
      if(icollision.eq.1) then
        call collision(g, lwf, rf, t1, t2, t3, nn1, nn2, tt1)
      else if(icollision.eq.2) then
        call collision_poisson(g, lwf, rf, t1, t2, t3, nn1, nn2, tt1)
      else
        call wrf_error_fatal("something wrong in icollision")
      end if
    end if
  end do
  end do

  ! breakup by Bleck (1970)
  do nn1 = 1, nkr
    ff_hydro1(nn1) = max(g(WATER,nn1),0.0)
  end do
  call breakup_water(ff_hydro1,dt)
  do nn1 = 1, nkr
    g(WATER,nn1) = max(ff_hydro1(nn1),0.0)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! other collisions

  drop_mass_init = sum(g(WATER,:))*dlnr

  if(tt1.lt.t00) then
    cs = 1
    ce = num_coll_type_lt0
  else
    cs = num_coll_type_lt0+1
    ce = num_coll_type_lt0+num_coll_type_ge0
  end if

  do n=cs, ce
    t1  = coll_table(n,1)
    t2  = coll_table(n,2)
    t31 = coll_table(n,3)
    t32 = coll_table(n,4)

    if(tt1.gt.t_icempl1.and.tt1.lt.t_icempl2.and.t1.eq.WATER.and.t2.eq.GRAUPEL) then
      drop_num_init = sum(g(WATER,kr_icempl:nkr)/mass_hydro(kr_icempl:nkr))*dlnr
    end if

    do nn1=1, nkr
    do nn2=1, nkr

      if(t1.eq.t2.and.nn1.eq.nn2) cycle  ! should be modified

      if(g(t1,nn1).gt.0.0.and.g(t2,nn2).gt.0.0.and.ckxx(t1,t2,nn1,nn2).gt.0.0) then
        if(nn2.gt.nn1) then
          if(icollision.eq.1) then
            call collision(g, lwf, rf, t1, t2, t32, nn1, nn2, tt1)
          else if(icollision.eq.2) then
            call collision_poisson(g, lwf, rf, t1, t2, t32, nn1, nn2, tt1)
          else
            call wrf_error_fatal("something wrong in icollision")
          end if
        else
          if(t1.ne.t2) then
            if(icollision.eq.1) then
              call collision(g, lwf, rf, t2, t1, t31, nn2, nn1, tt1)
            else if(icollision.eq.2) then
              call collision_poisson(g, lwf, rf, t2, t1, t31, nn2, nn1, tt1)
            else
              call wrf_error_fatal("something wrong in icollision")
            end if
          end if
        end if
      end if

    end do
    end do

    ! ice multiplication by Hallett and Mossop (1974)
    ! Note that the mass of hydrometeors does not conserved.
    if(tt1.gt.t_icempl1.and.tt1.lt.t_icempl2.and.t1.eq.WATER.and.t2.eq.GRAUPEL) then
      ff_hydro1 = g(WATER,:)
      ff_hydro2 = g(PLATE,:)
      call ice_mpl(ff_hydro1,ff_hydro2,drop_num_init,tt1)
      g(PLATE,:) = ff_hydro2
    end if

  end do

  drop_mass_finl = sum(g(WATER,:))*dlnr
  deldrop = drop_mass_finl - drop_mass_init

  if(deldrop.gt.0.0) then
    if(abs(deldrop/drop_mass_init).gt.tol) then
      write(6,*) "Water increases by collision.", drop_mass_init, drop_mass_finl, deldrop/drop_mass_init
    end if
  else
    dT_c(i,k) = -deldrop*Lm/rhoa1/Cp
    tt(i,k) = tt(i,k) + dT_c(i,k)
  end if

  do m=1, nhydro
    ff(m,:) = g(m,:)/mass_hydro*dlnr
  end do

  call check_validity(ff, lwf, rf, "collision")

  ff_hydro1 = ff(SNOW,:)
  ff_hydro2 = lwf(SNOW,:)
  ff_hydro3 = rf(SNOW,:)
  call snow_breakup(ff_hydro1, ff_hydro2, ff_hydro3)
   ff(SNOW,:) = ff_hydro1
  lwf(SNOW,:) = ff_hydro2
   rf(SNOW,:) = ff_hydro3

  call check_validity(ff, lwf, rf, "snow breakup")

  ! vapor diffusion ---------------------------------------------------------
  do iter=1, n_cond

    if(maxval(ff).le.0.0) exit

    tt1   = tt(i,k)
    qq1   = qq(i,k)
    pp1   = pp(i,k)
    rhoa1 = pp1/tt1/Rd

    ff_temp = ff/rhoa1 ! unit: # kg-1

    call diffuse(tt1, pp1, qq1, rhoa1, dt_cond, ff_temp, m_new, dT_ce1, dT_ds1)

    tt(i,k) = tt1
    qq(i,k) = qq1
    dT_ce(i,k) = dT_ce(i,k) + dT_ce1
    dT_ds(i,k) = dT_ds(i,k) + dT_ds1

    do m=1, nhydro
      if(maxval(ff(m,:)).gt.0.0) then
        m_temp    = m_new(m,:)
        ff_hydro1 = ff_temp(m,:)
        call remapping(m_temp,ff_hydro1)
        ff_temp(m,:) = ff_hydro1
      end if
    end do

    ff = ff_temp*rhoa1

  end do

  call check_validity(ff, lwf, rf, "diffusion")

  ! nucleation --------------------------------------------------------------
  T_old = tt(i,k)
  rccni = r_hydro(WATER,1)

  do n=nkr, 1, -1 ! large aerosol first

    fccni = fccn(n)

    if(fccni.gt.0.0) then
      tt1   = tt(i,k)
      qq1   = qq(i,k)
      pp1   = pp(i,k)
      rhoa1 = pp1/tt1/Rd
      ffwi  = ff(WATER,max(n-35,1))

      call water_nucl(tt1, pp1 ,qq1 ,rhoa1 ,fccni, rccni, ffwi)

      tt(i,k) = tt1
      qq(i,k) = qq1
      fccn(n) = fccni
      ff(WATER,max(n-35,1)) = ffwi
    end if

    rccni = rccni/exp(dlnr)

  end do

  call check_validity(ff, lwf, rf, "nucl_wat")

  dsupi = dsupi2d(i,k)
  dtemp = dtemp2d(i,k)

  if(dsupi.gt.0.0.or.dtemp.lt.0.0) then
    tt1   = tt(i,k)
    qq1   = qq(i,k)
    pp1   = pp(i,k)
    rhoa1 = pp1/tt1/Rd
    ffi1i = ff(COLUMN,1)
    ffi2i = ff(PLATE,1)
    ffi3i = ff(DENDRITE,1)

    call ice_nucl(tt1,pp1,qq1,rhoa1,dsupi,dtemp,ffi1i,ffi2i,ffi3i)

    tt(i,k) = tt1
    qq(i,k) = qq1
    ff(COLUMN,  1) = ffi1i
    ff(PLATE,   1) = ffi2i
    ff(DENDRITE,1) = ffi3i
  end if

  call check_validity(ff, lwf, rf, "nucl_ice")

  dT_n(i,k) = tt(i,k) - T_old

  ! type conversion ---------------------------------------------------------
  do m=COLUMN, HAIL
  do n=1, nkr
    if(ff(m,n).gt.0.0) then
      if(lwf(m,n).gt.lwf_upper_limit) then
        ff(WATER,n) = ff(WATER,n) + ff(m,n)
        ff(m,n) = 0.0
        lwf(m,n) = 0.0
        rf(m,n) = 0.0
      end if
    end if
  end do
  end do

  do n=kr_big_ice, nkr
    if(maxval(ff(COLUMN:DENDRITE,n)).gt.0.0) then
      lw_temp1 = ff(HAIL,n)*lwf(HAIL,n)
      lw_temp2 = ff(SNOW,n)*lwf(SNOW,n)

      ff(HAIL,n) = ff(HAIL,n) + ff(COLUMN,n)
      ff(COLUMN,n) = 0.0

      ff(HAIL,n) = ff(HAIL,n) + ff(PLATE,n)
      ff(PLATE,n) = 0.0

      ff(SNOW,n) = ff(SNOW,n) + ff(DENDRITE,n)
      ff(DENDRITE,n) = 0.0

      if(ff(HAIL,n).gt.0.0) then
        lwf(HAIL,n) = lw_temp1/ff(HAIL,n)
      else
        lwf(HAIL,n) = 0.0
      end if

      if(ff(SNOW,n).gt.0.0) then
        lwf(SNOW,n) = lw_temp2/ff(SNOW,n)
      else
        lwf(SNOW,n) = 0.0
      end if
    end if
  end do

  do n=1, nkr
    if(ff(SNOW,n).gt.0.0) then
      if(n.le.kr_snow) then
        if(rf(SNOW,n).gt.rf_upper_limit) then
          lw_temp1 = ff(HAIL,n)*lwf(HAIL,n)
          lw_temp2 = ff(SNOW,n)*lwf(SNOW,n)

          ff(HAIL,n) = ff(HAIL,n) + ff(SNOW,n)
          ff(SNOW,n) = 0.0
          lwf(SNOW,n) = 0.0
          rf(SNOW,n) = 0.0

          if(ff(HAIL,n).gt.0.0) then
            lwf(HAIL,n) = (lw_temp1+lw_temp2)/ff(HAIL,n)
          else
            lwf(HAIL,n) = 0.0
          end if
        end if
      else
        if(rho_hydro(SNOW,n).gt.rho_upper_limit) then
          lw_temp1 = ff(GRAUPEL,n)*lwf(GRAUPEL,n)
          lw_temp2 = ff(SNOW,n)*lwf(SNOW,n)

          ff(GRAUPEL,n) = ff(GRAUPEL,n) + ff(SNOW,n)
          ff(SNOW,n) = 0.0
          lwf(SNOW,n) = 0.0
          rf(SNOW,n) = 0.0

          if(ff(GRAUPEL,n).gt.0.0) then
            lwf(GRAUPEL,n) = (lw_temp1+lw_temp2)/ff(GRAUPEL,n)
          else
            lwf(GRAUPEL,n) = 0.0
          end if
        end if
      end if
    end if
  end do

  do n=kr_graupel, nkr
    if(ff(GRAUPEL,n).gt.0.0) then
      lw_temp1 = ff(GRAUPEL,n)*lwf(GRAUPEL,n)
      lw_temp2 = ff(HAIL,n)*lwf(HAIL,n)

      ff(HAIL,n) = ff(HAIL,n) + ff(GRAUPEL,n)
      ff(GRAUPEL,n) = 0.0

      if(ff(HAIL,n).gt.0.0) then
        lwf(HAIL,n) = (lw_temp1+lw_temp2)/ff(HAIL,n)
      else
        lwf(HAIL,n) = 0.0
      end if
    end if
  end do

  call check_validity(ff, lwf, rf, "conversion")

  ! aerosol replenishment by relaxation (Hyunho) ----------------------------
  island = xland1d(i)
  ivgtype = ivgtyp1d(i)
  alt = alt2d(i,k)

  call fccn_init(alt, island, ivgtype, isurban, ff_hydro1)

  do n=1, nkr
    fccn(n) = fccn(n) + (ff_hydro1(n) - fccn(n))/tau_ccn*dt
  end do

  ! unit conversion ---------------------------------------------------------
  rhoa1 = pp(i,k)/tt(i,k)/Rd
  do m=1, nhydro
    ffx(i,k,m,:) = ff(m,:)*mass_hydro/rhoa1
  end do
  lwx(i,k,:,:) = ffx(i,k,:,:)*lwf(:,:)
  rmx(i,k,:)   = ffx(i,k,SNOW,:)*(1.0-lwf(SNOW,:))*rf(SNOW,:)
  fccnx(i,k,:) = fccn/rhoa1

  th(i,k) = tt(i,k)/pii(i,k)

end do
end do

! sedimentation --------------------------------------------------------------
do i=1, ni

  rhoaz = pp(i,:)/tt(i,:)/Rd
  altz = alt2d(i,:)
  dzz = dz2d(i,:)

  do m=1, nhydro
  do n=1, nkr

    vtz = bulkvt(m,n)*sqrt(p1000hpa/pp(i,:)) ! why not vt_hydro?
    vt_sfc(m,n) = vtz(1)

    varz = ffx(i,:,m,n)
    call sediment(varz, vtz, rhoaz, altz, dzz, nk, dt)
    ffx(i,:,m,n) = varz

    if(m.ge.SNOW.and.m.le.HAIL) then
      varz = lwx(i,:,m,n)
      call sediment(varz, vtz, rhoaz, altz, dzz, nk, dt)
      lwx(i,:,m,n) = varz
    end if

    if(m.eq.SNOW) then
      varz = rmx(i,:,n)
      call sediment(varz, vtz, rhoaz, altz, dzz, nk, dt)
      rmx(i,:,n) = varz
    end if

  end do
  end do

  rhoa1 = pp(i,1)/tt(i,1)/Rd

  do m=1, nhydro
  do n=1, nkr

    precp(i) = ffx(i,1,m,n)/bulkrho(m,n)*vt_sfc(m,n)*rhoa1*dt*1000.

    if(m.eq.WATER) then
      rain_x(i) = rain_x(i) + precp(i)
    else if(m.eq.SNOW) then
      snow_x(i) = snow_x(i) + precp(i)
    else if(m.eq.GRAUPEL) then
      graupel_x(i) = graupel_x(i) + precp(i)
    else if(m.eq.HAIL) then
      hail_x(i) = hail_x(i) + precp(i)
    end if

    if(inan.eq.1) then
      if(ieee_is_nan(precp(i))) then
        write(6,*) "in precipitation", i, m, n
        write(6,*) ffx(i,1,m,n),bulkrho(m,n),vt_sfc(m,n)
        call wrf_error_fatal("nan detect")
      end if
    end if

  end do
  end do

end do

! WRF output ---------------------------------------------------------------
do i=1, ni
do k=1, nk

  qcloud  (i,k) = sum(ffx(i,k,WATER,1:kr_droplet))
  qrain   (i,k) = sum(ffx(i,k,WATER,kr_droplet+1:nkr))
  qice    (i,k) = sum(ffx(i,k,COLUMN:DENDRITE,:))
  qsnow   (i,k) = sum(ffx(i,k,SNOW,:))
  qgraupel(i,k) = sum(ffx(i,k,GRAUPEL,:))
  qhail   (i,k) = sum(ffx(i,k,HAIL,:))

  ncloud  (i,k) = sum(ffx(i,k,WATER,1:kr_droplet)/mass_hydro(1:kr_droplet))
  nrain   (i,k) = sum(ffx(i,k,WATER,kr_droplet+1:nkr) &
                     /mass_hydro(kr_droplet+1:nkr))
  nice    (i,k) = sum(ffx(i,k,COLUMN,  :)/mass_hydro) &
                + sum(ffx(i,k,PLATE,   :)/mass_hydro) &
                + sum(ffx(i,k,DENDRITE,:)/mass_hydro)
  nsnow   (i,k) = sum(ffx(i,k,SNOW,:)/mass_hydro)
  ngraupel(i,k) = sum(ffx(i,k,GRAUPEL,:)/mass_hydro)
  nhail   (i,k) = sum(ffx(i,k,HAIL,:)/mass_hydro)

  naerosol(i,k) = sum(fccnx(i,k,:))

end do
end do

end subroutine sbm2d

!==============================================================================
subroutine sediment(varz, vtz, rhoaz, altz, dzz, nz, dt)
implicit none

integer, intent(in) :: nz
real, intent(in) :: dt
real, dimension(nz), intent(in) :: vtz, rhoaz, altz, dzz
real, dimension(nz), intent(inout) :: varz

real :: tf, dt_fall
integer :: iter, it, k
real, dimension(nz)   :: dzz_dt
real, dimension(nz,2) :: flux

tf = minval(altz/(vtz+1.0e-10))
iter = int((2.0*dt)/tf)+1
dt_fall = dt/iter
dzz_dt = dt_fall/dzz

do it=1, iter
  ! hydrometeors that have terminal velocity
  do k=1, nz-1
    flux(k,1) = vtz(k)*dzz_dt(k)*varz(k)
    flux(k,2) = rhoaz(k+1)/rhoaz(k)*vtz(k+1)*dzz_dt(k)*varz(k+1)
  end do
  flux(nz,1) = vtz(nz)*dzz_dt(nz)*varz(nz)
  flux(nz,2) = 0.0

  do k=1, nz
    varz(k) = varz(k)-flux(k,1)+flux(k,2)
  end do
end do

end subroutine sediment

!==============================================================================
subroutine sbminit(ids, ide, jds, jde, kds, kde,      &
                   ims, ime, jms, jme, kms, kme,      &
                   its, ite, jts, jte, kts, kte,      &
                   xland, ivgtyp, isurban, hgt,       &
                   n_scalar, scalar, dt, restart      )
implicit none

real, intent(in) :: dt
logical, intent(in) :: restart
integer, intent(in) :: n_scalar
integer, intent(in) :: isurban,                      &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte
real, dimension(ims:ime,jms:jme), intent(in) :: xland
integer, dimension(ims:ime,jms:jme), intent(in) :: ivgtyp
real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: hgt
real, dimension(ims:ime,kms:kme,jms:jme,n_scalar), intent(inout) :: scalar

integer, parameter :: k0l_gl=16
integer, parameter :: k0g_gl=16

real, dimension(21) :: rl_turb
real, dimension(21,21,11) :: turbll_text1
real, dimension(13,13,11) :: turbll_text2
real, dimension(k0g_gl,k0l_gl) :: turbgl_text_400, turbgl_text_800
real, dimension(k0l_gl) :: rl_gl
real, dimension(k0g_gl) :: rg_gl

real, dimension(nkr) :: fccni
real, dimension(nkr,nhydro) :: temp
real, dimension(3,nkr,nkr) :: temp1
real, dimension(COLUMN:DENDRITE,3,nkr,nkr) :: temp2
real, dimension(COLUMN:DENDRITE,COLUMN:DENDRITE,3,nkr,nkr) :: temp3

integer :: hujisbm_unit1
logical :: opened
logical, external :: wrf_dm_on_monitor

integer :: i, j, k, kp, ip, jp, kq, jq, m1, m2, ivgtype
real :: x, y, mass0, alt, island

ywxx = 0.0

print *,'wrf_dm_on_monitor() =',wrf_dm_on_monitor()

!--- set timestep for diffusion
n_cond = max(nint(dt/dt_cond_ref),1)
dt_cond = dt/n_cond
write(6,*) "n_cond, dt_cond = ", n_cond, dt_cond

!--- read in various lookup tables
hujisbm_unit1=-1
if ( wrf_dm_on_monitor() ) then
  do i = 31,99
    inquire ( i , opened = opened )
    if ( .not. opened ) then
      hujisbm_unit1 = i
      exit
    endif
  enddo
endif

call wrf_dm_bcast_bytes ( hujisbm_unit1 , IWORDSIZE )
if ( hujisbm_unit1 < 0 ) then
  call wrf_error_fatal ( 'module_mp_sbm: cannot find unused unit to read in lookup table.')
endif

! model configuration parameters
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="model_config.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,*) n0ccn1
  read(hujisbm_unit1,*) kccn1
  read(hujisbm_unit1,*) n0ccn2
  read(hujisbm_unit1,*) kccn2
  read(hujisbm_unit1,*) urb_ccnf
  read(hujisbm_unit1,*) iremap
  read(hujisbm_unit1,*) icollision
  read(hujisbm_unit1,*) iturb
  read(hujisbm_unit1,*) ibreakup
  read(hujisbm_unit1,*) inan
  close(hujisbm_unit1)
  call wrf_message(' ***** file_config: succesfull *******')
endif
call wrf_dm_bcast_bytes ( n0ccn1 ,     RWORDSIZE )
call wrf_dm_bcast_bytes ( kccn1 ,      RWORDSIZE )
call wrf_dm_bcast_bytes ( n0ccn2 ,     RWORDSIZE )
call wrf_dm_bcast_bytes ( kccn2 ,      RWORDSIZE )
call wrf_dm_bcast_bytes ( urb_ccnf ,   RWORDSIZE )
call wrf_dm_bcast_bytes ( iremap ,     IWORDSIZE )
call wrf_dm_bcast_bytes ( icollision , IWORDSIZE )
call wrf_dm_bcast_bytes ( iturb ,      IWORDSIZE )
call wrf_dm_bcast_bytes ( ibreakup ,   IWORDSIZE )
call wrf_dm_bcast_bytes ( inan ,       IWORDSIZE )

write(6,*) "Read configuration from file"
write(6,*) "n0 over the land:       ", n0ccn1
write(6,*) "k over the land:        ", kccn1
write(6,*) "n0 over the ocean:      ", n0ccn2
write(6,*) "k over the ocean:       ", kccn2
write(6,*) "n0 multiplier for urban:", urb_ccnf
write(6,*) "iremap:    ", iremap
write(6,*) "icollision:", icollision
write(6,*) "iturb:     ", iturb
write(6,*) "ibreakup:  ", ibreakup
write(6,*) "inan:      ", inan

! capacitance
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="capacity43.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) temp
  close(hujisbm_unit1)
  call wrf_message(' ***** file1: succesfull *******')
  cap_hydro = transpose(temp)
end if
call wrf_dm_bcast_bytes ( cap_hydro , size ( cap_hydro ) * RWORDSIZE )

! masses :
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="masses43.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) mass_hydro
  close(hujisbm_unit1)
  call wrf_message(' ***** file2: succesfull *******')
endif
call wrf_dm_bcast_bytes ( mass_hydro , size ( mass_hydro ) * RWORDSIZE )

! terminal velosity :
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="termvels43.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) temp
  close(hujisbm_unit1)
  call wrf_message(' ***** file3: succesfull *******')
  bulkvt = transpose(temp)
endif
call wrf_dm_bcast_bytes ( bulkvt , size ( bulkvt ) * RWORDSIZE )

! bulkdensity :
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="bulkdens43.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) temp
  close(hujisbm_unit1)
  call wrf_message(' ***** file4: succesfull *******')
  bulkrho = transpose(temp)
end if
call wrf_dm_bcast_bytes ( bulkrho , size ( bulkrho ) * RWORDSIZE )

! bulkradius
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="bulkradii43.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,*) temp
  close(hujisbm_unit1)
  call wrf_message(' ***** file5: succesfull *******')
  r_hydro = transpose(temp)
end if
call wrf_dm_bcast_bytes ( r_hydro , size ( r_hydro ) * RWORDSIZE )

! kernels (depending on pressure) :
! water-water for 1000, 750, 500 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="ckll_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) (((ywxx(WATER,WATER,k,i,j),i=1,nkr),j=1,nkr),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file51: succesfull *******')
end if
call wrf_dm_bcast_bytes ( ywxx(WATER,WATER,:,:,:) , size ( temp1 ) * RWORDSIZE )

! water-ice for 750, 500, 300 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="ckli_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) ((((ywxx(WATER,m1,k,i,j),i=1,nkr),j=1,nkr),m1=COLUMN,DENDRITE),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file52: succesfull *******')
end if
call wrf_dm_bcast_bytes ( ywxx(WATER,COLUMN:DENDRITE,:,:,:) , size ( temp2 ) * RWORDSIZE )

! water-snow for 750, 500, 300 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="ckls_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) (((ywxx(WATER,SNOW,k,i,j),i=1,nkr),j=1,nkr),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file53: succesfull *******')
end if
call wrf_dm_bcast_bytes ( ywxx(WATER,SNOW,:,:,:) , size ( temp1 ) * RWORDSIZE )

! water-graupel for 750, 500, 300 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="cklg_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) (((ywxx(WATER,GRAUPEL,k,i,j),i=1,nkr),j=1,nkr),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file54: succesfull *******')
end if
call wrf_dm_bcast_bytes ( ywxx(WATER,GRAUPEL,:,:,:) , size ( temp1 ) * RWORDSIZE )

! water-hail for 750, 500, 300 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="cklh_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) (((ywxx(WATER,HAIL,k,i,j),i=1,nkr),j=1,nkr),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file55: succesfull *******')
end if
call wrf_dm_bcast_bytes ( ywxx(WATER,HAIL,:,:,:) , size ( temp1 ) * RWORDSIZE )

! ice-ice for 750, 500, 300 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="ckii_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) (((((ywxx(m1,m2,k,i,j),i=1,nkr),j=1,nkr),m2=COLUMN,DENDRITE),m1=COLUMN,DENDRITE),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file56: succesfull *******')
  end if
call wrf_dm_bcast_bytes ( ywxx(COLUMN:DENDRITE,COLUMN:DENDRITE,:,:,:) , size ( temp3 ) * RWORDSIZE )

! ice-snow for 750, 500, 300 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="ckis_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) ((((ywxx(m1,SNOW,k,i,j),i=1,nkr),j=1,nkr),m1=COLUMN,DENDRITE),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file57: succesfull *******')
end if
call wrf_dm_bcast_bytes ( ywxx(COLUMN:DENDRITE,SNOW,:,:,:) , size ( temp2 ) * RWORDSIZE )

! snow-snow for 750, 500, 300 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="ckss_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) (((ywxx(SNOW,SNOW,k,i,j),i=1,nkr),j=1,nkr),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file58: succesfull *******')
end if
call wrf_dm_bcast_bytes ( ywxx(SNOW,SNOW,:,:,:) , size ( temp1 ) * RWORDSIZE )

! snow-graupel for 750, 500, 300 hPa
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="cksg_z.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,900) (((ywxx(SNOW,GRAUPEL,k,i,j),i=1,nkr),j=1,nkr),k=1,3)
  close(hujisbm_unit1)
  call wrf_message(' ***** file59: succesfull *******')
end if
call wrf_dm_bcast_bytes ( ywxx(SNOW,GRAUPEL,:,:,:) , size ( temp1 ) * RWORDSIZE )

! collision enhancement by turbulence (water-graupel)
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="turb_gl_400.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,'(16f8.2)') ((turbgl_text_400(i,j),i=1,k0g_gl),j=1,k0l_gl)
  close(hujisbm_unit1)
  call wrf_message(' ***** file7: succesfull *******')
end if
call wrf_dm_bcast_bytes ( turbgl_text_400 , size ( turbgl_text_400 ) * RWORDSIZE )

if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="turb_gl_800.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,'(16f8.2)') ((turbgl_text_800(i,j),i=1,k0g_gl),j=1,k0l_gl)
  close(hujisbm_unit1)
  call wrf_message(' ***** file7: succesfull *******')
end if
call wrf_dm_bcast_bytes ( turbgl_text_800 , size ( turbgl_text_800 ) * RWORDSIZE )

! collision enhancement by turbulence (water-water)
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="turb_ll_part1.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,'(21e8.4)') turbll_text1
  close(hujisbm_unit1)
  call wrf_message(' ***** file8: succesfull *******')
end if
call wrf_dm_bcast_bytes ( turbll_text1 , size ( turbll_text1 ) * RWORDSIZE )

if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="turb_ll_part2.asc",status="old",action="read",err=2070)
  read(hujisbm_unit1,'(13e8.4)') turbll_text2
  close(hujisbm_unit1)
  call wrf_message(' ***** file9: succesfull *******')
end if
call wrf_dm_bcast_bytes ( turbll_text2 , size ( turbll_text2 ) * RWORDSIZE )

! breakup table
if ( wrf_dm_on_monitor() ) then
  open(unit=hujisbm_unit1,file="coeff_p43.asc",status="old",action="read",err=2070)
  do k=1,krbreak
  do i=1,krbreak
  do j=1,i
    read(hujisbm_unit1,*) kp,ip,jp,pkij(k,i,j) ! [pkij]=(g^3*cm^3)/s
  enddo
  enddo
  enddo
  close(hujisbm_unit1)

  open(unit=hujisbm_unit1,file="coeff_q43.asc",status="old",action="read",err=2070)
  do k=1,krbreak
  do j=1,krbreak
    read(hujisbm_unit1,*) kq,jq,qkj(k,j)       ! [qkj]=(g^3*cm^3)/s
  enddo
  enddo
  close(hujisbm_unit1)
end if
call wrf_dm_bcast_bytes ( pkij , size ( pkij ) * RWORDSIZE )
call wrf_dm_bcast_bytes ( qkj  , size ( qkj  ) * RWORDSIZE )

call wrf_message('******* end of reading files *******')

!------------------------------------------------------------------------------
!---- cgs to mks
 cap_hydro =  cap_hydro * 1.0e-2  ! cm -> m
mass_hydro = mass_hydro * 1.0e-3  ! g -> kg
   bulkvt  =    bulkvt  * 1.0e-2  ! cm s-1 -> m s-1
   bulkrho =    bulkrho * 1.0e+3  ! g cm-3 -> kg m-3
   r_hydro =    r_hydro * 1.0e-2  ! cm -> m
      pkij =       pkij * 1.0e-15 ! g3 cm3 s-1 -> kg3 m3 s-1
       qkj =        qkj * 1.0e-15 ! g3 cm3 s-1 -> kg3 m3 s-1

!---- initial aerosol distribution
if(.not.restart) then

  do i=its, ite
  do j=jts, jte
  do k=kts, kte

    island  = xland(i,j)
    ivgtype = ivgtyp(i,j)

    if(k.lt.kte) then
      alt = (hgt(i,k+1,j)+hgt(i,k,j))/2.0
    else
      alt = hgt(i,k,j)
    end if

    call fccn_init(alt, island, ivgtype, isurban, fccni) 
    scalar(i,k,j,1+nhydro*nkr+1:1+nhydro*nkr+nkr) = fccni
 
  end do
  end do
  end do

  if(wrf_dm_on_monitor()) then
    write(6,*) "initial ccn", xland(its,jts),ivgtyp(its,jts),isurban
    write(6,*) scalar(its,kts,jts,1+nhydro*nkr+1:1+(nhydro+1)*nkr)
    write(6,*) "initial height", kts, kte
    write(6,*) hgt(its,kts:kte,jts)
  end if

end if

!---- collision kernel
do m1=2, nhydro
do m2=1, m1-1
  do i=1, nkr
  do j=1, nkr
    ywxx(m1,m2,:,i,j) = ywxx(m2,m1,:,j,i)
  end do
  end do
end do
end do

ywxx = ywxx * 1.0e-6 * dlnr * dt ! cm-3 s-1 --> m-3 s-1

!---- turbulence-induced collision enhancement for water-water lookup table
do i=1, 21
  rl_turb(i)=real(i)*1.0e-6
end do

do k=1, 11
  do i=1, 13
    x=r_hydro(WATER,i)
    if(x.lt.rl_turb( 1)) x=rl_turb( 1)
    if(x.gt.rl_turb(21)) x=rl_turb(21)
    do j=1, 13
      y=r_hydro(WATER,j)
      if(y.lt.rl_turb( 1)) y=rl_turb( 1)
      if(y.gt.rl_turb(21)) y=rl_turb(21)
      turbll_grid(i,j,k)=f(x,y,rl_turb,rl_turb,turbll_text1(:,:,k),21,21)
      turbll_grid(i,j,k)=turbll_grid(i,j,k)*turbll_text2(i,j,k)
    end do
  end do
end do

!---- turbulence-induced collision enhancement for water-graupel lookup table
! in micro m
rl_gl( 1)=r_hydro(WATER,1)*1.e6
rl_gl( 2)=  8.0; rl_gl( 3)= 10.0; rl_gl( 4)= 16.0; rl_gl( 5)= 20.0
rl_gl( 6)= 30.0; rl_gl( 7)= 40.0; rl_gl( 8)= 50.0; rl_gl( 9)= 60.0
rl_gl(10)= 70.0; rl_gl(11)= 80.0; rl_gl(12)= 90.0; rl_gl(13)=100.0
rl_gl(14)=200.0; rl_gl(15)=300.0
rl_gl(16)=r_hydro(WATER,24)*1.e6

rg_gl( 1)=r_hydro(GRAUPEL,1)*1.e6
rg_gl( 2)=  30.0; rg_gl( 3)=  60.0; rg_gl( 4)= 100.0; rg_gl( 5)= 200.0
rg_gl( 6)= 300.0; rg_gl( 7)= 400.0; rg_gl( 8)= 500.0; rg_gl( 9)= 600.0
rg_gl(10)= 700.0; rg_gl(11)= 800.0; rg_gl(12)= 900.0; rg_gl(13)=1000.0
rg_gl(14)=2000.0; rg_gl(15)=3000.0
rg_gl(16)=r_hydro(GRAUPEL,33)*1.e6

do i=krming_gl,krmaxg_gl
  x=r_hydro(GRAUPEL,i)*1.0e6
  if(x.lt.rg_gl(     1)) x=rg_gl(     1)
  if(x.gt.rg_gl(k0g_gl)) x=rg_gl(k0g_gl)
  do j=krminl_gl,krmaxl_gl
    y=r_hydro(WATER,j)*1.0e6
    if(y.lt.rl_gl(     1)) y=rl_gl(     1)
    if(y.gt.rl_gl(k0l_gl)) y=rl_gl(k0l_gl)
    turbgl_grid_400(i,j)=max(f(x,y,rg_gl,rl_gl,turbgl_text_400,k0g_gl,k0l_gl)*0.5,1.0)
    turbgl_grid_800(i,j)=max(f(x,y,rg_gl,rl_gl,turbgl_text_800,k0g_gl,k0l_gl)*0.5,1.0)
  enddo
enddo

!---- collision table
coll_table(1:num_coll_type_lt0,:) =     &
  transpose(reshape([                   &
    WATER, SNOW, GRAUPEL, SNOW,         &
    WATER, GRAUPEL, GRAUPEL, GRAUPEL,   &
    WATER, HAIL, HAIL, HAIL,            &
    COLUMN, WATER, COLUMN, GRAUPEL,     &
    PLATE, WATER, PLATE, GRAUPEL,       &
    DENDRITE, WATER, DENDRITE, GRAUPEL, &
    COLUMN, COLUMN, SNOW, SNOW,         &
    COLUMN, PLATE, SNOW, SNOW,          &
    COLUMN, DENDRITE, SNOW, SNOW,       &
    COLUMN, SNOW, SNOW, SNOW,           &
    PLATE, PLATE, SNOW, SNOW,           &
    PLATE, DENDRITE, SNOW, SNOW,        &
    PLATE, SNOW, SNOW, SNOW,            &
    DENDRITE, DENDRITE, SNOW, SNOW,     &
    DENDRITE, SNOW, SNOW, SNOW,         &
    SNOW, SNOW, SNOW, SNOW,             &
    SNOW, GRAUPEL, SNOW, GRAUPEL      ],&
    [4,num_coll_type_lt0]))
coll_table(num_coll_type_lt0+1:num_coll_type_lt0+num_coll_type_ge0,:) = &
  transpose(reshape([                   &
    WATER, SNOW, GRAUPEL, SNOW,         &
    WATER, GRAUPEL, GRAUPEL, GRAUPEL,   &
    WATER, HAIL, HAIL, HAIL,            &
    SNOW, SNOW, SNOW, SNOW,             &
    SNOW, GRAUPEL, SNOW, GRAUPEL      ],&
    [4,num_coll_type_ge0]))

!---- courant number for bott collision
do i=1, nkr
do j=i, nkr
  mass0 = mass_hydro(i)+mass_hydro(j)
  do k=1, nkr-1
    if(mass0.ge.mass_hydro(k).and.mass0.lt.mass_hydro(k+1)) exit
  end do
  ima(i,j)=k
  if(k.lt.nkr) then
    cour(i,j)=log(mass0/mass_hydro(k))/log(mass_hydro(k+1)/mass_hydro(k))
  end if

  cour(j,i)=cour(i,j)
  ima(j,i)=ima(i,j)
enddo
enddo

! breakup collision table
do i=1, nkr
do j=1, nkr
  ywxx(WATER,WATER,:,i,j) = ywxx(WATER,WATER,:,i,j) &
                          * ecoalmass(mass_hydro(i),mass_hydro(j))
end do
end do

return

 900 format(6e13.5)
2070 continue

call wrf_error_fatal("error opening hujisbm_data")

end subroutine sbminit

!----------------------------------------------------------------------------
subroutine fccn_init(alt, island, ivgtype, isurban, fccni)
implicit none

real, intent(in) :: alt, island
integer, intent(in) :: isurban, ivgtype
real, dimension(nkr), intent(out) :: fccni

real :: a, b, a1, n0ccn, kccn, r0, s_kr, factor
integer :: kr

a=akoe/(t00+15.0)
b=bkoe
a1=2.0*(a/3.0)**1.5/sqrt(b)

if(alt.le.z0_ccn) then
  factor = 1.0
else
  factor = exp(-(alt-z0_ccn)/dz_ccn)
end if

if(island.gt.1.5) then ! ocean
  n0ccn = n0ccn2
  kccn  = kccn2
else if(ivgtype.eq.isurban) then ! urban
  n0ccn = n0ccn1 * urb_ccnf
  kccn  = kccn1
else ! rural
  n0ccn = n0ccn1
  kccn  = kccn1
end if

r0=r_hydro(WATER,1)*100. ! m -> cm

do kr=nkr, 1, -1
  s_kr=a1/(r0**1.5)
  if(r0.le.rccn_min.or.s_kr.ge.smax) then
    fccni(kr)=0.0
  else
    fccni(kr)=1.5*n0ccn*kccn*(100.*s_kr)**kccn*dlnr*factor
  end if
  r0=r0/exp(dlnr)
end do

end subroutine fccn_init
!---------------------------------!-------------------------------------------

real function f(x,y,x0,y0,table,k0,kk0)
implicit none

real,    intent(in) :: x, y
integer, intent(in) :: k0, kk0
real,    intent(in), dimension(k0    ) :: x0
real,    intent(in), dimension(kk0   ) :: y0
real,    intent(in), dimension(k0,kk0) :: table

integer :: k, ir, kk, iq
real :: p, q, ec, ek

do k=2,k0
  if(x.ge.x0(k-1).and.x.lt.x0(k)) then
    ir=k
  elseif(x.ge.x0(k0)) then
    ir=k0+1
  elseif(x.lt.x0(1)) then
    ir=1
  endif
enddo

do kk=2,kk0
  if(y.ge.y0(kk-1).and.y.lt.y0(kk)) then
    iq=kk
   elseif(y.ge.y0(kk0)) then
    iq=kk0+1
  elseif(y.lt.y0(1)) then
    iq=1
  endif
enddo

if(ir.lt.k0+1) then
  if(ir.ge.2) then
    p =(x-x0(ir-1))/(x0(ir)-x0(ir-1))
    q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
    ec=(1.d0-p)*(1.d0-q)*table(ir-1,iq-1)+ &
             p *(1.d0-q)*table(ir,  iq-1)+ &
             q *(1.d0-p)*table(ir-1,iq  )+ &
             p *      q *table(ir,  iq  )
  else
    q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
    ec=(1.d0-q)*table(1,iq-1)+q*table(1,iq)
  endif
else
  q =(y-y0(iq-1))/(y0(iq)-y0(iq-1))
  ec=(1.d0-q)*table(k0,iq-1)+q*table(k0,iq)
endif

f=ec
end function f

!----------------------------------------------------------------------------
subroutine freeze_lwf(dt, tt, pp, qq, rho_a,          &
                      mass, vt_w, vt_p, rho_p, ff, lwf)
implicit none

real, intent(in ) :: dt, pp, qq, rho_a, mass, vt_w, rho_p, ff
real, intent(inout) :: tt, lwf, vt_p

real :: e, Dv, Ka, etaa, Nsc, Npr, Nre
real :: fm_w, fm_i, v_i, v_w, v_p, r_p, r_w, vt_intp, fh, fv, Ts, delm

e  = pp*qq/(eps+qq)                      ! water vapor pressure
Dv = 2.11e-5*((tt/t00)**1.94)*(p00/pp)   ! diffusivity of water vapor
Ka = 4.184e-3*(5.69+0.017*(tt-t00))      ! heat conductivity of air
etaa = (1.718+0.0049*(tt-t00)-1.2e-5*(tt-t00)**2)*1.0e-5 ! dynamic viscosity of air (Physics and Dynamics of Clouds and Precipitation)
Nsc = etaa/(rho_a*Dv) ! Schmidt number
Npr = etaa/(Ka/cp)    ! Prandtl number

! start of freezing
fm_w = mass*lwf
fm_i = mass*(1.0-lwf)
v_i = fm_i/rho_p
v_w = fm_w/rho_water
v_p = v_i + v_w

r_w = (v_w/(4.0/3.0*pi))**(1.0/3.0)
r_p = (v_p/(4.0/3.0*pi))**(1.0/3.0)

vt_intp = vt_p*(1.0-lwf)*vt_w*lwf

Nre = vt_intp*2.0*rho_a*r_p/etaa ! Reynolds number
fh = ventilation_coef(Nre,Npr,HAIL)
fv = ventilation_coef(Nre,Nsc,HAIL)

! Ts and delm is given in Phillips et al. (JAS 2007)
!Ts   = surface_temp &
!       (e, tt, (Ls*Dv)/(Ka*Rv), fv/fh, (Ki*r_w)/(Ka*(r_p-r_w))/fh, t00)
Ts   = surface_temp &
       (e, tt, (Ls*Dv)/(Ka*Rv), fv/fh, 0.0, t00)
delm = dt*(4.0*pi*r_p/Lm)*(fh*Ka*(t00-tt)+(Dv*Ls/Rv)*fv*(esw(t00)/t00-e/tt))

! calculate tendency for temperature and lwf
if(delm.gt.0) then
  delm = min(delm, mass*lwf)
  delm = min(delm, (t00-tt)/ff/Lm*rho_a*Cp)
  tt  = tt  + delm*ff*Lm/rho_a/Cp
  lwf = lwf - delm/mass
end if

vt_p = vt_intp

end subroutine freeze_lwf

!===================================================================!
! freezing scheme (Vali, 1975; Biggs, 1953)                         !
!===================================================================!
subroutine freeze_drop(tin, rho, mass, ff1, ff_s, ff_l, lwf_l, dt, dtemp)
use ieee_arithmetic
implicit none

real, intent(in) :: rho, mass, dt, dtemp
real, intent(inout) :: tin, ff1, ff_s, ff_l, lwf_l

real :: dn, fliq_frac, liqm_bf, liqm_add, factor, tt_temp

if(inan.eq.1) then

if(ieee_is_nan(tin)) then
  write(6,*) "before freeze"
  write(6,*) tin, rho, mass
  write(6,*) dt, dtemp
  write(6,*) ff1, ff_s, ff_l, lwf_l
  write(6,*) dn, fliq_frac, liqm_bf, liqm_add
  call wrf_error_fatal("nan detected before freeze")
end if

end if

if(tin.ge.t00) then
  return
else if(tin-t00.gt.-30.0) then
  dn = min(max((-n0_frz*0.1*gam*(0.1*(t00-tin))**(gam-1.0))*dtemp,0.0),ff1)
else if(tin-t00.gt.-38.0) then
  dn = ff1*(1.0-exp(-mass*a_frz*exp(-b_frz*(tin-t00))*dt))
else
  dn = ff1
end if

! thermodynamic limits
if(mass.lt.mass_frz_bnd) then
  fliq_frac = 0.
else
  fliq_frac = fmin_liq_frac(tin-t00)
end if

tt_temp = tin + (Lm*dn*mass*(1.0-fliq_frac))/rho/Cp

if(tt_temp.gt.t00) then
  dn = (t00-tin)*rho*Cp/Lm/mass/(1.0-fliq_frac)
  tin = t00
else
  tin = tt_temp
end if

! change in number
ff1=ff1-dn

if(mass.lt.mass_frz_bnd) then
  ff_s      = ff_s+dn
else
  liqm_bf   = ff_l*lwf_l
  liqm_add  = dn  *fliq_frac
  ff_l      = ff_l+dn
  lwf_l     = max(min((liqm_bf+liqm_add)/ff_l,1.0),0.0)
endif

if(inan.eq.1) then
  if(ieee_is_nan(tin)) then
    write(6,*) "after freeze"
    write(6,*) tin, rho, mass
    write(6,*) dt, dtemp
    write(6,*) ff1, ff_s, ff_l, lwf_l
    write(6,*) dn, fliq_frac, liqm_bf, liqm_add
    call wrf_error_fatal("nan detected after freeze")
  end if
end if

end subroutine freeze_drop

!-------------------------------------------------------------------------
real function fmin_liq_frac(ttin_c)
implicit none

real, intent(in) :: ttin_c ! temperature in Celsius

real, parameter :: cw_zero = 1.000938
real, parameter :: c_i     = 0.503

real :: ttin, cw

ttin = max(ttin_c, -37.0)

if(ttin.lt.0.0) then
  cw = cw_zero -2.7052e-3*ttin -2.3235e-5*(ttin**2) +4.3778e-6*(ttin**3)
  fmin_liq_frac = min(max(1.0+((cw+cw_zero)/2.0)/((Lm/4184./ttin_c)+c_i),0.),1.)
else
  fmin_liq_frac = 1.0
endif

end function fmin_liq_frac
!----------------------------------------------------------------------------
subroutine melting(dt, tt, pp, qq, rho_a,                  &
                   mass, vt_w, vt_i, rho_i, ff, lwf, type, &
                   Nre, m_shed                             )
use ieee_arithmetic
implicit none

real, intent(in ) :: dt, pp, rho_a, mass, vt_w, rho_i, ff
real, intent(inout) :: tt, qq, lwf, vt_i
real, intent(out) :: Nre, m_shed
integer, intent(in) :: type

real :: ee, Dv, Ka, etaa, Nsc, Npr
real :: fm_i, fm_w, v_i, v_w, fh, fv, Ts, delm, delm_subl, delm_evap
real :: vt_intp, rho_intp, rho1, ar_i, ar_izero, cap_izero, cap, fl_star
real :: a_i, a_d, a_izero, a_crit, v_soakable
real :: Nre_dry, Nre_smooth, Nre_shed, Nre_fast, X_Best, X_Best_crit
real :: fm_w_soaked, fm_w_crit, frac_eqm, vt_start, vt_end

ee = pp*qq/(eps+qq)                      ! water vapor pressure
Dv = 2.11e-5*((tt/t00)**1.94)*(p00/pp)   ! diffusivity of water vapor
Ka = 4.184e-3*(5.69+0.017*(tt-t00))      ! heat conductivity of air
etaa = (1.718+0.0049*(tt-t00)-1.2e-5*(tt-t00)**2)*1.0e-5 ! dynamic viscosity of air (Physics and 
Nsc = etaa/(rho_a*Dv) ! Schmidt number
Npr = etaa/(Ka/cp)    ! Prandtl number

select case(type)

case(COLUMN,PLATE,DENDRITE,SNOW)

  rho_intp  = rho_i + lwf*(rho_water - rho_i)
  fm_w = mass*lwf
  fm_i = mass*(1.0-lwf)
  v_i = mass/rho_intp
  rho1 = fm_i/v_i

  ar_izero  = ar(mass, rho_i, type)
  ar_i      = ar_izero + lwf*(1.0-ar_izero)
  call cap_zero(mass, ar_i, rho1, cap_izero, fl_star, type)
  cap       = cap_izero*(0.8 + lwf*0.2)

  vt_intp   = vt_i*(1-chi_fra(lwf)) + vt_w*chi_fra(lwf)
  Nre       = vt_intp*fl_star*rho_a/etaa
  fv        = ventilation_coef(Nre, Nsc, type, arx=ar_i)
  fh        = fv

case(GRAUPEL,HAIL)

  fm_w = mass*lwf
  fm_i = mass*(1.0-lwf)
  v_i  = fm_i/rho_i
  v_w  = fm_w/rho_water

  if(rho_i < rho_crit) then
    v_soakable = v_i - fm_i/rho_ice
  else
    v_soakable = 0.0
  endif

  a_i     = (fm_i/(4.0/3.0*pi*rho_i))**(1.0/3.0)
  a_izero = (mass/(4.0/3.0*pi*rho_i))**(1.0/3.0)
  Nre_dry = vt_i*2.0*rho_a*a_izero/etaa
  X_Best  = 8.0*mass*rho_a*grav/(pi*etaa**2)
  Nre_smooth = fnre_sphere(X_Best)

  if(v_w < v_soakable) then

    a_d = a_i

    if(Nre_dry.lt.4000) then
      Nre = Nre_smooth
    else
      Nre = Nre_dry
    end if
    vt_intp = Nre/(2.0*a_d*rho_a)*etaa

  else

    a_d = ((v_i + (v_w - v_soakable))/(4.0/3.0*pi))**(1.0/3.0)
    fm_w_soaked = rho_water*v_soakable
    fm_w_crit = (0.268e-3 + (fm_i + fm_w_soaked)*0.1389)
    a_crit = ((v_i + fm_w_crit/rho_water)/(4.0/3.0*pi))**(1.0/3.0)
    Nre_shed = 4800.0 + 4831.5*(fm_i + fm_w_soaked)*1000.0
    frac_eqm  = min(max((fm_w-fm_w_soaked)/fm_w_crit,0.0),1.0)

    if(Nre_dry.lt.4000..or.rho_i.lt.800.) then
      vt_start = Nre_smooth*etaa/(2.0*a_i*rho_a)
!      vt_start = vt_low_density_trans  &
!                 (Nre_dry,Nre_smooth,vt_p,a_izero,etaa,rho_a,rho_p,mass)
    else
      vt_start = vt_i*a_izero/a_i
!      vt_start = vt_high_density_trans &
!                 (Nre_dry,Nre_smooth,vt_p,a_izero,etaa,rho_a)
    endif

    if(Nre_shed.lt.5000.) then
      vt_end = vt_w
    else if(Nre_shed.lt.25000.) then
      vt_end = 1.5e-5*Nre_shed/(2.0*a_crit)*sqrt(1.2*rho_a)
    else
      X_Best_crit = 8.0*(fm_i +fm_w_soaked + fm_w_crit)*rho_a*grav/(pi*etaa**2)
      Nre_fast    = sqrt(X_Best_crit/0.6)
      vt_end      = Nre_fast*etaa/(2.0*a_crit*rho_a)
    end if
    
!    vt_end  = vt_eqm &
!     (fm_i+fm_w_soaked,fm_w_crit,mass_hydro,vt_hydro(WATER,:),rho_a,etaa,a_crit)

    vt_intp = vt_start + (vt_end - vt_start)*frac_eqm
    Nre = vt_intp*(2.0*a_d*rho_a)/etaa

  endif

  fv = ventilation_coef(Nre,Nsc,HAIL)
  fh = ventilation_coef(Nre,Npr,HAIL)

  if(Nre < 6000.0) then
    cap = a_d
  else
    cap = a_i
  endif

end select

!Ts = surface_temp(ee, tt, (Ls*Dv)/(Ka*Rv), fv/fh, 0.0, t00)
!if(Ts.ge.t00) then
!  Ts = t00
!else
!  Nre = 0.0
!  m_shed = 0.0
!  return
!end if

! sublimation or evaporation
!if(lwf.eq.0.0) then
!  delm_subl = dt*(4.0*pi*cap*fv)*(ee/esi(Ts)-1.0)/((Ls/Rv/tt-1.0)*(Ls/Ka/tt)+Rv*tt/esi(Ts)/Dv)
!  qq = qq - delm_subl*ff/rho_a
!  tt = tt + delm_subl*ff/rho_a*Ls/Cp
!else
!  delm_evap = dt*(4.0*pi*cap*fv)*(Dv/Rv)*(ee/tt-esw(t00)/t00)
!  qq = qq + delm_evap*ff/rho_a
!  tt = tt + delm_evap*ff/rho_a*Lv/Cp
!end if

delm = dt*(4.0*pi*cap/Lm)*(Ka*(tt-t00)*fh+(Dv*Lv/Rv)*fv*(ee/tt-esw(t00)/t00))

! update temperature, lwf, and terminal velocity
if(delm.gt.0.0) then
  delm = min(delm, mass*(1.-lwf))
  delm = min(delm, (tt-t00)/ff/Lm*rho_a*Cp)
  tt   = tt  - delm*ff*Lm/rho_a/Cp
  lwf  = lwf + delm/mass
end if

vt_i = vt_intp

if(lwf.lt.-0.001.or.lwf.gt.1.001) then
  write(6,*) "during melting", lwf
  write(6,*) tt, delm, ff, mass
!  lwf = 0.0
end if

! diagnose whether shedding is needed or not
if(type.eq.GRAUPEL.or.type.eq.HAIL) then
  fm_i = mass*(1.0-lwf)
  fm_w = mass*lwf
  v_i  = fm_i/rho_i

  if(rho_i.gt.rho_crit) then
    v_soakable = 0.0
  else
    v_soakable = v_i - fm_i/rho_i
  end if
  fm_w_soaked = rho_water*v_soakable
  fm_w_crit   = 0.268e-3 + 0.1389*(fm_i + fm_w_soaked)

  if(Nre.gt.25000.) then
    m_shed = max(fm_w - (fm_w_soaked + fm_w_crit),0.0)
  else
    m_shed = max(fm_w - fm_w_soaked,0.0)
  end if
else
  m_shed = 0.0
end if

if(inan.eq.1) then
  if(ieee_is_nan(tt)) then
    write(6,*) "in melting", type
    write(6,*) rho_a, rho_i, mass, delm
    write(6,*) ar_izero, ar_i, lwf, fm_i
    write(6,*) cap_izero, fl_star, cap, Nre
    write(6,*) fv, Nsc, etaa, vt_i
    call wrf_error_fatal("nan in tt")
  end if
end if

end subroutine melting

!-------------------------------------------------------------------
subroutine shedding(Nre, m_shed, ff, lwf, mass, ff_w, ff_i, lwf_i, kr)
implicit none

integer, intent(in) :: kr
real, intent(in) :: Nre, ff, lwf, mass
real, intent(inout) :: m_shed
real, dimension(nkr), intent(inout) :: ff_w, ff_i, lwf_i

integer :: i
real :: d_drop, m_drop, n_drop, frac
double precision :: lw, lw1, m_ice, lwf_n

if(Nre.gt.15000.) then
  d_drop = 1.5e-3
else if(Nre.gt.10000.) then
  d_drop = 3.0e-3
else
  d_drop = 4.5e-3
end if

m_shed = min(m_shed, mass*lwf)

m_drop = min(m_shed*0.5,rho_water*pi/6.0*(d_drop**3))
n_drop = (m_shed/m_drop)*ff

if(m_drop.lt.mass_hydro(1)) then
  ff_w(1) = ff_w(1) + n_drop*m_drop/mass_hydro(1)
else
  do i=1, nkr-1
    if(m_drop.ge.mass_hydro(i).and.m_drop.lt.mass_hydro(i+1)) exit
  end do
  if(i.eq.nkr) then
    ff_w(nkr) = ff_w(nkr) + n_drop*m_drop/mass_hydro(nkr)
  else
    ff_w(i)   = ff_w(i)   + n_drop*(mass_hydro(i+1)-m_drop)/(mass_hydro(i+1)-mass_hydro(i))
    ff_w(i+1) = ff_w(i+1) + n_drop*(m_drop-mass_hydro(i)  )/(mass_hydro(i+1)-mass_hydro(i))
  end if
end if

m_ice = mass - m_shed
lw    = (mass*lwf-m_shed)*ff
lwf_n = lwf-m_shed/mass

 ff_i(kr) = 0.0
lwf_i(kr) = 0.0

if(m_ice.lt.mass_hydro(1)/10.0) then
  return
else if(m_ice.lt.mass_hydro(1)) then
  lw       = lw      + mass_hydro(1)*lwf_i(1)*ff_i(1)
  ff_i(1)  = ff_i(1) + ff*(m_ice/mass_hydro(1))
  lwf_i(1) = lw/mass_hydro(1)/ff_i(1)
else
  do i=1, kr-2
    if(m_ice.ge.mass_hydro(i).and.m_ice.lt.mass_hydro(i+1)) exit
  end do
  if(i.eq.nkr) then
    lw         = lw        + mass_hydro(nkr)*lwf_i(nkr)*ff_i(nkr)
    ff_i(nkr)  = ff_i(nkr) + ff*(m_ice/mass_hydro(nkr))
    lwf_i(nkr) = lw/mass_hydro(nkr)/ff_i(nkr)
  else
    frac     = (mass_hydro(i+1)-m_ice)/(mass_hydro(i+1)-mass_hydro(i))
    frac     = max(min(frac,1.0),0.0)
    lw1      = frac*(lw/m_ice*mass_hydro(i)) + mass_hydro(i)*lwf_i(i)*ff_i(i)
    ff_i(i)  = ff_i(i)  + ff*frac
    lwf_i(i) = lw1/mass_hydro(i)/ff_i(i)

    frac       = (m_ice-mass_hydro(i))/(mass_hydro(i+1)-mass_hydro(i))
    frac       = max(min(frac,1.0),0.0)
    lw1        = frac*(lw/m_ice*mass_hydro(i+1)) + mass_hydro(i+1)*lwf_i(i+1)*ff_i(i+1)
    ff_i(i+1)  = ff_i(i+1) + ff*frac
    lwf_i(i+1) = lw1/mass_hydro(i+1)/ff_i(i+1)
  end if
end if

end subroutine shedding

!----------------------------------------------------------------------------
real function esi(tt)
implicit none

real, intent(in) :: tt

esi = 611.*exp(a_i - b_i/tt)

end function esi

!----------------------------------------------------------------------------
real function esw(tt)
implicit none

real, intent(in) :: tt

esw = 611.*exp(a_w - b_w/tt)

end function esw

!----------------------------------------------------------------------------
real function ventilation_coef(fnrex, fnumber, type, arx)
implicit none

integer, intent(in) :: type
real,    intent(in) :: fnrex, fnumber
real,    intent(in), optional :: arx

real :: x_f, chi_fr, fx

x_f = (fnrex**(1.0/2.0))*(fnumber**(1.0/3.0))

select case (type)
  case(COLUMN)
    fx = 1.0                      &
       - 6.68e-3*  x_f/4.0        &
       + 2.39402*((x_f/4.0)**2.0) &
       + 0.73409*((x_f/4.0)**3.0) &
       - 0.73911*((x_f/4.0)**4.0)

  case(PLATE)
    fx = 1.0                       &
       - 0.06042*  x_f/10.0        &
       + 2.79820*((x_f/10.0)**2.0) &
       - 0.31933*((x_f/10.0)**3.0) &
       - 0.06247*((x_f/10.0)**4.0)

  case(DENDRITE)
    fx = 1.0                       &
       + 0.35463*  x_f/10.0        &
       + 3.55338*((x_f/10.0)**2.0)

  case(SNOW)
    if(.not.present(arx)) then
      write(6,*) "AR is needed when calculating snow vent coef."
    end if

    ! real snow is not spherical, so this should not be used
    if(abs(arx-1.0).lt.1.e-2) then
      if(x_f < 1.4) then
        fx = 1.0  + 0.108*x_f**2
      else
        fx = 0.78 + 0.308*x_f
      endif
    ! this is the correct formula for real snow
    else
      if(x_f <= 1.0) then
        fx = 1.0  + 0.14*x_f**2
      else
        fx = 0.86 + 0.28*x_f
      endif
    endif

  case(HAIL)
    if(fnrex < 6000.0) then
      if(x_f < 1.4) then
        fx = 1.0  + 0.108*x_f**2
      else
        fx = 0.78 + 0.308*x_f
      endif
      if(fnrex < 250.0) then
        fx = fx*2
      endif
    else
      if(fnrex < 20000.0) then
        chi_fr = 0.76
      else
        chi_fr = 0.57 + fnrex*9.e-6
      endif
      fx = chi_fr*x_f/2.0
    endif

  case default
    write(6,*) "somethig wrong in call ventilation_coef."
end select

fx = max(min(fx,150.0),1.0)

ventilation_coef = fx

end function ventilation_coef

!-----------------------------------------------------------------------------
real function surface_temp(eex, tempk, f_vap, fvofh, f_ice, t_liq_k)
implicit none

real, intent(in)::  eex, tempk, f_vap, fvofh, f_ice, t_liq_k

real :: tsx, tsxold, tsxnew, tdiff, f_tsx, f_tsxold
integer :: ilj

tsxold = 269.0
tsx = 270.0
tdiff = 1.0

ilj = 0

do while(tdiff > 1.0e-2)
  ilj = ilj + 1

  f_tsx   = (tempk-f_vap*fvofh*(esi(tsx)   /tsx   -eex/tempk+f_ice*t_liq_k)) &
          / (1.0 + f_ice) - tsx

  f_tsxold= (tempk-f_vap*fvofh*(esi(tsxold)/tsxold-eex/tempk+f_ice*t_liq_k)) &
          / (1.0 + f_ice) - tsxold

  tsxnew  = tsx - f_tsx*(tsx - tsxold)/(f_tsx - f_tsxold)

  tsxold = tsx
  tsx = tsxnew

  tdiff = abs(tsx-tsxold)

  if(ilj > 1.0e6) then
   write(6,*) 'not converging', eex, tempk, fvofh, f_ice
   exit
  endif
enddo

surface_temp = tsx

end function surface_temp

!-----------------------------------------------------------------------------
real function ar(fmassx, rhoix, type)
implicit none

real, intent(in) :: fmassx, rhoix
integer, intent(in) :: type

real :: a, b, m, d_equiv, shape, fl_i, d_i, h_i

select case (type)

case (COLUMN)

  m = 4.0*fmassx/pi/rhoix*1.0e9  ! unit of a: m3 -> mm3

  if(a.lt.0.002) then
    a = 0.5
    b = 1.0
  else
    a = 0.1973
    b = 0.414
  end if

  h_i = (m/(a**2))**(1.0/(2.0*b+1.0))
  d_i = a*h_i**b

  ar  = h_i/d_i 
  ar  = max(min(ar,5.0),2.0)

case (PLATE)

  d_i = ((fmassx*1.e-3)/0.0376)**(1.0/3.31)
  h_i = 0.0141*(d_i**0.474)

  ar = h_i/d_i

case (DENDRITE)

  d_i = ((fmassx*1.e-3)/0.00376)**(1.0/2.79)
  h_i = 0.00996*(d_i**0.415)

  ar = h_i/d_i

case (SNOW)

  ar = 0.3

end select

end function ar

!----------------------------------------------------------------------
subroutine cap_zero(fm_ice, ar_ice, rho_i, cap_izero, flstar, type)
implicit none

integer, intent(in) :: type
real, intent(in) :: fm_ice, ar_ice, rho_i
real, intent(out) :: cap_izero, flstar

real :: a_i, b_i, eps_i, omega_i, p_i

a_i   = (fm_ice/(4.0/3.0*pi*rho_i*ar_ice))**(1.0/3.0)

select case (type)

case (COLUMN)

  b_i   = ar_ice*a_i
  eps_i = sqrt(b_i**2 - a_i**2)
  cap_izero = eps_i/log((b_i+eps_i)/a_i)

  omega_i = 2.0*(pi*a_i**2) + 4.0*b_i*a_i
  p_i = 2.0*pi*a_i
  flstar = omega_i/p_i

case (PLATE,DENDRITE,SNOW)

  eps_i = sqrt(1.0 - ar_ice**2)

  if(eps_i.gt.0.001) then
    cap_izero = a_i*eps_i/asin(eps_i)

    omega_i = pi*a_i**2*(2.0+(ar_ice**2)/eps_i*log((1.0+eps_i)/(1.0-eps_i)))
    p_i = 2.0*pi*a_i
    flstar = omega_i/p_i
  else
    cap_izero = a_i
    flstar = 2.0*a_i
  endif

end select

end subroutine cap_zero

!-----------------------------------------------------------------------------
real function chi_fra(fra)
implicit none

real, intent(in) :: fra

integer, parameter :: ix_max = 14
real, dimension(ix_max) :: xxa, yya

integer :: ix

if(fra.le.0.0) then
  chi_fra = 0.0
  return
else if(fra.ge.1.0) then
  chi_fra = 1.0
  return
endif

xxa = (/0.,    0.1,   0.2,   0.3,   0.4,   0.5,   0.6,   &
        0.7,   0.75,  0.8,   0.85,  0.9,   0.95,  1.    /)
yya = (/0.,    0.0125,0.0312,0.05,  0.0812,0.1187,0.1749,&
        0.2436,0.2873,0.3498,0.4372,0.5684,0.7308,1.    /)

do ix=1, ix_max
  if(xxa(ix).gt.fra) then
    chi_fra=yya(ix-1)+(fra-xxa(ix-1))*(yya(ix)-yya(ix-1))/(xxa(ix)-xxa(ix-1))
    return
  end if
end do

write(6,*) "Something wrong in chi_fra"

end function chi_fra

!-----------------------------------------------------------------------
real function fnre_sphere(xd)
implicit none

real, intent(in) :: xd

real :: ww1, ww2, ww3

ww1 = log10(xd)
ww2 = ww1**2
ww3 = ww1**3

if(xd < 73.0) then
  fnre_sphere = xd/24.0
else if(xd < 562.0) then
  fnre_sphere = 10.0**(-1.7095+1.33438*ww1-0.11591*ww2)
else if(xd < 1830.) then
  fnre_sphere = 10.0**(-1.81391+1.34671*ww1-0.12427*ww2+6.3e-3*ww3)
else if(xd < 5.4e10) then
  fnre_sphere = 10.0**(-1.706026+1.225713*ww1-0.08962*ww2+3.567e-3*ww3)
else
  fnre_sphere = (xd/0.1)**0.5
endif

end function fnre_sphere

!-----------------------------------------------------------------------
subroutine kernel_pressure(pp)
use ieee_arithmetic
implicit none

real, intent(in) :: pp

real :: p_1, p_2, p_3
integer :: i, j, t1, t2

do t1=1, nhydro
do t2=1, nhydro
  if(t1.eq.WATER.and.t2.eq.WATER) then
    p_1 = 100000.; p_2 = 75000.; p_3 = 50000.
  else
    p_1 =  75000.; p_2 = 50000.; p_3 = 30000.
  end if
 
  do i=1, nkr
  do j=1, nkr
    if     (pp.ge.p_1) then
      ckxx(t1,t2,i,j)=ywxx(t1,t2,1,i,j)
    else if(pp.ge.p_2) then
      ckxx(t1,t2,i,j)=ywxx(t1,t2,2,i,j)*(p_1-pp)/(p_1-p_2) &
                     +ywxx(t1,t2,1,i,j)*(pp-p_2)/(p_1-p_2)
    else if(pp.ge.p_3) then
      ckxx(t1,t2,i,j)=ywxx(t1,t2,3,i,j)*(p_2-pp)/(p_2-p_3) &
                     +ywxx(t1,t2,2,i,j)*(pp-p_3)/(p_2-p_3)
    else
      ckxx(t1,t2,i,j)=ywxx(t1,t2,3,i,j)
    end if

    if(inan.eq.1) then
      if(.not.ieee_is_finite(ckxx(t1,t2,i,j))) then
        write(6,*) pp
        call wrf_error_fatal("inf of ckxx in kernel_pressure")
      end if
    end if

  end do
  end do
end do
end do

end subroutine kernel_pressure

!=================================================================
subroutine kernel_mod(tt, pp, qq, tke, n_sq, dx, dy, dz, rho, eps1, rey1)
use ieee_arithmetic
implicit none

real, intent(in) :: tt, pp, qq, tke, n_sq, dx, dy, dz, rho
real, intent(out) :: eps1, rey1

real :: factor, r1_new, r2_new, a_ij, b_ij, c_ij
integer :: i, j, k, m1, m2, n1, n2

call tice(tke, n_sq, dx, dy, dz, rho, eps1, rey1)

call onecwss(tt, pp, qq)

! update to reflect density
do m1=1, nhydro
do m2=1, nhydro
do n1=1, nkr
do n2=1, nkr

  if(rho_hydro(m1,n1).ne.bulkrho(m1,n1).or.rho_hydro(m2,n2).ne.bulkrho(m2,n2)) then
    r1_new = r_hydro(m1,n1)*(bulkrho(m1,n1)/rho_hydro(m1,n1))**(1.0/3.0)
    r2_new = r_hydro(m2,n2)*(bulkrho(m2,n2)/rho_hydro(m2,n2))**(1.0/3.0)

    a_ij = r1_new + r2_new
    b_ij = r_hydro(m1,n1) + r_hydro(m2,n2)
    c_ij = (a_ij/b_ij)**2

    ckxx(m1,m2,n1,n2)=ckxx(m1,m2,n1,n2)*c_ij
  end if

  if(inan.eq.1) then
    if(.not.ieee_is_finite(ckxx(m1,m2,n1,n2))) then
      write(6,*) a_ij, b_ij, c_ij
      call wrf_error_fatal("inf of ckxx in kernel_mod (density)")
    end if
  end if

end do
end do
end do
end do

! update to reflect terminal velocity
do m1=1, nhydro
do m2=1, nhydro
do n1=1, nkr
do n2=1, nkr

  if(vt_hydro(m1,n1).ne.vt_ref(m1,n1).or.vt_hydro(m2,n2).ne.vt_ref(m2,n2)) then
    if(abs(vt_ref(m1,n1)-vt_ref(m2,n2)).gt.1.0e-4) then
      factor = min(abs((vt_hydro(m1,n1)-vt_hydro(m2,n2))     &
                      /(vt_ref  (m1,n1)-vt_ref  (m2,n2))),5.0)
      ckxx(m1,m2,n1,n2)=ckxx(m1,m2,n1,n2)*factor
    end if 
  end if

  if(inan.eq.1) then
    if(.not.ieee_is_finite(ckxx(m1,m2,n1,n2))) then
      write(6,*) factor
      call wrf_error_fatal("inf of ckxx in kernel_mod (term_vel)")
    end if
  end if

end do
end do
end do
end do

end subroutine kernel_mod

!===================================================================
subroutine ice_mpl(g_d, g_i, dnc_old, tt)
implicit none

real, dimension(nkr), intent(in) :: g_d
real, dimension(nkr), intent(inout) :: g_i
real, intent(in) :: dnc_old, tt

real :: dnc_new, d_dnc, d_inc

dnc_new = sum(g_d(kr_icempl:nkr)/mass_hydro(kr_icempl:nkr))*dlnr
d_dnc = dnc_old-dnc_new

if(tt.le.t_icempl0) then
  d_inc=d_dnc*4.0e-3*(t_icempl1-tt)/(t_icempl1-t_icempl0)
else
  d_inc=d_dnc*4.0e-3*(t_icempl2-tt)/(t_icempl2-t_icempl0)
endif

if(d_inc.lt.0.0) then
  write(6,*) "error in ice_mpl"
  write(6,*) d_inc, dnc_old, dnc_new, d_dnc, tt
end if
g_i(1)=g_i(1)+d_inc*mass_hydro(1)/dlnr

end subroutine ice_mpl

!====================================================================
subroutine onecwss(tt, pp, qq)
use ieee_arithmetic
implicit none

real, intent(in) :: tt, pp, qq

real, parameter :: tc_max = -20.0
real, parameter :: tc_min = tt_coal - t00

real :: tc, satq2, temp, qs2, dele, t, p, factor
integer :: i, j, m1, m2

satq2(t,p)=3.80e3*(10**(9.76421-2667.1/t))/p
temp(t)=4.5185186e-5*t**3+3.4793e-3*t**2+9.31878e-2*t+8.8333e-1

tc = tt-t00

if(tc.le.0.0) then
  qs2  = satq2(tt,pp)
  dele = temp(tc)*qq*(eps+(1.0-eps)*qs2)/(eps+(1.0-eps)*qq)/qs2

  if     (tc.ge.-6.0) then
    factor = min(max(dele,0.5),1.0)
  else if(tc.ge.-12.5) then
    factor = 0.5
  else if(tc.ge.-17.0) then
    factor = 1.0
  else if(tc.ge.tc_max) then
    factor = 0.4
  else if(tc.ge.tc_min) then
    factor = 0.25*(tc-tc_min)/(tc_max-tc_min)
  else
    factor = 0.0
  endif

  ckxx(COLUMN:GRAUPEL,COLUMN:GRAUPEL,1:nkr,1:nkr) = &
  ckxx(COLUMN:GRAUPEL,COLUMN:GRAUPEL,1:nkr,1:nkr)*factor
 
  if(inan.eq.1) then
    do m1=COLUMN, GRAUPEL
    do m2=COLUMN, GRAUPEL
    do i=1, nkr
    do j=1, nkr
   
      if(.not.ieee_is_finite(ckxx(m1,m2,i,j))) then
        write(6,*) m1, m2, i, j, factor
        call wrf_error_fatal("inf of ckxx in onecwss1")
      end if

    end do
    end do
    end do
    end do
  end if

end if

if(factor.lt.0.0) then
  write(6,*) "in onecwss"
  write(6,*) tt, pp, qq, factor
  write(6,*) qs2, dele, tc_min, tc_max
  call wrf_error_fatal("negative factor in onecwss")
end if

end subroutine onecwss

!====================================================================
subroutine tice(tke, n2, dx, dy, dz, rho, eps1, rey1)
use ieee_arithmetic
!
! Reference : ARW technical note ver.3
!
implicit none

real, intent(in) :: tke, n2, dx, dy, dz, rho
real, intent(out) :: eps1, rey1

real, parameter :: ck = 0.1
real, parameter :: visd = 1.81e-05 ! dynamic viscosity

integer :: i, j
real :: visk, length, c, epsi, lambda, reyi, n, vel
real, dimension(13,13) :: array_ll
real, dimension(krming_gl:krmaxg_gl,krminl_gl:krmaxl_gl) :: array_gl

real :: extra_tice, eps0, rey0, rst, c1, c2, c3, c4, &
        len0, tketmp, ce1, ce2, coefc

! extrapolation of kernel using only 1 point
extra_tice(rst,epsi,eps0,reyi,rey0) = rst*(epsi/eps0)**1.5*(reyi/rey0)**0.25

len0 = (dx*dy*dz)**(1.0/3.0)
visk = visd/rho  ! kinematic viscosity. unit of rho should be MKS
tketmp = max(tke,1.0e-6)

! l = min[(dx*dy*dz)**(1/3), 0.76*sqrt(e)/N]
! c = 1.9*ck + (0.93-1.9*ck)*l/ds

length = len0
if(n2.gt.1.0e-6) then
  length = 0.76*sqrt(tketmp/n2)
  length = min(length,len0)
  length = max(length,0.001*len0)
end if

ce1 = ck*1.9
ce2 = max(0.0, 0.93-ce1)
!if(k.eq.1) then
!  coefc = 3.9
!else
  coefc = ce1 + ce2*length/len0
!end if
epsi = coefc*tketmp**1.5/length

reyi = sqrt(20.0/3/visk)*(tketmp/sqrt(epsi))

if(inan.eq.1) then
  if(ieee_is_nan(epsi).or.ieee_is_nan(reyi).or..not.ieee_is_finite(epsi).or..not.ieee_is_finite(reyi)) then
    write(6,*) "in tice", epsi, reyi
    write(6,*) tke, n2, n, c
    write(6,*) len0, length, vel, lambda
    call wrf_error_fatal("nan detected in tice")
  end if
end if

eps1 = epsi
rey1 = reyi

if(iturb.ne.0) then

! CGS
epsi = max(epsi*10000,  10.0)
reyi = max(reyi,      5000.0)

! classification

if(epsi.lt.50) then
  if(reyi.lt.10000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,1)
      c2 = turbll_grid(i,j,2)
      c3 = turbll_grid(i,j,3)
      c4 = turbll_grid(i,j,4)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,5000.,10000.,10.,50.)
    end do
    end do
  else if(reyi.lt.20000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,2)
      c2 = turbll_grid(i,j,2)
      c3 = turbll_grid(i,j,4)
      c4 = turbll_grid(i,j,4)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,10000.,20000.,10.,50.)
    end do
    end do
  else
    do i=1, 13
    do j=1, 13
      array_ll(i,j)=extra_tice(turbll_grid(i,j,2),epsi,10.,reyi,20000.)
    end do
    end do
  end if
else if(epsi.lt.100) then
  if(reyi.lt.10000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,3)
      c2 = turbll_grid(i,j,4)
      c3 = turbll_grid(i,j,5)
      c4 = turbll_grid(i,j,6)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,5000.,10000.,50.,100.)
    end do
    end do
  else if(reyi.lt.20000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,4)
      c2 = turbll_grid(i,j,4)
      c3 = turbll_grid(i,j,6)
      c4 = turbll_grid(i,j,6)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,10000.,20000.,50.,100.)
    end do
    end do
  else
    do i=1, 13
    do j=1, 13
      array_ll(i,j)=extra_tice(turbll_grid(i,j,4),epsi,50.,reyi,20000.)
    end do
    end do
  end if
else if(epsi.lt.200) then
  if(reyi.lt.10000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,5)
      c2 = turbll_grid(i,j,6)
      c3 = turbll_grid(i,j,7)
      c4 = turbll_grid(i,j,7)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,5000.,10000.,100.,200.)
    end do
    end do
  else if(reyi.lt.20000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,6)
      c2 = turbll_grid(i,j,6)
      c3 = turbll_grid(i,j,7)
      c4 = turbll_grid(i,j,8)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,10000.,20000.,100.,200.)
    end do
    end do
  else
    do i=1, 13
    do j=1, 13
      array_ll(i,j)=extra_tice(turbll_grid(i,j,6),epsi,100.,reyi,20000.)
    end do
    end do
  end if
else if(epsi.lt.500) then
  if(reyi.lt.10000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,7)
      c2 = turbll_grid(i,j,7)
      c3 = turbll_grid(i,j,9)
      c4 = turbll_grid(i,j,9)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,5000.,10000.,200.,500.)
    end do
    end do
  else if(reyi.lt.20000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,7)
      c2 = turbll_grid(i,j,8)
      c3 = turbll_grid(i,j,9)
      c4 = turbll_grid(i,j,10)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,10000.,20000.,200.,500.)
    end do
    end do
  else ! rey.ge.20000
    do i=1, 13
    do j=1, 13
      array_ll(i,j)=extra_tice(turbll_grid(i,j,8),epsi,200.,reyi,20000.)
    end do
    end do
  end if
else if(epsi.lt.1000) then
  if(reyi.lt.10000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,9)
      c2 = turbll_grid(i,j,9)
      c3 = turbll_grid(i,j,11)
      c4 = turbll_grid(i,j,11)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,5000.,10000.,500.,1000.)
    end do
    end do
  else if(reyi.lt.20000) then
    do i=1, 13
    do j=1, 13
      c1 = turbll_grid(i,j,9)
      c2 = turbll_grid(i,j,9)
      c3 = turbll_grid(i,j,10)
      c4 = turbll_grid(i,j,11)
      array_ll(i,j)=bilinear(c1,c2,c3,c4,reyi,epsi,10000.,20000.,500.,1000.)
    end do
    end do
  else ! rey.ge.20000
    do i=1, 13
    do j=1, 13
      array_ll(i,j)=extra_tice(turbll_grid(i,j,10),epsi,500.,reyi,20000.)
    end do
    end do
  end if
else ! eps.ge.1000
  if(reyi.lt.10000) then
    do i=1, 13
    do j=1, 13
      array_ll(i,j)=extra_tice(turbll_grid(i,j,11),epsi,1000.,reyi,10000.)
    end do
    end do
  else
    do i=1, 13
    do j=1, 13
      array_ll(i,j)=extra_tice(turbll_grid(i,j,11),epsi,1000.,reyi,20000.)
    end do
    end do
  end if
end if

if(inan.eq.1) then
  do i=1, 13
  do j=1, 13
    if(.not.ieee_is_finite(array_ll(i,j))) then
      write(6,*) "In tice"
      write(6,*) epsi, reyi, tke, n2
      write(6,*) length, c, lambda, vel
      write(6,*) dx, dz, rho, visk
      call wrf_error_fatal("tice inf")
    end if
  end do
  end do
end if

do i=1, 13
do j=1, 13
  ckxx(WATER,WATER,i,j)=ckxx(WATER,WATER,i,j)*max(array_ll(i,j),1.0)
end do
end do

! water-graupel collision
if(epsi.lt.100.0) then
  do i=krming_gl, krmaxg_gl
  do j=krminl_gl, krmaxl_gl
    array_gl(i,j)=1.0
  end do
  end do
else if(epsi.lt.400.0) then
  do i=krming_gl, krmaxg_gl
  do j=krminl_gl, krmaxl_gl
    array_gl(i,j)=1.0                 *(400.0-epsi)/(400.0-100.0) &
                 +turbgl_grid_400(i,j)*(epsi-100.0)/(400.0-100.0)
  end do
  end do
else if(epsi.lt.800.0) then
  do i=krming_gl, krmaxg_gl
  do j=krminl_gl, krmaxl_gl
    array_gl(i,j)=turbgl_grid_400(i,j)*(800.0-epsi)/(800.0-400.0) &
                 +turbgl_grid_800(i,j)*(epsi-400.0)/(800.0-400.0)
  end do
  end do
else if(epsi.lt.1600.0) then
  do i=krming_gl, krmaxg_gl
  do j=krminl_gl, krmaxl_gl
    array_gl(i,j)=turbgl_grid_800(i,j)*            &
                 (1.0*(1600.0-epsi)/(1600.0-800.0) &
                 +1.7*( epsi-800.0)/(1600.0-800.0))
  end do
  end do
else
  do i=krming_gl, krmaxg_gl
  do j=krminl_gl, krmaxl_gl
    array_gl(i,j)=turbgl_grid_800(i,j)*1.7
  end do
  end do
end if

do i=krming_gl, krmaxg_gl
do j=krminl_gl, krmaxl_gl
  ckxx(GRAUPEL,WATER,i,j)=ckxx(GRAUPEL,WATER,i,j)*array_gl(i,j)
  ckxx(WATER,GRAUPEL,j,i)=ckxx(WATER,GRAUPEL,j,i)*array_gl(i,j)

  if(inan.eq.1) then
    if(.not.ieee_is_finite(array_gl(i,j))) then
      write(6,*) i, j, epsi
      write(6,*) turbgl_grid_400(i,j), turbgl_grid_800(i,j)
      call wrf_error_fatal("inf of array_gl")
    end if
  end if

end do
end do

end if

end subroutine tice

!====================================================================
real function bilinear(f1,f2,f3,f4,x,y,x1,x2,y1,y2)
implicit none

real, intent(in) :: f1,f2,f3,f4,x,y,x1,x2,y1,y2

real :: f, p,q

p = (x-x1)/(x2-x1)
q = (y-y1)/(y2-y1)
f = (1.0-p)*(1.0-q)*f1+ &
         p *(1.0-q)*f2+ &
    (1.0-p)*     q *f3+ &
         p *     q *f4

bilinear = f

end function bilinear

!====================================================================
subroutine collision(g, lwf, rf, t1, t2, t3, i, j, tt)
use ieee_arithmetic
implicit none

integer, intent(in) :: t1, t2, t3, i, j
real, intent(in) :: tt
real, dimension(nhydro,nkr), intent(inout) :: g, lwf, rf

integer :: k
real :: g1, m1, g2, m2, g31, lwf31, rf31, g32, lwf32, rf32
real :: crn, beta, r_flux, r_rem, x01, dg1, dg2, dg, x1
real :: m_liq1, m_liq2, m_rim1, m_rim2

k  = ima(i,j)

g1 = g(t1,i)
m1 = mass_hydro(i)

g2 = g(t2,j)
m2 = mass_hydro(j)

if(inan.eq.1) then
  if(ieee_is_nan(g1).or.ieee_is_nan(g2)) then
    write(6,*) i, j, t1, t2
    write(6,*) g1, lwf(t1,i), rf(t1,i), m1
    write(6,*) g2, lwf(t2,j), rf(t2,j), m2
    call wrf_error_fatal("NaN is detected before collision.")
  end if
  if(.not.ieee_is_finite(ckxx(t1,t2,i,j))) then
    write(6,*) i, j, t1, t2
    write(6,*) ywxx(t1,t2,1,i,j),ywxx(t1,t2,2,i,j),ywxx(t1,t2,3,i,j)
    call wrf_error_fatal("Befoer collision, Inf")
  end if
  if(ieee_is_nan(ckxx(t1,t2,i,j))) then
    call wrf_error_fatal("nan of ckxx")
  end if
end if

if(ckxx(t1,t2,i,j).lt.0.0) then
  write(6,*) i, j, t1, t2
  write(6,*) ywxx(t1,t2,1,i,j),ywxx(t1,t2,2,i,j),ywxx(t1,t2,3,i,j)
  call wrf_error_fatal("ckxx is negative.")
end if

if(      g(t1,i).lt.0.0 &
   .or.lwf(t1,i).lt.0.0 &
   .or.lwf(t1,i).gt.1.0 &
   .or. rf(t1,i).lt.0.0 &
   .or. rf(t1,i).gt.1.0 &
   .or.  g(t2,j).lt.0.0 &
   .or.lwf(t2,j).lt.0.0 &
   .or.lwf(t2,j).gt.1.0 &
   .or. rf(t2,j).lt.0.0 &
   .or. rf(t2,j).gt.1.0 ) then
  write(6,*) "strange values are inputted in the collision."
  write(6,*) t1, i, g(t1,i), lwf(t1,i), rf(t1,i)
  write(6,*) t2, j, g(t2,j), lwf(t2,j), rf(t2,j)
end if

x01 = ckxx(t1,t2,i,j)*g1*g2
x01 = min(x01,g1*m2)
if((t2.ne.t3).or.(j.ne.k)) x01 = min(x01,g2*m1)

dg1 = x01/m2
dg2 = x01/m1

g(t1,i) = g(t1,i)-dg1
g(t2,j) = g(t2,j)-dg2

g(t1,i) = max(g(t1,i),0.0)
if((t2.ne.t3).or.(j.ne.k)) g(t2,j) = max(g(t2,j),0.0)

! start of redistribution (exponential flux method) by Bott (JAS 2000)
dg    = dg1+dg2
g31   = g  (t3,k)
lwf31 = lwf(t3,k)
rf31  = rf (t3,k)

if(dg.gt.0.0) then

  if(k.lt.nkr) then
    crn   = cour(i,j)

    g32   = g  (t3,k+1)
    lwf32 = lwf(t3,k+1)
    rf32  = rf (t3,k+1)

    x1 = g32/(g31+dg)
    if(x1.gt.1.0e-5.and.abs(x1-1.0).gt.1.0e-5.and.x1.lt.1.0e5) then
      beta = log(x1)
      r_flux = (exp(beta*0.5)-exp(beta*(0.5-crn)))/(exp(beta*0.5)-exp(-beta*0.5))
!     r_flux = (exp(beta*0.5)-exp(beta*(0.5-crn)))/beta
    else
      r_flux = (crn+crn**2)/2.0
    end if

    if(r_flux.gt.1.0.or.r_flux.lt.0.0) then
      r_flux = max(0.0,min(1.0,r_flux))
    end if

    r_flux = min(r_flux,1.0+g31/dg)
    g(t3,k+1) = g32 + dg*r_flux

  else

    r_flux = 0.0

  end if

  r_rem   = 1.0 - r_flux
  g(t3,k) = g31 + dg*r_rem
  g(t3,k) = max(g(t3,k),0.0)

  if(inan.eq.1) then
    if(g(t1,i).lt.0.or.g(t2,j).lt.0..or.g(t3,k).lt.0.or.g(t3,k+1).lt.0) then
      write(6,*) t1, t2, t3
      write(6,*) i, j, k
      write(6,*) g1, g2, g31, g32
      write(6,*) dg1, dg2, dg
      write(6,*) g(t1,i), g(t2,j), g(t3,k), g(t3,k+1)
      write(6,*) m1, m2, x01, ckxx(t1,t2,i,j)
      write(6,*) r_flux, r_rem
      call wrf_error_fatal("negative in collision")
    end if
  end if

! rimed fraction
  if(tt.lt.t00.and.t3.eq.SNOW) then

    if(t1.eq.WATER) then
      m_rim1 = dg1
    else if(t1.eq.SNOW) then
      m_rim1 = dg1*rf(t1,i)
    else
      m_rim1 = 0.0
    end if

    if(t2.eq.WATER) then
      m_rim2 = dg2
    else if(t2.eq.SNOW) then
      m_rim2 = dg2*rf(t2,j)
    else
      m_rim2 = 0.0
    end if

    if(k.lt.nkr) then
      if(g(t3,k+1).gt.0.0) then
        rf(t3,k+1) = max(min((g32*rf32 + (m_rim1+m_rim2)*r_flux)/g(t3,k+1),1.0),0.0)
      else
        rf(t3,k+1) = 0.0
      end if
    end if

    if(g(t3,k).gt.0.0) then
      rf(t3,k) = max(min((g31*rf31 + (m_rim1+m_rim2)*r_rem)/g(t3,k),1.0),0.0)
    else
      rf(t3,k) = 0.0
    end if

! liquid water fraction
  else if(tt.ge.t00.and.t3.ge.SNOW.and.t3.le.HAIL) then

    m_liq1 = dg1*lwf(t1,i)
    m_liq2 = dg2*lwf(t2,j)
    
    if(k.lt.nkr) then
      if(g(t3,k+1).gt.0.0) then
        lwf(t3,k+1) = max(min((g32*lwf32 + (m_liq1+m_liq2)*r_flux)/g(t3,k+1),1.0),0.0)
      else
        lwf(t3,k+1) = 0.0
      end if
    end if
    if(g(t3,k).gt.0.0) then
      lwf(t3,k) = max(min((g31*lwf31 + (m_liq1+m_liq2)*r_rem)/g(t3,k),1.0),0.0)
    else
      lwf(t3,k) = 0.0
    end if

  end if

end if
! end of redistribution

if(inan.eq.1) then
  if(ieee_is_nan(g(t1,i)).or.ieee_is_nan(lwf(t1,i)).or.ieee_is_nan(rf(t1,i)).or. &
     ieee_is_nan(g(t2,j)).or.ieee_is_nan(lwf(t2,j)).or.ieee_is_nan(rf(t2,j)).or. &
     ieee_is_nan(g(t3,k)).or.ieee_is_nan(lwf(t3,k)).or.ieee_is_nan(rf(t3,k)).or. &
     ieee_is_nan(g(t3,k+1)).or.ieee_is_nan(lwf(t3,k+1)).or.ieee_is_nan(rf(t3,k+1))) then

    write(6,*) i, j, k
    write(6,*) t1, t2, t3
    write(6,*) g(t1,i), g(t2,j), g(t3,k), g(t3,k+1)
    write(6,*) lwf(t1,i), lwf(t2,j), lwf(t3,k), lwf(t3,k+1)
    write(6,*) rf(t1,i), rf(t2,j), rf(t3,k), rf(t3,k+1)
    write(6,*) g1, g2, ckxx(t1,t2,i,j)
    write(6,*) ywxx(t1,t2,1,i,j),ywxx(t1,t2,2,i,j),ywxx(t1,t2,3,i,j)
    write(6,*) dg1, dg2, dg
    write(6,*) x1, beta, crn, r_flux
    call wrf_error_fatal("NaN in collision")
  end if
end if

end subroutine collision

!==========================================================================
subroutine collision_poisson(g, lwf, rf, t1, t2, t3, i, j, tt)
use ieee_arithmetic
implicit none

integer, intent(in) :: t1, t2, t3, i, j
real, intent(in) :: tt
real, dimension(nhydro,nkr), intent(inout) :: g, lwf, rf

double precision, dimension(:), allocatable :: p_it, g_temp, m_temp

double precision, parameter :: ln_p_min = log(1.0d-5)
double precision, parameter :: norm_range = 4.6 ! equiv. to 1.0e-5
double precision, parameter :: p_crit = 100.0

double precision :: m1, m2, n1, n2, p, ln_p, ln_p_it, ln_p_it0, del_n2, crn
double precision :: dg, dg1, dg2, g31, g32, lwf31, lwf32, rf31, rf32
double precision :: x1, beta, r_flux, r_rem
double precision :: m_liq1, m_liq2, m_rim1, m_rim2
integer :: it, k, k1, it_min, it_max, it_max2
logical :: flag

! note: 2 is larger than 1 !
m1 = mass_hydro(i)
m2 = mass_hydro(j)
n1 = g(t1,i)/m1*dlnr
n2 = g(t2,j)/m2*dlnr

if(inan.eq.1) then
  if(ieee_is_nan(n1).or.ieee_is_nan(n2)) then
    write(6,*) i, j, t1, t2
    write(6,*) n1, lwf(t1,i), rf(t1,i), m1
    write(6,*) n2, lwf(t2,j), rf(t2,j), m2
    call wrf_error_fatal("NaN is detected before collision.")
  end if
  if(.not.ieee_is_finite(ckxx(t1,t2,i,j))) then
    write(6,*) i, j, t1, t2
    write(6,*) ywxx(t1,t2,1,i,j),ywxx(t1,t2,2,i,j),ywxx(t1,t2,3,i,j)
    call wrf_error_fatal("Befoer collision, Inf")
  end if
  if(ieee_is_nan(ckxx(t1,t2,i,j))) then
    call wrf_error_fatal("nan of ckxx")
  end if
end if

if(ckxx(t1,t2,i,j).lt.0.0) then
  write(6,*) i, j, t1, t2
  write(6,*) ywxx(t1,t2,1,i,j),ywxx(t1,t2,2,i,j),ywxx(t1,t2,3,i,j)
  call wrf_error_fatal("ckxx is negative.")
end if

if(      g(t1,i).lt.0.0 &
   .or.lwf(t1,i).lt.0.0 &
   .or.lwf(t1,i).gt.1.0 &
   .or. rf(t1,i).lt.0.0 &
   .or. rf(t1,i).gt.1.0 &
   .or.  g(t2,j).lt.0.0 &
   .or.lwf(t2,j).lt.0.0 &
   .or.lwf(t2,j).gt.1.0 &
   .or. rf(t2,j).lt.0.0 &
   .or. rf(t2,j).gt.1.0 ) then
  write(6,*) "strange values are inputted in the collision."
  write(6,*) t1, i, g(t1,i), lwf(t1,i), rf(t1,i)
  write(6,*) t2, j, g(t2,j), lwf(t2,j), rf(t2,j)
end if

p = ckxx(t1,t2,i,j)/dlnr*n1
p = min(p, n1/n2)

if(p.gt.p_crit) then

  it_min = p - nint(norm_range*sqrt(p_crit))
  it_max = p + nint(norm_range*sqrt(p_crit))
  allocate(p_it(it_min:it_max))

  do it = it_min, it_max
    p_it(it) = 1.0d0/sqrt(2.0d0*pi*p_crit)*exp(-((it-p)**2/(2.0*p_crit)))
  end do

else

  it = 0
  flag = .false.
  ln_p = log(p)

  do
    if(it.eq.0) then
      ln_p_it = -p
    else
      ln_p_it = ln_p_it + ln_p - log(dble(it))
    end if

    if((flag.eqv..false.).and.(ln_p_it.gt.ln_p_min)) then
      flag = .true.
      it_min = it
      ln_p_it0 = ln_p_it
    end if

    if((flag.eqv..true.).and.(ln_p_it.lt.ln_p_min)) then
      it_max = it-1
      exit
    end if

    it = it + 1
  end do

  allocate(p_it(it_min:it_max))

  p_it(it_min) = exp(ln_p_it0)
  ln_p_it = ln_p_it0

  if(it_max.ge.it_min+1) then
    do it = it_min+1, it_max
      ln_p_it = ln_p_it + ln_p - log(dble(it))
      p_it(it) = exp(ln_p_it)
    end do
  end if

end if

it_min = max(it_min,1)

allocate(m_temp(it_min:it_max))
allocate(g_temp(it_min:it_max))

it_max2 = it_max

do it = it_min, it_max

  del_n2  = n2*p_it(it)
  n1 = n1 - it*del_n2

  if(n1.lt.0.0) then
    n1 = n1 + it*del_n2 ! restore
    del_n2 = n1/it
    n1 = 0.0
    p_it(it) = del_n2/n2
    it_max2 = it
    exit
  end if

  m_temp(it) = m2 + it*m1
  g_temp(it) = m_temp(it)*del_n2/dlnr

end do

n2 = n2*(1.0 - sum(p_it(it_min:it_max2)))

if(n1.lt.0.0.or.n2.lt.0.0) then
  write(*,*) n1, n2, "negative"
end if

if(ieee_is_nan(n1).or.ieee_is_nan(n2)) then
  write(*,*) i, j, m1, m2
  write(*,*) p, n1, n2, del_n2
  write(*,*) "NaN detected in poisson"
  stop
end if

g(t1,i) = n1*m1/dlnr
g(t2,j) = n2*m2/dlnr

! start of redistribution (exponential flux method) by Bott (JAS 2000)
do it = it_min, it_max2

  k = nkr
  do k1=1, nkr-1
    if(m_temp(it).ge.mass_hydro(k1).and.m_temp(it).lt.mass_hydro(k1+1)) then
      k = k1
      exit
    end if
  end do

  dg    = g_temp(it)
  dg1   = (it*m1/(it*m1+m2))*dg
  dg2   = (   m2/(it*m1+m2))*dg

  g31   = g  (t3,k)
  lwf31 = lwf(t3,k)
  rf31  = rf (t3,k)

  if(k.lt.nkr) then
    crn = log(m_temp(it)/mass_hydro(k))/(3.0*dlnr)

    g32   = g  (t3,k+1)
    lwf32 = lwf(t3,k+1)
    rf32  = rf (t3,k+1)

    x1  = g32/(g31+dg)
    if(x1.gt.1.0e-5.and.x1.lt.1.0e5.and.abs(x1-1.0).gt.1.0e-5) then
      beta = log(x1)
      r_flux = (exp(beta*0.5)-exp(beta*(0.5-crn)))/(exp(beta*0.5)-exp(-beta*0.5))
!     r_flux = (exp(beta*0.5)-exp(beta*(0.5-crn)))/beta
    else
      r_flux = (crn+crn**2)/2.0
    end if

    if(r_flux.gt.1.0.or.r_flux.lt.0.0) then
      r_flux = max(0.0d0,min(1.0d0,r_flux))
    end if
    r_flux = min(r_flux,1.0+g31/dg)

    g(t3,k+1) = g32 + dg*r_flux
  else
    r_flux = 0.0
  end if

  r_rem = 1.0-r_flux
  g(t3,k) = max(g31 + dg*r_rem,0.0)

  ! rimed fraction
  if(tt.lt.t00.and.t3.eq.SNOW) then

    if(t1.eq.WATER) then
      m_rim1 = dg1
    else if(t1.eq.SNOW) then
      m_rim1 = dg1*rf(t1,i)
    else
      m_rim1 = 0.0
    end if

    if(t2.eq.WATER) then
      m_rim2 = dg2
    else if(t2.eq.SNOW) then
      m_rim2 = dg2*rf(t2,j)
    else
      m_rim2 = 0.0
    end if

    if(k.lt.nkr) then
      if(g(t3,k+1).gt.0.0) then
        rf(t3,k+1) = (g32*rf32 + (m_rim1+m_rim2)*r_flux)/g(t3,k+1)
      else
        rf(t3,k+1) = 0.0
      end if
    end if

    if(g(t3,k).gt.0.0) then
      rf(t3,k) = (g31*rf31 + (m_rim1+m_rim2)*r_rem)/g(t3,k)
    else
      rf(t3,k) = 0.0
    end if

  ! liquid water fraction
  else if(tt.ge.t00.and.t3.ge.SNOW.and.t3.le.HAIL) then

    m_liq1 = dg1*lwf(t1,i)
    m_liq2 = dg2*lwf(t2,j)

    if(k.lt.nkr) then
      if(g(t3,k+1).gt.0.0) then
        lwf(t3,k+1) = (g32*lwf32 + (m_liq1+m_liq2)*r_flux)/g(t3,k+1)
      else
        lwf(t3,k+1) = 0.0
      end if
    end if
    if(g(t3,k).gt.0.0) then
      lwf(t3,k) = (g31*lwf31 + (m_liq1+m_liq2)*r_rem)/g(t3,k)
    else
      lwf(t3,k) = 0.0
    end if

  end if

end do

end subroutine collision_poisson

!==========================================================================
subroutine breakup_water(g1, dtime_coll)
implicit none

real, dimension(nkr), intent(inout) :: g1
real, intent(in) :: dtime_coll

real, dimension(nkr) :: g_init
real :: dropmass_before, dropmass_after, dropmass_delta, dtime_break
integer :: i, ndiv

dropmass_before = sum(g1)
g_init = g1
ndiv = 1

if(dropmass_before.eq.0.0) return

outer: do
  g1 = g_init
  dtime_break = dtime_coll/ndiv

  inner: do i=1, ndiv
    call coll_breakup(g1,dtime_break)

!! for implicit method
    if(ibreakup.eq.1) then
      dropmass_after = sum(g1)
      dropmass_delta = abs(dropmass_after-dropmass_before)
      if(dropmass_delta.gt.tol*dropmass_before) then
        ndiv = ndiv*2
        cycle outer
      end if
!! for explicit method
    else if(ibreakup.eq.2) then
      if(minval(g1).lt.0.0) then
        ndiv = ndiv*2
        cycle outer
      end if
    else
      call wrf_error_fatal("something wrong in ibreakup")
    end if

  end do inner

  exit outer

end do outer

end subroutine breakup_water

!==========================================================================
subroutine coll_breakup(g1,dt)
use ieee_arithmetic
implicit none

real, dimension(nkr), intent(inout) :: g1
real, intent(in) :: dt

integer, parameter :: offset = nkr-krbreak ! start of breakup

real, dimension(nkr) :: ft, gain, loss, weight, beta
integer :: i, j, k

gain = 0.0
loss = 0.0

if(inan.eq.1) then
  do i=1, nkr
    if(ieee_is_nan(g1(i))) then
      write(6,*) "in breakup", i
      call wrf_error_fatal("nan in breakup1")
    end if
  end do
end if

do i=1, nkr
  if(g1(i).lt.0.0) then
    write(6,*) "in breakup", i, g1(i)
    call wrf_error_fatal("negative in breakup")
  end if
end do

do k=1, nkr
  if(k.lt.nkr) then
    weight(k)=2.0/(mass_hydro(k+1)**2 - mass_hydro(k)**2)
  else
    weight(k)=2.0/((mass_hydro(k)*2)**2 - mass_hydro(k)**2)
  end if
end do

ft = g1/(3.0*mass_hydro**2)

do k=offset+1, nkr
  do i=offset+1, nkr
  do j=offset+1, i
    gain(k)=gain(k)+ft(i)*ft(j)*pkij(k-offset,i-offset,j-offset)
  end do
  end do

  do j=offset+1, nkr
    loss(k)=loss(k)+ft(j)*qkj(k-offset,j-offset)
  end do
end do

do k=offset+1, nkr
! analytic solution by Hyunho
!    if(loss(k).gt.0.0) then
!      beta(k) = exp(-weight(k)*loss(k)*dt)
!      ft(k) = gain(k)/loss(k)*(1.0-beta(k))+ft(k)*beta(k)
!    else
!      ft(k) = ft(k) + weight(k)*gain(k)*dt
!    end if
! implicit solution by Hyunho -> positive, stable
    if(ibreakup.eq.1) then
      ft(k) = (ft(k) + weight(k)*gain(k)*dt)/(1.0 + weight(k)*loss(k)*dt)
! explicit solution in the original HUCM
    else if(ibreakup.eq.2) then
      ft(k) = ft(k) + weight(k)*(gain(k) - ft(k)*loss(k))*dt
    else
      call wrf_error_fatal("something wrong in ibreakup")
    end if
  end do

g1 = ft*(3.0*mass_hydro**2)

if(inan.eq.1) then
  do i=1, nkr
    if(ieee_is_nan(g1(i))) then
      write(6,*) "in breakup", i
      write(6,*) gain(i), loss(i), weight(i)
      call wrf_error_fatal("nan in breakup2")
    end if
  end do
end if

end subroutine coll_breakup

!===================================================================
real function ecoalmass(x1,x2)
implicit none
! input is mks, but process is cgs

real, intent(in) :: x1, x2
real :: deta, dksi

! in cgs
deta = (x1/rho_water/pi*6.0)**(1.0/3.0)*1.0e2
dksi = (x2/rho_water/pi*6.0)**(1.0/3.0)*1.0e2

ecoalmass = ecoaldiam(deta,dksi)

end function ecoalmass

!====================================================================
real function ecoaldiam(deta,dksi)
implicit none

real, intent(in) :: deta, dksi

real, parameter :: qmin = 250.0e-4 ! 250 micro m
real, parameter :: qmax = 500.0e-4 ! 500 micro m

real :: dgr, dkl, q, e, x, sin1, cos1, e1, e2

dgr=max(deta,dksi)
dkl=min(deta,dksi)
q=(0.5*(dkl+dgr))*0.5 ! mean radius

if(dkl.lt.100.0e-4) then
  e = 1.0
else if(q.lt.qmin) then
  e = ecoalochs(dgr,dkl)
else if(q.ge.qmin.and.q.lt.qmax) then
  x=(q-qmin)/(qmax-qmin)
  sin1=sin(pi/2.0*x)
  cos1=cos(pi/2.0*x)
  e1=ecoalochs(dgr,dkl)
  e2=ecoallowlist(dgr,dkl)
  e=cos1**2*e1+sin1**2*e2
else
  e=ecoallowlist(dgr,dkl)
endif

ecoaldiam=max(min(1.0,e),1.0e-6)

end function ecoaldiam

!====================================================================
!coalescence efficiency (low & list)
real function ecoallowlist(dgr,dkl)
implicit none

real, intent(in) :: dgr, dkl
real, parameter :: sigma = 72.8   ! surface tension,[sigma]=g/s^2 (7.28e-2 n/m)
real, parameter :: aka = 0.778    ! empirical constant
real, parameter :: akb = 2.61e-4  ! empirical constant,[b]=2.61e6 m^2/j^2

real :: cke, st, sc, dstsc, et, qq1, qq2, ecl

call collenergy(dgr,dkl,cke,st,sc)

dstsc=st-sc         ! diff. of surf. energies   [dstsc] = g*cm^2/s^2
et=cke+dstsc        ! coal. energy,             [et]    =     "

if(et.lt.50.0) then    ! et < 5 uj (= 50 g*cm^2/s^2)
  qq1=aka/(1.0+(dkl/dgr))**2
  qq2=akb*sigma*(et**2)/sc
  ecl=qq1*exp(-qq2)
else
  ecl=0.0
endif

ecoallowlist=ecl

end function ecoallowlist

!====================================================================
!coalescence efficiency (beard and ochs)
real function ecoalochs(d_l,d_s)
implicit none

real, intent(in) :: d_l, d_s

real, parameter :: sigma=72.8  ! surface tension [sigma] = g/s^2 (7.28e-2 n/m)

real :: r_s, r_l, p, vtl, vts, dv, weber_number, p1, p2, p3, g, x

r_s=0.5*d_s
r_l=0.5*d_l
p=r_s/r_l

vtl=vtbeard(d_l)
vts=vtbeard(d_s)
dv=max(abs(vtl-vts),1.0e-30)

weber_number=r_s*dv**2/sigma

p1=1.0+p
p2=1.0+p**2
p3=1.0+p**3

x=sqrt(weber_number)*(2**1.5)/(6.0*pi)*p**4*p1/(p2*p3)
ecoalochs=0.767-10.14*x

end function ecoalochs

!====================================================================
!calculating the collision energy
subroutine collenergy(dgr,dkl,cke,st,sc)
implicit none

real, intent(in) :: dgr, dkl
real, intent(out) :: cke, st, sc

real, parameter :: rho=rho_water/1000.  ! water density,[rho]=g/cm^3
real, parameter :: sigma=72.8  ! surf. tension,(h2o,20°c)=7.28d-2 n/m [sigma]=g/s^2

real :: dgka2, dgka3, dgkb3, v1, v2, dv

dgka2=(dgr**2)+(dkl**2)
dgka3=(dgr**3)+(dkl**3)

if(dgr.ne.dkl) then
  v1=vtbeard(dgr)
  v2=vtbeard(dkl)
  dv=max((v1-v2)**2,1.0e-30)
  dgkb3=(dgr**3)*(dkl**3)
  cke=rho*pi/10.0*dv*dgkb3/dgka3  ! collision energy [cke]=g*cm^2/s^2
else
  cke = 0.0
endif

st=pi*sigma*dgka2             ! surf.energy (parent drop)  
sc=pi*sigma*dgka3**(2.0/3.0)  ! surf.energy (coal.system) 

!w1=cke/(sc+epsf)                  ! weber number 1
!w2=cke/(st+epsf)                  ! weber number 2
!dc=dgka3**(1.0d0/3.0d0)           ! diam. of coal. system

end subroutine collenergy

!====================================================================
!calculating terminal velocity (beard-formula)
real function vtbeard(diam)
implicit none

real, intent(in) :: diam

real :: aa
integer :: kr

aa = diam/2.0           ! radius in cm

if(aa.eq.r_hydro(WATER,1)) vtbeard=vt_hydro(WATER,1)

do kr=1,nkr-1
  if(aa.gt.r_hydro(WATER,kr).and.aa.le.r_hydro(WATER,kr+1)) then
    vtbeard=vt_hydro(WATER,kr+1)
  endif
enddo   

end function vtbeard

!====================================================================
subroutine diffuse(tt, pp, qq, rhoa, dt, ff, m_new, dT_ce, dT_ds)
use ieee_arithmetic
implicit none

real, dimension(nhydro,nkr), intent(in ) :: ff
real, dimension(nhydro,nkr), intent(out) :: m_new
real, intent(in) :: pp, rhoa, dt
real, intent(inout) :: tt, qq
real, intent(out) :: dT_ce, dT_ds

integer :: m, n
real, dimension(nhydro,nkr) :: dmdt
real :: ffi, capi, vti, dmdti, m23, m13, delm, a, Sw, Si, Sw_dt, Si_dt, &
        rw, ri, pw, pi, tto, qqo, swo, sio, sum_dmdt_w, sum_dmdt_i

! supersaturation
Sw = pp*qq/(eps+qq)/esw(tt)-1.0
Si = pp*qq/(eps+qq)/esi(tt)-1.0

tto = tt
qqo = qq
swo = Sw
sio = Si

do m=1, nhydro
do n=1, nkr
  capi =  r_hydro(m,n)
   vti = vt_hydro(m,n)

  call rate(pp, tt, rhoa, capi, vti, dmdti)

  dmdt(m,n) = dmdti
end do
end do

sum_dmdt_w = dot_product(dmdt(WATER,:),ff(WATER,:))
sum_dmdt_i = 0.0
do m=COLUMN, HAIL
  sum_dmdt_i = sum_dmdt_i + dot_product(dmdt(m,:),ff(m,:))
end do

if(qq.gt.1.0e-10) then
  a = eps/(qq*(eps+qq))
  rw = (a+lv*lv/Cp/Rv/(tt**2))*(1.0+Sw)*sum_dmdt_w
  ri = (a+lv*ls/Cp/Rv/(tt**2))*(1.0+Sw)*sum_dmdt_i
  pw = (a+ls*lv/Cp/Rv/(tt**2))*(1.0+Si)*sum_dmdt_w
  pi = (a+ls*ls/Cp/Rv/(tt**2))*(1.0+Si)*sum_dmdt_i
else
  rw = pp/(esw(tt)*eps)*sum_dmdt_w
  ri = pp/(esw(tt)*eps)*sum_dmdt_i
  pw = pp/(esi(tt)*eps)*sum_dmdt_w
  pi = pp/(esi(tt)*eps)*sum_dmdt_i
end if

call supsat(Sw, Si, Sw_dt, Si_dt, rw, ri, pw, pi, dt)

if(ieee_is_nan(Sw_dt).or.ieee_is_nan(Si_dt)) then
  write(6,*) Sw, Si, Swo, Sio
  write(6,*) rw, ri, pw, pi
  write(6,*) a, sum_dmdt_w, sum_dmdt_i
  call wrf_error_fatal("nan from supsat")
end if

do m=1, nhydro
do n=1, nkr
  if(m.eq.WATER) then
    m23 = mass_hydro(n)**(2.0/3.0)
    m13 = mass_hydro(n)**(1.0/3.0)
    delm = (2.0/3.0)*dmdt(m,n)*Sw_dt/m13
    a = m23+delm
    if(a.ge.0.0) then
      m_new(m,n)=a**(3.0/2.0)
    else
      m_new(m,n)=0.0
    endif
  else
    m_new(m,n)=mass_hydro(n)+dmdt(m,n)*Si_dt
  end if
end do
end do

!call S_to_tq(Sw, Si, tt, pp, qq)

delm = dot_product(m_new(WATER,:)-mass_hydro,ff(WATER,:))
dT_ce = delm*Lv/Cp
tt = tt + dT_ce
qq = qq - delm

delm = 0.0
do m=COLUMN, nhydro
  delm = delm + dot_product(m_new(m,:)-mass_hydro,ff(m,:))
end do
dT_ds = delm*Ls/Cp
tt = tt + dT_ds
qq = qq - delm

!if(abs(tt-tto).gt.5.0) then
!  write(6,*) "in diffusion"
!  write(*,*) tto, qqo, tt, qq
!  write(*,*) rhoa, pp
!  write(*,*) pw, rw, sum_dmdt_w, sum_dmdt_i
!  write(*,*) Sw, Si, swo, sio
!  write(*,*) Sw_dt, Si_dt
!  call wrf_error_fatal("too high temperature in diffusion")
!end if

end subroutine diffuse

!====================================================================
subroutine rate(pp, tt, rhoa, capi, vti, dmdt)
implicit none

real, intent(in) :: pp, tt, rhoa, capi, vti
real, intent(out) :: dmdt

real:: df_my, tc_my, vi_my, Nsh, Nre, reshm, ventplm, fd1, fk1

df_my=difv*(p00/pp)*(tt/t00)**1.94
tc_my=thca+tt*7.11e-5
vi_my=(visa*(291.15+120.0)/(tt+120.0)*(tt/291.15)**1.5)/rhoa

! Shumidt,Reynolds numbers and ventilation coeff. (PK p.541)
Nsh=vi_my/df_my
Nre=2.0*capi*vti/vi_my
reshm=Nsh**(1.0/3.0)*Nre**0.5
if(reshm.lt.1.4) then
  ventplm=1.0+0.108*reshm**2
else
  ventplm=0.78+0.308*reshm
endif
ventplm=min(ventplm,ventpl_max)

! Fd and Fk (water and ice have different values)
fd1=Rv*tt/df_my/esw(tt)
fk1=(Lv/Rv/tt-1.0)*Lv/tc_my/tt

dmdt=4.0*pi*capi*ventplm/(fk1+fd1)

end subroutine rate

!====================================================================
subroutine supsat(del1, del2, del1int, del2int, rw, ri, pw, pi, dt)
implicit none

real, intent(in   ) :: rw, ri, pw, pi, dt
real, intent(inout) :: del1, del2
real, intent(out  ) :: del1int, del2int
real, parameter :: epsilon = 1.0e-10
real, parameter :: ratio_threshold = 1.0e2

real :: expr, expp, expb, expg, alpha, beta, gamma, c11, c12, c21, c22, del10, del20

del10 = del1
del20 = del2

if(rw.le.epsilon.and.pi.le.epsilon) then

!  del1n = del1+dyn1*dt
!  del2n = del2+dyn2*dt
!  del1int = del1*dt+dyn1*dt*dt/2.0
!  del2int = del2*dt+dyn2*dt*dt/2.0

  del1 = del10
  del2 = del20
  del1int = del10*dt
  del2int = del20*dt

else if(rw.gt.ratio_threshold*pi) then

  expr=exp(-rw*dt)

!  del1n   =  del1*expr + (dyn1/rw)*(1.0-expr)
!  del1int = -del1*expr/rw + dyn1*dt/rw + dyn1*expr/(rw*rw) + del1/rw - dyn1/(rw*rw)
!  del2n   = pw*del1*expr/rw - pw*dyn1*dt/rw - pw*dyn1*expr/(rw*rw) + dyn2*dt + del2 - pw*del1/rw + pw*dyn1/(rw*rw)
!  del2int = pw*del1*expr/(-rw*rw) - pw*dyn1*dt*dt/(2.0*rw) + pw*dyn1*expr/(rw*rw*rw) + dyn2*dt*dt/2.0 + del2*dt - pw*del1*dt/rw + pw*dyn1*dt/(rw*rw) + pw*del1/(rw*rw) - pw*dyn1/(rw*rw*rw)

  del1    = del10*expr
  del1int = del10/rw*(1.0-expr)
  del2    = del20 - del10*(pw/rw)*(1.0-expr)
  del2int = del20*dt + pw*del10/(rw*rw)*(1.0-expr) - pw*del10*dt/rw

else if(pi.gt.ratio_threshold*rw) then

  expp=exp(-pi*dt)

!  del1n   =  ri*del2*expp/pi - ri*dyn2*dt/pi - ri*dyn2*expp/(pi*pi) + dyn1*dt + del1 - ri*del2/pi + ri*dyn2/(pi*pi)
!  del1int = -ri*del2*expp/(pi*pi) - ri*dyn2*dt*dt/(2.0*pi) + ri*dyn2*expp/(pi*pi*pi) + dyn1*dt*dt/2.0 + del1*dt - ri*del2*dt/pi + ri*dyn2*dt/(pi*pi) + ri*del2/(pi*pi) - ri*dyn2/(pi*pi*pi)
!  del2n   =  del2*expp + (dyn2/pi)*(1.0-expp)
!  del2int = -del2*expp/pi + dyn2*dt/pi + dyn2*expp/(pi*pi) + del2/pi - dyn2/(pi*pi)

  del2    = del20*expp
  del2int = del20/pi*(1.0-expp)
  del1    = del10 - del20*(ri/pi)*(1.0-expp)
  del1int = del10*dt + ri*del20/(pi*pi)*(1.0-expp) - ri*del20*dt/pi

else

  alpha = sqrt((rw-pi)**2+4.0*ri*pw)
  beta  = 0.5*(alpha+rw+pi)
  gamma = 0.5*(alpha-rw-pi)
  expb  = exp(-beta*dt)
  expg  = exp(gamma*dt)
!  g31 =  pi*dyn1 - ri*dyn2
!  g32 = -pw*dyn1 + rw*dyn2
!  g2  =  rw*pi-ri*pw

!  c11=(beta*del1 - rw*del1 - ri*del2 - beta*g31/g2 + dyn1)/alpha
!  c21=(gama*del1 + rw*del1 + ri*del2 - gama*g31/g2 - dyn1)/alpha
!  c12=(beta*del2 - pw*del1 - pi*del2 - beta*g32/g2 + dyn2)/alpha
!  c22=(gama*del2 + pw*del1 + pi*del2 - gama*g32/g2 - dyn2)/alpha

  c11=(beta *del10 - rw*del10 - ri*del20)/alpha
  c21=(gamma*del10 + rw*del10 + ri*del20)/alpha
  c12=(beta *del20 - pw*del10 - pi*del20)/alpha
  c22=(gamma*del20 + pw*del10 + pi*del20)/alpha

!  del1n = c11*expg + c21*expb + g31/g2
!  del2n = c12*expg + c22*expb + g32/g2
!  del1int = c11*expg/gama - c21*expb/beta + c21/beta - c11/gama + g31*dt/g2
!  del2int = c12*expg/gama - c22*expb/beta + c22/beta - c12/gama + g32*dt/g2

  del1 = c11*expg + c21*expb
  del2 = c12*expg + c22*expb
  if(abs(gamma).gt.1.0e-5) then
    del1int = c21/beta*(1.0-expb) - c11/gamma*(1.0-expg)
    del2int = c22/beta*(1.0-expb) - c12/gamma*(1.0-expg)
  else
    del1int = c21/beta*(1.0-expb) + c11*dt
    del2int = c22/beta*(1.0-expb) + c12*dt
  end if

end if

del1int = max(del1int,-dt)
del2int = max(del2int,-dt)

end subroutine supsat

!=============================================================================
subroutine S_to_tq(Sw, Si, tt, pp, qq)
implicit none

real, intent(in) :: Sw, Si, pp
real, intent(out) :: tt, qq

real :: es, e

tt = (b_w-b_i)/((a_w-a_i)+log((1.0+Sw)/(1.0+Si)))

es = esw(tt)
e  = es*(1.0+Sw)
qq = eps*e/(pp-e)

end subroutine S_to_tq

!====================================================================
subroutine remapping(rn,fi)
use ieee_arithmetic
implicit none

real, dimension(nkr), intent(in) :: rn
real, dimension(nkr), intent(inout) :: fi

integer :: i, kr
integer :: isign_diffusional_growth

real, dimension(nkr) :: psi,rrs
real :: rnt,rrt,rrp,rrm,gn1,gn2,gn3

psi = 0.0

! rrs(i), default masses before condensation or evaporation
! rn(i),  new     masses after  condensation or evaporation

isign_diffusional_growth = 0
rrs = mass_hydro

if(inan.eq.1) then

do kr=1, nkr
  if(ieee_is_nan(rn(kr))) then
    write(6,*) kr
    call wrf_error_fatal("nan in rn")
  end if
  if(ieee_is_nan(fi(kr))) then
    write(6,*) kr
    call wrf_error_fatal("nan in fi1")
  end if
end do

end if


do kr=1,nkr
  if(rn(kr).ne.rrs(kr)) then
    isign_diffusional_growth=1
    exit
  endif
enddo

if(isign_diffusional_growth.ne.0) then

! kovetz-olund method
  if(iremap.eq.1) then
    do kr=1,nkr
      if(rn(kr).gt.0.0.and.fi(kr).ne.0.0) then
        if(rn(kr).lt.rrs(1)) then
          psi(1  )=psi(1  )+fi(kr)*(rn(kr)/rrs(1  ))
        else if(rn(kr).ge.rrs(nkr)) then
          psi(nkr)=psi(nkr)+fi(kr)*(rn(kr)/rrs(nkr))
        else
          do i=1, nkr-1
            if(rn(kr).ge.rrs(i).and.rn(kr).lt.rrs(i+1)) exit
          enddo
          rnt=rn(kr  )
          rrt=rrs(i  )
          rrp=rrs(i+1)
          gn1=(rrp-rnt)/(rrp-rrt)
          gn2=(rnt-rrt)/(rrp-rrt)
          psi(i)  =psi(i)  +gn1*fi(kr)
          psi(i+1)=psi(i+1)+gn2*fi(kr)
        endif
      endif
    enddo

! 3point method
  else if(iremap.eq.2) then
    do kr=1,nkr
      if(rn(kr).gt.0.0.and.fi(kr).ne.0.0) then
        if(rn(kr).lt.rrs(1)) then
          psi(1  )=psi(1  )+fi(kr)*(rn(kr)/rrs(1  ))
        else if(rn(kr).ge.rrs(nkr)) then
          psi(nkr)=psi(nkr)+fi(kr)*(rn(kr)/rrs(nkr))
        else
          do i=1, nkr-1
            if((rn(kr).ge.rrs(i).and.rn(kr).lt.rrs(i+1))) exit
          enddo
          if(i.eq.1) then ! use KO method
            rnt=rn(kr  )
            rrt=rrs(i  )
            rrp=rrs(i+1)
            gn1=(rrp-rnt)/(rrp-rrt)
            gn2=(rnt-rrt)/(rrp-rrt)
            psi(i)  =psi(i)  +gn1*fi(kr)
            psi(i+1)=psi(i+1)+gn2*fi(kr)
          else
            ! Ref: Khain et al.(JAS 2008)
            ! rrm < rrt < rnt < rrp
            ! Note: gn1 is negative.
            rnt=rn(kr  )
            rrm=rrs(i-1)
            rrt=rrs(i  )
            rrp=rrs(i+1)

            gn1=(rrp-rnt)*(rrt-rnt)/(rrp-rrm)/(rrt-rrm)  ! + - + + = -
            gn2=(rrp-rnt)*(rrm-rnt)/(rrp-rrt)/(rrm-rrt)  ! + - + - = +
            gn3=(rrm-rnt)*(rrt-rnt)/(rrm-rrp)/(rrt-rrp)  ! - - - - = +

            if(psi(i-1).lt.-gn1*fi(kr)) then ! use KO method
              gn1=(rrp-rnt)/(rrp-rrt)
              gn2=(rnt-rrt)/(rrp-rrt)
              psi(i)  =psi(i)  +gn1*fi(kr)
              psi(i+1)=psi(i+1)+gn2*fi(kr)
            else ! use 3point method. Is it physical?
              psi(i-1)=psi(i-1)+gn1*fi(kr)
              psi(i)  =psi(i  )+gn2*fi(kr)
              psi(i+1)=psi(i+1)+gn3*fi(kr)
            endif
          endif
        endif
      endif
    enddo

  endif

else
  do kr=1,nkr
    psi(kr)=fi(kr)
  enddo

endif

! final output
do kr=1, nkr
  fi(kr)=psi(kr)
end do

if(inan.eq.1) then
  do kr=1, nkr
    if(ieee_is_nan(fi(kr))) then
      write(6,*) kr
      call wrf_error_fatal("nan in fi2")
    end if
  end do
end if

end subroutine remapping

!====================================================================!
subroutine ice_nucl(tt,pp,qq,rhoa,dsupi,dtemp,ffi1,ffi2,ffi3)
use ieee_arithmetic
implicit none

real, intent(in) :: pp, rhoa, dsupi, dtemp
real, intent(inout) :: tt, qq, ffi1, ffi2, ffi3

real, parameter :: tt_nucl1 =  -5.0
real, parameter :: tt_nucl2 =  -2.0
real, parameter :: tt_nucl3 = -20.0

real :: tpc, supi, n1, n2, dn1, dn2, dn, dmr

if(inan.eq.1) then
  if(ieee_is_nan(pp)) then
    write(6,*) "before ice nucl"
    write(6,*) tt, qq, pp, rhoa
    write(6,*) dsupi, dtemp
  end if
end if

supi = (pp*qq/(eps+qq))/esi(tt)-1.0
tpc = tt-t00
if(supi.le.0.0.or.tpc.gt.0.0.or.tpc.lt.tt_nucl_ice_min) return

if(dsupi.gt.0.0.and.tpc.lt.tt_nucl1) then
  n1=c1_mey*exp(a1_mey+b1_mey*min(100.0*supi,supi_max))
  dn1=b1_mey*n1*(min(100.0*dsupi,dsupi_max))
else
  dn1=0.0
endif

if(dtemp.lt.0.0.and.tpc.lt.tt_nucl2) then
  n2=c2_mey*exp(a2_mey-b2_mey*max(tpc,tt_nucl3))
  dn2=-b2_mey*n2*dtemp
else
  dn2=0.0
endif

dmr = min((dn1+dn2)*mass_hydro(1)/rhoa,qq)
dn  = dmr/mass_hydro(1)*rhoa

if     (tpc.gt. -4.0) then
  ffi2 = ffi2 + dn
else if(tpc.gt. -8.1) then
  ffi1 = ffi1 + dn
else if(tpc.gt.-12.7) then
  ffi2 = ffi2 + dn
else if(tpc.gt.-17.8) then
  ffi3 = ffi3 + dn
else if(tpc.gt.-22.4) then
  ffi2 = ffi2 + dn
else
  ffi1 = ffi1 + dn
endif

if(ffi1.lt.0.0.or.ffi2.lt.0.0.or.ffi3.lt.0.0) then
  write(6,*) tt, qq, dn1, dn2
  write(6,*) ffi1, ffi2, ffi3
  write(6,*) supi, dsupi, dtemp
  call wrf_error_fatal("negative in ice nucleation")
end if

if(inan.eq.1) then
  if(ieee_is_nan(ffi1).or.ieee_is_nan(ffi2).or.ieee_is_nan(ffi3)) then
    write(6,*) tt, qq, pp, rhoa
    write(6,*) supi, dsupi, dtemp
    write(6,*) n1, n2, dn1, dn2
    write(6,*) dn1, dn2
    write(6,*) ffi1, ffi2, ffi3
    call wrf_error_fatal("nan in ice nucleation")
  end if
end if

qq = qq - dmr
tt = tt + dmr*Ls/Cp

end subroutine ice_nucl

!===================================================================!
subroutine water_nucl(tt,pp,qq,rhoa,fccni,rccni,ffwi)
implicit none

real, intent(in) :: pp, rhoa, rccni
real, intent(inout) :: tt, qq, fccni, ffwi

real :: supw, akoe1, bkoe1, rcrit, rccni_p, dmr, qq1, tt1, factor

supw = (pp*qq/(eps+qq))/esw(tt)-1.0
if(supw.le.0.0.or.tt.lt.t00+tt_nucl_drop_min) return

akoe1=akoe/tt
bkoe1=bkoe
rcrit=(akoe1/3.0)*(4.0/bkoe1/(supw**2))**(1.0/3.0)*1.0e-2 ! cm -> m

rccni_p = rccni/exp(dlnr)

if(rcrit.lt.rccni_p) then
  factor = 1.0
else if(rcrit.lt.rccni) then
  factor = (log(rccni)-log(rcrit))/(log(rccni)-log(rccni_p))
else
  factor = 0.0
  return
end if

do
  dmr  = fccni*factor*mass_hydro(1)/rhoa
  qq1  = qq - dmr
  tt1  = tt + dmr*Lv/Cp
  supw = (pp*qq1/(eps+qq1))/esw(tt1)-1.0
  if(supw.ge.0.0) exit
  factor = factor * 0.95
end do

qq = qq1
tt = tt1

ffwi  = ffwi + fccni*factor
fccni = fccni*(1.0-factor)

end subroutine water_nucl

!----------------------------------------------------------------------------
subroutine snow_breakup(ffi, lwfi, rfi)
implicit none

real, dimension(nkr), intent(inout) :: ffi, lwfi, rfi
real, dimension(kr_snow1:kr_snow2), parameter :: &
        factor = (/0.05, 0.1, 0.2, 0.4/)
real, dimension(kr_snow1-1:nkr) :: m1, l1, r1
real :: mm, ll, rr
integer :: n

m1 = ffi(kr_snow1-1:nkr)*mass_hydro(kr_snow1-1:nkr)
l1 = m1*lwfi(kr_snow1-1:nkr)
r1 = (m1-l1)*rfi(kr_snow1-1:nkr)

! zero out from kr_snow2+1 to nkr
mm = sum(m1(kr_snow2:nkr))
ll = sum(l1(kr_snow2:nkr))
rr = sum(r1(kr_snow2:nkr))

m1(kr_snow2) = mm
l1(kr_snow2) = ll
r1(kr_snow2) = rr

 ffi(kr_snow2+1:nkr) = 0.0
lwfi(kr_snow2+1:nkr) = 0.0
 rfi(kr_snow2+1:nkr) = 0.0

! from kr_snow2 to kr_snow1 by -1
do n=kr_snow2, kr_snow1, -1
  m1(n-1) = m1(n-1)+factor(n)*m1(n)
  l1(n-1) = l1(n-1)+factor(n)*l1(n)
  r1(n-1) = r1(n-1)+factor(n)*r1(n)

  m1(n) = (1.0-factor(n))*m1(n)
  l1(n) = (1.0-factor(n))*l1(n)
  r1(n) = (1.0-factor(n))*r1(n)

  ffi(n) = m1(n)/mass_hydro(n)
  if(m1(n).gt.0.0) then
    lwfi(n) = l1(n)/m1(n)
    lwfi(n) = min(max(lwfi(n),0.0),1.0)
  else
    lwfi(n) = 0.0
  end if
  if(m1(n).gt.0.0.and.l1(n).lt.m1(n)) then
    rfi(n) = r1(n)/(m1(n)-l1(n))
    rfi(n) = min(max(rfi(n),0.0),1.0)
  else
    rfi(n) = 0.0
  end if
end do

end subroutine snow_breakup

!----------------------------------------------------------------------------
subroutine vt_pressure(ppi)
implicit none

real, intent(in) :: ppi

integer :: m, n

do m=1, nhydro
do n=1, nkr

  vt_ref(m,n) = bulkvt(m,n)*sqrt(p1000hpa/ppi)

  vt_hydro(m,n) = vt_ref(m,n)

end do
end do

end subroutine vt_pressure

!----------------------------------------------------------------------------
subroutine vt_wet(lwf)
implicit none

real, dimension(nhydro, nkr), intent(in) :: lwf
integer :: m, n

do m=SNOW, HAIL
do n=1, nkr

  vt_hydro(m,n) = vt_hydro(m,n) + lwf(m,n)*(vt_hydro(WATER,n)-vt_hydro(m,n))

end do
end do

end subroutine vt_wet

!----------------------------------------------------------------------------
subroutine hydrometeor_density(lwf, rf)
use ieee_arithmetic
implicit none

real, dimension(nhydro, nkr), intent(in) :: lwf, rf

real :: lwfi, rfi, dry_snow, rho_agg, dm
integer :: n1, n2

do n1=1, nhydro
do n2=1, nkr

  rho_hydro(n1,n2) = bulkrho(n1,n2)

end do
end do

! snow density considering lwf and rf
do n1=1, nkr

  lwfi = lwf(SNOW,n1)
   rfi =  rf(SNOW,n1)

  dry_snow = mass_hydro(n1)*(1.0-lwfi)*(1.0-rfi)

  if(dry_snow.lt.mass_hydro(1)) then
    rho_agg = bulkrho(SNOW,1)
  else if(dry_snow.ge.mass_hydro(nkr)) then
    rho_agg = bulkrho(SNOW,nkr)
  else
    do n2=1, nkr-1
      if(dry_snow.ge.mass_hydro(n2).and.dry_snow.lt.mass_hydro(n2+1)) then
        dm = mass_hydro(n2+1)-mass_hydro(n2)
        rho_agg = bulkrho(SNOW,n2  )*(mass_hydro(n2+1)-dry_snow)/dm &
                 +bulkrho(SNOW,n2+1)*(dry_snow-mass_hydro(n2)  )/dm
        exit
      end if
    end do
  end if

  rho_hydro(SNOW,n1) = (rho_water*rho_rim*rho_agg) &
                      /(rho_rim*rho_agg*lwfi &
                       +rho_agg*rho_water*(1.0-lwfi)*rfi &
                       +rho_water*rho_rim*(1.0-lwfi)*(1.0-rfi))

  rho_ice_snow(n1) = (rho_rim*rho_agg)/(rho_rim*(1.0-rfi)+rho_agg*rfi)

  if(inan.eq.1) then
    if(ieee_is_nan(rho_hydro(SNOW,n1)).or.ieee_is_nan(rho_ice_snow(n1))) then
      write(6,*) n1, lwfi, rfi
      write(6,*) rho_hydro(SNOW,n1), rho_ice_snow(n1)
      call wrf_error_fatal("nan detected in hydrometeor_density")
    end if
  end if

end do

end subroutine hydrometeor_density

!----------------------------------------------------------------------------
subroutine check_validity(ff, lwf, rf, subname)
use ieee_arithmetic
implicit none

real, dimension(nhydro,nkr), intent(inout) :: ff, lwf, rf
character(len=*), intent(in) :: subname

integer :: m, n

! NaN test
if(inan.eq.1) then
  do m=1, nhydro
  do n=1, nkr
    if(ieee_is_nan(ff(m,n))) then
      write(6,*) m, n, trim(subname)
      call wrf_error_fatal("NaN detected in ff")
    end if

    if(ieee_is_nan(lwf(m,n))) then
      write(6,*) m, n, trim(subname)
      call wrf_error_fatal("NaN detected in lwf")
    end if

    if(ieee_is_nan(rf(m,n))) then
      write(6,*) m, n, trim(subname)
      call wrf_error_fatal("NaN detected in rf")
    end if
  end do
  end do
end if

do m=1, nhydro
do n=1, nkr
  if(ff(m,n).lt.0.0) then
    if(ff(m,n).lt.-0.0001) then
      write(6,*) m, n, ff(m,n), trim(subname)
      call wrf_error_fatal("negative ff detect")
    else
      ff(m,n) = 0.0
    end if
  end if

  if(lwf(m,n).lt.0.0) then
    if(lwf(m,n).lt.-0.0001) then
      write(6,*) m, n, lwf(m,n), trim(subname)
      call wrf_error_fatal("negative lwf detected during check")
    else
      lwf(m,n) = 0.0
    end if
  end if

  if(lwf(m,n).gt.1.0) then
    if(lwf(m,n).gt.1.0001) then
      write(6,*) m, n, lwf(m,n), trim(subname)
      call wrf_error_fatal("large lwf detected during check")
    else
      lwf(m,n) = 1.0
    end if
  end if

  if(rf(m,n).lt.0.0) then
    if(rf(m,n).lt.-0.0001) then
      write(6,*) m, n, rf(m,n), trim(subname)
      call wrf_error_fatal("negative rf detected during check")
    else
      rf(m,n) = 0.0
    end if
  end if

  if(rf(m,n).gt.1.0) then
    if(rf(m,n).gt.1.0001) then
      write(6,*) m, n, rf(m,n), trim(subname)
      call wrf_error_fatal("large rf detected during check")
    else
      rf(m,n) = 1.0
    end if
  end if
end do
end do

end subroutine check_validity

!===================================================================!
end module module_mp_sbm
