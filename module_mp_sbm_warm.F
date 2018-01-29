module module_mp_sbm_warm
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
integer, private, save :: isubsidence
integer, private, save :: iremap       ! 1 for KO, 2 for new
integer, private, save :: ibreakup     ! 1 for implicit, 2 for explicit
integer, private, save :: ikernel      ! 1 for Pinsky, 2 for Hall
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

integer, parameter, private :: scal = 1
real,    parameter, private :: dlnr = log(2.0)/3.0/scal

real, dimension(nkr), save, private :: mass_hydro
real, dimension(nkr), save, private :: vt_hydro
real, dimension(nkr), save, private :: rho_hydro
real, dimension(nkr), save, private :: r_hydro
real, dimension(nkr), save, private :: bulkvt
real, dimension(nkr,nkr), save, private :: ckxx
real, dimension(3,nkr,nkr), save, private :: ywxx

real, parameter, private :: rccn_min = 0.003e-4 ! cm
real, parameter, private :: smax     = 0.20     ! 120% of RH
real, parameter, private :: z0_ccn   = 2000.    ! 2.0 km
real, parameter, private :: dz_ccn   = 2000.    ! 2.0 km

real, parameter, private :: t00      = 273.15
real, parameter, private :: p00      = 101325.
real, parameter, private :: p1000hpa = 100000.

real, parameter, private :: Rd  = 287.
real, parameter, private :: Rv  = 461.51
real, parameter, private :: eps = Rd/Rv
real, parameter, private :: Cp  = 1004.7
real, parameter, private :: Ls  = 2.8345e6      ! latent heat for sublimation
real, parameter, private :: Lv  = 2.5008e6      ! latent heat for evaporation
real, parameter, private :: b_w = Lv/Rv
real, parameter, private :: a_w = b_w/t00
real, parameter, private :: b_i = Ls/Rv
real, parameter, private :: a_i = b_i/t00

real, parameter, private :: pi    = 3.141593
real, parameter, private :: grav  = 9.81

real, parameter, private :: rho_water = 1000.

real,    parameter, private :: tt_coal = t00 - 40.0
integer, dimension(nkr,nkr), save, private :: ima
real,    dimension(nkr,nkr), save, private :: cour

integer, parameter, private :: krbreak = 28
real,    dimension(krbreak,krbreak,krbreak), save, private :: pkij
real,    dimension(krbreak,krbreak),         save, private :: qkj

real,    parameter, private :: ventpl_max = 5.0

real, parameter, private :: ro_solute = 2.16
real, parameter, private :: akoe = 3.3e-05
real, parameter, private :: bkoe = 2.0*4.3/(22.9+35.5)*(4.0/3.0)*pi*ro_solute

real, parameter, private :: tt_nucl_drop_min = -45.0

real, parameter, private :: tau_ccn = 3600.0 ! relaxation time is i h

real, parameter, private :: dt_cond_ref = 1.0
real, save,      private :: dt_cond
integer, save,   private :: n_cond

real, parameter, private :: Dsub = 2.2124e-6    ! 3.75e-6 s-1
real, parameter, private :: vsub_max = 5.0e-3 ! 0.5 cm s-1

real, save, private :: turbll_grid(13,13,11)

contains
!--------------------------------------------------------------------------
subroutine sbm_warm(th_phy, qv, pi_phy, p_phy,    &
                    dx, dy, dz8w,                 &
                    tke_opt, tke3d, n2_3d,        &
                    xland, ivgtyp, isurban,       &
                    scalar, n_scalar,             &
                    qc, qr, nc, nr, na,           &
                    rainnc, rainncv,              &
                    dTn3d, dTce3d,                &
                    eps_phy, rey_phy,             &
                    its, ite, jts, jte, kts, kte, &
                    ims, ime, jms, jme, kms, kme, &
                    ids, ide, jds, jde, kds, kde, &
                    p_start, dt, itimestep        )
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
      dz8w                                                 &
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
    , nc                                                   &
    , nr                                                   &
    , na                                                   &
    , dTn3d                                                &
    , dTce3d                                               &
    , eps_phy                                              &
    , rey_phy
real, dimension(ims:ime,kms:kme,jms:jme,n_scalar), intent(inout) :: &
      scalar
real, dimension(ims:ime,jms:jme), intent(in   ) :: &
      xland
real, dimension(ims:ime,jms:jme), intent(inout) :: &
      rainnc
real, dimension(ims:ime,jms:jme), intent(out  ) :: &
      rainncv
integer, dimension(ims:ime,jms:jme), intent(in) :: &
      ivgtyp

! Local variables
integer :: i, j, k, m, n, ii, kk
real,    dimension(ite-its+1) :: rain1d, xland1d
integer, dimension(ite-its+1) :: ivgtyp1d
real,    dimension(ite-its+1,kte-kts+1) ::         &
      th2d, qv2d, pp2d, pi2d, hgt2d, dz2d          &
    , qc2d, qr2d                                   &
    , nc2d, nr2d, na2d                             &
    , dTn2d, dTce2d                                &
    , tke2d, n2_2d, eps2d, rey2d
real,    dimension(ite-its+1,kte-kts+1,nkr) :: fccn
real,    dimension(ite-its+1,kte-kts+1,nkr) :: ffx

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

! From 3-D to 2-D
do j=jts, jte

  ! set model height
  do i=its, ite

    ii = i-its+1
    dz2d (ii,1) = dz8w(i,kts,j)
    hgt2d(ii,1) = dz2d(ii,1)/2.0

    do k=kts+1, kte

      kk = k-kts+1
      dz2d (ii,kk) = dz8w(i,k,j)
      hgt2d(ii,kk) = hgt2d(ii,kk-1) + (dz2d(ii,kk-1)+dz2d(ii,kk))/2.0 

    end do
  end do

  do k=kts, kte
  do i=its, ite

    ii = i-its+1
    kk = k-kts+1

    th2d(ii,kk) = th_phy(i,k,j)
    qv2d(ii,kk) = max(qv(i,k,j),0.)
    pp2d(ii,kk) = p_phy(i,k,j)
    pi2d(ii,kk) = pi_phy(i,k,j)

    tke2d(ii,kk) = tke3d(i,k,j)
    n2_2d(ii,kk) = n2_3d(i,k,j)

    xland1d (ii) = xland (i,j)
    ivgtyp1d(ii) = ivgtyp(i,j)

! ff for hydrometeors (kg kg-1)
    do n=1, nkr
      ffx(ii,kk,n) = scalar(i,k,j,(p_start-1)+n)
    end do

! ff for aerosol (# kg-1)
    do n=1, nkr
      fccn(ii,kk,n) = scalar(i,k,j,(p_start-1)+nkr+n)
    end do

    if(inan.eq.1) then
      if(ieee_is_nan(th2d(ii,kk)).or. &
         ieee_is_nan(qv2d(ii,kk)).or. &
         ieee_is_nan(dz2d(ii,kk)).or. &
         ieee_is_nan(sum(ffx(ii,kk,:))).or. &
         ieee_is_nan(sum(fccn(ii,kk,:)))) then

         write(6,*) i, k, j
         write(6,*) th2d(ii,kk), qv2d(ii,kk), dz2d(ii,kk)
         write(6,*) sum(ffx(ii,kk,:)), sum(fccn(ii,kk,:))
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
  call sbm2d_warm(ite-its+1, kte-kts+1, dx, dy, dt,         &
                  th2d, qv2d, pp2d, pi2d, hgt2d, dz2d,      &
                  xland1d, ivgtyp1d, isurban,               &
                  ffx, fccn,                                &
                  qc2d, qr2d,                               &
                  nc2d, nr2d, na2d,                         &
                  rain1d,                                   &
                  dTn2d, dTce2d,                            &
                  tke2d, n2_2d, eps2d, rey2d,               &
                  its, j, itimestep                         )
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
         ieee_is_nan(nc2d(ii,kk)).or. &
         ieee_is_nan(nr2d(ii,kk)).or. &
         ieee_is_nan(na2d(ii,kk)).or. &
         ieee_is_nan(dTn2d(ii,kk)).or. &
         ieee_is_nan(dTce2d(ii,kk)).or. &
         ieee_is_nan(sum(ffx(ii,kk,:))).or. &
         ieee_is_nan(sum(fccn(ii,kk,:))).or. &
         ieee_is_nan(rain1d(ii))) then

         write(6,*) i, k, j
         write(6,*) th2d(ii,kk), qv2d(ii,kk), qc2d(ii,kk), qr2d(ii,kk)
         write(6,*) nc2d(ii,kk), nr2d(ii,kk), na2d(ii,kk)
         write(6,*) dTn2d(ii,kk), dTce2d(ii,kk)
         write(6,*) sum(ffx(ii,kk,:)), sum(fccn(ii,kk,:))
         write(6,*) rain1d(ii)
         call wrf_error_fatal("nan from mp.")
      end if
    end if

    th_phy(i,k,j) = th2d(ii,kk)

    qv(i,k,j) = qv2d(ii,kk)
    qc(i,k,j) = qc2d(ii,kk)
    qr(i,k,j) = qr2d(ii,kk)
    nc(i,k,j) = nc2d(ii,kk)
    nr(i,k,j) = nr2d(ii,kk)
    na(i,k,j) = na2d(ii,kk)

    dTn3d(i,k,j)  = dTn2d(ii,kk)
    dTce3d(i,k,j) = dTce2d(ii,kk)

    eps_phy(i,k,j) = eps2d(ii,kk)
    rey_phy(i,k,j) = rey2d(ii,kk)

    do n=1, nkr
      scalar(i,k,j,(p_start-1)+n) = ffx(ii,kk,n)
    end do

    do n=1, nkr
      scalar(i,k,j,(p_start-1)+nkr+n) = fccn(ii,kk,n)
    end do

    if(k==kts) then
      rainncv(i,j) = rain1d(ii)
      rainnc(i,j)  = rainnc(i,j) + rainncv(i,j)
    end if

  end do
  end do

end do

end subroutine sbm_warm

!--------------------------------------------------------------------------
subroutine sbm2d_warm(ni, nk, dx, dy, dt,           &
                      th, qq, pp, pii, alt2d, dz2d, &
                      xland1d, ivgtyp1d, isurban,   &
                      ffx, fccnx,                   &
                      qcloud, qrain,                &
                      ncloud, nrain,                &
                      naerosol,                     &
                      rain_x,                       &
                      dT_n, dT_ce,                  &
                      tke, n2, eps2d, rey2d,        &
                      i_id, j_id, itimestep         )
use ieee_arithmetic
implicit none

integer, intent(in) :: ni, nk, i_id, j_id, isurban, itimestep
real,    intent(in) :: dx, dy, dt

real,    dimension(ni    ),      intent(out  ) :: rain_x
real,    dimension(ni    ),      intent(in   ) :: xland1d
integer, dimension(ni    ),      intent(in   ) :: ivgtyp1d
real,    dimension(ni, nk),      intent(in   ) :: pp, pii, tke, n2, alt2d, dz2d
real,    dimension(ni, nk),      intent(inout) :: th, qq
real,    dimension(ni, nk),      intent(out  ) ::                 &
   qcloud, qrain, ncloud, nrain, naerosol, dT_n, dT_ce,           &
   eps2d, rey2d
real,    dimension(ni, nk, nkr), intent(inout) :: ffx
real,    dimension(ni, nk, nkr), intent(inout) :: fccnx

! local variables
integer :: i, j, k, m, n, nn1, nn2, nn1s, nn1e, nn2s, nn2e, iter, ivgtype
real, dimension(ni,nk) :: tt
real, dimension(nkr) :: ff, g, dg, fccn, m_new, ff_hydro1, vt_sfc
real, dimension(ni) :: precp
real, dimension(nk) :: rhoaz, altz, dzz, vsub, varz, vtz
real :: tt1, pp1, qq1, rhoa1, tke1, n2_1, eps1, rey1, &
        ffi, massi, vti, capi,                        &
        fccni, rccni, ffwi, m_newi,                   &
        rain_mass_init, rain_mass_finl, sum_temp,     &
        T_old, mass_ref, ak, dz, island, alt

! initialization of output array
rain_x = 0.
dT_n = 0.; dT_ce = 0.

do i=1, ni
do k=1, nk

  tt(i,k) = th(i,k)*pii(i,k)
  rhoa1 = pp(i,k)/tt(i,k)/Rd

  ! unit conversion ---------------------------------------------------------
  ! ffx : kg kg-1, ff : # m-3
  ! fccnx : # kg-1, fccn : # m-3

  ! caution: at the initialization, the unit of fccn is # m-3, not # kg-1.
  ff = max(ffx(i,k,:)/mass_hydro*rhoa1,0.0)
  if(itimestep.gt.1) then
    fccn = max(fccnx(i,k,:)*rhoa1,0.0)
  else
    fccn = max(fccnx(i,k,:),0.0)
  end if

  if(inan.eq.1) then
    if(ieee_is_nan(tt(i,k)).or. &
       ieee_is_nan(pp(i,k)).or. &
       ieee_is_nan(qq(i,k))) then
      write(6,*) i_id+i-1, j_id, k
      write(6,*) tt(i,k), pp(i,k), qq(i,k)
      call wrf_error_fatal("nan in tt before integration")
    end if
  end if

  call check_validity(ff, "start")

  call vt_pressure(pp(i,k))

  ! collision ---------------------------------------------------------------
  tt1   = tt(i,k)
  pp1   = pp(i,k)
  rhoa1 = pp1/tt1/Rd
  tke1  = tke(i,k)
  n2_1  = n2(i,k)
  dz    = dz2d(i,k)

  call kernel_pressure(pp1)

  call kernel_mod(tke1, n2_1, dx, dy, dz, rhoa1, eps1, rey1)

  eps2d(i,k) = eps1
  rey2d(i,k) = rey1

  g = ff*mass_hydro/dlnr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! water-water: sepatated collision of c-c, c-r, r-r

! cloud-cloud
  rain_mass_init = sum(g(kr_droplet+1:))*dlnr
  do nn1 =     1, kr_droplet
  do nn2 = nn1+1, kr_droplet
    if(g(nn1).gt.0.0.and.g(nn2).gt.0.0.and.ckxx(nn1,nn2).gt.0.0) then
      if(icollision.eq.1) then
        call collision(g, nn1, nn2)
      else if(icollision.eq.2) then
        call collision_poisson(g, nn1, nn2)
      else
        call wrf_error_fatal("something wrong in icollision")
      end if
    end if
  end do
  end do
  rain_mass_finl = sum(g(kr_droplet+1:))*dlnr

! cloud-rain
  rain_mass_init = rain_mass_finl
  do nn1 =            1, kr_droplet
  do nn2 = kr_droplet+1, nkr
    if(g(nn1).gt.0.0.and.g(nn2).gt.0.0.and.ckxx(nn1,nn2).gt.0.0) then
      if(icollision.eq.1) then
        call collision(g, nn1, nn2)
      else if(icollision.eq.2) then
        call collision_poisson(g, nn1, nn2)
      else
        call wrf_error_fatal("something wrong in icollision")
      end if
    end if
  end do
  end do
  rain_mass_finl = sum(g(kr_droplet+1:))*dlnr

! rain-rain
  do nn1 = kr_droplet+1, nkr
  do nn2 =        nn1+1, nkr
    if(g(nn1).gt.0.0.and.g(nn2).gt.0.0.and.ckxx(nn1,nn2).gt.0.0) then
      if(icollision.eq.1) then
        call collision(g, nn1, nn2)
      else if(icollision.eq.2) then
        call collision_poisson(g, nn1, nn2)
      else
        call wrf_error_fatal("something wrong in icollision")
      end if
    end if
  end do
  end do

  ! breakup by Bleck (1970)
  ff_hydro1 = g
  call breakup_water(ff_hydro1,dt)
  g = ff_hydro1

  ff = g/mass_hydro*dlnr

  call check_validity(ff, "collision")

  ! vapor diffusion ---------------------------------------------------------
  do iter=1, n_cond

    if(maxval(ff).le.0.0) exit

    tt1   = tt(i,k)
    qq1   = qq(i,k)
    pp1   = pp(i,k)
    rhoa1 = pp1/tt1/Rd

    ff_hydro1 = ff/rhoa1

    call diffuse(tt1, pp1, qq1, rhoa1, dt_cond, ff_hydro1, m_new)

    tt(i,k) = tt1
    qq(i,k) = qq1
    dT_ce(i,k) = dT_ce(i,k) + (tt1 - tt(i,k))

    if(maxval(ff).gt.0.0) then
      call remapping(m_new,ff_hydro1)
    end if

    ff = ff_hydro1*rhoa1

  end do

  call check_validity(ff, "diffusion")

  ! nucleation --------------------------------------------------------------
  T_old = tt(i,k)
  rccni = r_hydro(1)

  do n=nkr, 1, -1 ! large aerosol first

    fccni = fccn(n)

    if(fccni.gt.0.0) then
      tt1   = tt(i,k)
      qq1   = qq(i,k)
      pp1   = pp(i,k)
      rhoa1 = pp(i,k)/tt(i,k)/Rd
      ffwi  = ff(max(n-35,1))

      call water_nucl(tt1, pp1 ,qq1 ,rhoa1 ,fccni, rccni, ffwi)

      tt(i,k) = tt1
      qq(i,k) = qq1

      fccn(n) = fccni
      ff(max(n-35,1)) = ffwi
    end if

    rccni = rccni/exp(dlnr)

  end do

  dT_n(i,k) = tt(i,k) - T_old

  call check_validity(ff, "nucleation")

  ! aerosol replenishment by relaxation (Hyunho) ----------------------------
  island = xland1d(i)
  ivgtype = ivgtyp1d(i)
  alt = alt2d(i,k)

  call fccn_init(alt, island, ivgtype, isurban, ff_hydro1)

  do n=1, nkr
    fccn(n) = fccn(n) + (ff_hydro1(n) - fccn(n))/tau_ccn * dt
  end do

  ! unit conversion ---------------------------------------------------------
  rhoa1        = pp(i,k)/tt(i,k)/Rd
  ffx(i,k,:)   = ff*mass_hydro/rhoa1
  fccnx(i,k,:) = fccn/rhoa1

  th(i,k) = tt(i,k)/pii(i,k)

end do
end do

! sedimentation -----------------------------------------------------------
do i=1, ni

  rhoaz = pp(i,:)/tt(i,:)/Rd
  altz = alt2d(i,:)
  dzz = dz2d(i,:)

  if(isubsidence.eq.1) then
    vsub = min(Dsub*altz,vsub_max)
  else
    vsub = 0.0
  end if

  do n=1, nkr
    vtz = bulkvt(n)*sqrt(p1000hpa/pp(i,:))
    vt_sfc(n) = vtz(1)

    varz = ffx(i,:,n)
    call sediment(varz, 'h', vtz, vsub, rhoaz, altz, dzz, nk, dt)
    ffx(i,:,n) = varz
  end do

  if(isubsidence.eq.1) then
    vtz = 0.0
    do n=1, nkr
      varz = fccnx(i,:,n)
      call sediment(varz, 's', vtz, vsub, rhoaz, altz, dzz, nk, dt)
      fccnx(i,:,n) = varz
    end do

    varz = th(i,:)
    call sediment(varz, 's', vtz, vsub, rhoaz, altz, dzz, nk, dt)
    th(i,:) = varz

    varz = qq(i,:)
    call sediment(varz, 's', vtz, vsub, rhoaz, altz, dzz, nk, dt)
    qq(i,:) = varz
  end if

  rhoa1 = pp(i,1)/tt(i,1)/Rd

  do n=1, nkr
    precp(i) = ffx(i,1,n)/rho_hydro(n)*vt_sfc(n)*rhoa1*dt*1000.
    rain_x(i) = rain_x(i) + precp(i)
  end do

  if(inan.eq.1) then
    if(ieee_is_nan(precp(i))) then
      write(6,*) "in precipitation", i, n
      write(6,*) ffx(i,1,n),rho_hydro(n),vt_hydro(n)
      call wrf_error_fatal("nan detect")
    end if
  end if

end do

! large-scale forcing from RICO (vanZanten et al. 2011)
!do i=1, ni
!do k=1, nk
!
!  alt = alt2d(i,k)
!
!  th(i,k) = th(i,k) - 2.5*(dt/86400.0)
!  qq(i,k) = qq(i,k) + min((1.3456/2980.0*alt)-1.0,0.3456)*1.0e-3*(dt/86400.0)
!
!end do
!end do

! WRF output
do i=1, ni
do k=1, nk

  qcloud  (i,k) = sum(ffx(i,k,1:kr_droplet))
  qrain   (i,k) = sum(ffx(i,k,kr_droplet+1:nkr))

  ncloud  (i,k) = sum(ffx(i,k,1:kr_droplet)/mass_hydro(1:kr_droplet))
  nrain   (i,k) = sum(ffx(i,k,kr_droplet+1:nkr)/mass_hydro(kr_droplet+1:nkr))

  naerosol(i,k) = sum(fccnx(i,k,:))

end do
end do

end subroutine sbm2d_warm

!==============================================================================
subroutine sediment(varz, itype, vtz, vsub, rhoaz, altz, dzz, nz, dt)
implicit none

integer, intent(in) :: nz
real, intent(in) :: dt
character(len=1), intent(in) :: itype
real, dimension(nz), intent(in) :: vtz, vsub, rhoaz, altz, dzz
real, dimension(nz), intent(inout) :: varz

real :: tf, dt_fall
integer :: iter, it, k
real, dimension(nz)   :: dzz_dt
real, dimension(nz,2) :: courant, flux

tf = minval(altz/(vtz+vsub+1.0e-10))
iter = int((2.0*dt)/tf)+1
dt_fall = dt/iter
dzz_dt = dt_fall/dzz

do k=1, nz-1
  courant(k,1) =                     vsub(k)  *dzz_dt(k)
  courant(k,2) = rhoaz(k+1)/rhoaz(k)*vsub(k+1)*dzz_dt(k)
end do
courant(nz,:) = vsub(nz)*dzz_dt(nz)

do it=1, iter
  ! hydrometeors that have terminal velocity
  if(itype.eq.'h') then
    do k=1, nz-1
      flux(k,1) =                     (vtz(k)  +vsub(k)  )*dzz_dt(k)*varz(k)
      flux(k,2) = rhoaz(k+1)/rhoaz(k)*(vtz(k+1)+vsub(k+1))*dzz_dt(k)*varz(k+1)
    end do
    flux(nz,1) = (vtz(nz)+vsub(nz))*dzz_dt(nz)*varz(nz)
    flux(nz,2) = 0.0

  ! scalars that hang aloft
  else if(itype.eq.'s') then
    do k=1, nz-1
      flux(k,1) = courant(k,1)*varz(k)
      flux(k,2) = courant(k,2)*varz(k+1)
    end do
    flux(nz,:) = 0.0

  ! error handle
  else
    call wrf_error_fatal("something wrong in sedimentation")
  end if

  do k=1, nz
    varz(k) = (varz(k)-flux(k,1)+flux(k,2))/(1.0-courant(k,1)+courant(k,2))
  end do
end do

end subroutine sediment

!==============================================================================
subroutine sbminit_warm(ids, ide, jds, jde, kds, kde,      &
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

! local parameters for terminal velocity (Beard 1976)
real, parameter :: l0 = 6.62e-6
real, parameter :: h0 = 1.818e-4
real, parameter :: rhoa = 1.189e-3
real, parameter :: drho = 1.0-rhoa
real, parameter :: sigma = 72.8
real, parameter ::       &
      b10 = -0.318657e+1 &
    , b11 =  0.992696    &
    , b12 = -0.153193e-2 &
    , b13 = -0.987059e-3 &
    , b14 = -0.578878e-3 &
    , b15 =  0.855176e-4 &
    , b16 = -0.327815e-5
real, parameter ::       &
      b20 = -0.500015e+1 &
    , b21 =  0.523778e+1 &
    , b22 = -0.204914e+1 &
    , b23 =  0.475294    &
    , b24 = -0.542819e-1 &
    , b25 =  0.238449e-2
real, parameter :: c1 = drho*(grav*100.)/(18.*h0)
real, parameter :: c2 = 4.*rhoa*drho*(grav*100.)/(3.*h0**2)
real, parameter :: c3 = 4.*drho*(grav*100.)/(3.*sigma)
real, parameter :: Np = ((sigma**3)*(rhoa**2)/((h0**4)*drho*(grav*100.)))**(1.0/6.0)

real, dimension(21) :: rl_turb
real, dimension(21,21,11) :: turbll_text1
real, dimension(13,13,11) :: turbll_text2

! Hall kernel
real, dimension(nkr,nkr) :: e_int, kernel_raw
real, dimension(nkr) :: r_ref
real, dimension(15) :: r0
real, dimension(21) :: rat
real, dimension(15,21) :: ecoll
real :: p, q, rq
integer :: ir, iq, kk

real, dimension(nkr) :: fccni

integer :: hujisbm_unit1
logical :: opened
logical, external :: wrf_dm_on_monitor

integer :: i, j, k, kr, kp, ip, jp, kq, jq, in, jn, m1, m2, ivgtype
real :: x, y, Csc, Nre, d_cm, mass0, alt, island

ywxx = 0.0

print *,'wrf_dm_on_monitor() =',wrf_dm_on_monitor()

!--- set timestep for diffusion
n_cond = max(nint(dt/dt_cond_ref),1)
dt_cond = dt/n_cond
write(6,*) "n_cond, dt_cond = ", n_cond, dt_cond

!--- mass, density, and radius of water
   r_hydro(1) = 2.0e-6 ! 2 mm
 rho_hydro(1) = 1.0e+3 ! 1 g cm-3
mass_hydro(1) = rho_hydro(1)*4.0/3.0*pi*(r_hydro(1)**3)

do kr=2, nkr
  mass_hydro(kr) = mass_hydro(kr-1)*(2**(1.0/scal)) ! mass doubled every scal bins
   rho_hydro(kr) = 1.0e+3
     r_hydro(kr) = (mass_hydro(kr)/rho_hydro(kr)*3.0/4.0/pi)**(1.0/3.0)
end do

!--- terminal velocity (Beard 1976)
do kr=1, nkr
  d_cm = r_hydro(kr)*1.0e2*2.0 ! diameter in cm
  Csc  = 1.0+2.51*l0/d_cm

  if(d_cm.lt.0.0019) then     ! diameter less than 19 micro m
    bulkvt(kr) = c1*Csc*(d_cm**2)
  else if(d_cm.lt.0.107) then ! diameter less than 1.07 mm
    x = log(c2*d_cm**3)
    y = b10 + b11*x + b12*x**2 + b13*x**3 + b14*x**4 + b15*x**5 + b16*x**6
    Nre = Csc*exp(y)
    bulkvt(kr) = h0*Nre/rhoa/d_cm
  else if(d_cm.lt.0.7) then   ! diameter less than 7 mm
    x = log(c3*d_cm**2*Np)
    y = b20 + b21*x + b22*x**2 + b23*x**3 + b24*x**4 + b25*x**5
    Nre = Np*exp(y)
    bulkvt(kr) = h0*Nre/rhoa/d_cm
  else                        ! diameter larger than 7 mm
    bulkvt(kr) = 920.
  end if

  bulkvt(kr) = bulkvt(kr)*1.0e-2 ! cgs to mks
end do

if ( wrf_dm_on_monitor() ) then
  write(6,*) "kr     mass        radius         rho       term_vel"
  do kr=1, nkr
    write(6,'(I2,1x,4E13.5)') kr, mass_hydro(kr), r_hydro(kr), rho_hydro(kr), bulkvt(kr)
  end do
end if

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
  read(hujisbm_unit1,*) ikernel
  read(hujisbm_unit1,*) iturb
  read(hujisbm_unit1,*) ibreakup
  read(hujisbm_unit1,*) inan
  read(hujisbm_unit1,*) icollision
  read(hujisbm_unit1,*) isubsidence
  close(hujisbm_unit1)
  call wrf_message(' ***** file_config: succesfull *******')
endif
call wrf_dm_bcast_bytes ( n0ccn1 ,      RWORDSIZE )
call wrf_dm_bcast_bytes ( kccn1 ,       RWORDSIZE )
call wrf_dm_bcast_bytes ( n0ccn2 ,      RWORDSIZE )
call wrf_dm_bcast_bytes ( kccn2 ,       RWORDSIZE )
call wrf_dm_bcast_bytes ( urb_ccnf ,    RWORDSIZE )
call wrf_dm_bcast_bytes ( iremap ,      IWORDSIZE )
call wrf_dm_bcast_bytes ( icollision ,  IWORDSIZE )
call wrf_dm_bcast_bytes ( ikernel ,     IWORDSIZE )
call wrf_dm_bcast_bytes ( iturb ,       IWORDSIZE )
call wrf_dm_bcast_bytes ( ibreakup ,    IWORDSIZE )
call wrf_dm_bcast_bytes ( isubsidence , IWORDSIZE )
call wrf_dm_bcast_bytes ( inan ,        IWORDSIZE )

write(6,*) "Read configuration from file"
write(6,*) "n0 over the land:       ", n0ccn1
write(6,*) "k over the land:        ", kccn1
write(6,*) "n0 over the ocean:      ", n0ccn2
write(6,*) "k over the ocean:       ", kccn2
write(6,*) "n0 multiplier for urban:", urb_ccnf
write(6,*) "iremap:     ", iremap
write(6,*) "icollision: ", icollision
write(6,*) "ikernel:    ", ikernel
write(6,*) "iturb:      ", iturb
write(6,*) "ibreakup:   ", ibreakup
write(6,*) "isubsidence:", isubsidence
write(6,*) "inan:       ", inan

! kernels (depending on pressure) :
! water-water for 1000, 750, 500 hPa
if ( wrf_dm_on_monitor() ) then

  if(ikernel.eq.1) then
! Pinsky kernel
    open(unit=hujisbm_unit1,file="ckll_z.asc",status="old",action="read",err=2070)
    read(hujisbm_unit1,900) (((ywxx(k,i,j),i=1,nkr),j=1,nkr),k=1,3)
    close(hujisbm_unit1)
    call wrf_message(' ***** file51: succesfull *******')
  else
! Hall kernel
    open(unit=hujisbm_unit1,file="eff_hall.txt",status="old",action="read",err=2070)
    read(hujisbm_unit1,*) ecoll

    close(hujisbm_unit1)
    call wrf_message(' ***** file51: succesfull *******')

    r0 =  (/ 6.0,  8.0, 10.0,  15.0,  20.0,  25.0,  30.0, 40.0, &
            50.0, 60.0, 70.0, 100.0, 150.0, 200.0, 300.0/)
    rat = (/0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, &
            0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0/)

    r_ref = r_hydro*1.0e6

    do j=1, nkr
    do i=1, j
      do k=2, 15
        if (r_ref(j).le.r0(k).and.r_ref(j).ge.r0(k-1)) then
          ir=k
        else if (r_ref(j).gt.r0(15)) then
          ir=16
        else if (r_ref(j).lt.r0(1)) then
          ir=1
        endif
      enddo

      rq=r_ref(i)/r_ref(j)

      do kk=2,21
        if (rq.le.rat(kk).and.rq.gt.rat(kk-1)) iq=kk
      enddo

      if (ir.lt.16) then
        if (ir.ge.2) then
          p = (r_ref(j)-r0(ir-1))/(r0(ir)-r0(ir-1))
          q = (rq-rat(iq-1))/(rat(iq)-rat(iq-1))
          e_int(j,i)=(1.-p)*(1.-q)*ecoll(ir-1,iq-1) &
                    +p     *(1.-q)*ecoll(ir,  iq-1) &
                    +(1.-p)*q     *ecoll(ir-1,iq  ) &
                    +p     *q     *ecoll(ir,  iq  )
        else
          q = (rq-rat(iq-1))/(rat(iq)-rat(iq-1))
          e_int(j,i) = (1.-q)*ecoll(1,iq-1)+q*ecoll(1,iq)
        endif
      else
        q = (rq-rat(iq-1))/(rat(iq)-rat(iq-1))
        e_int(j,i) = min((1.-q)*ecoll(15,iq-1)+q*ecoll(15,iq),1.0)
      endif
      e_int(i,j) = e_int(j,i)
    enddo
    enddo

    do i=1, nkr
    do j=1, nkr
      kernel_raw(i,j) = pi*(r_hydro(i)+r_hydro(j))**2 &
                          *abs(bulkvt(i)-bulkvt(j))   &
                          *e_int(i,j)                 &
                          *1.0e6                          ! mks to cgs
    end do
    end do

    do i=1, nkr
    do j=1, nkr
      ip = max(i-1,1  )
      in = min(i+1,nkr)
      jp = max(j-1,1  )
      jn = min(j+1,nkr)
      ywxx(:,i,j) = ( kernel_raw(ip,j )     &
                    + kernel_raw(i ,jp)     &
                    + kernel_raw(i ,j )*4.0 &
                    + kernel_raw(in,j )     &
                    + kernel_raw(i ,jn) )/8.0

      if(i.eq.j) ywxx(:,i,j) = 0.5*ywxx(:,i,j)
    end do
    end do

  end if

end if
call wrf_dm_bcast_bytes ( ywxx , size ( ywxx ) * RWORDSIZE )

! collision enhancement by turbulence
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
pkij = pkij * 1.0e-15 ! g3 cm3 s-1 -> kg3 m3 s-1
 qkj =  qkj * 1.0e-15 ! g3 cm3 s-1 -> kg3 m3 s-1

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
    scalar(i,k,j,1+nkr+1:1+nkr+nkr) = fccni

  end do
  end do
  end do
 
  if(wrf_dm_on_monitor()) then
    write(6,*) "initial ccn", xland(its,jts),ivgtyp(its,jts),isurban
    write(6,*) scalar(its,kts,jts,1+nkr+1:1+2*nkr)
    write(6,*) "initial height", kts, kte
    write(6,*) hgt(its,kts:kte,jts)
  end if

end if

ywxx = ywxx * 1.0e-6 * dlnr * dt ! cm-3 s-1 --> m-3 s-1

!---- turbulence-induced collision enhancement for water-water lookup table
do i=1, 21
  rl_turb(i)=float(i)*1.0e-6
end do

do k=1, 11
  do i=1, 13
    x=r_hydro(i)
    if(x.lt.rl_turb( 1)) x=rl_turb( 1)
    if(x.gt.rl_turb(21)) x=rl_turb(21)
    do j=1, 13
      y=r_hydro(j)
      if(y.lt.rl_turb( 1)) y=rl_turb( 1)
      if(y.gt.rl_turb(21)) y=rl_turb(21)
      turbll_grid(i,j,k)=f(x,y,rl_turb,rl_turb,turbll_text1(:,:,k),21,21)
      turbll_grid(i,j,k)=turbll_grid(i,j,k)*turbll_text2(i,j,k)
    end do
  end do
end do

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
end do
end do

! breakup collision table
do i=1, nkr
do j=1, nkr
  ywxx(:,i,j) = ywxx(:,i,j) * ecoalmass(mass_hydro(i),mass_hydro(j))
end do
end do

return

 900 format(6e13.5)
2070 continue

call wrf_error_fatal("error opening hujisbm_data")

end subroutine sbminit_warm

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

r0=r_hydro(1)*100. ! m -> cm

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

!----------------------------------------------------------------------------
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
real function esw(tt)
implicit none

real, intent(in) :: tt

esw = 611.*exp(a_w - b_w/tt)

end function esw

!----------------------------------------------------------------------------
real function esi(tt)
implicit none

real, intent(in) :: tt

esi = 611.*exp(a_i - b_i/tt)

end function esi

!===================================================================
subroutine kernel_pressure(pp)
implicit none

real, intent(in) :: pp

real :: p_1, p_2, p_3
integer :: i, j

p_1 = 100000.; p_2 = 75000.; p_3 = 50000.
 
do i=1, nkr
do j=1, nkr
  if     (pp.ge.p_1) then
    ckxx(i,j)=ywxx(1,i,j)
  else if(pp.ge.p_2) then
    ckxx(i,j)=ywxx(2,i,j)*(p_1-pp)/(p_1-p_2) &
             +ywxx(1,i,j)*(pp-p_2)/(p_1-p_2)
  else if(pp.ge.p_3) then
    ckxx(i,j)=ywxx(3,i,j)*(p_2-pp)/(p_2-p_3) &
             +ywxx(2,i,j)*(pp-p_3)/(p_2-p_3)
  else
    ckxx(i,j)=ywxx(3,i,j)
  end if

end do
end do

end subroutine kernel_pressure

!=================================================================
subroutine kernel_mod(tke, n2, dx, dy, dz, rho, eps1, rey1)
implicit none

real, intent(in) :: tke, n2, dx, dy, dz, rho
real, intent(out) :: eps1, rey1

call tice(tke, n2, dx, dy, dz, rho, eps1, rey1)

end subroutine kernel_mod

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

real :: extra_tice, eps0, rey0, rst, c1, c2, c3, c4, &
        len0, tketmp, ce1, ce2, coefc

! extrapolation of kernel using only 1 point
extra_tice(rst,epsi,eps0,reyi,rey0) = rst*(epsi/eps0)**1.5*(reyi/rey0)**0.25

len0   = (dx*dy*dz)**(1.0/3.0)
visk   = visd/rho  ! kinematic viscosity. unit of rho should be MKS
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
  ckxx(i,j)=ckxx(i,j)*max(array_ll(i,j),1.0)
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
subroutine collision(g, i, j)
use ieee_arithmetic
implicit none

integer, intent(in) :: i, j
real, dimension(nkr), intent(inout) :: g

integer :: k
real :: g1, m1, g2, m2, g3, g3n, g31, g32
real :: crn, beta, r_flux, r_rem, x01, dg1, dg2, dg, x1

k = ima(i,j)

g1 = g(i)
m1 = mass_hydro(i)

g2 = g(j)
m2 = mass_hydro(j)

if(inan.eq.1) then
  if(ieee_is_nan(g1).or.ieee_is_nan(g2)) then
    write(6,*) i, j, g1, g2
    call wrf_error_fatal("NaN is detected before collision.")
  end if
  if(.not.ieee_is_finite(ckxx(i,j))) then
    call wrf_error_fatal("Befoer collision, Inf")
  end if
  if(ieee_is_nan(ckxx(i,j))) then
    call wrf_error_fatal("nan of ckxx")
  end if
  if(ckxx(i,j).lt.0.0) then
    write(6,*) i, j
    write(6,*) ywxx(1,i,j),ywxx(2,i,j),ywxx(3,i,j)
    call wrf_error_fatal("ckxx is negative.")
  end if
  if(g1.lt.0.0.or.g2.lt.0.0) then
    write(6,*) "strange values are inputted in the collision."
    write(6,*) i, g1, j, g2
  end if
end if

x01 = ckxx(i,j)*g1*g2
x01 = min(x01,g1*m2)
if(k.ne.j) x01 = min(x01,g2*m1)
dg1 = x01/m2
dg2 = x01/m1

g(i) = g(i)-dg1
g(j) = g(j)-dg2

! start of redistribution (exponential flux method) by Bott (JAS 2000)
dg  = dg1+dg2
g31 = g(k)

if(dg.gt.0.0) then

  if(k.lt.nkr) then
    crn = cour(i,j)
    g32 = g(k+1)
    x1  = g32/(g31+dg)

    if(x1.gt.1.0e-3.and.x1.lt.1.0e3.and.abs(x1-1.0).gt.1.0e-2) then
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
    g(k+1) = g32 + dg*r_flux

  else

    r_flux = 0.0

  end if

  r_rem = 1.0-r_flux
  g(k) = g31 + dg*r_rem

end if
! end of redistribution

if(inan.eq.1) then
  if(ieee_is_nan(g(i)).or.ieee_is_nan(g(j)).or.ieee_is_nan(g(k))) then
    write(6,*) i, j, k
    write(6,*) g(i), g(j), g(k), g(k+1)
    write(6,*) g1, g2, ckxx(i,j)
    write(6,*) ywxx(1,i,j),ywxx(2,i,j),ywxx(3,i,j)
    write(6,*) dg1, dg2, dg
    write(6,*) x1, beta, crn, r_flux
    call wrf_error_fatal("NaN in collision")
  end if
end if

end subroutine collision

!==========================================================================
subroutine collision_poisson(g, i, j)
use ieee_arithmetic
implicit none

integer, intent(in) :: i, j
real, dimension(nkr), intent(inout) :: g

double precision, dimension(:), allocatable :: p_it, g_temp, m_temp

double precision, parameter :: ln_p_min = log(1.0d-5)
double precision, parameter :: norm_range = 4.6 ! equiv. to 1.0e-5
double precision, parameter :: p_crit = 100.0

double precision :: m1, m2, n1, n2, p, ln_p, ln_p_it, ln_p_it0, del_n2, crn
double precision :: dg, dg1, dg2, g31, g32
double precision :: x1, beta, r_flux, r_rem
double precision :: mass_init, mass_final
integer :: it, k, k1, it_min, it_max, it_max2
logical :: flag

! note: 2 is larger than 1 !
m1 = mass_hydro(i)
m2 = mass_hydro(j)
n1 = g(i)/m1*dlnr
n2 = g(j)/m2*dlnr

p = ckxx(i,j)/dlnr*n1
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

g(i) = n1*m1/dlnr
g(j) = n2*m2/dlnr

! start of redistribution (exponential flux method) by Bott (JAS 2000)
do it = it_min, it_max2

  k = nkr
  do k1=1, nkr-1
    if(m_temp(it).ge.mass_hydro(k1).and.m_temp(it).lt.mass_hydro(k1+1)) then
      k = k1
      exit
    end if
  end do

  dg  = g_temp(it)
  g31 = g(k)

  if(k.lt.nkr) then
    crn = log(m_temp(it)/mass_hydro(k))/log(mass_hydro(k+1)/mass_hydro(k))
    g32 = g(k+1)
    x1  = g32/(g31+dg)
    if(x1.gt.1.0e-3.and.x1.lt.1.0e3.and.abs(x1-1.0).gt.1.0e-2) then
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

    g(k+1) = g32 + dg*r_flux
  else
    r_flux = 0.0
  end if

  r_rem = 1.0-r_flux
  g(k) = g31 + dg*r_rem

end do

end subroutine collision_poisson

!==========================================================================
subroutine breakup_water(g1, dtime_coal)
implicit none

real, dimension(nkr), intent(inout) :: g1
real, intent(in) :: dtime_coal

real, dimension(nkr) :: g1_coll
real :: dropmass_before, dropmass_after, dropmass_delta, dtime_break
integer :: i, ndiv

dropmass_before = sum(g1)
g1_coll = g1
ndiv = 1

outer: do
  g1 = g1_coll
  dtime_break = dtime_coal/ndiv

  inner: do i=1, ndiv
    call coll_breakup(g1,dtime_break)

!! for implicit method
    if(ibreakup.eq.1) then
      dropmass_after = sum(g1)
      dropmass_delta=abs(dropmass_after-dropmass_before)
      if(dropmass_delta.gt.tol*dropmass_before) then
        ndiv = ndiv*2
        cycle outer
      end if
!! for explicit method
    else if(ibreakup.eq.2) then
      if(minval(g1).lt.0.0) then
        ndiv = ndiv*2
        cycle outer
      endif
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

ft = g1/(3.0*mass_hydro**2)

!breakup: Bleck (1970) first order method
!pkij: gain coefficients
!qkj : loss coefficients
do k=offset+1, nkr

  do i=offset+1, nkr
  do j=offset+1, i
    gain(k)=gain(k)+ft(i)*ft(j)*pkij(k-offset,i-offset,j-offset)
  enddo
  enddo

  do j=offset+1, nkr
    loss(k)=loss(k)+ft(j)*qkj(k-offset,j-offset)
  enddo

  if(k.lt.nkr) then
    weight(k)=2.0/(mass_hydro(k+1)**2-mass_hydro(k)**2)
  else
    weight(k)=2.0/((mass_hydro(k)*2)**2-mass_hydro(k)**2)
  end if

enddo

do k=offset+1, nkr
! analytic solution by Hyunho
!  if(loss(k).gt.0.0) then
!    beta(k)=exp(-weight(k)*loss(k)*dt)
!    ft(k)=gain(k)/loss(k)*(1.0-beta(k))+ft(k)*beta(k)
!    else
!      ft(k) = ft(k) + weight(k)*gain(k)*dt
!  end if
! implicit solution by Hyunho -> positive, stable
  if(ibreakup.eq.1) then
    ft(k) = (ft(k) + weight(k)*gain(k)*dt)/(1.0 + weight(k)*loss(k)*dt)
! explict solution in the original version
  else if(ibreakup.eq.2) then
    ft(k)=ft(k) + weight(k)*(gain(k) - loss(k)*ft(k))*dt
  else
    call wrf_error_fatal("something wrong in ibreakup")
  end if
end do

g1 = ft*(3.0*mass_hydro**2)

if(inan.eq.1) then
  do i=1, nkr
    if(ieee_is_nan(g1(i))) then
      write(6,*) "in breakup", i, gain, loss
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

! multiply by 1.0e2 to convert mks to cgs
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

if(aa.le.r_hydro(1))   vtbeard=bulkvt(1)
if(aa.gt.r_hydro(nkr)) vtbeard=bulkvt(nkr)

do kr=1,nkr-1
  if(aa.gt.r_hydro(kr).and.aa.le.r_hydro(kr+1)) then
    vtbeard=bulkvt(kr+1)
  endif
enddo   

end function vtbeard

!====================================================================
subroutine diffuse(tt, pp, qq, rhoa, dt, ff, m_new)
use ieee_arithmetic
implicit none

real, dimension(nkr), intent(in ) :: ff
real, dimension(nkr), intent(out) :: m_new
real, intent(in) :: pp, rhoa, dt
real, intent(inout) :: tt, qq

integer :: n
real, dimension(nkr) :: dmdt
real :: ffi, capi, vti, dmdti, sum_dmdt, m23, m13, delm, a, Sw, Si, Sw_dt, Si_dt, &
        rw, pw, tto, qqo, swo, sio

! supersaturation
Sw = pp*qq/(eps+qq)/esw(tt)-1.0
Si = pp*qq/(eps+qq)/esi(tt)-1.0

if(abs(Sw).lt.1.0e-4.or.qq.le.0.0) then
  m_new = mass_hydro
  return
end if

tto = tt
qqo = qq
swo = Sw
sio = Si

if(inan.eq.1) then
  if(Sw.gt.0.5) then
    write(6,*) "in diffusion"
    write(6,*) pp, qq, tt, (Sw+1.0)*100.0
    call wrf_message("too high rh.")
  end if
end if

do n=1, nkr
  capi =  r_hydro(n)
   vti = vt_hydro(n)

  call rate(pp, tt, rhoa, capi, vti, dmdti)

  dmdt(n) = dmdti
end do

sum_dmdt = dot_product(dmdt,ff)

if(qq.gt.1.0e-10) then
  a = eps/(qq*(eps+(1.0-eps)*qq))
  rw = (a+lv*lv/Cp/Rv/(tt**2))*(1.0+Sw)*sum_dmdt
  pw = (a+lv*ls/Cp/Rv/(tt**2))*(1.0+Si)*sum_dmdt
else
  rw = pp/(esw(tt)*eps)*sum_dmdt
  pw = pp/(esi(tt)*eps)*sum_dmdt
end if

call supsat(Sw, Si, Sw_dt, Si_dt, rw, 0.0, pw, 0.0, dt)

do n=1, nkr
  m23 = mass_hydro(n)**(2.0/3.0)
  m13 = mass_hydro(n)**(1.0/3.0)
  delm = (2.0/3.0)*dmdt(n)*Sw_dt/m13
  a = m23+delm
  if(a.ge.0.0) then
    m_new(n)=a**(3.0/2.0)
  else
    m_new(n)=0.0
  endif
end do

!call S_to_tq(Sw, Si, tt, pp, qq)

delm = dot_product(m_new-mass_hydro,ff)
tt   = tt + delm*Lv/Cp
qq   = qq - delm

! validity test
if(inan.eq.1) then
  Sw = pp*qq/(eps+qq)/esw(tt)-1.0
  Si = pp*qq/(eps+qq)/esi(tt)-1.0

  if(ieee_is_nan(tt).or.ieee_is_nan(qq)) then
    write(6,*) "in diffusion"
    write(6,*) tt, pp, qq, rhoa
    write(6,*) ff
    write(6,*) m_new
    write(6,*) dmdt
    write(6,*) m23, m13, delm, a
    call wrf_error_fatal("nan detected in diffuse")
  end if
  do n=1, nkr
    if(ieee_is_nan(m_new(n))) then
      write(6,*) "in diffusion"
      write(6,*) tt, pp, qq, rhoa
      write(6,*) ff
      write(6,*) m_new
      write(6,*) dmdt
      write(6,*) m23, m13, delm, a
      call wrf_error_fatal("nan detected in diffuse")
    end if
  end do
  if(Sw*swo.lt.0) then
    write(6,*) "in diffusion", Sw, swo, Sw_dt
    write(6,*) "sign of supersaturation was changed."
  end if
  if(qq.lt.0) then
    write(6,*) "in diffusion", Sw, swo, Sw_dt
    call wrf_error_fatal("negative vapor.")
  end if
  if(abs(tt-tto).gt.5.0) then
    write(6,*) "in diffusion"
    write(6,*) tto, qqo, tt, qq
    write(6,*) rhoa, pp
    write(6,*) pw, rw, sum_dmdt
    write(6,*) Sw, Si, swo, sio
    write(6,*) Sw_dt, Si_dt
!    write(6,*) ffi, dm
    call wrf_error_fatal("too high temperature in diffusion")
  end if
end if

end subroutine diffuse

!====================================================================
subroutine rate(pp, tt, rhoa, capi, vti, dmdt)
implicit none

real, parameter :: difv = 2.11e-5    ! diffusivity of vapor
real, parameter :: thca = 4.38e-3    ! thermal condunctivity of air
real, parameter :: visa = 1.827e-5   ! dynamic viscosity of air

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

if(rw.le.epsilon.and.ri.le.epsilon) then

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

  alpha = sqrt((rw-pi)**2+4.0*pw*ri)
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

!==============================================================================
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

!==============================================================================
subroutine remapping(rn,fi)
implicit none

real, dimension(nkr), intent(in) :: rn
real, dimension(nkr), intent(inout) :: fi

integer :: i, kr
integer :: isign_diffusional_growth

real, dimension(nkr) :: psi,rrs
real :: rnt,rrt,rrp,rrm,gn1,gn2,gn3

psi = 0.0

! rrs(i), g - default masses before condensation or evaporation
! rn(i),  g - new     masses after  condensation or evaporation

isign_diffusional_growth = 0
rrs = mass_hydro

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
            else ! use 3point method.
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

return
end subroutine remapping

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
subroutine vt_pressure(ppi)
implicit none

real, intent(in) :: ppi

integer :: n

do n=1, nkr

  vt_hydro(n) = bulkvt(n)*sqrt(p1000hpa/ppi)

end do

end subroutine vt_pressure

!----------------------------------------------------------------------------
subroutine check_validity(ff, subname)
use ieee_arithmetic
implicit none

real, dimension(nkr), intent(inout) :: ff
character(len=*), intent(in) :: subname

integer :: n

! NaN test

if(inan.eq.1) then
  do n=1, nkr
    if(ieee_is_nan(ff(n))) then
      write(6,*) n, trim(subname)
      call wrf_error_fatal("NaN detected in ff")
    end if
  end do
end if

do n=1, nkr
  if(ff(n).lt.0.0) then
    if(ff(n).lt.-0.0001) then
      write(6,*) n, ff(n), trim(subname)
      call wrf_error_fatal("negative ff detect")
    else
      ff(n) = 0.0
    end if
  end if
end do

end subroutine check_validity

!===================================================================!
end module module_mp_sbm_warm
