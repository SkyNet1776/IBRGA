      program ibrga
      common nsll, kpr, fracsl(10), dsdxsl(10), surfsl(10), nslp(10),
     1 tsl(10), pbrch, pbase, pmean, bbr(10), abr(10), deltat, y(20),
     1 igrad
      character cutfil * 10, bdfile * 10, style * 10
      character title(15) * 4, vsn * 4
      dimension br(10), trav(10), rp(10), tr(10), forcp(10), tempp(10),
     1 covp( 10)
      dimension chwp(10), rhop(10), gamap(10), nperfs(10), glenp(10),
     1 pdp(10 ), gdiap(10), alpha(10, 10), beta(10, 10), pres(10, 10)
      dimension a(4), b(4), ak(4), d(20), p(20), z(20), frac(10),
     1 surf(10 ), volp(10), dsdx(10), nbr(10), ibo(10), tbo(10),
     1 d2xdt2(10), tng(10)
      real lambda, jlzp, j2zp, j3zp, j4zp, jlzb, j2zb, j3zb, j4zb
      real 11,12,13
      dimension chdist(6), chdiam(6), bint(10), projtr(20), projms(20)
      dimension nsl(10), surfo(10), dsdxn(10)
      data pi/3.14159/vsn/'4hboat'/ c
c
c     USER'S MANUAL FOR IBRGA
c
C
c     IBRGA relies on an input database consisting of all
c     numerical parameters essential for running the code. Values may
c     be in metric units or in Imperial units, but must be consistent
c     throughout a dataset. Below is a compilation of a typical data
c     base showing the name and location of each parameter. The names
c     for the numerical values are prefixed with an alphabetical
c     designator corresponding to the position at which the data is to
c     appear, that is, from left to right. The data may be separated
c     by blanks or commas. Measurement units, if any, are shown to the
c     right of each input. In general, metric units of weight and mass
c     are the meter and kilogram, respectively; corresponding Imperial
c     units are the inch and pound. The only exceptions are Imperial
c     units of propellant impetus, which are foot-pounds per pound mass.
c
c
c     title card - up to 60 characters of title and identification
c
c     parameter information and placement:
c
c     A B C D E F G H I J K
c
c     record 1 Metric Imperial
c     A. - chamber volume (when canister (m**3) (in**3)
c     model is used, this will be the
c     volume after canister bursts)
c     B. - groove diameter (i) (in)
c     C. - land diameter (m) (in)
c     D. - groove/land ratio (land, groove,
c     and groove/land ratio used to
c     calculate the tube bore area)
c     E. - twist (units are turns/caliber)
c     F. - projectile travel (m) (in)
c     G. - gradient switch (integer value
c     designating the gradient equation
c     (1 = Lagrange, 2 = Chambrage,
c     pg 57
c     3 = Two-phase, 4 - RGA,
c     5 = Lagrange w/ bt, 6 = Cham. w/ bt)
c     H. - variable projectile mass switch
c     (0=no, 1=yes)
c     I. - igniter canister model switch
c     (0=no, l=yes)
c     J. - friction factor (normally 1 for
c     granular, 0.01 for stick and
c     0.1 for partially cut propellant;
c     only used when gradient = 3 or 4)
c
c     record la (read if and only if gradient - 5 or 6)
c     A. - boattail diameter (M) (in)
c     B. - boattail length (i) (in)
c
c     record lb (Read if and only if gradient - 2 or 4 or 6)
c     A. - number of point pairs to describe
c     chamber geometry, integer I <- 5
c     B. - initial distance from breech (i) (in)
c     (must be 0.0)
c     C. - diameter at initial distance (M) (in)
c
c
c
c -   Ith distance from breech (i) (in)
c     (initial position of the base
c     of the projectile)
c     - Ith diameter at Ith distance (W) (in)
c     (used to calculate bore area -
c     overrides record 1 groove and
c     land diameter specifications)
c     (Note: chamber geometry is used
c     to calculate the chamber volume
c     which overrides record 1 chamber
c     volume description.)
c
c
c     record 2
c     A. - projectile mass (kg) (lb)
c     B. - switch to calculate energy lost
c     to air resistance, an integer
c     either 0 = no loss, or 1 = loss
c     C. - fraction of bore resistance work
c     used to heat tube (0.0<=f<=1.0)
c     D. - gas pressure ahead of projectile (MPa) (psi)
c
c     record 2A (Read if and only if variable projectile mass
c     switch is 1)
c     A. - number of point pairs to describe
c     variable projectile mass =< 20
c     B. - initial projectile travel (i) (in)
c     (conceptually should be 0.0)
c     C. - initial projectile mass (kg) (lb)
c     (overrides value from record 2)
c     D. - projectile travel at which first (i) (in)
c     mass change occurs
c     E. - new projectile mass at travel D. (kg) (ib)
c     .
c     Page 58
c
c     x - i-th projectile travel where mass (i) (in)
c     change occurs
c     y - i-th new projectile mass value (kg) (ib)
c
c     record 2B (Read if and only if i;niter
c     canister model switcn = 1)
c     A. - pressure at which the igniter (MPa) (psi)
c     canister will burst
c     B. - volume of igniter canister (m**3) (in**3)
c     (used as chamber volume until
c     burst pressure achieved)
c     C. - canister diameter (assumes a (W) (in)
c     right circular cylinder)
c
c     record 3
c     A. - number of pairs of barrel
c     resistance points (integer <= 10)
c     B. - bore resistance (MPa) (psi)
c     C. - travel (i) (in)
c
c
c
c     - Jth bore resistance (MPa) (psi)
c     - Jth travel (i) (in)
c
c     record 4
c     A. - mass of recoiling parts (kg) (lb)
c     B. - number of recoil point pairs
c     (must be an integer = 2)
c     C. - recoil force (force to overcome (N) (lb)
c     before recoil start - rod preload)
c     D. - time of rod preload (must be 0.0) (s) (s)
c     E. - recoil force (constant resistive (N) (lb)
c     force after rise time)
c     F. - rise time (time to go from recoil (s) (s)
c     start to constant resistive
c     recoil force)
c
c     record 5
c     A. - free convective heat transfer (W/m**2/K) (in-lb/in**2
c     coefficient /s/K)
c     B. - chamber wall thickness (wall (i) (in)
c     depth which is heated uniformly)
c     C. - heat capacity of chamber wall (J/kg/K) (in-lb/lb/K)
c     D. - initial temperature of tube and (K) (K)
c     chamber walls
c     E. - heat loss coefficient (usually 1.,
c     but may be set to 0.0 in order to
c     eliminate heat loss)
c     F. - density of chamber wal' (kg/m**3) (lb/in**3)
c
c     record 6
c     A. - impetus of igniter (J/kg) (ft-lb/lb)
c     B. - adiabatic flame temperature of (K) (K)
c     igniter material
c     C. - covolume of igniter (m**3/kg) (in**3/lb)
c     D. - ratio of specific heats of igniter
c     Page 59
C     E. - mass of igniter (kg) (lb)
c
c     record 7
c     A. - number of propellants
c     (integer <= 10)
c
c     record 8
c     A. - impetus of propellant (J/kg) (ft-lb/ilb)
C     B. - adiabatic flame temperature (K) (K)
c     C. - covolume of propellant (m**3/kg) (in**3/lb)
c     D. - ratio of specific heats
c     E. - mass of propellant kg) (lb)
c     F. - density of propellant (kg/m**3) (lb/in**3)
c     G. - propellant form function indicator
c     (integer; may be one of:
c     0 solid cylindrical grain
c     I single-perf cylindrical grain
c     2 spherical grain
c     7 seven-perf cylindrical grain
c     15 nineteen-perf hexagonal grain
c     19 nineteen-perf cylindrical grain)
c     H. - length of propellant grain (M) (in)
c     I. - diameter of perforations in the (W) (in)
c     propellant grains (ignored if not
c     required, but must be present)
c     J. - outside diameter of propellant (i) (in)
c     grain (for the hexagonal grain
c     it is the distance between
c     rounded corners)
c
c     (Record 8 repeated for each propellant)
c
c     record 9
c     A. - number of burning rate triplet points
c     (integer J <= 10)
c     B. - exponent
c     C. - coefficient (m/s-MPa**e) (in/s-psi**e)
c     D. - pressure (upper pressure limit (MPa) (psi)
c     for which the previous exponent
c     and coefficient are valid) c
c
c
c     -Jth exponent
c     -Jth coefficient (m/s-MPa**e) (in/s-psi**e)
c     Jth pressure (if pressure should (MPa) (psi)
c     exceed this limit, then this
c     burning rate equation is used
c     for all higher pressures)
c
c     (Record 9 repeated for each propellant)
c
c     record 10
c     A. - integration time increment (ms) (ms)
c     B. - print increment (Ms) (ms)
c     C. - upper limit on integration time (ms) (ms)
c     for stopping calculation
c*
c     Page 60
c*    conversion factors for imperial units => metric units
c*
c     length inches * .0254 => meters
c     mass lb * .45359237 => kilograms
c     area in^2 * .00064516 => m^2
c     volume in^3 * 0.000016387064 => m^3
c     pressure lb/in^2 * 6894.757 => pascals
c     velocity ft/s / 3.28083 => m/s
c     in/s * .0254 => m/s
c     energy ft-lb * 1.3558179 => joules
c     in-lb * 0.1129848 => joules
c     density lb/in^3 * 27679.9 => kg/m^3
c     force/mass (ft-lb)/ib * 2.989067 => j/kg
c     covolume in^A3/lb / 27679.9 => m^3/kg
c
c     source: engineering design handbook metric conversion guide
c     darcom pamphlet 706-470, july 1976
c
c     density grams/cc * 1000. => kg/m^3
c
c     call gettim(ihr,imin,isec,ihuns)
c
      write( *, 830)
      read( * 840)bdfile
      open(unit = 2, err = 810, file = bdfile, status = 'old', iostat =
     1 ios)
      nzp=0
      rewind 2
      write( *, 850)
      read( *, 840)outfil
      open(unit = 6, err = 820, file = outfil)
      do 10 i = 1, 20
      p(20) = 0.
      y(20) = 0.
      z(20) = 0.
      d(20) = 0.
      10 continue
      write( *, 870)
      read( *, 840 )style
      mode = 0
      if(style(1:1).eq.'m' .or. style(1:1).eq.'M') mode = 1
      if(style(1:1).eq.'e' .or. style(1:1).eqo'E') mode = 2
      if(mode.eq.0) write( *, 880)
      if(mode.eq.0) stop
      read(2, 885)title
      write(6, 1236)title,vsn
      write(6, 860)bdfile
      read(2, *, end = 790, err = 800)cham, grve, aland, glr, twst,
     1 travp, igrad, ivpm, ihl, fs0
      if(ijrad.gt.l)go to 20
      write(6, 890)
      igrad = 1
      go to 140
c
c     define chambrage assumes nchpts=number of points to define
c     chamber > or = 2 < or = 5 (?),chdiam(i) defines chamber diameter
c     at chdist (i) chamber distance. chdiam(nchpts) is assumed to be
c     the bore diameter and chdist(i) is assumed to be 0, i.e. at the
c     Page 61
c     breech. assumes truncated cones.
c
20    if(igrad.eq.3)go to 130
      if(igrad.eq.4)go to 30
      if(iqrad.GE.5)go to 25
      write(6, 900, err 800)
      go to 40
25    read (2 , ,end =790,err - 800) btdia,btlen
      if (mode.eq.1) then
      btvol=pi*btdia*btdia/4 . *btlj
      write (6, 955)btdia,btlen,btvol
955   format(/,lx,'b~attail diameter m ',el6.6/lx,'boattail length
     &m ',el6.6/lx,'boattail volume m**3 ',e16.6,/)
      else
      btvol=pi*btdia*btdia/4. *btlen
      write (6, 965,err=800)btdia,btlen,btvol
965   format(/,lx,'boattail diameter in ',el6.6/1x,'boattail length
     &in ',e16.6/lx,'boattaiJ. volume in**3 ',e16.6,/)
      btdia=btdia*0 .0254
      btlen=btlen*0 .0254
      btvol=btvol*0 .0254"'*3
      cham=cham* . 6387064e-5
      endif
35    if(igrad.eq.5)then
      write (6, 975)
975   format(lx,'using lagrange with boattail gradient')
      go to 140
      end if
      if(igrad.eq.6)then
      write (6, 995)
995   format(lx,'using chambraqe with boattail gradient')
      go to 40
      endif
30    write(6, 910)
      go to 40
40    read(2, *, end - 790, err - 800)nchpts, (chdist(i), chdiam(i), i
      1 = 1, nchpts)
      if (mode.eq. l)then
      write(6, 920, err = 800) (chdist(i), cridiam(i), i - 1, nchpts)
      goto 60
      else
      write(6, 925, err - 800) (chdist(i), chdiam(i), i - 1, nchpts)
      do 50 i = 1, nchpts
      chdist(i) = chdist(i) * 0.0254
      chdiam(i) = chdiam(i) * 0.0254
50    continue
      endif
c
c     calculate chamber integrals and volume
c
60    if(nchpts.gt.5) write(6, 930, err - 800)
      if(nchpts.gt.5)nchpts - 5
      bore = chdiam(nchpts)
      if(chdist(l).ne.0.0)write(6, 940, err = 800)
      chdist(l) =0.0
c     do 54 I=i,nchpts
c     chdist(I)=0.01*chdist(I)
c     Page 62
c54   chdiam(I)=0.0l*chdiam(I)
c     calculate chamber integrals and volume
      if(nchpts.gt.5) write(6,44,err=30)
44    format(lx,'use first 5 points')
      if (nchpts .gt .5) nchpts=5
      bore=chdiam (nchpts)
      if(chdist (1).ne.0.O)write(6,45,err-3O)
45    format(lx,' # points ? '
      chdist (1)=0. .0
      ptl=chdist (nchpts)
      btd=btdia
      btl=btlen
      call jint(btd,btl,ptl,ptl,nchptS,chdist,chdiam,bint,bvol)
41    cham=bvol+btvol
c     write(6,47,err=30)bint(l),bint(3),bint(4 )
c     format(lx,'bint 1 = ',e14.6,' bint 3 = ',e14.6,' bint. 4 =',el4.
c    &6)
      chmlen=chdist (nchpts)
      go to 140
130   write(6, 950)
140   if(mode.eq.l)then
      write(6, 960, err -800)charn, grve, aland, glr,
     & twst, travp, igrad, ivpm, ihl, fs0
      cham=cham-btvol
      endif
      if (mode .eq. 2) then
      cham--cham/ . 6387064e-5
      write(6, 970, err = 800)cham, grve, aland, glr, twst, travp,
     1 igrad, ivprn, ihl, fs0
      chain = chain * 1.6387064e - 5
      cham=cham-btvol
      grve =grve * 0.0254
      aland =aland * 0.0254
      travp =travp * 0.0254
      endif
      read(2, *, end = 790, err = 800)prwt0, iair, htfr, pgas0
      if (iode.eq.l)then
      prwt = prwtO
      pgas = pgas0 * 1.0e6
      elseif (mode.eq.2)then
      prwt = prwt0 * 0.45359237
      pgas = pgas0 * 6894.757
      endif
      if (ivpm.eq.1) then
      read(2, * )nvpmp, (projtr(i), projms(i), i = 1, nvpmp)
      if(mode.eq.1)write(6, 980)nvpmp, (projtr(i),
     1 projms(i), i = 1, nvpmp)
      if (mode .eq. 2) then
      write(6, 985)nvpmp, (projtr(i), projms(i), i =1, nvpmp)
      do 150 i = 1, nvpmp
      projtr(i) = projtr(i) * 0.0254
      projms(i) = prcjms(i) * 0.45359237
150   continue
      endif
      prwt =projms(1)
      prwt0 prwt
      if(mode.eq.2) prwt0 = prwt / 0.45359237
      endif
c     Page 63
      write(6, 1050)
      if(ihl.eq.1)then
      read(2, * )burstp, highv, highd
      if(mode.eq.1)write(6, 990)burstp, highv, highd
      if (mode .eq. 2) then
      write(6, 1000)burstp, highv, highd
      burstp =burstp * 0.006894757
      highv highv * 1.6387064e - 5
      highd =highd * 0.0254
      endif
      burstp =burstp * 1.e6
      areaw =4. *highv / highd + pi * highd *highd /2.
      endif
      read(2, *~end =790, err = 800)npts, (br(i), trav(i), i = ,npts)
      read(2, *~end =790, err = 800)rcwt, nrp, (rp(i), tr(i), i=1,nrp)
      read(2, *,end =790, err = 800)ho, tshl, cshl, twal, hi, rhocs
      read(2, *~end =790, err = 800)forciq, tempi, covi, gamai, chwi
      read(2, *)nprop
      read(2, *,end = 790, err = 800) (forcp(i), tempp(i), covp(i), 1 gamap(i), chwp(i), rhop(i), nperfs(i), glenp(i), pdp(i), gdiap(i)
     1 , i= 1, nprop)
      if (mode .eq.1) then
      write(6, 1010, err = 800)prwt0, iair, htfr, pgas0
      write(6, 1030, err = 800)npts, (br(i), trav(i), i = 1, npts) write(6, 1060, err = 800)rcwt, nrp, (rp(i), tr(i), i = 1, nrp) write(6, 1080, err = 800)ho, tshl, cshl, twal, hi, rhocs
            write(6, 1100, err - 800)forcig, tempi, covi, gamai, chwi
            write(6, 1236)title,vsn
            write(6, 1120)nprop
            write(6, 1130, err = 800) (i, forcp(i), tempp(i), covp(i),
     1 gamap(i), chwp(i), rhop(i), nperfs(i), glenp(i), pdp(i),
     1 gdiap(i), i = 1, nprop)
      endif
      if (mode .eq. 2) then
            write(6, 1020, err = 800)prwt0, iair, htfr, pgasO
            write(6, 1040, err = 800)npts, (br(i), trav(i), i = 1, npts)
            write(6, 1070, err = 800)rcwt, nrp, (rp(i), tr(i), i - 1, nrp)
            write(6, 1090, err = 800)ho, tshi, cshl, twal, hi, rhocs
            write(6, 1110, err = 800)forcig, tempi, ccvi, gamai, chwi
            write(6, 1236)title,vsn
            write(6, 1120)nprop
            write(6, 1140, err = 800) (i, forcp(i), tempp(i), covp(i),
     1 gamap(i), chwp(i), rhop(i), nperfs(i), glenp(i), pdp(i),
     1 gdiap(i), i = 1, nprop)
      endif
      do 170 j = 1, nprop
      read(2, *,end= 790, err = 800)nbr(j), (alpha(j, i), beta(j,i)
     1 , pres(j, i), i = 1, nbr(j))
      if (mode.eq.1)write(6, 1160) nbr(j)
      if (mode .eq. 2) write (6, 1170) nbr (j)
      do 160 i = 1, nbr(j)
      if(mode.eq.1)write(6, 1180) aipha(j, 1), beta(j, i),
     1 pres(j, i)
      if (mode .eq. 2) then
      write(6, 1180) alpha(j, i), beta(j, i), pres(j, i)
      rate = beta(j, i) * pres(j, i) ** alpha(j, i)
      pres(j, i) = pres(j, i) * 0.006894757
      beta(j, i) = 0.0254 * rate / pres(j, i) ** alpha(j,i)
c     Page 64
      endif
160   continue
170   continue
c
c     convert units to program requirements
c
      do 180 i = 1, npts
      if(mode.eq.l)br(i) = br(i) * 1.e6
      if (mode.eq.2) then
      br(i) = br(i) * 6894.757
      trav(i) = trav(i) * 0.0254
      endif
180   continue
      do 200 j = 1, nprop
      if (mode.eq.2)then
      forcp(j) = forcp(j) * 2.989067
      covp(j) = covp(j) / 27679.9
      chwp(j) = chwp(j) * 0.45359237
      rhop(j) = rhop(j) * 27679.9
      glenp.(j) = glenp(j) * 0.0254
      pdp(j) - pdp(j) * 0.0254
      gdiap(j) = gdiap(j) * 0.0254
      endif
      do 190 i = 1, nbr(j)
      pres(j, i) - pres(j, i) * 1.e6
190   continue
200   continue
      if (mode.eq.2) then
      do 210 i = 1, nrp
      rp(i) = rp(i) * 0.1129848
      tr(i) = tr(i) * 0.0254
210   continue
c
c     conversion factor for free convective heat transfer coeff
c     w/m**2 * (0.00064516 m**2/in**2) * (1.0/1.3558179 ft-lb-s/w) *
c     (12.0 in/ft) = 0.005710147 in-lb/in**2-s
c
      ho - ho / 0.005710147
      tshl = tshl * 0.0254
      cshl = cshl * 2.989067 / 12.0
      rhocs = rhocs * 27679.9
      icwt = rcwt * 0.45359237
      forcig = forcig * 2.989067
      covi = covi / 27679.9
      chwi = chwi * 0.45359237
      endif
      tmpi = 0.0
      do 220 i = 1, nprop
      tmpi = tmpi + chwp(i)
      kpr = i
      call prf7l0(pdp(i), gdiap(i), glenp(i), nperfs(i), 0., frac(i)
     1 , volp(i), surf(i), dsdx(i))
      tng(i) = chwp(i) / rhop(i) / volp(i)
      surfo(i) = surf(i)
      write(6, 1150)i, tng(i)
220   continue
      tmpi = tmpi + chwi
      write (6, 1050)
c     Page 65
      read(2, *, end = 790, err = 800)deltat, deltap, tstop,nzpi
      write(6, 1190, err = 800)deltat, deltap, tstop
      write( *~1200)
      deltat =deitat * 0.001
      deltap =deltap * 0.001
      tstop = stop * .001
      if(igrad.eq.2.or.igrad.eq.4 .or.igrad.eq.6)go to 230
      bore = (gIr * grve *grve + aland * aland) / (gir + 1.)
      bore = sqrt (bore)
230   areab = pi * bore *bore I4.
      areaa=pi* (btdia/2.) **2
      areaba-areab-areaa
      lambda =1. / ((13.2 + 4. *log10(100. *bore)) **2)
      iplot =0
      pltdt =deltat
      Pitt =0.
      pmaxm =0.0
      pmaxbr = 0.0
      pmaxba = 0.0
      tpmaxm = 0.0
      tpmxbr = 0.0
      tpmxba = 0.0
      tpmax 0.0
      a(l) =0.5
      a(2) =1. - sqrt(2.) / 2.
      a(3) =1. + sqrt(2.) / 2.
      a(4) =1. / 6.
      b(l) =2.
      b(2) =1.
      b(3) =1.
      b(4) =2.
      ak(l) =0.5
      ak(2) =a(2)
      ak(3) =a(3)
      ak(4) 0.5
      vp0 = 0.0
      trO = 0.0
      tcw = 0.0
      if(igrad.eq.3)chmlen =chain areab
      if (igrad. eq. 5) chmlen= (cham+btvol) /areab
      zb = chmlen
      zp = chmlen
      grien = 0.
      grdiam = 0.
      egama = 0.
      do 240 i =1, nprop
      grien =grien + chwp(i) * glenp(i)
      grdiam = grdiam + chwp(i) * gdiap(i)
      ibo(i) = 0
      egama = egama + chwp(i) * gamap(i)
      nsl(i) = 0
      vp0 = chwp(i) / rhop(i) + vp0
240   continue
      volgi = chain vp0 - chwi * covi
      grien = grien /(tinpi - chwi)
      grdiam = grdiam / (tinpi - chwi)
      egamna = (egama + chwi *gamai) / tinpi
      ismn = 0
c     Page 66
      odlnr = 0.
      vfO 0 chain - vp0
      epsO = 1. - vp0 / chain
      eps = eps0
      gasden = chwi /vf 0
      prden = tinpi IvpO
      ug =0.
      up =0.
      pinean - forcig * chwi / volgi
      if(ihl.eq.1)pinean =forcig * chwi/
     1 (highv - vp0 chwi *covi)
      pbase - pmean
      pbrch = pmean
      opbase = pmean
      voig = volgi
      volgi = volgi + vp0
      wailt = twal
      tgas = temnpi
      told = 0.
      tgaso = tgas
      dtgaso = 0.
      covi = covi
      t = 0.
      ptime = 0.0
      ibrp = 12
      z(3) = 1.
      nde = ibrp + nprop
      if(mode.eq.1)write(6, 1210)areab, pinean, vpO, volgi
      if (mode .eq. 2) then
      argl = areab /0.00064516
      arg2 = pinean /6894.757
      arg3 = vp0 / 0.000016387064
      arg4 = volgi / 0.000016387064
      write(6, 1220)argl, arg2, arg3, arg4
      endif
      write(6, 1236)title,vsn
      write(6, 1230)
      if(mode.eq.1)write(6, 1232)
      if(inode.eq.2)write(6, 1234)
      lines = 4
      linmax = 62
      iswl = 0
      prwt0 = prwt
250   continue
      do 690 j = 1, 4
c     find barrel resistance
c
      if(ivpin.ne.1)go to 270
      do 260 k = 2, nvpmp
      if(y(2) + y(7).lt.projtr(k))go to 270
      prwt = projms(k)
260   continue
270   if(ihl.eq.1)go to 300
      do 280 k = 2, npts
      if(y(2) + y(7).ge.trav(k))go to 280
      go to 290
280   continue
c     67
      k = npts
290   resp = (trav(k) - y(2) - y(7)) / (trav(k) - trav(k-1)
      resp = br(k) - resp * (br(k) - br(k - 1))
c
c     find mass fraction burned
c
300   do 320 k - 1, nprop
      kpr = k
      if(ibo(k) .eq.1)goto320
      nsll - 0
      call prf7lO(pdp(k), qdiap(k), glenp(k), nperfs(k),
1     y(ibrp + k), frac(k), volp(k), surf(k), dsdx(k))
      nsl(k) = nsll
      if(nsl(k) .eq.0)goto 310
      if(nslp(k).eq.1)go to 310
      write(6, 1240)k
      lines = lines + 1
      nslp(k) 1
      tsl(k) =y( 3)
      ism = 1
310   continue
      if(frac(k).lt..9999) go to 320
      frac(k) =1.
      tbo(k) =y(3)
      ibo(k) =1
      ism = 1
      write(6, 1250)k
      lines = lines + 1
320   continue
      if(ihl.eq.1)goto370
c
c     energy loss to projectile translation
c
      elpt = y(11)
c
c     elpt-prwt*y(1)*y(1)/2.
c
      eptdot =prwt * y (1) * z (1)
      z(11) =eptdot
c
c     energy loss due to projectile rotation
c
      elpr = y(12)
c
c     elpr=pi*pi*prwt*((y(l)+y(6))**2)/4.*twst*twst
c
      eprdot = pi * pi * prwt * (y(l) + y(6)) *(z(1) + z(6))
     1 /2. * twst * twst
      z(12) = eprdot
c
c     energy loss due to gas and propellant motion
c
      if(igrad.eq.l)go to 340
      if(igrad.eq.3)go to 350
      if(igrad.eq.4)go to 330
      if(igrad.eq.5)go to 352
      if(igrad.eq.6)go to 355
      pt = y(2) + y(7)
c     Page 68
      vzp =bvol + areab * pt
      j4zp bint(4) + ((bvol + areab * pt) ** 3 -bvol **3) /3.
     1 / areab / areab
      elgpm tmpi * y(l) * y(1) * areab *areab *j4zp I2. /vzp
     1 /vzp/vzp
      go to 360
330   pb =y(7) + y(lO)
      vzb bvol + areab * pb
      j4zb =bint(4) + (vzb ** 3 - bvol **3) / 3. / areab / areab
      elgpm =(1. -eps) * up *up * areab ** 2 * prden *j4zb + eps
     1 * ug *ug *areab **2 *gasden * j4zb
      elgpm =elgpm / 2. Ivzb /vzb + gasden *areab * ullen /6.
     1 * (3. y(l) * y(l) + 3. *y(l) * ullen *dlnrho + ullen **2
     1 * dlnrho ** 2)
c
c     approximate epdot
c
      epdot = tmpi * y(l) * z(1) / 3.
      go to 360
340   elgpm = tmpi * (y(J.) * y(l) - y(l) * y(6) + y(6) *y(6)) /6.
      go to 360
350   elgpm =areab *zb / 6. *(eps *gasden *ug *ug + (.-eps)
     1 * prden * up *up)
      elgpm = elgpm + gasden *areab *ullen /6. *(3. *y(1)*
     1 yMl + 3. * yMl * ulleri * dlnrho + ullen **2 *dlnrho **2)
c
c     approximate epdot
c
      epdot = tmpi * y (1) * z (1) / 3.
      go to 360
352   zp=y(2)+y(7)+chmlen
      za=zp-btlen
      vzp=zp*areab-btvol
      vza=za*areab
c     write(6,2050)tmpi,vzpj,areaa,y(1)
c     elgpm=tmpi*y (1)*y (1)/vzpl2.
c     elgpm=elgpm* (areab*areab/3. /vzp/vzp* (areab*za**3
c    & +((areab*za+areaba*(zp-za))**3-(areab*za)**3)/areaba/areaba)
c    & -areaa*areab/vzp/areaba/areaba* ((areab*za+areaba* (zp-za))
c    & **2-(areab*za)**2)+areaa*areaa/areaba*(zp-za))
      elgpml=tmpi*areab**2*y (1)**2/2 ./vzp/vzplvzp
      taql= (areab*za**3/3 .+
     & (areab*za+(zp-za)*areaba) **3/3./areaba/areaba.-(areab*za)**3/
     & 3./areaba/areaba)
      elgpm2=tmpi*areab*areaa*y (1)**2/vzp/vzp
      taq2= ((areab*za+ (zp-za)
     & *areaba) **2/2./areaba/areaba-~(areab*za) **2/2./areaba/areaba)
      elgpm3=tmpi*y (1)**2*areaa**2/2. /vzp
      taq3= (zp-za) /areaba
      elgpm=elgpml *taql-elgpm2*taq2+elgpm3*taq3
c     write(6,2050)tmpi,vzp,areaa,y(1)
c 2050 format(' tmpi',e17.10,' vzp',e17.i0,' areaa',e17.10,' yMl)
c & ,el7.10)
      go to 360
c     Page 69
355   pt-=y(2)+y(7)
      ptl-pt+chmnleri
      call jint(btdia,btlen,ptl,ptl,nchpts,chdist,chdiam,bflt~bvolzp)
      vzp=bvolzp
      qlzp-bint (10)
      q2zp=biit (2)
      q6zp=bint (4)
      q8zp~bint (8)
      q3zp=bint (5)
      q4zp=bint (6)
      q5zp-bint (7)
      q9zp=bint (9)
c     write(6 ,7 6 )qlzp,q2zp,q3zp,q4zp,qszp,q6zp~q7zpq8zp~ptlpt2 76 format(lx,l0ell.4)
      delta=l.
      ptl-chmlen+y (2) +y (7)
      pt2=chmlen+y (2) +y (7) -btlen
      zp~=ptl
      z a =Pt 2
      calljint(btdia,btlen,ptl,pt2,nchpts~chdist~chdiam~bint~bvolza)
      vza = bvolza
      qlza=biit (1)
      q2za=bint (2)
      q4za-bint (6)
      q5za=bint (7)
      q6za=bint (4)
      q7za=bint (3)
      q3zpza=q3zp
      q9zpza=ql zp
c     write (6, 7 7 )qlza~q2 za,q3za,q4zaiqsza,q6za,q7za~q8za~ptl~pt2
77    -format(Ix,10ell.4)
c     elgpm-tmpi*y(1)*y(1)*areab**2*q6zp/2./vzp**3 - tmpi*areaa*areab* "o 
     &y(l)*y(1)*q9zpza/vzp/vzp + tmpi*y(1)*y(1)*areaa**2*q3zpza/ "C &2./vzp
      elgpml=tmpi*y(l)*y(1) *areab**2*q6zp/2./vzp**3
      elgpm2=tmpi*areaa*areab*y (1)*y(l) *qlzp/vzp/vzp elgpm3=tmpi*y (1)*y (1)*areaa**2*q3zp12 . vzp elgpm=e lgpml-elgpm2+elgpm3
      go to 360
c
c     energy loss from bore resistance
c
360   elbr = y(4)
      z(4) = areab * resp * (y(l) +y()
      ebrdot = z(4) c
c     energy loss due to recoil
c
      elrc =rcwt *y(6) *y(6) / 2.
      erdot =rcwt * y(6) z z(6)
c
c     energy loss due to heat loss
c
      areaw -cham /areab *pi *bore +2. areab +pi bore*
     1 (y(2) +y()
370   avden = 0.0
      avc = 0.0
c     Page 70
      avcp =0. 0
      Z18 =0
      Z19 =0
      do 380 k =1, nprop
      z18 = forcp(k) * gamap(k) * chwp(k) *frac(k) / (gaxnap(k)
     1 - 1.) / tempp(k) + z18
      Z19 = chwp(k) * frac(k) + z19
      avden = avden + chwp(k) *frac(k)
380   continue
      avcp = (z18 + forcig * gamai *chwi / (gamai - 1.) /temapi)/
     1 (Z19 + chwi)
      avden = (avden + chwi) / (voig + covi)
      avvel = .5 * (y(l) + y(6))
      htns = lambda *avcp * avden *avvel + ho
      z(5) = areaw *htns * (tgas -walit) * hi
      elht = y(5)
      ehdot =z(5)
      wailt = (elht + htfr * elbr) Icshl / rhocs /areaw /tshl +
     1 twal
c
c     energy loss due to air resistance in tube
c     (assume no drag resistance on air/tube interface)
c
      if(ihl.eq.1)goto 410
      air =iair
      z(8) y y(1) *pgas * air
      elar =areab *Y(8)
      eddot z z(8) *areab
c
c     recoil
c
      z(6) =0.0
      if(pbrch.le.rp(l) / areab)go to 400
      rfor = rp(2)
      if(y(3) - tr0.ge.tr(2))go to 390
      rifor = (tr (2) -(y (3) - trO)) (tr (2) -tr (1))
      rfor = rp(2) -rfor * (rp(2) -rp(l))
390   z(6) = areab /rcwt * (pbrch -rfor /areab - resp)
      if(y(6) .lt.0.0)y(6) = 0.0
      Z(7) = y(6)
      goto 410
400   trO = y( 3 )
410   continue
c
c     calculate gas temperature
c
      eprop = 0.0
      rprop = 0.0
      dmfogt = 0.0
      dmfog = 0.0
      do 420 k =1, nprop
      eprop =eprop + forcp(k) * chwp(k) * frac(k) /(gamap(k)
      rprop =rprop + forcp(k) * chwp(k) * frac(k) /(gamap(k)
     1 - 1.) /tempp(k)
      dmfogt =dmfogt + forcp(k) *rhop(k) * tng(k) *surf(k)*
     1 z(ibrp + k) / ((gamap(k) -1.) * tempp(k))
      drnfog = drnfog + forcp(k) * rhop(k) * tng(k) * surf(k)*
c     Page 71
     1 z(ibrp + k) / (gamap(k)-1.
      420 continue
      tenerg = elpt + elpr + elgpm + elbr + elrc + elht + elar
      tgas = (eprop + forcig * chwi / (gamai - 1.) - elpt - elpr -
     1 elgpm - eibr - elrc - elht - elar) / (rprop + forcig * chwi
     1 / (gamai - 1.) / tempi)
      tedot = epdot + eprdot + eddot + ebrdot + erdot + ehdot+eptdot
      dtgas = (dxnfog - tedot - tgas * dmfogt) / (rprop + forcig*
     1 chwi / (gamai - 1.) / tempi)
c
c     find free volume
c
      V1 = 0.0
      covi = 0.0
      do 430 k - 1, nprop
      vi = chwp(k) * (1. - frac(k)) / rhop(k) + v1
      covi covi + chwp(k) * covp(k) *frac(k)
430   continue
      voig = volgi + areab * (y(2) + y(7)) -v1 - covl
      if(ihl.eq.i)volg = highv - vi - covi
c
c     calculate mean pressure
c
      ri 0.0
      do 440 k = 1, nprop
      ri = ri + forcp(k) * chwp(k) * frac(k) /tempp(k)
440   continue
      pmean = tgas / voig * (ri + forcig *chwi /tempi)
      if(ihl.eq.i)go to 640
      resp = resp + pgas * air
      if(igrad.eq.i)go to 590
      if(igrad.eq.2)go to 450
      if(igrad.eq.3)go to 470
      if(igrad.eq.4)go to 540
      if(igrad.eq.5)go to 582
      if(igrad.eq.6)go to 585
450   if(iswi.ne.0)go to 460
      pbase = pmean
      pbrch = pmean
      if(pbase.gt.resp + 1.)iswi 1
      go to 620
c
c     use chambrage pressure gradient equation
c
460   jizp = bint(i) + (bvoi * pt + areab /2. *pt *pt) /areab
      j2zp = (bvol + areab * pt) ** 2 /areab /areab
      j3zp = bint(3) + areab *bint(i) *pt + bvol * pt * pt /2. +
     1 areab pt *pt * pt /6.
      a2t = -tmpi *areab *areab / prwt / vzp /vzp
      alf = 1. - a2t * jizp
      alt = tmpi * areab * (areab * y (i) y y(1) /vzp + areab *resp
     1 / prwt) / vzp / vzp
      bt = - tmpi * y(i) * y(i) *areab *areab /2. / vzp / vzp/vzp
      bata =-alt *jizp - bt *j2zp
      gamma =aif + a2t * j3zp /vzp
c     Page 72
      delta =bata + alt * j3zp / vzp + bt *j4zp /vzp
c
c     calculate base pfessure
c
      pbase = (pmean - delta) / gamma
c
c     calculate breech pressure
c
      pbrch = alf * pbase + bata
      go to 610
c
c     use 2 phase gradient equation
c
470   if(iswl.ne0O)goto 480
      pbase = pmean
      pbrch = pmean
      if(pbase.gt.resp + 1)iswl = 1
      go to 620
480   if(iswl.eq.2)go to 580
      vzp = chain + areab * (y(2) +y()
      vzb = chain + areab * (y(10) +y()
      phi = (
      phidot = 0.
      dmorho = 0.
      dmcov = 0.
      dmromw = 0.
      rmomw = 0.
      vfree = vzp - vi
      do 490 k 1, nprop
      rmoinw =rmomw + chwp(k) * frac(k) forcp(k) /tempp(k)
      phi = chwp(k) * frac(k) + phi
      if(ibo(k).eq.l)go to 490
      dmorho = dmorho + tng(k) surf(k) *z(ibrp + k)
      phidot = rhop~k W tng (k) * surf~k M z (ibrp + k) + phidot
      dmcov =rhop (k) *tng (k) * surf (k) z z(ibrp + k) * covp (k)
     1 + dmcov
      dmroinw = dniroinw + rhop (k) * tng (k) *surf (k * z (ibrp + k)
     1 * forcp(k) / tempp(k)
      490 continue
      rmoinw =rmoinw + chwi * forcig Itempi
      gasmas =phi + chwi
      gasden =gasmas / vfree
      phi = (phi + chwi) / tinpi
      if (phi.gt.0.999) then
      iswl = 2
      rbin = phase /pinean
      rbrm = pbrch Ipmean
      if(phi.ge.l.)go to 580
      endif
      dmdt = phidot
      phidot = phidot / tinpi
      vdotov = (dxnorho + areab * y(l)) /vfree
      dlnrho = dmdt / gasinas - vdotov
      dvoldt = dmorho + areab * y(l) - dmcov
c
c     get time derivative of m-ean pressure
c
      dpindt = (dmroinw * tgas - pinean * dvoldt + dtgas *rinomw) /volg
c     Page 73
      voiprp -0.
      ef fdia =0.
      dxndidt = 0.
      dxndmor = 0.
      avelen = 0.
      avedia =0.
      do 500 k -1, nprop
      if(ibo(k).eq.1)go to 500
      volprp = volprp + (1. - frac(k)) *chwp(k) /rhop(k)
      dmdrndt = dmdrndt + rhop(k) *tng(k) *dsdx(k) *z(ibrp + k)
     1 * z(ibrp + k)
      dmdmdt - dnidxndt + rhop(k) *tng(k) *surf(k) *d2xdt2(k)
      drndmror = dmdmor + (dsdx(k) *z(ibrp + k) '~2 surf(k)*
     1 d2xdt2(k)) * tng(k)
      effdia - effdia + 6. * voip(k) / surf(k) *(1. - frac(k))
     1 * chwp(k)
500 continue
      cit = dxndrndt / gasxnas - dni mor / ",free + vdotov ** 2 -(dmdt
     1 I gasmas) ** 2
      d2lnr = cit - areab ** 2 * pbase / vfree / prwt
      d2lnr = d2inr + areab * areab * resp / vfree / prwt
      zp = chmlen + y(2) + y(7)
      zb = chrnien + y(10) + y(7)
      ulien =zp -zb
      cnow =tmpi -gasrnas
      vp = y(l)
      effdia =effdia / cnow
      prden =cnow /volprp
      up = Y(9)
      phistr = phi -gasden * areab * ulien / tmpi
      ulldot = vp -up
      dphist = phidot -gasden * areab / trnpi * (ulidot + ullen*
     1 dinrho)
      eps = 1. - (1. -phi) * tmpi / prden / vzb
      epsdot =phidot *trnpi /prden /vzb + (1. - phi) *trnpi *up
     1 * areab / prden /vzb /vzb
      ug = up + (vp + ullen * dinrho -up) / eps
      alam = (1.5 * grien / grdiamn) **.666666667
      alam = (0.5 + grien / grdiarn) /alam
      alam= alam ** 2.17
c
c     vis kg/s/rn
c
      vis = .00007
      ren = gasden / vis * effdia * abs(ug - up)
      if(ren.lt.1.)ren =1.
      fsrg = 2.5 * alam /ren ** .081 * ((l. - eps) /(1. -epsO)*
     1 epsO / eps) ** .45
      fsc = fsrg *fs0
      phi2 = 1. -phi - phistr * (1. - eps) / eps
      philp =dphist * ug -phidot *up - phistr * epsdot /eps/
     1 eps *(vp + ullen *dinrho - up) + phistr * ulidot *dlnrho
     1 /eps +2. * phistr *ugl/ zb* (ug -up)
      philp =philp + phi2 * asden /effdia / prden*
     1 (ug-up) ** 2*fsc
      ak2 = 1. / (1. -phi2 *trnpi /prden / vzb)
c     Page 74
      Phil = Philp + phistr * z(1) / eps + ullen * phistr * d2lnr/
     1 eps
c
c acceleration of forward boundary of propellant bed
c
      z(9) = gasden *(ug - up) ** 2 *fsc / prden / effdia + tmpi
      1 * Phil *ak2 /vzb / prden
      Z(10) =y(9)
      e = phistr /eps * (1. - ullen *areab /vfree) * areab / prwt
      dd = ullen *phistr * cit /eps
      akil = tmpi * e- * ak2 / zb Ivzb
      akl2 =tmpi * ak2 * (Philp + dd) / zb /vzb -akil * resp
      pbase = pmean - akl2 * zb *zb /2. + gasden *ullen * resp*
     1 areab / prwt
      pbase = pbase + akl2 * zb *zb *(zb I3. + ullen) /2. /zp
      pbase = pbase - gasden * ullen **2 *areab * resp /2. /zp
     1 / prwt
      pbase = pbase - gasden * ullen **2 /2. * (1. - 2. *ullen/
     1 3. / zp) * (cit - dlnrho **2)
      pbase =pbase - areab ** 2 *gasden *ullen ** 2 * (1. -2.*
     1 ullen/3. / zp) * resp /prwt /vfree /2.
      deno = -akil * zb ** 3 / 6. /zp - ullen * akil * zb *zb/
     1 2. / zp
      deno = deno + gasden * ullen *areab / prwt - areab ** 2
     1 gasden *ullen ** 2 * (1. - 2. * ullen / 3. / zp) /2./
     1 vfree /prwt
      deno = deno -gasden * ullen w* 2 * areab /2. /zp /prwt +
     1 1. + akl * zb * zb / 2.
      pbase = pbase / deno
      if (ism. eq.0) goto530
      if (ism. eq.1) goto5l0
      goto520
510   ism =2
      tss =sqrt(egama * rmomw / gasmas * tgas)
      write(6, * )tss
      tss = ullen / (ullen * odlnr + tss)
      tso =y(3)
      write(6, * )tss,, tso
520   coefbp = (tss + tso - y(3) - deltat) / tss
      if(coefbp.gt.l.)coefbp = 1.
      if (coefbp .le. 0.) then
      coefbp = 0.
      ism = 0
      endif
      pbase = coefbp * opbase + (1. - coefbp) * pbase
      if(mode.eq.1)write(6, * )coefbp, opbase, pbase, ism
      if (mode .eq.2) then
      argi opbase /6894.757
      arg2 = pbase /6894.757
      write(6, * )coefbp, argi, arg2, ism
      endif
530   odlnr = dlnrho
      opbase pbase
      pbrch = pbase *(1. + akil zb * zb / 2. + gasden * ullen
     1 areab /orwt -areab ** 2 *gasden * ullen ** 2 / 2. / vfree
     1 /prwt)
c     Page 75
      pbrch = pbrch + akl2 * zb * zb /2. -gasden * ullen * areab
     1 * resp / prwt
      pbrch =pbrch + gasden *ullen **2 /2. * (cit - dlnrho **2)
      pbrch - pbrch + areab **2 * gasden *ullen **2 *resp I2.
     1 / vfree / prwt
      go to 610
c
c     using rga gradient c
540   if(iswi.ne.0)go to 550
      pbase =pmean
      pbrch =pmean
      if(pbase.gt.resp + 1.)iswi 1
      go to 620
550   if(iswl.eq.2)go to 580
      vzp =chain + areab * (y(2) +y()
      vzb = chain + areab * (y(10) +y7)
      jlzb = bint(1) + (bvol * pb + areab /2. *pb *pb) /areab
      j2zb = (bvo'. + areab * pb) ** 2 Iareab /areab
      j3zb = bint(3) + areab *bint(i) *pb + bvol *pb *pb /2. +
     1 areab / 6. * pb ** 3
      phi = 0.
      phidot -0.
      dmorho -0.
      dmcov = 0.
      dinromw =0.
      rmomw = 0.
      vfree = vzp - vi
      do 560 k =1, nprop
      rmoinw =rmomw + chwp(k) *frac(k) forcp(k) /tempp(k)
      phi = chwp(k) * frac(k) + phi
      if(ibo(k).eq.i)go to 560
      dmorho = dmorho + tng(k) *surf(k) *z(ibrp + k) phidot = rhop (k) *tng (k) * surf (k) z z(ibrp + k) + phidot
      dmcov = rhop~k W tng W) * surf (k z z(ibrp + k) * covp (k)
     1 + dmtcov
      dmr omw = dmromw + r hop (k) * t ng (k) s sur f (k) z z(ib rp +t k)
     1 * forcp(k) / tempp(k)
      560 continue
      rmomw =rmomw + chwi * forcig /tempi
      gasmas =phi + chwi
      gasden gasrnas / vfree
      phi = (phi + chwi) / tinpi
      if (phi.gt.0.99) then
      iswi = 2
      rbm = pbase /pinean
      rbrm = pbrch /pinean
      if(phi.ge.i.)go to 580
      endif
      dxndt = phidot
      phidot -phidot / tinpi
      vdotov = (dinorho + areab * yMi) /vfree
      dlnrho =dmdt / gasmas - vdotov
      dvoldt =drnorho + areab * yMi - drncov
c
c     get time derivative of mean pressure
c     Page 76
c
      dpmdt =(dm~romw * tgas - pmean * dvoldt + dtgas * rmornw) /volg
      volprp =0.
      effdia 0.
      dxndxdt =0.
      dmdmor =0.
      avelen 0.
      avedia =0.
      do 570 k = 1, nprop
      if(ibo(k).eq.1)go to 570
      volprp = volprp + (1. - frac(k)) * chwp(k) / rhop(k)
      dxndrdt = dmdmdt + rhop (k) *tng (k) *dsdx (k) z z(ibrp + k)
     1 * z(ibrp + k)
      dxndmrdt = drndrdt + rhop(k) *tng(k) *surf(k) *d2xdt2(k)
      dmdmor = dxndmor + (dsdx(k) *z(ibrp + k) **2 + surf (k)*
      d2xdt2(k)) * tng(k)
      effdia - effdia + 6. * volp(k) / surf(k) *(1. - frac(k))
     1 * chwp(k)
      570 continue
      clt =dmdrndt / gasmas - dmdrnor / vfree + vdotov ** 2 -(dmdt
     1 / gasmas) ** 2
      d2lnr = cit - areab ** 2 * pbase / vfree / prwt
      d2lnr =d2lnr + areab * areab * resp / vfree / prwt
      zp = chmlen + y (2) + y (7)
      zb = chmlen + y (10) + y (7)
      ullen zp -zb
      cnow =tmpi -gasmas
      vp = y(l)
      effdia =effdia / cnow
      prden cnow /volprp
      up = yý'9)
      phistr = phi -gasden * areab * ullen / trnpi
      ulidot = vp -up
      dphist = phidot -gasden * areab / tmpi * (ulidot + ullen*
      dlnrho)
      eps = 1. - (1. -phi) * tmpi / prden / vzb
      epsdot = phidot *trnpi /prden /vzb + (1. - phi) *trnpi *up
     1 * areab / prden /vzb /vzb
      ug = up + (vp + ullen * dlnrho -up) / eps
      alam = (1.5 * grien /grdiam) **.666666667
      alarn = (0.5 + grien /grdiarn) /alarn
      alam =alarn ** 2.17
c
c     vis kg/s/rn
c
      vis =.00007
      ren =gasden / vis * effdia * abs(ug - up)
      if(ren.1t.1.)ren =1.
      fsrg = 2.5 * alam /ren **.081 * ((1. - eps) /(1. -epsO)*
      epsO / eps) ** .45
      fsc =fsrg *fsQ
      phi2 = 1. -phi - phistr *(1. - eps) / eps
      philp dphist * ug - phidot * up - phistr * epsdot /eps/
      eps *(vp + ullen * dlnrho - up) + phistr * ulidot *dlnrho
      /eps +2. * areab *phistr * ug /vzb* (ug -up)
      philp = philp + phi2 * gasden / effdia /prden*
c     Page 77
     1 (ug -up) **2 *fsc
      ak2 =1. / (1. -phi2 * tmpi /prden / vzb)
      phil =philp + phistr * z(1) /eps + ullen * phistr * d2lnr/
     1 eps
c     acceleration of forward boundary of propellant bed
c
      z(9) = gasden * (ug -up) ** 2 * fsc /prden / effdia + tmpi
     1 * phil *ak2 / vzb /prden
      Z(10) =Y( 9)
      phi3 phistr * ug *ug + (1. - phi) *up * up
      e= 1 . - ullen * areab / vfree
      dd =ullen * phistr *cit /eps
      alt =tmpi * areab Ivzb /vzb * (phi3 * areab /vzb - (philp
     1 + dd -e * phistr *areab *resp / eps / prwt) *ak2)
      a2t = (- tmpi * e *phistr *areab **2 / vzb /vzb / eps/
     1 prwt) *ak2
      bt =- tmpi * phi3 *areab* 2/2. /vzb /vzb/ vzb
      pbase =pmean -gasden * ullen **2 /2. * (cit -dlnrho **2
     1 + areab ** 2 *resp /prwt /vfree) *(1. - 2. *areab*
     1 ullen /3. / vzp)
      pbase =pbase - alt *j3zb /vzp - bt *j4zb / vzp - areab*
     1 ullen *alt * jlzb /vzp
      pbase =pbase - areab * bt *ullen * j2zb /vzp - gasden*
     1 areab ** 2 *ullen ** 2 * resp / 2. / vzp /prwt + alt*
     1 jlzb + bt *j2zb + areab * gasden * ullen *resp / prwt
      deno =1. -* areab *ullen * a2t * jlzb / vzp - gasden * areab
     1 **2 *ullen ** 2 /2. / vzp / prwt + a2t *j3zb / vzp - a2t
     1 *jlzb + gasden * ullen * areab /prwt
      demo = deno - gasden * ullen ** 2 *areab **2 /2. /vfree/
     1 prwt + gasden *areab ** 3 * ullen ** 3 /3. /vzp /vfree/
     1 prwt
      pbase = pbase /deno
      pbrch = pbase *(1. - a2t * jlzb + gasden *ullen * areab/
     1 prwt - gasden *ullen **2 * areab ** 2 / 2. / vfree / prwt)
     1 + gasden * ullen ** 2 /2. * (cit. - dlnrho ** 2 + areab **2
     1 * resp Iprwt /vfree) - alt, * jlzb -bt * j2zb -areab*
     1 gasden *ullen *resp Iprwt
      go to 610
580   pbase = rbm *pmean
      pbrch = rbrm *pmean
      go to 610
582   areazm=areab
      areaza=areazm- pi * btdia**2/4.
c     writek(G,*)areazm, areaza
delta=l.
      qlzp= (vzp**2- (areab*za) **2) /2. /areaba**2
      q2zp=vzp* *2/areaba/areaba
      q3zp= (zp-za) /areaba
      q4zp=vzp/areaba/areaba
      q5zp=l ./areaba/areaba
      q6zp=(areab*za**3/3.+(areab*za+areaba*(zp-za))**3/3./areaba
     &**2-(areab*za)**3I3.Iareaba**2)
      q8zp= (zp-za) **212.
      q9zp=(vzp**3-(areab*za)**3)/6./areaba**2
     & -(areab*za)**2*(zp-za)/2./areaba
      qlza=za**2/2.
c     Page 78
      q2za-vza* *2 /areab/areab
      q7za-areab*za**3/6.
      vp-y (1)
      bt= -tmpi * y(l)*y(l)*areab*areab/2./vzp/vzp/vzp
      hl= tmpi*areab*areaa*y(l) *y (1)!(vzp*vzp)
      akl=-tmpi *areab*vp*vp*vza*areaa/vzp/vzp/areaza* *2
      akl=akl+tmpi*areaa**2*vp**2/vzp/areaza**2/2.
      akl=akl+tmpi*vp*vp*areab**2*vza**2/2 ./vzp**3/areaza**2
      ajl=-tmpi*areaa**2*vp*vp/2. /vzp
      fk= -2. *tmpi*vp**2*areazm*areaa/areaza/vzp*
     & (areab*vza/vzp/areazm-1l 4**2
      fk= fk/ (areaza+areazm)
      rrl= 1. + tmpi*areab*areaa*qlza/prwt/vzp**2
      tal-tmpi*areab**2*y (1)**2/vzp**3
      ta4=tmpi*areab**2*resp/prwt/vzp**2
      zal= (bt*q2za + (tal+ta4)*qlza)/rrl
      zaO= 1./rnl
      za2= - (tmpi*areab*areaba*qlza/prwt/vzp**2 ) /rrl
      a3 = tmpi*areab**2*y(l)**2/vzp**3 - tmpi*areab*areaa*zal/prwt/
     & vzp**2 + tmpi*areab**2*resp/vzp**2/prwt
      a4l= -tmpi*areab*areaa*za2/vzp**2/prwt
      a42= -tmpi*areab*areaba/vzp**2/prwt
      a4 = a41+a42
      c a4 - -tmpi*areab*areaa*za2/vzp**2/prwt - tmpi*areab*areaba/
      c & vzp**2/prwt
      a5 = - tmpi*areab*areaa*zaO/prwt/vzp**2
      c3 = tmpi*areaa**2*zal/prwt/vzp - tmpi*areaa*areab*
     & resp/prwt/vzp - tmpi*areaa*areab*y(l)*y(l)/vzp/vzp
      c41 = tmpi*areaa*areaa*za2/vzp/prwt
      c42 = tmpi*areaa*areaba/vzp/prwt
      c4 = c41+c42
c     c4 = delta*tmpi*areaa*areaa*za2/vzp/prwt + delta*tmpi*areaa*
c     & areaba/vzp/prwt
      c5 = tmpi*areaa*areaa*zaO/vzp/prwt
      11 = zal+fk+a3*ql zp+c3*q3zp+bt*q2zp+hl*q4zp+ajl*qszp+akl
      12 = 1-za2-a4*qlzp - c4*q3zp
      13 = zaO + a5*qlzp +c5*q3zp
      bi = (a3*q7za+a3*q9zp+bt*q6zp+zal*(vzp-vza)+fk*(vzp-vza)
     1 *+c3*q8zp+h1*qlzp+ajl*q3zp+akl*(vzp-vza))/vzp
      b2 = (a4*q7za+a4*q9zp+za2* (vzp-vza) +c4*q8zp) /vzp
      b3 = (vza+a5*q7za+a5*q9zp+zaO* (vzp-vza) +c5*q8zp) /vzp
c     calculate base pressure
      pbase=(pmean/b3 - bl/b3 + 11/13)/( 12/13 + b2/b3)
c     calculate breech pressure
      pbrch=pmean/b3 - bl/b3 - pbase*b2/b3
      pza=zal + za2*pbase + zaO*pbrch
c     write(6,7)y(3),z(l) ,y(l),y(2),pmean,pbase,pbrch
c     write(6,7) bl,b2,b3,ll,12,13
c     write(6,7) a3,a4,a5,c3,c4,c5
      z(l)= (areaa*pza+areaba*pbase-areab*resp) /prwt
      go to 615
585   If(za.gt.chdist(nchpts))goto 275
      Do 269 I=2,nchpts
      If(za.lt.chdist(I).and.za.qt.(chdist(I-l)))goto 274
269   continue
274   diam =( za - chdist(I-l)) /(chdist(I) - chdist(I-l))*
     & (chdiam(I) - chdiam(I-l)) + chdiam(I-l)
      areazm = pi * diam**2/4.
c     Page 79
      areaza - areazm - pi * btdia**2/4.
275   continue
      vp-y(1)
      bt- -tmpi * y (1)*y(l) *areab*areab/2 ./vzp/vzp/vzp
      hl= tmpi*areab*areaa*y(1) *y(l) /(vzp*vzp)
      ak1=-tmpi*areab*vp*vp*vza*areaa/vzp/vzp/areaza* *2
      akl=akl+tmpi*areaa**2*vp* *2/vzp/areaza**2/2.
      akl=akl+tmpi*vp*vp*areab**2*vza**2/2./vzp**3/areaza**2
      ajl=-tmpi*areaa**2*vp*vp/2. /vzp
      fk= -2. *tmpi*vp**2*areazm*areaa/areaza/vzp*
     & (areab*vza/vzp/areazm-1.) **2
      fk= fkl (areaza+areazm)
      rrl= 1. + tmpi*areab*areaa*qlza/prwt/vzp**2
      tal=tmpi*areab**2*y (1)**2/vzp**3
      ta4=tmpi*areab**2*resp/prwt/vzp**2
      za1l= (bt*q2za + (tal+ta4)*qlza)/rrl
      zaO= 1./rrn
      za2= -(tmpi*areab*areaba*qlza/prwt/vzp**2 ) /rrl
      a3 - tmpi*areab**2*y(1)**2/vzp**3 - tmpi*areab*areaa*zal/prwt/
     & vzp**2 + tmpi*areab**2*resp/vzp**2/prwt
      a41= -tmpi*areab*areaa*za2/vzp**2/prwt
      a42- -tmpi*areab*aieaba/vzp**2/prwt
      a4 = a4l+a42
c     a4 - -tmpi*areab*areaa*za2/vzp**2/prwt - tmpi*areab*areaba/
c    & vzp**2/prwt
      a5 = - tmpi*areab*areaa*zaO/prwt/vzp**2
      c3 = tmpi*areaa**2*zal/prwt/vzp - tmpi*areaa*areab*
     & resp/prwt/vzp - tmpi*areaa*areab*y(l)*y(1)/vzp/vzp
      c4l = tmpi*areaa*areaa*za2/vzp/prwt
      c42 = tmpi*areaa*areaba/vzp/prwt
      c4 = c43+c42
c     c4 = delta*tmpi*areaa*areaa*za2/vzp/prwt + delta*tmpi*areaa*
c    & areaba/vzp/prwt
      c5 = trnpi*areaa*areaa*zaO/vzp/prwt
      11 = zal+fk+a3*qlzp+c3*q3zp+bt*q2zp+hl*q4zp+ajl*qszp+ak1
      12 = 1-za2-a4*qlzp - c4*q3zp
      13 = zaO + a5*qlzp +c5*q3zp
      bi = (a3*q7za+a3*q9zp+bt*q6zp+zal* (vzp-vza) +fk* (vzp-vza)
     *+c3*q8zp+hl*qlzp+ajl*q3zp+akl* (vzp-vza) )/vzp
      b2 = (a4*q7za+a4*q9zp+za2*(vzp-vza)+c4*q8zp)/vzp
      b3 = (vza+a5*q7za+a5*q9zp+zaO* (vzp-vza) +c5*q8zp) /vzp
c     calculate base pressure
      pbase=(pmean/b3 - bl/b3 + 11/13)/( 12/13 + b2/b3)
c     calculate breech pressure
      pbrch=pmean/b3 - bl/b3 - pbase*b2/b3
      pza=zal + za2*pbase + zaO*pbrch
c     write (6,7) y(3),z (1), y(1), y(2) ,pmean,pbase,pbrch
c     write(6,7) bl,b2,b3,ll,.12,13
c     write(6,7) a3,a4,a5,c3,c4,c5
      z (1) =(areaa*pza+areaba*pbase-areab*resp) /prwt
      go to 615
c     use lagrange pressure gradient equation
c
590   if(iswl.ne.0)go to 600
      if(pmean.lt.resp)resp = pmean
c
c     calctilate base pressure C.
c     Page 80
600 tmp2 = 1.0 + tmpi. / 2.0 / prwt
      tmr2 = 1.0 + tmpi / 2.0 / rcwt
      tmr3 = 1.0 + tmpi / 3.0 / rcwt
      tmr4 = rfor / areab + resp -pgas *air
      pbase = pmean / tmr3
      1 - tmpi / 2.0 /tmr2 * (tmr4 /rcwt - resp / prwt)
      2 + tmpi / 3.0 Itmr3 * (tmr4 /rcwt - resp / prwt / 2.0)
      pbase = pbase I(tmp2 / tmr2 - tmpi / tmr3 / prwt /6.0)
c
      if(pbase.gt.resp + 1.)iswl =1
c     calculate breech pessure
c
      pbrch =pbase * tnip2 / tmr2 + trpi / 2.0 /trnr2*
     1 (tmr4 / rcwt - resp / prwt)
c
c     calculate projectile acceleration
c
610   z(1) - areab * (pbase - resp) /prwt
615   if(z(1).lt.0.0)go to 620
      go to 630
620   if(iswl.eq.0)z(1) =0.0
630   if (y (1) . lt. 0. 0)y (l) = 0. 0
      z(2) = y(l)
c
c     get burning rate
c
640   do 670 m = 1, nprop
      z(ibrp + mn) = 0.0
      d2xdt2(m) = 0.0
      if(ibo(m).eq.1) goto 670
      do 650 k = 1, nbr(m)
      if(pmean.gt.pres(m, k))go to 650
      go to 660
      650 continue
      k = nbr(m)
      660 pmix = pinean
      if(igrad.eq.3)pmix = pbrch - (akil * pbase + akl2) /6.*
     1 zb * zb
      if(igrad.eq.4)pmix = pbrch + (alt + a2t * phase) * j3zb/
     1 vzb + bt * j4zb / vzb
      if(pmix.lt. .99 * pmean)pmix = pmean
      z(ibrp + in) =beta(m, k) * (pinix * i.e - 6) **alpha(in, k)
      abr(in) = alpha(m, k)
      bbr(m) = beta (m, k)
      d2xdt2(in) = beta(in, k) * alpha(m, k) * (pinix * .e -6) *
     1 (alpha(m, k) -1.) * dpindt * i.e - 6
670   continue
      do 680 i = 1, nde
      d(i) = (z(i) -b(j) * p(i)) * a(j)
      y(i) = deltat *d(i) + y(i)
      P(i) = 3. * d(i) - ak(j) * z(i) + p(i)
680   continue
690   continue
      nzp=nzp+1
      if(igrad.ne.6)goto 2003
      if(nzp.ne.nzpi)goto 2003
      dzp=zp/50.
c     Page 81
      open (unit=11, file=' dist .dat')
      open(unit=12,file='press.dat')
      do 2001 i=1,50
      ddzp-i*dzp
      call jint(btd,btl,zp,ddzp,nchpts,chdist,chdiarn,bint,bvol)
      if(ddzp.lt.za)then
      pz-pbrch.-tmpi*areab*z (1)*bint (1)/vzp/vzp
     & +tal*bint(1)+bt*bint(2)
      pzl-pbrch+ (a3+a5*pbrch+a4*pbase) *bint (1)+bt*bint (2)
      else
      pz=pza +fk-tmpi*areab*z(l)*bint(10)/vzp**2
     &+tmpi*areaa*z (1)*bint (5) /vzp+bt*bint (2) +akl
     & +h1*bint(6)-~bt*2.*bint(10)-h1*bint(5)+aj1*bint (7)
      pz1=~za0*pbrch+zal+za2*pbase+fk+a3*bint (10)
     &+a4*pbase*bint(10)+a5*pbrch*bint(l0)
     &+c3*bint (5) +c4*pbase*bint (5)+c5*pbrch*bint (5)
     &+bt*bint (2) +h1*bint (6) +ajl*bint (7) +akl
      endif
      write(6,*)vza ',za,' ddzp ',ddzp,' pza ',pza
      write (11,*) ddzp
      write (12, *)pz/1 .e6
2001  continue
2003  if(prwt0.ne.prwt)then
      if(mode.eq.1)write(6, 1450)prwt
      if(mode.eq.2)then
      argl = prwt / 0.45359237
      write(6, 1450)argl
      endif
      prwtO =prwt
      lines = lines + 1
      endif
      t = t + deltat
      told = y(3 )
      if(ihl.eq.1 .and. pmean.gt.burstp)then
      write(6, 1440)
      ihl = 2
      lines = lines + 1
      endif
      if(pmaxm.gt.pmean)go to 700
      pmaxrn pmean
      tpmaxrn y(3)
700   if(pmaxba.gt.pbase)go to 710
      pmaxba = pbase
      tpmxba = y(3)
710   if(pmaxbr.gt.pbrch)go to 720
      pmaxbr = pbrch
      tpmxbr = y(3)
720   continue
      if(y(3).lt.ptime)go to 730
      ptime = ptime + deltap
      pjt =y( 2 ) + y(7 )
      argO =y(3) * 1000.
      if(mode.eq.1)write(6, 1270)arg0, z(1), y(l), pjt, pmean, pbase,
     1 pbrch
      if(nzp.ne.nzpi)goto 2004
      write(6,*)'dpza ',dpza,' pza ',pza,' dl ',dl,' gi ',gl
c     write(6, *) areaa*pza,areaba*pbase,areaa*pza+areaba*pbase
c     Page 82
      write(6,*)Iareaa ',areaa,' areab ',areab,' areaba ',areaba
      write(6,*)' vza ',vza,' vzp ',vzp
      write (6, *) 'qiza-qiOza'
      write (6, *) qiza, q2za, q3za, q4za
      write (6, *) q5za, q6za, q7za, q8za
      write (6, *) q9za, qiOza
      write (6, *)' qlzp-q9zp'
      write(6, *)qlzp,q2zp,q3zp,q4zp
      write(6, *)qszp,q6zp,q7zp,q8zp
      write (6, *)q9zp
      write(6, *)
      write(6,*)Izao,zal,za2,resp ',zaO,zal,za2,resp
      write (6, *) 'b1,b2,b3, 11 ',bl,b2,b3, 11
      write(6,*) 'l2,l3,a3,a4 ',l2,l3,a3,a4
      write(6,*)I a5,c3,c4,c5 ',a5,c3,c4,c5
      write(6,*)vbt,hl,akl,ajl,fk ',bt,h1,akl,ajl,fk
2004  if(mode.eq.2)then
      argi - y(l) /0.0254
      arg2 - pjt /0.0254
      arg3 = pmean /6894.757
      arg4 = pbase /6894.757
      arg5 = pbrch /6894.757
      arg6 =z(1) /0.0254
      write(6, 1270)arg0, arg6, argi, arg2, arg3, arg4, arg5
      endif
      lines = lines + 1
      if (igrad.gt .2)then
      pjt = y(2) + y(7)
      prt = y(10) + y(7)
      if(mode.eq.1)write(6, 1280)prt, pjt
      if (mode .eq. 2) then
      argi = prt /3.28083
      arg2 = pjt /3.28083
      write(6, 1280)argl, arg2
      endif
      lines = lines + 1
      endif
730   continue
      if (lines.gt .linmax) then
      write(6, 1236) title,vsn
      write(6, 1230)
      if(mode.eq.l)write(6, 1232)
      if(mode.eq.2)write(6, 1234)
      lines = 4
      endif
      if(t.gt.tstop)goto 740
      if(y(2) + y(7).gt.travp)go to 740
      rmvelo = y(l)
      tmvelo = y(23)
      rcvelo = y(6)
      disto = y(2) + y(7)
      go to 250
740   if(lines.gt.linmax-nprop-25) write(6, 1236) title,vsn
      write(6, 1290)t, y(3)
      if (mode. eq. 1) write (6, 1300) -maxm, tpmaxm, pmaxba, tpmxba, ?maxbr,
      1 tpmxbr
      if (mode.eq. 2) then
      pmaxm = pmaxm / 6894.757
c     Page 83
      pinaxba -pmaxba / 6894.757
      pmaxbr = pmaxbr / 6894.757
      write(6, l3lO)pmaxm, tpmaxrn, pmaxba, tpmxba, pmaxbr, tpmxbr
      endif
      if(y(2) + y(7).le.travp)goto 750
      dfract =(travp - disto) /(y(2) + y(7) - disto)
      rmvel =(y(l) - rmvelo) *dfract + rmvelo
      tmvel (y(3) - tmfvelo) *dfract + tmvelo
      rcvel -(y(6) - rcvelo) *dfract + rcvelo
      if(mode.eq.1)write(6, 1320)rmvel, tmvel, rcvel
      if (mode.eq.2)then
      rmvel - rmvel * 3.28083
      rcvel = rcvel * 3.28083
      write(6, 1330)rmvel, tmvel, rcvel
      endif
      go to 760
750   if(mode.eq.1)write(6, l3 40)y(l), y(3)
      if (mode .eq. 2) then
      argi = y(1) * 3.28083
      write(6, 1350)arql, y(3)
      endif
760   efi = chwi * forcig / (gamai - 1.)
      efp = 0.0
      do 770 i 1, nprop
      efp =efp + chwp(i) * forcp(i) /(gaxnap(i) -1.0)
770   continue
      tenerg =efi + efp
      if(mode.eq.1)write(6, 1360)tenerg
      if(mode.eq.2)tenerg - tenerg / 0.1129848
      if(mode.eq.2)write(6, 1370)tenerg
      tengas = chwi * forcig *tgas / (gamai - 1.) / tempi
      do 780 i - 1, nprop
      tengas =(frac(i) *chwp(i) * forcp(i) * tgas /tempp(i)/
      1 (gamap(i) - 1.)) + tengas
780   continue
      write(6, 1380) (i, frac(i), tbo(i), i =1, nprop)
      if(mode.eq.1)write(6, 1390)
      if (mode.eq.2)then
      tengas = tengas / 0.1129848
      elpt = elpt /0.1129848
      elpr = elpr /0.1129848
      elgpm = elgpm / 0.1129848
      elbr = elbr /0.1129848
      elrc = elrc /0.1129848
      elht = elht /0.1129848
      31ar = elar /0.1129848
      write (6, 1400)
      endif
      pctenl - tengas / tenerg * 100.0
      pcten2 = elpt /tenerg * 100.0
      pcten3 = elpr /tenerg * 100.0
      pcten4 = elgpm / tenerg * 100.0
      pcten5 = elbr /tenerg * 100.0
      pcten6 = elrc /tenerg * 100.0
      pcten7 = elht /tenerg * 100.0
      pcten8 = elar /tenerg * 100.0
      write(6, 1410)tengas, pctenl, elpt, pcten2, elpr, pcten3, 1 elgprn, pcten4, elbr, pcten5, elrc, pcten6, elht, pcten7,
c     Page 84
      2 elar, pcten8
      stop
790   write( *, 1420)
      stop
800   write( *, 1430)
810   continue
820   continue
      stop
830   format (' input name of data file to be used as input ')
840   format (alO)
850   format (' input name of output file ')
860   format (/' the input file is ',alO/)
870   format (' input data units - "metric" or "english"')
880   format (' must use quoted "im" or "e" as first character of',
     1 ' input file'!' to specify metric or english input units')
885   format (15a4)
890   format (lx,'using lagrange pressure gradient')
900   format (lx,'using chambrage pressure gradient')
910   format (lx,'using rga gradient')
920   format (///,' chamber distance m chamber diameter m',/
     1 (5x,e14.6,5x,e14.6))
925   format (///,' chamber distance in chamber diameter in',/
     1 (5x,e14.6,6x,e14.6))
930   format (lx,'use first 5 points')
940   format (lx,' # points ? ')
950   format (lx,'using 2 phase gradient equation')
960   format (/" chamber volume m**3 ',e16.6/
     1 ' groove diam m ",e16.6/" land diam m ',e16.6/
     1 ' groove/land ratio ',e16.6/" twist turns/caliber ",e16.6/
     1 ' projectile travel m',e16.6//' gradient # ',i3,//
     1 ' variable mass switch',i3/' container switch',i7!
     1 ' friction factor ',e18.6/)
970   format (/' chamber volume in**3 ",e16.6/
     1 ' groove diam in ',e16.6/' land diam in ",e16.6/
     1 ' groove/land ratio ',e16.6/' twist turns/caliber ',e16.6/
     1 ' projectile travel in',el6.6//' gradient # ",i3,//
     1 ' variable mass switch',i3/' container switch',i7/
     1 ' friction factor ',e18.6/)
980   format (lx,'number of variable projectile mass points ',i2,/
     1 lx,' travel (m) projectile mass (kg)'/
     1 (1x,e14.6,e14.6))
985   format (lx,'number of variable projectile mass points ',i2,/
     1 lx,' travel (in) projectile mass (lb)'/
     1 (lx,e14.6,e14.6))
990   format (ix,'canister burst pressure (mpa)',el4.6/
     1 Ix,'canister volume (m^3) ',e14.6/
     1 lx,'canister diameter (m) ',e14.6)
1000  format (lx,'canister burst pressure (psi)',el4.6/
     1 lx,'canister volume (inA3) ',e14.6/
     1 lx,'canister diameter (in) ',e14.6)
1010  format (/' projectile mass kg',34x,e14.6/ 1 ' switch to calculate energy lost to air resistance ',i3/
     1 ' fraction of work against bore used to heat the tube',el4.6/
     1 ' gas pressure mpa ',e14.6)
1020  format (W' projectile mass lb',34x,e14.6/
     1 ' switch to calculate energy lost to air resistance ',i3/ 1 ' fraction of work against bore used to heat the tube',el4.6/
     1 ' gas pressure psi ',e14.6)
c     Page 85
1030  format (I' number barrel resistance points ',i2/
     1 ' bore resistance mpa - travel m '/(3x,e14.6,8x,el4.6))
1040  format (/' number barrel resistance points ',i2/
     1 ' bore resistance psi - travel inches '/(3x,e14.6,e22.6))
1050  format (1x)
1060  format (U
     1 ' mass of recoiling parts kg ',e14.6/
     1 ' number of recoil point pairs ',i3/
     1 ' recoil force J',' recoil time sec'/, (lx,e14.6,3x,e14.6))
1070  format (/
     1 ' mass of recoiling parts lb ',e14.6/
     1 ' number of recoil point pairs ',i3/
     1 ' recoil force in-lb',' recoil time sec'/
     1 (1x, e14.6,3x, e14.6))
1080  format (U
     1 ' free convective heat transfer coefficient w/m^2 k ',e14.6/
     1 ' chamber wall thickness m ' e14.6/
     1 ' heat capacity of steel of chamber wall j/kg k ',e14.6/
     1 ' initial temperature of chamber wall k ',e14.6/
     1 ' heat loss coefficient ' ,e14.6/
     1 ' density of chamber wall steel kg/m^3 ',e14.6/)
1090  format (U
     1 ' free convective heat transfer coef in-lb/inA2-s-k ',e14.6/
     1 ' chamber wall thickness (inches) ',e14.6/
     1 ' heat capacity of steel of chamber wall in-lb/lb-k ',e14.6/
     1 ' initial temperature of chamber wall k ',e14.6/
     1 ' heat loss coefficient ',e14.6/
     1 ' density of chamber wall steel lb/in^3 ',e14.6/)
1100  format (
     1 ' impetus of igniter propellant j/kg ',e14.6/
     1 ' adiabatic flame temperature of igniter propellant k',e14.6/
     1 ' covolume of igniter mA3/kg ',e14.6/
     1 ' ratio of specific heats for igniter ',e14.6/
     1 ' initial mass of igniter kg ',e14.6/)
1110  format
     1 ' impetus of igniter propellant ft-lb/lb ',e14.6/
     1 ' adiabatic flame temperature of igniter propellant k',e14.6/
     1 ' covolume of igniter ftA3/lb ', e14.6/
     1 ' ratio of specific heats for igniter ',e14.6/
     1 ' initial mass of igniter lb ',e14.6/)
1120  format (/P there are ',i2,' propellants'//)
1130  format ((' for propellant number',i2//
     1 ' impetus of propellant j/kg ',e14.6/
     1 ' adiabatic temperature of propellant k ',e14.6/
     1 ' covolume of propellant mA3/kg ',e14.6/
     1 ' ratio of specific heats for propellant ',e14.6/
     1 ' initial mass of propellant kg ',e14.6/
     1 ' density of propellant kg/mA3 ',e14.6/
     1 ' number of perforations of propellant ',i3/
     1 ' length of propellant grain m ',e14.6/
     1 ' perforation diameter m ',e14.6/
     1 ' outside diameter of propellant grain m ',e14.6/)/)
1140  format ((' for propellant number',i2//
     1 ' impetus of propellant ft-lb/lb ',e14.6/
     1 ' adiabatic temperature of propellant k ',e14.6/
     1 ' covolume of propellant in^3/lb ',e14.6/
     1 ' ratio of specific heats for propellant ',e14.6/
     1 ' initial mass of propellant lb ',e14.6/
c     Page 86
     1 ' density of propellant lb/in^3 ,e14.6/
     1 ' number of perforations of propellant ',i3/
     1 ' length of propellant grain in ',e14.6/
     1 ' perforation diameter in ',e14.6/
     1 ' outside diameter of propellant grain in',el4.6/)/) 1150 format (P' for propellant ',i2,' the total number of ains'
     1 ,' is ',e14.6)
1160  format (' number of burning rate points',i2/
     1 ' exponent coefficient pressure'/
     1 ' m/sec-mpa-ai mpa') 1170 format (' number or burning rate points',i2/
     1 ' exponent coefficient pressure'/
     1 ' in/sec-psi^ai psi') 1180 format (lx,e14.6,5x,e14.6,4x,e14.6)
1190  format (' time increment msec ',e14.6/
     1 ' print increment msec ',e14.6/
     1 ' time to stop calculation msec',el4.6) 1200 format (lx,'end input data -- i.b. calculation start')
1210  format (W' area bore m^2 ',e27.6/' pressure from i-7n pa',e21.6/
     1 ,' volume of unburnt prop m^3 ',e14.6/
     1 ' init cham vol-cov ign mA3 ',e15.6) 1220 format (I' area bore in^2 ',e29.6/' pressure from ign psi',e23.6/
     1 ,' volume of unburnt prop in^3 ',e16.6/
     1 ' init cham vol-cov ign in^3 ',e17.6) 1230 format (I' time accel velocity distance pr mean',
     1 ' pr base pr brch')
1232  format ( ' (ms) (m/s/s) (mIs) (M) (Pa) ',
     1 ' (Pa) (Pa)')
1234  format ( ' (ms) (in/s/s) (in/s) (in) (psi) ',
     1 ' (psi) (psi)') 1236 format (lhl,3x,15a4,' rga.',a4)
1240  format (' propellant',i2,' has slivered') 1250 format (' propellant',i2,' has burned out')
1270  format (lx,lp7ell.4)
1280  format (lx,'prop travel',ell.4,'proj travel',ell.4) 1290 format (/lx,' deltat t', e14.6, ' intg t',e14.6)
1300  format (P' pmaxmean pa ',Ipel4.7,' time at pmaxmean sec ',
     1 Opel6.6/' pmaxbase pa ',ipe14.7,' time at pmaxbase sec ',
     1 Opel6.6/' pmaxbreech pa ',lpel4.7,' time at pmaxbreech sec ',
     1 Opel4.6) 1310 format (W' pmaxmean psi',fl14.3,' time at pmaxmean sec ',
     1 Opel6.6/' pmaxbase psi',fl4.3,' time at pmaxbase sec ',
     1 Opel6.6/' pmaxbreech psi',fl4.3,' time at pmaxbreech sec ',
     1 Opel4.6) 1320 format (/lx,'muzzle velocity m/s ',e14.6,' time of muzzle exit
     1 e14.6,' sec'//Ix,'recoil velocity m/s ',e14.6) 1330 format (/lx,'muzzle velocity ft/s ',e14.6,' time of muzzle exit ',
     1 e14.6,' sec'//lx,'recoil velocity ft/s ',e14.6) 1340 format (I ' velocity of projectile m/s ',e14.6,' at time ',e14.6,
     1 ' sec') 1350 format (U ' velocity of projectile ft/s ',e14.6,' at time ',e14.6,
     1 ' sec')
1360  format (/lx,'total initial energy available j = ',e14.6/) 1370 format (/lx,'total initial energy available in-lb = ',e14.6/) 1380   format (' propellant mass fraction burnt time (sec)'/
     1 (4x,i2,9x,e14.6,5x,e14.6))
1390  format (/' ** energy summary **',23x,'joules',1lx,'Q') 1400 format (/' ** energy summary **',23x,' in-lb',llx,'%')
c     Page 87
1410  format (lx,'total eziirgy remain'ng in gas',llx,' = ',e14.6,fll.4
     1 /lx,'energy loss from projecuile translation = ',el4.6,fl1.4
     1 /lx,'t-nergy loss frcn projectile rotation = ',e14.6,fll.4
     1 /lx,'energy lost to gas and propellant motion - ',el4.6,fll.4
     1 /Ilx,'energy lost to bore resistance -',el4.6,fl1.4
     1 /lx,'energy lost to recoil = ',el4.6,fll.4
     1 /l,'energy loss from heat transfer = ',eI4.6,fll.4
     1 /ix,'eiergy lost to air resistance = ',e14.6,fll.4)
1420  format (Jx,'end of file encountcr')
1430  format l1x,'read or write error')
1440  format (' canister burst pressute achieved')
1450  forrat (' projectile mass ti;nsition point - new mass =
     1 lpell.4)
      end
      subroutine prt7lO(pd, gd, 1l, np, x, frac, vol, surf, dsdx)
      common nsl, kpr, fracsl(lO), dsdxsl(10), surfsl(10), nslplO),
     1 tsl(10), pbrch, pbase, pmean, bbr(10), abr(l0), deltat, yar(20),
     1 igrad
      dimension ts(10), coef(10)
      pi = 3.141593
      nsl = 0
c
c     pd=perforation diameter
c     gd=outer dia
c     gl=grain length
c     np=number of perfs
c
c     surf=output surface area
c     frac=output mass fraction of propellant burned
c
c     w = web = distance between perforation edges
c     = distance between outside perf edge and edge of grain
c
c     p = distance between perforation centers
c
c     xl = distance to inner sliver burnout
c
c     x2 = distance to outer sliver burnout (frac=l)
c
      if(np.eq.0) go to 70
      if(np.eq.l)go to 90
      if(np.eq.2)go to 210
      if(np.eq.7)go to 10
      if(np.eq.19)go to 110
      if(np.eq.15)go to 180
      write(6, 220)
      stop
      10 w = (gd - 3. * pd) / 4.
      d = w + pd
      sqr3 = sqrt (3.)
      xl = d / sqr3 - pd / 2.
      x2 = (14. - 3. * sqr3) * d / 13. - pd / 2.
      vO = pi/ 4. *gl * (gd *gd - 7. *pd * pd)
      sO = 2. *v0 /gl + pi * gl * (gd + 7. * pd)
      if (x.gt.w / 2. + .0000001) goto 20
      vol = pi / 4. * (gl - 2. * x) * ((gd - 2. * x) ** 2- 7. * (pd +
     1 2. * x) ** 2)
      surf = 2. * vol / (gl - 2. * x) + pi * (gl - 2. * x) * ((gd - 2.
c     Page 88
     1 * x) + 7. * (pd + 2. * x))
      frac - 1. - vol / v0
      dsdx - - 4 * pi* (gd + 7. *pd - 3. * gl + 18. *x)
      dsdxsl(kpr) = dsdx
      fracsl(kpr) = frac
      surfsl(kpr) = surf
      return
      20 nsl = 1
      coef(kpr) = 0.
      if(igrad.eq.l.or.igrad.eq.2)go to 40
      if(nslp(kpr) .eq.l)goto 30
      tsl(kpr) = yar(3)
      ts(kpr) = w / 2. * ( - 1. + (pbrch / pmean) ** abr(kpr)) /
     1 (bbr(kpr) * (pbase * l.e - 6) ** abr(kpr))
30    continue
      coef(kpr) (ts(kpr) + tsl(kpr) - (deltat + yar(3))) I ts(kpr)
      if(coef(kpr).gt.l.)coef(kpr) = 1.
      if(coef(kpr).lt.0.)coef(kpr) = 0.
40    if(x.ge.x2)goto 60
      sl = 0.
      s2 = 0.
      v1 = 0.
      v2 = 0.
      dsldx = 0.
      ds2dx = 0.
      y = sqrt((pd + 2. * x) ** 2 - d *d)
      theta = atan(y / d)
      al = theta / 4. * (pd + 2. *x) **2 - d/ 4. *y
      if(x.ge.xl)goto 50
      vl = 3. / 4. * (gl - 2. * x)
      vl = vl * (2. * sqr3 * d * d - pi * (pd + 2. * x) 4, 2 + 24. * al)
      sl = 2. * vl / (gl - 2. * x)
      sl = sl + 3. * (gl - 2. * x) * (pi - 6. * theta) * (pd + 2. * x)
50    yl = sqrt((gd - 2. * x) ** 2 - (5. * d - 2. * (pd + 2. * x)) ** 2)
      chi = atan(yl / (5. * d - 2. * (pd + 2. * x)))
      y2 = sqrt((pd + 2. * x) ** 2 - (3. * d - 2. * (pd + 2. * x)) ** 2)
      phi = atan(y2 / (3. * d - 2. * (pd + 2. * x)))
      a2 = phi * (pd + 2. * x) ** 2 - chi * (gd - 2. * x) ** 2
      a2 = (a2 + 2. * sqr3 * d * sqrt((3. * d - pd - 2. * x) * (3. * d
     1 - gd + 2. * x))) / 8.
      v2 = pi * (gd - 2. * x) ** 2 - 6. * sqr3 * d * d - 4. * pi * (pd
     1 + 2. * x) ** 2
      v2 = (v2 + 24. * (al + 2. * a2)) * (gl - 2. * x) / 4.
      s2 = 2. * v2 / (gl - 2. * x)
      s2 = s2 + (gl- 2. * x) * ((pi - 6. * chi) * (gd- 2. * x) + 2. *
     1 (2. * pi - 3. * phi - 3. * theta ) * (pd + 2. * x))
      vol = vl + v2
      surf = sl + s2
      frac = 1. - vol / vO
      dsdx = - surf / (x2 - x)
      dsdx = coef(kpr) * dsdxsl(kpr) + (1. - coef(kpr)) * dsdx
      dsdxsl(kpr) = dsdx
      frac = coef(kpr) * fracsl(kpr) + (1. - coef(kpr)) * frac
      fracsl(kpr) = frac
      surf = coef(kpr) * surfsl(kpr) + (1. - coef(kpr)) * surf
c     Page 89
      surfsl(kpr) = surf
      return
60    vol = 0.
      surf = 0.
      frac = fracsl(kpr) * coef(kpr) + 1. - coef(kpr)
      fracsl(kpr) = frac
      if(frac.gt..9999) frac = 1.
      if(frac.gt..9999)return
      dsdx = 0.
      dsdx = dsdxsl(kpr) * coef(kpr)
      dsdxsl(kpr) = dsdx
      if(abs(dsdx).lt.1.)dsdx = 0.
      surf = surfsl(kpr) * coef(kpr)
      surfsl(kpr) = surf
      return
c
c     zero perf calculations start here.
c
70    if(gd - 2. * x.le.0.0) go to 80
      v0 = pi * gd * gd / 4. * gl
      vol = pi * (gd - 2. * x) ** 2 / 4. * (gl - 2. * x)
      frac = 1. - vol / vO
      surf = pi / 2. * (gd- 2. * x) ** 2 + pi * (gd- 2. * x) * (gl1 2. *x)
      dsdx = - 2. * pi * (gd + gl - 6. * x)
      return
c
c     grain completely burned
c
80    surf = 0.
      frac = 1.0
      vol = 0.
      dsdx = 0.
      nsl = 1
      return
c
c     one perf calculation starts here
c
90    if(gd - pd - 4. * x.le.0.0) goto 80
      v0 = pi / 4. * (gd * gd - pd * pd) * gl
      vol -pi / 4. * ((gd - 2. * x) ** 2- (pd + 2. * x) ** 2) * (gl1 2. * x)
      frac = 1. - vol / vO
      surf = pi / 2. * ((gd - 2. *x) **2 - (pd + 2. *x) **2)
      surf = surf + pi * (gd - 2. * x) * (gl - 2. * x)
      surf = surf + pi * (pd + 2. * x) * (gl - 2. * x)
      dsdx = - 4. * pi * (gd + pd)
      return
c
c     below is the calculation for the cylindrical 19 perf grain.
c     programmer:fred robbins
c     input
c
c     p = perf diameter
c     d = grain diameter
c     gl = grain length
c     x = distance burnt
c
c     Page 90
c     output
c
c     vol - the volume of one grain at x.
c     surf = the surface area of one grain at x.
c     frac = the fraction of grain burnt at x.
c
c     w=web
c
110   p = pd
      d - gd
      w (d - 5. * p) / 6.
      pi 3.141592654
      sqrt3 = sqrt (3.)
      sqrt5 = sqrt(5.)
      sqrt6 - sqrt(6.)
      sqrtlO = sqrt(1l.)
c
c     initial volume and surface area
c
      v0 - pi! 4. *gl* (d * d 19. *p *p)
      s0 = 2. * vO / gl + pi * gl * (d + 19. * p)
c
c     xl = distance to inner sliverr burnout
c     x2 = distance to outer sliver burnout
c     dbc = distance between perforation centers
c     assumes burnout does not occur in longitudinal direction
c     wl = secondary web
c
      dbc = w + p
      wl = 0.5 * (d - p -2. * sqrt3 * dbc) xl = dbc / sqrt3 - p / 2.
      x2 = 0.25 * (dbc * (6. - sqrtlO) - 2. * p) if(x.gt.w / 2.)go to 120
c
c     not slivered yet c
      vol = pi / 4. * (gl - 2. * x) * ((d - 2. * x) ** 2 - 19. * (p + 2.
     1 * x) ** 2)
      surf = 2. * vol / (gl - 2. * x) + pi * (gl - 2. * x) * (d - 2. *
     1 x + 19. * (p + 2. * x))
      dsdx = pi * ( - 4 * d + 36 * gl - 76 * p - 216 * x)
      frac = 1. - vol / vO
      dsdxsl(kpr) = dsdx
      fracsl(kpr) = frac
      surfsl(kpr) = surf
      return
c
c     vl=total volume of inner sliver, v2=total volume of outer sliver
c     sl=total surface area of inner slivers, s2=total surface area of
c     outer slivers
c
120   vl = 0.
      v2 = 0.
      sl = 0.
      s2 = 0.
      delta = 0.
      chi = 0.
c     Page 91
      nsl = 1
      coef(kpr) = 0.
      if(igrad.eq.l.or.igrad.eq.2)go to 140
      if(nslp(kpr).eq.1)goto 130
      tsl(kpr) = yar(3)
      ts(kpr) = w / 2. * ( - 1. + (pbrch / pmean) ** abr(kpr)) /
     1 (bbr(kpr) * (pbase * l.e - 6) ** abr(kpr))
130   continue
      coef(kpr) = (ts(kpr) + tsl(kpr) - (deltat + yar(3))) / ts(kpr)
      if(coef(kpr).gt.i.)coef(kpr) = 1.
      if(coef(kpr).lt.0.)coef(kpr) = 0.
140   a3 = 0.
      if(x.ge.x2)go to 170
      theta = acos(dbc / (p + 2. * x)) al = theta /4. * (p + 2. * x) ** 2 - dbc / 4. * sqrt((p + 2. * x)
     1 ** 2 - dbc * dbc) if(x.gt.xl)go to 150
      vl = 3. * (gl - 2. * x) * (2. * sqrt3 * dbc * dbc - 4i + (p + 2.
     1 * x) *- 2 + 24 * al)
      sl = 2. * vl / (gl - 2. * x) + 12. * (gl - 2. * x) * (pi -6. *
     1 theta) * (p + 2. * x)
150   phi = acos((5. * d - 13. * p - 36. * x) / (12. * (p + 2. * x)))
      xi = acos((13. * d - 5. * p - 36. * x) / (12. * (d - 2. * x))) if(x.le.wl / 2.)go to 160
      delta = acos((2. * d - p - 6. * x) / (sqrt3 * (d - 2. * x)))
      chi = acos((d - 2. * p - 6. * x) / (sqrt3 * (p + 2. * x)))
      a3 = .125 * (chi * (p + 2. * x) ** 2 - delta * (d - 2. * x) ** 2
     1 + 2. * sqrt6 * dbc * sqrt(6. * dbc * (p + 2. * x - dbc) - (p + 2.
     1 * x) ** 2))
160   a2 = .125 * (phi * (p + 2. * x) ** 2 - xi * (d - 2. * x) ** 2 + 2.
     1 * sqrt5 * dbc * sqrt((5. * dbc - p - 2. * x) * (5. * dbc - d + 2.
     1 * x)
      v2 = .25 * (gl- 2. * x) * (pi * (d- 2. * x) ** 2 - 7. * pi *
     1 (p+2. * x) ** 2 - 24. * sqrt3 * dbc * dbc + 48. * (al + a2 + a3))
      s2 = 2. * v2/ (gl- 2. * x) + (gl- 2. * x) * ((d- 2. * x) * (pi 1 - 6. * (xi + delta)) + (p + 2. * x) * (7. * pi - 6. * (2. * theta
     1 + chi + phi)))
      vol . vl + v2
      surf = sl + s2
      dsdx = - surf / (x2 - x)
      frac 1. - vol / vO
      dsdx = coef(kpr) * dsdxsl(kpr) + (I. - coef(kpr)) * dsdx
      dsdxsl(kpr) = dsdx
      frac = coef(kpr) * fracsl(kpr) + (1. - coef(kpr)) * frac
      fracsl(kpr) = frac
      surf = coef(kpr) * surfsl(kpr) + (1. - coef(kpr)) * surf
      surfsl(kpr) = surf
      return
170   vol = 0.
      surf = 0.
      frac = fracsl(kpr) * coef(kpr) + 1. - coef(kpr)
      fracsl(kpr) = frac
      if(frac.gt..9999) frac = 1.
      if(frac.gt..9999)return
      dsdx = 0.
c     Page 92
      dsdx - dsdxsl(kpr) * coef(kpr)
      dsdxsl(kpr) = dsdx
      if(abs(dsdx).lt.l.)dsdx = 0.
      surf - surfsl(kpr) * coef(kpr)
      surfsl(kpr) - surf
      return c
c     below is the calculation for the 19 perf hex grain.
c     programmer:karen a. cieslewicz<std.cont.>
c
c     translation of the input values.
c     p= perf diameter
c     d= grain diameter
c     gl= grain length
c     x- distance burnt
c
c     translation of the output values.
c     vol= volume of one grain at x.
c     surf= surface area of one grain at x.
c     frac= mass fraction of the grain burnt at x.
c
c     assignment statement for pi.
      180 pi = 3.141592654
      sqrt3 = sqrt(3.)
      p = pd
      d = gd
c
c     d=6w + 5p is the statement for the grain diameter which will be
c     used to calculate the web.
c
c     to calculate the web.
      w - (d - 5. * p) / 6.
c
c     below is the equation to calculate the distance between the perf cen-
c     ters.
      dpc = p + w
c     to calculate the grain diameter between the flats.
      f - 2. * (sqrt3 * dpc + p / 2. + w)
c
c      to calculate the distance burnt.
      xl = dpc / sqrt3 - p / 2.
x2 = 0.125 * (5. * dpc - 4. * p)
c
c     to calculate the area.
      a - sqrt3/3. * ((w + p / 2.) ** 2) - pi / 6. * ((w + p / 2.) **2)
c
c     to calculate the initial volume of the sharp corner grain.
      vs = gl / 4. * (2. * sqrt3 * f ** 2 - 19. * pi * p ** 2)
c
c     to calculate the volume that will be removed from the grain.
      vr = 6. * a * gl
c
c     to calculate the initial volume for the 19hex grain with rounded
c     corners. 
      vo = vs - vr
c
c     to calculate the initial surface area of the sharp corner grain.
      ss = 2. *vs / gl + gl * (2. * sqrt3 * f + 19. * pi * p)
c     Page 93
c
c     to calculate the surface area that will be removed from the grain.
      sr - 12. * a + gl * (w + p / 2.) * (4. * sqrt3 - 2. * pi)
c
c     to calculate the initial surface area for the 19hex grian with rounded
c     corners.
      so - ss - sr
c
c     to calculate the unknows of the grain under the condition x.le.5*w.
      if(O.le.x.and.x.le.w / 2.) then
      a = sqrt3 / 3. * (w - 2. * x + (p + 2. * x) / 2.) ** 2 - pi /
     1 6. * (w - 2. * x + (p + 2. * x) / 2.) ** 2
c     to calculate the volume that will be removed from the sharp corner gra
      vr = 6. * a * (gl - 2. * x)
c     to calculate the volume for the sharp corner grain at some distance bu
      vn = .25 * (gl - 2. * x) * (2. * sqrt3 * (f - 2. * x) ** 2. -
     1 19. * pi* (p + 2. * x) **2.)
c
c     to calculate the volume for the 19hex grain with rounded corners.
      v = vn - vr
c
c     to calculate the surface area that will be removed from the sharp
c     corner grain.
      sr = 12. * a + (gl - 2. * x) * (w - 2 * x + (p + 2. * x) /2.)
     1 * (4. * sqrt3 - 2. * pi)
c     to calculate the surface area for the sharp corner grain.
      sn = 2. * v / (gl -2. *x) + (gl - 2. * x) * (2. * sqrt3* (f
     1 - 2. * x) + 19. *pi *(p + 2. * x))
c
c     to calculate the surface area for 19hex grain with rounded corners.
      s = sn - sr
c     to calculate the mass fraction.
      frac = 1 - v / vo
      dsdx = - 8. * sqrt3 * (f - 2. * x) - 76. * pi * (p + 2. * x)
     1 + (gl - 2. * x) * ( - 4. * sqrt3 + 38. * pi) + 16 * sqrt3 *
     1 (w + p /2. - x) - 8. *pi* (w + p / 2. -x) + (gl -2. * x)
     1 * (4. *sqrt3 - 2. * pi)
      surf = s
      vol = v
      dsdxsl(kpr) = dsdx
      fracsl(kpr) = frac
      surfsl(kpr) = surf
      return
      endif
c
c     due to the cross section at the sliver point x-.5*w there will be 24
c     identical inner slivers,12 identical side slivers, after slivering th
c     surface area and the volume function become more complex. each type
c
c     sliver will be treated separately and later the volumes will be combin
c     to complete the function.
c
c     to calculate the 12 identical side slivers for the grain x-.5/w.
      nsl - 1
      coef(kpr) = 0.
      if(igrad.eq.l.or.igrad.eq.2)go to 200
      if(nslp(kpr) .eq.l)goto 190
c     Page 94
      tsi(kpr) = yar(3)
      ts(kpr) =w I2. * .+(bc ~al *abr(kpr))/
     1 (bbr(kpr) *+pae ~ 6)rc abrep)) /*
190   continue pae le-6 abkr)
      coef(kpr) = (ts(kpr) + tsl(kpr) -(deltat + yar(3))) / ts(kpr)
      if(coef(kpr).gt.1.)coef(kpr) - 1.
      if(coef(kpr).lt.0.)coef(kpr) = 0.
200   if(w / 2.lt.x.and.x.lt.xi.and.x.lt.x2) then
c
c     to calculate the areas of the grain.
      a = sqrt3 /3. *(w - 2. * x + (p + 2. *x) / 2.) ** 2 - pi/
     1 6. * (w -2. *x + (p + 2. * x) /2.) **2
      theta = acos(dpc /(p + 2. *x))
      al =theta/4.* (p +2.*x) *2 -dpc/4. *sqrt((p +2.
     1 * x) **2 -dpc ** 2)
      omega =acos(2. *dpc / (p + 2. *x . a2=0.125 * (p +2. * x) * ((p + 2 . *x) *(omega +sin(omega
     1 )) -2. * dpc * sin(omega))
      c to calculate the volumes of the grain.
      vi = 3. * (gi - 2. * x) * (2. * sqrt3 * dpc ** 2 - pi * (p +
     1 2. *x) ** 2 + 24. * al)
      v2 = 6. * (gl - 2. * x) * (2. * dpc **2 - dpc *(p + 2. *x)
     1 - pi / 4. * (p + 2. * x) ** 2 + 2. *al + 4. *a2)
      c to calculate the surface areas of the grain.
      s1 =2. * vl / (gl-2. *x) + 12. *(gi-2. *x) *(pi-6.
     1 * theta) * (p + 2 . x)
      s2 =2. * v2 / (gl-2. *x) + 12. *(gl-2. *x) *(dpc + (p
     1 + 2. * x) * (pi / 2. - omega - theta - sin(omega)))
c     to calculate the total volume and total surface area.
      vf =v1 + v2
      sf = si + s2
c     to calculate the mass fraction.
      frac = 1. - vf /vo
      surf = sf
      dsdx = - surf /(x2 - x)
      Vol = vf
      dsdx = coef(kpr) * dsdxsl(kpr) + (1. - coef(kpr)) * dsdx
      dsdxsl(kpr) = dsdx
      frac = coef(kpr) * fracsl(kpr) + (1. - coef(kpr)) * frac
      fracsl(kpr) =frac
      surf = coef(kpr) * surfsl(kp~r) + (1. - coef(kpr)) * surf
      surfsl(kpr) = surf
      return
      endif
      if(x.gt.xl.and.x.lt.x2)then
c     to calculate the area of the grain.
      a = sqrt3 /3. *(w - 2. * x + (p + 2. *x) / 2.) ** 2 - piI
     1 6. * (w-2. x+ (p+2. * x) /2.) *2
      theta = acos(dpc /(p + 2. *x))
      al1=theta/4.*(p +2.*x)**2 -dpc/4. *sqrt((p +2.
     1 * x) **2 -dpc ** 2)
      omega =acos(2. *dpc / (p + 2. *x)-1)
      a2 = 0.125 * (p + 2. * x) *H(p + 2. *x) *(omega +
     1 sin(omega)) - 2. * dpc * sin(omega))
c     to calculate the volume of the grain.
      v2 =6. * (g1 -2. * x) * (2. * dpc *2 .dpc *(p +2. * x)
     1 - pi / 4. * (p + 2. * x) ** 2 + 2. *al + 4. *a2)
c     to calculate the surface area of the grain.
c     Page 95
      s2 = 2. * v2 / (gI - 2. * x) + 12. * (gI - 2. * x) * (dpc + (p
     1 + 2. * x) * (pi / 2. - omega - theta - sin(omega)))
c     to calculate the volume and the surface area.
      vf = v2
      sf = s2
c     to calculate the the mass fraction.
      frac = 1 - vf /vo
      surf = sf
      dsdx = - surf /(x2 - x)
      Vol =vf
      dsdx = coef(kpr) * dsdxsl(kpr) + (1. - coef(kpr)) * dsdx
      dsdxsl(kpr) = dsdx
      frac = coef(kpr) * fracsl(kpr) + (1. - coef(kpr)) * frac
      fracsl(kpr) = frac
      surf = coef(kpr) * surfsl(kpr) + (1. - coef(kpr)) * surf
      surfsl(kpr) = surf
      return
      endif
      if (x.gt .x2) then
      dsdx = 0.
      surf = 0.
      Vol =0.
      frac = fracsl(kpr) * coef(kpr) + 1. -coef(kpr)
      fracsl(kpr) = frac
      if(frac.gt. .9999) frac = 1.
      if(frac.gt. .9999)return
      dsdx = 0.
      dsdx = dsdxsl(kpr) * coef(kpr)
      dsdxsl(kpr) = dsdx
      if(abs(dsdx).lt.l.)dsdx = 0.
      surf = surfsl(kpr) * coef(kpr)
      surfsl(kpr) = surf
      return
      endif
      stop
c
c     spherical grain calculations start here
c
210   if (gd .le. 2.*x) goto 80
      vol =pi /6. * (gd -2. * x) ** 3
      surf pi *(gd - 2. *x) ** 2
      frac ((gd - 2. * x) /gd) ** 3
      dsdx =pi * 4. * (2. *x - gd)
      return
c
220   format (lx,'unacceptable granulation')
      end
      subroutine jint(btdia,btlen,x,y,nchpts,chdist,chdiam,bint,bvol)
      dimension bint(l0),chdist(6),chdiam(6)
      pi=3. 141593
      areaa-pi*btdia*btdia/4.
      distbp-x-bt len
      point s-10.
      step=y/points
      zz=0.
      bvoll-O.
      i chg= 0
      do 1 j=1,10
c     Page 96
      bint (j) =0.
      2. continue
      if(step.1t. 1.e-1O)then
      bint (7) =1./ (pi* (chdiam(l) **2/4.)) **2
      return
      endif
      nt sw=0
      do 2 j=2,nchpts
      if(y.gt.chdist(j)) go to 2
      nchp=j
      i chg= 1
      holddt-chdist (j)
      holddm=chdiam (j)
      diam- (y-chdist (j-1) )/ (chdist (j) -chdist (j-1))
      chdiam(j)=chdiam(j-l)+diam*(chdiaxn(j)-~chdiam(j-1))
      chdist (j)=y
      go to 3
      2 continue
      nchp=nchpts+l
      chdist (nchp)=y
      chdiam(nchp) =chdiam(nchpts)
c     write(6,*)chdist(nchp-1) ,chdist(nchp)
c     910 format(lx,2e11.4)
3     continue
      areal=chdiam(1)**2 * pi / 4.
      bint5o=0 .0
      do 58 il=l,nchp-1
      npt=int ((chdist (11+2)-chdist (ii)) /step)
      if (npt.le.0)npt=l
c     write(6,*)npt,chdist(il),chdist(il+l),,step
      stepl= (chdist (il+1) -chdist (ii)) /npt
c     write(6,912)chdiam(nchp-1),chdiam(nchp)
c     912 format(lx,2el1.4)
      rl=chdiam(il) *5
c     areal=pi*rl*rl
      do 57 i~1,npt
      zz~zz+stepl
      if P(zz .gt .chdist (nchp) )zz=chdist (nchp)
      diam= (zz-chdist (ii))! (chdist (il+1ý)-chdist (ii))
      diam=chdiam(i1)+cliam* (chciiam(iI+1)-chdiam(il))
      r2=C.5*diam
      area2=pi*r2*r2
      bvol2=bvoll+stepl*pi/3. *(rl*rl+rl*r2+r2*r2)
      if (zz .gt.distbp)then
      if (intsw.eq.0)then
      diam=(distbp-chdist(ii))/(chdi'st(i14-1)-chdiz-t(ilflI
      diarn=chdiam(i1)+diam*(chdiam(i1+l)-chdiam(i1))
      r2a=0.5*diam
      area2=pi*r2a*r2a
      stepla=stepl- (zz-distbp)
      bvol2=bvoll+stepla*pi/3.* (rl*rl+rl*r2a+r2a*r2a)
      bintl=bint (1)+0. 5*stepla* (bvoll/areal+bvol2/area2)
      bint (2) =bvol2*bvol2/area2/area2
      bint (3) =bint (3) +0. 5*stepla* (bint (1)*areal+bintl*area2)
      bint(4)=bint(4)+.5*stepla*(bvoll*bvoll/areal+bvol2*bvol2/area2)
c     bint(5)=bi-nt(5)+.5*stepla*(1./areal+1./area2)
      bint (6)=bvol2/area2/area2
      bint (7)"1 ./area2/area2
c     Page 97
      bint(1) =bintl
      bvoll=bvol2
      areal-area2-areaa
      area2=pi *r2 *r2
      area2-area2-areaa
      stepla=zz-distbp
      rl=r2a
      bvol2=bvoll+stepla*pi/3.* (rl*rl+rl*r2+r2*r2)
      bvo12=bvol2-stepla*areaa
      bintlO=bint (10)+0.5*stepla* (bvoll/areal+bvol2/area2)
      bint5-barnt(5)+0.5*stepla*(1./areal +1./area2)
      bint (2) =bvol2*bvol2/area2/area2
c     bint(3)=bint(3)+O.5*stepla*(bint(l)*areal+bintl*area2)
      bint(4)=bint(4)+.5*stepla*(bvoll*bvoll/areal+bvol2*bvol2/area2)
c     bint(5)=bint(5)+.5*stepla*(l./areal+1./area2)
      bint (6) =bvol2/area2/area2
      bint (7) =1./area2/area2
c     bint5a=bint5o+.5*stepla*(1./areal + l./area2)
      bint(8)=bint(8)+.5*stepla*(areal*bint(5) + bint5*area2)
      bint(9)=bint(9)+.5*stepla*(areal*bint(10)+area2*bintlO)
      bint (1O)=bintl0
      bint (l)=bintl
      bint (5)=bint5
c     bint5o=bint5a
      areal~area2
      bvoll=bvol2
      rl=r2
      intsw=1
      go to 57
      else
      area2=area2-areaa
      bvol2=bvol2-stepl *areaa
      bintl0=bint(10)+0.5*stepl*(bvoll/areal+bvol2/area2)
      bint (2) =bvol2*bvol2/area2/larea2
c     bint(3)=bint(3)+O.5*stepl*(bint(l)*areal+bintl*area2)
      bint (4) =bint (4) +0. 5*stepl* (bvoll*bvoll/areal+bvol2*bvol2/area2)
c     bint(5)=bint(5)+.5*stepl*(l./areal+1,./area2)
      bint (6) =bvol2/area2/area2
      bint (7) =1./area2/area2
      bint5=bint(5)+.5*stepl*(1./areal + 1./area2)
      bint(8)=bint(8)+.5*stepl*(areal*bint(5) + bint5*area2)
      bint (9) =bint (9) +. 5*stepl* (areal*bint (10) +area2*bintl0)
      bint (10) =bint 10
c     bint(l)=bintl
      bint (5) =bint5
c     bint5o=birit5a
      areal=area2
      rl=r2
      bvoll=bvol2
      go to 57
      endif
      endif
      bintl=bint (1)+0. 5*stepl* (bvoll/areal+bvol2/area2)
      bint (2) =bvol2*bvol2/area2/area2
      bint (3) =bint (3) +0. 5*stepl* (bint (1)*areal+bintl*area2)
      bint (4) =bint (4) +0. 5*stepl* (bvoll*bvoll/areal+bvol2*bvo12/area2)
c     bint(5)=bint(5)+.5*stepl*(l./areal+l./area2)
      bint (6) =bvol2/area2/area2
c     Page 98
      bint (7) 1./area2/area2
      bint (1) =bintl
      areal-~area2
      rl-sr2
      bvol1=bvol2
57    continue
58    continue
      bvol=bvol2
      if (ichg.eq.1) then
      chdiam(nchp) =holddm
      chdist (nchp) =holddt
      endif
c     write (6, 915 )
c     915 forrnat('1',1x,'Leaving Jint')
      return
end
