      program ibrgac
      character bdfile*10,outfil*10
      dimension br(10),trav(10),rp(10),tr(10),forcp(10),tempp(10),covp(
     &10
      dimension chwp(10),rhop(10),gamap(10),nperfs(10),glenp(10),pdpi(10
     &,p dpo(10),gdiap(10),dbpcp(10),alpha(10,10),beta(10,10),
     &pres(10,10)
      dimension a(4),b(4),ak(4),d(20),y(20),p(20),z(20),frac(10),surf(10
     &),nbr(10),ibo(10)
      real lambda,j1zp,j2zp,j3zp,j4zp
      dimension chdist(5),chdiam(5),bint(4)
c     call gettim(ihr,imin,isec,ihuns)
      pi=3.14159
      write(*,15)
15    format(' input name of data file to be used as input ')
      read(*,10)bdfile
10    format(a10)
      open(unit=2,err=999,file=bdfile,status='old',iostat=ios)
      rewind 2
      write(*,25)
25    format(' input name of output file ')
      read(*,10)outfil
      open(unit=3,err=998,file=outfil,status='new')
      write(3,16)bdfile
16    format(' THE INPUT FILE IS ',a10)
      read(2,*,end=20,err=30)cham,grve,aland,glr,twst,travp,igrad
      if(igrad.gt.1)go to 51
      write(3,55)
55    format(1x,'using Lagrange pressure gradient')
      go to 52
c     define chambrage assumes nchpts=number of points to define
c     chamber > or = 2 < or = 5 (?),chdiam(I) defines chamber diameter 
c     at chdist (I) chamber distance. chidiam(nchpts) is assumed to be
c     the bore diameter and chdist(i) is assumed to be -, i.e. at the
c     breech. Assumes truncated cones.
51    write(3,47,err=30)
47    format(1x,'Using chambrage pressure gradient')
      read(2,*,end=20,err=30)nchpts,(chdist(I),chdiam(I),I=1,nchpts)
      write(3,53,err=30)(chdist(I),chdiam(I),I=1,nchpts)
53    format(///,'   chamber distance cm   chamber diameter cm',/(5x,e14
     &.6,5x,e14.6))
      do 54 I=1,nchpts
      chdist(I)=0.01*chdist(I)
54    chdiam(I)=0.01*chdist(I)
c     calculate chamber integrals and volume
      if(nchpts.gt.5) write(3,44,err=30)
44    format(1x,'use first 5 points')
      if(nchpts.gt.5)nchpts=5
      bore=chdiam(nchpts)
      if(chdist(2).ne.0.0)write(3,45,err=30)
45    format(1x,' # points ? ')
c     page 15
      chdist(1)=0.0
      pi3=pi/3.0
      b1=0.0
      b2=0.0
      b3=0.0
      b4=0.0
      points=25.0
56    points=points+points
      step=chdist(nchpts)/points
      zz=0.0
      bint(1)=0.0
      bint(3)=0.0
      bint(4)=0.0
      bvol=0.0
      r2=0.5*chdiam(1)
      k=1
      j=int(points+0.5)
      do 57 I=1,j	!is that right? 1 not l?
      zz=zz+step
      if(k.eq.nchpts-1)go to 46
      do 58 I1=k,nchpts-1
      if (zz.gt.chdist(I1).and. zz.lt.chdist(I1+1)go to 59
58    continue
      I1=nchpts-1
      k=I1
46    diam=(zz-chdist(k))/(chdist(k+1)-chdist(k))
      diam=chdiam(k)+diam*(chdiam(k+1)-chdiam(k))
      r1=0.5*diam
      area=pi*(r1+r2)*(r1+r3)/4.
      bvol=bvol+step*pi3*(r1*r1+r1*r2+r2*r2)
      bint(1)=bint(1)+step*bvol/area
      bint(3)=bint(3)+step*area*bint(1)
      bint(4)=bint(4)+step*bvol*bvol/area
57    r2=r1
      temp=abs(1.0-b1/bint(1))
      if(abs(1.0-b3/bint(3)).gt.temp)temp=abs(1.0-b3/bint(3))
      if(abs(1.0-b4.bint(4)).gt.temp)temp=abs(1.0-b4/bint(4))
      if(temp.le.0.001)go to 41
      b1=bint(1)
      b3=bint(3)
      b4=bint(4)
      go to 56
41    cham=bvol*1.e6
c     write(3,47,err=30)bint(1),bint(3),bint(4)
c     format(1x,'bint 1 = ',e14.6,' bint 3 = ',e14.6,' bint 4 = ',e14.
c    &6)
      chmlen=chdist(nchpts)
52    write(3,40,err=30)cham,grve,aland,glr,twst,travp
40    format(1x,'chamber volume cm**3',e.14.6,/' groove diam cm',14.6,/
     &' land diam cm',e14.6,/' groove/land ratio',e14.6,/' twist turns
     &/caliber ',e14.6,/' projectile travel cm', e14.6///)
      cham-cham*1.e-6
      grve=grve*1.e-2
      aland=aland*1.e-2
      travp=travp*1.e-2
c     page 16
      read(2,*,end=20,err=30)prwt,iair,htfr,pgas
      write(3,50,err=30)prwt,iair,htfr,pgas
50    format(1x,'projectile mass kg',e14.6./' switch to calculate energ
     &y lost to air resistance J',i2,/' fraction of work against bore u
     &sed to heat the tube',e14.6/1x,' gas pressure Pa' ,e14.6)
     
      read(2,*,end=20,err=30)npts,(br(i),trav(i),i=1,npts) !i=1 or i=l?
      write(3,60,err=30)npts,(br(i),trav(i),i=1,npts)
60    format(1x,'number barrel resistance points',i2,/' bore resistance
     & MPa - travel cm'/(1x,e14.6,e14.6))
      write(3,65)
      do 62 i=1,npts
      br(i)=br(i)*1.e6
      trav(i)=trav(i)*1.e-2
62    continue
65    format(1x)
      read(2,*,end=20,err=30)rcwt,nrp,(rp(i),tr(i),i=1,nrp)
      write(3,70,err=30)rcwt,nrp,(rp(i),tr(i),i=1,nrp)
70    format(1x,' mass of recoiling parts kg',e14.6,/' number of recoi
     &l point pairs',i2,/' recoil force N',' recoil time sec'/,(1x,e14
     &.6,3x,e14.6))
      write(3,65)
      read(2,*,end=20,err=30)ho,tsh,cshl,twal,hl,rhocs
      write(3,75,err=30)ho,tshl,cshl,twal,hl,rhocs
75    format(1x,' free convective heat transfer coefficient wcm**2 K',
     &e14.6,/' chamber wall thickness cm',e14.6,/' heat capacity of st
     &eel of chamber wall J/g K',e14.6,/' initial temperature of chambe
     &r wall K',e14.6,/' heat loss coefficient',e14.6,/' density of ch
     &amber wall steel g/cm**3',e14.6//)
      ho=ho/1.e-4
      tshl=tshl*1.e-2
      cshl=cshl*1.e+3
      rhocs=rhocs*1.e-3/1.e-6
      read(2,*,end=20,err=30)forcig,covi,tempi,chwi,gamai
      write(3,85,err=30)forcig,covi,tempi,chwi,gamai
85    format(1x,' impetus of igniter propellant J/g',e14.6,/' covolume
     & of igniter cm**3/g',e14.6,/' adiabatic flame temperature of igni 
     &ter propellant K',e14.6,/' initial mass of igniter kg',e14.6,/' r
     &atio of specific heats for igniter',e14.6//)
      forcig=forcig*1.e+3
      covi=covi*1.e-6/1.e-3
      read(2,*,end=20,err=30)nprop,(forcp(i),tempp(i),covp(i),chwp(i),
     &rhop(i),gamap(i),nperfs(i)glenp(i),pdpi(i),pdpo(i),gdiap(i),dbpc
     &(i),i=1,nprop)
      write(3,95,err=30)(i,forcp(i),tempp(i),covp(i),chwp(i),
     &,rhop(i),gamap(i),nperfs(i),glenp(i),pdpi(i),pdpo(i),gdiap(i),dbpc
     &p(i),I=1,nprop)
95    format((' for propellant number',i2,/' impetus of propellant J/g
     &',e14.6,/' adiabatic temperature of propellant K',e14.6,/' covol
     &ume of propellant cm**3/g',e14.6/' initial mass of propellant kg'
     &,e14.6,/' density of propellant g/cm**#',e14.6' ratio of specifi
     &c heats for propellant',e14.6/' number of perforations of propell
     &ant',i2/' length of propellant grain cm',e14.6/' diameteter of inn
     &er perforation in propellant grains cm',e14.6/' diameter of outer
     &perforation of propellant grains cm',e14.6/' outside diameter of
c     page 17
