    program ibrgac
      character bdfile*10,outfil*10
      dimension br(10),trav(10),rp(10),tr(10),forcp(10),tempp(10),covp(&
     &10)
      dimension chwp(10),rhop(10),gamap(10),nperfs(10),glenp(10),pdpi(10&
     &,p dpo(10),gdiap(10),dbpcp(10),alpha(10,10),beta(10,10),&
     &pres(10,10)
      dimension a(4),b(4),ak(4),d(20),y(20),p(20),z(20),frac(10),surf(10&
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
53    format(///,'   chamber distance cm   chamber diameter cm',/(5x,e14&
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
c     format(1x,'bint 1 = ',e14.6,' bint 3 = ',e14.6,' bint 4 = ',e14.&
c    &6)
      chmlen=chdist(nchpts)
52    write(3,40,err=30)cham,grve,aland,glr,twst,travp
40    format(1x,'chamber volume cm**3',e.14.6,/' groove diam cm',14.6,/&
     &' land diam cm',e14.6,/' groove/land ratio',e14.6,/' twist turns&
     &/caliber ',e14.6,/' projectile travel cm', e14.6///)
      cham-cham*1.e-6
      grve=grve*1.e-2
      aland=aland*1.e-2
      travp=travp*1.e-2
c     page 16
      read(2,*,end=20,err=30)prwt,iair,htfr,pgas
      write(3,50,err=30)prwt,iair,htfr,pgas
50    format(1x,'projectile mass kg',e14.6./' switch to calculate energ&
     &y lost to air resistance J',i2,/' fraction of work against bore u&
     &sed to heat the tube',e14.6/1x,' gas pressure Pa' ,e14.6)
     
      read(2,*,end=20,err=30)npts,(br(i),trav(i),i=1,npts) !i=1 or i=l?
      write(3,60,err=30)npts,(br(i),trav(i),i=1,npts)
60    format(1x,'number barrel resistance points',i2,/' bore resistance&
     & MPa - travel cm'/(1x,e14.6,e14.6))
      write(3,65)
      do 62 i=1,npts
      br(i)=br(i)*1.e6
      trav(i)=trav(i)*1.e-2
62    continue
65    format(1x)
      read(2,*,end=20,err=30)rcwt,nrp,(rp(i),tr(i),i=1,nrp)
      write(3,70,err=30)rcwt,nrp,(rp(i),tr(i),i=1,nrp)
70    format(1x,' mass of recoiling parts kg',e14.6,/' number of recoi&
     &l point pairs',i2,/' recoil force N',' recoil time sec'/,(1x,e14&
     &.6,3x,e14.6))
      write(3,65)
      read(2,*,end=20,err=30)ho,tsh,cshl,twal,hl,rhocs
      write(3,75,err=30)ho,tshl,cshl,twal,hl,rhocs
75    format(1x,' free convective heat transfer coefficient wcm**2 K',&
     &e14.6,/' chamber wall thickness cm',e14.6,/' heat capacity of st&
     &eel of chamber wall J/g K',e14.6,/' initial temperature of chambe&
     &r wall K',e14.6,/' heat loss coefficient',e14.6,/' density of ch&
     &amber wall steel g/cm**3',e14.6//)
      ho=ho/1.e-4
      tshl=tshl*1.e-2
      cshl=cshl*1.e+3
      rhocs=rhocs*1.e-3/1.e-6
      read(2,*,end=20,err=30)forcig,covi,tempi,chwi,gamai
      write(3,85,err=30)forcig,covi,tempi,chwi,gamai
85    format(1x,' impetus of igniter propellant J/g',e14.6,/' covolume&
     & of igniter cm**3/g',e14.6,/' adiabatic flame temperature of igni& 
     &ter propellant K',e14.6,/' initial mass of igniter kg',e14.6,/' r&
     &atio of specific heats for igniter',e14.6//)
      forcig=forcig*1.e+3
      covi=covi*1.e-6/1.e-3
      read(2,*,end=20,err=30)nprop,(forcp(i),tempp(i),covp(i),chwp(i),&
     &rhop(i),gamap(i),nperfs(i)glenp(i),pdpi(i),pdpo(i),gdiap(i),dbpc&
     &(i),i=1,nprop)
      write(3,95,err=30)(i,forcp(i),tempp(i),covp(i),chwp(i),&
     &,rhop(i),gamap(i),nperfs(i),glenp(i),pdpi(i),pdpo(i),gdiap(i),dbpc&
     &p(i),I=1,nprop)
95    format((' for propellant number',i2,/' impetus of propellant J/g&
     &',e14.6,/' adiabatic temperature of propellant K',e14.6,/' covol&
     &ume of propellant cm**3/g',e14.6/' initial mass of propellant kg'&
     &,e14.6,/' density of propellant g/cm**#',e14.6' ratio of specifi&
     &c heats for propellant',e14.6/' number of perforations of propell&
     &ant',i2/' length of propellant grain cm',e14.6/' diameteter of inn&
     &er perforation in propellant grains cm',e14.6/' diameter of outer&
     &perforation of propellant grains cm',e14.6/' outside diameter of&
c     page 17
     &propellant graincm',e14.6/' distance between perf centers cm',e1&
     &4.6)//)
      tmpi=0.0
      do 96 i=1,nprop
      forcp(i)=forcp(i)*1.e+3
      covp(i)=covp(i)*1.e-6/1.e-3
      rhop(i)=rhop(i)*1.e-3/1.e-6
      glenp(i)=glenp(i)*0.01
      pdpi(i)=pdpi(i)*0.01
      pdpo(i)=pdpo(i)*0.01
      gdiap(i)=gpiap(i)*0.01
      dbpcp(i)=dbpcp(i)*0.01
      tmpi=tmpi+chwp(i)
96    continue
      tmpi=tmpi+chwi
      do 97 j=1,nprop
      read(2,*,end=20,err=30)nbr(j),(alpha(j,i),beta(j,i),pres(j,i),&
     &i=1,nbr(j))
110   format(1x,'number of burning rate points',i2/3x,' exponent',8x,'&
     & coefficient',10x,' pressure'/5x,'-',15x,'cm/sec-MPa**ai',10x,'MP&
     &a',/(1x,e14.6,5x,e14.6,15x,14.6))
      do 112 i=1,nbr(j)
      beta(j,i)=beta(j,i)*1.e-2
      pres(j,i)=pres(j,i)*1.e6
112   continue
97    continue
      write(3,65)
      read(2,*,end=20,err=30)deltat,deltap,tstop
      write(3,120,err=30)deltat,deltap,tstop
120   format(1x,'time increment msec',14.6' print increment msec',e14&
     &.6/1x,'time to stop calculation msec ',e14.6)
      write(*,130)
      deltat=deltat*0.001
      deltap=deltap*0.001
      tstop=stop*.001
130   format(1x,'the data has been read')
      if(igrad.gt.1)go to 131
      bore=(glr*grve*grve+aland*aland)/(glr+1.)
      bore=sqrt(bore)
131   areab=pi*bore*bore/4.
      lambda=1./((13.2+4.*log10(100.*bore))**2)
      pmaxm=0.0
      pmaxbr=0.0
      pmaxba=0.0
      tpmaxm=0.0
      tpmaxbr=0.0
      tpmaxba=0.0
      tpmax=0.0
      a(1)=0.5
      a(2)=1.-sqrt(2.)/2.
      a(3)=1.+sqrt(2.)/2.
      a(4)=1./6.
      b(1)=2.
c     Page 18
      b(2)=1.
      b(3)=1.
      b(4)=2.
      ak(1)=0.5
      ak(2)=a(2)
      ak(3)=a(3)
      ak(4)=0.5
      vp0=0.0
      tr0=0.0
      tcw=0.0
      do 5 i=1,nprop
      ibo(i)=0
5     vp0=chwp(i)/rhop(i)+vp0
      volgi=cham-vp0-chwi*covi
      pmean=forcig*chwi/volgi
      volg=volgi
      volgi=volgi+vp0
      wallt=twal
      ptime=0.0
      ibrp=8
      z(3)=1.
      nde=ibrp+nprop
      write#,132)areab,pmean,vp0,volgi
132   format(1x,'area bore m^2 ',e16.6,' pressure from ign Pa',e16.6,/&
     &1x,' volume of unburnt prop m^3 ',e16.6,' init cham vol-cov ign m&
     &^3 ',16.6)
      write(3,6)
6     format(1x,'      time       acc      vel      dis      mpress      &
     &       pbase       pbrch       ')
      iswl=0
19    continue
      do l1 J=1,4
c     FIND BARREL RESISTANCE
      do 201 k=2, npts
      if(y(2)+y(7).ge.trav(k))go to 201
      go to 203
201   continue
      k=npts
203   resp=(trav(k))-y(2)-y(7))/trav(k)-trav(k-1))
      resp=br(k)-resp*(br(k)-br(k-1))
c     FIND MASS FRACTION BURNING RATE
      do 211 k=1, nprop
      if(ibo(k).eq.1)goto211
      call prf017(pdpo(k),pdpi(k),gdiap(k),dbpcp(k),glenp(k),surf(k),fra&
     &c(k),y(ibrp+k),nperfs(k),u)
      if(surf(k).lt.1.e-10) ibo(k)=1
211   continue
      k=nprop
c     ENERGY LOSS TO PROJECTILE TRANSLATION
      elpt=prwt*y(1)*y*(1)/2.
c     ENERGY LOSS DUE TO PROJECTILE ROTATION
      if(igrad.le.1)go to 214
      pt=y(2)+y(7)
c     Page 19
      vzp=bvol+areab*pt
      j4zp=bint(4)+((bvol+areab*pt)**3-bvol**3)/3./areab/areab
      elgpm=tmpi*y(1)*y(1)*areab*areab*j4zp/2./vzp/vzp/vzp
      go to 216
214   elgpm=tmpi*y(1)*y(1)/6.
c     ENERGY LOSS FROM BORE RESISTANCE
216   elbr=y(4)
      z(4)=areab*resp*y(1)
c     ENERGY LOSS DUE TO RECOIL
      elrc=rcwt*y(6)*y(6)/2
c     ENERGY LOSS DUE TO HEAT LOSS
      areaw=cham/areab*pi*bore+2.*areab+pi*bore*(y(2)+y(7))
      avden=0.0
      avc=0.0
      avcp=0.0
      z18=0
      z19=0
      do 213 k=1,nprop
      z18=forcp(k)*gamap(k)*chwp(k)*frac(k)/(gamap(k)-1.)/tempp(k)+z18
      z19=chwp(k)*frac(k)+z19
      avden=avden+chwp(k)*frac(k)
      continue
      avcp=(z18+forcig*gamai*chwi/(gamai-1.)/tempi)/(z19+chwi)
      avden=(avden+chwi)/(volg+cov1)
      avvel=.5*y(1)
      htns=lambda*avcp*avden*avvel+ho
      z(5)=areaw*htns*(tgas=wallt)*hl
      elht=y(5)
      wallt=(elht+htfr*elbr)/cshl/rhocs/areaw/tshl+twal
c     write(3,*)lambda,avcp,avden,avvel,ho,areaw,htns,tgas,wallt,hl,z(5)&
c    &,elht
c     ENERGY LOSS DUE TO AIR RESISTANCE
      air=iair
      z(8)=y(1)*pgas*air
      elar=areab*y(8)
c     RECOIL
      z(6)=0.0
      if(pbrch.le/rp(1)/areab)go to 221
      rfor=rp(2)
      if(y(3)-tr0.ge.tr(2))go to 222
      rfor=(tr(2)-(y(3)-tr0))/(tr(2)-tr(1))
      rfor=rp(2)-rfor*(rp(2)-rp(1))
222   z(6)=areab/rcwt*(pbrch-rfor/areab-resp)
      if(y(6).lt.0.0)y(6)=0.0
      z(7)=y(6)
      goto 223
221   tr0=y(3)
223   continue
c     CALCULATE GAS TEMPERATURE
      eprop=0.0
      rprop=0.0
      do 231 k=1,nprop
      eprop=eprop+forcp(k)*chwp(k)*frac(k)/(gamap(k)-1.)
      rprop=rprop+forcp(k)*chwp(k)*frac(k)/(gamap(k-1.)/tempp(k)
231   continue
c     Page 20
      tenergy=elpt+elpr+elgpm+elbr+elrc+elht+elar
      tgas=(eprop+forcig*chwi/(gamai-1.)-elpt-elpr-elgpm-elbr-elrc-elht&
     &-elar)/(rprop+forcig*chwi/(gamai-1.)/tempi)
c     FIND FREE VOLUME
      v1=0.0
      cov1=0.0
      do 241 k=1, nprop
      v1=chwp(k)*(1.-frac(k))/rhop(k)+v1
      cov1=cov1+chwp(k)*covp(k)*frac(k)
241   continue
      volg=volgi+areab*(y(2)+y(7))-v1-cov1
c     CALCULATE MEAN PRESSURE
      r1=0.0
      do 251 k=1,nprop
      r1=r1+forcp(k_*chwp(k)*frac(k)/tempp(k)
251   continue
      pmean=tgas/volg*(r1+forcig*chwi/tempi)
      resp=resp+pgas*air
      if(igrad.le.1)go to 252
      if(iswl.ne.0)go to 253
      pbase=pmean
      pbrch=pmean
      if(pbase.gt.resp=1.)iswl=1
      go to 257
c     USE CHAMBRAGE PRESSURE GRADIENT EQUATION
253   j1zp=bint(1)+(bvol*pt+areab/2.*pt*pt)/areab
      j2zp=(bvol+areab*pt)**2/areab/areab
      j3zp=bint(3)+areab*bint(1)*pt+bvol*pt*pt/2,+areab*pt*pt*pt/6.
      a2t=-tmpi*areab*areab/prwt/vzp/vzp
      alf=1.-a2t*j1zp
      alt=tmpi*areab*(areab*y(1)*y(1)/vzp+areab*resp/prwt)/vzp/vzp
      bt=-tmpi*y(1)*y(1)*areab*areab/2./vzp/vzp/vzp
      bata=-alt*jlzp-bt*j2zp
      gamma=alf+a2t*j3zp/vzp
      delta=bata+alt*j3zp/vzp+bt*j4zp/vzp
c     calculate base pressure
      pbase=(pmean-delta)/gamma
c     calculate breech pressure
      pbrch=alf*pbase+bata
      go to 254
c     USE LAGRANGE PRESSURE GRADIENT EQUATION
252   if(isw1.ne.0)go to 256
      if(pmean.lt.resp)resp=pmean
c     CALCULATE BASE PRESSURE
256   pbase=(pmean+tmpi*resp/3./prwt)/(1.+tmpi/3./prwt)  
      if(pbase.gt.resp+1.)isw1=1
c     CALCULATE BREECH PRESSURE
      pbrch=pbase+tmpi/2./prwt*(pbase-resp)
c     CALCULATE PROJECTILE ACCELERATION
254   z(1)=areab*(pbase-resp)/prwt
      if(z(1).lt.0.0)go to 257
      go to 258
257   if(isw1.eq.0)z(1)=0.0
258   if(y(1).lt.0.0)y(1)=0.0
      z(2)=y(1)
c     Page 21
c     GET BURNING RATE
      do 264 m=1,nprop
      z(ibrp+m)=0.0
      if(ibo(m).eq.1) got 264
      do 262 k=1,nbr(m)
      if(pmean.gt.pres(m,k))go to 262
      go to 263
262   continue
      k=nbr(m)
263   z(ibrp+m)=beta(m,k)*(pmean*1.e-6)**alpha(m,k)
264   continue
      do 21 i=1,nde
      d(i)=(z(i)-b(j)*p(i))*a(j)
      y(i)=deltat*d(i)+y(i)
      p(i)=3.*d(i)-ak(j)*z(i)+p(i)
21    continue
11    continue
      t=t+deltat
      if(pmaxm.gt.pmean)go to 281
      pmaxm=pmean
      tpmaxm=y(3)
281   if(pmaxba.gt.pbase)go to 282
      pmaxba=pbase
      tpmaxba=y(3)
282   if(pmaxbr.gt.pbrch)go to 283
      pmaxbr=pbrch
      tpmaxbr=y(3)
283   continue
      if(y(3).lt.ptime)go to 272
      ptime=ptime+deltap
      write(3,7)y(3),z(1),y(1),y(2),pmean,pbase,pbrch
7     format(1x,7e11.4)
316   format(1x)
272   continue
      if(t.gt.tstop)goto 200
      if(y(2).gt.travp)go to 200
      rmvelo=y(1)
      tmvelo=y(3)
      disto=y(2)
      go to 19
200   write(3,311)t,y(3)
311   format(1x,' deltat t', e14.6,' intg t',e14.6)
      write(3,312)pmaxm,tpmaxm
312   format(1x,'PMAXMEAN P ',e14.6,' time at PMAXMEAN sec ',e14.6)
      write(3,313)pmaxba,tpmaxba
313   format(1x,'PMAXBASE Pa ',e14.6,' time at PMAXBASE sec ',e14.6)
      write(3,314)pmaxbr,tpmaxbr
314   format(1x,'PMAXBREECH Pa ',e14.6,' time at PMAXBREECH sec ',e14.6)
      if(y(2).le.travp)go to 303
      dfract=(travp=disto)/(y(2)-disto)
      rmvel=(y(1)-rmvelo)*dfract+rmvelo
      tmvel=(y(3)-tmvelo)*drafct+tmvelo
      write(3,318)rmvel,tmvel
318   format(1x,'muzzle velocity m/s ',e14.6,' time of muzzle velocity s&
     &ec ',e14.6)
c     Page 22
      goto 319
303   write(3,327)y(1),y(3)
327   format(1x,'veloicty of projectile m/s ',e14.6,' at this time msec&
     &',e14.6)
319   efi=chwi*forcig/(gamai-1.)
      efp=0.0
      do 315 i=1,nprop
      efp=efp+chwp(i)*forcp(i)/(gamap(i)-1.0)
315   continue
      tenerg=efi+efp
      write(3,317)tenerg
317   format(1x,'total initial energy available J = ',e14.6)
      tengas=chwi*forcig*tgas/(gamai-1.)/tempi
      do 135 i=1,nprop
      tengas=(frac(i)*chwp(i)*forcp(i)*tgas/tempp(i)/(gamap(i)-1.))+teng&
     &as
      write(3,328)i,frac(i)
328   format(' FOR PROPELLANT ',I2,' MASSFRACT BURNT IS ',e14.6)
135   continue
      write(3,136)tengas
136   format(1x,'total energy remaining in gas J= ',e14.6)
      write(3,320)elpt
320   format(1x,'energy loss from projectiletranslation J= ',e14.6)
      write(3,321)elpr
321   format(1x,'energyloss from projectile roration J= ',e14.6)
      write(3,322)elgpm
322   format(1x,'energy lost to gas and propellant motion J= ',e14.6)
      write(3,323)elbr
323   format(1x,'energy lost to bore resistance J ',e14.6)
      write(3,324)elrc
324   format(1x,'energy lost to recoil J= ',e14.6)
      write(3,325)elht
325   format(1x,'energy loss from heat transfer J= ',e14.6)
      write(3,326)elar
326   format(1x,'energy lost to air resistance J= ',e14.6)
c     call gettim(ihro,imino,iseco,ihunso)
c     time=(ihro-ihr)*3600.+(imino-imin)*60.+(iseco-isec)+(ihunso-inhuns)&
c    &/100.
c     write(3,*)time
      stop
20    write(*,140)
140   format(1x,'end of file encounter')
      stop
30    write(*,150)
999   continue
998   continue
150   format(1x,'read or write error')
      stop
      end
      SUBROUTINE PRF017(P,P1,D,D1,L,SURF,MASSF,X,NP,u)
      IMPLICIT REAL*4(A-Z)
C     
C     P=OUTER PERF DIA
C     P1=INNER PERF DIA
C     D=OUTER DIA
c     Page 23
C     D1=DISTANCE BETWEEN PERF CENTRES
C     L=GRAIN LENGTH
C     NP=NUMBER OF PERFS
C
C     SURF=OUTPUT SURFACE AREA
C     MASSF=OUTPUT MASS FRACTION OF PROPELLANT BURNER
C
C     W=WEB BETWEEN OUTER PERFS
C     W0=OUTER WEB
C     W1=WEB BETWEEN OUTER AND INNER PERFS
C     W4=MINIMUM WEB
      INTEGER ITYM,NP
      DATA PI,SQRT3/3.14159,1.732051/,ITYM/0/
      DATA HAFPI,PIFOR,TWOPI/1.570796,.785398,6.283185/
C
      IF(ITYM.GT.0)GO TO 10
      P1SQ=P1*P1
      DISQ=D1*D1
      PSQ=P*P
      DSQ=D*D
      D1SQ3=DI*SQRT3
      D2SQ3=D1SQ*SQRT3
      IF(NP.EQ.0)GO TO 2000
      IF(NP.EQ.1)GO TO 3000
      IF(NP.NE.7)GO TO 60
      IF(P1.GT.(P+D1*(SQRT3-1))) GO TO 60
      IF(D.GE.D1*(SQRT3+1.)-P)GO TO 130
60    WRITE(6,90)
90    FORMAT(1X,'UNACCEPTABLE GRANULATION')
      STOP
130   W=D1-P
      IF(W.LT.0)GO TO 60
      W0=(D-P-2.*D1)/2.
      IF(W.LT.0.)GO TO 60
      W1=(2.*D1-P-P1)/2.
      IF(W1.LT.0.)GO TO 60
      X1=(P1SQ-PSQ+4.*D1SQ-2.*P1*D1SQ3)/4./(D1SQ3+P-P1)
      X2=(4.*D1SQ+D*D-2.*D*D1SQ3-PSQ)/4./(-D1SQ3+P+D)
      A=PI*L*(D+P1+6.*P)+HAFPI*(DSQ-P1SQ-6.*PSQ)      !image obscured, guess
      U=PI*L/4.*(DSQ-P1SQ-6.*PSQ)
      W4=AMIN1(W,W0,W1)
10    MASSF=0.
      TWOX=X+X
      XSQ=X*X
      P1P2X=P1+TWOX
      PP2X=P+TWOX
      DM2X=D-TWOX
      LM2X=L-TWOX
      P12XSQ=P1P2X*P1P2X
      PP2XSQ=PP2X*PP2X
      IF(NP.EQ.0)GO TO 2000
      IF(NP.EQ.1)GO TO 3000
      IF(IM2X.GT.0)GO TO 340
c     page 24
      SURF=0
      V=0
      GO TO 850
340   S0=PI*LM2X*LM2X*(D+P1+6.*P+12.*X)+HAFPI*(DM2X*DM2X&
     & -P1P2X*P1P2X-6.*PP2X*PP2X)
      V0=PIFOR*LM2X*(DM2X*DM2X-P1P2X*P1P2X-6.*PP2X*PP2X)
      IF(X.GT.W4/2.)GO TO 360
      MASSF=-TWOX/L/(DSQ-P1SQ-6.*PSQ)
      MASSF=MASSF*(24.*XSQ+(24.*P+4.*P1+4.*D-12.*L)*+P1SQ&
     & +6.*PSQ-2.*L*D-2.*P1*L-12.*L*P-DSQ)
      SURF=S0
      RETURN
360   IF(X.GT.W1/2.)GO TO 390
      F2=0.
      L2=0.
      A3=0.
      A4=0.
      GO TO 460
390   Z=(2.*D1+P+P1+4.*X)/4
      B3=((P1-P)*(P1+P+4.*X)+4.*D1SQ)/4./D1/P1P2X
      A3=ATAN(SQRT(1.-B3*B3)/B3)
      B4=((P-P1)*(P+P1+4.*X)+4.*D1SQ/4./D1/PP2X
      A4=ATAN(SQRT(1.-B4*B4)/B4)
      F2=AR/4.*P12XSQ+A4/4.*PP2XSQ&
     & -SQRT(Z*(Z-D1)*(2.*Z-P-TWOX)*(2.*Z-P1-TWOX))
      L2=LM2X*(A4*PP2X+A3*P1P2X)
460   IF(X.GT.W/2.)GO TO 490
      F3=0.
      L3=0.
      A5=0.
      GO TO 530
490   B5=D1/PP2X
      A5=ATAN(SQRT(1.-B5*B5)/B5)
      F3=(A5*PP2XSQ-D1*SQRT(PP2XSQ-D1SQ))/2
      L3=2.*A5*LM2X*PP2X
530   IF(X.GT.W0/2.)GO TO 560
      F1=0.
      L1=0.
      A1=0.
      A2=0.
      GO TO 650
560   Y=(2.*D1+P+D)/4.
      B1=((D+P)*(D-P-4.*X)-4.*D1SQ)/4./D1/PP2X
      A1=ATAN(SQRT(1.-B1*B1)/B1)
      IF(A1.GT.0.)GO TO 610
      A1=PI+A1
610   B2=((D+P)*(D-P-4.*X)+4.*D1SQ)/4./D1/DM2X
      A2=ATAN(SQRT(1.-B2*B2)/B2)
      F1=A1/4.*PP2XSQ-A2/4.*DM2XSQ+SQRT(Y*(Y-D1)&
     & *(2.*Y-P-TWOX)*(2.*Y-D+TWOX))
      L1=LM2X*(A1*PP2X+A2*DM2X)
650   IF(X.GT.W/2.)GO TO 690
      SURF-S0+12.*(F1+F2+F3)-6.*(L1+L2+L3)
      V=V0+6.*(F1+F2+F3)*LM2X
      GO TO 850
c     Page 25
690   IF(X.LT.X1)GO TO 730
      S1=0.0
      V1=0.0
      GO TO 760
730   S1=3.*D2SQ3-PI*PP2XSQ-HAFPI*P12XSQ&
     & +6.*F3+12.*F2
      S1=S1+LM2X*(2.*(PI-3.*A5-3.*A4)*PP2X+(PI-6.*A3)&
     & *P1P2X)
      V1=LM2X/2.*(3.*D2SQ3-PI*PP2XSQ)&
     & -HAFPI*P12XSQ+.6*F3+12.*F2)
760   IF(X.LT.X2) GO TO 800
      S2=0.0
      V2=0.0
      GO TO 830
800   S2=HAFPI*DM2XSQ-3.*D2SQ3-TWOPI*PP2XSQ&
     & +12.*F1+6.*F3\
      S2=S2+LM2X*((PI-6.*A2)*DM2X+2.*(TWOPI-3.*A1-3.*A5)&
     & *PP2X)
      V2=LM2X/2.*(HAFPI*DM2XSQ-3.*D2SQ3-TWOPI&
     & *PP2XSQ+12.*F1+6.*F3)
830   SURF=S1+S2
      V=V1+V2
850   MASSF=1.-V/U
      RETURN
C
C     ZERO PERF CALCULATIONS START HERE
C
2000  if(d-2*x.le.0.0) go to 2001
      twox=x+x
      xsq=x*x
      MASSF=TWOX*(DSQ+2.*L*D-4.*X*D-TWOX*L+4.*XSQ)/(DSQ*L)
      u=dsq*l*pi/4.
      SURF=PI*(DSQ/2.-4.*D*X-TWOX*L+D*L+6.*XSQ)
      RETURN
2001  surf=0.0
      massf=1.0
      u=dsq*l*pi/4
      return
C
C     ONE PERF CALCULATIONS START HERE
C
3000  (if(d-p-4.*x.le.0.) go to 3001
      twox=x+x
      MASSF=TWOX*(DSQ+2.*L*D-4.*X*D-PSQ+2.*P*L-4.*P*X)&
     & /(DSQ*L-PSQ*L)
      u=dsq*l*pi/4.-psq*l*pi/4.
3001  surf=0.0
      massf=1.0
      u=dsq*l*pi/4.-psq*l*pi/3
      return
      END
