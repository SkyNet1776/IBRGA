      program ibrga
      common nsl I kpr,fracsl(1 0),dsdxsl(10),surfsl( 10),
     & nslp(I0),tsl(10),pbrch,pbase,pmean,bbr( 10),abr(10),
     & deltat,y(20),igrad
      character outfll* 10
      character bdflle* 10
      dimension br(10),trav(10),rp(10),tr(10),forcp(10),tempp(10),covp(
     &10)
      dimension chwp( 10),rhop(10),gamap(10),nperfs(10),glenp(1O),pdp(1O
     &),gdiap(10),alpha(10, 10),beta(10, 10),pres(10, 10)
      dimension a(4),b(4),ak(4),d(20),p(20),z(20),frac(10),surf(10
     &),gdiap(10),alpha(10,10),beta(10,10),pres(10,10)
     &),volp(10),dsdx(10),nbr(10),ibo(10),tbo(10),d2xdt2(10),tng(10)
      real lambdaj 1zp,32zp~j3zp j4zp~jlzbj2zbj3zbj4zb
      dimension chdist(5),chdiam(5),bint(4)
      dimension nsl(10),surfo(10),dsdxn(10)
      L call gettim(ihr~imin,isec,ihuns)
      pi=3. 14159
      v-wite(* 15)
15 format(' input name of data file to te used as input '
      r-ead(*, 10)bdfile
10 format(a10)
      open(unit--2,err--999,file=bdflle,status='old' ,iostat=ios)
      rewind 2
      write(*,25)
25 format(' input name of output file '
      read(*, 10)outfil
      open(unit=-6,err--998,fl~e=outfi1)
      write(6, 16)bdflle
16 format(' the input file is ',alO)
do 9 i=1,20
p(i)=0O.
y(i)=0.
Z(i)=0O.
d(i)=0.
9 continue
readI(2,*,end=20,err-=30)cham,grvc,aland,glr,twst,travp,igrad
&,fs0
if(igrad.gt.l)go to 51
write(6,55)
55 format(lx,'using Lagrange pressure gradient')
igrad=l
go to 52
c define chambrage assumes nchpts;=number of points to define
c chamber > or = 2 < or = 5 (?),chdiam(I) defines chamber diameter
c at chdist (1) chamber distance. chdiam(nchpts) is assumed to be
c the bore diameter and chdist(i) is assumed to be 0, i.e. at the
c breech. Assumes trnmcated cones.
51 if(-igrad.eq.3)go to 401
if(igrad.eq.4)go to 434
write(6,47,err--30)
47 fonnat(lx,'Using chambrage pressure gradient')
go to 436
95
434 write(6,437)
437 format(lx,'using rga gradient')
go to 436
436 read(2,*,end=20,err-=30)nchpts,(chdist(l),chdiam(I),I=l nchpts)
write(6,53,err=-30)(chdist(I),chdiaxn(I),I= 1,nchpts)
53 format(///,' chamber distance cm c~hamnber diameter cm'J(5x~el4
&.6,5x,el4.6))
do 54 1=1 ,nchpts
chdist(I)=O.Ol *chdistqI)
54 chdiam(I)=O.OlI *chdiam(I)
c calculate chamber integrals and volume
if(nchpts.gt.5) write(6,44,err-=30)
44 formnat(lx,'use first 5 points')
if(nechpts.gt.5)nchpts=-5
bore=chdiam(nchpts)
if(chdist(l).ne.O.O)write(6,45,err--3O)
45 format(lx,' # points? '
chdist(l)=-O.O
pi3=piI3.O
bl=-O.O
b2=0.0
b3=0.0
b4=0.0
points=25.O
56 points--points+points
step=chdist(nchpts)/points
ZZ=O.O
bint(l)=-O.O
bint(3)=-O.O
bint(4)=-O.O
bvol=O.O
r2=O.5*chdiam(l)
k=1
j=int(points+O.5)
do 57 I=lj
zz=zz+step
if(k.eq.nzhpts-l)go to 46
do 58 Il=k~nchpts-1
if(zz.gt.chdist(I l).and. zz.ltchdist(I 1+l))go to 59
58 continue
Il=nchpts-I
59 k=II
46 diam=(zz-chdist(k))/(chdist(k+l)-chdist(k))
diam=chdiam(k)+diam*(chdiam(k+l)-chdiam(k))
rl=O.5*diam
areapi*(rl+r2)*(rl+jd2),'4.
bvol=bvol+steplbpi3*(rl *rl+rl *r2+r2*r2)
bint(l)=bint(l)+step~'bvol/area
bint(3)=bint(3)+step*area*bint(l)
bint(4)=bint(4)+step*bvol*bvolIatea
57 r2=rl
temp=abs(l .0-b l/bint( 1))
if(abs( 1.O-b3/bint(3)).gt.temp)temp=abs(l .O-b3/bint(3))
96
if(abs(1 .0-b4/bint(4)).gt.temp)temp=abs(1 .0-b4/bint(4))
if(temp.le.0.001)go to 41
bl=bint(1)
b3=bint(3)
b4=bint(4)
go to 56
41 chiam=bvol*l .e6
c write(6,47,err-=30)bint(1),bint(3),bint(4)
c fonnat(lx,'bint I = ',e14.6,' bint 3 = ',e14.6,' bint 4 ='4e4.
c &M)
chmlen--chdist(nichpts)
go to 52
401 write(6,402)
402 format(lx,'using 2 phase gradient equation')
goto 52
52 write(6,4O,err--30)cham,grve~aland,glr~twst~travp,igrad,fsO
40 fonnat~lx,'chamber volume cm**3',e14.6X/ groove diam cm',el4.6,/
W' land diam cmn',el4.6,/' groove/land ratio',e14.6,/' twist turns
&/caliber ',el4.6,I' projectile travel cm',e14.6
&/' gradient # ',i3X/ friction factor ',e14.611/)
cham=cham*l1.e-6
grve=grve*1.e-2
aland=aland*lI.e-2
travp=travp* 1 e-2
read(2,*,end=20,err-=30)pn.it~iair,htfr,pgas
write(6,50,err=-30)prv-t,iair,htfr,pgas
50 format(lx,'projectile mass kg',e14.6,/' switch to calculate energ
&y lost to air resistance J',i2,/' fraction of work against bore u
&sed to heat the tube',el4.611x,' gas pressure Pa' ,el4.6)
read(2,*,end=20,enr-30)npts,(br(-*,trav(i),i=l npts)
write(6,60,err-=30)npts,(br(i),trav(i),i=1,npts)
60 format(lx,'number barrel resistance points',i21' bore resistance
& M~a - travel cm'/(lx~el4.6,e14.6))
write(6,65)
do 62 i= l,npts
br(i)=br(i)* 1.e6
trav(i)=-trv(i)* 1.e-2
62 continue
65 format(J x)
read(2,*,end=20,err--3o)rcwt,nrp,(rp(i),tr(i),i= 1 JnP)
write(6,70,err-=30)rcwt,nip,(rp(i),tr(i),i=I ,nip)
70 fonnat(lx,' mass of recoiling parts kg',e14.6j' number of recoi
&I point pairs',i2XI recoil force N',' recoil time sec'/,(lx,e14
&.6,3x,014.6))
write(6,65)
read(2,*,end=20,erI--30)ho,tshl,cshl,twal,hl,rhocs
write(6,75,err-=30)ho,tshl,cshl,twal,hIlitocs
75 format(lx,' free convective heat transfer coefficient wlcm**2 K',
&el4.6,f chamber wall thickness cm',e14.6X/ heat capacity of st
&eel of chamber wall Jig K',e14.6,/' initial temperature of chambe
&r wall K',el4.6,/' heat loss coefficient',e14.6,/' density of ch
&amber wail steel g/cm**3'4e4.6/f)
ho=ho/1 .e-4
tshl=tshl* 1 .e-2
97
cshl=cshl*1 .e+3
rhocs--rhocs* 1 .e-31l .e-6
read(2,*,end=20,err-=30)forcig,covi,tempi,chwi,gamai
write(6,85,err--3O)forcig,covi,tempi,chwi,gamnai
85 format(lx,' impetus of igniter propellant J/g',e14.6,/' covolume
& of igniter cm**31g',e14.6,I' adiabatic flame temperature of igni
&ter propellant K',e14.6j' initial mass of igniter kg',el4.6,/' r
&atio of specific heats for igniter',e14.6//)
forcig--forcig* I e+3
coy j-covj"' Ie-611 .e-3
read(2,*)nprop
write(6,98)nprop
98 fornat(' there are ',i2,' propellants')
read(2,*,end=2O,err--3O)(forcp(i),tempp(i),covp(i),chwp(i),
&rhop(i),gamap(i),nperfs(i),glenp(i),pdp(i),gdiap(i),i=1 ,nprop)
write(6,95,err-=30)(i,forcp(i),tempp(i),covp(i),chwp(i)
&,rhop(i),gamnap(i),nperfs(i),glenp(i),pdp(i),gdiap(i),I=1,nprop)
95 format((' for propellant number',i2,I' impetus of propellant Jig
&',e14.6,/' adiabatic temperature of propellant K',el4.6,/' covol
&ume of propellant cm**31g',e14.6/' initial mass of propellant kg'
&,e14.6/' density of propellant glcm**3',e14.6/' ratio of specifi
&c heats for propellant',e14.6/' number of perforations of propell
&ant',i21' length of propellant grain cm',e14.6/' diameter of per
&foration in propellant grains cm',e14.6/' outside diameter of pro
&pellant grain cm',e14.61)If)
tmpi=-O.0
do 96 i=1I,nprop
forcp(i)=forcpQi)* 1.e+3
covp(i)--covp(i)*1 .e-611 .e-3
rhop(i)=rhiop(i)* I.e-311 .e-6
glenp(i)-=glenp(i)*0.0l
pdp(i)--pdp(i)*0.01
gdiap(i)--gdiap(i)*0.O1
tmnpi=tmnpi4-chwp(i)
kpr--i
call prf7lO(pdp(i),gdiap(i),glenp(i),nperfs(i),0.,
& frac(i),volp(i),surf(i),dsdx(i))
tng(i)=chwp(i)/rhop(i)/volp(i)
surfo(i)=-surf(i)
write(6,408)i,tng(i)
408 format(' for propellant ',i2,' the total number of grains'
& ,' is ',el4.6)
96 continue
tmnpi=tmnpi+chwi
do 97 j=1,nprop
read(2,*,end=2O,err-=3O)nbr(j),(alphaoj,i),betaoj,i),presoi,i),
&i=1 ,nbroj))
write(6,1 I 0,err=-30)nbr(j),(alphaoj,i),betaoj,i),presoj,i),
&i-l ,nbrOj))
110 fonnat(lx,'number of burning rate points',i2/3x,' exponent',8x,'
& coefficient', lOx,' pressure'I5x,'-',l5x,'cm/sec-MPa**ai',10x,'MOP
&a'/( lx,el4.6,5x,e14.6,15x,e14.6))
do 112 i=1,nbr(j)
beta~ji)=beta~ji)* 1.e-2
98
presOj,i)=presOj,i)* 1 .e6
112 continue
97 continue
write(6,65)
read(2,*,end=20,err-=30)deltat,deltap,tstop
wiite(6, 120,err-=30)deltat,deltap,tstop
120 format(lx,'time increment msec',e14.6,' print increment msec',e14
&.6/lx,'time to stop calculation msec '4e4.6)
write(*, 130)
deltat=deltat*0.001
deltap=deltap*0.001
tstop=tstop*.O01
130 format(lx,'end input data -- I.B. calculation start')
if(igrad.eq.2.or.igrad.eq.4)go to 131
bore=(glr*grve*grve+aland*aland)/(glr+ 1.)
bore=sqrt(bore)
131 areab=pi*bore*bore/4.
lambda-- 1./((l13.2+4. *logl1O(l 00. *bore))**2)1
iplot=O
pltdt=deltat
pltt=0.
pmaxm=-O.0
pmaxbr=-O.0
pmaxba=0O.0
tpmaxm=-O.0
tpmxbr=0O.0
tpmxba=0O.0
tpmax=O.0
a(1)=0.5
a(2)=1.-sqrt(2.)/2.
a(3)= 1 .+sqrt(2.)/2.
a(4)-1 .16.
b(1)--2.
b(2)=1.
b(3)--1.
b(4)--2.
ak(1)=0.5
ak(2)=a(2)
ak(3)=a(3)
ak(4)=0.5
VP0=0.0
t,0=0.0
tcw=0O.0
if(igrad.eq.3)chmlen=cham/areab
zb=chmlen
zp=chmlen
grlen=0O.
grdiam=0.
egama=O.
do 5 i=1I,nprop
grlen--grlen+chwp(i)*glenp(i)
grdiam=grdiamn+chwp(i)*gdiap(i)
ibo(i)=0O
egama=egama+chwp(i)*gamap(i)
99
nsl(i)=O
5 vp0=chwp(i)/rhopQl)+vp0
volgi=cham-vpO-chwi*covi
grlen7-grien/(timpi-chwi)
grdiam=grdiam/(tmpi-chwi)
egama=<egama+chwi*gamai)/tmnpi
ism_=O
odlnr=-0.
vfD=cham-vp0
epsO=1.-vpO/cham
eps=epsO
gasden--chwi/vlO
prden--tmilp0o
ug=0O.
up=o.
pmean--forcig*chwi/volgi
pbase=pmean
pbrch=pmean
opbase=pmean
volg--volgi
volgi=volgi+vp0
walit--twal
tgas=tempi
told=-0.
tgaso=tgas
dtgaso=0.
COY 1=covi
t=-O.
ptime=0.0
ibrp=I0
.z(3)= 1.
nde=ibtp+nprop
write(6, 132)areab,pmean~vpO,volgi
132 format(1x,'area bore mA2 ',el6.6,/' pressure from ign pa',e16.6,/
&,' volume of unbumnt prop m'A3 ',e16.6,/
&' init chain vol-coy ign mA3 ',e16.6)
wiite(6,6)
6 fonnat(lx,' time acc vel dis inpress
& pbase pbrch '
iswl=0
19 continue
do I1I J=1,4
c FIND BARREL RESISTANCE
do 201 k=2,npts
if(y(2)+y(7).ge.trav(k))go to 201
go to 203
201 continue
k=npts
203 resp=(trav(k)-y(2)-y(7))/(trav(k)-trav(k-1))
resp=br(k)-resp*(br(k)-br(k-1))
c FIND MASS FRACTION BURNED
do 211 k=1,nprop
kpr--k
if(ibo(k).eq. 1)goto2l 1
100
nsll=O
call prf7 1O(pdp(k),gdiap(k),glenp(k),nperfs(k),y(ibrp+k)
&,frac(k),volp(k),surf(k),dsdx(k))
nsl(k)=nsll
if(nsl(k).eq.O)goto 212
if(nslp(k).eq.l)go to 212
write(6,2213)k
2213 format(' propeliantJi2,' has slivered')
nslp(k)--1
tsl(k)=y(3)
ism=1
212 continue
if(frac(k).1L.9999) go to 211
frac(k)= 1.
tbo(k)--y(3)
ibo(k)=1
ism= 1
write(6,456)k
456 format(' propellant',i2,' has burned out')
211 continue
c ENERGY LOSS TO PROJECTILE TRANSLATION
elpt--prwt*y(1 )*y(1)12.
eptdot=prwt*y(l)*z(1)
c ENERGY LOSS DUE TO PROJECTILE ROTATION
elpr--pi*pi*prwt*y(1)*y(1)/4.*twst*twst
eprdot=pi*pi~'prwt*y(l)*z(1)/2.*twst*twst
c ENERGY LOSS DUE TO GAS AND PROPELLANT MOTION
if(igrad.eq.1)go to 214
if(igrad-eq.3)goto 217
iffigrad-eq.4)go to 438
pt=y(2)+y(7)
vzp=bvol+areab*pt
j4zp=bint(4).i((bvo1+areab*pt)**3-bvol**3)/3./area~b/areab
elgpm=tmpi*y(1)*y(1 )*areab*areab*j4zp/2./vzp/vzp/vzp
go to 216
438 pb=y(7)+y(1O)
vzb=bvol+areab*pb
j4zb=bint(4)+(vzb**3-bvol**3)/3./areatbIareab
elgpm=(1 .-eps)*up*up*areab**2*prden*j4zb+
& eps*ug*ug*areab**2*gasden*j4zb
elgpm=elgpml2./vzb/Vzb+gasden*areab*ullen/6.*
& (3-*y(1)*y(1)+3.*y(1)*u11en*diho+ujJlen**2*dlnrjj**2)
c approximate epdot
epdot=trnpi*y(l)*z(1)t3.
go to 216
214 elgpm--tmpi*y(1)*y(1)16.
go to 216
217 elgprn-areab*zb/6.*(eps*gasden*ug*ug+(I .-eps)*prden*up*up)
eIgpm=elgpm4~gasden*areab*ullen/6.*(3.*y( l)*y(1)+ & 3.*Y(1)*ujllen*dlnrhozIlen**2*dphuo**2)
c approximate epdot
epdot=tmpi*y(1)*z( 1)13.
c ENERGY LOSS FROM BORE RESISTANCE
216 elbr--y(4)
101
z(4)=areab*resp*y(l)
ebrdot--z(4)
c ENERGY LOSS DUE TO RECOIL
elrc=rcwt*y(6)*y(6)f2.
erdot--rcwt*y(6)*z(6)
c ENERGY LOSS DUE TO HEAT LOSS
areaw=cham/areab*pi*bore.,2.*areab+pi*bore*(y(2)+y(7))
avden-=O.O
avc=O.O
avcp=O.O
z18=0
Z19=0
do 213 k=l,nprop
z18=forcp(k)*gamap(k)*chwp(k)*frac(k)/(gamap(k)-1 .)Itempp(k)+z18
zl9=-chwp(k)*frac~k)+z19
avdewn-avden+chwp(k)*frac(k)
213 continue
avcp=(zl8+forcig*gamnai*chwi/(gamai-1 .)/tempi)/(z19+chwi)
avdenr-(avden+chwi)/(volg+covl)
avvel=.5*y(l)
huis=lambda*avcp*avden*avvel+ho
z(5)=areaw*htns*(tgas-wallt)*hl
elht--y(5)
ehdot~z(5)
wallt=(elht+htfr*elbr)/cshl/rhocs/areaw/tshl+twaI
c write(6,*)amnbda~avcp,avden,avve1,ho,areaw,hitns~tgas,wal1t~h1,z(5)
c &,elht
c ENERGY LOSS DUE TO AIR RESISTANCE
air--iair
z(8)=y(l)*pgas*air
elar--areab*y(8)
eddot=z(8)*areab
c RECOIL
z(6)=O.O
if(pbrch.le.rp(l)/areab)go to 221
rfor--rp(2)
if(y(3)-ti.ge.tr(2))go to 222
rfor--(tr(2)-(y(3)-tio))I(tr(2)-tr(l))
rfbr-*p2)-rfor*(rp(2)-rp(1))
222 z(6)=areab/rcwt*(pbrch-rfor/areab-resp)
if(y(6).lt.O.O)y(6)=O.O
z(7)=y(6)
goto 223
221 tI%)=y(3)
223 continue
c CALCULATE GAS TEMPERATURE
eprop=0.O
rprop=O.O
dmfogt=-O.O
dmfog=O0.O
do 231 k=1,nprop
eprop=eprop+forcp(k)*chwp(k)*frac(k)/(gamap(k)- I.)
wpo~po+forcp(k)*chwp(k)*frac(k)I(gamnap(k)- I .)Itempp(k) dmfogt=dmfogt+forcp(k)*rhop(k)*tng(k)*surf(k)*z(ibrp+k)/
102
& ((gamap(k)-1 )*tempp(k))
dmfog=dmfog+forcp(k)*rhop(k)*tng(k)*surf(k)*z(ibrp+k)/
& (gamap(k)-l.)
231 continue
ternerg=elpt+elpr+elgpm+elbr+elrc+elht+elar
tgas=(eprop+forcig*chwi/(ganiai-I1.)-elpt-clpr-elgpm-elbr-elrc-elht
&-elar)/(rprop+forcig*chwiI(gamai- 1.)Itempi)
tedot=epdot+eprdot+eddot+ebrdot+erdot+ehdot+eptdot
dtgas=(dmfog-tedot-tgas*dmfogt)/(rprop+forcig*chwi/
& (gamai-I.)/tempi)
c FIND FREE VOLUME
vl=0O.0
covl=O.0
do 241 k=l,nprop
v 1=chwp(k)*(1.-frac(k))/rhop(k Wv 1
coy 1=cov 1+chwp(k)*covp(k)*frac(k)
241 continue
volg--volgi+areab*(y(2)+y(7))-vl1-coy I
c CALCULATE MEAN PRESSURE
rl=0O.0
do 251 k=I,nprop
rl=rl+forcp(k)*chwp(k)*frac(k)/tempp(k)
251 continue
pmean--tgaslvolg*(rI+forcig*chwi/tempi)
259 represp+pgas*air
if(igrad.eq.1)go to 252
if(igrad.eq.2)goto 403
if(igrad.eq.3)go to 404
if(i.grad.eq.4)go to 441
403 if(iswl.ne.0)go to 253
pbase=pmean
pbrch=pmean
if(pbase.gt.resp+1 .)iswl=1
go to 257
c USE CHAMBRAGE PRESSURE GRADIENT EQUATION
253 j lzp=bint(l)+(bvol*pt+areabf2.*pt*pt`)/areab
j2zp=(bvol+areab*pt)**2/areab/areab
j3zp=bint(3)+areab*bint(I)*pt+bvol*pt*pt/2.+areab*pt*pt*pt/6.
a2t---tmpi*areab*areao/prwt/vzp/vzp
alf=1.-a2t*jlzp
alt=tmnpi*areab*(alvab*y(1)*y(1)/vzp+areab*resp/prwt)/vzp/vzp
bt---tmpi*y(l)*y(1)*areab*areab/2./Vzp/vzp/vzp
bata=--alt*j lzp-bt*j2zp
gamma=alf+a2t*j3zp/vzp
delta=bata+alt*j3zplvzp+bt*j4zp/vzp
c calculate base pressure
pbase=(pmean-delta)/gwmma.
c calculate breech pressure
pbrch=alf~pbase+bata
go to 254
c USE 2 PHASE GRADIENT EQUATION
404 IF(ISWI.NE.0)GOTO 407
pbase=pmean
pbrch=pmean
103
if(pbase.gtamsp+l)iswl11
go to 257
407 if(iswl.eq.2)go to 411
vzp=cham+aw~ab*(Y(2W+y(7))
vzb=cham+areab*(y(1O)+y(7))
phi=O.
phidot=0.
dmorho=0.
dmcov=0O.
dmromw=O.
rrnaomw=O.
vfime=vzp-vl
do 405 k=1,nprop
nnomw=nnomw+chwp(k)*frac(k)*forcp(k)ItemPP(k;'
phi=chwp(k)*frac(k)+phi
if(ibo(k).eq.1)go to 405
dmorho--dmorho+tng(k)*surf(k)*z(ibrp+k)
phfidot7-rhop(k)*tng(k)*surf(k)*z(ibrp+k)+phidot
dmcov=rhop(k)*tng(k)*surf(k)*z(ibrp+k)*covp(k)+dmcov
dmromw=dmromw+lbop(k)*tng(k)*surf(k)*z(ibrp+k)*
& foivp(k)Itempp(k)
405 continue
rmomw=rinomw+chwi*forcigttempi
gasmas=phi+chwi
gasdenr-gasmnaslvfree
phi=(phi+chwi)/tmPi
if (phi.gt.O.999) then
iswl=2
rbm=pbaselpmean
rbrm=pbrcb/pmean
if(phi.ge.l.)go to 411
endif
dmdt=phidot
phidot--phidotltmpi
vdotov=(dmorho+areab*y(l))Ivfree
dlnrho=dmtdtlgasmas-vdotov
dvoldt=dmorho+areab*y(l)-dmcov
c GET TIME DERIVATIVE OF MEAN PRESSURE
dpmndt=(dmromw*tgas-pmean*dvoldt+dtgas*rmomw)Ivolg
volpqp=O.
effdia=0O.
dmdmdt=0.
dmdmor=0O.
avelen=0.
avedia=0.
do 406 k=lnprop
if(ibo(k).eq.l)go to 406
volprp~volprp+(l.-frac(k))*chwp(k)Irhp(k)
dmdmdt=dmdmdt+rhop(k)*tng(k)*dsdx(k)*z(ibrp+k)*z(ibrp+k)
dmdmdt=dmdmdt+rhop(k)*tng(k)*surf(k)*d2xdt2(k)
dmdmor--dmdmor+(dsdx(k)*z(ibrp+k)**2+surf(k)*d2xdt2(k))*tng(k)
effdia=effdia+6.*volp(k)/surf(k)*(1 .-frc(k))*chwp(k)
406 continue
clt=dmdmdt/gasmas-dmdmortvfree+vdotov**2-(dmdt/gasmas)**2
104
d2Inr--c It-areab* *2 *pbase/v free/prwt
d2 Inr--d21nr+ areab* areab* re sp/v free/prwt
zp=chnilen+y(2)+y(7)
zb=chmlen+y(1O)+y(7)
ullen=zp-zb
cnow=tmpi-gasmas
vp=y(1)
effdia=effdia/cnow
prden--cnow/volprp
up=Y(9)
phistr--phi-gasden*areab*ullen/tmpi
ulldot--vp-up
dphist--phidot-gasden*areab/tlnpi*(uJ~dot+ullen*dlnrho)
eps=l.-(1.-phi)*tmnpi/prden/vzb
epsdot=phidot*tmpi/prden/vzb+(1 .-phi)*tmpi*up*areab/
& prden/tvzb/vzb
ug--up+(vp+ullen~dlnrho-up)/eps
alam=(1 .5*grlje/grdialn)**.666666667
alam=(O.5+grleri/grdiam)/Iaam
alt~am**2.17
c VIS kgs/sin
vis=.00007
ren--gasden/vis*effdia*abs(ug-up)
if(ren.1t. 1.)ren=1.
fsrg--2.5*alaln/ren**.08 I *((1 .-eps)/I1.-epsO)*eps0leps)**.45
fsc=fsrg*fsO
phi2-l .-phi-phistr"'(1.-eps)Ieps
phi lp=dphist*ug-phidot*up-phistr*epsdot/eps/eps
"& *(vpullen*dlnrho..up)+phit*ujldot*dlndio/eps
"& +2.*phistr*ug/zb*(ug-up)
philp=phi 1p+phi2*gasden/effdia/prden*(ug-up)**2*fsc
ak2=l1.(1 .-phi2*tmpi/prden/vzb)
phi 1=phi lp+phistr*z(1)/eps+u11en*phistr*d2Inr/eps
c ACCELERATION OF FORWARD BOUNDARY OF PROPELLANT BED
Z(9)1=a~r2qdon*(Ua-Up)**2*fs.ctprdenx/effdia+tmpi*phi1*ak2
&tvzb/prden
z(1O)=Y(9)
e=phistr/eps*( 1.-ullen*areab/vfree)*area~b/Prwt
dd=ullen*phistr*clt/epr,
akI 1=tmnpi*e*ak2lzb/vzb
akl2=tmpi*ak2*(phi lp+dd)/zb/vzb-akl I *rep
pbase=pmean-akl2*zb*zb/2.+gasden*ullen*resp*aea~b/Prwt
pbase=pbase+akl2*zb*zb*(zbI3.+ullen)t2./zp
pbase=pbase-gasden*ullen**2*areab*resp/2./zp/prwt
pbase=pbase-gasden*ullen**2I2.*(1 .-2.*u1Ien/3./Zp)*
& (ci t-dlnrho**2)
pbase=pbase-arcab**2*gasden*ullen**2*
&(1 .-2.*ullen./zp)*1esp/prwt/vfrect2.
deno=-akl 1 *zb**3/6 /Zp..ullen*akl 1 *zb*zb/2I/Zp
deno=deno+gasden*ullen*areab/prwt-areab**2*gasden*ullen**2
&*(1 .-2.*ullen/3./zp)t2./vfree/prwt
deno=deno-gasden*u11en**2*areab/2./Zp/prwt+1 .+akl I *zb*zb/2.
pbase=pbase/deno
if(ism.eq.O)goto453
105
if(ism.cq. 1)goto451
goto452
451 ism=2
ts~qtean~mm/aia~gs
write(6,*)tss
tss~ullen/(ullen*odlnr+tss)
tso=y(3)
write(6,*)tss,tso
452 coefbp=(tss+tso-y(3)-deltat)/t~s
if(coefbp.gt. I .)coefbp= 1.
if(coefbp.le.O.)then
coefbp=O.
ism=-O
endif
pbase=coefbp*opbase+(1.-coefbp)*pbase
write(6,*)coefbp,opbase,pbase,ism
453 odlnr--dlnrho
opbase=pbase
pbrch=pbase*(1 .+akl I1*zb*zbI2.+gasden*ullen*areablprwt
& -areab**2*ga~den*ullen**2/2.IvfreeIprwt)
pbrch=pbrch+ak,.2*zb*zb/2.-gasden*ullen*areab*respp4rwt
pbrch~pbrch+gasden*ullen**2/2.*(clt~dflnt1o**2)
pbrch=pbrch+areab**2*gasdefl*u~lef**2*resp/2./vfree/prwt
go to 254
C USING RGA GRADIENT
441 if(iswl.ne.O)go to 444
pbase=pmnean
pbrch=pmean
if(pbase.gt~iesp+1 .)iswl=1
go to 257
444 if(iswl.eq.2)go to 411
V7.p=cham+areab*(y(2)+y(7))
vzb=chan1+areab*(y(1O)+y(7))
j Ilzb=bint(1)+(bvol*pb+areabl2.*pb*pb)/areab
j2zb=(bvol+areab*pb)**2/areab/areab
j3zb=bint(3)+areab*bint(1)*pb+bvol*pb*pbt2.+area~b/6.*pb**3
phi=O.
phidot=-O.
dmorho=O.
dmcov=O.
dmromw=O.
rmomw=O.
vfree=vzp-v I
do 442 k=1,nprop
rmomw=rmomw+chwp(k)*frac(k)*forcp(k)Itempp(k)
phi=chwp(k)*frac(k)+phli
ii(it*o~).eq.1)go to 442
dmorho=dmorho+tng(k)*sur(k)*z(ibrp+k)
phidot=rhop(k)*tng(k)*surf(k)*z(ibrp+k)+phidot
dmcov=rhop(k)*tng(k)*surf(k)*z(ibrp+k)*covp(k)+dmcov
dmromw=dmromw+rhop(k)*tng(k)*surf(k)*z(ibip+k)*
& forcp(k)Itempp(k)
442 continue
zrmomw--rmomw+chwi*fbrcig/tempi
106
gasmas=phi+chwi
gasden--gasmaslvfree
phi=(phi+chwi)Itmpi
if (phi.gt.0.99) then
iswl=2
rbmn=pbase/pmean
rbrm=pbrch/pmean
if(phi.ge.l.)go 9 411
endif
dmdt--phidot
phidot=phidot/tmpi
vdotov=(dmorho+areab*y(l))Ivfree
dhlnrhodmdt/gasmas-vdotov
dvoldt=dmorho+artab*y(1).dmcov
c get time derivative of mean pressure
dpmdt=(dmmmw*tgas-pmean*dvoldt+dtgas*rmomw)Ivolg
volprp=O.
effdia=-O.
dmdmdt=0O.
dmdmor=O.
aveleng0.
avedia=O.
do 443 k=1,nprop
if(ibo(k).eq.1)go to 443
volprp=volprp+(1 .-frc(k))*chwp(k)Irhop(k)
dmdmdt--dmdmdt+rhop(k)*tng(k)*dsdx(k)*z(ibtp+k)*z(iblp+k)
dmdmndt~~dmdmdt+rhop(k)*tng(k)*surf(k)*d2xdt2(k)
dmdmor--dmdmor+(dsdx(k)*z(ibrp+k)**2+surf(k)*d2xdt2(k))*tng(k)
effdia=effdia+6.*volp(k)/surf(k)*(1.-frac(k))*chwp(k)
443 continue
clt=dmdmdt/gasmas-dmdmnortvfree+vdotov**2-(dmdt/gasmas)**2
d21nr--clt-amzl*2'*pbase/vfree-/Prwt
d2Inr--d2Inr+areab*areab*resptvfree/prwt
zp=chmlen+y(2)+y(7)
zb=chmlen+y(10)+y(7)
ullen7-zp-zb
cnow--tmpi-gasmas
vp=y(1)
effdia=.4'fdia/cnow
piden--cnowlvolprp
up=Y(9)
phistr--phi-gasden*areab*u11en/tmpi
ulldot=vp-up
dphist=phidot-gasden*areab~txpi*(uUdot+ullen*dnh~o)
eps=l.-(1.-phi)*tmnpi/prden/vzb
epsdot=phiidot*tmpi/prden/vzb+(1 .-phi)*tmpi*up*area~b/
& prden/vzb/vzb
ug--up+(vp+uUen*dlnrho-up)/eps
alam=(1.5*grlen/grdiam)**.666666667
alam=(0.5+grlen/grdiam)Ialamn
alamaJa**2.17
c VIS kg/s/rn vis=.00007
ren=gasden/vis*effdia*abs(ug-up)
107
if(ren.lt.l1.)ren--1.
fsrg=2.5*alam/ren**.08 I *((1 .-eps)/(1 .-epsO)*cpsO/eps)**.45
fsc=fsrg*fsO
phi2=1.-phi-phistr*(1.-eps)/eps
phi 1p=dphist*ug-phidot*up-phistr*epsdot/eps/eps
& *(vp~ulen*dlnrho.up)+phistr*ulldoL*dlprho/eps
& +2.*areab*phistr*ug/vzb*(ug-up)
phi lp=phi 1p+phi2*gasden/effdia/prden*(ug-up)**2*fsc
ak2=1 .1(1 .phi2*tmpi/prden/vzb)
phi 1=phi lp+phistr*z( 1)/eps+ullen*phistr*d2lnr/eps
c ACCELERATION OF FORWARD BOUNDARY OF PROPELLANT BED
z(9)=gasden*(ug-up)**2*fsc/prdent/effdia.ItmpiV-phi 1*ak2
&Ivzblprden
z(1O)=Y(9)
phi3=phistr*ug*ug+(1 .-phi)*up*up
e=l1.-ullen*areabtvfree
dd=ullen*phistr*cltleps
alt---mpi*areab/vzb/vzb*(phi3*areab/vzb-(phi lp+dd-e
& *phistr~areab*resp/eps/prwt)*ak2)
a2t=(-tinpi*e*phistr*areab**2/vzb/vzb/eps/prwt)*ak2
bt---tmpi*phi3*areab**2/2./Vzb/vzb/vzb
pbase=pmean-gasden*ullen**2/2.*(clt-dlnrho**2
$ +areab**2*resp/prwt/vfree)*(1 .-2.*areab*u11en/3.Ivzp)
pbase=pbase-alt*j3zb/vzp-bt*j4zb/vzp-areab*ullen*alt*j lzb/vzp
pbase=pbase-areab*bt*ullen*j2zbtvzp-gasden*areab**2*ullen**2
&*respf2./Vzp/prwt+alt*jlIzb+bt*j2zb+areab*gasden*ullen*resp/prwt
deno=1 .+areab*ullen*a2t*j izb/vzp-gasden*areab**2*ullen**2
& t2./vzplprwt+a2t*j3zb/vzp-a2t*j lzb+gasden*ullen*areab/prwt
deno=deno-gasden*ullen**2*areab**2t2./vfree/Prwt
& +gasden*areab**3*ullen**3/3./Vzp/vfirecýprwt
pbase=pbase/deno
pbrch=pbase*(1 ..a2t*j lzb+gasden*ullen*areab/prwt
& -gasden*ullen**2*areab**2/2./vfree/prwt)
&+gasden*ullen**2t2.*(clIt-dlrn-tio**2+areab**2*resp/prwt/vfree)
& -alt*j 1zb-bt*j2zb-areab*gasden*ullen*resp/prwt
go to 254
411 pbase=rbm*pmean
pbrch=rbrin*pmean
go to 254
c USE LAGRANGE PRESSURE GRADIENT EQUATION
252 if(iswl1.ne.O)go to 256
if(pmean.lLresp)resp=pmean
c CALCULATE BASE PRESSURE
256 pbase=(pmean+tmnpi*resp/3./prwt)/(1 .+tmnpil3./Prwt)
if(pbase.gt.resp+1 .)iswl=1
c CALCULATE BREECH PESSURE
pbrch=-pbase+tinpi/2./prwt*(pbase-resp)
c CALCULATE PROJECTILE ACCELERATION
254 z(1 )=areab*(pbase-resp)/prwt
if(z(1).lt.O.O)go to 257
go to 258
257 if(iswl .eq.O)z(l)=O.O
258 if(y(l).lt.O.O)y(1)=-O.O
z(2)=y(1)
108
c GET BURNING RATE
do 264 m=1,nprop
z(ibrp+m)=O.O
d2xdt2(m)=O.O
if(ibo(m).eq.1) goto 264
do 262 k=l,nbr(m)
if(pmean.gt.pres(m,k))go t~o 262
go to 263
262 continue
k=nbr(m)
263 pmix=pinean
if(igrad.eq.3)pmix=pbrch-(akl 1 *pbase+akl2)/6.*zb*zb
if(igrad.eq.4)pmix--pbrch+(alt+a2t*pbase)*j3zb/vzb+bt*j4zb/vzb
if(Ipmix.ILt.99*pmean)ptnix=pmean
z(ibrp+m)--beta(m,k)*(pmix*lI.e-6)**alpha(m~k)
abr(m)=alpha(m,k)
bbr(m)=beta(m,k)
d2xdt2(m)=beta(m,k)*alpha(m,k)*(pmix*lI.e-6)
& **(alpha(m~k)-1 )*dpmdt*1 .e-6
264 continue
do 21 i=1,nde
y(i)=deltat*d(i)+y(i)
21 continue
11 continue
t=t+deltat
told=y(3)
if(pmaxm.gLpmean)go to 281
pmnaxm=pmean
tpmaxm=y(3)
281 if(pmaxba.gt.pbase)go to 282
pmaxba--pbase
tpmxba7-y(3)
282 if(pmaxbr.gt.pbrch)go to 283
pmnaxbr--pbrch
tpmxbr--y(3)
283 continue
if(y(3)JIt.ptime)go to 272
ptime=ptime+deltap
pjt=y(2)+y(7)
write(6,7)y(3),z(I),y(1),pjt,pmean,pbase,pbrch
7 format(lx,7el 1.4)
if(igrad.gt.2)then
pjt=y(2)+y(7)
prt=y(1O)+y(7)
wiite(6,427)prt,pjt
427 forlnat(lx,'prop travel',ell1.4,'proj travel',el 1.4)
endif
272 continue
if(t.gt.tstop)goto 200
if(y(2)+y(7).gt.travp)go to 200
rmvelo=y(l)
tinvelo=y(3)
109
disto=y(2)+y(7)
go to 19
200 write(6,31 1)t,y(3)
311 format(lx,' deltat t', e14.6, ' intg t',e14.6)
write(6,3 12)pmaxm,tpmaxm
312 format(1x,'PMAXMEAN Pa. ',e14.6,' time at PMAXMEAN sec ',e14.6)
wnite(6,3 13)pmaxba,tpmxba
313 fonnat(lx,'PMAXBASE Pa. ',e14.6,' time at PMAXBASE sec ',e14.6)
write(6,3 14)pmaxbr,tpmxbr
314 format(lx,'PMAXBREECH Pa ',e14.6,' time at PMAXBREECH sec ',el4.6)
if(y(2)+y(7).le.travp)go to 303
dfract=(travp-disto)/(y(2)+y(7)-disto)
rmvel=(y(l)-rmvelo)*dfract+rinvelo
tmnve1=(y(3)-tmvelo)*dfract+tmnvelo
write(6,3 18)rmvel,tmnvel
318 format(lx,'muzzle VELOCITY rn/s ',e14.6,' time of muzzle velocity s
&ec ',e 14.6)
goto 319
303 write(6,327)y(1),y(3)
327 format(lx,'velocity of projectile rn/s ',e14.6,' at this time msec
&',e14.6)
319 efi=chwi*fbrcigI(gamai-l.)
efp=O.0
do 315 i=1,np~rop
efp=efp+chwp(i)*forcp(i)/(gaznap(i)-1 .0)
315 continue
tenerg=efi+efp
write(6,3 17)tenerg
317 format(lx,'total initial energy available J =',el4.6)
tengas=chwi*forcig*tgas/(gamai I .)/tempi
do 135 i=1,nprop
tengas=(frac(i)*chwp(i)*forcp(i)*tgas/tempp(i)/(gamap(i)-l .))+teng
&as
write(6,328)i,frac(i),tbo(i)
328 format(' FOR PROPELLANT ',12,' MASSFRACT BURNT IS ',e14.6
&,' at time in sec ',el4.6)
135 continue
write(6, 136)tengas
136 format(lx,'total energy remaining in gas J= ',e14.6)
write(6,320)elpt
320 format(lx,'energy loss from projectile translation J= ',e14.6)
write(6,32 1)elpr
321 format(lx,'energy loss from projectile rotation J= ',e14.6)
write(6,322)elgpm
322 format(lx,'energy lost to gas and propellant motion J= ',e14.6)
write(6,323)elbr
323 format(lx,'energy lost to bore resistance J= ',e14.6)
write(6,324)elrc
324 format(lx,'energy lost to recoil J= ',el4.6)
write(6,325)elht
325 format(lx,'energy loss from heat transfer J= ',e14.6)
write(6,326)elar
326 format(lx,'energy lost to air resistance J= ',e14.6)
stop
110
20 write(*,140)
140 format(lx,'end of file encounter')
stop
30 write(*,150)
999 continue
998 continue
150 format(1 x,'read or write error')
stop
end
subroutine prf~l1O(pd,gd,gl,np,x,frac~vol,surf~dsdx)
common nsl~kpr~fracsl(10),dsdxsl(10),surfsl(10),
& nslp(10),tsl(10),pbrch,pbase,pmean,bbr(10),abr(IO),
& deltat~yar(20) 'grad
dimension ts(10),coef( 10)
pi=3. 141593
nsl=0
C
C pd=perforation diameter
C gd=OUTER DIA
C gl=GRAIN LENGTH
C NP=-NUMBER OF PERFS
C
C SURF=OUTPUT SURFACE AREA
C frac=OUTPUT MASS FRACTION OF PROPELLANT BURNED
C
C w =web = distance between perforation edges
C =distance between outside perf edge and edge of grain
C
C p =distance between perforation centers
C
C xI = distance to inner sliver burnout
C
C x2 = distance to outer sliver burnout (frac=1)
C
if(np.eq.0) go to 2000
IF(NP.EQ.1)GO TO 3000
IF(NF.eq.7)GO TO 61
if(np.eq.19)go to 4000
if(np.eq.15)go to 5000
60 WRITE(6,90)
90 FORMAT(1X, 'UNACCEPTABLE GRANULATION')
STOP
61 w=(gd-3.*pd)I4.
d=w+pd
sqr3=sqrt(3.)
xl=d/sqr3-pd/2.
x2=(14.-3.*sqr3)*d/13.-pd/2.
v0=piI4.*g1*(gd*gd-7.*pd*pd)
sO=2.*vO/gl+pi*gl*(gd+7.*pd)
if (x.gtwt2.+.0000001) goto 20
vol=pi/4.*(gl.2.*x)*((gd-2.*x)**2-7.*(pd+2.*x)**2)
surf=2.*volI(gl-2.*x)+pi*(g1-2.*x)*((gd-2.*x)+
& 7.*(pd+2.*x))
frac=I.-vol/vO
dsdx=-4*pi*(gd+7.*pd-3.*g1+1 8.*x)
dsdxsl(kpr)=dsdx
fracsl~kpr)=frac
surfsl(kpr)--surf
return
20 nsl=1
coef(kpr)=0O.
if(igrad.eq.l1.or.igrad.eq.2)go to 726
if(nslp(kpr).eq. I)goto 26
tsl(kpr)--yar(3)
ts(kpr)=w/2.*(- 1 .+(pbrcb4,mean)**abr(kpr))/
& (bbr(kpr)*(pbase*1 .e-6)**abr(kpr))
26 continue
coef(kpr)=-(ts~kpr)-itsl(kpr)-(deltat+yar(3)))Its(kpr)
if(coef(kpr).gt. 1 .)cocf(kpr)= 1.
if(coef(kpr).lt.0.)coef(kpr)=O.
726 if(x.ge.x2)goto 30
sl=O.
s2=0.
vl=0.
v2=0.
dsldx=0O.
ds2dx=0O.
y=sqrt((pd+2.*x)**2-d*d)
theta=atan(y/d)
a1=theta/4.*(pd+2.*x)**2-d/4.*y
if(x.ge.xl)goto 25
v 1=3./4.*(g1-2.*x)
vl=v1*(2.*sqr3*d*d-pi*(pd+2.*x)**2+24.*a1)
sl=2.*v 1I(g1-2.*x)
s1=s1+3.*(g1-2.*x)*(pi-6.*thieta)*(pxi+2.*x)
25 yl=sqrt((gd-2.*x)**2-(5.*d-2.*(pd+2.*x))**2)
chi=atan(y1I(5.*d-2.*(Pd+2.*x)))
y2=sqrt((pd+2.*x)**2-(3.*d-2.*(pd+2.*x))**2)
Phi=atan(y2/(3.*d-2.*(Pd+2.*x)))
a2=phi*(pd+2.*x)**2-chij*(gd-.2*x)**2
a2=(a2+2.*sqr3*d*sqrt((3.*d-pd-2.*x)*(3.*d-gd+2.*x)))I8.
v2=pi*(gd-2.*x)**2-6.*sqr3*d*d~4.*pi*(pd+2.*x)**2
v2=(v2+24.*(al+2.*a2))*(g1-2.*x)/4.
s2=2.*v2/(g1-2.*x)
s2s+g-.x*(i6*h)(d-.x+.(.p-.pi3*ht
& )*(pd+2.*x))
vol--vl+v2
surf--sl+s2
frac=l.-vol/vO
dsdx=--surf/(x2-x)
dsdx=coef(kpr)*dsdxsl(kpr)+(1 .-coef(kpr))*dsdx
dsdxsl(kpr)=dsdx
frac=coef(kpr)*fracsl(kpr)+(1 .-coef(kpr))*frac
fracsl(kpr)=frac
surf--cof(kpr)*surfsl(kpr)+(1 .-coef(kpr))*surf
surfsl(kpr)=surf
return
30 vol=0.
112
surf=0O.
frac=fracsl(kpr)*coef(kpr)+1 .-coef(kpr)
hracsl(kpr)=frac
if(frac.gt..9999) frac=l.
if(frac.gt..9999)return
dsdx=g.
dsdx=dsdxsl(kpr)*coef(kpr)
dsdxsl(kpr)=dsdx
if(abs(dsdx).ILt 1 .)dsdx=0.
surf-surifsl(kpr)*coef(kpr)
surfsl(kpr>)surf
return
C
C ZERO PERF CALC4ULATIONS START HERE.
C
2000 if(gd-2.*x.le.0.0) go to 2001
vO=pi*gd*gd/4.*gl
vo1l-pi*(gd-2.*x)**2I4.*(g1.2.*x)
frac=1.-voJ/vO
suf-i2*g-.x*2p*g-.x*g-.x
dsdx=-2.*Pi*(gd+g1.6.*x)
return
2001 surf=0.
frac=1 .0
Vo1=0.
dsdx=0.
nsl=l
return
C
c one perf calculation starts here
C
3000 if(gd-pd-4.*x.le.0.0) goto 3001
v0=pi/4.*(gd*gd-pd*pd)*gl
vo1l~pi/4.*((gd-2.*x)**24(pd+2.*x)**2)*(gl.2.*x)
frac=l.-vol/vO
surzf=PJ2.*((gd2.2*x)**2-pjj+2.*x)**2)
surf--surf+Pi*(gd-2.*x)*(gl-2.*x)
surf--surf+pi*(pd+2.*x)*(gl-2.*x)
dsdx--4.*pi*(gd+pd)
return
3001 surf=0.
frac= 1.0
vol=O.
dsdx=0.
nsl=l
return
C
C Below is the calculation for the cylindrical 19 perf grain.
C
C INPUT
C
C P =PERF DIAM1ETER
C D =GRAIN DIAMETER
C GL = GRAIN LENGTH
113
C X = DISTANCE BURNT
C
C OUT`PUT
C
C VOL = THE VOLUME OF ONE GRAIN AT X.
C SURF = THE SURFACE AREA OF ONE GRAIN AT X.
C FRAC = THE FRACTION OF GRAIN BURNT AT X.
C
C W=WEB
C
4000 p=pd
d=gd
W=(D-5.*P)/6.
P1=3.141592654
SQRT3=SQRT(3.)
SQRT5=SQRT(5.)
SQRT6--SQRT(6.)
SQRT1O=SQRT( 10.)
C
C INITIAL VOLUME AND SURFACE AREA
C
V0=PI/4.*GL*(D*D-I9.*P*P)
S0=2.*VO/GL+PI*GL*(D+ 19.*P)
C
C X1 = DISTANCE TO INNER SLIVERR. BURNOUT
C X2 = DISTANCE TO OUTER SLIVER BURNOUT
C DBC = DISTANCE BETWEEN PERFORATION CENTERS
C ASSUMES BURNOUT DOES NOT OCCUR IN LONGITUDINAL DIRECTION
C WI = SECONDARY WEB
C
DBC=-W+P
Wl1=0.5*(D-P-2.*SQRT3*DBC)
X1=DBC/SQRT3-Pt2.
X2=0.25*(DBC*(6.-SQRT10)-2.*P)
IF(X.GT.W/2.)GO TO 110
C
C NOT SLIVERED YET
C
VOL=-P1/4.*(GL-2.*X)*((D-2.*X)**2- 19.*(P+2.*X)**2)
SURF=2.*VOL/(GL-2.*X)+PI*(GL-2.*X)*(D-2.*X+ 19.*(P+2.*X))
dsdx=pi*(~4*D+36*GL-76*P-216*x)
FRAC=I.-VOL/VO
dsdxsl(kpr)--dsdx
fracsl(kpr)=frac
surfsl(kpr)=-surf
RETURN
C
C VI=TOTAL VOLUME OF INNER SLIVER, V2=TOTAL VOLUME OF OUTER SLIVER
C SI=TOTAL SURFACE AREA OF INNER SLIVERS, S2-=TOTAL SURFACE AREA OF
C OUTER SLIVERS
C
110 V1=0o.
V2=0.
51=0o.
114
S2=0.
DELTA=0O.
CHI=-O.
NSL--l
cocf(kpr)=0O.
if(igrad.eq.l1.or.igrad.eq.2)go to 727
if(nslp(kpr).eq. 1)goto 728
tsl1k pr)=-yar(3)
ts(kpr)=-w/2.*( I .+(pbrcli4,mean)**abr(kpr))/
& (bbr(kpr)*(pbase* 1.e-6)**abr(kpr))
728 continue
coef~kpr)=-(ts(kpr)+tsl(kpr)-(deltat+yar(3)))/ts(kpr)
if(coef(kpr).gt. 1 .)coef(kpr)= 1.
if(coef(kpr).1LO0.)cocf(kpr)=0.
727 A3=0O.
IF(X.GE.X2)GO TO 130
THETA=ACOS(DBC/(P+2.*X))
A1=THETAI4.*(P+2.*X)**2-DBC/4.*SQRT((P+2.*X)**2-DBC*DBC)
IF(X.GT.X1)GO TO 120
V1=3.*(GL-2.*X)*(2.*SQRT3*DBC*DBC.PI*(P+2.*X)**2+24*Al)
SF=2.*V1/(GL-2.*X)+12.*(GL-2.*X)*(PI-6.*THETA)*(P+2.*X)
120 PHI=ACOS((5.*D-1 3.*P.36.*XQ/(12.*(P+2.*X)))
XI=ACOS((1 3*D-.5*P.36.*X)/(12.*(D2.2*X)))
IF(X.LE.W112.)GO TO 125
DELTA=ACOS((2.*D-P-6.*X)/(SQRT3*(D-2.*X)))
CHI=ACOS((D-2.*P-6.*X)/(SQRT3*(P+2.*X)))
A3=. 125 *(CHI*(P+2.*X)**2..DELTA*(D-.2*X)**2
& +2.*SQRT6*DBC*SQRT(6.*DBC*(P+2.*X-DBC)-(P+2.*X)**2))
125 A2=. 125*(PHI*(P+2.*X)**2-XI*(D-2.*X)**2
& +2.*SQRT5*DBC*SQRT((5.*DBC-P-2.*X)*(5.*DBC-D+2.*X)))
V2=.25*(GL-2.*X)*(PI*(D-2.*X)**2-7.*PI*(P+2.*X)**2
& -24.*SQRT3*DBC*DBC+48.*(Al+A2+A3))
S22*2(L2*)(L2*)(D2*)(I6*X+ET)
& +(P+2.*X)*(7.*PI-6.*(2.*THETA+CHI+PHI)))
VOL=-VI+V2
SURF=S 1+52
DSDX=-SURFI(X2-X)
FRAC=1.-VOL/VO
dsdx=coef(kpr)*dsdxsI(kpr)+(1 .-coef(kpr))*dsdx
dsdxsl(kpr)=dsdx
frac=coef(kpr)*fracsl(kpr)+(1 .-coef(kpr))*frac
fracsl(kpr)-frac
surf=cocf(kpr)*surfsl(kpr)+(1 .-coef(kpr))*surf
surfsl(kpr)--surf
RETURN
130 VOL,=0.
SURF=0.
frac=fracsl(kpr)*coef(kpr)+1.-coef(kpr)
fracsl(kpr)=-frac
if(frac.gt..9999) frac=1.
if(frac.gt..9999)retum
dsdx=0.
dsdx--dsdxsl(kpr)*coef(kpr)
dsdxsl(kpr)=dsdx
115
if(abs(dsdx).1L 1.)dsdx=O.
surf=surfsl(kpr)*coef(kpr)
surfsl(kpr)=surf
RETURN
C
C Below is the calculation for the 19 perf hex grain.
C
C
C Translation of the input values.
C p= perf diameter
C d= grain diameter
C gl= grain length
C x= distance burnt
C
C Translation of the output values.
C vol= volume of one grain at x.
C surf= surface area of one grain at x.
C frac= mass fraction of the grain burnt at x.
C
C Assignment statement for pi.
5000 pi=3.141592654
sqrt3=sqrt(3.)
p=pd
d=gd
C
C d=6w + 5p is the statement for the grain diameter which will be
C used to calculate the web.
C
C To calculate the web.
w= (d-5.*p)/6.
C
C Below is the equation to calculate the distance between the perf cenC ters.
dpc= p + w
C To calculate the grain diameter between the flats.
f= 2.*(sqrt3*dpc + p/2. + w)
C
C To calculate the distance burnt
Xl=dpc/sqrt3-p12.
X2= 0.125*(5.*dpc-4.*p)
C
C To calculate the area.
A=sqrt3/3.*((w+p/2.)**2)-pi/6.*((w+p/2.)**2)
C To calculate the initial volume of the sharp comer grain.
Vs=gl/4.*(2.*sqrt3*f**2-19.*pi*p**2)
C
C To calculate the volume that will be removed from the grain.
Vr=6.*A*gl
C
C To calculate the initial volume for the 19hex grain with rounded
C comers.
Vo= Vs - Vr
C
C To calculate the initial surface area of the sharp comer grain.
116
Ss=2.*Vs/gl+gl*(2.*sqrt3*f+19.*pi*p)
C
C To calculate the surface area that will be removed from the grain.
Sr=-12.*A+gl*(w+p/2.)*(4.*sqrt3-2.*pi)
C
C To calculate the initial surface area for the 19hex grian with rounded
C comers.
So= Ss-Sr
C
C To calculate the unknows of the grain under the condition x.le.5*w.
if(O.le.x.and.x.le.w/2.) then
A=sqrt3/3.*(w-2.*x+(p+2.*x)/2.)**2-pi/6.*
&(w-2.*x+(p+2.*x)/2.)**2
C To calculate the volume that will be removed from the sharp comer grain.
Vr=6.*A*(gl-2.*x)
C To calculate the volume for the sharp comer grain at some distance burnt.
Vn=.25*(gl-2.*x)*(2.*sqrt3*(f-2.*x)**2.
& -19.*pi*(p+2.*x)**2.)
C
C To calculate the volume for the 19hex grain with rounded comers.
V= Vn-Vr
C
C To calculate the surface area that will be removed from the sharp
C comer grain.
Sr=12.*A+(gl-2.*x)*(w-2*x+(p+2.*x)/2.)*(4.*sqrt3-2.*pi)
C To calculate the surface area for the sharp comer grain.
Sn=2.*V/(gl-2.*x)+(gl-2.*x)*(2.*sqrt3*(f-2.*x)+
& 19.*pi*(p+2.*x))
C
C To calculate the surface area for 19hex grain with rounded comers.
S= Sn-Sr
C To calculate the mass fraction.
frac= 1-V/Vo
dsdx=-8.*sqrt3*(f-2.*x)-76.*pi*(p+2.*x)+(gl-2.*x)*(-4.*sqrt3+38.
"& *pi)+ 16*sqrt3*(w+p/2.-x)-8.*pi*(w+p/2.-x)+(gl-2.*x)*
"& (4.*sqrt3-2.*pi)
surf=S
vol=V
dsdxsl(kpr)=dsdx
fracsl(kpr)=frac
surfsl(kpr)=surf
return
endif
C
C Due to the cross section at the sliver point x=.5*w there will be 24
C identical inner slivers,12 identical side slivers. After slivering the
C surface area and the volume function become more complex. Each type of
C sliver will be treated seperately and later the volumes will be combined
C to complete the function.
C
C To calculate the 12 identical side slivers for the grain x=.5/w.
nsl=l
coef(kpr)=O.
if(igrad.eq.l.or.igrad.eq.2)go to 729
117
if(nslp(kpr).eq.l)goto 730
tsl(kpr)=yar(3)
ts(kpr)=w/2.*(- 1 .+(pbrch4,mean)**abr(kpr))/
& (bbr(kpr)*(pbase*lI.e-6)**abr(kpr))
730 continue
coef(kpr)=(ts(kpr)+tsl(kpr)-(deltat+yar(3)))/ts(kpr)
if(coef(kpr).gt.l1.)coef(kpr)= 1.
if(coef(kpr).lt.0.)coef(kpr)=0.
729 if(w/2.lt.x.and.x.lt.XlI.and.x.lt.X2) then
C
C To calculate the areas of the grain.
A=sqrt3/3.*(w-2.*x+(p+2.*x)t2.)**2-pi/6.*
&(w-2.*x+(p+2.*x)t2.)**2
theta--acos(dpcl(p+2.*x))
Al=theta/4.*(p+2.*x)**2-dpc/4.*sqrt((p+2.*x)**2-dpc**2)
omega--acos(2.*dpc/(p+2.*x)- 1.)
A2=O0.125*(p+i2.*x)*((p+2.*x)*(omega+sin(omega))-2.*dpc*sin(omega))
C To calculate the volumes of the grain.
Vl=3.*(gl-2.*x)*(2.*sqrt3*dpc**2-pi*(p+2.*x)**2+24.*A1)
V2=6.*(gl-2.*x)*(2.*dpc**2-dpc*(p+2.*x)-pi/4.*(p+2.*x)**2
&+2.*A1+4.*A2)
C To calculate the surface areas of the grain.
S 1=2.*Vl/(gl-2.*x)+ 12.*(gl-2.*x)*(pi-6.*theta)*(p+2.*x)
S2=2.*V2/(gl-2.*x)+12.*(gl-2.*x)*(dpc+i(p+2.*x)*(pi/2.-omega
&-theta-sin(omega)))
C To calculate the total volume and total surface area.
Vf=Vl+V2
Sf=S1+S2
C To calculate the mass fraction.
frac=1I.- Vf/o
surf=Sf
dsdx---surf/(x2-x)
vol=Vf
dsdx=coef(kpr)*dsdxsl(kpr)+( 1.-coef(kpr))*dsdx
dsdxsl(kpr)=dsdx
frac=coef(kpr)*fracsl(kpr)+(l.-coef(kpr))*frac
fracsl(kpr)=-frac
surf=coef(kpr)*surfsl(kpr)+( .-coef(kpr))*surf
surfsl(kpr)--surf
return
endif
if(x.gLXl .and.x.lt.X2)then
C To calculate the area of the grain.
A=sqrt3/3.*(w-2.*x+(P+2.*x)f2.)**2-pi/6.*
&(w-2.*x+(p+2.*x)/2.)**2
theta=acos(dpc/(p+2.*x))
Al=theta/4.*(p+2.*x)**2-dpc/4.*sqrt((p+2.*x)**2-dpc**2)
omega--acos(2.*dpc/(p+2.*x)-l.)
A2=0. l
2 5*(p+2.*x)*((p+2.*x)*(omega+sin(omega))-2.*dpc*sin(omega))
C To calculate the volume of the grain.
V2=6.*(gl-2.*x)*(2.*dpc**2-dpc*(p+2.*x)-pi/4.*(P+2.*x)**2
&+2.*Al+4.*A2)
C To calculate the surface area of the grain.
S2=2.*V2/(gI.2.*x)+ l2.*(gl-2.*x)*(dpc+(p+2.*x)*(pi/2.-omega
118
&-tlieta-sin(omega)))
C To calculate the volume and the surface area.
Vf=V2
Sf=S2
C To calculate the the mass fraction.
frac=l-Vf/Vo
surf.-Sf
dsdx=--surfl(x2-x)
vol=Vf
dsdx=coef(kpr)*dsdxsl(kpr)+(l .- maf(kpr))*dsdx
dsdxsl(kpr)=dsdx
*ft-ac -,coe f(kp r) *frac sl(kpr)+( 1. -coe f(kpr)) *frac
fracsl(kpr)=frac
surf=coef(kprY)'surfsl(kpr)+(l .-coef(kpr))*surf
surfsl(kpr)=surf
return
endif
if(x.gLX2)then
dsdx=O.
surf=-O.
vol=O.
hu -- fracsl(kpr) *coef(kpr)+ 1. -coef(k pr)
fracl(kpr)=frac
if(frac.gt..9999) frac=l.
if(frac.gt..9999)Tetum
dsdx=O.
dsdx=dsdxsl(kpr)*coef(kpr)
dsdxsl(kpr)=dsdx
if(abs(dsdx).lL I.)dsdx=O.
surff-surfsl(kpr)*coef(kpr)
surfsl(kpr)=surf
return
endif
end
