c                       ***** PROGRAM MINEOS *****
c
c   This program needs one input file and creates two output files.
c   It also needs a few screen inputs for program control.
c   The input file contains the spherical reference model.
c   The output files are
c
c        1) an ASCII file with model listing + summaries of the
c           properties of the normal modes, and
c        2) a binary file with model listing and the summaries and
c           eigenfunctions of the normal modes.
c
c   The program also needs a few screen inputs for program control:
c
c       name:  name of the input reference model file
c       namo1: name of the ASCII output file
c       namo2: name of the binary output file
c       jcom,eps,wgrav
c           jcom=1 radial modes, 2 toroidal modes, 3 spheroidal modes,
c           4 inner core toroidal modes.
c           eps controls the accuracy of the Runge-Kutta integration scheme.
c           The relative accuracy of an eigenfrequency will be 2-3 x eps.
c           It also controls the precision with which a root is found and
c           the minimum relative separation of two roots with the same angular
c           order. It is safe to set eps=1.d-7.
c           wgrav is the frequency in millihertz above which gravitational terms
c           are neglected -- this gives about a factor of 3 increase in speed.
c           (L.Z. -- However, this can change the eigenfrequencies and
c           eigenfunctions of the low-frequency normal modes. With the powerful
c           workstations nowadays, wgrav should always be set large enough,
c           e.g. 1000., so that gravity is always effectively accounted for.
c       lmin,lmax,wmin,wmax,nbran
c           lmin - lmax defines the range of angular orders to be computed.
c           if jcom=1 this is ignored. wmin - wmax defines the frequency range
c           to be computed (in millihertz).
c           nbran specifies the number of mode branches to be computed (0 will
c           cause all branches to be computed, otherwise the nbran lowest
c           frequency branches above wmin will be computed)
c
c   Structure of the input reference model file:
c
c       Line 1: title (80 chars)                (20a4 format)
c       Line 2: ifanis,tref,ifdeck              (unformatted)
c           ifanis=1 for anisotropic model, 0 for isotropic
c           tref=ref period(secs) of model for dispersion correction.
c           if tref is .le. 0. no correction is made.
c           ifdeck=1 for card deck model, 0 for polynomial model.
c
c             *** card deck model ***
c
c       Line 3: n,nic,noc                       (unformatted)
c           n=no of levels,nic=index of solid side of icb,noc=index of
c           fluid side of mcb.note that n must be .le. 223.
c       Line 4-n: r,rho,vpv,vsv,qkappa,qshear,vph,vsh,eta
c               (f8.0,3f9.2,2f9.1,2f9.2,f9.5)
c           if isotropic vp=vpv,vs=vsv and vph,vsh,eta are unspecified.
c           if the q model is not specified then no disp. correction is made.
c           S.I. units are used,e.g. radius in metres and velocities in m/s.
c
c             *** polynomial model ***
c
c       Line 3: nreg,nic,noc,rx                 (unformatted)
c           nreg is the number of regions in the model,nic and noc as before
c           and rx is the normalising radius for the polynomials.
c           rx is given in kms and is usually 6371.
c       Line 4:   nlay,r1,r2                      (unformatted)
c           nlay is the number of levels to be used in the region extending
c           from radius r1 to r2 (in kms).
c       Line 5-9 (iso) or Line 5-12 (ani): coefs     (5f9.5)
c           5 sets of coefficients are required for each region of an isotropic
c           model and 8 sets for an anisotropic model. Each polynomial can be
c           up to a quartic in r/rx (5 coefficients) and are ordered thusly:
c           rho,vpv,vsv,qkappa,qshear,vph,vsh,eta. The coeffs are given in the
c           usual mixed seismological units (rho in g/cc, vel in km/s etc.)
c           conversion to S.I. is done by the program.
c           Lines 4-9 (iso) or 4-12 (ani) are repeated for each region of the
c           model
c
c   Structure of the ASCII output file:
c
c       This file lists the reference model and the normal-mode properties, i.e.
c       frequency in rad/sec and millihz, period in sec, group velocity in km/s,
c       Q and a parameter which is the ratio of kinetic to potential energy
c       minus 1. This parameter should be small (of order eps) if the
c       eigenfunction is accurate and if there are enough radial knots in the
c       reference model specification to allow quodratures to be done accurately.
c       You will probably see some degradation in this parameter for strongly
c       exponential modes such as Stoneley modes as well as higher-frequency
c       fundamental modes.
c
c   Structure of the binary output file
c     
c       This is a sequential access binary file. The first record contains the
c       screen-input controlling parameters, followed by the reference model
c       list. Then the property summaries and eigenfunctions of the normal modes
c       are stored in the same order as they are in the ASCII output file. The
c       writing statement for each mode is
c
c           write(ioeig) (abuf(i),i=1,nvec)
c
c       where abuf is a real*4 array and nvec is 8+6*n for spheroidal modes,
c       8+2*(n-noc) for toroidal modes, 8+2*n for radial modes, and 8+2*nic for
c       inner-core toroidal modes. The first eight 4-byte slots of abuf are n,l,
c       both of 4-byte integers,  frequecy, Q and group velocity, all of 8-byte
c       floating-point numbers. The rest of abuf contains w(1...n) and
c       dw/dr(1...n) for toroidal modes and u(1...np), du/dr(1...n), v(1...n),
c       dv/dr(1...n), phi(1...n) and dphi/dr(1...n) for spheroidal modes, where
c       phi is the coefficient in the spherical harmonic expansion for the
c       change in gravitational potential due to the free oscillations (see
c       Dahlen & Tromp 1998, eq.8.18). The eigenfunctions are normalized in
c       such a way that

c           frequency**2 times integral:[rho*w*w*r*r] dr and
c           frequency**2 times integral:[rho*(u*u+v*v)*r*r] dr

c       are 1 for toroidal and spheroical modes, respectively. The eigenfunctions
c       are defined by the spherical harmonic expension, eq.(C.135) in Dahlen &
c       Tromp (1998) with the spherical harmonics Y_lm defined therein.
c
c   The computation is done using a non-dimensionalized system in which the
c   non-dimensionalization constants (dnorm) are 5515 mg/m**3 for density, 
c   radius of the earth (rn=6371000) for length or radius, and
c   sqrt(pi*dnorm*g) for frequency or 1/sqrt(pi*dnorm*g) for time, where g is
c   Newtons gravity constant (6.6723e-11). These normalizations result in
c
c       acceleration non-dimensionalization = pi*dnorm*g*rn = an
c       velocity non-dimensionalization     = rn*sqrt(pi*dnorm*g) = vn
c       
c   H.Y. Yang 2008
c   In order to reconcile the radial output format with spheroidal one, 
c   the other four eigenfunctions (V=0,dV/dr=0,phi=[eqn.8.148 in D.&T., 1998],
c   dphi/dr=-4*rho(r)*U(r))[here is normalized density]
c   are calculated and write out. [phi is normalized by "pigrho"=pi*G*dnorm 
c   which is inconsistent with U and V]
c   Phase velocity correction: l+0.5 --> sfl3
c--------------------------------------------------------------------------

	character*128 name,namo1,namo2
	real*8 eps,wgrav,wmin,wmax,wwminn
     
	common/counter/wwminn,ijkl,ijk
	common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
	
	data iin,iout,ioeig/7,8,9/

c   Screen-input parameters for program control.
      
	write(6,*) 'Enter spherical-earth reference model file name:'
	call flush(6)
	read(*,'(a)') name
      
	write(6,*) 'Output file name for modal summaries (ASCII)'
	call flush(6)
	read(*,'(a)') namo1
      
	write(6,*) 'Output file name for modal eigenfunctions (binary)'
	call flush(6)
	read(*,'(a)') namo2
     
10	write(6,*) 'Enter jcom (1=rad; 2=tor; 3=sph; 4=ictor), eps and wgrav'
	call flush(6)
	read(*,*) jcom,eps,wgrav
	if((jcom.ge.1).and.(jcom.le.4)) then
	else
	  write(6,*) 'Incorrect mode-type specification: must be 1<=jcom<=4.'
	  call flush(6)
	  goto 10
	endif

	write(6,*) 'Enter lmin,lmax,wmin,wmax,nbran'
	call flush(6)
	read(*,*) lmin,lmax,wmino,wmaxo,nbran

c   Initialize the total count of the number of modes.
	
	iter=0
11	ijk=0
	iter=iter+1
	
	wmin=wmino
	wmax=wmaxo
	open(iout,file=namo1,form='formatted',status='unknown')
	open(ioeig,file=namo2,form='unformatted',access='sequential')
 	write(ioeig) jcom,real(wmin),real(wmax),lmin,lmax,real(wgrav)
	call flush(ioeig)
	
c   Read in the reference earth model parameters and write them in the two
c   output files.
         
	call model(iin,name,iout,ioeig,jcom)     

c   Start mode calculations.
     
	call wtable(iout,ioeig,jcom,eps,wgrav,lmin,lmax,wmin,wmax,nmode,nbran)
	
     	if(ijk.gt.0) then
	  close(iout)
	  close(ioeig)
	  write(6,'(/i6,a/)') ijk,' modes in total.'
	  call flush(6)
	else
	  close(iout,status='delete')
	  close(ioeig,status='delete')
c	  write(6,'(/a/)') '  No modes found. Output files deleted.'
c	  call flush(6)
	endif
	
c   In case no modes are found while nmode > 0, adjust wmin and recalculate.

	if((ijk.eq.0).and.(nmode.gt.0)) then
	  wmino=wmino+(wmaxo-wmino)/100.d0
	  if(iter.le.50) then
	    write(6,'(/a/)') '  No modes found. Adjust wmin and recalculate.'
	    call flush(6)
	    goto 11
	  else
	    write(6,'(/a/)') '  No modes found after 50 iterations. Output files deleted.'
	    call flush(6)
	  endif
	endif
	
	stop
	end
c
c
	subroutine wtable(iout,ioeig,jjcom,epss,wgravs,lmin,lmax,wmin,wmax,nmode,nbran)
	
c   Makes up table of frequencies.

	include'mineos.inc'
	implicit real*8(a-h,o-z)
	dimension wt(level),fmxx(11),nmxx(11)
	
	common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     1	  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
	common/shanks/b(46),c(10),step(8)
	common/shank/dx,stepf,maxo,in
	common/mtab/we(lev1),de(lev1),ke(lev1),wtry(lev5),bm(lev5),um(lev5)
	common/counter/wwminn,ijkl,ijk

	data inss/5/     
	data fmxx/52.d0,82.d0,102.d0,202.d0,252.d0,502.d0,1002.d0,1502.d0,
     1	  2002.d0,2502.d0,3002.d0/
c	data nmxx/200,300,350,650,800,1600,3100,4500,5900,7500,9000/
	data nmxx/400,600,800,1600,3100,4500,5900,7500,9000,11000,13000/

	nmx=10000
	do ii=11,1,-1
	  if(wmax.le.fmxx(ii)) nmx=nmxx(ii)
	end do
	
	jcom=jjcom
	cmhz=pi/500.d0
	stepf=1.d0
	eps=epss
	eps1=eps
	eps2=eps
	wgrav=wgravs
	wgrav=wgrav*cmhz
	
c	write(iout,100) eps,eps1,wgrav
c	call flush(iout)
100	format(/,'Integration precision =',g12.4,'  Root precision =',
     1	  g12.4,'  Gravity cut off =',g12.4,' rad/s',///,6x,'mode',
     2	  8x,'w(rad/s)',7x,'w(mhz)',10x,'t(secs)',6x,'grp vel(km/s)',
     3	  8x,'q',13x,'Raylquo',/)
     
	call steps(eps)
	wmin=wmin*cmhz
	wmax=wmax*cmhz
	wwminn=wmin
	wmin0=wmin
	if(nbran.le.0.or.nbran.gt.nmx) nbran=nmx
	if(lmin.le.0) lmin=1
	if(jcom.ne.1) goto 6
	lmin=0
	lmax=0
6	l=lmin-1
10	l=l+1
	i1i2=0
	i11=-100
	i22=-100

c   The following five lines are moved here from the beginning of the subroutine.
c   They seem to belong here to initialize the arrays at the start of computation
c   for each l.

	do i=1,lev5
	  bm(i)=0.d0
	  wtry(i)=0.d0
	end do
	nev=2
	
c   The following few lines are added to make the minimum frequency for each l
c   adjustable according to the lowest computed eigenfrequency for the previous
c   value of l.

c	if(ijk.le.0) then
c	  wt(1)=wwminn
c	else
c	  if((wwminn/cmhz).le.2.d0) then
c	  wt(1)=wwminn/3.d0
c	  elseif((wwminn/cmhz).le.10.d0) then
c	  wt(1)=wwminn-1.d0*cmhz
c	  elseif((wwminn/cmhz).le.50.d0) then
c	  wt(1)=wwminn-2.d0*cmhz
c	  else
c	  wt(1)=wwminn-3.d0*cmhz
c	  endif
c	endif
	wt(1)=wmin0
	wt(2)=wmax
	
	if(l.gt.lmax.or.wt(1).ge.wmax) return
	ijkl=0
	knsw=1
	maxo=inss
	fl=l
	fl1=fl+1.d0
	fl2=fl+fl1
	fl3=fl*fl1
	sfl3=dsqrt(fl3)
	
c   Determine mode count.

	we(1)=wt(1)
	call detqn(we(1),ke(1),de(1),0,ratiomax)
	imax=2*nmx
	do i=2,imax
	  we(i)=wmax
	  ke(i)=-10
	end do
	do i=2,nev
	  call entry(wt(i),imax,kei)
	  nmode=kei-ke(1)
	  kkk=i
	  if(nmode.ge.nbran) goto 25
	end do
	
25	write(6,900) l,wt(1)/cmhz,wt(kkk)/cmhz,nmode,nbran
	call flush(6)
900	format(/' Degree l=',i6,': count between ',f8.3,
     1	  ' mHz and ',f8.3,' mHz =',i5,'. MAX=',i6/)
     
	if(nmode.le.0) goto 10
c	if(nmode.le.0) return
	if(nmode.gt.nbran) nmode=nbran
	imax=2*nmode

c   Fill up table using bisection.

	indx=2

c   ke(1)+1 is the branck number of the first mode. ke(2*nmode)+1  
c   is the branch number of the last (nmode-th) mode.
c   ichk below is the branch number of the first mode.

	ichk=ke(1)+1
35	if(ke(indx).ne.ichk) goto 40
	indx=indx+2
	if(indx.gt.imax) goto 60
	ichk=ichk+1
	goto 35
40	i1=indx-1
45	indx=indx+2
	if(indx.ge.imax) goto 46
	ichk=ichk+1
	if(ke(indx).ne.ichk) goto 45
46	i2=min0(indx,imax)
	wtst=0.5d0*(we(i2)+we(i1))
	
c   The following few lines are added to avoid infinite loops which happens
c   when the consecutive values of i1 and i2 do not change.

	if((i11.eq.i1).and.(i22.eq.i2)) then
	  i1i2=i1i2+1
	else
	  i1i2=0
	endif
	if(i1i2.gt.100) then
	  i1i2=0
	  goto 50
	endif
	i11=i1
	i22=i2
	
	if((we(i2)-we(i1))/wtst.lt.eps2) goto 50
	call entry(wtst,imax,ktst)
	indx=i1+1
	ichk=ke(i1)+1
	goto 35
50	continue
c50	write(*,901) we(i1),ke(i1),we(i2),ke(i2)
901	format('Problem in table: ',2(g16.8,i5))
	j1=i1+1
	j2=i2-1
	do i=j1,j2
	  de(i)=1.d0
	end do
	indx=i2+2
	if(indx.ge.imax) goto 60
	ichk=ke(i2)+1
	goto 35
60	nev=imax
c	write(*,902) (i,we(i),ke(i),de(i),i=1,nev)
902	format(i5,1pd17.7,i6,1pd17.7)

c   Find roots.

	knsw=0
	maxo=8
	call rotspl(nev,eps1,wmin,wmax,wt,iout,ioeig)
	goto 10
	
	return
	end
c
c
	subroutine entry(w,imax,kei)
      
	include'mineos.inc'
	implicit real*8(a-h,o-z)
	
	common/mtab/we(lev1),de(lev1),ke(lev1),wtry(lev5),bm(lev5),um(lev5)
	
	call detqn(w,kei,dei,0,ratiomax)
	indx=min0(max0(2*(kei-ke(1)),1),imax)
	if(indx.eq.1.and.we(1).lt.w) goto 10
	if(indx.eq.imax.and.we(imax).gt.w) goto 10
	if(kei.ne.ke(indx)) goto 5
	if(we(indx).gt.w) goto 10
	indx=indx+1
	if(we(indx).lt.w) goto 10
	return
5	we(indx)=w
	ke(indx)=kei
	de(indx)=dei
	indx=indx+1
10	we(indx)=w
	ke(indx)=kei
	de(indx)=dei
	
	return
	end
c
c
	subroutine rotspl(nev,eps1,wmin,wmax,wt,iout,ioeig)
	
c   Find roots by spline interpolation.

	include'mineos.inc'
	implicit real*8(a-h,o-z)
      
	dimension x(60),det(60),qx(3,60),wrk(60),wt(1),ichar(4)
	
	common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     1	  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
	common/mtab/we(lev1),de(lev1),ke(lev1),wtry(lev5),bm(lev5),um(lev5)
	common/counter/wwminn,ijkl,ijk

cxxxxxxx tol was 1.d-9, itmax was 15.
	
	data tol,itmax /1.d-18,100/
        CHARACTER*4 :: ichar
        ichar = ( /' S',' T',' S',' C'/)
	
	knev=1
	nmo=nev/2
	do 100 i=1,nmo
	k1=2*i
	k=k1-1
	if(de(k)*de(k1).gt.0.d0) goto 100
	nord=ke(k1)
	if(l.eq.1) nord=nord+1
	nind=nord+1

c   The following line is added to avoid numerical problem when the value
c   of nind is less than 1, which happens when ke(k1) is less than 0.

	if(nind.le.0) goto 100
		
	det(1)=de(k)
	det(2)=de(k1)
	x(1)=we(k)
	x(2)=we(k1)
	
c	if(wtry(nind).le.x(1)) goto 5
c	if(wtry(nind).ge.x(2)) goto 5
c	call detqn(wtry(nind),knt,ftry,0,ratiomax)
c	if(ftry*det(1).lt.0.d0) goto 5
c	x(1)=wtry(nind)
c	det(1)=ftry
c5	c=x(1)
c	if(dabs(det(1)).gt.dabs(det(2))) c=x(2)
c	write(*,910) x(1),det(1),x(2),det(2),c
910	format(6g18.10)

	iter=0
10	continue

c	b=0.5d0*(x(1)+x(2))

	b=x(1)+(0.d0-det(1))*(x(2)-x(1))/(det(2)-det(1))
	iter=iter+1

c	j=1
c	m=2
c	ntry=2
c15	t=dabs(b*eps1)
c	if(dabs(b-c).lt.t) goto 65
c	call detqn(b,knt,fb,0,ratiomax)
c	ind=1
c	do 20 m=2,ntry
c	ind=ind+1
c	if(b.lt.x(m)) goto 25
c20	continue
c25	ntry=ntry+1
c	j2=ntry
c30	j1=j2-1
c	x(j2)=x(j1)
c	det(j2)=det(j1)
c	if(j1.eq.ind) goto 35
c	j2=j2-1
c	goto 30
c35	continue
c	x(ind)=b
c	det(ind)=fb
c	idn=0
c	do 40 m=2,ntry
c	idn=idn+1
c	iup=idn+1
c40	if(det(idn)*det(iup).le.0.d0) goto 45
c45	ind=iup
c	if(dabs(det(idn)).lt.dabs(det(iup))) ind=idn
c	c=x(ind)
c	if(ntry.ge.itmax) goto 60
c	call dsplin(ntry,x,det,qx,wrk)
c	del=-det(ind)/qx(1,ind)
c50	delx=-det(ind)/(qx(1,ind)+del*qx(2,ind))
c	if(dabs(delx-del).lt.tol) goto 55
c	if(del*delx.lt.0.d0) goto 60
c	del=delx
c	goto 50
c55	b=c+delx
c	if(b.ge.x(idn).and.b.le.x(iup)) goto 15
c60	continue
c	x(1)=x(idn)
c	x(2)=x(iup)
c	det(1)=det(idn)
c	det(2)=det(iup)
c	goto 10
	
	
c   Write out frequencies.

65	continue

	call detqn(b,knt,fb,1,ratiomax)
c	write(*,900) b,fb,x(1),det(1),x(2),det(2)
900	format(6g20.12)

	if((fb*det(1)).gt.0.d0) then
	  x(1)=b
	  det(1)=fb
	else
	  x(2)=b
	  det(2)=fb
	endif
			
c	write(*,900) b,fb,dabs(x(2)-x(1)),tol
	if((dabs(x(2)-x(1)).ge.tol).and.(iter.le.itmax)) goto 10
c	if((dabs(x(2)-x(1)).ge.tol).and.(iter.le.itmax)) then
c	  write(*,*) '1 stop -- 2 continue'
c	  read(*,*) istop
c	  if(istop.eq.2) goto 10
c	endif
	
	wdiff=(b-wray*wn)/b
	tcom=2.d0*pi/b
	wmhz=1000.d0/tcom
	gcom=vn*cg/1000.d0
	qmod=0.d0
	if(qinv.gt.0.d0) qmod=1.d0/qinv

c   The maximum value of the group velocity is forced to be 999.999999. The
c   unreasonably large values for group velocity usually happens to the
c   inner-core shear modes which are of no practical interests in seismology.
	
	if(gcom.ge.1000.d0) gcom=999.999999d0
	write(6,201) nord,ichar(jcom),l,wmhz,tcom,gcom,qmod,wdiff,ratiomax
	call flush(6)
c	write(*,202) nord,ichar(jcom),l,b,qmod,wmhz,tcom,gcom,wdiff
	write(iout,200) nord,ichar(jcom),l,b,wmhz,tcom,qmod,gcom,wdiff,ratiomax
	call flush(iout)
	ijkl=ijkl+1
	if(ijkl.eq.1) wwminn=b
200	format(i5,a2,i5,7g16.7)
201	format(i5,a2,1x,i6,1x,f20.15,1x,f11.5,1x,f10.6,1x,f13.6,1x,e11.4,1x,e11.4)
202	format(i5,a2,i5,2(1x,f20.15),4g16.7)

	call modout(b,qmod,gcom,ratiomax,ioeig)
	
	if(nord.eq.0) wmin=b
	if(bm(nind).le.0.d0) goto 70
	blin=b+cg*wn
	bp2=5.d0*bm(nind)-4.d0*b+2.d0*(um(nind)+2.d0*cg*wn)
	bdiff=dmax1(dabs(blin-bp2),eps1*bp2)
	bup=bp2+bdiff
	bdwn=bp2-bdiff
	wtry(nind)=bdwn
	goto 71
70	bup=b+2.d0*cg*wn
	wtry(nind)=b
71	bm(nind)=b
	um(nind)=cg*wn
	if(bup.ge.wmax) goto 100
	knev=knev+1
	wt(knev)=bup
100	continue
	nev=knev+1
	wt(1)=wmin
	wt(nev)=wmax
	
	return
	end
c
c
	subroutine modout(wcom,qmod,gcom,ratiomax,ioeig)
	
	include'mineos.inc'
	implicit real*8(a-h,o-z)
	integer krec(ndata)
	dimension abuf(lev3)
	
	common/counter/wwminn,ijkl,ijk
	common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     1	  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
	common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
	common/eifx/a(14,ndata),ndum(ndata)
	common/eigsav/krec,nrec

        common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     1    qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     2    qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     3    fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     4    nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     5    aspl(3,ndata)
        REAL*8 pU(ndata),SO(ndata),dS,sumall
        REAL rrec(ndata)
        INTEGER j1
      
	if(jcom.eq.1) then
	  nvec=7+6*nrec
	elseif(jcom.eq.2) then
	  nvec=7+2*nrec
	elseif(jcom.eq.3) then
	  nvec=7+6*nrec
	elseif(jcom.eq.4) then
	  nvec=7+2*nrec
	else
	endif
	
	if(jcom.eq.2) then
	  do i=1,noc
	    a(1,i)=0.d0
	    a(2,i)=0.d0
	  end do
	elseif(jcom.eq.4) then
	  do i=nic+1,n
	    a(1,i)=0.d0
	    a(2,i)=0.d0
	  end do
	else
	endif

	abuf(1)=dble(nord)
	abuf(2)=dble(l)
	abuf(3)=dble(wcom)
	abuf(4)=dble(qmod)
	abuf(5)=dble(gcom)
	abuf(6)=dble(wdiff)
	abuf(7)=ratiomax
	jst=7
	j=0
	if(jcom.eq.1) then
          SO=0.d0
          DO i=1,n
            IF(krec(i).EQ.1) THEN
              j=j+1
              pU(j)=rho(i)*a(1,i)
              rrec(j)=r(i)
c           Using trapzoidal method to get integrand
              dS = 0.d0
              IF ( j>=2 ) THEN
                j1=j-1
                dS=(pU(j)+pU(j1))*(rrec(j)-rrec(j1))/2.d0
                SO(j)=SO(j1)+dS
              ENDIF
            ENDIF
          END DO
          sumall=SO(j)
          factor=4.d0
          j=0
          do i=1,n
            if(krec(i).eq.1) then
              j=j+1
              abuf(j+jst)=a(1,i)
              abuf(j+nrec+jst)=a(2,i)
              abuf(j+2*nrec+jst)=0.d0
              abuf(j+3*nrec+jst)=0.d0
              abuf(j+4*nrec+jst)=factor*(sumall-SO(j))
              abuf(j+5*nrec+jst)=-factor*pU(j)
            endif
          end do
	elseif(jcom.eq.2) then
	  do i=noc+1,n
	    if(krec(i).eq.1) then
	      j=j+1
	      abuf(j+jst)=a(1,i)
	      abuf(j+nrec+jst)=a(2,i)
	    endif
	  end do
	elseif(jcom.eq.3) then
	  do i=1,n
	    if(krec(i).eq.1) then
	      j=j+1
	      abuf(j+jst)=a(1,i)
	      abuf(j+nrec+jst)=a(2,i)
	      abuf(j+2*nrec+jst)=a(3,i)
	      abuf(j+3*nrec+jst)=a(4,i)
	      abuf(j+4*nrec+jst)=a(5,i)
	      abuf(j+5*nrec+jst)=a(6,i)
	    endif
	  end do
	elseif(jcom.eq.4) then
	  do i=1,nic
	    if(krec(i).eq.1) then
	      j=j+1
	      abuf(j+jst)=a(1,i)
	      abuf(j+nrec+jst)=a(2,i)
	    endif
	  end do
	endif
	ijk=ijk+1
	
c	WRITING EIGENFUNCTIONSS

        write(ioeig) (abuf(i),i=1,nvec)
	call flush(ioeig)
     	ijkl=ijkl+1
	
	return
	end
c
c
	subroutine model(iin,name,iout,ioeig,jjcom)
	
	include'mineos.inc'
	implicit real*8(a-h,o-z)
	character*128 name
	integer*4 ititle(20)
	integer krec(ndata)
	real*8 lcon,ncon,lspl,nspl,wrk(8)
	real*4 rrec(ndata),rhorec(ndata),vpvrec(ndata),vsvrec(ndata),vphrec(ndata)
	real*4 vshrec(ndata),qkapparec(ndata),qshearrec(ndata),etarec(ndata)
	
	common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     1	  qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     2	  qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     3	  fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     4	  nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     5	  aspl(3,ndata)
	common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     1	  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
	common/eifx/vpv(ndata),vph(ndata),vsv(ndata),vsh(ndata),
     1	  eta(ndata),nnnn(19*ndata)
	common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
	common/eigsav/krec,nrec
	
	data bigg,tau,rhobar/6.6723d-11,1.d3,5515.d0/
	

	pi=4.d0*datan(1.d0)

	do i=1,ndata
	  krec(i)=0
	end do

 	open(iin,file=name,status='old',form='formatted')
	read(iin,'(20a4)') (ititle(i),i=1,20)
	read(iin,*) ifanis,tref,ifdeck
	if(ifdeck.eq.0) goto 1000
	
c       *** card deck model ***

	read(iin,*) n,nic,noc,lrec
c	read(iin,1055) (r(i),rho(i),vpv(i),vsv(i),
c     1    qkappa(i),qshear(i),vph(i),vsh(i),eta(i),i=1,n)
 	read(iin,*) (r(i),rho(i),vpv(i),vsv(i),
     1    qkappa(i),qshear(i),vph(i),vsh(i),eta(i),i=1,n)
c     	do i=1,nic
c     	  rho(i)=12850.d0
c	  vpv(i)=11160.d0
c	  vph(i)=11160.d0
c	  vsv(i)=3580.d0
c	  vsh(i)=3580.d0
c	end do
c     	do i=nic+1,noc
c     	  rho(i)=11000.d0
c	  vpv(i)=9150.d0
c	  vph(i)=9150.d0
c	  vsv(i)=0.d0
c	  vsh(i)=0.d0
c	end do
c     	do i=noc+1,n
c     	  rho(i)=4500.d0
c	  vpv(i)=10000.d0
c	  vph(i)=10000.d0
c	  vsv(i)=5200.d0
c	  vsh(i)=5200.d0
c	end do
    	irstart=1
	if(jjcom.eq.2) irstart=noc+1
     	if(lrec.eq.0) then
     	  do i=irstart,n
	    krec(i)=1
	  end do
	elseif(lrec.eq.1) then
	  i1=irstart
	  
c   Below 6330km radius, results are stored at 1sample/km.

	  do j=0,6330
	    rcur=real(j)*1.e3
	    xdiff=1.e20
	    do i=i1,n
	      if(xdiff.gt.abs(rcur-real(r(i)))) then
	        xdiff=abs(rcur-real(r(i)))
	        ki=i
	      endif
	    end do
	    krec(ki)=1
	    i1=ki
	  end do
	  
c   Both sides of a discontinuity are always saved. 

	  do i=irstart,n-1
	    if(r(i).eq.r(i+1)) then
	      krec(i)=1
	      krec(i+1)=1
	    endif
	  end do
	  krec(irstart)=1
	  krec(n)=1

c   Above 6330km radius, all nodes are saved. 

	  do i=irstart,n
	    if(r(i).ge.6330000.) krec(i)=1
	  end do
	  
	endif
105	format(f8.0,3f9.2,2f9.1,2f9.2,f9.5,1x,i1)
1055	format(f8.0,3f9.2,2f9.1,2f9.2,f9.5)
	goto 2000
	
c       *** polynomial model ***

1000	read(iin,*) nreg,nic,noc,rx
	rx=rx*tau
	n=0
	knt=0
	jj=5
	if(ifanis.ne.0) jj=8
	do 10 nn=1,nreg
	read(iin,*) nlay,r1,r2,ndec
	r1=r1*tau
	r2=r2*tau
	dr=(r2-r1)/dble(nlay-1)
	do i=1,nlay
	  n=n+1
	  r(n)=r1+dr*dble(i-1)
	  if(mod(i-1,ndec).eq.0) then
	    krec(n)=1
	  endif
	end do
	do j=1,jj
	  read(iin,'(5f9.5)') (wrk(i),i=1,5)
	  do i=1,nlay
	    ind=knt+i
	    rt=r(ind)/rx
	    val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
	    if(j.eq.1) rho(ind)=val*tau
	    if(j.eq.2) vpv(ind)=val*tau
	    if(j.eq.3) vsv(ind)=val*tau
	    if(j.eq.4) qkappa(ind)=val
	    if(j.eq.5) qshear(ind)=val
	    if(ifanis.eq.0) then
	    else
	      if(j.eq.6) vph(ind)=val*tau
	      if(j.eq.7) vsh(ind)=val*tau
	      if(j.eq.8) eta(ind)=val
	    endif
	  end do
	end do
	knt=knt+nlay
	if(mod(nlay-1,ndec).eq.0) goto 10
	krec(n)=1
10	continue
2000	if(ifanis.ne.0) goto 3000
	do i=1,n
	  vph(i)=vpv(i)
	  vsh(i)=vsv(i)
	  eta(i)=1.d0
	end do
3000	if(iout.lt.0) goto 30

c   Write model to the two output files.

	j=0
	nrecic=0
	nrecoc=0
	do i=1,n
	  if(krec(i).eq.1) then
	    j=j+1
	    rrec(j)=real(r(i))
	    rhorec(j)=real(rho(i))
	    vpvrec(j)=real(vpv(i))
	    vsvrec(j)=real(vsv(i))
	    vphrec(j)=real(vph(i))
	    vshrec(j)=real(vsh(i))
	    qkapparec(j)=real(qkappa(i))
	    qshearrec(j)=real(qshear(i))
	    etarec(j)=real(eta(i))
	    if(i.eq.nic) nrecic=j
	    if(i.eq.noc) nrecoc=j
	  endif
	end do
	nrec=j
	if(jjcom.eq.1) then
	  nvec=7+2*nrec
	elseif(jjcom.eq.2) then
	  nvec=7+2*nrec
	elseif(jjcom.eq.3) then
	  nvec=7+6*nrec
	elseif(jjcom.eq.4) then
	  nvec=7+2*nrec
	else
	endif
	xmem=real(nvec)*8./1024.
	
	write(6,'(/a,i5,a,f8.3,a)') ' Eigenfunctions are saved for ',nrec,
     1	  ' radial nodes.',xmem,' KB for each mode.'
     	call flush(6)

	write(ioeig) nrec,nrecic,nrecoc,ifanis,real(tref),(rrec(i),i=1,nrec),
     1	  (rhorec(i),i=1,nrec),(vpvrec(i),i=1,nrec),
     2	  (vsvrec(i),i=1,nrec),(vphrec(i),i=1,nrec),
     3	  (vshrec(i),i=1,nrec),(qkapparec(i),i=1,nrec),
     4	  (qshearrec(i),i=1,nrec),(etarec(i),i=1,nrec)
     	call flush(ioeig)
     
c	write(iout,900) (ititle(k),k=1,20),tref
c     	call flush(iout)
c	write(iout,905) (i,r(i),rho(i),vpv(i),vph(i),vsv(i),vsh(i),
c     1	  eta(i),qshear(i),qkappa(i),krec(i),i=1,n)
c     	call flush(iout)
900	format(1x,20a4,' ref per =',f6.1,' secs',///,2x,'level',
     1	  4x,'radius',8x,'rho',9x,'vpv',9x,'vph',9x,'vsv',
     2	  9x,'vsh',9x,'eta',9x,'qmu ',8x,'qkap',/)
905	format(i6,f12.1,5f12.2,f12.5,2f12.2,1x,i1)

c   Non-dimensionalize and TI parameters.

30	rn=r(n)
	gn=pi*bigg*rhobar*rn
	vn2=gn*rn
	vn=dsqrt(vn2)
	wn=vn/rn
	do i=1,n
	  r(i)=r(i)/rn
	  if(i.gt.1) then
	    if(dabs(r(i)-r(i-1)).lt.1.d-7) r(i)=r(i-1)
	  endif
	  if(qshear(i).gt.0.d0) qshear(i)=1.d0/qshear(i)
	  if(qkappa(i).gt.0.d0) qkappa(i)=1.d0/qkappa(i)
	  rho(i)=rho(i)/rhobar
	  acon(i)=rho(i)*vph(i)*vph(i)/vn2
	  ccon(i)=rho(i)*vpv(i)*vpv(i)/vn2
	  lcon(i)=rho(i)*vsv(i)*vsv(i)/vn2
	  ncon(i)=rho(i)*vsh(i)*vsh(i)/vn2
	  fcon(i)=eta(i)*(acon(i)-2.d0*lcon(i))
	  fmu(i)=(acon(i)+ccon(i)-2.d0*fcon(i)+5.d0*ncon(i)+
     1	    6.d0*lcon(i))/15.d0
	  flam(i)=(4.d0*(acon(i)+fcon(i)-ncon(i))+ccon(i))/9.d0
     1	    -2.d0*fmu(i)/3.d0
	  rat=4.d0*fmu(i)/(3.d0*(flam(i)+2.d0*fmu(i)))
	  xlam(i)=((1.d0-rat)*qkappa(i)-.5d0*rat*qshear(i))/(1.d0-1.5d0*rat)
	  xa2(i)=(1.d0-rat)*qkappa(i)+rat*qshear(i)
	end do
	call drspln(1,n,r,rho,qro,wrk)
	
c   Compute gravity functions and spline.

	call grav(g,rho,qro,r,n)
	call drspln(1,n,r,g,qg,wrk)
	call drspln(1,n,r,fcon,fspl,wrk)
	call drspln(1,n,r,lcon,lspl,wrk)
	if(ifanis.eq.0) goto 60
	call drspln(1,n,r,acon,aspl,wrk)
	call drspln(1,n,r,ccon,cspl,wrk)
	call drspln(1,n,r,ncon,nspl,wrk)
60	nsl=n
	if(vsv(nsl).gt.0.d0) goto 70
65	nsl=nsl-1
	if(vsv(nsl).le.0.d0) goto 65
70	nicp1=nic+1
	nocp1=noc+1
	nslp1=nsl+1
	tref=0.5d0*tref/pi
	
	close(iin)
	
	return
	end
c
c
      subroutine detqn(wdim,knt,det,ifeif,ratiomax)
      
c**** Supevises the integration of the equations,it returns the value
c**** of the secular determinant as det and the count of zero crossings.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),fspl(3,ndata),
     +lcon(ndata),lspl(3,ndata),ncon(ndata),nspl(3,ndata),
     +ccon(ndata),cspl(3,ndata),acon(ndata),aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,ndata),nnnn(ndata)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension ass(14),vf(ndata),zi(4)
      
      iback=0
      w=wdim/wn
      wsq=w*w
      iexp=0
      kount=0
      kg=0
      fct=0.d0
      if(tref.gt.0.d0) fct=2.d0*dlog(tref*wdim)/pi
      goto (2,3,1,3),jcom
    1 if(wdim.le.wgrav) kg=1
      nvefm=2+kg*3
      nvesm=5+kg*9
      call sdepth(wdim,ls)
      if(ls.ge.nocp1) goto 25
      if(ls.ge.nicp1) goto 20
      if(ls.ge.2) goto 15
      r10=4.5d-4*(fl+0.5d0)/wdim
      if(r10.ge.r(2)) goto 15
      r(1)=r10
      g(1)=rho(1)*r(1)*1.333333333333333d0
c*** Form the spheroidal-mode starting solution in the solid inner core.
   15 call spsm(ls,nvesm,ass)
c*** Propagate through inner core ***
      call sprpmn(ls,nic,ass,vf,nvesm,iexp)
      r(1)=0.d0
      g(1)=0.d0
c*** Obtain the solution on the fluid side of the s/f ICB.
      call sfbm(ass,kg,iback)
   20 is=max0(ls,nicp1)
c*** If the first level ls is in the outer core, form the starting solution
c*** in the fluid outer core.
      if(is.eq.ls) call fpsm(ls,nvefm,ass)
c*** Propagate through outer core ***
	call fprpmn(is,noc,ass,vf,nvefm,iexp)
c*** Obtain the solution on the solid side of the f/s CMB.
      call fsbm(ass,kg,iback)
   25 is=max0(ls,nocp1)
c*** If the first level ls is above the CMB, form the starting solution
c*** in the fluid outer core.
      if(is.eq.ls) call spsm(ls,nvesm,ass)
c*** Propagate through mantle ***
      call sprpmn(is,nsl,ass,vf,nvesm,iexp)
      if(nsl.ne.n) goto 40
      dnorm=a(1,nsl)*a(1,nsl)
      do 26 i=2,nvesm
   26 dnorm=dnorm+a(i,nsl)*a(i,nsl)
      det=a(5,nsl)/dsqrt(dnorm)
      goto 45
   40 call sfbm(ass,kg,iback)
c*** Propagate through ocean ***
      call fprpmn(nslp1,n,ass,vf,nvefm,iexp)
      if(kg.eq.0) det=a(2,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
      if(kg.ne.0) det=a(5,n)/dsqrt(a(1,n)**2+a(2,n)**2+a(3,n)**2+
     +   a(4,n)**2+a(5,n)**2)
   45 if(ls.gt.noc) det=-det
      if(knsw.ne.1) goto 50
      if(ls.gt.noc) kount=kount-2
      irem=mod(kount,2)
      if(irem.eq.0.and.det.lt.0.d0) kount=kount+1
      if(irem.ne.0.and.det.gt.0.d0) kount=kount+1
      knt=kount
c      write(*,*) knt
   50 if(ifeif.eq.0) goto 100
c*** This does eigenfunction calculation for spheroidal modes ***
      iback=1
      jexp=0
      nbakf=1+kg*3
      nbaks=4+kg*10
      do 55 i=1,nbaks
   55 ass(i)=0.d0
      if(n.eq.nsl) goto 65
      if(kg.ne.0) goto 75
c*** Minor at the free (ocean) surface *** 
      ass(1)=dsign(1.d0,a(1,n))
      goto 80
   65 if(kg.eq.0) goto 75
c*** Minor at the free (solid) surface *** 
      asi1=a(3,n)*a(3,n)+a(12,n)*a(12,n)
      asi2=a(4,n)*a(4,n)+a(11,n)*a(11,n)
      if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
      if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
      goto 85
   75 asi1=a(3,n)*a(3,n)
      asi2=a(4,n)*a(4,n)
      if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
      if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
      if(n.eq.nsl) goto 85
c*** Propagate through ocean (downward) ***
   80 call fprpmn(n,nslp1,ass,vf,nbakf,jexp)
      call fsbm(ass,kg,iback)
   85 nto=max0(ls,nocp1)
c*** Propagate through mantle (downward) ***
      call sprpmn(nsl,nto,ass,vf,nbaks,jexp)
      if(nto.eq.ls) goto 90
c*** Obtain the solution on the fluid side of the s/f CMB.
      call sfbm(ass,kg,iback)
      nto=max0(ls,nicp1)
c*** Propagate through outer core (downward) ***
      call fprpmn(noc,nto,ass,vf,nbakf,jexp)
      if(nto.eq.ls) goto 90
c*** Obtain the solution on the solid side of the f/s ICB.
      call fsbm(ass,kg,iback)
      nto=max0(ls,2)
c*** Propagate through inner core (downward) ***
      call sprpmn(nic,nto,ass,vf,nbaks,jexp)
c   90 if(dabs(det).gt.1.d-4) call remedy(ls)
   90 continue
      call eifout(ls,ratiomax)
      goto 100
c*** radial modes ***
    2 ls=2
      call rps(ls,ass)
      call rprop(ls,n,ass)
      det=a(2,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
      knt=kount-1
      if(ifeif.eq.0) goto 100
      a(1,1)=0.d0
      a(2,1)=0.d0
      do 205 i=ls,n
      ff=fcon(i)*(1.d0+xlam(i)*fct)
      cc=ccon(i)*(1.d0+xa2(i)*fct)
  205 a(2,i)=(a(2,i)-2.d0*ff*a(1,i)/r(i))/cc
      zi(1)=0.d0
      zi(2)=0.d0
      zi(3)=0.d0
      do 210 i=ls,n
      im=i-1
  210 if(r(i).ne.r(im)) call gauslv(r(im),r(i),im,zi,3)
      rnrm=1.d0/(w*dsqrt(zi(1)))
      cg=0.d0
      qinv=zi(2)/(wsq*zi(1))
      wray=dsqrt(zi(3)/zi(1))
      do 215 i=2,n
      do 215 j=1,2
  215 a(j,i)=a(j,i)*rnrm
      goto 100
c*** toroidal modes ***
    3 nb=nocp1
      n2=nsl
      ass(1)=1.d0
      ass(2)=0.d0
      if(jcom.eq.2) goto 300
      nb=2
      a(1,1)=0.d0
      a(2,1)=0.d0
      n2=nic
  300 q=0.d0
      ls=nb
      IF(r(ls)==0.d0) ls=ls+1
      call startl(ls,n2,fmu,ls,q)
      if(ls.ne.nocp1) call tps(ls,ass)
      call tprop(ls,n2,ass)
      det=a(2,n2)/dsqrt(a(1,n2)*a(1,n2)+a(2,n2)*a(2,n2))
      if(ifeif.eq.0) goto 335
      do 305 i=ls,n2
  305 a(2,i)=a(1,i)/r(i)+a(2,i)/(lcon(i)*(1.d0+qshear(i)*fct))
      if(ls.eq.nb) goto 315
      ls1=ls-1
      do 310 i=nb,ls1
      a(1,i)=0.d0
  310 a(2,i)=0.d0
  315 do 320 i=1,4
  320 zi(i)=0.d0
      do 325 i=ls,n2
      im=i-1
  325 if(r(i).ne.r(im)) call gauslv(r(im),r(i),im,zi,4)
      rnrm=1.d0/(w*dsqrt(zi(1)))
c      cg=(fl+0.5d0)*zi(2)/(w*zi(1))
      cg=sfl3*zi(2)/(w*zi(1))
      qinv=zi(3)/(wsq*zi(1))
      wray=dsqrt(zi(4)/zi(1))
      do 330 i=ls,n2
      do 330 j=1,2
  330 a(j,i)=a(j,i)*rnrm
      goto 100
  335 if(knsw.ne.1) goto 100
      knt=kount-1
      if(jcom.eq.4.or.l.eq.1) goto 100
      irem=mod(knt,2)
      if(irem.eq.0.and.det.lt.0.d0) goto 100
      if(irem.ne.0.and.det.gt.0.d0) goto 100
      knt=knt+1

100	continue
      
      return
      end
c
c
      subroutine sprpmn(jf,jl,f,h,nvesm,iexp)
      
c*** propagate a minor vector in a solid region from level jf to jl ***

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,ndata),inorm(ndata)
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,stepf,maxo,in
      dimension f(1),h(nvesm,1),s(14),fp(14),rne(6)
      data econst/1048576.d0/
      
      maxo1=maxo-1
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 45
   10 x=y
      y=r(i)
      if(y.eq.x) goto 45
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=0.5d0*(alfsq+betasq+delsq)-fl3/(x*x)
      qt=dsqrt(dabs(fksq))+dsqrt(dabs(fksq-delsq))+2.d0/zs
      q=(qt+dble(kg)*sfl3/x)/stepf
      del=dble(jud)*step(maxo)/q
      dxs=0.d0
   15 y=x+del
      if(dble(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,maxo1)
      dxs=dx
      do 30 j=1,nvesm
   30 s(j)=f(j)
      do 35 ni=1,in
      z=x+c(ni)
      call derms(iq,z,f,h(1,ni),0,qff,qll,qaa)
   35 call rkdot(f,s,h,nvesm,ni)
      if(knsw.ne.1) goto 40
      call derms(iq,y,f,fp,1,qff,qll,qaa)
      call zknt(s,h,f,fp,x,y,1)
   40 x=y
      if(y.ne.r(i)) goto 15
   45 size=dabs(f(1))
      do 50 j=2,nvesm
   50 size=dmax1(size,dabs(f(j)))
   55 if(size.lt.1024.d0) goto 65
      do 60 j=1,nvesm
   60 f(j)=f(j)/econst
      size=size/econst
      iexp=iexp+20
      goto 55
   65 if(iback.eq.0) goto 85
      inorm(i)=inorm(i)+iexp
      if(kg.eq.0) goto 70
      t1=f(4)+f(8)
      t2=t1+f(4)
      t1=t1+f(8)
      t3=f(8)-f(4)
      rne(1)=ar(6,i)*f(10)-ar(14,i)*f(9)+ar(13,i)*t3
     1      -ar(1,i)*f(7)-ar(7,i)*f(6)+ar(8,i)*f(5)
     2      +ar(12,i)*f(3)-ar(2,i)*f(2)+ar(3,i)*f(1)
      rne(2)=ar(6,i)*f(13)+ar(14,i)*t2+ar(13,i)*f(12)
     1      -ar(1,i)*f(11)-ar(9,i)*f(6)-ar(7,i)*f(5)
     2      +ar(11,i)*f(3)-ar(4,i)*f(2)-ar(2,i)*f(1)
      rne(3)=ar(6,i)*f(14)-ar(7,i)*t1-ar(8,i)*f(12)
     1      +ar(13,i)*f(11)-ar(9,i)*f(9)+ar(14,i)*f(7)
     2      +ar(10,i)*f(3)+ar(11,i)*f(2)+ar(12,i)*f(1)
      rne(4)=ar(14,i)*f(14)+ar(7,i)*f(13)+ar(12,i)*f(12)
     1      -ar(2,i)*f(11)-ar(9,i)*f(10)-ar(11,i)*t3
     2      +ar(4,i)*f(7)+ar(10,i)*f(5)+ar(5,i)*f(1)
      rne(5)=ar(13,i)*f(14)+ar(8,i)*f(13)-ar(12,i)*t2
     1      -ar(3,i)*f(11)+ar(7,i)*f(10)-ar(11,i)*f(9)
     2      -ar(2,i)*f(7)+ar(10,i)*f(6)+ar(5,i)*f(2)
      rne(6)=ar(1,i)*f(14)+ar(13,i)*f(13)-ar(2,i)*t1
     1      -ar(3,i)*f(12)+ar(14,i)*f(10)-ar(4,i)*f(9)
     2      -ar(11,i)*f(6)-ar(12,i)*f(5)+ar(5,i)*f(3)
      goto 75
   70 rne(1)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
      rne(2)=-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
      rne(3)=-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
      rne(4)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
   75 do 80 jj=1,6
   80 ar(jj,i)=rne(jj)
      goto 95
   85 inorm(i)=iexp
      do 90 j=1,nvesm
   90 ar(j,i)=f(j)
   95 if(i.eq.jl) goto 100
      i=i+jud
      goto 10
      
100	continue
	
      return
      end
c
c
      subroutine fprpmn(jf,jl,f,h,nvefm,iexp)
      
c*** Propagate the minor vector in a fluid region from level jf to jl ***

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,ndata),inorm(ndata)
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,stepf,maxo,in
      dimension f(1),h(nvefm,1),s(5),fp(5)
      data econst/1048576.d0/

      if(nvefm.eq.1) goto 85
      maxo1=maxo-1
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      goto 45
   10 x=y
      y=r(i)
      if(y.eq.x) goto 45
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=(dsqrt(dabs(alfsq-fl3/(x*x)))+1.d0/zs+dble(kg)*sfl3/x)/stepf
      del=dble(jud)*step(maxo)/q
      dxs=0.d0
   15 y=x+del
      if(dble(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,maxo1)
      dxs=dx
      do 30 j=1,nvefm
   30 s(j)=f(j)
      do 35 ni=1,in
      z=x+c(ni)
      call dermf(iq,z,f,h(1,ni),0,qff)
   35 call rkdot(f,s,h,nvefm,ni)
      if(knsw.ne.1) goto 40
      call dermf(iq,y,f,fp,1,qff)
      call zknt(s,h,f,fp,x,y,0)
   40 x=y
      if(y.ne.r(i)) go to 15
   45 size=dabs(f(1))
      do 50 j=2,nvefm
   50 size=dmax1(size,dabs(f(j)))
   55 if(size.lt.1024.d0) goto 65
      do 60 j=1,nvefm
   60 f(j)=f(j)/econst
      size=size/econst
      iexp=iexp+20
      goto 55
   65 if(iback.eq.0) goto 70
      inorm(i)=inorm(i)+iexp
      rne2   =-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
      ar(1,i)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
      rne3   =-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
      ar(4,i)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
      ar(2,i)=rne2
      ar(3,i)=rne3
      goto 80
   70 inorm(i)=iexp
      do 75 j=1,nvefm
   75 ar(j,i)=f(j)
   80 if(i.eq.jl) goto 100
      i=i+jud
      goto 10
   85 do 90 i=jl,jf
      inorm(i)=inorm(i)+iexp
      do 91 j=1,2
   91 ar(j,i)=ar(j,i)*f(1)
   90 continue

100	continue
   
      return
      end
c
c
      subroutine derms(iq,z,f,fp,iknt,qff,qll,qaa)
      
c*** calculates minor vector derivative (fp) in a solid ***

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 nn,ll,lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension f(1),fp(1)
      
      if(iknt.ne.0) goto 19
      t=z-r(iq)
      if(t.ne.0.d0) goto 5
      ro=rho(iq)
      gr=g(iq)
      ff=fcon(iq)*qff
      ll=lcon(iq)*qll
      nn=ncon(iq)*qll
      cc=ccon(iq)*qaa
      aa=acon(iq)*qaa
      goto 15
    5 ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      if(ifanis.ne.0) goto 10
      nn=ll
      cc=ff+ll+ll
      aa=cc
      goto 15
   10 nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
      aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
   15 zr=1.d0/z
      sfl3z=sfl3*zr
      rogr=ro*gr
      c11=1.d0/cc
      c22=1.d0/ll
      dmg=aa-nn-ff*ff*c11
      zdmg=zr*dmg
      t11=-2.d0*ff*zr*c11+zr
      t12=sfl3z*ff*c11
      t21=-sfl3z
      t22=zr+zr
      s22=-ro*wsq
      s11=s22+4.d0*zr*(zdmg-rogr)
      s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn)
      s12=sfl3z*(rogr-zdmg-zdmg)
      if(kg.ne.0) goto 25
      s11=s11+4.d0*ro*ro
      if(iback.eq.1) goto 20
      b11=t11+t22
      b33=t11-t22
      fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
      fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
   19 if(kg.ne.0) goto 29
      fp(3)=s22*f(1)-2.d0*t12*f(2)+b33*f(3)+c11*f(5)
      fp(4)=-s11*f(1)+2.d0*t21*f(2)-b33*f(4)-c22*f(5)
      fp(5)=-2.d0*s12*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
      goto 100
   20 fp(1)=t22*f(1)-t21*f(2)-c22*f(3)
      fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
      fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
      fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
      goto 100
   25 t31=-4.d0*ro
      t33=-fl*zr
      s13=-fl1*zr*ro
      s23=ro*sfl3z
      if(iback.eq.1) goto 30
      b11=t11+t22-t33
      b33=t11-t22-t33
      b44=t22-t11-t33
      b55=-t11-t22-t33
      b32=-t12-t12
      b42=t21+t21
      b52=-s12-s12
      b313=-s23-s23
      b414=s13+s13
      b914=t31+t31
      fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
      fp(2)=s12*f(1)-t33*f(2)-t21*f(3)+t12*f(4)-s13*f(13)-s23*f(14)
      fp(6)=4.d0*f(1)-b55*f(6)+c22*f(8)-c11*f(9)
      fp(7)=4.d0*f(2)+s12*f(6)+t33*f(7)-t21*f(8)+t12*f(9)-t31*f(13)
      fp(8)=4.d0*f(3)+s22*f(6)+b32*f(7)-b44*f(8)+c11*f(10)
      fp(9)=4.d0*f(4)-s11*f(6)+b42*f(7)-b33*f(9)-c22*f(10)+b914*f(14)
      fp(10)=4.d0*f(5)+b52*f(7)+s11*f(8)-s22*f(9)-b11*f(10)+b914*f(12)
      fp(11)=-t31*f(2)+s13*f(7)+s23*f(9)-t11*f(11)+t21*f(12)
     +      -s11*f(13)+s12*f(14)
      fp(12)=t31*f(3)+s23*f(7)-s13*f(8)+t12*f(11)-t22*f(12)
     +      +s12*f(13)-s22*f(14)
      fp(13)=s23*f(6)-c11*f(11)+t11*f(13)-t12*f(14)
      fp(14)=-t31*f(1)+s13*f(6)-c22*f(12)-t21*f(13)+t22*f(14)
   29 fp(3)=s22*f(1)+b32*f(2)+b33*f(3)+c11*f(5)+b313*f(13)
      fp(4)=-s11*f(1)+b42*f(2)+b44*f(4)-c22*f(5)+b414*f(14)
      fp(5)=b52*f(2)+s11*f(3)-s22*f(4)+b55*f(5)-b313*f(11)+b414*f(12)
      goto 100
   30 b11=t22+t33
      b22=t11+t33
      b33=t11+t22
      b55=t22-t33
      b66=t11-t33
      b99=t11-t22
      t4=f(4)+f(8)
      t5=t4+f(8)
      t4=t4+f(4)
      fp(1)=b11*f(1)-t21*f(2)-t31*f(3)-4.d0*f(5)+c22*f(7)
      fp(2)=-t12*f(1)+b22*f(2)-4.d0*f(6)+c11*f(11)
      fp(3)=b33*f(3)-c22*f(9)+c11*f(12)
      fp(4)=-s23*f(1)+s13*f(2)+t31*f(6)
      fp(5)=s13*f(3)+b55*f(5)-t21*f(6)-c22*f(10)
      fp(6)=s23*f(3)-t12*f(5)+b66*f(6)-c11*f(13)
      fp(7)=s22*f(1)-s12*f(2)-b55*f(7)+t31*f(9)+4.d0*f(10)+t12*f(11)
      fp(8)=s23*f(1)-s12*f(3)-t21*f(9)+t12*f(12)
      fp(9)=s23*f(2)-s22*f(3)-t12*t5+b99*f(9)-c11*f(14)
      fp(10)=s23*(f(4)-f(8))-s22*f(5)+s12*f(6)+s13*f(9)-b11*f(10)
     1      +t12*f(13)
      fp(11)=-s12*f(1)+s11*f(2)-t4*t31+t21*f(7)-b66*f(11)+4.d0*f(13)
      fp(12)=-s13*f(1)+s11*f(3)+t21*t5-t31*f(5)-b99*f(12)+c22*f(14)
      fp(13)=-t4*s13+s12*f(5)-s11*f(6)+t21*f(10)-s23*f(12)-b22*f(13)
      fp(14)=s12*t5-s13*f(7)-s11*f(9)+t31*f(10)-s23*f(11)+s22*f(12)
     1      -b33*f(14)
     
100	continue

      return
      end
c
c
      subroutine dermf(iq,z,f,fp,iknt,qff)
      
c*** Calculates minor vector derivative (fp) in a fluid ***

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension f(1),fp(1)
      
      if(iknt.ne.0) goto 14
      t=z-r(iq)
      if(t.ne.0.d0) goto 5
      ro=rho(iq)
      flu=fcon(iq)*qff
      gr=g(iq)
      goto 10
    5 ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      flu=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
   10 t21=-4.d0*ro
      zr=1.d0/z
      t12=fl3*zr*zr/wsq
      t11=gr*t12-zr
      s11=ro*(gr*gr*t12-wsq)+t21*gr*zr
      c11=-t12/ro+1.d0/flu
   14 if(kg.ne.0) goto 15
      fp(1)=t11*f(1)+c11*f(2)
      fp(2)=(s11-t21*ro)*f(1)-t11*f(2)
      goto 100
   15 if(iknt.ne.0) goto 19
      t22=-fl*zr
      s22=ro*t12
      b11=t11+t22
      s12=ro*b11
      if(iback.eq.1) goto 20
      b33=t11-t22
      fp(1)=b11*f(1)+4.d0*f(3)-c11*f(4)
      fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
   19 fp(3)=s22*f(1)-(t12+t12)*f(2)+b33*f(3)+c11*f(5)
      fp(4)=-s11*f(1)+(t21+t21)*f(2)-b33*f(4)-4.d0*f(5)
      fp(5)=-(s12+s12)*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
      goto 100
   20 fp(1)=t22*f(1)-t21*f(2)-4.d0*f(3)
      fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
      fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
      fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
      
100	continue

      return
      end
c
c
      subroutine eifout(lsmin,rmax)
      
c*** Massages spheroidal mode eigenfunctions before output ***

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 ll,lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,ndata),inorm(ndata)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension zi(4)
      
      i1=min0(nic,max0(2,lsmin))
      i2=nic
    5 if(i1.eq.i2) goto 20
      do 10 iq=i1,i2
      ff=fcon(iq)*(1.d0+xlam(iq)*fct)
      ll=lcon(iq)*(1.d0+qshear(iq)*fct)
      zr=1.d0/r(iq)
      sfl3z=sfl3*zr
      d=1.d0/(ccon(iq)*(1.d0+xa2(iq)*fct))
      v=a(2,iq)
      if(kg.ne.0) goto 15
      a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(3,iq)
      a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(4,iq)/ll
      a(5,iq)=0.d0
      a(6,iq)=0.d0
      goto 10
   15 a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(4,iq)
      a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(5,iq)/ll
      a(5,iq)=a(3,iq)
      a(6,iq)=4.d0*(a(6,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
   10 a(3,iq)=v
   20 if(i2.eq.nsl) goto 25
      i1=min0(nsl,max0(lsmin,nocp1))
      i2=nsl
      goto 5
   25 i1=min0(noc,max0(lsmin,nicp1))
      i2=noc
   30 if(i1.eq.i2) goto 50
      do 35 iq=i1,i2
      zr=1.d0/r(iq)
      sfl3z=sfl3*zr
      ffi=1.d0/(flam(iq)*(1.d0+xlam(iq)*fct))
      if(kg.ne.0) goto 40
      p=a(2,iq)
      a(5,iq)=0.d0
      a(6,iq)=0.d0
      goto 45
   40 p=a(3,iq)
      a(5,iq)=a(2,iq)
      a(6,iq)=4.d0*(a(4,iq)-rho(iq)*a(1,iq))-fl*zr*a(5,iq)
   45 a(3,iq)=sfl3z*(g(iq)*a(1,iq)-p/rho(iq)+a(5,iq))/wsq
      a(2,iq)=sfl3z*a(3,iq)-a(1,iq)*zr+p*ffi
   35 a(4,iq)=sfl3z*(a(1,iq)+p*(qro(1,iq)/(rho(iq)**2)+g(iq)*ffi)/wsq)
   50 if(n.eq.nsl.or.i2.eq.n) goto 55
      i1=nslp1
      i2=n
      goto 30
   55 imax=0
      do 60 iq=lsmin,n
   60 imax=max0(inorm(iq),imax)
      do 65 iq=lsmin,n
      iexp=inorm(iq)-imax
      al=0.d0
      if(iexp.ge.-120) al=2.d0**iexp
      do 65 j=1,6
   65 a(j,iq)=a(j,iq)*al
      lsm1=max0(1,lsmin-1)
      do 70 i=1,lsm1
      do 70 j=1,6
   70 a(j,i)=0.d0
      if(l.gt.1.or.lsmin.gt.2) goto 75
      a(2,1)=1.5d0*a(1,2)/r(2)-0.5d0*a(2,2)
      a(4,1)=1.5d0*a(3,2)/r(2)-0.5d0*a(4,2)
   75 do 80 j=1,4
   80 zi(j)=0.d0
      i1=max0(lsmin,2)
      do 85 iq=i1,n
      ip=iq-1
   85 if(r(iq).ne.r(ip)) call gauslv(r(ip),r(iq),ip,zi,4)
      cg=zi(2)/(w*zi(1))
      wray=dsqrt(2.d0*zi(4)/zi(1))
      qinv=2.d0*zi(3)/(wsq*zi(1))
      rnorm=1.d0/(w*dsqrt(zi(1)))
      do 90 iq=i1,n
      zr=1.d0/r(iq)
      a(1,iq)=a(1,iq)*zr
      a(2,iq)=(a(2,iq)-a(1,iq))*zr
      a(3,iq)=a(3,iq)*zr
      a(4,iq)=(a(4,iq)-a(3,iq))*zr
      a(5,iq)=a(5,iq)*zr
      a(6,iq)=(a(6,iq)-a(5,iq))*zr
      a(1,iq)=a(1,iq)*rnorm
      a(2,iq)=a(2,iq)*rnorm
      a(3,iq)=a(3,iq)*rnorm
      a(4,iq)=a(4,iq)*rnorm
      a(5,iq)=a(5,iq)*rnorm
   90 a(6,iq)=a(6,iq)*rnorm
      if(lsmin.gt.2.or.l.gt.2) goto 100
      if(l.eq.2) goto 95
      a(1,1)=a(1,2)-0.5d0*a(2,2)*r(2)
      a(2,1)=0.d0
      a(3,1)=a(3,2)-.5d0*a(4,2)*r(2)
      a(4,1)=0.d0
      a(6,1)=1.5d0*a(5,2)/r(2)-0.5d0*a(6,2)
      goto 100
   95 a(2,1)=1.5d0*a(1,2)/r(2)-0.5d0*a(2,2)
      a(4,1)=1.5d0*a(3,2)/r(2)-0.5d0*a(4,2)
      
100	continue

	aintic=0.d0
	bintic=0.d0
	cintic=0.d0
	dintic=0.d0
	do i=1,nic
	  aintic=aintic+dabs(a(1,i))/dble(nic)
	  bintic=bintic+dabs(a(2,i))/dble(nic)
	  cintic=cintic+dabs(a(3,i))/dble(nic)
	  dintic=dintic+dabs(a(4,i))/dble(nic)
	end do
	aintoc=0.d0
	bintoc=0.d0
	cintoc=0.d0
	dintoc=0.d0
	do i=nicp1,noc
	  aintoc=aintoc+dabs(a(1,i))/dble(noc-nic)
	  bintoc=bintoc+dabs(a(2,i))/dble(noc-nic)
	  cintoc=cintoc+dabs(a(3,i))/dble(noc-nic)
	  dintoc=dintoc+dabs(a(4,i))/dble(noc-nic)
	end do
	aintmt=0.d0
	bintmt=0.d0
	cintmt=0.d0
	dintmt=0.d0
	do i=nocp1,n
	  aintmt=aintmt+dabs(a(1,i))/dble(n-noc)
	  bintmt=bintmt+dabs(a(2,i))/dble(n-noc)
	  cintmt=cintmt+dabs(a(3,i))/dble(n-noc)
	  dintmt=dintmt+dabs(a(4,i))/dble(n-noc)
	end do
	
	aratio=1.d98
	if(aintmt.gt.0.d0) aratio=(aintic+aintoc)/aintmt
	bratio=1.d98
	if(bintmt.gt.0.d0) bratio=(bintic+bintoc)/bintmt
	cratio=1.d98
	if(cintmt.gt.0.d0) cratio=(cintic+cintoc)/cintmt
	dratio=1.d98
	if(dintmt.gt.0.d0) dratio=(dintic+dintoc)/dintmt
	abmax=dmax1(aratio,bratio)
	abcmax=dmax1(abmax,cratio)
	rmax=dmax1(abcmax,dratio)
	
c	write(*,'(4(1x,e15.8))') aratio,bratio,cratio,dratio

      return
      end
c
c
      subroutine gauslv(r1,r2,iq,fint,nint)
      
c*** Fifth order Gauss-Legendre integration ***

      implicit real*8(a-h,o-z)
      dimension fint(1),vals(4),vals1(4),sum(4),w(2),x(2)
      data w,x/0.478628670499366d0,0.236926885056189d0,
     +         0.538469310105683d0,0.906179845938664d0/
     
      y1=0.5d0*(r2+r1)
      y2=0.5d0*(r2-r1)
      call intgds(y1,iq,vals)
      do 5 j=1,nint
    5 sum(j)=0.568888888888889d0*vals(j)
      do 10 i=1,2
      t1=x(i)*y2
      call intgds(y1+t1,iq,vals)
      call intgds(y1-t1,iq,vals1)
      do 10 j=1,nint
   10 sum(j)=sum(j)+w(i)*(vals(j)+vals1(j))
      do 15 j=1,nint
   15 fint(j)=fint(j)+y2*sum(j)
   
      return
      end
c
c
      subroutine sdepth(wdim,ls)
      
c*** finds starting level,ls, for a given l and w ***

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      data aw,bw,dw/-2.d-3,2.25d-3,1.28d-3/
      
      q=0.d0
      w=wdim/wn
      wsoc=aw+dw*fl
c     if(wdim.gt.wsoc) goto 10
      call startl(nocp1,nsl,fmu,ls,q)
      if(ls.eq.nsl) ls=ls-1
      if(ls.gt.nocp1) goto 100
   10 wsic=aw+bw*fl
      if(wdim.gt.wsic) goto 20
      call startl(nicp1,noc,flam,ls,q)
      if(ls.eq.noc) ls=ls-1
      if(ls.gt.nicp1) goto 100
   20 call startl(2,nic,fmu,ls,q)
      if(ls.eq.nic) ls=ls-1
      
100	continue

      return
      end
c
c
      subroutine sfbm(ass,kg,iback)
      
c*** Convert minor vector at a solid/fluid boundary ***

      implicit real*8(a-h,o-z)
      dimension ass(14),as(14)
      
      do 10 j=1,14
      as(j)=ass(j)
   10 ass(j)=0.d0
      if(iback.eq.1) goto 30
      if(kg.ne.0) goto 20
      ass(1)=as(3)
      ass(2)=as(5)
      goto 100
   20 ass(1)=as(8)
      ass(2)=-as(12)
      ass(3)=as(3)
      ass(4)=-as(10)
      ass(5)=as(5)
      goto 100
   30 if(kg.ne.0) goto 40
      ass(1)=-as(3)
      goto 100
   40 ass(1)=as(7)
      ass(2)=-as(9)
      ass(3)=-as(10)
      ass(4)=-as(14)

100	continue
      
      return
      end
c
c
      subroutine fsbm(ass,kg,iback)
      
c*** Convert minor vector at a fluid/solid boundary ***

      implicit real*8(a-h,o-z)
      dimension ass(14),as(14)
      
      do 10 j=1,14
      as(j)=ass(j)
   10 ass(j)=0.d0
      if(iback.eq.1) goto 30
      if(kg.ne.0) goto 20
      ass(1)=as(1)
      ass(4)=-as(2)
      goto 100
   20 ass(6)=as(1)
      ass(14)=as(2)
      ass(1)=as(3)
      ass(9)=as(4)
      ass(4)=-as(5)
      goto 100
   30 if(kg.ne.0) goto 40
      ass(1)=-as(1)
      goto 100
   40 ass(1)=-as(1)
      ass(3)=-as(2)
      ass(5)=-as(3)
      ass(12)=as(4)
      
100	continue

      return
      end
c
c
      subroutine zknt(s,sp,f,fp,x,y,ifsol)
      
c*** Given minor vector and derivs, constructs mode count ***

      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension s(1),sp(1),f(1),fp(1),xs(4),val(4)

c   ifsol=1 for solid. kg=0 no gravity effect.

	do i=1      
      if(ifsol.eq.0.and.kg.eq.0) goto 5
      y1=s(5)
      y2=f(5)
      y1p=sp(5)
      y2p=fp(5)
      t1=s(3)-s(4)
      t2=f(3)-f(4)
      t1p=sp(3)-sp(4)
      t2p=fp(3)-fp(4)
      goto 10
    5 y1=s(2)
      y2=f(2)
      y1p=sp(2)
      y2p=fp(2)
      t1=s(1)
      t2=f(1)
      t1p=sp(1)
      t2p=fp(1)
   10 h=y-x
      ns=0
      if(kount.ne.0) goto 15
      a1=y2-y1
      a2=0.d0
      a3=0.d0
      a22=0.d0
      a33=0.d0
      goto 50
   15 a1=h*y1p
      a2=-h*(2.d0*y1p+y2p)+3.d0*(y2-y1)
      a3=h*(y1p+y2p)-2.d0*(y2-y1)
      a33=3.d0*a3
      a22=2.d0*a2
      if(a3.ne.0.d0) goto 20
      if(a2.eq.0.d0) goto 50
      xs(2)=-a1/a22
      if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
      goto 50
   20 disc=a2*a2-a1*a33
      if(disc) 50,25,30
   25 xs(2)=-a2/a33
      if(xs(2).ge.0.d0.and.xs(2).le.1.d0) ns=1
      goto 50
   30 disc=dsqrt(disc)
      tr1=(-a2+disc)/a33
      tr2=(-a2-disc)/a33
      if(dabs(a33).gt.dabs(a1)) goto 35
      fac=a1/a33
      tr1=fac/tr1
      tr2=fac/tr2
   35 if(tr1.lt.0.d0.or.tr1.gt.1.d0) goto 40
      xs(2)=tr1
      ns=1
   40 if(tr2.lt.0.d0.or.tr2.gt.1.d0) goto 50
      ns=ns+1
      xs(ns+1)=tr2
      if(ns.lt.2) goto 50
      if(tr2.ge.tr1) goto 50
      xs(2)=tr2
      xs(3)=tr1
   50 val(1)=y1
      xs(1)=0.d0
      ns2=ns+2
      val(ns2)=y2
      xs(ns2)=1.d0
      if(ns.eq.0) goto 60
      ns1=ns+1
      do 55 j=2,ns1
      t=xs(j)
   55 val(j)=y1+t*(a1+t*(a2+t*a3))
   60 ift=0
      do 100 j=2,ns2
      if(val(j-1)*val(j).gt.0.d0) goto 100
      if(val(j-1).ne.0.d0) goto 65
      tes=t1*a1
      goto 90
   65 rt1=0.5d0*(xs(j-1)+xs(j))
      rt=rt1
      do 70 i=1,5
      v=y1+rt*(a1+rt*(a2+rt*a3))
      vp=a1+rt*(a22+rt*a33)
      add=-v/vp
      rt=rt+add
      if(dabs(add).lt.1.d-5) goto 75
      if(dabs(rt-rt1).le.0.5d0) goto 70
      rt=rt1
      goto 75
   70 continue
   75 if(ift.ne.0) goto 85
      if(kount.ne.0) goto 80
      b1=t2-t1
      b2=0.d0
      b3=0.d0
      goto 85
   80 b1=h*t1p
      b2=-h*(2.d0*t1p+t2p)+3.d0*(t2-t1)
      b3=h*(t1p+t2p)-2.d0*(t2-t1)
      ift=1
   85 tes=t1+rt*(b1+rt*(b2+rt*b3))
      vp=a1+rt*(a22+rt*a33)
      tes=tes*vp
   90 if(tes.lt.0.d0) kount=1+kount
      if(tes.gt.0.d0) kount=kount-1
  100 continue
  
      return
      end
c
c
      subroutine baylis(q,maxo1)
      
c    baylis returns the coefficients for rks integration.
c    See E. Baylis Shanks(1966 A. M. S.) and references therein for the
c    coefficients. The eight Runge-Kutta-Shanks formulae are (1-1), (2-2)
c    (3-3), (4-4), (5-5), (6-6), (7-7), and (8-10). For orders greater than 4 the
c    formulae are approximate rather than exact so incurring less roundoff.

      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,stepf,maxo,i
      
      ds=q*dabs(dx)
      do 10 j=1,maxo1
      if(ds.gt.step(j)) go to 10
      i=j
      go to 15
   10 continue
      i=maxo
   15 c(1)=0.d0
      goto (1,2,3,4,5,6,7,8),i
    1 b(1)=dx
      goto 100
    2 c(2)=dx
      b(1)=dx
      b(2)=0.5d0*dx
      b(3)=1.d0
      goto 100
    3 c(2)=0.5d0*dx
      c(3)=dx
      b(1)=c(2)
      b(2)=-dx
      b(3)=-2.d0
      b(4)=0.16666666666667d0*dx
      b(5)=4.d0
      b(6)=1.d0
      goto 100
    4 c(2)=0.01d0*dx
      c(3)=0.6d0*dx
      c(4)=dx
      b(1)=c(2)
      b(2)=-0.17461224489790d+02*dx
      b(3)=-0.10343618513324d+01
      b(4)=0.59691275167780d+02*dx
      b(5)=-0.10140620414448d+01
      b(6)=0.30814908546230d-01
      b(7)=-0.25555555555556d+01*dx
      b(8)=-0.11165449632656d+01
      b(9)=-0.22568165070006d+00
      b(10)=-0.49077733860351d-01
      goto 100
    5 c(2)=1.1111111111111d-04*dx
      c(3)=3.0d-01*dx
      c(4)=7.5d-01*dx
      c(5)=dx
      b(1)=c(2)
      b(2)=-0.40470000000000d+03*dx
      b(3)=-0.10007412898443d+01
      b(4)=0.25301250000000d+04*dx
      b(5)=-0.10004446420631d+01
      b(6)=0.74107010523195d-03
      b(7)=-0.11494333333333d+05*dx
      b(8)=-0.10004929965491d+01
      b(9)=0.52629261224803d-03
      b(10)=-0.12029545422812d-03
      b(11)=0.92592592592593d-01*dx
      b(12)=0.00000000000000d+00
      b(13)=0.47619047619048d+01
      b(14)=0.42666666666667d+01
      b(15)=0.77142857142857d+00
      goto 100
    6 c(2)=3.3333333333333d-03*dx
      c(3)=0.2d0*dx
      c(4)=0.6d0*dx
      c(5)=9.3333333333333d-01*dx
      c(6)=dx
      b(1)=c(2)
      b(2)=-0.58000000000000d+01*dx
      b(3)=-0.10344827586207d+01
      b(4)=0.64600000000000d+02*dx
      b(5)=-0.10216718266254d+01
      b(6)=0.30959752321982d-01
      b(7)=-0.62975802469136d+03*dx
      b(8)=-0.10226149961576d+01
      b(9)=0.24906685695466d-01
      b(10)=-0.37737402568887d-02
      b(11)=-0.54275714285714d+04*dx
      b(12)=-0.10225567867765d+01
      b(13)=0.25375487829097d-01
      b(14)=-0.31321559234596d-02
      b(15)=0.12921040478749d-03
      b(16)=0.53571428571429d-01*dx
      b(17)=0.00000000000000d+00
      b(18)=0.61868686868687d+01
      b(19)=0.77777777777778d+01
      b(20)=0.40909090909091d+01
      b(21)=-0.38888888888889d+00
      goto 100
    7 c(2)=5.2083333333333d-03*dx
      c(3)=1.6666666666667d-01*dx
      c(4)=0.5d0*dx
      c(5)=dx
      c(6)=8.3333333333333d-01*dx
      c(7)=dx
      b(1)=c(2)
      b(2)=-0.25000000000000d+01*dx
      b(3)=-0.10666666666667d+01
      b(4)=0.26166666666667d+02*dx
      b(5)=-0.10421204027121d+01
      b(6)=0.61228682966918d-01
      b(7)=-0.64500000000000d+03*dx
      b(8)=-0.10450612653163d+01
      b(9)=0.51262815703925d-01
      b(10)=-0.77519379844961d-02
      b(11)=-0.93549382716049d+02*dx
      b(12)=-0.10450293206756d+01
      b(13)=0.48394546673620d-01
      b(14)=-0.11877268228307d-01
      b(15)=-0.39590894094358d-03
      b(16)=0.35111904761905d+03*dx
      b(17)=-0.10446476812124d+01
      b(18)=0.52479782656724d-01
      b(19)=-0.71200922221468d-02
      b(20)=-0.61029361904114d-03
      b(21)=0.27463212856852d-02
      b(22)=0.46666666666667d-01*dx
      b(23)=0.57857142857143d+01
      b(24)=0.78571428571429d+01
      b(25)=0.00000000000000d+00
      b(26)=b(23)
      b(27)=0.10000000000000d+01
      goto 100
    8 c(2)=0.14814814814815d0*dx
      c(3)=0.22222222222222d0*dx
      c(4)=0.33333333333333d0*dx
      c(5)=0.5d0*dx
      c(6)=0.66666666666667d0*dx
      c(7)=0.16666666666667d0*dx
      c(8)=dx
      c(9)=0.83333333333333d0*dx
      c(10)=dx
      b(1)=c(2)
      b(2)=0.55555555555556d-01*dx
      b(3)=0.30000000000000d+01
      b(4)=0.83333333333333d-01*dx
      b(5)=0.00000000000000d+00
      b(6)=0.30000000000000d+01
      b(7)=0.12500000000000d+00*dx
      b(8)=0.00000000000000d+00
      b(9)=0.00000000000000d+00
      b(10)=0.30000000000000d+01
      b(11)=0.24074074074074d+00*dx
      b(12)=0.00000000000000d+00
      b(13)=-0.20769230769231d+01
      b(14)=0.32307692307692d+01
      b(15)=0.61538461538461d+00
      b(16)=0.90046296296295d-01*dx
      b(17)=0.00000000000000d+00
      b(18)=-0.13881748071980d+00
      b(19)=0.24832904884319d+01
      b(20)=-0.21182519280206d+01
      b(21)=0.62467866323908d+00
      b(22)=-0.11550000000000d+02*dx
      b(23)=-0.35064935064935d+00
      b(24)=0.50389610389610d+01
      b(25)=-0.28398268398268d+01
      b(26)=0.52813852813853d+00
      b(27)=-0.34632034632035d+01
      b(28)=-0.44097222222222d+00*dx
      b(29)=-0.14173228346457d+00
      b(30)=0.53385826771654d+01
      b(31)=-0.35905511811023d+01
      b(32)=0.70866141732284d-01
      b(33)=-0.45354330708661d+01
      b(34)=-0.31496062992126d-01
      b(35)=0.18060975609756d+01*dx
      b(36)=-0.54692775151925d-01
      b(37)=0.47967589466576d+01
      b(38)=-0.22795408507765d+01
      b(39)=0.48615800135044d-01
      b(40)=-0.34031060094530d+01
      b(41)=-0.40513166779204d-01
      b(42)=0.48615800135044d+00
      b(43)=0.48809523809524d-01*dx
      b(44)=0.65853658536585d+00
      b(45)=0.66341463414634d+01
      b(46)=0.52682926829268d+01
      i=10

100	continue
      
      return
      end
c
c
      subroutine grav(g,rho,qro,r,n)
      
c*** given rho and spline coeffs,computes gravity ***

      implicit real*8(a-h,o-z)
      dimension g(1),rho(1),qro(3,1),r(1)
      
      g(1)=0.d0
      do 10 i=2,n
      im1=i-1
      del=r(i)-r(im1)
      rn2=r(im1)*r(im1)
      trn=2.d0*r(im1)
      c1=rho(im1)*rn2
      c2=(qro(1,im1)*rn2+trn*rho(im1))*0.5d0
      c3=(qro(2,im1)*rn2+trn*qro(1,im1)+rho(im1))/3.d0
      c4=(qro(3,im1)*rn2+trn*qro(2,im1)+qro(1,im1))*.25d0
      c5=(trn*qro(3,im1)+qro(2,im1))*0.2d0
   10 g(i)=(g(im1)*rn2+4.d0*del*(c1+del*(c2+del*(c3+del*(c4+del*
     +    (c5+del*qro(3,im1)/6.d0))))))/(r(i)*r(i))
     
      return
      end
c
c
      subroutine startl(jf,jl,v,ls,q)
      
c*** Finds start level between jf and jl using velocityv and ang. ord. l.
c*** Upon entry q is the value of the exponent at r(jf) or at the turning
c*** point(q=0) depending on previous calls to startl. Upon exit q is the
c*** value of the exponent at the starting level ls.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      dimension rrlog(ndata),p(ndata),v(1)
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      data ifirst/1/
      
      if(ifirst.ne.1) goto 5
      ifirst=0
      vertno=-dlog(eps)
      do 1 i=3,n
    1 rrlog(i)=0.5d0*dlog(r(i)/r(i-1))
    5 do 10 j=jf,jl
      pp=fl3-wsq*r(j)*r(j)*rho(j)/v(j)
      if(pp.le.0.d0) goto 15
   10 p(j)=dsqrt(pp)
   15 p(j)=0.d0
   20 k=j
      j=j-1
      if(j.le.jf) go to 25
      q=q+rrlog(dble(k))*(p(j)+p(k))
      if(q.lt.vertno) go to 20
      ls=j
      goto 100
   25 ls=jf

100	continue
	
      return
      end
c
c
      subroutine steps(eps)
      
c*** computes 8 dimensionless step sizes for rks integration

      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,stepf,maxo,in
      
      ps=dlog(eps)
      fac=1.d0
      do 2 n=1,8
      fn=dble(n+1)
      fac=fac*fn
      x=(dlog(fac)+ps)/fn
      x=dexp(x)
      s=x
      do 1 i=1,n
    1 s=x*dexp(-s/fn)
    2 step(n)=s
    
      return
      end
c
c
      subroutine drspln(i1,i2,x,y,q,f)
      
c   rspln computes cubic spline interpolation coefficients
c   for y(x) between grid points i1 and i2 saving them in q.  The
c   interpolation is continuous with continuous first and second
c   derivitives.  It agrees exactly with y at grid points and with the
c   three point first derivitives at both end points (i1 and i2).
c   x must be monotonic but if two successive values of x are equal
c   a discontinuity is assumed and seperate interpolation is done on
c   each strictly monotonic segment.  The arrays must be dimensioned at
c   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
c   for rspln.
c                                                     -RPB
      implicit real*8(a-h,o-z)
      dimension x(1),y(1),q(3,1),f(3,1),yy(3)
      
      equivalence (yy(1),y0)
      
      data yy/3*0.d0/
      
      j1=i1+1
      y0=0.d0
c   bail out if there are less than two points total.
      if(i2-i1)13,17,8
 8    a0=x(j1-1)
c   search for discontinuities.
      do 3 i=j1,i2
      b0=a0
      a0=x(i)
      if(a0-b0)3,4,3
 3    continue
 17   j1=j1-1
      j2=i2-2
      go to 5
 4    j1=j1-1
      j2=i-3
c   see if there are enough points to interpolate (at least three).
 5    if(j2+1-j1)9,10,11
c   only two points.  use linear interpolation.
 10   j2=j2+2
      y0=(y(j2)-y(j1))/(x(j2)-x(j1))
      do 15 j=1,3
      q(j,j1)=yy(j)
 15   q(j,j2)=yy(j)
      go to 12
c   more than two points.  do spline interpolation.
 11   a0=0.d0
      h=x(j1+1)-x(j1)
      h2=x(j1+2)-x(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
c   calculate derivitive at near end.
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      b1=b0
c   explicitly reduce banded matrix to an upper banded matrix.
      do 1 i=j1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
 1    b0=f(3,i)
c   take care of last two rows.
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
c   calculate derivitive at far end.
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
c   solve upper banded matrix by reverse iteration.
      do 2 j=j1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
 2    i=k
      q(1,i)=b1
c   fill in the last point with a linear extrapolation.
 9    j2=j2+2
      do 14 j=1,3
 14   q(j,j2)=yy(j)
c   see if this discontinuity is the last.
 12   if(j2-i2)6,13,13
c   no.  go back for more.
 6    j1=j2+2
      if(j1-i2)8,8,7
c   there is only one point left after the latest discontinuity.
 7    do 16 j=1,3
 16   q(j,i2)=yy(j)
c   fini.

13	continue   
 
      return
      end
c
c
      subroutine dsplin(n,x,y,q,f)
      
      implicit real*8(a-h,o-z)
      dimension x(1),y(1),q(3,1),f(3,1),yy(3)
      
      equivalence (yy(1),y0)
      
      data yy/3*0.d0/
      
      a0=0.d0
      j2=n-2
      h=x(2)-x(1)
      h2=x(3)-x(1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
      b0=(y(1)*(h-h2)+y(2)*h2-y(3)*h)/y0
      b1=b0
      do 5 i=1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.d0*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
    5 b0=f(3,i)
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
      do 10 j=1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
   10 i=k
      q(1,i)=b1
      do 15 j=1,3
   15 q(j,n)=yy(j)
   
      return
      end
c
c
      subroutine trknt(y1,y1p,y2,y2p,x,y)
      
c*** toroidal and radial mode counter ***

      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension xs(2),val(4)
      
      ns=0
      if(kount.eq.0) goto 60
      h=y-x
      a1=h*y1p
      a2=-h*(2.d0*y1p+y2p)+3.d0*(y2-y1)
      a3=h*(y1p+y2p)-2.d0*(y2-y1)
      a33=3.d0*a3
      a22=2.d0*a2
      if(a3.ne.0.d0) goto 20
      if(a2.eq.0.d0) goto 50
      xs(1)=-a1/a22
      if(xs(1).ge.0.d0.and.xs(1).le.1.d0) ns=1
      goto 50
   20 disc=a2*a2-a1*a33
      if(disc) 50,25,30
   25 xs(1)=-a2/a33
      if(xs(1).ge.0.d0.and.xs(1).le.1.d0) ns=1
      goto 50
   30 disc=dsqrt(disc)
      tr1=(-a2+disc)/a33
      tr2=(-a2-disc)/a33
      if(dabs(a33).gt.dabs(a1)) goto 35
      fac=a1/a33
      tr1=fac/tr1
      tr2=fac/tr2
   35 if(tr1.lt.0.d0.or.tr1.gt.1.d0) goto 40
      xs(1)=tr1
      ns=1
   40 if(tr2.lt.0.d0.or.tr2.gt.1.d0) goto 50
      ns=ns+1
      xs(ns)=tr2
      if(ns.lt.2) goto 50
      if(tr2.ge.tr1) goto 50
      xs(1)=tr2
      xs(2)=tr1
   50 if(ns.eq.0) goto 60
      ns1=ns+1
      do 55 j=2,ns1
      t=xs(j-1)
   55 val(j)=y1+t*(a1+t*(a2+t*a3))
   60 val(1)=y1
      ns2=ns+2
      val(ns2)=y2
      do 100 j=2,ns2
  100 if(val(j-1)*val(j).le.0.d0) kount=kount+1
      if(val(1).eq.0.d0) kount=kount-1
      
      return
      end
c
c
      subroutine rprop(jf,jl,f)
      
c*** propagates soln ,f, for radial modes from jf to jl ***

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,ndata),nnnn(ndata)
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,stepf,maxo,in
      dimension h(2,10),s(2),f(2)
      
      maxo1=maxo-1
      y=r(jf)
      vy=dsqrt((flam(jf)+2.d0*fmu(jf))/rho(jf))
      i=jf
      go to 50
   10 iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 50
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      vx=vy
      vy=dsqrt((flam(i)+2.d0*fmu(i))/rho(i))
      q=dmax1(w/vx+1.d0/x,w/vy+1.d0/y)
      del=step(maxo)/q
      dxs=0.d0
   15 y=x+del
      if(y.gt.r(i)) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,maxo1)
      dxs=dx
      s(1)=f(1)
      s(2)=f(2)
      do 40 ni=1,in
      z=x+c(ni)
      t=z-r(iq)
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      if(ifanis.ne.0) goto 30
      nn=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      cc=ff+nn+nn
      aa=cc
      goto 35
   30 nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
      aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
   35 z=1.d0/z
      a21=-ro*wsq+4.d0*z*(z*(aa-nn-ff*ff/cc)-ro*gr)
      h(1,ni)=(f(2)-2.d0*ff*z*f(1))/cc
      h(2,ni)=a21*f(1)+2.d0*z*f(2)*(ff/cc-1.d0)
   40 call rkdot(f,s,h,2,ni)
      if(knsw.ne.1) goto 45
      fp=a21*f(1)+2.d0*z*f(2)*(ff/cc-1.d0)
      call trknt(s(2),h(2,1),f(2),fp,x,y)
   45 x=y
      if(y.ne.r(i)) go to 15
   50 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.ne.jl) go to 10
      
      return
      end
c
c
      subroutine tprop(jf,jl,f)
      
c*** propagates f from jf to jl - toroidal modes ***

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/a(14,ndata),nnnn(ndata)
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,stepf,maxo,in
      dimension h(2,10),s(2),f(2)
      
      fl3m2=fl3-2.d0
      maxo1=maxo-1
      y=r(jf)
      vy=fmu(jf)/rho(jf)
      i=jf
      go to 50
   10 iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 50
      qll=1.d0+qshear(iq)*fct
      vx=vy
      vy=fmu(i)/rho(i)
      qx=1.d0/x+dsqrt(dabs(wsq/(vx)-fl3/(x*x)))
      qy=1.d0/y+dsqrt(dabs(wsq/(vy)-fl3/(y*y)))
      q=dmax1(qx,qy)
      del=step(maxo)/q
      dxs=0.d0
   15 y=x+del
      if(y.gt.r(i)) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,maxo1)
      dxs=dx
      s(1)=f(1)
      s(2)=f(2)
      do 40 ni=1,in
      z=x+c(ni)
      t=z-r(iq)
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      nn=ll
      if(ifanis.ne.0) nn=(ncon(iq)+
     +    t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      z=1.d0/z
      h(1,ni)=z*f(1)+f(2)/ll
      h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
   40 call rkdot(f,s,h,2,ni)
      if(knsw.ne.1) goto 45
      fp=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
      call trknt(s(2),h(2,1),f(2),fp,x,y)
   45 x=y
      if(y.ne.r(i)) goto 15
   50 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.ne.jl) go to 10
      
100	continue

      return
      end
c
c
      subroutine rkdot(f,s,h,nvec,ni)
      
c*** performs dot product with rks coefficients ***

      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,stepf,maxo,in
      dimension s(1),f(1),h(nvec,1)
      
      goto (1,2,3,4,5,6,7,8,9,10),ni
    1 do 21 j=1,nvec
   21 f(j)=s(j)+b(1)*h(j,1)
      goto 100
    2 do 22 j=1,nvec
   22 f(j)=s(j)+b(2)*(h(j,1)+b(3)*h(j,2))
      goto 100
    3 do 23 j=1,nvec
   23 f(j)=s(j)+b(4)*(h(j,1)+b(5)*h(j,2)+b(6)*h(j,3))
      goto 100
    4 do 24 j=1,nvec
   24 f(j)=s(j)+b(7)*(h(j,1)+b(8)*h(j,2)+b(9)*h(j,3)+b(10)*h(j,4))
      goto 100
    5 do 25 j=1,nvec
   25 f(j)=s(j)+b(11)*(h(j,1)+b(12)*h(j,2)+b(13)*h(j,3)+b(14)*h(j,4)+
     +b(15)*h(j,5))
      goto 100
    6 do 26 j=1,nvec
   26 f(j)=s(j)+b(16)*(h(j,1)+b(17)*h(j,2)+b(18)*h(j,3)+b(19)*h(j,4)+
     +b(20)*h(j,5)+b(21)*h(j,6))
      goto 100
    7 do 27 j=1,nvec
   27 f(j)=s(j)+b(22)*(h(j,1)+b(23)*h(j,3)+b(24)*h(j,4)+b(25)*h(j,5)+
     +b(26)*h(j,6)+b(27)*h(j,7))
      goto 100
    8 do 28 j=1,nvec
   28 f(j)=s(j)+b(28)*(h(j,1)+b(29)*h(j,3)+b(30)*h(j,4)+b(31)*h(j,5)+
     +b(32)*h(j,6)+b(33)*h(j,7)+b(34)*h(j,8))
      goto 100
    9 do 29 j=1,nvec
   29 f(j)=s(j)+b(35)*(h(j,1)+b(36)*h(j,3)+b(37)*h(j,4)+b(38)*h(j,5)+
     +b(39)*h(j,6)+b(40)*h(j,7)+b(41)*h(j,8)+b(42)*h(j,9))
      goto 100
   10 do 30 j=1,nvec
   30 f(j)=s(j)+b(43)*(h(j,1)+h(j,10)+b(44)*(h(j,4)+h(j,6))+
     +b(45)*h(j,5)+b(46)*(h(j,7)+h(j,9)))

100	continue
     
      return
      end
c
c
      subroutine intgds(rr,iq,vals)
      
c*** interpolates integrands for normalisation,cg,q etc..for use with gauslv.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,ndata),nnnn(ndata)
      dimension q(3),qp(3),vals(1)
      data d1,d2,d3,d4,d5,d6,d7/0.111111111111111d0,
     + 0.066666666666667d0,0.666666666666667d0,1.333333333333333d0,
     + 2.666666666666667d0,3.333333333333333d0,5.333333333333333d0/
     
      t=rr-r(iq)
      hn=1.d0/(r(iq+1)-r(iq))
      hsq=hn*hn
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      iq1=iq+1
      ifun=3
      if(jcom.ne.3) ifun=1
      do 10 i=1,ifun
      i2=2*i
      i1=i2-1
      a=((ar(i2,iq)+ar(i2,iq1))+2.d0*hn*(ar(i1,iq)-ar(i1,iq1)))*hsq
      b=-(2.d0*ar(i2,iq)+ar(i2,iq1))*hn-3.d0*(ar(i1,iq)-ar(i1,iq1))*hsq
      q(i)=(ar(i1,iq)+t*(ar(i2,iq)+t*(b+t*a)))/rr
   10 qp(i)=ar(i2,iq)+t*(2.d0*b+t*3.d0*a)
      rro=(rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq))))*rr
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      if(ifanis.ne.0) goto 15
      nn=ll
      cc=ff+ll+ll
      aa=cc
      goto 20
   15 qaa=1.d0+xa2(iq)*fct
      nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
      aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
   20 qrka=d1*(4.d0*(aa+ff-nn)+cc)
     1     *(qkappa(iq)+t*hn*(qkappa(iq1)-qkappa(iq)))
      qrmu=d2*(aa+cc-2.d0*ff+5.d0*nn+6.d0*ll)
     1     *(qshear(iq)+t*hn*(qshear(iq1)-qshear(iq)))
      if(jcom.ne.3) goto 25
      q1sq=q(1)*q(1)
      q2sq=q(2)*q(2)
      vals(1)=rr*rro*(q1sq+q2sq)
c      fac=(fl+.5d0)/sfl3
      fac=sfl3/sfl3
      vals(2)=(sfl3*(ll*q1sq+aa*q2sq)+q(2)*((rro*gr+2.d0*(nn-aa-ll)+ff)
     +   *q(1)+rro*q(3)-ff*qp(1))+ll*qp(2)*q(1))*fac
     +   +0.25d0*q(3)*(qp(3)+fl*q(3))
      t2=qrka+d7*qrmu
      t3=qrka+d4*qrmu
      t4=qrka+d6*qrmu
      t5=qrka-d5*qrmu
      t6=qrka-d3*qrmu
      vals(3)=0.5d0*((fl3*qrmu+t2)*q1sq+(2.d0*qrmu+fl3*t3)*q2sq)
     1 -q(1)*sfl3*t4*q(2)+q(1)*(t5*qp(1)+sfl3*qrmu*qp(2))+q(2)*(-2.d0*
     2 qrmu*qp(2)-sfl3*t6*qp(1))+0.5d0*(t3*qp(1)*qp(1)+qrmu*qp(2)*qp(2))
      vals(4)=0.5d0*((fl3*ll+4.d0*(rro*(rro-gr)+aa-nn-ff)+cc)*q1sq+
     +(4.d0*ll-nn-nn+fl3*aa)*q2sq +fl*fl*0.25d0*q(3)*q(3)+cc*qp(1)*qp(1)+
     +ll*qp(2)*qp(2)+0.25d0*qp(3)*qp(3))+q(3)*(rro*sfl3*q(2)+fl*0.25d0*qp
     +(3))+q(1)*(sfl3*(rro*gr+2.d0*(nn-aa-ll)+ff)*q(2)+rro*(qp(3)-q(3))+
     +(ff+ff-cc)*qp(1)+sfl3*ll*qp(2))-q(2)*(sfl3*ff*qp(1)+(ll+ll)*qp(2))
      return
   25 q(1)=q(1)*rr
      vals(1)=rr*rro*q(1)*q(1)
      if(jcom.eq.1) goto 30
      vals(2)=nn*q(1)*q(1)
      t1=(rr*qp(1)-q(1))**2
      t2=(fl3-2.d0)*q(1)*q(1)
      vals(3)=(t1+t2)*qrmu
      vals(4)=t1*ll+t2*nn
      return
   30 t1=(rr*qp(1)+2.d0*q(1))**2
      t2=d4*(rr*qp(1)-q(1))**2
      vals(2)=t1*qrka+t2*qrmu
      vals(3)=rr*qp(1)*(cc*rr*qp(1)+4.d0*ff*q(1))+4.d0*q(1)*q(1)
     +    *(aa-nn-rro*gr)
     
      return
      end
c
c
      subroutine fpsm(ls,nvefm,ass)
      
c*** Spheroidal mode start solution in a fluid region using sph. Bessel fns.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension ass(1)
      
      x=r(ls)
      fla=flam(ls)*(1.d0+xlam(ls)*fct)
      vpsq=fla/rho(ls)
      xi=g(ls)/x
      qsq=(wsq+dble(kg)*4.d0*rho(ls)+xi-fl3*xi*xi/wsq)/vpsq
      zsq=qsq*x*x
      call bfs(l,zsq,eps,fp)
      if(kg.eq.0) goto 20
      u=(fl-fp)/qsq
      c1=fl*g(ls)-wsq*x
      c2=fl2*c1*0.25d0/x-rho(ls)*fl
      ass(1)=-x*fl*vpsq-c1*u
      ass(2)=-x*fl*fla
      ass(3)=-fl*fl2*vpsq*0.25d0-u*c2
      ass(4)=x*fla*c1
      ass(5)=-x*fla*c2
      goto 25
   20 ass(1)=-(fl3*xi/wsq+fp)/qsq
      ass(2)=x*fla
   25 sum=ass(1)*ass(1)
      do 30 i=2,nvefm
   30 sum=sum+ass(i)*ass(i)
      sum=1.d0/dsqrt(sum)
      if(ass(nvefm).lt.0.d0) sum=-sum
      do 35 i=1,nvefm
   35 ass(i)=ass(i)*sum
   
      return
      end
c
c
      subroutine spsm(ls,nvesm,ass)
      
c*** Spheroidal mode start solution in a solid region using sph. Bessel fns.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(6,2),e(15),ass(1)
      
      x=r(ls)
      ro=rho(ls)
      fu=fmu(ls)*(1.d0+qshear(ls)*fct)
      flu=flam(ls)*(1.d0+xlam(ls)*fct)+2.d0*fu
      vssq=fu/ro
      vpsq=flu/ro
      zeta=4.d0*ro
      xi=g(ls)/x
      alfsq=(wsq+dble(kg)*zeta+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vpsq*vssq))
      fksq=0.5d0*(alfsq+betasq+delsq)
      qsq=fksq-delsq
      zsq=qsq*x*x
      b=xi/(vssq*(betasq-qsq))
      k=1
    5 call bfs(l,zsq,eps,fp)
      a(1,k)=fl3*b+fp
      a(2,k)=1.d0+b+b*fp
      a(3,k)=-zsq
      a(4,k)=b*a(3,k)
      a(5,k)=1.d0
      a(6,k)=fl1-fl3*b
      if(k.eq.2) goto 10
      zsq=fksq*x*x
      b=-flu/(fu*fl3*b)
      k=2
      goto 5
   10 jj=3+2*kg
      kk=jj+1
      ll=0
      do 15 i=1,jj
      i1=i+1
      do 15 j=i1,kk
      ll=ll+1
   15 e(ll)=a(i,1)*a(j,2)-a(j,1)*a(i,2)
      if(kg.ne.0) goto 20
      ass(1)=x*x*e(1)
      ass(2)=fu*x*sfl3*(2.d0*e(1)-e(5))
      ass(3)=fu*x*(e(3)-2.d0*e(1))
      ass(4)=x*(flu*e(4)+4.d0*fu*e(1))
      ass(5)=fu*(flu*(e(6)+2.d0*e(4))+4.d0*fu*(fl3*(e(5)-e(1))
     +     -e(3)+2.d0*e(1)))
      goto 25
   20 c0=wsq-xi*fl
      c1=ro*fl+0.25d0*fl2*c0
      c2=2.d0*fu/x
      c3=c2*(fl-1.d0)
      ass(6)=x*x*(c0*e(1)-zeta*(fl*e(8)-e(4)))
      ass(14)=flu*(fl*e(6)-e(2))
      ass(13)=fu*sfl3*(fl*e(7)-e(3))
      ass(1)=x*(c1*e(1)-ro*(fl*e(9)-e(5)))
      ass(7)=x*flu*(c0*e(2)-zeta*fl*e(11))/sfl3+c2*sfl3*ass(6)
      ass(8)=x*fu*(c0*e(3)-zeta*fl*e(13))-c2*ass(6)
      ass(12)=(flu*fl*e(10)+2.d0*(ass(14)+sfl3*ass(13)))*fu/x
      ass(2)=flu*(c1*e(2)-ro*fl*e(12))/sfl3+c2*sfl3*ass(1)
      ass(3)=fu*(c1*e(3)-ro*fl*e(14))-c2*ass(1)
      ass(9)=(x*c0*ass(14)+sfl3*ass(7)-c3*fl*ass(6))/fl
      ass(11)=(sfl3*ass(12)+c3*(sfl3*ass(14)-fl*ass(13)))/fl
      ass(4)=(c1*ass(14)+sfl3*ass(2)-c3*fl*ass(1))/fl
      ass(10)=(x*c0*ass(11)-c3*(sfl3*ass(9)+fl*ass(7)))/sfl3
      ass(5)=(c1*ass(11)-c3*(sfl3*ass(4)+fl*ass(2)))/sfl3
   25 sum=ass(1)*ass(1)
      do 30 i=2,nvesm
   30 sum=sum+ass(i)*ass(i)
      sum=1.d0/dsqrt(sum)
      if(ass(5).lt.0.d0) sum=-sum
      do 35 i=1,nvesm
   35 ass(i)=ass(i)*sum
   
      return
      end
c
c
      subroutine rps(i,a)
      
c*** radial mode start soln using sph bessel fns.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(2)
      
      fla=flam(i)*(1.d0+xlam(i)*fct)
      sig=fla+2.d0*fmu(i)*(1.d0+qshear(i)*fct)
      zsq=r(i)*r(i)*rho(i)*(wsq+4.d0*g(i)/r(i))/sig
      call bfs(1,zsq,eps,fp)
      a(1)=r(i)
      a(2)=sig*fp+2.d0*fla
      
      return
      end
c
c
      subroutine tps(i,a)
      
c*** toroidal mode start soln using sph bessel fns.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      dimension a(2)
      
      fu=fmu(i)*(1.d0+qshear(i)*fct)
      zsq=r(i)*r(i)*wsq*rho(i)/fu
      call bfs(l,zsq,eps,fp)
      a(1)=r(i)
      a(2)=fu*(fp-1.d0)
      
      return
      end
c
c
      subroutine bfs(l,xsq,eps,fp)
      
c  This routine calculates spherical bessel function of the ist kind.
c  fp is equivalent to (r*dj/dr)/j
c  where r is radius and j is the sbf of order l and argument x=k*r
c  the technique employs the continued fraction approach
c  described in W. Lentz's article in applied qptics, vol.15, #3, 1976

      implicit real*8(a-h,o-z)
      real*8 numer,nu
      
      if(xsq.le.0.d0) goto 10
      x=dsqrt(xsq)
      lp1=l+1
      rx=2.0d0/x
      nu=dble(lp1)-0.5d0
      rj=nu*rx
      rx=-rx
      denom=(nu+1.d0)*rx
      numer=denom+1.0d0/rj
      rj=rj*numer/denom
      nm1=1
    2 nm1=nm1+1
      rx=-rx
      a3=(nu+dble(nm1))*rx
      denom=a3+1.d0/denom
      numer=a3+1.d0/numer
      ratio=numer/denom
      rj=rj*ratio
      if(dabs(dabs(ratio)-1.d0).gt.eps) goto 2
      fp=rj*x-dble(lp1)
      goto 100
c  series solution
   10 f=1.d0
      fp=dble(l)
      a=1.d0
      b=dble(l+l)+1.d0
      c=2.d0
      d=dble(l)+2.d0
   15 a=-a*xsq/(c*(b+c))
      f=f+a
      fp=fp+a*d
      if(dabs(a*d).lt.eps) goto 20
      c=c+2.d0
      d=d+2.d0
      goto 15
   20 fp=fp/f

100	continue
   
      return
      end
c
c
      subroutine remedy(ls)
      
c    Obtains the eigenfunction of an awkward spheroidal mode by
c    integrating to the icb or the mcb.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,ndata),inorm(ndata)
      common/arem/a(6,3,ndata)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension af(4,2),as(6,3),afr(4)
      
      rms1=0.d0
      rms2=0.d0
c      write(*,900) ls,noc
      if(ls.gt.noc) goto 100
      iexp=0
      do 10 k=1,2
      do 10 j=1,4
   10 af(j,k)=0.d0
      af(1,1)=1.d0
      if(kg.eq.1) af(2,2)=1.d0
      if(nsl.eq.n) goto 5
      do 6 i=nslp1,n
      do 6 k=1,3
      do 6 j=1,6
    6 a(j,k,i)=0.d0
      call fprop(n,nslp1,af,iexp)
    5 call fsbdry(af,as,kg)
      do 7 k=1,3
      do 7 j=1,6
    7 a(j,k,nsl)=as(j,k)
      if(n.ne.nsl) call ortho(n,nsl,as,kg)
      call sprop(n,nsl,nocp1,as,iexp)
      call sfbdry(n,nocp1,as,af,kg)
      imtch=noc
      do 11 i=1,4
   11 afr(i)=ar(i,noc)
      if(ls.gt.nic) goto 15
      icomp=0
      call match(n,noc,kg,af,afr,icomp,rms1)
      if(icomp.eq.0) goto 100
      call fprop(noc,nicp1,af,iexp)
      imtch=nic
      do 12 i=1,4
   12 afr(i)=ar(i,nicp1)
   15 icomp=-1
      call match(n,imtch,kg,af,afr,icomp,rms2)
100	continue
c      write(6,900) ls,rms1,rms2
      call flush(6)
900   format('    Remedy with start level : ',i6,2(1x,e12.5))

      return
      end
c
c
      subroutine fprop(jf,jl,f,iexp)
      
c    fprop propagates the fundamental matrix f from jf to jl (a fluid region)

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,ndata),inorm(ndata)
      common/arem/a(6,3,ndata)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,sdum,idum,in
      dimension f(4,2),s(4,2),h(4,2,10)
      data econst/1048576.d0/
      
      kk=kg+1
      jj=2*kk
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 80
   10 x=y
      y=r(i)
      if(y.eq.x) goto 80
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=dmax1(sfl3/x,dsqrt(dabs(alfsq-fl3/(x*x)))+1.d0/zs)
      del=dble(jud)*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(dble(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,7)
      dxs=dx
      do 20 k=1,kk
      do 20 j=1,jj
   20 s(j,k)=f(j,k)
      d=fl3/wsq
      do 40 ni=1,in
      z=x+c(ni)
      t=z-zs
      zr=1.d0/z
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      flu=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      gr=(g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq))))*zr
      t21=-4.d0*ro
      t12=d*zr*zr
      t11=(gr*d-1.d0)*zr
      s11=-ro*(wsq+4.d0*gr-gr*gr*d)
      c11=-t12/ro+1.d0/flu
      if(kg.eq.0) s11=s11-t21*ro
      if(kg.eq.0) goto 25
      t22=-fl*zr
      s22=ro*t12
      s12=ro*(t11+t22)
   25 do 70 k=1,kk
      if(kg.ne.0) goto 30
      h(1,k,ni)=t11*f(1,k)+c11*f(2,k)
      h(2,k,ni)=s11*f(1,k)-t11*f(2,k)
      goto 35
   30 h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(3,k)
      h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+4.d0*f(4,k)
      h(3,k,ni)=s11*f(1,k)+s12*f(2,k)-t11*f(3,k)-t21*f(4,k)
      h(4,k,ni)=s12*f(1,k)+s22*f(2,k)-t12*f(3,k)-t22*f(4,k)
   35 do 70 j=1,jj
      go to (701,702,703,704,705,706,707,708,709,710),ni
  701 f(j,k)=s(j,k)+b(1)*h(j,k,1)
      go to 70
  702 f(j,k)=s(j,k)+b(2)*(h(j,k,1)+b(3)*h(j,k,2))
      go to 70
  703 f(j,k)=s(j,k)+b(4)*(h(j,k,1)+b(5)*h(j,k,2)+b(6)*h(j,k,3))
      go to 70
  704 f(j,k)=s(j,k)+b(7)*(h(j,k,1)+b(8)*h(j,k,2)+b(9)*h(j,k,3)+
     +b(10)*h(j,k,4))
      go to 70
  705 f(j,k)=s(j,k)+b(11)*(h(j,k,1)+b(12)*h(j,k,2)+b(13)*h(j,k,3)+
     +b(14)*h(j,k,4)+b(15)*h(j,k,5))
      go to 70
  706 f(j,k)=s(j,k)+b(16)*(h(j,k,1)+b(17)*h(j,k,2)+b(18)*h(j,k,3)+
     +b(19)*h(j,k,4)+b(20)*h(j,k,5)+b(21)*h(j,k,6))
      go to 70
  707 f(j,k)=s(j,k)+b(22)*(h(j,k,1)+b(23)*h(j,k,3)+b(24)*h(j,k,4)+
     +b(25)*h(j,k,5)+b(26)*h(j,k,6)+b(27)*h(j,k,7))
      go to 70
  708 f(j,k)=s(j,k)+b(28)*(h(j,k,1)+b(29)*h(j,k,3)+b(30)*h(j,k,4)+
     +b(31)*h(j,k,5)+b(32)*h(j,k,6)+b(33)*h(j,k,7)+b(34)*h(j,k,8))
      go to 70
  709 f(j,k)=s(j,k)+b(35)*(h(j,k,1)+b(36)*h(j,k,3)+b(37)*h(j,k,4)+
     +b(38)*h(j,k,5)+b(39)*h(j,k,6)+b(40)*h(j,k,7)+b(41)*h(j,k,8)+
     +b(42)*h(j,k,9))
      go to 70
  710 f(j,k)=s(j,k)+b(43)*(h(j,k,1)+h(j,k,10)+b(45)*h(j,k,5)+
     +b(44)*(h(j,k,4)+h(j,k,6))+b(46)*(h(j,k,7)+h(j,k,9)))
   70 continue
   40 continue
      x=y
      if(y.ne.r(i)) go to 15
   80 size=0.d0
      do 81 k=1,kk
      do 81 j=1,jj
   81 size=dmax1(size,dabs(f(j,k)))
   82 if(size.lt.1024.d0) goto 84
      do 83 k=1,kk
      do 83 j=1,jj
   83 f(j,k)=f(j,k)/econst
      size=size/econst
      iexp=iexp+20
      goto 82
   84 inorm(i)=iexp
      do 85 k=1,kk
      do 85 j=1,jj
   85 a(j,k,i)=f(j,k)
      if(i.eq.jl) goto 100
      i=i+jud
      go to 10

100	continue
      
      return
      end
c
c
      subroutine sprop(li,jf,jl,f,iexp)
      
c    sprop propagates the fundamental matrix f from jf to jl (a solid region)
c    if iorth=1 the columns of f are orthogonalized at each level
c    except in regions of oscillatory p and s.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(ndata),fmu(ndata),flam(ndata),qshear(ndata),
     +qkappa(ndata),xa2(ndata),xlam(ndata),rho(ndata),
     +qro(3,ndata),g(ndata),qg(3,ndata),fcon(ndata),
     +fspl(3,ndata),lcon(ndata),lspl(3,ndata),ncon(ndata),
     +nspl(3,ndata),ccon(ndata),cspl(3,ndata),acon(ndata),
     +aspl(3,ndata)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,wgrav,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,wdiff,jcom,nord,l,kg,kount,knsw,ifanis,iback
      common/eifx/ar(14,ndata),inorm(ndata)
      common/arem/a(6,3,ndata)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/shanks/b(46),c(10),step(8)
      common/shank/dx,sdum,idum,in
      dimension f(6,3),s(6,3),h(6,3,10)
      data econst/1048576.d0/
      
      kk=kg+2
      jj=2*kk
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
      go to 80
   10 x=y
      y=r(i)
      if(x.eq.y) goto 80
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=0.5d0*(alfsq+betasq+delsq)
      al=fl3/(x*x)
      jorth=1
      aq=fksq-delsq-al
      if(aq.gt.0.d0) jorth=0
      qs=dsqrt(dabs(fksq-al))+1.d0/zs
      qf=dsqrt(dabs(aq))+1.d0/zs
      q=dmax1(sfl3/x,qs,qf)
      del=dble(jud)*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(dble(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,7)
      dxs=dx
      do 20 k=1,kk
      do 20 j=1,jj
   20 s(j,k)=f(j,k)
      do 50 ni=1,in
      z=x+c(ni)
      t=z-zs
      ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
      gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
      ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      if(ifanis.ne.0) goto 25
      nn=ll
      cc=ff+ll+ll
      aa=cc
      goto 30
   25 nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
      cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
      aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
   30 zr=1.d0/z
      sfl3z=sfl3*zr
      rogr=ro*gr
      c11=1.d0/cc
      c22=1.d0/ll
      dmg=aa-nn-ff*ff*c11
      zdmg=zr*dmg
      t11=-2.d0*ff*zr*c11+zr
      t12=sfl3z*ff*c11
      t21=-sfl3z
      t22=zr+zr
      s22=-ro*wsq
      s11=s22+4.d0*zr*(zdmg-rogr)
      s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn)
      s12=sfl3z*(rogr-zdmg-zdmg)
      if(kg.eq.0) s11=s11+4.d0*ro*ro
      if(kg.eq.0) goto 35
      t31=-4.d0*ro
      t33=-fl*zr
      s13=-fl1*zr*ro
      s23=ro*sfl3z
   35 do 70 k=1,kk
      if(kg.eq.1) goto 40
      h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(3,k)
      h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+c22*f(4,k)
      h(3,k,ni)=s11*f(1,k)+s12*f(2,k)-t11*f(3,k)-t21*f(4,k)
      h(4,k,ni)=s12*f(1,k)+s22*f(2,k)-t12*f(3,k)-t22*f(4,k)
      goto 45
   40 h(1,k,ni)=t11*f(1,k)+t12*f(2,k)+c11*f(4,k)
      h(2,k,ni)=t21*f(1,k)+t22*f(2,k)+c22*f(5,k)
      h(3,k,ni)=t31*f(1,k)+t33*f(3,k)+4.d0*f(6,k)
      h(4,k,ni)=s11*f(1,k)+s12*f(2,k)+s13*f(3,k)-t11*f(4,k)-t21*f(5,k)
     +    -t31*f(6,k)
      h(5,k,ni)=s12*f(1,k)+s22*f(2,k)+s23*f(3,k)-t12*f(4,k)-t22*f(5,k)
      h(6,k,ni)=s13*f(1,k)+s23*f(2,k)-t33*f(6,k)
   45 do 70 j=1,jj
      go to (701,702,703,704,705,706,707,708,709,710),ni
  701 f(j,k)=s(j,k)+b(1)*h(j,k,1)
      go to 70
  702 f(j,k)=s(j,k)+b(2)*(h(j,k,1)+b(3)*h(j,k,2))
      go to 70
  703 f(j,k)=s(j,k)+b(4)*(h(j,k,1)+b(5)*h(j,k,2)+b(6)*h(j,k,3))
      go to 70
  704 f(j,k)=s(j,k)+b(7)*(h(j,k,1)+b(8)*h(j,k,2)+b(9)*h(j,k,3)+
     +b(10)*h(j,k,4))
      go to 70
  705 f(j,k)=s(j,k)+b(11)*(h(j,k,1)+b(12)*h(j,k,2)+b(13)*h(j,k,3)+
     +b(14)*h(j,k,4)+b(15)*h(j,k,5))
      go to 70
  706 f(j,k)=s(j,k)+b(16)*(h(j,k,1)+b(17)*h(j,k,2)+b(18)*h(j,k,3)+
     +b(19)*h(j,k,4)+b(20)*h(j,k,5)+b(21)*h(j,k,6))
      go to 70
  707 f(j,k)=s(j,k)+b(22)*(h(j,k,1)+b(23)*h(j,k,3)+b(24)*h(j,k,4)+
     +b(25)*h(j,k,5)+b(26)*h(j,k,6)+b(27)*h(j,k,7))
      go to 70
  708 f(j,k)=s(j,k)+b(28)*(h(j,k,1)+b(29)*h(j,k,3)+b(30)*h(j,k,4)+
     +b(31)*h(j,k,5)+b(32)*h(j,k,6)+b(33)*h(j,k,7)+b(34)*h(j,k,8))
      go to 70
  709 f(j,k)=s(j,k)+b(35)*(h(j,k,1)+b(36)*h(j,k,3)+b(37)*h(j,k,4)+
     +b(38)*h(j,k,5)+b(39)*h(j,k,6)+b(40)*h(j,k,7)+b(41)*h(j,k,8)+
     +b(42)*h(j,k,9))
      go to 70
  710 f(j,k)=s(j,k)+b(43)*(h(j,k,1)+h(j,k,10)+b(45)*h(j,k,5)+
     +b(44)*(h(j,k,4)+h(j,k,6))+b(46)*(h(j,k,7)+h(j,k,9)))
   70 continue
   50 continue
      x=y
      if(y.ne.r(i)) go to 15
   80 size=0.d0
      do 81 k=1,kk
      do 81 j=1,jj
   81 size=dmax1(size,dabs(f(j,k)))
   82 if(size.lt.1024.d0) goto 84
      do 83 k=1,kk
      do 83 j=1,jj
   83 f(j,k)=f(j,k)/econst
      size=size/econst
      iexp=iexp+20
      goto 82
   84 inorm(i)=iexp
      do 85 k=1,kk
      do 85 j=1,jj
   85 a(j,k,i)=f(j,k)
      if(jorth.eq.1) call ortho(li,i,f,kg)
      if(i.eq.jl) goto 100
      i=i+jud
      go to 10

100	continue
      
      return
      end
c
c
      subroutine sfbdry(jf,jl,as,af,kg)
      
c*** The tangential traction scalar is forced to vanish at the solid
c*** side of a s/f boundary(level jl). a(j,3,i) is elliminated for
c*** i=jf...jl and af is loaded from a at level jl.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      common/arem/a(6,3,ndata)
      dimension as(6,1),af(4,1)
      
      n1=min0(jf,jl)
      n2=max0(jf,jl)
      if(kg.ne.0) goto 25
      i1=1
      i2=2
      if(dabs(as(4,2)).gt.dabs(as(4,1))) goto 10
      i1=2
      i2=1
   10 rat=-as(4,i1)/as(4,i2)
      do 15 i=n1,n2
      do 15 j=1,4
   15 a(j,1,i)=a(j,i1,i)+rat*a(j,i2,i)
      af(1,1)=a(1,1,jl)
      af(2,1)=a(3,1,jl)
      goto 100
   25 ab53=dabs(as(5,3))
      do 30 k=1,2
      i1=k
      i2=3
      if(ab53.gt.dabs(as(5,k))) goto 35
      i1=3
      i2=k
   35 rat=-as(5,i1)/as(5,i2)
      do 40 i=n1,n2
      do 40 j=1,6
   40 a(j,k,i)=a(j,i1,i)+rat*a(j,i2,i)
      af(1,k)=a(1,k,jl)
      af(2,k)=a(3,k,jl)
      af(3,k)=a(4,k,jl)
   30 af(4,k)=a(6,k,jl)
   
100	continue

      return
      end
c
c
      subroutine fsbdry(af,as,kg)
      
c    fsbdry creates solid fundamental matrix as from fluid fundamental matrix
c    af. It is presumed that fsbdry is used to cross a f/s boundary.

      implicit real*8(a-h,o-z)
      dimension af(4,1),as(6,1)
      
      do 10 i=1,3
      do 10 j=1,6
   10 as(j,i)=0.d0
      if(kg.ne.0) goto 20
      as(1,1)=af(1,1)
      as(3,1)=af(2,1)
      as(2,2)=1.d0
      goto 100
   20 do 25 k=1,2
      as(1,k)=af(1,k)
      as(3,k)=af(2,k)
      as(4,k)=af(3,k)
   25 as(6,k)=af(4,k)
      as(2,3)=1.d0
      
100	continue

      return
      end
c
c
      subroutine match(n,j,kg,af,afr,icomp,rms)
      
      include'mineos.inc'
      implicit real*8(a-h,o-z)
      common/eifx/ar(14,ndata),inorm(ndata)
      common/arem/a(6,3,ndata)
      dimension af(4,1),afr(1),afi(4)
      
      k=j+2
      rms=0.d0
      fnor=0.d0
      if(kg.eq.1) go to 20
      c=(af(1,1)*afr(1)+af(2,1)*afr(2))/(af(1,1)**2+af(2,1)**2)
      do 5 i=1,2
      afi(i)=af(i,1)*c
      rms=rms+(afi(i)-afr(i))**2
    5 fnor=fnor+afr(i)*afr(i)
      rms=dsqrt(rms/fnor)
      if(icomp.lt.0) goto 6
cxxxxxxxxx rms was 1.d-3.
      if(rms.lt.1.d-6) goto 6
      icomp=1
      goto 100
    6 idiff=inorm(j)-inorm(j+1)
      inorm(j+1)=inorm(j)
  999 format(4g20.10)
      do 10 i=k,n
      inorm(i)=inorm(i)+idiff
      do 10 jj=1,4
   10 ar(jj,i)=c*a(jj,1,i)
      goto 100
   20 a2=(af(3,1)*afr(1)-af(1,1)*afr(3))/(af(1,2)*af(3,1)-af(1,1)
     +   *af(3,2))
      a1=(af(3,2)*afr(1)-af(1,2)*afr(3))/(af(1,1)*af(3,2)-af(3,1)
     +   *af(1,2))
      do 21 i=1,4
      afi(i)=a1*af(i,1)+a2*af(i,2)
      rms=rms+(afi(i)-afr(i))**2
   21 fnor=fnor+afr(i)*afr(i)
      rms=dsqrt(rms/fnor)
c     print 999,rms
      if(icomp.lt.0) goto 22
      if(rms.lt.1.d-15) goto 22
      icomp=1
      goto 100
   22 idiff=inorm(j)-inorm(j+1)
      inorm(j+1)=inorm(j)
      do 25 i=k,n
      inorm(i)=inorm(i)+idiff
      do 25 jj=1,6
   25 ar(jj,i)=a1*a(jj,1,i)+a2*a(jj,2,i)
   
100	continue

      return
      end
c
c
      subroutine ortho(li,lc,b,kg)
      
c    Finds the orthogonal matrix v such that the columns of b*v are orthogonal
c    the array a is replaced by a*v for levels li - lc. Array b is replaced
c    by b*v and is then ready fo entry to sprop at level lc. This is intended
c    to diminish the onset of degeneracy caused by rapid exponential growth
c    in the mantle for modes with deeply turning s and shallowly turning p.

      include'mineos.inc'
      implicit real*8(a-h,o-z)
      common/arem/a(6,3,ndata)
      dimension b(6,1),as(6,3)
      
      i1=min0(lc,li)
      i2=max0(lc,li)
      nc=kg+2
      nr=2*nc
      call svd(b,nr,nc)
      do 25 i=i1,i2
      do 20 j=1,nc
      do 20 k=1,nr
      as(k,j)=0.d0
      do 20 l=1,nc
   20 as(k,j)=as(k,j)+a(k,l,i)*b(l,j)
      do 25 j=1,nc
      do 25 k=1,nr
   25 a(k,j,i)=as(k,j)
      do 35 j=1,nc
      do 35 k=1,nr
   35 b(k,j)=a(k,j,lc)
   
100	continue

      return
      end
c
c
      subroutine svd(a,mrow,ncol)
      
c   Section I Chapter 10 Wilkenson and Reinsch (1971 ,Springer).
c   The matrix a is overwritten with v(ncol,ncol), the right side orthogonal
c   matrix in the svd decomposition. For use only in EOS subs as , to reduce
c   branching points, I have used the fact that ncol is lt mrow.

      implicit real*8(a-h,o-z)
      dimension a(6,1),e(10),q(10)
      
cxxxxxxxxxxxxxxxx eps was 1.5d-14, tol was 1.d-293.
c      eps=1.5d-14
c      tol=1.d-293
      eps=1.5d-18
      tol=1.d-293
      g=0.d0
      x=0.d0
      do 60 i=1,ncol
      l=i+1
      e(i)=g
      s=0.d0
      do 10 j=i,mrow
   10 s=s+a(j,i)*a(j,i)
      if(s.gt.tol) go to 15
      q(i)=0.d0
      if(l.gt.ncol) goto 60
      go to 30
   15 q(i)=dsign(dsqrt(s),-a(i,i))
      h=a(i,i)*q(i)-s
      a(i,i)=a(i,i)-q(i)
      if(l.gt.ncol) go to 60
      do 25 j=l,ncol
      s=0.d0
      do 20 k=i,mrow
   20 s=s+a(k,i)*a(k,j)
      f=s/h
      do 25 k=i,mrow
   25 a(k,j)=a(k,j)+f*a(k,i)
   30 s=0.d0
      do 35 j=l,ncol
   35 s=s+a(i,j)*a(i,j)
      if(s.ge.tol)go to 40
      g=0.d0
      go to 60
   40 g=dsign(dsqrt(s),-a(i,l))
      h=a(i,l)*g-s
      a(i,l)=a(i,l)-g
      do 45 j=l,ncol
   45 e(j)=a(i,j)/h
      do 55 j=l,mrow
      s=0.d0
      do 50 k=l,ncol
   50 s=s+a(j,k)*a(i,k)
      do 55 k=l,ncol
   55 a(j,k)=a(j,k)+s*e(k)
   60 x=dmax1(dabs(q(i))+dabs(e(i)),x)
      goto 100
   75 if(g.eq.0.d0)go to 91
      h=a(i,l)*g
      do 80 j=l,ncol
   80 a(j,i)=a(i,j)/h
      do 90 j=l,ncol
      s=0.d0
      do 85 k=l,ncol
   85 s=s+a(i,k)*a(k,j)
      do 90 k=l,ncol
   90 a(k,j)=a(k,j)+s*a(k,i)
   91 do 95 j=l,ncol
      a(i,j)=0.d0
   95 a(j,i)=0.d0
  100 a(i,i)=1.d0
      g=e(i)
      l=i
      i=i-1
      if(i.ge.1)go to 75
      ep=eps*x
      k=ncol
  105 l=k
  110 if(dabs(e(l)).le.ep)go to 125
      if(dabs(q(l-1)).le.ep) go to 115
      l=l-1
      if(l.ge.1)go to 110
  115 c=0.d0
      s=1.d0
      do 120 i=l,k
      f=s*e(i)
      e(i)=c*e(i)
      if(dabs(f).le.ep)go to 125
      g=q(i)
      h=dsqrt(f*f+g*g)
      c=g/h
      s=-f/h
  120 q(i)=h
  125 z=q(k)
      if(l.eq.k)go to 145
      x=q(l)
      y=q(k-1)
      g=e(k-1)
      h=e(k)
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.d0*h*y)
      g=dsqrt(f*f+1.d0)
      f=((x-z)*(x+z)+h*(y/(f+dsign(g,f))-h))/x
      c=1.d0
      s=1.d0
      lp1=l+1
      do 140 i=lp1,k
      g=e(i)
      y=q(i)
      h=s*g
      g=c*g
      z=dsqrt(f*f+h*h)
      im1=i-1
      e(im1)=z
      c=f/z
      s=h/z
      f=s*g+c*x
      g=c*g-s*x
      h=s*y
      y=c*y
      do 130 j=1,ncol
      x=a(j,im1)
      z=a(j,i)
      a(j,im1)=c*x+s*z
  130 a(j,i)=c*z-s*x
      z=dsqrt(f*f+h*h)
      q(im1)=z
      c=f/z
      s=h/z
      f=s*y+c*g
  140 x=c*y-s*g
      e(l)=0.d0
      e(k)=f
      q(k)=x
      go to 105
  145 if(z.ge.0.d0)go to 155
      q(k)=-z
      do 150 j=1,ncol
  150 a(j,k)=-a(j,k)
  155 k=k-1
      if(k.ge.1)go to 105
      
      return
      end
