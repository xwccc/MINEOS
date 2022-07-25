        program mineos_strip

c   This program takes the merged MINEOS binary output file for normal-mode eigenfunctions,
c   strips everything below a given depth rstrip and writes the eigenfunctions U, dU/dr,
c   V and dV/dr for spheroidal modes and W and dW/dr for toroidal modes from rstrip to
c   the free surface into kstrip binary files. Each of these files contains the
c   eigenfunctions of all the modes at one radial node. The radius of the radial
c   node (in m) is used as part of the file name. There is also a single binary index
c   file which contains the info on the reference model and the modal properties.
c
c   Version 1.0 -- January 2006. (C) Li Zhao, IESAS, Taipei.

        parameter(nfmax=200,nfmineos=20000)
        character*120 funmineos(nfmineos),fotpath,fileot,fileoth,modelid,ftmp,probfile,adum*7
        double precision rn,pi,bigg,rhobar,wn,phi,cv,dl,omega,alpha,qmod
        double precision ur,dur,vr,dvr,xr,yr,xl,tol,errmax
        double precision rstrip,rsprt,rcmb,ricb,rk1
 	double precision, dimension (:), allocatable :: radius,dn,pv,sv,ph,sh,qk,qm,eta,eig,r
	real, dimension (:,:), allocatable :: tdum
 	integer, dimension (:), allocatable :: kistrip
        integer::rlen
        data bigg,rhobar,rstrip,tol/6.6723d-11,5515.0d0,800.d3,1.d-3/

        pi=4.d0*datan(1.d0)
        wn=dsqrt(pi*bigg*rhobar)
	if3=0

c   Obtain path and file names of input and output files.

        write(*,'(/a)') '    Enter the 1st MINEOS binary output filename (.fun)'
        read(*,'(a)') funmineos(1)
	
c   Min and max of angular degrees. Note: the angular degree must be continuous 
c   from lmin to lmax.

        write(*,'(/a)') '    Enter the min and max angular degrees:'
        read(*,*) lmin,lmax
        write(*,'(/a)') '    Enter the min and max frequency (mHz):'
        read(*,*) fmin,fmax
        write(*,'(/a)') '    Enter the path for the stripped files'
        read(*,'(a)') fotpath
        write(*,'(/a)') '    Enter the reference model ID (PREM)'
        read(*,'(a)') modelid
        write(*,'(/a)') '    Maximum depth of extraction in km (0. for CMB):'
        read(*,*) rstrip
        write(*,'(/a)') '    Sampling radius interval in km:'
        read(*,*) rsprt
        rsprt=rsprt*1.d3
	write(*,'(/a)') '    Tolerance for Rayleigh quotient (1.d-2):'
	read(*,*) errmax
	
c   Generate the MINEOS binary file names. 

	nmineos=lmax-lmin+1
	ftmp=funmineos(1)
	lf=index(ftmp,' ')-1
	write(adum,'(i7.7)') lmin
	ftmp(lf-10:lf-4)=adum(1:7)
	if(ftmp(1:lf).ne.funmineos(1)(1:lf)) then
	  write(*,'(a,a,a,i5,a)') '    1st MINEOS binary file ',funmineos(1)(1:lf),
     1	    ' and lmin ',lmin,' mismatch. Stop.'
     	  stop
	endif
	do l=lmin+1,lmax
	  ftmp=funmineos(1)
	  lf=index(ftmp,' ')-1
	  write(adum,'(i7.7)') l
	  ftmp(lf-10:lf-4)=adum
	  funmineos(l-lmin+1)(1:lf)=ftmp(1:lf)
	end do

c   Open and read the header of the MINEOS binary file.

	lf0=index(funmineos(1),' ')-1
        open(1,file=funmineos(1)(1:lf0),status='old',form='unformatted',access='sequential')
        read(1) jcom
        read(1) nr,nic,noc
        close(1)

c   Number of saved eigenfunction values according to mode types.
	
        if(jcom.eq.3) then
          mrec=6*nr+7
        elseif(jcom.eq.1) then
          mrec=2*nr+7
        else
          mrec=2*(nr-noc)+7
        endif
	
c   Declaring arrays.
	
	allocate(tdum(nr,9))
	allocate(eig(mrec),r(nr),kistrip(nr))
	allocate(radius(nr),dn(nr),pv(nr),sv(nr),ph(nr),sh(nr),qk(nr),qm(nr),eta(nr))

c   Now read the model parameters from the MINEOS binary file.
	
        open(1,file=funmineos(1)(1:lf0),status='old',form='unformatted',access='sequential')
        read(1) jcom,wmin1,wmax1,llmin,llmax,wgrav
        read(1) nr,nic,noc,ifanis,tref,((tdum(ii,jj),ii=1,nr),jj=1,9)
        close(1)
        do i=1,nr
          radius(i)=dble(tdum(i,1))
          dn(i)=dble(tdum(i,2))
          pv(i)=dble(tdum(i,3))
          sv(i)=dble(tdum(i,4))
          qk(i)=dble(tdum(i,5))
          qm(i)=dble(tdum(i,6))
          ph(i)=dble(tdum(i,7))
          sh(i)=dble(tdum(i,8))
          eta(i)=dble(tdum(i,9))
        end do
	
        write(*,'(/a,i1,a,i6,a,i6)') '    Mode type: ',jcom,'    lmin = ',lmin,'  lmax = ', lmax
        write(*,'(a,f7.2,a,f7.2)') '    fmin = ', fmin,'  fmax = ',fmax

        rn=radius(nr)
        if(jcom.eq.2) then
          rcmb=radius(1)
          ricb=0.d0
        else
          rcmb=radius(noc)
          ricb=radius(nic)
        endif

c   Find out the node of the ocean bottom. Surface if the model is continental.

        do i=nr,1,-1
          nocb=i
          if(sv(i).gt.0.d0) goto 1
        end do

1       continue

c   Set the radius where the stripping starts. This bottom radius can be 800-km depth 
c   if only synthetics are needed. It is set to rcmb for toroidal modes and 0 for 
c   spheroidal and radial modes when kernels are needed. 

	if(rstrip.eq.0.d0) rstrip=rn-rcmb

c   Strip all the knots deeper than rstrip. kstrip is the number of saved radial nodes.
c   Between rstrip and 6300-km depth, the samples are decided according to the radial 
c   sampling interval rsprt. From 6300-km depth to the surface, every node in the 
c   original MINEOS binary file is saved. 

c   Find the index ks2 in radius near 6300-km depth deviding two regimes of 
c   different sampling intervals.

        do i=1,nr
          if(radius(i).lt.6300.d3) then
          else
            ks2=i
	    goto 5
          endif
        end do

c   Sample from rstrip to radius(ks2).

5       kstrip=0
	do i=nr,1,-1
	  ist1=i
	  if((rn-radius(i)).ge.rstrip) goto 51
	end do
51	continue
        do i=ist1,nr
          if(i.ge.ks2) goto 6
          if((rn-radius(i)).le.rstrip) then
            if((i.eq.ist1).or.(i.eq.nr)) then
              kstrip=kstrip+1
              r(kstrip)=radius(i)/rn
              kistrip(kstrip)=i
            elseif(dabs(radius(i)-radius(i+1)).lt.tol) then
	      kstrip=kstrip+1
	      r(kstrip)=radius(i)/rn
	      kistrip(kstrip)=i
	      kstrip=kstrip+1
	      r(kstrip)=radius(i+1)/rn
	      kistrip(kstrip)=i+1
	    else
	      rk1=0.d0
	      if(kstrip.gt.0) rk1=r(kstrip)*rn
              if(dabs(radius(i)-rk1).lt.rsprt) then
              else
                kstrip=kstrip+1
                r(kstrip)=radius(i)/rn
                kistrip(kstrip)=i
              endif
            endif
          endif
        end do

c   Sample from radius(ks2) to radius(nocb).

6       continue
        do i=ks2,nr
          kstrip=kstrip+1
          r(kstrip)=radius(i)/rn
          kistrip(kstrip)=i
        end do
	
        write(*,'(/a,i6)') '    Total radial nodes:    ',nr
        write(*,'(a,i6)') '    Level of ocean bottom: ',nocb
        write(*,'(a,i9,a,i6,a,f9.4,a)') '    Nodes after strip = ',kstrip,
     1    '.  Starting from ',kistrip(1),' at ',r(1)*rn/1.d3,' km.'

c   Since there is a limit on the maximum number of concurrently opened files by a
c   single process, the writing is divided into nloop pieces, each of which involve
c   nfmax files, hence nfmax radial nodes. The 1st piece corresponds to the first
c   nfmax radial nodes starting from rstrip, the 2nd piece corresponds to the next
c   nfmax radial nodes, and so on, until the nloop-th piece which corresponds to
c   the top kstrip-(nloop-1)*nfmax radial nodes. The nloop pieces can be run either
c   in parallel or consecutively depending on iloop0 entered below.
c   If iloop0=0, then a single run of this program will write all the piece by piece.
c   If iloop0=i, where 1<=i<=nloop, then the current run of this program only writes
c                the files corresponding to the i-th piece.
c   Thus, a maximum of nloop copies of this program can be run simultaneously.

	if(kstrip.le.nfmax) then
	  nloop=1
	  iloop0=0
	elseif(mod(kstrip,nfmax).eq.0) then
	  nloop=int(kstrip/nfmax)
	else
	  nloop=int(kstrip/nfmax)+1
	endif
	
        write(*,'(/a,i5,a,i2,a)') '    A total of ',kstrip,
     1    ' files need to be processed in ',nloop,' loops.'
     
     	if(nloop.gt.1) then
          write(*,'(/a)') '    Enter loop number (0 for all):'
          read(*,*) iloop0
          if((iloop0.lt.0).or.(iloop0.gt.nloop)) then
            write(*,'(a,i2)') '    Stop. iloop must be from 0 to ',nloop
            goto 1000
          endif
	endif

	kf=index(fotpath,' ')-1

c   Create output file directory.

	ftmp=''
	ftmp='mkdir -p '//fotpath(1:kf)
	kff=kf+9
	if(ftmp(kff:kff).ne.'/') then
	  kff=kff+1
	  ftmp(kff:kff)='/'
	endif
	iop=10
300	iop=iop+1
	if(iop.gt.kff) goto 303
	if(ftmp(iop:iop).ne.'/') goto 300
	call system(ftmp(1:iop))
	goto 300
303	continue

	mkf=index(modelid,' ')-1
        fileot=fotpath(1:kf)//'/'//modelid(1:mkf)//'_'
	kf=index(fileot,' ')-1
        fileoth(1:kf-1)=fileot(1:kf-1)
        fileoth(kf:kf+3)='.sph'
        fileot(kf+1:kf+12)='00000000.sph'
        rlen=5*2

        if(jcom.eq.2) then
          fileoth(kf+1:kf+3)='tor'
          fileot(kf+10:kf+12)='tor'
          rlen=2*2
	elseif(jcom.eq.1) then
          fileoth(kf+1:kf+3)='rad'
          fileot(kf+10:kf+12)='rad'
          rlen=2*2
        endif

c   Start looping over the nloop pieces.

	mr=0
        iloop=iloop0
100     continue
        if(iloop0.eq.0) iloop=iloop+1

        kr1=nfmax*(iloop-1)+1
        kr2=nfmax*iloop
        if(kr2.ge.kstrip) kr2=kstrip
	mr=mr+(kr2-kr1)+1

        write(*,'(/a,i2,a,f8.3,a,f8.3,a/16x,a,i5,10x,a,i5/)') '    Loop # ',iloop,'.  Rmin: ',
     1    r(kr1)*rn/1.d3,'(km). Rmax: ',r(kr2)*rn/1.d3,'(km).','kr1: ',kr1,'kr2: ',kr2

c   Generate the names of the stripped files according to the radius of the radial nodes.
c   Then open the files for writing.

        do i=kr1,kr2
          ir=int(r(i)*rn)
	  write(adum,'(i7.7)') ir
	  fileot(kf+1:kf+7)=adum(1:7)
          fileot(kf+8:kf+8)='0'
          if(i.eq.1) then
          else
            if(dabs(rn*(r(i)-r(i-1))).lt.tol) then
              fileot(kf+8:kf+8)='1'
            endif
          endif
          if((mod(i,50).eq.0).or.(i.eq.kr1).or.(i.eq.kr2))
     1      write(*,'(a,a)') '    Opening stripped file: ',fileot(1:kf+12)
          open(100+i,file=fileot(1:kf+12),status='unknown',form='unformatted',
     1	    access='direct',recl=rlen)
        end do
	
        ntmode=0
	nmin=100000
	nmax=-1
	ntot=0
	
c   Write the header of the index file (only in 1st loop). Open the problem file.

	if(iloop.eq.1) then
          write(*,'(/a,a)') '    Opening index file: ',fileoth(1:kf+3)
          open(2,file=fileoth(1:kf+3)//'.tmp',status='unknown',form='unformatted',access='sequential')
          write(2) jcom,fmin,fmax,lmin,lmax,wgrav
          write(2) rn,wn,wn,ifanis,tref
          write(2) kstrip,(r(i1),i1=1,kstrip)
          write(2) nr,(radius(i1),i1=1,nr),(dn(i1),i1=1,nr),(pv(i1),i1=1,nr),
     1      (sv(i1),i1=1,nr),(ph(i1),i1=1,nr),(sh(i1),i1=1,nr),(eta(i1),i1=1,nr)
          write(2) nr,(radius(i1),i1=1,nr),(qk(i1),i1=1,nr),(qm(i1),i1=1,nr)
	  probfile=fileoth(1:kf+3)//'.prob'
          write(*,'(/a,a)') '    Opening problem mode file: ',probfile(1:kf+8)
	  open(3,file=probfile(1:kf+8),status='unknown')
	  if3=1
	endif
	
c   Loop over all the MINEOS binary files.

	do j1=1,nmineos
	
          open(1,file=funmineos(j1),status='old',form='unformatted',access='sequential')
          read(1) jcom,wmin1,wmax1,llmin,llmax,wgrav
          read(1) nnr,nnic,nnoc,ifanis,tref,((tdum(ii,jj),ii=1,nnr),jj=1,9)

c   Start reading the MINEOS binary file mode by mode.

20        read(1,end=25) eig
	  ntot=ntot+1
	  if(dabs(eig(6)).gt.errmax) then
	    if(iloop.eq.1) write(3,'(a,2(1x,i5),5(1x,g16.7))') 'Problem mode:    ',
     1	      int(eig(1)),int(eig(2)),eig(3),eig(4),eig(5),eig(6),eig(7)
	    goto 20
	  endif
          ntmode=ntmode+1
	  nn=int(eig(1))
	  if(nmin.gt.nn) nmin=nn
	  if(nmax.lt.nn) nmax=nn
          ll=int(eig(2))
	  dl=0.d0
          if(ll.gt.0) dl=1.d0/dsqrt(dble(ll)*(dble(ll)+1.d0))

c   Write the record for the current mode in the index file (only in 1st loop). 

          if(iloop.eq.1) then
            cv=(rn/1.d3)*eig(3)*dl
            alpha=0.5d0*eig(3)/eig(4)
            if((jcom.eq.3).or.(jcom.eq.1)) then
              phi=eig(7+4*nr+nocb)
            else
              phi=0.d0
            endif
            write(2) nn,ll,eig(3),eig(4),alpha,phi,cv,real(eig(5)),real(eig(6))
          endif

          omega=eig(3)/wn
          xl=dsqrt(eig(2)*(eig(2)+1.d0))

c   Renormalizing the eigenfunctions according to Dahlen & Tromp (1998).

          do i1=8,mrec
            eig(i1)=eig(i1)*omega
          end do

c   Note: V, dV/dr, W and dW/dr are renormalized.

	  do j=kr1,kr2
            i=kistrip(j)
            if(jcom.eq.3) then
              ur=eig(7+i)
              dur=eig(7+nr+i)
              vr=eig(7+2*nr+i)
              dvr=eig(7+3*nr+i)
              xr=2.d0*ur-xl*vr
              yr=r(j)*dvr-vr+xl*ur
              write(100+j,rec=ntmode) ur,dur,vr,xr,yr
            elseif(jcom.eq.2) then
              ur=eig(7+i)
              dur=r(j)*eig(7+nr+i)-ur
              write(100+j,rec=ntmode) ur,dur
            elseif(jcom.eq.1) then
              ur=eig(7+i)
              dur=eig(7+nr+i)
              write(100+j,rec=ntmode) ur,dur
           endif
	  end do

          if(mod(ntmode,20000).eq.0) then
            write(*,'(a,i2,a,i7,1x,i4,1x,i5,4(1x,e11.4))') '    Loop # ',iloop,' Mode # ',ntmode,
     1        int(eig(1)),ll,1.d3*eig(3)/(2.d0*pi),eig(4),eig(5),eig(6)
          endif

          goto 20
25        close(1)

	end do
	
        if(iloop.eq.1) then
          close(3)
          open(3,file=fileoth(1:kf+3),status='unknown',form='unformatted',access='sequential')
	  rewind(2)
          read(2) jcom,fmin,fmax,lmin,lmax,wgrav
	  write(3) jcom,fmin,fmax,ntmode,nmin,nmax,lmin,lmax,wgrav
          read(2) rn,wn,wn,ifanis,tref
          write(3) rn,wn,wn,ifanis,tref
          read(2) kstrip,(r(i1),i1=1,kstrip)
          write(3) kstrip,(r(i1),i1=1,kstrip)
          read(2) nr,(radius(i1),i1=1,nr),(dn(i1),i1=1,nr),(pv(i1),i1=1,nr),
     1      (sv(i1),i1=1,nr),(ph(i1),i1=1,nr),(sh(i1),i1=1,nr),(eta(i1),i1=1,nr)
          write(3) nr,(radius(i1),i1=1,nr),(dn(i1),i1=1,nr),(pv(i1),i1=1,nr),
     1      (sv(i1),i1=1,nr),(ph(i1),i1=1,nr),(sh(i1),i1=1,nr),(eta(i1),i1=1,nr)
          read(2) nr,(radius(i1),i1=1,nr),(qk(i1),i1=1,nr),(qm(i1),i1=1,nr)
          write(3) nr,(radius(i1),i1=1,nr),(qk(i1),i1=1,nr),(qm(i1),i1=1,nr)
	  do i=1,ntmode
            read(2) nn,ll,omega,qmod,alpha,phi,cv,x1,x2
            write(3) nn,ll,omega,qmod,alpha,phi,cv,x1,x2
	  end do
	  close(2,status='delete')
          close(3)
	endif

        do i=kr1,kr2
          close(100+i)
        end do

        if((iloop0.eq.0).and.(kr2.lt.kstrip)) goto 100

	write(*,'(/a,i9,a,i5,a,i5,a,e8.1,a)') '    Done.  ',ntmode,' modes at ',
     1	  mr,' radial nodes are saved.',ntot-ntmode,' bad modes (ERR:',errmax,').'
     	if(if3.eq.1) write(*,'(/a,a,a/)') '    See ',probfile(1:kf+8),' for the list of bad modes.'	
     
1000    continue

        stop
        end
