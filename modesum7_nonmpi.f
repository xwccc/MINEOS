        module modesummodule
	
	double precision pi,grav,rnorm,dnorm	
	parameter(nfmax=50,nstmax=10000)
        parameter(grav=6.6723d-11,pi=3.14159265358979d0)
        parameter(rnorm=6371.d3,dnorm=5515.d0)
        character(len=120) paths(nfmax),files(nfmax),modepath,modefile,opath
	character(len=120) eventid,stat,amode*3,inname(nstmax),innamex
        double precision, dimension (:,:), allocatable :: uksr,duksr,vksr
	double precision, dimension (:,:), allocatable :: xksr,yksr,ak
	double precision, dimension (:), allocatable :: fx,atten
        double precision, dimension (:,:,:), allocatable :: glf0
	double precision, dimension (:), allocatable :: rknots,rmod,sv
        double precision bpf(4),mg(3,3),mp(3,3),xm(6),xm0,t1,t2,dt,tmshift,hdur
        double precision evtlat,evtlon,depths0,stlat,stlon,stelv,rs,rr
	real, dimension (:), allocatable :: weight
	integer, dimension (:,:), allocatable :: icount
	integer, dimension (:), allocatable :: lk,ntbl
	integer nmodes,nm,nmax1,lmax1,kf,nrept,nfiles,nfilet,nfiler
	integer knots,knrad,nrmax
	integer nst
	parameter(nchmax=nstmax*120)
	
	end module modesummodule
c
c
c   This program computes the complete synthetic seismograms in SNREI 
c   Earth models by the normal-mode summation method. It uses the 
c   normal-mode eigensolutions calculated by the MINEOS code and 
c   post-processed by mineos_strip.f which creates a summary file 
c   containing the eigenfrequencies and other modal properties of all 
c   the modes and a set of files each containing the eigenfunctions 
c   of all the modes at a specifc radius.
c
c   An input file modesum.in provides the event and station information 
c   as well as the component and time and frequency ranges of the 
c   synthetic seismogram.
c
c   The code here is written based on unpublished research notes 
c   on the summation of normal modes in SNREI models. 
c
c   (C) Li Zhao, 2006; 2007; 2009; 2010.
c       zhaol@earth.sinica.edu.tw
c       Institute of Earth Sciences (IES)
c       Academia Sinica (AS)
c       128 Academia Road Sec.2
c       Nankang, Taipei 115
c       Taiwan.
c
c   Version 1.0 - Nov. 2006, IESAS.
c   Version 2.0 - Nov. 2006, IESAS. Enable multiple mode table files.
c                 Compute Z and R components simultaneously.
c   Version 2.1 - Jan. 2007, IESAS. Merge spheroidal and toroidal sums 
c                 so that spheroidal-mode contribution to T component 
c                 and toroidal-mode to R component are also calculated. 
c                 All three component are now computed in the same run.
c   Version 3.0 - Aug. 2007, IESAS. Correct mistake discovered by 
c                 Hsin-Ying Yang in transforming the moment tensors 
c                 from GCS to PCS.
c                 Add function to input source fault-plane solutions.
c   Version 4.0 - March 2007, IESAS. Adapt to new stripped mode table 
c                 format which incorporates spheroidal eigenfunctions  
c                 in the outer and inner cores.
c   Version 5.0 - February 2009, IESAS. Allow sources and receivers 
c                 in the ocean. Adopt a more complete time function 
c                 ftime, although its effect is extremely small.
c   Version 6.0 - October 2010, IESAS. Added radial modes in summation.
c                 Add a choice of applying sine-square frequency tapers 
c                 in normal-mode summation to reduce Gibbs oscillations.
c   Version 7.0 - November 2010, IESAS. MPI-enabled. Also added source-time 
c                 function convolution and time shift following the 
c                 Global CMT convention (Dziewonski & Woodhouse 1983).
c
c   In case you use this code, please cite the following article:
c
c       Yang, H.-Y., Zhao, L. and Hung, S.-H., Synthetic seismograms by normal-mode
c           summation: a new derivation and numerical examples, Geophys. J. Int., 
c           183, 1613-1632.
c
        program modesum7_mpi

	use modesummodule
        character(len=120) comp*1
        character(len=40) datex,timex,datex1,datex2,timex1,timex2        
        double precision, dimension (:,:), allocatable :: u
        double precision, dimension (:), allocatable :: fx2p,fx2m
        double precision, dimension (:), allocatable :: t,capl,capl2,utmp
        double precision phis,phir,thetas,thetar,depths,depthr
	double precision phirs,srs2,crs2,srs,crs,s0,s1r,s1i,s2r,s2i
        double precision vnorm,tnorm,fnorm,anorm,mnorm
        double precision xl,att,ftime,xmm,factor
	double precision slat,slon,rlat,rlon,delta,azim,bazim
	double precision fxf,wfilter,bpf12,bpf34,rdum
	double precision alpha,alphatau,expon,source
	
c   Non-dimensionalization parameters for time, velocity, frequency, 
c   acceleration as well as moment tensor. Constants dnorm and rnorm 
c   are for density and radius. The two basic quantities are rnorm and
c   dnorm usually taking values of the radius of the Earth (6371000 m)
c   and 5515 kg/m^3 as hardwired in MINEOS. If the normal modes are
c   computed with software other than MINEOS, these two values can be
c   defined so as to make most of the functions close to unity.
 
	tnorm=1.d0/dsqrt(pi*dnorm*grav)
        vnorm=rnorm/tnorm
        fnorm=1.d0/tnorm
        anorm=vnorm/tnorm
        mnorm=dnorm*vnorm**2*rnorm**3
	alpha=1.628d0
	
c   Record the date and time at the beginning of the job.

	call date_and_time(datex,timex)
        write(*,*) ' '
        write(*,'(a,a4,a,a2,a,a2,a,a2,a,a2,a,a4)') 
     1    '    Starting date and time:                     ',
     2    datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',
     3    timex(1:2),':',timex(3:4),':',timex(5:8)

c   Read input file.

        call scrinp
	
	do 1000 iloop=1,nst
	
	innamex=''	
	innamex=inname(iloop)
	
	call readinp
	
	nt0=int((t2-t1)/dt)
	
	alphatau=alpha/hdur
	nj=int(1.5d0*hdur/dt)+1
     
        write(*,'(a,2x,f6.2,2x,f7.2,2x,f7.1)') 
     1    'Source coordinates (lat, long, depth): ',
     1    evtlat,evtlon,depths0
        write(*,'(a,2x,f6.2,2x,f7.2,2x,f8.2)') 
     1	  'Station coordinates (lat, long, depth): ',
     2    stlat,stlon,stelv
        
	slat=evtlat
	slon=evtlon
	rlat=stlat
	rlon=stlon
	
c   Compute the geocentric source-receiver distance. Note: azim is the 
c   angle at the source from North to the source->receiver line measured  
c   clockwise. bazim is the angle at the receiver from North to the 
c   receiver->source line measured clockwise. 
     	
	call azimth(0,slat,slon,rlat,rlon,delta,azim,bazim)
               
	phirs=pi-azim*pi/180.d0
	srs2=dsin(2.d0*phirs)
	crs2=dcos(2.d0*phirs)
	srs=dsin(phirs)
	crs=dcos(phirs)
	
        write(*,'(a,4(2x,f8.4))') 'Distance and azimuths (deg): ',
     1	  delta,azim,bazim
     
c   Normalize the moment tensor elements.
     
	do i=1,3
          do j=1,3
            mp(i,j)=mg(i,j)/mnorm
          end do
        end do

c   Set moment tensor elements of extremely small absolute value to zero.
	
	xmm=0.d0
	do i=1,3
	  do j=1,3
	    xmm=xmm+mp(i,j)**2
	  end do
	end do
	xmm=dsqrt(xmm)
	do i=1,3
	  do j=1,3
	    if(dabs(mp(i,j)/xmm).lt.1.d-13) mp(i,j)=0.d0
	  end do
	end do
	
	write(*,'(a,6(1x,e10.3))') 'Normalized moment tensor: ',
     1	  mp(1,1),mp(2,2),mp(3,3),mp(1,2),mp(1,3),mp(2,3)  
           
c   Find out lmax by reading the summary mode files.

	lmin1=1000000
	lmax1=-1
	nm=0
	nmax1=-1
	nrmax=-1
	do i=1,nfiles+nfilet+nfiler
	  amode='sph'
	  if(i.gt.nfiles) amode='tor'
	  if(i.gt.(nfiles+nfilet)) amode='rad'
	  modepath=''
	  modepath=paths(i)
	  modefile=''
	  modefile=files(i)
	  k=index(modepath,' ')-1
	  k1=index(modefile,' ')-1
	  if(modepath(k:k).ne.'/') then
	    k=k+1
	    modepath(k:k)='/'
	  endif
	  write(*,'(/a,a)') '    Examining: ',modepath(1:k)//modefile(1:k1)
     1	    //'.'//amode
          open(10,file=modepath(1:k)//modefile(1:k1)//'.'//amode,status='old',
     1      form='unformatted',access='sequential')
          read(10) jcom,wwmin,wwmax,nmodes,nnmin,nnmax,llmin,llmax,wgrav
	  read(10) rdum
          read(10) knots
          read(10) knrad
	  close(10)
	  if(lmin1.gt.llmin) lmin1=llmin
	  if(lmax1.lt.llmax) lmax1=llmax
	  if(i.le.nfiles) nm=nm+nmodes
	  if(nmax1.lt.nnmax) nmax1=nnmax
	  if(nrmax.lt.knots) nrmax=knots
	  if(nrmax.lt.knrad) nrmax=knrad
          write(*,'(a,2(1x,i6,1x,i6,1x,f7.2))') 
     1	    '    min-max l in mode table:  ',llmin,llmax,wwmax
     	  if(dabs(bpf(4)).lt.wwmax) bpf(4)=wwmax
	end do
	
	write(*,'(/a,i6,1x,i6)') 
     1	  '    Total range of l:                   ',lmin1,lmax1
	write(*,'(/a,i6,1x,i6)') 
     1	  '    Maximum n and size of modal arrays: ',nmax1,nm
	
	if(iloop.eq.1) then
	  allocate(rknots(nrmax))
	  allocate(rmod(nrmax),sv(nrmax))
	  allocate(fx(nm),atten(nm),ak(3,nm),fx2p(nm),fx2m(nm))
	  allocate(ntbl(nm))
	  allocate(uksr(nm,2),duksr(nm,2),vksr(nm,2))
	  allocate(xksr(nm,2),yksr(nm,2))
	  allocate(icount(0:nmax1,0:lmax1))
	  allocate(lk(nm),weight(nm))
	  allocate(capl(0:lmax1),capl2(lmax1))
	  allocate(glf0(lmax1,5,5))
	endif
	allocate(t(0:nt0),u(3,0:nt0),utmp(0:nt0))
	
     	bpf12=(pi/2.d0)/(bpf(2)-bpf(1))
	bpf34=(pi/2.d0)/(bpf(4)-bpf(3))
	if(bpf(2).gt.0.d0) then
	  write(*,'(/a,4(1x,f7.2))') 
     1	    '    Frequencies of sine-square tapers applied (mHz): ',
     2	    (bpf(ibp),ibp=1,4)
	endif
	
        do l=0,lmax1
          capl(l)=dsqrt((2.d0*dble(l)+1.d0)/(4.d0*pi))
          if(l.ge.1) capl2(l)=0.25d0*dsqrt((dble(l)-1.d0)*(dble(l)+2.d0))
        end do
          
c   Calculating the generalized Legendre functions for theta=delta from 
c   l=1 to l=lmax1. 

        write(*,'(/a,i6/)') '    Calculating the GLF up to l=',lmax1
        call glegendre(delta)

c   Loop over mode table files and conduct mode summation.

	write(*,'(a,i3)') '    Starting mode summation. Number of mode files: ',
     1	  nfiles+nfilet+nfiler
	
	do it=0,nt0
	  t(it)=t1+dble(it)*dt
	end do
	
	write(*,'(a,2(2x,f8.3))') 
     1	  '    Summation time window (begin, end): ',t(0),t(nt0)

	u=0.d0
	
	nmodest=0
	nsumtot=0
	nrept=0
	do 50 kf=1,nfiles+nfilet+nfiler
	
	if((kf.eq.1).or.(kf.eq.(nfiles+1))) icount=0

	modepath=paths(kf)
	modefile=files(kf)
	
	amode='sph'
	if(kf.gt.nfiles) amode='tor'
	if(kf.gt.(nfiles+nfilet)) amode='rad'
	  
c   Read in the nomal mode eigenfrequencies and group velocities for 
c   the background earth model.

	nmodes=0
        call modect

	if(nmodes.eq.0) goto 50
          
c   Compute the time-independent part of the summation coefficients for
c   individual modes.

	write(*,'(/a)') '    Computing time-independent mode-sum coefficients.'

	do k=1,nmodes
	  if(weight(k).ne.0.) then
	    l=lk(k)
	    
	    fx2p(k)=fnorm**2/(fx(k)**2+atten(k)**2)
	    fx2m(k)=(fx(k)**2-atten(k)**2)/fnorm**2

c   Apply a sine-square filter to normal-mode summation. When bpf(2)<=0, 
c   no filter is applied.
	    
	    wfilter=1.d0
	    if(bpf(2).gt.0.d0) then
	      fxf=500.d0*fx(k)/pi
	      if((fxf.le.bpf(1)).or.(fxf.ge.bpf(4))) then
	        wfilter=0.d0
	      elseif((fxf.ge.bpf(2)).and.(fxf.le.bpf(3))) then
	      elseif((fxf.gt.bpf(1)).and.(fxf.lt.bpf(2))) then
	        wfilter=(dsin(bpf12*(fxf-bpf(1))))**2
	      elseif((fxf.gt.bpf(3)).and.(fxf.lt.bpf(4))) then
	        wfilter=(dsin(bpf34*(bpf(4)-fxf)))**2
	      endif
	    endif
	    
c   Spheroidal modes. 

     	    if(kf.le.nfiles) then
	      
	      s0=mp(1,1)*duksr(k,1)+0.5d0*(mp(2,2)+mp(3,3))*xksr(k,1)/rs
	      s1r=0.5d0*mp(1,2)*yksr(k,1)/rs
	      s1i=0.5d0*mp(1,3)*yksr(k,1)/rs
	      s2r=capl2(l)*(mp(2,2)-mp(3,3))*vksr(k,1)/rs
	      s2i=2.d0*capl2(l)*mp(2,3)*vksr(k,1)/rs
	      
	      ak(1,k)=wfilter*capl(l)*uksr(k,2)*(2.d0*(s2r*crs2+s2i*srs2)
     1	        *glf0(l,1,3)+2.d0*(s1r*crs+s1i*srs)*glf0(l,2,3)+s0*glf0(l,3,3))
	      ak(2,k)=wfilter*capl(l)*vksr(k,2)*((s2r*crs2+s2i*srs2)
     1	        *(glf0(l,1,2)+glf0(l,5,2))+(s1r*crs+s1i*srs)
     2	        *(glf0(l,2,2)-glf0(l,4,2))+s0*glf0(l,3,2))
	      ak(3,k)=-wfilter*capl(l)*vksr(k,2)*((s2r*srs2-s2i*crs2)
     1	        *(glf0(l,1,2)-glf0(l,5,2))+(s1r*srs-s1i*crs)
     2	        *(glf0(l,2,2)+glf0(l,4,2)))
     
	    elseif((kf.gt.nfiles).and.(kf.le.(nfiles+nfilet))) then
	    
c   Toroidal modes. 

	      s0=0.d0
	      s1r=-0.5d0*mp(1,3)*duksr(k,1)/rs
	      s1i=0.5d0*mp(1,2)*duksr(k,1)/rs
	      s2r=-2.d0*capl2(l)*mp(2,3)*uksr(k,1)/rs
	      s2i=capl2(l)*(mp(2,2)-mp(3,3))*uksr(k,1)/rs
	      
	      ak(1,k)=0.d0
	      ak(2,k)=-wfilter*capl(l)*uksr(k,2)*((s2r*srs2-s2i*crs2)
     1	        *(glf0(l,1,2)-glf0(l,5,2))+(s1r*srs-s1i*crs)
     2	        *(glf0(l,2,2)+glf0(l,4,2)))
	      ak(3,k)=-wfilter*capl(l)*uksr(k,2)*((s2r*crs2+s2i*srs2)
     1	        *(glf0(l,1,2)+glf0(l,5,2))+(s1r*crs+s1i*srs)
     2	        *(glf0(l,2,2)-glf0(l,4,2)))
     
c   Radial modes. 

     	    elseif(kf.gt.(nfiles+nfilet)) then
	      
	      s0=mp(1,1)*duksr(k,1)+(mp(2,2)+mp(3,3))*uksr(k,1)/rs
	      
	      ak(1,k)=wfilter*capl(l)*uksr(k,2)*s0
	      ak(2,k)=0.d0
	      ak(3,k)=0.d0
     
	    endif
	  endif
	end do

c   Loop over time samples and sum the modes.

	write(*,'(/a/)') '    Computing the waveform time series.'
	mt=nt0/15
	nt1=0
	if(t(nt1).le.0.d0) nt1=1
	do it=nt1,nt0
	  if(mod(it,mt).eq.0) write(*,'(a,i7,a,i7)') 
     1	    '    Time step: ',it, ' out of ',nt0
          do k=1,nmodes
	    if(it.eq.nt1) nmodest=nmodest+1
	    if(weight(k).ne.0.) then
	      if(it.eq.nt1) nsumtot=nsumtot+1
              att=dexp(-atten(k)*t(it))
	      ftime=fx2p(k)*(fx2m(k)*fx2p(k)*(1.d0-att*dcos(fx(k)*t(it)))-
     1	        2.d0*(fx(k)*atten(k)/fnorm**2)*fx2p(k)*att*dsin(fx(k)*t(it)))
	      do j=1,3
                u(j,it)=u(j,it)+ak(j,k)*ftime
	      end do
	    endif
          end do
	end do

50	continue

c   Convolve with source time function if hdur>0.

	if(hdur.gt.0.d0) then
	  do jc=1,3
	    utmp=0.d0
	    do it=0,nt0
	      do j=-nj,nj
	        if((it.gt.j).and.((it-j).le.nt0)) then
	          source=0.d0
	          expon=alphatau**2*(dble(j)*dt)**2
	          if(expon.le.50.d0) source=alphatau*dexp(-expon)/dsqrt(pi)
	          utmp(it)=utmp(it)+source*u(jc,it-j)*dt
	        endif
	      end do
	    end do
	    u(jc,0:nt0)=utmp(0:nt0)
	  end do
	endif
	
c   Shift the seismograms by time shift in Global CMT solution. 

	nshift=int(tmshift/dt)
	do jc=1,3
	  if(nshift.gt.0) then
	    do it=nt0-nshift,0,-1
	      jt=it+nshift
	      u(jc,jt)=u(jc,it)
	    end do
	    do it=0,nshift-1
	      u(jc,it)=0.d0
	    end do
	  elseif(nshift.lt.0) then
	    do it=-nshift,nt0
	      jt=it+nshift
	      u(jc,jt)=u(jc,it)
	    end do
	    do it=nt0+nshift+1,nt0
	      u(jc,it)=0.d0
	    end do
	  endif
	end do
	
c   Open output file and write the header. Note: T compmnent is reversed here 
c   to conform with the sign convention in SAC. 

	ko=index(opath,' ')-1
        if(opath(ko:ko).ne.'/') then
          ko=ko+1
          opath(ko:ko)='/'
        endif
	k=index(stat,' ')-1
	ke=index(eventid,' ')-1
        open(1,file=opath(1:ko)//eventid(1:ke)//'.'//stat(1:k)//'.ASC',
     1	  status='unknown')
        write(1,'(a)') eventid(1:ke)
        if(xm0.le.10.) write(1,'(f9.5,1x,f10.5,1x,f6.2,1x,f3.0,3(1x,f7.2),2(1x,f6.2))') 
     1	  evtlat,evtlon,depths0,xm0,(xm(i),i=1,3),tmshift,hdur
        if(xm0.gt.10.) write(1,'(f9.5,1x,f10.5,1x,f6.2,1x,f3.0,6(1x,f7.3),2(1x,f6.2))') 
     1	  evtlat,evtlon,depths0,xm0,(xm(i),i=1,6),tmshift,hdur
        write(1,'(a)') stat(1:k)
        write(1,'(3(1x,f10.5),2(1x,f8.4))') stlat,stlon,stelv,delta,azim
	write(1,'(i1,1x,a3)') 3,'ZRT'
        write(1,'(i6,2x,f8.3)') nt0+1,t(0)
	do it=0,nt0
          write(1,'(f9.3,3(1x,e15.8))') t(it),(u(j,it)*rnorm,j=1,2),-u(3,it)*rnorm
	end do
        close(1)
	
	write(*,'(/a,i9,a,i9,a,i9)') '    Summation done. Total mode number: ',
     1	  nmodest,'. Repeated: ',nrept,'. Summed: ',nsumtot
	write(*,'(a,a)') '    Output filename: ',
     1	  opath(1:ko)//eventid(1:ke)//'.'//stat(1:k)//'.ASC'

c   Record the date and time at the beginning of the job.

        call date_and_time(datex,timex)
        write(*,'(/a,a4,a,a2,a,a2,a,a2,a,a2,a,a4/)') 
     1    '    Finishing date and time:                    ',
     2    datex(1:4),'-',datex(5:6),'-',datex(7:8),'.  ',
     3    timex(1:2),':',timex(3:4),':',timex(5:8)

	deallocate(t,u,utmp)
	
1000	continue

        stop
        end
c
c
        subroutine scrinp
        
c   This subroutine reads in the list of input file names.

        use modesummodule
	character(len=120) listfile
	
        call getarg(1,listfile)
	open(1,file=listfile,status='old')
	nst=0
1	read(1,'(a)',end=2) inname(nst+1)
	nst=nst+1
	goto 1
2	close(1)

	return
	end
c
c	
        subroutine readinp
        
c   This subroutine reads in a few basic input parameters from the keyboard.
c   They include the event id, usually an eight-digit integer specifying 
c   event date, event location and moment tensor, station code and location, 
c   component (T, R or Z), path for the mode table file, mode table file id, 
c   and time and frequency ranges of the synthetic seismogram.

        use modesummodule
	double precision strike,dip,rake,m11,m22,m33,m12,m13,m23
	double precision mbar,xm00
	
	bpf=0.d0
	
	open(1,file=innamex,status='old')
        read(1,'(a)') eventid
	read(1,*) evtlat,evtlon,depths0
	read(1,*) xm0,(xm(i),i=1,6),tmshift,hdur
	read(1,'(a)') stat
	read(1,*) stlat,stlon,stelv
	read(1,*) t1,t2,dt
	read(1,*) nfiles
	do i=1,nfiles
	  read(1,'(a)') paths(i)
	  read(1,'(a)') files(i)
	end do
	read(1,*) nfilet
	do i=1,nfilet
	  read(1,'(a)') paths(i+nfiles)
	  read(1,'(a)') files(i+nfiles)
	end do
	read(1,*) nfiler
	do i=1,nfiler
	  read(1,'(a)') paths(i+nfiles+nfilet)
	  read(1,'(a)') files(i+nfiles+nfilet)
	end do
	read(1,'(a)') opath
	read(1,*) bpf(2),bpf(3)
	close(1)
	
	bpf=bpf*1000.d0
	bpf(1)=bpf(2)/2.d0
	
c   In Global CMT convention M0 has the unit of dyne*cm. Here 7 is subtractedfrom 
c   it to convert to N.m.

	xm00=10.d0**(xm0-7.d0)
	mg(1,1)=xm(1)
	mg(2,2)=xm(2)
	mg(3,3)=xm(3)
	mg(1,2)=xm(4)
	mg(1,3)=xm(5)
	mg(2,3)=xm(6)
	mg(2,1)=mg(1,2)
	mg(3,1)=mg(1,3)
	mg(3,2)=mg(2,3)
	if(xm0.gt.10.) goto 1
	
c   If xm0 <= 10, source mechanism is given in terms of strike/dip/rake 
c   instead of moment tensor. The fault-plane solution is converted to 
c   moment tensor here. strike=xm(1), dip=xm(2), rake=xm(3).
c   First calculate the moment tensor in XYZ system at the source in which
c   X - East; Y - North; Z - Up.

	strike=xm(1)*pi/180.d0
	dip=xm(2)*pi/180.d0
	rake=xm(3)*pi/180.d0
	m11=dsin(dip)*dcos(rake)*dsin(2.d0*strike)-dsin(2.d0*dip)*dsin(rake)
     1	  *(dcos(strike)**2)
	m22=-dsin(dip)*dcos(rake)*dsin(2.d0*strike)-dsin(2.d0*dip)*dsin(rake)
     1	  *(dsin(strike)**2)
	m33=dsin(2.d0*dip)*dsin(rake)
	m12=dsin(dip)*dcos(rake)*dcos(2.d0*strike)+0.5d0*dsin(2.d0*dip)
     1	  *dsin(rake)*dsin(2.d0*strike)
	m13=dcos(dip)*dcos(rake)*dsin(strike)
     1	  -dcos(2.d0*dip)*dsin(rake)*dcos(strike)
	m23=dcos(dip)*dcos(rake)*dcos(strike)
     1	  +dcos(2.d0*dip)*dsin(rake)*dsin(strike)
     
c   xm0 is assumed to be moment magnitude and should be converted to scalar moment.  
c   The conversion assumes the relation log(M0)=1.5Mw+16.1. This equation also 
c   assumes the units of dyne and cm, hence 7 is subtracted to convert to N and m.

     	mbar=dsqrt(m11*m11+m22*m22+m33*m33+2.d0*(m12*m12+m13*m13+m23*m23))/dsqrt(2.d0)
	xm00=10.d0**(1.5d0*xm0+16.1d0-7.d0)/mbar
	tmshift=0.d0
	hdur=2.4d0*xm00**(1.d0/3.d0)/1.d6
	
	write(*,'(a,f5.2,a,3(1x,f7.2),2(1x,f6.2))') '    Mw: ',xm0,
     1	  '. Fault-plane, tmshift and hdur: ',xm(1),xm(2),xm(3),tmshift,hdur
	
c   Then obtain the moment tensor elements in the geographic spherical polar 
c   system at the source in which r - up; theta - south; phi - east.

	mg(1,1)=m33
	mg(2,2)=m22
	mg(3,3)=m11
	mg(1,2)=-m23
	mg(1,3)=m13
	mg(2,3)=-m12
	mg(2,1)=mg(1,2)
	mg(3,1)=mg(1,3)
	mg(3,2)=mg(2,3)

1	do i=1,3
	  do j=1,3
	    mg(i,j)=mg(i,j)*xm00
	  end do
	end do
	
	write(*,'(a,6(1x,e10.3))') '    Moment tensor in GCS: ',
     1	  mg(1,1),mg(2,2),mg(3,3),mg(1,2),mg(1,3),mg(2,3)  
	
	return        
	end         
c
c
        subroutine glegendre(distan)

c   This subroutine computes the generalized Legendre function (GLF) denoted 
c   by P_lm^N(cos_theta) as defined in Dahlen & Tromp (1998). It is calculated 
c   for theta=distan and for l from 1 to lmax0. The seismologically relavent
c   functions are those with -2<=m<=2 and -2<=N<=2. The function is 
c   represented by a 3-dimensional array with the first index indicating 
c   l from 1 to lmax0, the second index from 1 to 5 indicating m from -2 to 2, 
c   and the third index from 1 to 5 indicating N from -2 to 2.
c
c   The explicit formulae are used in calculating P_lm^N for l<=2 as well as
c   P_lm^N for l=3 and m=2. All the other functions are computed using the
c   recursive in my notes on GLFs and standard symmetry relations.
c
c   Note: when the recursive relation is used, the round-off error quickly
c   accumulates and the result is unreliable for large l. In this case, the
c   symmetry relations must be used as extensively as possible to reduce the
c   usage of recursive relation.

        use modesummodule
        double precision theta,stheta,ctheta,recur1,recur2,recur3,distan,xl
        
c   Initializing the arrays for generalized Legendre functions.

        p0=0.d0
            
c   The generalized Legendre function is computed for the angle of the
c   epicentral distance distan.

        theta=distan*pi/180.d0
        stheta=dsin(theta)
        ctheta=dcos(theta)          

c   For functions P_lm^N when l=1.

        glf0(1,2,2)=(1.d0+ctheta)/2.d0
        glf0(1,2,3)=stheta/dsqrt(2.d0)
        glf0(1,2,4)=(1.d0-ctheta)/2.d0
        glf0(1,3,2)=-stheta/dsqrt(2.d0)
        glf0(1,3,3)=ctheta
        glf0(1,3,4)=stheta/dsqrt(2.d0)
        glf0(1,4,2)=(1.d0-ctheta)/2.d0
        glf0(1,4,3)=-stheta/dsqrt(2.d0)
        glf0(1,4,4)=(1.d0+ctheta)/2.d0
          
c   For functions P_lm^N when l=2.

        glf0(2,1,1)=(1.d0+ctheta)**2/4.d0
        glf0(2,1,2)=stheta*(1.d0+ctheta)/2.d0
        glf0(2,1,3)=dsqrt(6.d0)*stheta**2/4.d0
        glf0(2,1,4)=stheta*(1.d0-ctheta)/2.d0
        glf0(2,1,5)=(1.d0-ctheta)**2/4.d0
        glf0(2,2,1)=-stheta*(1.d0+ctheta)/2.d0
        glf0(2,2,2)=(1.d0+ctheta)*(2.d0*ctheta-1.d0)/2.d0
        glf0(2,2,3)=dsqrt(6.d0)*stheta*ctheta/2.d0
        glf0(2,2,4)=(1.d0-ctheta)*(2.d0*ctheta+1.d0)/2.d0
        glf0(2,2,5)=stheta*(1.d0-ctheta)/2.d0
        glf0(2,3,1)=dsqrt(6.d0)*stheta**2/4.d0
        glf0(2,3,2)=-dsqrt(6.d0)*stheta*ctheta/2.d0
        glf0(2,3,3)=(3.d0*ctheta**2-1.d0)/2.d0
        glf0(2,3,4)=dsqrt(6.d0)*stheta*ctheta/2.d0
        glf0(2,3,5)=dsqrt(6.d0)*stheta**2/4.d0
        glf0(2,4,1)=-stheta*(1.d0-ctheta)/2.d0
        glf0(2,4,2)=(1.d0-ctheta)*(2.d0*ctheta+1.d0)/2.d0
        glf0(2,4,3)=-dsqrt(6.d0)*stheta*ctheta/2.d0
        glf0(2,4,4)=(1.d0+ctheta)*(2.d0*ctheta-1.d0)/2.d0
        glf0(2,4,5)=stheta*(1.d0+ctheta)/2.d0
        glf0(2,5,1)=(1.d0-ctheta)**2/4.d0
        glf0(2,5,2)=-stheta*(1.d0-ctheta)/2.d0
        glf0(2,5,3)=dsqrt(6.d0)*stheta**2/4.d0
        glf0(2,5,4)=-stheta*(1.d0+ctheta)/2.d0
        glf0(2,5,5)=(1.d0+ctheta)**2/4.d0
              
c   For functions P_lm^N when l=3 and m=2.

        glf0(3,5,1)=(3.d0*ctheta+2.d0)*(1.d0-ctheta)**2/4.d0
        glf0(3,5,2)=-dsqrt(10.d0)*stheta*(3.d0*ctheta+1.d0)*
     1    (1.d0-ctheta)/8.d0
        glf0(3,5,3)=dsqrt(30.d0)*stheta**2*ctheta/4.d0
        glf0(3,5,4)=-dsqrt(10.d0)*stheta*(3.d0*ctheta-1.d0)*
     1    (1.d0+ctheta)/8.d0
        glf0(3,5,5)=(3.d0*ctheta-2.d0)*(1.d0+ctheta)**2/4.d0
          
c   The symmetry relations are used to obtain P_lm^N when l=3 and m=-2.

        glf0(3,1,1)=glf0(3,5,5)
        glf0(3,1,2)=-glf0(3,5,4)
        glf0(3,1,3)=glf0(3,5,3)
        glf0(3,1,4)=-glf0(3,5,2)
        glf0(3,1,5)=glf0(3,5,1)
          
c   Recursive relation is used to obtain P_lm^N for (3,0,-1), (3,0,0),
c   (3,0,1), (3,1,-1) and (3,1,1).
          
        l=2
        m=3
        m1=m-3
        do n=2,4
          n1=n-3
          recur1=dble(l+1)*dsqrt(dble(l-m1))/
     1      dsqrt(dble(l+n1+1)*dble(l-n1+1)*dble(l-m1+1))
          recur2=dsqrt(dble(l+m1+1)/dble(l-m1))
          recur3=dble(n1)/dble(l+1)
          glf0(l+1,m,n)=recur1*(stheta*glf0(l,m+1,n)+
     1      recur2*(ctheta-recur3)*glf0(l,m,n))
        end do
              
        m=4
        m1=m-3
        do n=2,4,2
          n1=n-3
          recur1=dble(l+1)*dsqrt(dble(l-m1))/
     1      dsqrt(dble(l+n1+1)*dble(l-n1+1)*dble(l-m1+1))
          recur2=dsqrt(dble(l+m1+1)/dble(l-m1))
          recur3=dble(n1)/dble(l+1)
          glf0(l+1,m,n)=recur1*(stheta*glf0(l,m+1,n)+
     1      recur2*(ctheta-recur3)*glf0(l,m,n))
        end do

c   Symmetry relations are used to obtain other functions P_lm^N.

        glf0(l+1,4,3)=-glf0(l+1,3,4)
        glf0(l+1,4,1)=glf0(l+1,5,2)
        glf0(l+1,4,5)=glf0(l+1,1,2)
        glf0(l+1,3,1)=glf0(l+1,1,3)
        glf0(l+1,3,5)=glf0(l+1,5,3)
        glf0(l+1,2,1)=glf0(l+1,5,4)
        glf0(l+1,2,2)=glf0(l+1,4,4)
        glf0(l+1,2,3)=glf0(l+1,3,4)
        glf0(l+1,2,4)=glf0(l+1,4,2)
        glf0(l+1,2,5)=glf0(l+1,1,4)
            
c   Now ready to calculate functions P_lm^N with l>3.

        do l=3,lmax1-1
	
c   First calculate functions with m=-2 using recursive relation.

	  if(mod(l,5000).eq.0) write(*,'(a,i7)') '    Processing l=',l

          m=1
          m1=m-3
          do n=1,5
            n1=n-3
            recur1=dble(l+1)*dsqrt(dble(l-m1))/
     1        dsqrt(dble(l+n1+1)*dble(l-n1+1)*dble(l-m1+1))
            recur2=dsqrt(dble(l+m1+1)/dble(l-m1))
            recur3=dble(n1)/dble(l+1)
            glf0(l+1,m,n)=recur1*(stheta*glf0(l,m+1,n)+
     1        recur2*(ctheta-recur3)*glf0(l,m,n))
          end do
                 
c   For m=2, the functions are obtained by symmetry relation.

          glf0(l+1,5,1)=glf0(l+1,1,5)
          glf0(l+1,5,2)=-glf0(l+1,1,4)
          glf0(l+1,5,3)=glf0(l+1,1,3)
          glf0(l+1,5,4)=-glf0(l+1,1,2)
          glf0(l+1,5,5)=glf0(l+1,1,1)
            
c   Again, recursive relation is used to obtain P_lm^N for (l,0,-1), (l,0,0),
c   (l,0,1), (l,1,-1) and (l,1,1).

          m=3
          m1=m-3
          do n=2,4
            n1=n-3
            recur1=dble(l+1)*dsqrt(dble(l-m1))/
     1        dsqrt(dble(l+n1+1)*dble(l-n1+1)*dble(l-m1+1))
            recur2=dsqrt(dble(l+m1+1)/dble(l-m1))
            recur3=dble(n1)/dble(l+1)
            glf0(l+1,m,n)=recur1*(stheta*glf0(l,m+1,n)+
     1        recur2*(ctheta-recur3)*glf0(l,m,n))
          end do
              
          m=4
          m1=m-3
          do n=2,4,2
            n1=n-3
            recur1=dble(l+1)*dsqrt(dble(l-m1))/
     1        dsqrt(dble(l+n1+1)*dble(l-n1+1)*dble(l-m1+1))
            recur2=dsqrt(dble(l+m1+1)/dble(l-m1))
            recur3=dble(n1)/dble(l+1)
            glf0(l+1,m,n)=recur1*(stheta*glf0(l,m+1,n)+
     1        recur2*(ctheta-recur3)*glf0(l,m,n))
          end do

c   Other functions are obtained by symmetry relations.
              
          glf0(l+1,4,3)=-glf0(l+1,3,4)
          glf0(l+1,4,1)=glf0(l+1,5,2)
          glf0(l+1,4,5)=glf0(l+1,1,2)
          glf0(l+1,3,1)=glf0(l+1,1,3)
          glf0(l+1,3,5)=glf0(l+1,5,3)
          glf0(l+1,2,1)=glf0(l+1,5,4)
          glf0(l+1,2,2)=glf0(l+1,4,4)
          glf0(l+1,2,3)=glf0(l+1,3,4)
          glf0(l+1,2,4)=glf0(l+1,4,2)
          glf0(l+1,2,5)=glf0(l+1,1,4)

c   All P_lm^N are computed for the current l.

        end do

c   End of loop over l. All glf are computed for the current thetasq.

c   Converting P_lm^N to X_lm^N.

        do l=1,lmax1
          xl=dsqrt((2.d0*dble(l)+1.d0)/(4.d0*pi))
          do m=1,5
            do n=1,5
              glf0(l,m,n)=xl*glf0(l,m,n)
            end do
          end do
        end do
        
        return
        end
c
c
        subroutine modect

c   This subroutine reads the normal-mode summary file and the stripped
c   eigenfunction files for the source and receiver radii.

        use modesummodule
        character(len=8) adep*8
        integer rlen,iczero,icone
        double precision qtbl,ptbl,ctbl,dn,eta,pv,ph,sh,rq,qa,qb
	double precision xdiff,r1x,rs0,rr0,rn
	double precision depths,depthr,rob
	double precision ffmin,ffmax

        rlen=5*2
        if(amode.ne.'sph') then
          rlen=2*2
        endif
	
	weight=0.
	iczero=0
	icone=1
	
c   First open mode table header file and read header

	k=index(modepath,' ')-1
	k1=index(modefile,' ')-1
	if(modepath(k:k).ne.'/') then
	  k=k+1
	  modepath(k:k)='/'
	endif
	write(*,'(/a,a)') '    Reading: ',
     1	  modepath(1:k)//modefile(1:k1)//'.'//amode
        open(20,file=modepath(1:k)//modefile(1:k1)//'.'//amode,status='old',
     1    form='unformatted',access='sequential')
        read(20) jcom,wwmin,wwmax,nmodes0,nnmin,nnmax,llmin,llmax,wgrav	
        read(20) rn

c   r() contains the radial grids at which the eigenfunctions are stored 
c   during stripping.

        read(20) knots,(rknots(i1),i1=1,knots)
	
c   rad() contains the radial grids given in the reference model file. For 
c   spheroidal modes, rad() starts from the center of the earth. For toroidal 
c   modes, it starts at CMB+.

        read(20) nrad,(rmod(i1),i1=1,nrad),(dn,i1=1,nrad),
     1	  (pv,i1=1,nrad),(sv(i1),i1=1,nrad),(ph,i1=1,nrad),
     2	  (sh,i1=1,nrad),(eta,i1=1,nrad)
        read(20) nqk,(rq,i1=1,nqk),(qa,i1=1,nqk),(qb,i1=1,nqk)
	if((kf.eq.1).or.(kf.eq.(nfiles+1)).or.(kf.eq.(nfiles+nfilet+1))) then
          do i=1,nm
            read(20,end=100) ntbl(i),lk(i),fx(i),qtbl,atten(i),ptbl,
     1	      ctbl,gtbl,etbl
	    
c   For the first file being read, set the weight of every mode to 1.

	    weight(i)=1.
	    
c   When a mode is read (and used in summation), set icount(n,l)=1.

	    icount(ntbl(i),lk(i))=icone
          end do
	else
          do i=1,nm
            read(20,end=100) ntbl(i),lk(i),fx(i),qtbl,atten(i),ptbl,ctbl,
     1	      gtbl,etbl
	    
c   If mode (n,l) has icount(n,l)=0, then set its weight to 1 so that it's  
c   used in the summation.

	    if(icount(ntbl(i),lk(i)).eq.iczero) weight(i)=1.
	    
c   If mode (n,l) has icount(n,l)=1, which means the mode has already been 
c   use in the summation (the mode is also contained in a previous mode 
c   table file), its weight is set to 0 so that it won't be used again.

	    icount(ntbl(i),lk(i))=icone
          end do
	endif
100     close(20)

	nmodes=i-1
	nrep=0
	do i=1,nmodes
	  if(weight(i).eq.0.) nrep=nrep+1
	end do
	nrept=nrept+nrep
	
	ffmin=1.d10
	ffmax=-1.d0
	do i=1,nmodes
	  if(ffmin.gt.fx(i)) ffmin=fx(i)
	  if(ffmax.lt.fx(i)) ffmax=fx(i)
	end do
	
        write(*,'(a,1x,i9,2(1x,i6,1x,i6))')
     1    '    Mode number, min-max n and l in mode table:  ',
     2	  nmodes,nnmin,nnmax,llmin,llmax
 	write(*,'(a,i9)') '    Number of repeated modes: ',nrep
	write(*,'(a,2(2x,f9.6))') '    Min. and Max. frequencies (Hz): ',
     1	  ffmin/(2.d0*pi),ffmax/(2.d0*pi)	
     
c   Find out the node numbers of ICB and CMB.
	
	do i=1,nrad
	  licb=i
	  if(sv(i).lt.1.d-5) goto 1
	end do
1	do i=licb,nrad
	  lcmb=i
	  if(sv(i).gt.1.d-5) goto 2
	end do
2	continue
	if(jcom.eq.2) licb=0
	if(jcom.eq.2) lcmb=1
	write(*,'(/a,3(2x,i6))') 
     1	  '    Total, ICB and CMB nodes: ',nrad,licb,lcmb
	
c   Find seafloor or surface boundary.

        lomb=nrad
        do i=nrad,lcmb,-1
          if(sv(i).gt.1.d-10) then
            lomb=i
            goto 55
          end if
        end do
        stop '  ***PROBLEM IN FINDING SEAFLOOR OR SURFACE BOUNDARY***'
55      continue         
        rob=rmod(lomb)
	
        depthr=-stelv/1.d3
	if(depthr.lt.0.d0) depthr=0.d0
        depths=depths0
	if(rob.lt.rnorm) then
	  write(*,'(a,f9.4)') '    Radius of sea floor: ',rob/1.d3
          if(depths0.le.((rnorm-rob)/1.d3)) then
            write(*,'(/a,3x,f7.3,a)') 
     1	    '    *** Source depth in the ocean: ',real(depths),'km'
          else
          endif
          if(depthr.le.((rnorm-rob)/1.d3)) then
            write(*,'(/a,3x,f7.3,a)') 
     1      '    *** Receiver depth in the ocean: ',real(depthr),'km'
          else
          endif
        else
	  write(*,'(a)') '    This is a continental model.'
        endif
	
c   Convert source and receiver depths to radii. Note: In input file source 
c   depth depths0 is given in km (downward from 0 at surface), and station 
c   elevation stelv is given in m (upward from 0 at surface).

        rs0=rnorm-depths*1.d3
        rr0=rnorm-depthr*1.d3
	
	write(*,'(a,2(4x,f9.1))') '    Final source and station radii:   ',
     1	  rs0,rr0

	if((jcom.eq.3).or.(jcom.eq.1)) goto 151

        if((rs0.lt.rmod(1)).or.(rr0.lt.rmod(1))) then
          nmodes=0
          write(*,'(a,3(1x,f10.5))') 
     1	    '    Source or receiver below the CMB. Toroidal modes ignored.',
     2      rmod(1)/1.d3,rs0/1.d3,rr0/1.d3
          return
        endif

151     continue

        xdiff=1.d20
	r1x=rs0/rnorm
        do 15 i=1,knots
          if(dabs(rknots(i)-r1x).le.xdiff) then
	    xdiff=dabs(rknots(i)-r1x)
            rs=rknots(i)
            js=i
          endif
15      continue
        ir=int(rs*rnorm)
	write(adep(1:7),'(i7.7)') ir
        adep(8:8)='0'
        if(js.eq.1) then
        elseif((rnorm*dabs(rknots(js)-rknots(js-1))).lt.1.d-2) then
          if(rs0.ge.(rs*rnorm)) adep(8:8)='1'
          if(rs0.lt.(rs*rnorm)) adep(8:8)='0'
        endif

c   Open the mode table - this assumes one using mineos_table.
c   kmode determines the record length: 4 bytes for every mode.

	write(*,'(/a,a)') '    S: ',
     1	  modepath(1:k)//modefile(1:k1)//'_'//adep//'.'//amode
        open(30,file=modepath(1:k)//modefile(1:k1)//'_'//adep//'.'//amode,
     1	  status='old',form='unformatted',access='direct',recl=rlen)
        do 10 j=1,nmodes
          if(jcom.eq.3) then
            read(30,rec=j) uksr(j,1),duksr(j,1),vksr(j,1),xksr(j,1),yksr(j,1)
          else
            read(30,rec=j) uksr(j,1),duksr(j,1)
          endif
10      continue
        close(30)

        xdiff=1.d20
	r1x=rr0/rnorm
        do 25 i=1,knots
          if(dabs(rknots(i)-r1x).le.xdiff) then
            xdiff=dabs(rknots(i)-r1x)
            rr=rknots(i)
            jr=i
          endif
25      continue
        ir=int(rr*rnorm)
	write(adep(1:7),'(i7.7)') ir
        adep(8:8)='0'
        if(jr.eq.1) then
        elseif((rnorm*dabs(rknots(jr)-rknots(jr-1))).lt.1.d-2) then
          if(rr0.ge.(rr*rnorm)) adep(8:8)='1'
          if(rr0.lt.(rr*rnorm)) adep(8:8)='0'
        endif

c   Open the mode table - this assumes one using mineos_table.
c   kmode determines the record length: 4 bytes for every mode.

	write(*,'(a,a)') '    R: ',
     1	  modepath(1:k)//modefile(1:k1)//'_'//adep//'.'//amode
        open(30,file=modepath(1:k)//modefile(1:k1)//'_'//adep//'.'//amode,
     1	  status='old',form='unformatted',access='direct',recl=rlen)
        do 40 j=1,nmodes
          if(jcom.eq.3) then
            read(30,rec=j) uksr(j,2),duksr(j,2),vksr(j,2),xksr(j,2),yksr(j,2)
          else
            read(30,rec=j) uksr(j,2),duksr(j,2)
          endif
40      continue
        close(30)
	
        return
        end
c
c
	subroutine azimth(ellips,slat,slon,rlat,rlon,delta,azim,bazim)

c   This routine uses Euler angles to find the geocentric distance,
c   azimuth, and back azimuth for a source-reciever pair.
c
c   Input
c
c     slat  - source geographic latitude in decimal degrees
c     slon  - source longitude in decimal degrees
c     rlat  - reciever geographic latitude in decimal degrees
c     rlon  - reciever longitude in decimal degrees
c
c   Output
c
c     delta - geocentric source-reciever distance in decimal degrees of arc
c     azim  - geocentric azimuth from the source to the reciever
c     bazim - geocentric back azimuth from the reciever to the source
c
c   The distance calculated here delta is always between 0 and 180 degrees. 
c   Accordingly, the azimuth and back azimuth are defined for the minor 
c   arc between (slat,slon) and (rlat,rlon).

        implicit real*8 (a-h,o-z)
	integer ellips

        data flt/298.25d0/

	pi=4.d0*datan(1.d0)
	dtor=pi/180.d0
	
	if(ellips.ne.0) then
          e=1.d0/flt
	else
	  e=0.d0
	endif
	
c   Convert to geocentric coordinates and from latitude to colatitude.

        slatra=dtor*slat
        w=dsin(slatra)
        s=((2.d0-e)*w+4.d0*e*(w**3))*e*dcos(slatra)
        scolat=pi/2.d0-slatra+s
        rlatra=dtor*rlat
        w=dsin(rlatra)
        s=((2.d0-e)*w+4.d0*e*(w**3))*e*dcos(rlatra)
        rcolat=pi/2.d0-rlatra+s

        slonra=slon*dtor
        rlonra=rlon*dtor
        c2=dcos(scolat)
        s2=dsin(scolat)
        c1=dcos(slonra)
        s1=dsin(slonra)
        slatrc=dsin(rcolat)

c   Find the azimuth and distance by rotating the source to the north pole.

        x0=slatrc*dcos(rlonra)
        y0=slatrc*dsin(rlonra)
        z0=dcos(rcolat)
        x1=c1*x0+s1*y0

        z0=dcos(rcolat)
        x1=c1*x0+s1*y0
        y1=-s1*x0+c1*y0
        z1=z0
        x2=c2*x1-s2*z1
        y2=y1
        z2=c2*z1+s2*x1
        call angles(x2,y2,z2,delta,azim)
        azim=180.d0-azim

c   Find the back azimuth by rotating the receiver to the north pole.

        c2=dcos(rcolat)
        s2=dsin(rcolat)
        c1=dcos(rlonra)
        s1=dsin(rlonra)
        slatrc=dsin(scolat)
        x0=slatrc*dcos(slonra)
        y0=slatrc*dsin(slonra)
        z0=dcos(scolat)
        x1=c1*x0+s1*y0
        y1=-s1*x0+c1*y0
        z1=z0
        x2=c2*x1-s2*z1
        y2=y1
        z2=c2*z1+s2*x1
        call angles(x2,y2,z2,delta,bazim)
        bazim=180.d0-bazim

        return
        end
c
c
        subroutine angles(x,y,z,theta,phi)

c   Finds the angles theta and phi of a spherical polar coordinate
c   system from the cartesion coordinates x, y, and z.

        implicit real*8 (a-h,o-z)

        parameter(eps=1.d-14)

	pi=4.d0*datan(1.d0)

        rtod=180.d0/pi
        arg1=dsqrt(x*x+y*y)
        theta=datan2(arg1,z)
        if(dabs(x).le.eps.and.dabs(y).le.eps) then
          phi=0.d0
        else
          phi=datan2(y,x)
        endif
        phi=phi*rtod
        theta=theta*rtod

        return
        end
