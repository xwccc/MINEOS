c   This program merges the binary normal-mode eigenfunction files computed by a series of
c   MINEOS runs into a single binary file with the same format.
c   
c   Input files:   MINEOS generated normal mode eigenfunctions files (binary).
c
c   Output file:   Merged eigenfunction file (binary).

	parameter(nfmax=100,nmax=2000)
	character*120 frefile,funfile(nfmax),outfun,outfre,mtype*2
	integer nn(nmax)
        real, dimension (:,:), allocatable :: tdum
	double precision, dimension (:), allocatable :: buf
	double precision, dimension (:,:), allocatable :: eig
	double precision ndum,fmhz,pi,period

	pi=4.d0*datan(1.d0)
	
c   Obtain names of input and output files.

	write(*,'(/a)') '    Number of input files:'
	read(*,*) nfile
	
	write(*,'(/a)') '    Enter the filenames (.fun):'
	do i=1,nfile
	  read(*,'(a)') funfile(i)
	end do
	write(*,'(/a)') '    Output filename (.fun):'
	read(*,'(a)') outfun
	lo=index(outfun,' ')-1
	outfre(1:lo-3)=outfun(1:lo-3)
	outfre(lo-2:lo)='fre'
	lr=index(funfile(1),' ')-1
	frefile(1:lr-3)=funfile(1)(1:lr-3)
	frefile(lr-2:lr)='fre'

c   Read one input file to get the array dimensions.

	open(1,file=funfile(1),status='old',form='unformatted',access='sequential')
        read(1) jcom,wmin1,wmax1,llmin,llmax,wgrav
        read(1) nr,nic,noc
	close(1)

c   Number of saved eigenfunction values according to mode types.
	
        if(jcom.eq.3) then
          mrec=6*nr+7
	  mtype=' S'
        elseif(jcom.eq.1) then
          mrec=2*nr+7
	  mtype=' R'
        else
          mrec=2*(nr-noc)+7
	  mtype=' T'
        endif

c   Read the input file and find out the range of n.

	minn=1000000
	maxn=-1
	wmin0=1.e10
	wmax0=-1.
	do i=1,nfile
	  open(100,file=funfile(i),status='old',form='unformatted',
     1	    access='sequential')
          read(100) jcom,wmin1,wmax1,llmin,llmax,wgrav
	  if(wmin0.gt.wmin1) wmin0=wmin1
	  if(wmax0.lt.wmax1) wmax0=wmax1
          read(100) nr
6	  read(100,end=5) ndum
	  nx=int(ndum)
	  if(minn.gt.nx) minn=nx
	  if(maxn.lt.nx) maxn=nx
	  goto 6
5	  close(100)
	end do
	
	write(*,'(/a,2(1x,i4),a,2(1x,i5),a,2(1x,f8.3))') 
     1	  '    n:',minn,maxn,';  l:',llmin,llmax,';  f(mHz):',wmin0,wmax0
     
c   Declaring arrays.
	
	allocate(tdum(nr,9))
	allocate(buf(mrec))
	allocate(eig(0:(maxn-minn),mrec))
	
c   Reading the MINEOS eigenfunctions files and in case of repeating modes 
c   the one comes later replaces the previous one. 

	nm=0
	do i=1,nfile
	  open(100,file=funfile(i),status='old',form='unformatted',
     1	    access='sequential')
          read(100) jcom,wmin1,wmax1,llmin,llmax,wgrav
          read(100) nr,nic,noc,ifanis,trefl,((tdum(ii,jj),ii=1,nr),jj=1,9)
1	  read(100,end=2) buf
	  nx=int(buf(1))
	  if(nm.eq.0) then
	    nm=nm+1
	    nn(nm)=nx
	    do jj=1,mrec
	      eig(nx-minn,jj)=buf(jj)
	    end do
	    goto 1
	  else
	    do j=1,nm
	      if(nx.eq.nn(j)) then
		do jj=1,mrec
		  eig(nx-minn,jj)=buf(jj)
		end do
		goto 1
	      endif
	    end do
	    nm=nm+1
	    nn(nm)=nx
	    do jj=1,mrec
	      eig(nx-minn,jj)=buf(jj)
	    end do
	  endif
	  goto 1
2	  close(100)
	end do
	
c   Write the output.
	
	open(1,file=outfun,status='unknown',form='unformatted',access='sequential')
        write(1) jcom,wmin0,wmax0,llmin,llmax,wgrav
        write(1) nr,nic,noc,ifanis,trefl,((tdum(ii,jj),ii=1,nr),jj=1,9)
	do n=minn,maxn
	  write(1) (eig(n-minn,jj),jj=1,mrec)
	end do
	close(1)
	open(1,file=outfre,status='unknown')
	do n=minn,maxn
	  n1=int(eig(n-minn,1))
	  l1=int(eig(n-minn,2))
	  fmhz=eig(n-minn,3)/(2.d0*pi)
	  period=1.d0/fmhz
	  fmhz=1000.d0*fmhz
	  write(1,'(i5,a2,i5,7g16.7)') n1,mtype,l1,eig(n-minn,3),fmhz,period,
     1	    eig(n-minn,4),eig(n-minn,5),eig(n-minn,6),eig(n-minn,7)
	end do
	close(1)
	
     	lf=index(funfile(1),' ')-1
	lo=index(outfun,' ')-1
     	write(*,'(a,a/a,a/)') '    Replace ',funfile(1)(1:lf),
     1	  '    with    ',outfun(1:lo)
     	write(*,'(a,a/a,a/)') '    Replace ',frefile(1:lf),
     1	  '    with    ',outfre(1:lo)
		
	stop
	end
