c+--------------------------------------------------------------------------------+
c   This program reads the normal-mode eigenfunction files computed by MINEOS.
c   For toroidal modes, the output functions are W and dW/dr. For spheroidal
c   modes, the output functions are U, dU/dr, V and dV/dr. All these functions
c   are defined and normalized in the same way as in Dahlen & Tromp (1998).
c
c   Input:  sequential-access binary file created by MINEOS containing the
c           normal-mode eigenfunctions U, V, dU/dr, dV/dr, Phi1 and Phi2
c           for spheroidal modes and W and dW/dr for toroidal modes.
c
c   Keyboard input: binary filename from MINEOS
c                   output filename
c                   n and l of the mode to be read.
c+--------------------------------------------------------------------------------+
c
        program read_mineos

        parameter(nrmax=100000,nrmax10=10*nrmax)
        character*200 funfile,outfile
        real*4 tdum(nrmax,9),r(nrmax)
	real*8 eig(nrmax10)
        real*8 u(nrmax),du(nrmax),v(nrmax),dv(nrmax)
        real*8 pi,wrad,gv,qmod,err
        real*8 dnorm,tnorm,fnorm,gravity
        data dnorm,gravity/5515.d0,6.6723d-11/

        pi=4.d0*datan(1.d0)
        tnorm=1.d0/dsqrt(pi*dnorm*gravity)
        fnorm=1.d0/tnorm

        write(*,'(/a)') '  Input binary eigenfunction file (from MINEOS):'
        read(*,'(a)') funfile

        write(*,'(/a)') '  Output ASCII eigenfunction file:'
        read(*,'(a)') outfile

        write(*,'(/a)') '  n and l of the mode to be read:'
        read(*,*) n0,l0

c   Reading the header of the eigenfunction file for info on reference model.

        open(1,file=funfile,status='old',
     1    form='unformatted',access='sequential')

        read(1) jcom,wmin,wmax,lmin,lmax,wgrav
        read(1) n,nic,noc,ifanis,trefl,((tdum(i,j),i=1,n),j=1,9)
	
        lbottom=0
        if(jcom.eq.3) then
          mrec=6*n+7
          nr=n
        elseif(jcom.eq.1) then
          mrec=2*n+7
          nr=n
        else
          mrec=2*(n-noc)+7
          nr=n-noc
          lbottom=noc
        endif

c   Start reading and processing the eigenfunctions mode by mode.

2       read(1,err=5,end=5) (eig(j),j=1,mrec)
	nn=int(eig(1))
	ll=int(eig(2))
	wrad=eig(3)
	qmod=eig(4)
	gv=eig(5)
	err=eig(6)
		
        if((nn.ne.n0).or.(ll.ne.l0)) goto 2
        close(1)

c   Renormalizing the eigenfunctions according to Dahlen & Tromp (1998).

        do i=8,mrec
          eig(i-7)=eig(i)*(wrad/fnorm)
        end do

c   Processing the eigenfunctions.

        if(jcom.eq.3) then
          do i=1,nr
            r(i)=tdum(lbottom+i,1)
            u(i)=dble(eig(i))
            du(i)=dble(eig(n+i))
            v(i)=dble(eig(2*n+i))
            dv(i)=dble(eig(3*n+i))

c   Setting numerically nonzero but insignificantly small functional values
c   to exactly zero to avoid later numerical problems.

            if(dabs(u(i)).le.1.d-37) u(i)=0.d0
            if(dabs(du(i)).le.1.d-37) du(i)=0.d0
            if(dabs(v(i)).le.1.d-37) v(i)=0.d0
            if(dabs(dv(i)).le.1.d-37) dv(i)=0.d0
          end do

          write(*,'(a,i4,1x,a,1x,i4,1x,f20.14,2(1x,f10.3),1x,e15.8)') 'Mode found: ',
     1	    nn,'S',ll,wrad,gv,qmod,err
          open(2,file=outfile,status='unknown')
          write(2,'(a,i4,1x,a,1x,i4,1x,f20.14,2(1x,f10.3),1x,e15.8)') 'Mode found: ',
     1	    nn,'S',ll,wrad,gv,qmod,err
          do i=1,nr
            write(2,'(f12.4,5(1x,e15.8))') r(i),u(i),du(i),v(i),dv(i)
          end do
          close(2)
          goto 15

        elseif(jcom.eq.1) then
          do i=1,nr
            r(i)=tdum(lbottom+i,1)
            u(i)=dble(eig(i))
            du(i)=dble(eig(nr+i))

c   Setting numerically nonzero but insignificantly small functional values
c   to exactly zero to avoid later numerical problems.

            if(dabs(u(i)).le.1.d-37) u(i)=0.d0
            if(dabs(du(i)).le.1.d-37) du(i)=0.d0
          end do

          write(*,'(a,i4,1x,a,1x,i4,1x,f20.14,2(1x,f10.3),1x,e15.8)') 'Mode found: ',
     1	    nn,'R',ll,wrad,gv,qmod,err
          open(2,file=outfile,status='unknown')
          write(2,'(a,i4,1x,a,1x,i4,1x,f20.14,2(1x,f10.3),1x,e15.8)') 'Mode found: ',
     1	    nn,'R',ll,wrad,gv,qmod,err
          do i=1,nr
            write(2,'(f12.4,5(1x,e15.8))') r(i),u(i),du(i)
          end do
          close(2)
          goto 15

        else

          do i=1,nr
            r(i)=tdum(lbottom+i,1)
            u(i)=dble(eig(i))
            du(i)=dble(eig(nr+i))

c   Setting numerically nonzero but insignificantly small functional values
c   to exactly zero to avoid later numerical problems.

            if(dabs(u(i)).le.1.d-37) u(i)=0.d0
            if(dabs(du(i)).le.1.d-37) du(i)=0.d0
          end do

          write(*,'(a,i4,1x,a,1x,i4,1x,f20.14,2(1x,f10.3),1x,e15.8)') 'Mode found: ',
     1	    nn,'T',ll,wrad,gv,qmod,err
          open(2,file=outfile,status='unknown')
          write(2,'(a,i4,1x,a,1x,i4,1x,f20.14,2(1x,f10.3),1x,e15.8)') 'Mode found: ',
     1	    nn,'T',ll,wrad,gv,qmod,err
          do i=1,nr
            write(2,'(f12.4,5(1x,e15.8))') r(i),u(i),du(i)
          end do
          close(2)
          goto 15

        endif

5       continue

        write(*,'(/a/)') '   Requested mode not found.'
        stop

15      continue

        stop
        end
