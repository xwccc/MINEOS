program reformat
implicit none
character(100)::type,fre,pre,deg,filenm
integer::l,lmax,nmax,err,n,length,lf,nn
logical::alive
real::w,f
write(*,'(/a)') '    Enter the 1st MINEOS ascii output filename (.fre)'
        read(*,'(a)') fre
! get all modes
nmax=0
err=0
l=1
open(4,file='freq.txt',status='replace')
open(3,file=fre,status='old')
do while (.true.) 
   read(3,*,iostat=err) n,type,l,w,f
   if (err .ne. 0) then
        exit
   endif
   write(4,*) n,l,f
end do
close(3)
nmax=n
l=1
lf=index(fre,' ')-1
pre=fre(1:lf-8)
!write(*,*) pre
do while (.true.)
	l=l+1
	write(deg,"(i4.4)") l
!	write(*,*)deg
        filenm=trim(pre)//trim(deg)//'.fre'
!	write(*,*)filenm
	inquire(file=filenm, exist=alive)
	if (alive .eqv. .false.) then
		exit
    endif
    open(3,file=filenm,status='old')
    err=0
	do while (err==0) 
	   read(3,*,iostat=err) n,type,l,w,f
           if (err .ne. 0) then
            exit
           endif
	   write(4,*) n,l,f
	end do
	close(3)
enddo
close(4) 
lmax=l-1
! reshape
open(5,file='freqsortn.txt',status='replace')
do nn = 0,nmax
	open(4,file='freq.txt',status='old')
 	err=0
	do while (.true.) 
	   read(4,*,iostat=err) n,l,f
           if (err .ne. 0) then
            exit
           endif
	   if (n .eq. nn) then
	   write(5,*) n,l,f
           endif
	end do
	close(4)
enddo
close(5)
end
