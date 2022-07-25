c   This program checks the normal-mode eigenfrequency file created by MINEOS
c   to see if there are following problems:
c
c      1. Whether there are any modes missing including
c
c         -- if the lowest mode 1S1 is calculated;
c         -- if all the fundamental modes are calculated;
c         -- if there is any gap in n numbers for each l.
c
c      2. Whether there are any modes with negative n.
c
c   Note: Even after this check, there could still be problems in the results,
c   e.g. sometimes everything seems fine for a mode except that its group
c   velocity and/or Q values are NaNs.
c  2006/4/24, 2006/11/17  slightly modified by Ysinging
c  2007/3/30 read eigenfreq and check if there is any freq jump [fixed l]
c  2007/10/24 capable to read new format (ftype=3) and check error
c  2007/11/8 output mode with error n 
c  2008/12/23 capable to read new format (ftype=4)

	program checkmode
	dimension n(550000),l(550000),fmhz(550000),err(550000),qmod(550000),gmod(550000)
        INTEGER, PARAMETER :: nlayer=50000, fid=20, fidf=21
        INTEGER :: ftype, Ier, choix, ifmt
        REAL :: tmp, rayp, mhz2w
        CHARACTER(LEN=128) :: filename, outname, outnamef, cdum*10
        LOGICAL :: iswrong1, iswrong2
	double precision crit_err

        iswrong1=.FALSE.
        iswrong2=.FALSE.
c        crit_err = 1.0e-5
        
c   Get the eigenfrequency file name from the keyboard.

c        WRITE(*,*) 'Type of input file: (before merge[1],
c    &              after merge[2], another[3], new fre format[4]):'
        READ (*,*)  ftype 
c       write(*,*) 'input filename (the immediate MINEOS output file *.fre):'
	read(*,'(a)') filename
c       WRITE(*,*) 'output filename : '
        READ(*,'(A)') outname
c	WRITE(*,*) 'filename for modes with large Rayleigh quotients'
	READ(*,'(A)') outnamef
c       WRITE(*,'(A/A)') "Starting program...","Check the size of modes...>>>"

c   write out the missing mode for fixed n but jumped l [the last block of if-statment]
c       WRITE(*,*) 'Wish to write part of the missing modes...' 
c       WRITE(*,*) 'if freq. starts from very low freq. - recommand[0]... else [1]'
        READ(*,*) choix 

	READ(*,*) crit_err

c   Open the eigenfrequency file and read the header.

	k=-1
        OPEN(1,FILE=filename, STATUS='old',err=3)	
        IF ( ftype == 1 ) THEN
	  do i=1, nlayer
	    read(1,'(a10)') cdum
	    if(cdum(7:10).eq.'MODE'.or.cdum(7:10).eq.'mode') goto 1
	  end do
1	  read(1,'(a10)') cdum
	  write(*,'(a)') cdum, 'Blank is right!!'
        END IF

c   Start reading the mode branch numbers and angular degrees one by one.
        IF ( ftype == 1 ) THEN
          i = 1
	  DO
            READ(1,100,IOSTAT=Ier) n(i), cdum, l(i), tmp, fmhz(i)
            IF ( Ier < 0 ) THEN
c             WRITE(*,*) "<<< Check successfully"
              EXIT
            END IF   
            i = i + 1
	  END DO
        ELSE IF ( ftype == 2 ) THEN
          i = 1
          DO
            READ(1,*,IOSTAT=Ier) n(i), l(i), tmp, fmhz(i)
            IF ( Ier < 0 ) THEN
c             WRITE(*,*) "<<< Check successfully"
              EXIT
            END IF
            i = i + 1
          END DO
        ELSE IF (ftype == 3) THEN
          i = 1
          DO
            READ(1,*,IOSTAT=Ier) n(i), l(i), fmhz(i), tmp, tmp, err(i)
            IF ( Ier < 0 ) THEN
c             WRITE(*,*) "<<< Check successfully"
              EXIT
            END IF
            i = i + 1
          END DO
        ELSE IF (ftype ==4) THEN
          i=1
          DO
            READ(1,*,IOSTAT=Ier) n(i), cdum, l(i), tmp, fmhz(i), tmp, qmod(i), gmod(i), err(i)
            WRITE(*,*) n(i), l(i), fmhz(i), err(i)
            IF (Ier<0) THEN
c             WRITE(*,*) "<<< Check successfully"
              EXIT
            END IF
            i=i+1
          END DO
        ELSE
          STOP 'check type of input file'
        END IF
2	close(1)

        k=i-1
        WRITE(*,*) '# of data=', k
        IF(k==0) THEN
	  lenf=index(filename,' ')-1
	  write(*,*) 'No modes in specified file: ',filename(1:lenf)
	  STOP
	ENDIF

3	continue
	IF(k<0) THEN
	  lenf=index(filename,' ')-1
	  write(*,*) 'Specified file: ',filename(1:lenf),' does not exist.'
	  STOP
	ENDIF

        OPEN(UNIT=fid,FILE=outname)
	OPEN(UNIT=fidf,FILE=outnamef)
        WRITE(UNIT=fid,FMT='(A)') '   Checking:  '//TRIM(filename)
c   Write out the problem mode
     

       IF ((ftype == 3).OR.(ftype == 4)) THEN
         mhz2w=2.0*3.1415926/1000.0

cccccccc frequency error 
         DO i = 1, k
           IF (ABS(err(i)) > crit_err ) THEN 
             iswrong1=.TRUE.
             rayp = (l(i)+0.5)/(mhz2w*fmhz(i))
	     fmhz1=0.
	     fmhz2=0.
	     if(i.gt.1) fmhz1=fmhz(i-1)
	     if(i.lt.k) fmhz2=fmhz(i+1)
             WRITE(UNIT=fidf,FMT='(A,I5,1X,I5,1X,F10.5,1X,E11.4,2(1X,F10.5),3(1X,F9.2))')
     &       'Inaccurate: ', n(i), l(i), fmhz(i), err(i), fmhz1, fmhz2, rayp, gmod(i), qmod(i) 
           END IF
         END DO
ccccccccc problem for mode finding 
         IF (Ier > 0 ) STOP '***open wrongn file error***'
         DO i = 2, k
           IF ( (l(i) == l(i-1)) .AND. (n(i) == n(i-1)+1) .AND. (fmhz(i) == fmhz(i-1)) ) THEN
             iswrong2=.TRUE.
             WRITE(UNIT=fid,FMT='(A,/,I5,2X,I5,2X,F12.2)') 
     &        'possible wrong n:',  n(i), l(i), fmhz(i)
           END IF
         END DO     
       END IF

c100	format(i5,a2,i5)
c101     FORMAT(I5,1X,I5)
100     FORMAT(I5,A2,I5,1X,F11.7,4X,F12.5)
	iflag=0
c   Reading is done. There are a total of k modes in the file.
c   First check if 1S1 is missing.
	
	if((n(1).ge.2).and.(l(1).eq.1)) then
	  write(fid,*) '   Missing first mode: ',1,1
	  iflag=iflag+1
	endif

	if((n(1).ge.1).and.(l(1).ge.2)) then
	  write(fid,*) '  Missing fundamental mode: ',0,l(1)
	  iflag=iflag+1
	endif
        
	nmax=-50000
	nmin=50000
	lmax=-50000
	lmin=50000
	do i=2,k
	  if((n(i)-n(i-1)).eq.1.and.l(i).eq.l(i-1)) then
	  elseif(n(i).lt.0) then
	    n(i)=n(i-1)+1
	    write(fid,*) '   Bad mode: ',n(i),l(i)
	    iflag=iflag+1
	  elseif(l(i).eq.l(i-1)) then
	    do nn=n(i-1)+1,n(i)-1
	      write(fid,*) '   Missing mode: ',nn,l(i)
	      iflag=iflag+1
	    end do
	  else
	    if(n(i).eq.0.and.(l(i)-l(i-1)).eq.1) then
	    elseif((l(i)-l(i-1)).gt.1) then
	      do ll=l(i-1)+1,l(i)-1
	        write(fid,*) '   Missing degree: ',ll
	        iflag=iflag+1
	      end do
	    else
              IF ( choix == 0 ) THEN
	      do nn=0,n(i)-1
	        write(fid,*) '  Possible  Missing mode: ',
     &                   0,l(i),' to ',nn,l(i)     
	        iflag=iflag+1
	      end do
              END IF
	    endif
	  endif
	  if(n(i).ge.0.and.n(i).lt.nmin) nmin=n(i)
	  if(n(i).ge.0.and.n(i).gt.nmax) nmax=n(i)
	  if(l(i).lt.lmin) lmin=l(i)
	  if(l(i).gt.lmax) lmax=l(i)
	  if(n(i).lt.0.or.l(i).lt.0) then
	    write(fid,*) '   Problematic mode: ',n(i),l(i)
	    write(fid,*) '   Previous mode: ',n(i-1),l(i-1)
	  endif

c  Check if there is any eigenfreqeuncy jump
          IF ((n(i) > n(i-1)) .AND. (l(i) == l(i-1)) .AND. (fmhz(i) < fmhz(i-1))) THEN
            WRITE(fid,'(A,I5,A,I5,I8)') '    Eigenfreq wrong: ',n(i-2),'-',n(i), l(i)  
          ENDIF           
	end do

	IF (lmin > l(1)) lmin = l(1)
        IF (nmin > n(1)) nmin = n(1)
        IF (iflag >= 1) THEN	
       	  write(fid,*) '   Mode checking is done. Branch and degree bounds: '
	  write(fid,*) '        n: ',nmin,nmax,'  l: ',lmin,lmax
	  write(fid,*) '   Number of missing modes: ',iflag
          CLOSE(fid)
c	  stop
        ELSE
c         IF ((iswrong1 .eqv. .FALSE.).AND.(iswrong2 .eqv. .FALSE.)) THEN 
	  IF(iswrong2 .eqv. .FALSE.) THEN
            CLOSE(UNIT=fid,STATUS='DELETE')
c            STOP
          END IF
        END IF

	if(iswrong1 .eqv. .FALSE.) THEN
	  close(fidf,status='delete')
	else
	  close(fidf)
	endif

	stop
	end
