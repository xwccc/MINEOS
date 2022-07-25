	program mineos_drv

c   Adapted from the MPI version mineos_mpi.
c
c   This programe computes a mode catalogue down to 3.3 seconds
c   in a few minutes on exactly 300 processor cores.
c
c   Author: Dimitri Komatitsch, University of Pau, France, August 2008
c
c   Adapted from MPI version by Li Zhao, CNRS-Toulouse, France. June 2009.
c
c   Modified by Li Zhao, Academia Sinica, Taipei, Taiwan. September 2008.
c
c     -- Number of processors made flexible.
c     -- Checkmode output files only for runs with missing modes.
c     -- MINEOS input parameters changed to command-line arguments.
c     -- Ranges of frequency and angular degree determined 
c        automatically via Jeans relation.
c
c   Modified by H.Y. Yang, Dec. 2008.
c
c     -- Add radial modes.
c     -- Separate different tasks into individual subroutines

	implicit none
	character(len=120) :: string,modelpath,modelname,modelfile,tablepath
	character(len=120) :: buffer(20)
	integer, parameter :: IOUT=35,IIN=36
	integer :: myrank,l,sizeprocs,ier,i,j,ll,skip,ndeg,ib,ie
	integer :: ncpu,nloop
	integer :: min_ang_deg,max_ang_deg,max_ang_deg1
	double precision :: grav,eps,fmhzmin0,fmhzmin1,fmhzmin,fmhzmax,crit_err
	double precision :: time
        CHARACTER(LEN=128) :: frename, funname, dummy
	integer :: ifanis,ifdeck,nmodel,nic,noc,lrec
	double precision :: tref,rmodel,dmodel,pmodel,smodel,pmax

	sizeprocs=1
	myrank=0

	ncpu=sizeprocs
	
c   Get the command-line arguments.
	
	do i=1,10
	  call getarg(i,buffer(i))
	end do
	read(buffer(1),'(a)') modelpath
	read(buffer(2),'(a)') modelname
	read(buffer(3),'(a)') tablepath
	read(buffer(4),*) grav
	read(buffer(5),*) min_ang_deg	
	read(buffer(6),*) max_ang_deg
	read(buffer(7),*) eps
	read(buffer(8),*) fmhzmin0
	read(buffer(9),*) fmhzmax
	read(buffer(10),*) crit_err
	
        CALL check_dir(modelpath)
	modelfile=TRIM(modelpath)//TRIM(modelname)
	CALL check_dir(tablepath)
	
c   Create the output directory

	call system('mkdir -p MINEOS_INPUT_FILES')
	call system('mkdir -p MINEOS_SCROUTPUT_FILES')
	call system('mkdir -p CHECKMODE_INPUT_FILES')
	call system('mkdir -p CHECKMODE_OUTPUT_FILES')
	call system('mkdir -p CHECKMODE_FREQERR_FILES')
	string='mkdir -p '//TRIM(tablepath)
	call system(string)

c   Master process reads model file and determine the maximum ray parameter.

	if(myrank.eq.0) then
	  open(IIN,file=modelfile,status='old')
	  read(IIN,'(a)') dummy
	  read(IIN,*) ifanis,tref,ifdeck
	  read(IIN,*) nmodel,nic,noc,lrec
	  do i=1,noc+1
	    read(IIN,'(a)') dummy
	  end do
	  do i=noc+2,nmodel
	    read(IIN,*) rmodel,dmodel,pmodel,smodel
	    if(smodel.ne.0.d0) pmax=rmodel/smodel
	  end do
	  close(IIN)
	endif

        string=TRIM(tablepath)//TRIM(modelname)
	if((min_ang_deg.eq.0).and.(max_ang_deg.eq.0)) goto 2000
	
c   Determine the maximum angular degree from maximum ray parameter by Jeans 
c   relation and further increase it by 10%.
	
	max_ang_deg1=int((fmhzmax/500.d0)*3.14159d0*pmax*1.1d0)
	if(max_ang_deg.gt.max_ang_deg1) max_ang_deg=max_ang_deg1
	
c   Compute for angular degree between l=min_ang_deg and l=max_ang_deg.
c   Remember that degree l starts at 1, not 0.
c   Compute l=1, ncpu+1, 2*ncpu+1, 3*ncpu+1... on process 0;
c   Compute l=2, ncpu+2, 2*ncpu+2, 3*ncpu+2... on process 1;
c   etc... in order to balance the computation load because computing 
c   low l degrees is far more expensive than computing high l degrees.

	nloop=int(max_ang_deg/ncpu)+1
	if(myrank.eq.0) then
	  write(6,'(i5,a,i3,a,i5,a,i5)') 
     1	    nloop,' loops on ',ncpu,' CPUs for l from ',
     2	    min_ang_deg,' to ',max_ang_deg
     	  call flush(6)
	endif

	ndeg=max_ang_deg-min_ang_deg
        ib=min_ang_deg
	DO j=myrank, ndeg, ncpu
	  l=min_ang_deg+j
          ie=ib+ncpu-1
          IF (ie > max_ang_deg) ie = max_ang_deg
          ll=ib+ie-l
          string=TRIM(tablepath)//TRIM(modelname)
	  
c   Determine the minimum frequency from l and maximum ray parameter.

	  fmhzmin1=500.d0*(dble(l)+0.5d0)/pmax/3.14159265d0
	  
c   Further reduce minimum frequency from Jean's prediction by a factor of 8.

	  fmhzmin1=fmhzmin1/8.d0
	  if(fmhzmin1.gt.fmhzmin0) fmhzmin0=fmhzmin1
	  IF(fmhzmin0.ge.fmhzmax) EXIT 	  
         
c   Spheroidal.

	  CALL frange(l,fmhzmin0,fmhzmin)

c   Run mineos for spheroidal.

          CALL xmineos(3,l,eps,fmhzmin,fmhzmax,grav,TRIM(modelfile),
     1	    TRIM(string),frename)
	          
c   Run checkmode for spheroidal

          CALL xcheckmode(3,l,TRIM(frename),crit_err)

c   Toroidal.

	  CALL frange(ll,fmhzmin0,fmhzmin)
	  
c   Run mineos for toroidal.

          CALL xmineos(2,ll,eps,fmhzmin,fmhzmax,grav,TRIM(modelfile),
     1	    TRIM(string),frename)
	  
c   Run checkmode for toroidal 

          CALL xcheckmode(2,ll,TRIM(frename),crit_err)       
 
	  ib = ie+1
	end do

2000	continue

        IF (myrank==ncpu-1) THEN
	
c   Radial.

          CALL frange(0,fmhzmin0,fmhzmin)
	  
c   Run mineos for radial.

          CALL xmineos(1,0,eps,fmhzmin,fmhzmax,grav,TRIM(modelfile),
     1	    TRIM(string),frename)
	  
c   Run checkmode for radial.

          CALL xcheckmode(1,0,TRIM(frename),crit_err)
	  
        END IF
c
c
        CONTAINS
	
	SUBROUTINE check_dir(dirr)
	
          IMPLICIT NONE
          CHARACTER(LEN=*), INTENT(INOUT) :: dirr

 	  INTEGER :: k
	  k=LEN_TRIM(dirr)
	  IF (dirr(k:k) /= '/') THEN
	    k=k+1
	    dirr(k:k)='/'
  	  END IF 
	  
   	END SUBROUTINE check_dir
c
c
        SUBROUTINE frange(l,fmin0,fmin)
	
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: l
          DOUBLE PRECISION, INTENT(IN) :: fmin0
          DOUBLE PRECISION, INTENT(OUT):: fmin
          
c   fmhzbrack are rough estimates of minimum frequencies for ranges of 
c   l that can be used by mineos so that fundamental modes won't be missed.
          
          INTEGER, DIMENSION(11), PARAMETER ::
     &     lbrack=(/6,21,51,101,151,201,251,301,401,501,9999999/)
          DOUBLE PRECISION, DIMENSION(11), PARAMETER ::
     &     fmhzbrack=(/0.005,0.08,2.,5.,7.,10.,12.,15.,20.,25.,30./)

	  fmin=fmhzbrack(1)
	  if(l.eq.0) return

c   Set the minimum frequency. If fmhzmin<0, then use fmhzbrack according to l.

          if(fmin0.gt.0.d0) then
            fmin=fmin0
          else
            fmin=fmhzbrack(1)
            do i=1,10
              if((l.ge.lbrack(i)).and.(l.lt.lbrack(i+1))) fmin=fmhzbrack(i+1)
            end do
          endif
	  
        END SUBROUTINE frange
c
c
        SUBROUTINE xmineos(jcom,l,eps,fmhzmin,fmhzmax,grav,modelfile,
     1	  eigname,frename)
     
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: jcom, l
          DOUBLE PRECISION, INTENT(IN) :: fmhzmin, fmhzmax, grav
          CHARACTER(LEN=*), INTENT(IN) :: modelfile, eigname
          CHARACTER(LEN=128),INTENT(OUT) :: frename

          DOUBLE PRECISION :: eps
          CHARACTER(LEN=128) :: finp, funname, ftmp, mtype*3, id*8, string*256
          INTEGER :: ie, ib      
	  
          IF (jcom==1) THEN
            mtype='rad'
          ELSEIF (jcom==2) THEN
            mtype='tor'
          ELSEIF (jcom==3) THEN
            mtype='sph'
          END IF
	  
c    Make mineos input file

          id='_0000000'
          ie=LEN_TRIM(id)
          ib=ie-6
          WRITE(id(ib:ie),'(I7.7)') l
          finp='MINEOS_INPUT_FILES/mineos.inp'//TRIM(mtype)//TRIM(id)
          ftmp=TRIM(eigname)//'_'//TRIM(mtype)//TRIM(id)
          frename=TRIM(ftmp)//'.fre' 
          funname=TRIM(ftmp)//'.fun'
          open(unit=IOUT,file=finp,action='write')
          write(IOUT,'((A)/(A)/(A)/,I3,2(1X,E21.14)/,I5,1X,I5,2(1X,F8.3),A)')
     &      TRIM(modelfile), TRIM(frename), TRIM(funname),
     &      jcom,eps,grav,
     &      l,l,fmhzmin,fmhzmax,' 0'
          call flush(IOUT)
          close(IOUT)
	  
c   Run mineos.

          write(6,'(a,i5,a,f8.3,a,f8.3)')
     1      '  Running mineos for '//TRIM(mtype)//
     2	    ' for l=',l,'. fmhzmin=',fmhzmin,', fmhzmax=',fmhzmax
          call flush(6)
          string='./mineos < '//TRIM(finp)//
     &           ' > MINEOS_SCROUTPUT_FILES/mineos.out'//TRIM(mtype)//TRIM(id)
          call system(string)
          call system('sleep 5')
	  
        END SUBROUTINE xmineos
c
c
        SUBROUTINE xcheckmode(jcom,l,fname,crit_err)
	
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: jcom, l
          CHARACTER(LEN=*) :: fname 
 
          INTEGER, PARAMETER :: IOUT=35
          INTEGER :: ib, ie
          CHARACTER(LEN=128) :: finp, fout, fhead, ferr, mtype*3, string*256
	  double precision crit_err

c    Make checkmode input file 
       
          IF (jcom==1) THEN
            mtype='rad'
          ELSEIF (jcom==2) THEN
            mtype='tor'
          ELSEIF (jcom==3) THEN
            mtype='sph'
          END IF
          finp='CHECKMODE_INPUT_FILES/checkmode.'//TRIM(mtype)//'inp_0000000'
          fout='CHECKMODE_OUTPUT_FILES/checkmode.'//TRIM(mtype)//'out_0000000'
	  ferr='CHECKMODE_FREQERR_FILES/checkmode.'//TRIM(mtype)//'err_0000000'
          ie=LEN_TRIM(finp)          
          ib=ie-6
          WRITE(finp(ib:ie),'(I7.7)') l
          ie=LEN_TRIM(fout)
          ib=ie-6
          WRITE(fout(ib:ie),'(I7.7)') l
	  ie=LEN_TRIM(ferr)
	  ib=ie-6
	  WRITE(ferr(ib:ie),'(I7.7)') l
          open(unit=IOUT,file=finp,action='write')
          WRITE(IOUT,'((A)/(A)/(A)/(A)/(A)/(E11.4))') '4', TRIM(fname),
     1	    TRIM(fout),TRIM(ferr),'0',crit_err
          call flush(IOUT)
          close(IOUT)
	  
c    Run checkmode

          write(6,'(a,i3,a,i5)')
     1      '  Running checkmode for '//
     2	    TRIM(mtype)//' for l=',l
          call flush(6)
          write(fhead,"('checkmode.out',(A),'_',i7.7)") mtype, l
          string='./checkmode < '//TRIM(finp)//
     &          ' > MINEOS_SCROUTPUT_FILES/'//TRIM(fhead)
          call system(string)
	  
c    Wait to make sure the files have time to be written to the
c    shared file system

          call system('sleep 5')
	  
        END SUBROUTINE xcheckmode

	end
