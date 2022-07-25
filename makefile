#
OFLAGS = ifort -O4 -c -132 -zero -save -o
PFLAGS = ifort -O4    -132 -zero -save -o
MPIOFLAGS = mpiifort -O4 -c -132 -zero -save -o
MPIPFLAGS = mpiifort -O4    -132 -zero -save -o
#
#
OBJECTMP = modesum7_mpi.o
OBJECTP = mineos_mpi.o
OBJECTD = mineos_drv.o
OBJECTG = mineos_merge.o
OBJECT  = mineos.o	
OBJECTC = checkmode.o
OBJECTS = mineos_strip.o
OBJECTM = modesum7_nonmpi.o
OBJECTR = read_mineos.o
#
modesum7_mpi:    $(OBJECTMP)
	$(MPIPFLAGS) modesum7_mpi    $(OBJECTMP)
#
modesum7_mpi.o: modesum7_mpi.f
	$(MPIOFLAGS) modesum7_mpi.o  modesum7_mpi.f
#
mineos_mpi:    $(OBJECTP)
	$(MPIPFLAGS) mineos_mpi    $(OBJECTP)
#
mineos_mpi.o: mineos_mpi.f
	$(MPIOFLAGS) mineos_mpi.o  mineos_mpi.f
#
mineos_drv:     $(OBJECTD)
	$(PFLAGS) mineos_drv    $(OBJECTD)
#
mineos_drv.o: mineos_drv.f
	$(OFLAGS) mineos_drv.o  mineos_drv.f
#
mineos_merge:     $(OBJECTG)
	$(PFLAGS) mineos_merge    $(OBJECTG)
#
mineos_merge.o: mineos_merge.f
	$(OFLAGS) mineos_merge.o  mineos_merge.f
#
mineos:         $(OBJECT)
	$(PFLAGS) mineos      $(OBJECT)
#
mineos.o: mineos.f
	$(OFLAGS) mineos.o    mineos.f
#
checkmode:      $(OBJECTC)
	$(PFLAGS) checkmode      $(OBJECTC)
#
checkmode.o: checkmode.f
	$(OFLAGS) checkmode.o    checkmode.f
#
mineos_strip:   $(OBJECTS)
	$(PFLAGS) mineos_strip  $(OBJECTS)
#
mineos_strip.o: mineos_strip.f
	$(OFLAGS) mineos_strip.o mineos_strip.f
#
modesum7_nonmpi:   $(OBJECTM)
	$(PFLAGS) modesum7_nonmpi  $(OBJECTM)
#
modesum7_nonmpi.o: mineos_strip.f
	$(OFLAGS) modesum7_nonmpi.o modesum7_nonmpi.f
#
read_mineos:   $(OBJECTR)
	$(PFLAGS) read_mineos  $(OBJECTR)
#
read_mineos.o: read_mineos.f
	$(OFLAGS) read_mineos.o read_mineos.f
#
all_nonmpi: mineos_drv mineos_merge mineos checkmode mineos_strip modesum7_nonmpi read_mineos
#
all_mpi: modesum7_mpi mineos_mpi
#
all: modesum7_mpi mineos_mpi mineos_drv mineos_merge mineos checkmode mineos_strip modesum7_nonmpi read_mineos
#
clean:
	rm *.o *.mod modesum7_mpi mineos_mpi mineos_drv mineos_merge mineos checkmode mineos_strip modesum7_nonmpi read_mineos
