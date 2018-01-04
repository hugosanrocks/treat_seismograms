include $(HOME)/.makefilerc
#
EXEC=../bin/
LIBS=
INC=

all : arrtime 

.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(F90_OPTION) $(INC) -c $*.f90

.SUFFIXES: .f .o
.f.o:
	$(F77) $(F77_OPTION) $(INC) -c $*.f

.SUFFIXES: .c .o
.c.o:
	$(CC) $(C_OPTION) $(INC) -c $*.c

clean:
	$(RM) *.o *~ *.mod \#*

#
SUB = sub_hs.o time_3d.o
#

arrtime : $(SUB) arrtime.o
	$(F90) $(F90_OPTION) arrtime.o $(SUB) $(LIBS) -o $(EXEC)arrtime

