########################
#       TARGET          #
#########################

TARGET= krigingestimate.a

SUF90=f90
SUF77=f
.SUFFIXES: .f90 .f .o

#########################
#      COMPILATION      #
#########################

FAD	= mpif90
F90	= mpif90
F77	= mpif77

#FFLAGS  = -mp -convert big_endian -fpconstant -zero -c
FFLAGS  = -r8 -O4 -openmp 

LIBS    = 

SRCS =  dimKrig.o main.o functions.o\
        threebarcost.o latin.o mpi.o\
        read_set.o Dutch.o Dutchgeninterp.o\
        read_sample.o check_sample.o \
        make_krig.o reduce_data.o tool.o eva_sample.o \
	matrixR.o petitR.o \
        search_krig.o meta.o diff.o output_des.o \
        correct.o LUroutines.o\
        higher.o variance.o\
        indirect.o trust.o trustool.o \
        rank.o vrange.o extraroutine.o\
        update.o DynamicPointSelection.o\
        ludim.o lusol2.o \
        bfgs.o bfgs_routines.o eva_bfgs.o \
        dimGA.o ga.o gatool.o \
        post.o monaco.o make_sample.o
        
OBJS =  ${SRCS:.$(SUF)=.o}

all:  $(TARGET)

$(TARGET): $(OBJS) 
	   ar rvs $@ $(OBJS)
	@echo " ----------- ${TARGET} created ----------- "
.$(SUF90).o:
	$(F90) $(FFLAGS) -c $<
.$(SUF77).o:
	$(F77) $(FFLAGS) -c $<

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clean:
	rm -f $(OBJS) *.L *.msg *.???~ *~
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

