#########################
#       TARGET          #
#########################

TARGET= kriging

SUF90=f90
SUF77=f
.SUFFIXES: .f90 .f .o

GSL_prefix = /usr/local

#########################
#      COMPILATION      #
#########################

FAD	= mpif90
F90	= mpif90
F77	= mpif77

#FFLAGS  = -mp -convert big_endian -fpconstant -zero -c

FFLAGS  = -O3 -fopenmp  -ffree-line-length-10000 #-fbounds-check -fdefault-integer-8

#LIBS    = -L/usr/local/lib -L/Softwareinstall/gsl-1.15/.libs/lgsl -L/Softwareinstall/gsl-1.15/.libs/lgslcblas -lm #-L/Softwareinstall/gsl-1.15/.libs/lmir
LFLAGS =  -L$(GSL_prefix)/lib -lgsl -lgslcblas -lm
LIBS = -ldl -lstdc++

%.mod : %.o
	@if [! -f $@ ]; then \
	rm $< \
	$(MAKE) $< \
	fi

SRCS =  dimKrig.o Timer.o main.o functions.o \
        latin.o mpi.o optimize.o threebarcost.o \
        read_set.o Dutch.o Dutchgeninterp.o\
        read_sample.o check_sample.o \
        make_krig.o reduce_data.o tool.o eva_sample.o \
	matrixR.o petitR.o \
        search_krig.o meta.o diff.o output_des.o \
        correct.o LUroutines.o\
        higher.o variance.o\
        indirect.o trust.o trustool.o \
        rank.o vrange.o nieder.o \
        update.o DynamicPointSelection.o\
        ludim.o lusol2.o halton.o \
        bfgs.o bfgs_routines.o eva_bfgs.o \
        dimGA.o ga.o gatool.o \
        scf.o post.o monaco.o make_sample.o\
        scf_df.o scf_df_df.o scf_df_df_db.o \
        scf_db.o scf_db_db.o scf_db_db_df.o \
	scf_db_df.o scf_df_db.o scf_df_df_db_db.o \
	hammersley.o sobol.o faure.o rmsebound.o

OBJS =  ${SRCS:.$(SUF)=.o}

all:  $(TARGET)

$(TARGET): $(OBJS) libmir.a tapenade.a
	$(F90) $(FFLAGS) -o $(TARGET) $(OBJS) libmir.a tapenade.a $(LFLAGS) -Wl,-rpath=.
	@echo " ----------- ${TARGET} created ----------- "
.$(SUF90).o:
	$(F90) $(FFLAGS) -c $<
.$(SUF77).o:
	$(F77) $(FFLAGS) -c $<

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clean:
	rm -f $(OBJS) *.L *.msg *.???~ *~
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

