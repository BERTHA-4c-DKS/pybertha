# redo

# for BlueGene Q
FORBGQ=no

# debug yes or no
DEBUG=no

# profile yes or no
PROFILE=no

#is used only by serial
#use Intel compiler
USEINTEL=no

#LIBXC
LIBXC=no

BERTHAROOT=/home/redo/Sources/bertha_ng

###
## NO BLUEGENE
###
ifeq ($(FORBGQ),no)
  ifeq ($(USEINTEL),yes)
    FC = ifort
    CC = icc

    ifeq ($(PROFILE),yes)
      FFLAGS = -pg
      CFLAGS = -pg
      LINKFLAGS = -pg
    else
      FFLAGS =
      CFLAGS =
    endif

    # intel 
    #BLASLAPACK = -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
    #      $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a \
    #      -Wl,--end-group -lpthread -lm
    # BLASLAPACK = -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_rt.so  -Wl,--end-group -lpthread -lm -ldl
    # SCALDIR=/usr/local/SCALAPACK_intel_14_0_0_ompi2
    # BLACSDIR=/usr/local/BLACS_intel_14_0_0_ompi2
    # SCALAPACK=-L$(SCALDIR) -lscalapack
    # SCALAPACK=  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential \
    #     -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl
    #BLACS=$(BLACSDIR)/LIB/blacs_MPI-LINUX-0.a $(BLACSDIR)/LIB/blacsF77init_MPI-LINUX-0.a \
    # 	$(BLACSDIR)/LIB/blacs_MPI-LINUX-0.a
    
    #dynamic 
    #BLASLAPACK =   ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
    #	-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

    #static 
    BLASLAPACK =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a \
    	$(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread

    #BLASLAPACK = -llapack -lblas
    
    FFLAGS+= -I${MKLROOT}/include/intel64/lp64 -mkl=parallel 

    INCLUDE = 

    ifeq ($(DEBUG),yes)
      FFLAGS += -r8 -check all -check noarg_temp_created -traceback -warn all -O0 -g -132
      CFLAGS += -D_FILE_OFFSET_BITS=64 -O0 -g
    else
      #FFLAGS += -r8 -check all -check noarg_temp_created -traceback -warn all -O2 -132 
      FFLAGS += -r8 -warn all -O3 -132
      #FFLAGS += -C -O0 -r8 -warn all -132 -I./$(MODIR)
      CFLAGS += -D_FILE_OFFSET_BITS=64 -O3
    endif

    LIBS += $(BLASLAPACK)
  else
    FC = gfortran
    CC = gcc
    FOPT = 
    INCLUDE = 
    
    # gnu standard
    BLASLAPACK = -llapack -lblas
    SCALAPACK=-L/usr/lib64/openmpi/lib/ -lscalapack 
    BLACS=-L/usr/lib64/openmpi/lib/ -lmpiblacs

    # gnu custom
    #BLACSDIR=/home/mat/local/lib
    #BLASLAPACK = -L/home/mat/local/lib -ltmg -lreflapack -lrefblas 
    #SCALAPACK=-L/home/mat/local/lib -lscalapack 
    #BLACS=$(BLACSDIR)/libscalapack.a
    #SCALAPACK=-L/opt/share/scalapack_gnu -lscalapack

    ifeq ($(PROFILE),yes)
      FFLAGS = -pg
      CFLAGS = -pg
      LINKFLAGS = -pg
    else
      FFLAGS =
      CFLAGS =
    endif

    ifeq ($(DEBUG),yes)
      FFLAGS += -finit-local-zero -fdefault-double-8 -fdefault-real-8 -O0 -ffixed-line-length-132 -fbacktrace -ffpe-trap=zero,overflow,underflow -g -W -Wall -I./$(MODIR)
      CFLAGS += -D_FILE_OFFSET_BITS=64 -O0 -g -W -Wall
    else
      #FFLAGS += -finit-local-zero -fdefault-double-8 -fdefault-real-8 -O2 -I./$(MODIR) -W -Wall -ffixed-line-length-132
      FFLAGS +=  -fdefault-double-8 -fdefault-real-8 -O3 -I./$(MODIR) -W -Wall -ffixed-line-length-132
      CFLAGS += -D_FILE_OFFSET_BITS=64 -O3 -W -Wall
    endif

    LIBS += $(BLASLAPACK)
  endif
  
  CFLAGS += -W -Wall -DUSE_UNDER
else
  ## for the node
  #FC = bgxlf_r
  #CC = bgxlc_r
  # provo O2 per avere uno strict nella compilazion
  #FBASICF = -qrealsize=8 -O2 -qarch=qp -qtune=qp -qnostaticlink -qstackprotect -qhalt=w -w -qcheck -qflttrap
  #FFLAGS = -qfixed $(FBASICF)
  #CFLAGS = -O2 -qarch=qp -qtune=qp -qnostaticlink  -qstackprotect -qhalt=w -w -qcheck -qflttrap
  #LIBS = $(BGQLAPACK)

  # for the FE
  FC = xlf_r
  CC = xlc_r
  
  # blue gene 
  BGQLAPACK = -L${LAPACK_LIB} -llapack -L$(ESSL_LIB) -lesslbg 
  BGQFELAPACK = -L${LAPACK_LIB} -llapack -L$(ESSL_LIB) -lessl

  FBASICF = -qrealsize=8 -O2 -q64 -qnostaticlink -qstackprotect -qhalt=w -w -qcheck -qflttrap
  FFLAGS = -qfixed $(FBASICF)
  CFLAGS = -O2 -q64 -qnostaticlink -qstackprotect -qhalt=w -w -qcheck -qflttrap
  LIBS = $(BGQFELAPACK)
endif

FFLAGS += -I../common 

CFLAGS += -fPIC
FFLAGS += -fPIC

ifeq ($(LIBXC),yes)
  # Use libxc of a distribution DIRLIBXC to be set
  DIRLIBXC = /usr/lib/x86_64-linux-gnu
  #DIRLIBXC = /usr/local/libxc
  CFLAGS += -DLIBXC 
  FFLAGS += -DLIBXC 
  INCLUDE += -I$(DIRLIBXC)/include
  LIBS += -L$(DIRLIBXC) -lxc -lxcf90 
endif

ifeq ($(USEINTEL),yes)
  CFLAGS += -DUSEINTELCMP
  FFLAGS += -DUSEINTELCMP
endif

FFLAGS += -DDUMPFOCKMTX

MAKE = make

.SUFFIXES:

%.o:	%.c
	$(CC) $(CFLAGS) $(COPT) $(INCLUDE) -c $< 

%.o:	%.F
	$(FC) $(FFLAGS) $(FOPT) $(INCLUDE) -o $@ -c $< 

%.o:	%.f
	$(FC) $(FFLAGS) $(FOPT) $(INCLUDE) -o $@ -c $< 

%.o:	%.f90
	$(FC) $(FFLAGS) $(FOPT) $(INCLUDE) -o $@ -c $< 

%.o:	%.F90
	$(FC) $(FFLAGS) $(FOPT) $(INCLUDE) -o $@ -c $< 

ifeq ($(FORBGQ),no)
%.o:    %.f90
	$(FC) $(FFLAGS) $(FOPT) $(INCLUDE) -o $@ -c $< 
else
%.o : %.f90
	$(FC) $(FBASICF) -o $@ -c $<
endif
