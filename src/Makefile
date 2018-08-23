include ../../config.mk

BERTHAWLIB = bertha_wrapper.so

OBJ = bertha_wrapper.o \
      c_wrapper.o \
      cubefile.o

all : $(BERTHAWLIB)

FFLAGS+= -fPIC -I$(ROOTDIR)/common -I$(ROOTDIR)/serial
CFLAGS+= -fPIC 

$(BERTHAWLIB): $(OBJ)
	$(FC) -shared $(OBJ) -o $(BERTHAWLIB) -L$(ROOTDIR)/lib -lbertha -lberthaserial $(LIBS)
	cp $(BERTHAWLIB) $(ROOTDIR)/lib

clean:
	rm -f *.o *.mod *__genmod.f90 $(BERTHAWLIB) $(ROOTDIR)/lib/$(BERTHAWLIB)
