include ../config.mk

BERTHAWLIB = bertha_wrapper.so

OBJ = bertha_wrapper.o \
      c_wrapper.o \
      cubefile.o

all : $(BERTHAWLIB)

FFLAGS+= -fPIC -I${BERTHAROOT}/common -I${BERTHAROOT}/serial
CFLAGS+= -fPIC 

$(BERTHAWLIB): $(OBJ)
	$(FC) -shared $(OBJ) -o $(BERTHAWLIB) -L${BERTHAROOT}/lib -lbertha -lberthaserial $(LIBS)
	cp $(BERTHAWLIB) ../lib

clean:
	rm -f *.o *.mod *__genmod.f90 $(BERTHAWLIB) ../lib/$(BERTHAWLIB)
