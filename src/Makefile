include ../config.mk

BERTHAWLIB = bertha_wrapper.so

OBJ = \
      c_wrapper.o \
      cubefile.o \
      bertha_wrapper.o

all : $(BERTHAWLIB)

FFLAGS+= -fPIC -I${BERTHAROOT}/common -I${BERTHAROOT}/serial
CFLAGS+= -fPIC 

$(BERTHAWLIB): $(OBJ)
	cp ${BERTHAROOT}/lib/* ../lib
	$(FC) -shared $(OBJ) -o $(BERTHAWLIB) -L../lib -lbertha -lberthaserial $(LIBS)
	strip $(BERTHAWLIB)
	cp $(BERTHAWLIB) ../lib

clean:
	rm -f *.o *.mod *__genmod.f90 $(BERTHAWLIB) ../lib/$(BERTHAWLIB)
