include ../config.mk

BERTHAWLIB = bertha_wrapper.so

OBJ = \
      c_wrapper.o \
      cubefile.o \
      bertha_wrapper.o

all : $(BERTHAWLIB)

FFLAGS+= -I${BERTHAROOT}/common -I${BERTHAROOT}/serial

$(BERTHAWLIB): $(OBJ)
	cp ${BERTHAROOT}/lib/* ../lib
	$(FC) -shared $(LINKFLAGS) $(OBJ) -o $(BERTHAWLIB) -L../lib -lbertha -lberthaserial $(LIBS)
	strip $(BERTHAWLIB)
	cp $(BERTHAWLIB) ../lib

clean:
	rm -f *.o *.mod *__genmod.f90 $(BERTHAWLIB) ../lib/$(BERTHAWLIB)
