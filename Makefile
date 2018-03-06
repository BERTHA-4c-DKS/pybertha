include ../config.mk

BERTHAWLIB = bertha_wrapper.so

OBJ = bertha_wrapper.o

all : $(BERTHAWLIB)

FFLAGS+= -fPIC -I../common -I../serial

$(BERTHAWLIB): $(OBJ)
	$(FC) -shared $(OBJ) -o $(BERTHAWLIB) -L../lib -lbertha -lberthaserial $(LIBS)

clean:
	rm -f *.o *.mod $(BERTHAWLIB)
