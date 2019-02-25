# redo

include ./config.mk

all:
	make -C ./src

clean:
	make -C ./src clean
