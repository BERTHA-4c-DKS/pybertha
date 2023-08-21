include ../config.mk

BERTHAWLIB = bertha_wrapper.so

#version
MAJOR=0
MINOR=11
PATCH=0

CFLAGS+=-DMAJOR=$(MAJOR) -DMINOR=$(MINOR) -DPATCH=$(PATCH)

COMMONDIR = ${BERTHAROOT}/common
COMMONMODDIR = ${BERTHAROOT}/common/MOD
SERIALMODIR = ${BERTHAROOT}/serial/MOD
SERIALDIR = ${BERTHAROOT}/serial

OBJ = \
	$(COMMONMODDIR)/rlight.o \
	$(COMMONMODDIR)/constants.o \
	$(COMMONMODDIR)/limits.o \
	$(COMMONMODDIR)/shiftr.o \
	$(COMMONMODDIR)/spec.o \
	$(COMMONMODDIR)/opensh.o \
	$(COMMONMODDIR)/clabel.o \
	$(COMMONMODDIR)/energy.o \
	$(COMMONMODDIR)/timetot.o \
	$(COMMONMODDIR)/looptime.o \
	$(COMMONMODDIR)/gammam.o \
	$(COMMONMODDIR)/funfm.o \
	$(COMMONMODDIR)/fittingdata2.o \
	$(COMMONMODDIR)/fittingdata3.o \
	$(COMMONMODDIR)/mod_f_matrix.o \
	$(COMMONMODDIR)/labels.o \
	$(COMMONMODDIR)/dtrngle.o \
	$(COMMONMODDIR)/locgeo.o \
	$(COMMONMODDIR)/indsys.o \
	$(COMMONMODDIR)/cont.o \
	$(COMMONMODDIR)/nuc.o \
	$(COMMONMODDIR)/temp.o \
	$(COMMONMODDIR)/gen.o \
	$(COMMONMODDIR)/rnorm.o \
	$(COMMONMODDIR)/periodictable.o \
	$(COMMONMODDIR)/time.o \
	$(COMMONMODDIR)/gammat.o \
	$(COMMONMODDIR)/gammaf2.o \
	$(COMMONMODDIR)/s3j.o \
	$(COMMONMODDIR)/ang.o \
	$(COMMONMODDIR)/ints.o \
	$(COMMONMODDIR)/consts.o \
	$(COMMONMODDIR)/klc.o \
	$(COMMONMODDIR)/indx.o \
	$(COMMONMODDIR)/ijc.o \
	$(COMMONMODDIR)/libxc_funcs.o \
	$(COMMONMODDIR)/libxc_master.o \
        $(COMMONDIR)/cmatrix.o \
        $(COMMONDIR)/rmake.o \
        $(COMMONDIR)/stepn.o \
        $(COMMONDIR)/ess0.o \
        $(COMMONDIR)/utilf.o \
        $(COMMONDIR)/meminfo.o \
        $(COMMONDIR)/readstatm.o \
        $(COMMONDIR)/tempo.o \
        $(COMMONDIR)/util.o \
        $(COMMONDIR)/angmat.o \
        $(COMMONDIR)/aufbau.o \
        $(COMMONDIR)/vrs.o \
        $(COMMONDIR)/chi.o \
        $(COMMONDIR)/chigrad_last.o \
        $(COMMONDIR)/c_rks_p86.o \
        $(COMMONDIR)/cutoff.o \
        $(COMMONDIR)/ell0.o \
        $(COMMONDIR)/rnlarge.o \
        $(COMMONDIR)/esgtf.o \
        $(COMMONDIR)/emakess.o \
        $(COMMONDIR)/emakell.o \
        $(COMMONDIR)/factrl.o \
        $(COMMONDIR)/factrl3.o \
        $(COMMONDIR)/factrl4.o \
        $(COMMONDIR)/gen_oh.o \
        $(COMMONDIR)/gfinit.o \
        $(COMMONDIR)/gradhgfamplitude.o \
        $(COMMONDIR)/gradhgflap.o \
        $(COMMONDIR)/hgfamplitude.o \
        $(COMMONDIR)/hgfkinetic.o \
        $(COMMONDIR)/hgflap.o \
        $(COMMONDIR)/hgfoverlap.o \
        $(COMMONDIR)/hgtf.o \
        $(COMMONDIR)/indtuv.o \
        $(COMMONDIR)/jfitenergy.o \
        $(COMMONDIR)/klset.o \
        $(COMMONDIR)/lcase.o \
        $(COMMONDIR)/lebedev.o \
        $(COMMONDIR)/oemakell.o \
        $(COMMONDIR)/oemakess.o \
        $(COMMONDIR)/overlapmake.o \
        $(COMMONDIR)/p_block.o \
        $(COMMONDIR)/plegdr.o \
        $(COMMONDIR)/rnorma.o \
        $(COMMONDIR)/rnormf.o \
        $(COMMONDIR)/rnsmall.o \
        $(COMMONDIR)/s_block.o \
        $(COMMONDIR)/slmake.o \
        $(COMMONDIR)/ssmake.o \
        $(COMMONDIR)/stepl.o \
        $(COMMONDIR)/steplm.o \
        $(COMMONDIR)/xc_b88pw86.o \
        $(COMMONDIR)/xc_rks_hcth.o \
        $(COMMONDIR)/ecart.o \
        $(COMMONDIR)/errmsg.o \
	$(SERIALMODIR)/dcoeff.o \
	$(SERIALMODIR)/coeff.o \
	$(SERIALMODIR)/fmat.o \
	$(SERIALMODIR)/denhgf.o \
	$(SERIALMODIR)/realtime.o \
	$(SERIALMODIR)/openaccinterfaces.o \
        $(SERIALDIR)/diagonalizer.o \
        $(SERIALDIR)/initfitting.o \
        $(SERIALDIR)/flabel.o \
        $(SERIALDIR)/dens.o \
        $(SERIALDIR)/hden.o \
        $(SERIALDIR)/gdensh.o \
        $(SERIALDIR)/densh.o \
        $(SERIALDIR)/v_calculation.o \
        $(SERIALDIR)/factlu.o \
        $(SERIALDIR)/jmatfitt.o \
        $(SERIALDIR)/dftinit_fit.o \
        $(SERIALDIR)/dftinit_fit_allext.o \
        $(SERIALDIR)/dftinit_fit_minmem.o \
        $(SERIALDIR)/dftinit_fit_quadpregrid.o \
        $(SERIALDIR)/rjmat.o \
        $(SERIALDIR)/dftint.o \
        $(SERIALDIR)/oneel.o \
        $(SERIALDIR)/shftlv.o \
        $(SERIALDIR)/atomic.o \
        $(SERIALDIR)/smatrix.o \
        $(SERIALDIR)/pmatrix.o \
        $(SERIALDIR)/vtvec.o \
        $(SERIALDIR)/vavec.o \
        $(SERIALDIR)/jmatfittt.o \
        $(SERIALDIR)/jmatfitta.o \
        $(SERIALDIR)/gridset.o \
        $(SERIALDIR)/rlebdv.o \
        $(SERIALDIR)/fit_excpot.o \
        $(SERIALDIR)/gradfu.o \
        $(SERIALDIR)/lda1.o \
        $(SERIALDIR)/rel1e.o \
        $(SERIALDIR)/fock.o \
        $(SERIALDIR)/ilmake.o \
        $(SERIALDIR)/ismake.o \
        $(SERIALDIR)/fit_dens.o \
        $(SERIALDIR)/gradfit_dens.o \
        $(SERIALDIR)/xc_blyp.o \
        $(SERIALDIR)/rint.o \
        $(SERIALDIR)/klinit.o \
        $(SERIALDIR)/e2c1.o \
        $(SERIALDIR)/x_rks_b88.o \
        $(SERIALDIR)/getmem.o \
        $(SERIALDIR)/density_vxc.o \
        $(SERIALDIR)/dftint_Exc.o \
        $(SERIALDIR)/check_Exc.o \
        $(SERIALDIR)/density.o \
        $(SERIALDIR)/rnormnr.o \
        $(SERIALDIR)/gencubefile.o \
        $(SERIALDIR)/dipolematrix.o \
        $(SERIALDIR)/expmat.o \
        $(SERIALDIR)/epscalc.o \
        $(SERIALDIR)/fit_basis_xyz.o \
        $(SERIALDIR)/onel_pot_from_grid.o \

ifeq ($(USEOPENMP),yes)
        OBJ += \
	  $(SERIALDIR)/jmatfitttopenmp.o
endif

ifeq ($(USECUDA),yes)
        OBJ += \
	  $(SERIALDIR)/jmatfitttopenacc.o
endif

OBJ += \
      c_wrapper.o \
      cubefile.o \
      bertha_wrapper.o

all : $(BERTHAWLIB)

FFLAGS+= -I${BERTHAROOT}/common -I${BERTHAROOT}/serial

$(BERTHAWLIB): $(OBJ)
	$(FC) -shared $(LINKFLAGS) $(OBJ) -o $(BERTHAWLIB) $(LIBS)
	cp $(BERTHAWLIB) ../lib

clean:
	rm -f *.o *.mod *__genmod.f90 $(BERTHAWLIB) ../lib/$(BERTHAWLIB) rm -f $(SERIALDIR)/*.o $(SERIALDIR)/*.mod $(SERIALDIR)/*__genmod.f90 \
	  $(SERIALMODIR)/*.o $(SERIALMODIR)/*.mod $(SERIALMODIR)/*__genmod.f90  $(COMMONMODDIR)/*.o $(COMMONMODDIR)/*.mod $(COMMONMODDIR)/*__genmod.f90 \
	  $(COMMONDIR)/*.o $(COMMONDIR)/*.mod $(COMMONDIR)/*__genmod.f90
