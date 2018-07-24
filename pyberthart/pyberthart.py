import argparse
import os.path
import ctypes
import numpy
import sys
import re

import scipy.linalg as scila
from numpy.linalg import eigvalsh
from scipy.linalg import eigh

sys.path.insert(0, '../src/')
import berthamod
import rtutil

listpulses = ""
for key in rtutil.funcswitcher:
    listpulses += key
    listpulses += " "

parser = argparse.ArgumentParser()

parser.add_argument("-f","--inputfile", help="Specify BERTHA input file (default: input.inp)", required=False, 
        type=str, default="input.inp")
parser.add_argument("-t","--fittfile", help="Specify BERTHA fitting input file (default: fitt2.inp)", required=False, 
        type=str, default="fitt2.inp")
parser.add_argument("-c","--fitcoefffile", help="Specify BERTHA fitcoeff output file (default: fitcoeff.txt)",
        required=False, type=str, default="fitcoeff.txt")
parser.add_argument("-e","--vctfile", help="Specify BERTHA vct output file (default: vct.txt)", required=False, 
        type=str, default="vct.txt")
parser.add_argument("-p","--ovapfile", help="Specify BERTHA ovap output file (default: ovap.txt)", required=False, 
        type=str, default="ovap.txt")
parser.add_argument("-s", "--dumpfiles", help="Dumpfile on, default is off", required=False,
        default=False, action="store_true")
parser.add_argument("-m", "--dt", help="Specify dt to be used (default: 0.1)", required=False,
        default=0.1, type=numpy.float64)
parser.add_argument("-T", "--totaltime", help="Specify total time )deaful: 1.0)", required=False,
        default=1.0, type=numpy.float64)
parser.add_argument("--pulse", help="Specify the pulse to use [" + listpulses + "] default kick", required=False, 
        type=str, default="kick")
parser.add_argument("--pulseFmax", help="Specify the pulse Fmax value (default: 0.0001)", 
        default=0.0001, type=numpy.float64)
parser.add_argument("--pulsew", help="Specify the pulse w value if needed (default: 0.0)", 
        default=0.0, type=numpy.float64)
parser.add_argument("-d", "--debug", help="Debug on, prints debug info to debug_info.txt", required=False, 
        default=False, action="store_true")
parser.add_argument("-v", "--verbosity", help="Verbosity level 0 = minim, -1 = print iteration info, " + 
        "1 = maximum (defaul -1)", required=False, default=-1, type=int)
parser.add_argument("--iterations", help="Use iteration number instead of progressbar",
        required=False, default=False, action="store_true")

args = parser.parse_args()

bertha = berthamod.pybertha("../../lib/bertha_wrapper.so")

fittcoefffname = args.fitcoefffile
vctfilename = args.vctfile
ovapfilename = args.ovapfile

fnameinput = args.inputfile
if not os.path.isfile(fnameinput):
    print "File ", fnameinput, " does not exist"
    exit(1)

fittfname = args.fittfile
if not os.path.isfile(fittfname):
    print "File ", fittfname , " does not exist"
    exit(1)

verbosity = args.verbosity
dumpfiles = int(args.dumpfiles)

bertha.set_fittcoefffname(fittcoefffname)
bertha.set_ovapfilename(ovapfilename)
bertha.set_vctfilename(vctfilename)
bertha.set_fnameinput(fnameinput)
bertha.set_fittfname(fittfname)

bertha.set_verbosity(verbosity)
bertha.set_dumpfiles(dumpfiles)

bertha.init()

ndim = bertha.get_ndim()
nshift = bertha.get_nshift()
nocc = bertha.get_nocc()
sfact = bertha.get_sfact()
nopen = bertha.get_nopen()

print "Verbosity       : ", verbosity
print "Dumpfiles       : ", dumpfiles
print ""
print "Matrix dimension: ", ndim
print "            nocc: ", nocc
print "          nshift: ", nshift
print "           nopen: ", nopen
print "     level shift: ", sfact
print ""

ovapm, eigem, fockm, eigen = bertha.run()
if (fockm is None) or (eigen is None) or (fockm is None) \
        or (eigen is None):
    print "Error in bertha run"
    exit(-1)

print ""
print "Final results "
sum=0.0
for i in range(nocc+nopen):
    print "eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact)
    sum+=eigen[i+nshift]-sfact
print "      lumo       %20.8f"%(eigen[i+nshift+1])
print "      Sum of eigen : ", sum
erep = bertha.get_erep()
etotal = bertha.get_etotal()

print ""
print "total electronic energy  = %20.8f"%(etotal-(sfact*nocc))
print "nuclear repulsion energy = %20.8f"%(erep)
print "total energy             = %20.8f"%(etotal+erep-(sfact*nocc))
print " "

bertha.realtime_init()

print "Start RT"

debug = args.debug
dt = args.dt
t_int = args.totaltime
niter = int(t_int/dt)

print "Debug: ", debug
print "dt : ", dt
print "Total time  : ", t_int
print "Number of iterations: ", niter

ene_list = []
dip_list = []
imp_list = []
Enuc_list = []

C = eigem
D_0 = numpy.zeros((ndim,ndim), dtype=numpy.complex128)
for num in range(nocc):
    D_0[num+nshift,num+nshift]=1.0+0.0j

fo = sys.stderr
if debug:
    fo = open("debug_info.txt", "w")

#print type(eigem)
#C_inv used to backtransform D(AO)

try: 
    C_inv = numpy.linalg.inv(eigem)
except LinAlgError:
    print "Error in numpy.linalg.inv of eigem" 
    exit(1)

if debug:
  test = numpy.matmul(C_inv, eigem)
  fo.write("Check if atol is 1.e-14 for inversion of C: %s\n"% \
        numpy.allclose(numpy.eye((ndim),dtype=numpy.complex128), \
        test,atol=1.e-14))

if debug:
  diff = test - numpy.eye((ndim),dtype=numpy.complex128)
  mdiff = numpy.max(diff)
  fo.write("  maxdiff is: %.12e %.12e\n"%(mdiff.real, mdiff.imag))

if debug:
  test = numpy.matmul(numpy.conjugate(C.T),numpy.matmul(ovapm,C))
  fo.write("Check orthonormal orbitals (atol = 1.e-14): %s\n"% \
        numpy.allclose(numpy.eye((ndim),dtype=numpy.complex128), \
        test, atol=1.e-14))
  diff = test - numpy.eye((ndim),dtype=numpy.complex128)
  mdiff = numpy.max(diff)
  fo.write("  maxdiff is: %.12e %.12e\n"%(mdiff.real, mdiff.imag))

#build density in ao basis

occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)
iocc = 0

for i in range(ndim):
    if i >= nshift and iocc < nocc:
        for j in range(ndim):
            occeigv[j, iocc] = eigem[j, i]
        iocc = iocc + 1

Da = numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))

if debug:
  #check trace(S Da)
  trace_ds = numpy.trace(numpy.matmul(Da,ovapm))
  trace_dsfock = numpy.trace(numpy.matmul(Da,fockm))
  fo.write("Density matrix trace at  t0: %.12e %.12e \n"%(trace_ds.real,trace_ds.imag))
  fo.write("Trace of fock*density at t0: %.12e %.12e \n"%(trace_dsfock.real, trace_dsfock.imag))

direction = 2
normalise = 1
dipz_mat = bertha.get_realtime_dipolematrix (direction, normalise)

if debug:
  fockmh = numpy.conjugate(fockm.T)
  diff_fockmh = fockm-fockmh
  mdiff = numpy.max(diff_fockmh)
  fo.write("Check max diff fockm-fockmh: %.12e %.12e\n"%\
        (mdiff.real, mdiff.imag))
  fo.write("Fockm (t=0) is hermitian: %s \n"%numpy.allclose(fockm,fockmh,atol=1.e-15))

print "Start first mo_fock_mid_forwd_eval "

fock_mid_init = rtutil.mo_fock_mid_forwd_eval(bertha,Da,fockm,0,numpy.float_(dt),\
	dipz_mat,C,C_inv,ovapm,ndim, debug, fo, args.pulse, args.pulseFmax, args.pulsew)

if (fock_mid_init is None):
    print "Error accurs in mo_fock_mid_forwd_eval"
    exit(1)

if debug:
  fock_mid_h=numpy.conjugate(fock_mid_init.T)
  diff_fock_mid_h=fock_mid_init-fock_mid_h
  fo.write("Max diff fock_mid_init-fock_mid_h"%(numpy.max(diff_fock_mid_h)))
  fo.write('Fockm (t=1/2) is hermitian: %s\n'% 
        numpy.allclose(fock_mid_init,fock_mid_h,atol=1.e-14))

fockp_mid_init=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(fock_mid_init,C))
u=rtutil.exp_opmat(fockp_mid_init,numpy.float_(dt))
#u=rtutil.exp_opmat(fockp_mid_init,numpy.float_(dt))
#u=scila.expm(-1.j*fockp_mid_init*dt)
temp=numpy.matmul(D_0,numpy.conjugate(u.T))
Dp_t1=numpy.matmul(u,temp)

#check u if unitary
if debug:
  test_u=numpy.matmul(u,numpy.conjugate(u.T))
  fo.write('U is unitary : %s' % 
        numpy.allclose(test_u,numpy.eye(u.shape[0]),atol=1.e-14))

#backtrasform Dp_t1
D_t1=numpy.matmul(C,numpy.matmul(Dp_t1,numpy.conjugate(C.T)))

if debug:
  diff = D_t1 - Da 
  mdiff = numpy.max(diff)
  fo.write("Max diff density: %.12e %.12e \n"%(mdiff.real, mdiff.imag))

dip_list.append(numpy.trace(numpy.matmul(Da,dipz_mat)))
dip_list.append(numpy.trace(numpy.matmul(D_t1,dipz_mat)))
if debug:
  fo.write("Dipole: %.12e\n"%(numpy.trace(numpy.matmul(Da,dipz_mat))).real)
  fo.write("Dipole: %.12e\n"%(numpy.trace(numpy.matmul(D_t1,dipz_mat))).real)

if debug:
  tfock = numpy.trace(numpy.matmul(D_t1,fockm))
  fo.write("Trace fock*density t1: %.12e, %.12e\n"%(tfock.real, tfock.imag))
  trace_ds=numpy.trace(numpy.matmul(D_t1,ovapm))
  fo.write(" Traceds: %.12e %.12ei\n" % (trace_ds.real,trace_ds.imag))

Ndip_z=0.0
#estrarre le 3 componenti del dipolo nucleare

D_ti=D_t1
Dp_ti=Dp_t1
#aggiungere repulsione nucleare
#Enuc_list.append(-func_t0*Ndip_z+Nuc_rep) #just in case of non-zero nuclear dipole
fockm_ti=bertha.get_realtime_fock(D_ti.T)
ene_list.append(numpy.trace(numpy.matmul(Da,fockm)))
ene_list.append(numpy.trace(numpy.matmul(D_ti,fockm_ti)))

print "Starting iterations ..."
print ""

fock_mid_backwd=numpy.copy(fock_mid_init)
for j in range(1,niter):

    if (args.iterations):
        print "Iteration ", j, " of ", niter-1
    else:
        rtutil.progress_bar(j, niter-1)

    fock_mid_tmp = rtutil.mo_fock_mid_forwd_eval(bertha,numpy.copy(D_ti), \
            fock_mid_backwd,j,numpy.float_(dt),dipz_mat,C,C_inv,ovapm,ndim,\
                    debug, fo, args.pulse, args.pulseFmax, args.pulsew)
    
    if (fock_mid_tmp is None):
        print "Error accurs in mo_fock_mid_forwd_eval"
        exit(1)

    if debug:
      fo.write('%.8f\n' % numpy.trace(numpy.matmul(ovapm,D_ti)).real)
    Ah = numpy.conjugate(fock_mid_tmp.T)
    if debug:
      fo.write('Fock_mid hermitian: %s\n' % numpy.allclose(fock_mid_tmp,Ah,atol=1.e-14))
    # transform fock_mid_init in MO basis
    fockp_mid_tmp = numpy.matmul(numpy.conjugate(C.T),numpy.matmul(fock_mid_tmp,C))
    u = rtutil.exp_opmat(numpy.copy(fockp_mid_tmp),numpy.float_(dt))
    # u=scila.expm(-1.0j*fockp_mid_tmp*dt)
    # check u is unitary
    test_u = numpy.matmul(u,numpy.conjugate(u.T))
    
    if (not numpy.allclose(numpy.eye(u.shape[0]),test_u,atol=1.e-14)):
        print("U is not unitary\n")
        exit(1)

    #evolve the density in orthonormal basis
    #check the trace of density to evolve
    if debug:
      fo.write('trace of density to evolve: %.8f\n' % numpy.trace(Dp_ti).real)
    temp=numpy.matmul(Dp_ti,numpy.conjugate(u.T))
    Dp_ti_dt = numpy.matmul(u,temp)
    #backtransform Dp_ti_dt
    D_ti_dt=numpy.matmul(C,numpy.matmul(Dp_ti_dt,numpy.conjugate(C.T)))
    if debug:
      fo.write('  Trace of D_ti_dt %.8f\n' % numpy.trace(Dp_ti_dt).real)
    #dipole expectation for D_ti_dt
    dip_list.append(numpy.trace(numpy.matmul(dipz_mat,D_ti_dt)))
    if debug:
      fo.write("Dipole: %.12e\n"%(numpy.trace(numpy.matmul(dipz_mat,D_ti_dt)).real))

    #Energy expectation value at t = t_i_dt 
    fockm_ti_dt=bertha.get_realtime_fock(D_ti_dt.T)

    ene_list.append(numpy.trace(numpy.matmul(D_ti_dt,fockm_ti_dt)))
    
    # update D_ti and Dp_ti for the next step
    # message for debug
    # fo.write('here I update the matrices Dp_ti and D_ti\n')
    D_ti = numpy.copy(D_ti_dt)
    Dp_ti = numpy.copy(Dp_ti_dt)
  
    if debug:
      fo.write('  Trace of Dp_ti %.8f\n' % numpy.trace(Dp_ti).real)
      fo.write('  Trace of ovapm * D_ti  %.8f\n' % numpy.trace(numpy.matmul(ovapm,D_ti)).real)
      fo.flush()

    #update fock_mid_backwd for the next step
    fock_mid_backwd=numpy.copy(fock_mid_tmp)

print ""
print ""
print "Done"

if debug:
    fo.close()

t_point = numpy.linspace(0.0, niter*dt, niter+1)
numpy.savetxt('dipole.txt',numpy.c_[t_point.real,numpy.array(dip_list).real], fmt='%.12e')
numpy.savetxt('ene.txt',numpy.c_[t_point.real,numpy.array(ene_list).real], fmt='%.12e')

bertha.finalize()
