import ctypes
import numpy
import sys
import re
import os.path

import rtutil
import scipy.linalg as scila
from numpy.linalg import eigvalsh
from scipy.linalg import eigh

sys.path.insert(0, '../src/')
import berthamod

bertha = berthamod.pybertha("../../lib/bertha_wrapper.so")

fittcoefffname = "fitcoeff.txt"
vctfilename = "vct.txt" 
ovapfilename = "ovap.txt"
fnameinput = "input.inp"
fittfname = "fitt2.inp"

verbosity = -1
dumpfiles = 0


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
print "fock[0,10] = ", fockm[0,10]
if (fockm is None) or (eigen is None) or (fockm is None) \
        or (eigen is None):
    print "Error in bertha run"
    exit(-1)

print ""
print "Final results "
for i in range(nocc+nopen):
    print "eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact)
    
print "      lumo       %20.8f"%(eigen[i+nshift+1])

erep = bertha.get_erep()
etotal = bertha.get_etotal()

print ""
print "total electronic energy  = %20.8f"%(etotal-(sfact*nocc))
print "nuclear repulsion energy = %20.8f"%(erep)
print "total energy             = %20.8f"%(etotal+erep-(sfact*nocc))

bertha.realtime_init()

D_0=numpy.zeros((ndim,ndim),dtype=numpy.complex128)
for num in range(nocc):
    D_0[num,num]=1.0+0.0j
dt=0.1
#containers
ene_list = []
dip_list = []
imp_list = []
Enuc_list = []
C = eigem
print type(eigem)
#C_inv used to backtransform D(AO)
try: 
    C_inv = numpy.linalg.inv(eigem)
except LinAlgError:
    print "error" 
test=numpy.matmul(C_inv,eigem)
test1=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(ovapm,C))
#print test1
print numpy.allclose(numpy.eye((ndim),dtype=numpy.complex128),test1)
#build density in ao basis
occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)

iocc = 0

for i in range(ndim):
    if i >= nshift and iocc < nocc:
        for j in range(ndim):
            occeigv[j, iocc] = eigem[j, i]
        iocc = iocc + 1

Da=numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))
#check trace(S Da)
print "traccia prima della propagazione"
trace_ds=numpy.trace(numpy.matmul(Da,ovapm))
print('%.14f %.14fi\n' % (trace_ds.real,trace_ds.imag))
direction = 2
normalise = 1

dipz_mat = bertha.get_realtime_dipolematrix(direction, normalise)

fockmh=numpy.conjugate(fockm.T)
diff_fockmh=fockm-fockmh
print "max diff fockm-fockmh"
print numpy.max(diff_fockmh)
print('Fockm (t=0) is hermitian: %s\n' % numpy.allclose(fockm,fockmh,atol=1.e-15))


fock_mid_init = rtutil.mo_fock_mid_forwd_eval(bertha,Da,fockm,0,numpy.float_(dt),\
	dipz_mat,C,C_inv,ovapm,ndim)

fockmh=numpy.conjugate(fock_mid_init.T)
diff_fockmh=fock_mid_init-fockmh
print "max diff fock_mid_init-fockmh"
print numpy.max(diff_fockmh)
print('Fockm (t=1/2) is hermitian: %s\n' % numpy.allclose(fock_mid_init,fockmh,atol=1.e-15))
diff =fock_mid_init-fockm
print "massimo valore di fock-fock_1/2"
print numpy.max(diff)
Ah=numpy.conjugate(fock_mid_init.T)
print('Fock_mid hermitian: %s\n' % numpy.allclose(fock_mid_init,Ah,atol=1.e-15))

fockmh=numpy.conjugate(fockm.T)
diff_fockmh=fockm-fockmh
print "max diff fockm-fockmh"
print numpy.max(diff_fockmh)
print('Fockm (t=0) is hermitian: %s\n' % numpy.allclose(fockm,fockmh,atol=1.e-15))
fockp_mid_init=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(fock_mid_init,C))
fockmp=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(fockm,C))
print fockmp[0:3,0:3]
u=rtutil.exp_opmat(fockmp,numpy.float_(dt))
#u=rtutil.exp_opmat(fockp_mid_init,numpy.float_(dt))
#u=scila.expm(-1.j*fockp_mid_init*dt)
temp=numpy.matmul(D_0,numpy.conjugate(u.T))
Dp_t1=numpy.matmul(u,temp)
#check u if unitary
test_u=numpy.matmul(u,numpy.conjugate(u.T))
print('U is unitary :%s' % numpy.allclose(test_u,numpy.eye(u.shape[0]),atol=1.e-15))
#backtrasform Dp_t1
D_t1=numpy.matmul(C,numpy.matmul(Dp_t1,numpy.conjugate(C.T)))
diff = D_t1 - Da 
print "diff density: ", numpy.max(diff)
dip_list.append(numpy.trace(numpy.matmul(Da,dipz_mat)))
dip_list.append(numpy.trace(numpy.matmul(D_t1,dipz_mat)))

print "trace fock*density t0: ", numpy.trace(numpy.matmul(Da,fockm))
print "trace fock*density t1: ", numpy.trace(numpy.matmul(D_t1,fockm))

trace_ds=numpy.trace(numpy.matmul(D_t1,ovapm))
print('%.14f %.14fi\n' % (trace_ds.real,trace_ds.imag))

exit()

Ndip_z=0.0
#estrarre le 3 componenti del dipolo nucleare
D_ti=D_t1
Dp_ti=Dp_t1
#aggiungere repulsione nucleare
#Enuc_list.append(-func_t0*Ndip_z+Nuc_rep) #just in case of non-zero nuclear dipole
fockm_ti=bertha.get_realtime_fock(D_ti.T)

ene_list.append(numpy.trace(numpy.matmul(D_ti,fockm_ti)))

#trace of Dp_t1
trace_ds=numpy.trace(numpy.matmul(D_ti,ovapm))
print('%.14f %.14fi\n' % (trace_ds.real,trace_ds.imag))
niter=10
fo = open("err.txt","w")

fock_mid_backwd=numpy.copy(fock_mid_init)
for j in range(1,1):
    fock_mid_tmp=rtutil.mo_fock_mid_forwd_eval(bertha,numpy.copy(D_ti),fock_mid_backwd,j,numpy.float_(dt),\
        dipz_mat,C,C_inv,ovapm,ndim)

    fo.write('%.8f\n' % numpy.trace(numpy.matmul(S,D_ti)).real)
    Ah=numpy.conjugate(fock_mid_tmp.T)
    fo.write('Fock_mid hermitian: %s\n' % numpy.allclose(fock_mid_tmp,Ah))
    #transform fock_mid_init in MO basis
    fockp_mid_tmp=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(fock_mid_tmp,C))
    u=rtutil.exp_opmat(numpy.copy(fockp_mid_tmp),numpy.float_(dt))
    #u=scila.expm(-1.0j*fockp_mid_tmp*dt)
    #check u is unitary
    test_u=numpy.matmul(u,numpy.conjugate(u.T))
    if (not numpy.allclose(numpy.eye(u.shape[0]),test_u)):
        print('U is not unitary\n')
    #   else:
    #       print('U is unitary')
    #evolve the density in orthonormal basis
    #check the trace of density to evolve
    fo.write('tr of density to evolve: %.8f\n' % numpy.trace(Dp_ti).real)
    temp=numpy.matmul(Dp_ti,numpy.conjugate(u.T))
    Dp_ti_dt=numpy.matmul(u,temp)
    #backtransform Dp_ti_dt
    D_ti_dt=numpy.matmul(C,numpy.matmul(Dp_ti_dt,numpy.conjugate(C.T)))
    fo.write('%.8f\n' % numpy.trace(Dp_ti_dt).real)
    #dipole expectation for D_ti_dt
    dip_list.append(numpy.trace(numpy.matmul(dipz_mat,D_ti_dt)))
    #Energy expectation value at t = t_i_dt 
    fockm_ti=bertha.get_realtime_fock(D_ti_dt.T)

    ene_list.append(numpy.trace(numpy.matmul(D_ti_dt,fockm_ti_dt)))
    #####################################    
    #update D_ti and Dp_ti for the next step
    #message for debug
    fo.write('here I update the matrices Dp_ti and D_ti\n')
    D_ti=numpy.copy(D_ti_dt)
    Dp_ti=numpy.copy(Dp_ti_dt)
    fo.write('%.8f\n' % numpy.trace(Dp_ti).real)
    fo.write('%.8f\n' % numpy.trace(numpy.matmul(S,D_ti)).real)
    #update fock_mid_backwd for the next step
    fock_mid_backwd=numpy.copy(fock_mid_tmp)

fo.close







bertha.finalize()
