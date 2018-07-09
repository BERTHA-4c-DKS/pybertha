import numpy 
import sys
import scipy.linalg as scila

def exp_opmat(mat,dt):
#first find eigenvector of a  matrix defined as -iF*dt, F beig hermitian
# and take the exponential of the diagonal
#for our purpose mat=-iF, F being hermitian
    try: 
       w,v=numpy.linalg.eigh(dt*mat)
    except LinAlgError:
       print "error"
    diag=numpy.exp(-1.j*w)
#build the diagonal matrix
# use numpy.diagflat(w)
#   dmat=numpy.zeros(mat.shape,dtype=float)
#   for num in range(diag.shape[0]):
#       dmat[num,num]=diag[num]
    dmat=numpy.diagflat(diag)
# for a general matrix Diag = M^(-1) A M
# M is v 
    try:
       v_i=numpy.linalg.inv(v)
    except LinAlgError:
       print "error"
#   v_i=numpy.conjugate(v.T)
#transform back
#matmul introduced in numpy 1.10 is preferred with respect numpy.dot 
    tmp=numpy.matmul(dmat,v_i)
#for real matrices conjugation=transposition
    mat_exp=numpy.matmul(v,tmp)
    return mat_exp
#################################################

def kick(Fmax,t):
    func = 0.0
    if t >0:
      func=0.0
    elif (t == 0.0):
      func = Fmax
    else:
     print('t out of range')

    return func

#################################################

def mo_fock_mid_forwd_eval(bertha,D_ti,fock_mid_ti_backwd,i,delta_t,dipole_z,C,C_inv,S,ndim):

   fock_inter=numpy.zeros((ndim,ndim),dtype=numpy.complex128)   
   
   #D_ti is in AO basis 
   #transform in the MO ref basis
   Dp_ti= numpy.matmul(C_inv,numpy.matmul(D_ti,numpy.conjugate(C_inv.T)))
   k=1
   t_arg=numpy.float_(i)*numpy.float_(delta_t)
   fockmtx = bertha.get_realtime_fock(D_ti.T)
   fock_ti_ao=fockmtx+dipole_z*kick(0.0001,t_arg)
   #dipole matrix null for test
   dens_test=numpy.zeros((ndim,ndim),dtype=numpy.complex128)
   fock_guess = 2.00*fock_ti_ao - fock_mid_ti_backwd
   while True:
        fockp_guess=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(fock_guess,C))
        u=exp_opmat(fockp_guess,delta_t)
        #u=scila.expm(-1.j*fockp_guess*delta_t)
        test_u = numpy.matmul(u,numpy.conjugate(u.T))
        diff = test_u - numpy.eye(u.shape[0])
        print "max diff: ", numpy.max(diff)
        print('U is unitary(fock_mid) :%s' % numpy.allclose(test_u,numpy.eye(u.shape[0]),atol=1.e-14))
        tmpd=numpy.matmul(Dp_ti,numpy.conjugate(u.T))
        Dp_ti_dt=numpy.matmul(u,tmpd)
    #backtrasform Dp_ti_dt
        D_ti_dt=numpy.matmul(C,numpy.matmul(Dp_ti_dt,numpy.conjugate(C.T)))
    #build the correspondig Fock , fock_ti+dt
        
        fock_ti_dt_ao=bertha.get_realtime_fock(D_ti_dt.T)+dipole_z*kick(0.0001,t_arg+delta_t)
        fock_inter= 0.5*fock_ti_ao + 0.5*fock_ti_dt_ao
    #update fock_guess
        fock_guess=numpy.copy(fock_inter)
        if k >1:
        #test on the norm: compare the density at current step and previous step
        #calc frobenius of the difference D_ti_dt_mo_new-D_ti_dt_mo
            diff=D_ti_dt-dens_test
            norm_f=numpy.linalg.norm(diff,'fro')
            print norm_f
            if norm_f<(1e-6):
                tr_dt=numpy.trace(numpy.matmul(S,D_ti_dt))
                break
        dens_test=numpy.copy(D_ti_dt)
        k+=1
   return fock_inter
