import scipy.linalg as scila
import numpy 
import sys

#######################################################################

def progress_bar (count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush() 

#######################################################################

def exp_opmat(mat,dt):
    # first find eigenvector of a  matrix defined as -iF*dt, F beig hermitian
    # and take the exponential of the diagonal
    # for our purpose mat=-iF, F being hermitian
    
    try: 
       w,v=numpy.linalg.eigh(dt*mat)
    except LinAlgError:
       return None 

    diag=numpy.exp(-1.j*w)
    # build the diagonal matrix
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
       return None 
       
    # v_i=numpy.conjugate(v.T)
    # transform back
    # matmul introduced in numpy 1.10 is preferred with respect 
    # numpy.dot 
    tmp = numpy.matmul(dmat,v_i)
    
    #for real matrices conjugation=transposition
    
    mat_exp = numpy.matmul(v,tmp)
    
    return mat_exp

#######################################################################

def kick (Fmax, w, t):

    w = 0.0
    func = 0.0
    if t > 0:
      func = 0.0
    elif (t == 0.0):
      func = Fmax
    else:
      return None 
      # t out of range

    return func

#######################################################################

def gauss_env (k, w, t, t0=3.0, s=0.2):
    
    #s is the pulse width
    #when used the envelope must be multiplied by sin(wt)
    # (typically a few time step
    w = 0.0
    
    func=k*numpy.exp(-(t-t0)**2.0/(2.0*s**2.0))
    
    return func

#######################################################################

def envelope (Fmax, w, t):
   
   if (t >= 0.0 and t<= 2.00*numpy.pi/w):
      Amp =(w*t/(2.00*numpy.pi))*Fmax
   elif (t > 2.00*numpy.pi/w and t < 4.00*numpy.pi/w):
      Amp = Fmax
   elif ( t >= 4.00*numpy.pi/w and t <= 6.00*numpy.pi/w):
      Amp = (3.00 -w*t/(2.00*numpy.pi))*Fmax
   elif ( t > 6.00*numpy.pi/w):
      Amp = 0.0
   else :
      Amp = 0.0

   return Amp

#######################################################################

def sin_env (Fmax, w, t):
   
   # 1-oscillation 
   if (t >= 0.0 and t<= 2.00*numpy.pi/w):
      Amp = Fmax
   else:
      Amp = 0.0

   return Amp

#######################################################################

funcswitcher = {
    "kick": kick,
    "gauss_env": gauss_env,
    "envelope": envelope,
    "sin_env": sin_env
     }
   
#######################################################################

def mo_fock_mid_forwd_eval(bertha, D_ti, fock_mid_ti_backwd, i, delta_t,
    dipole_z, C, C_inv, S, ndim, debug=False, odbg=sys.stderr, 
    impulsefunc="kick", fmax=0.0001, w=0.0): 

   func = funcswitcher.get(impulsefunc, lambda: kick)

   fock_inter = numpy.zeros((ndim,ndim),dtype=numpy.complex128)   
   
   # D_ti is in AO basis 
   # transform in the MO ref basis
   
   Dp_ti = numpy.matmul(C_inv,numpy.matmul(D_ti,numpy.conjugate(C_inv.T)))
   k = 1
   t_arg = numpy.float_(i) * numpy.float_ (delta_t)
   fockmtx = bertha.get_realtime_fock(D_ti.T)
   
   pulse = func(fmax, w, t_arg)
   if debug: 
       odbg.write("Pulse: %.8f\n"%(pulse))

   if pulse is None:
     return None 

   fock_ti_ao = fockmtx - (dipole_z * pulse)
   # dipole matrix null for test
   dens_test = numpy.zeros((ndim,ndim),dtype=numpy.complex128)
   fock_guess = 2.00*fock_ti_ao - fock_mid_ti_backwd
   while True:
        fockp_guess = numpy.matmul(numpy.conjugate(C.T), \
                numpy.matmul(fock_guess,C))

        u = exp_opmat(fockp_guess,delta_t)

        if u is None:
          return None 

        #u=scila.expm(-1.j*fockp_guess*delta_t)

        if debug:
          test_u = numpy.matmul(u,numpy.conjugate(u.T))
          diff = test_u - numpy.eye(u.shape[0])
          maxdiff = numpy.max(diff)
          odbg.write("U is unitary(fock_mid) : %s"% 
                numpy.allclose(test_u,numpy.eye(u.shape[0]),atol=1.e-14))
          odbg.write("  max diff: %.12e %.12e \n"%
                (maxdiff.real, maxdiff.imag))

        tmpd = numpy.matmul(Dp_ti,numpy.conjugate(u.T))
        Dp_ti_dt = numpy.matmul(u,tmpd)
        #backtrasform Dp_ti_dt
        D_ti_dt = numpy.matmul(C,numpy.matmul(Dp_ti_dt,numpy.conjugate(C.T)))
        #build the correspondig Fock , fock_ti+dt
        
        pulse = func (fmax, w, t_arg + delta_t)
        if debug: 
            odbg.write("Pulse: %.8f\n"%(pulse))

        if pulse is None:
          return None 

        fock_ti_dt_ao=bertha.get_realtime_fock(D_ti_dt.T)-(dipole_z*pulse)
        fock_inter = 0.5*fock_ti_ao + 0.5*fock_ti_dt_ao
        #update fock_guess
        fock_guess = numpy.copy(fock_inter)
        if k > 1:
            # test on the norm: compare the density at current step and previous step
            # calc frobenius of the difference D_ti_dt_mo_new-D_ti_dt_mo
            diff = D_ti_dt-dens_test
            norm_f = numpy.linalg.norm(diff,'fro')
            if debug:
              odbg.write("Norm: %.10e\n" %(norm_f))
              odbg.flush()
            if norm_f < (1e-6):
                if debug:
                  odbg.write(" Converged after %i interpolations\n" % (k))
                  odbg.write("   i = %i" % i)
                  odbg.write("   norm of D_ti_dt_(%i)-D_ti_dt(%i) : %.8f\n" % (k,k-1,norm_f))
                tr_dt = numpy.trace(numpy.matmul(S,D_ti_dt))
                break

        dens_test = numpy.copy(D_ti_dt)
        k += 1

   return fock_inter
