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

def exp_opmat(mat,dt,debug=False,odbg=sys.stderr):
    # first find eigenvectors and eigvals of F (hermitian)
    # and take the exponential of the -i*w*dt, w being eigvals of F
    #for debug
    if False :
       mat_h = numpy.conjugate(mat.T)
       diff_mat_h=mat-mat_h
       odbg.write("Max diff mat-mat_h: (%.15f, %.15fi)\n"%(numpy.max(diff_mat_h.real),numpy.max(diff_mat_h.imag)))
       odbg.write('Inputted  matrix is hermitian: %s\n'% 
        numpy.allclose(mat,mat_h,atol=1.e-14))
       if ( not numpy.allclose(mat,mat_h,atol=1.e-14)):
          rdiff = diff_mat_h.real
          maxrdiff = numpy.max(rdiff)
          idiff = diff_mat_h.imag
          maxidiff = numpy.max(idiff)
          #absvalue
          absmaxrdiff = numpy.abs(numpy.max(rdiff))
          absmaxidiff = numpy.abs(numpy.max(idiff))
          odbg.write('Abs values of max(rdiff), max(idiff): %.15f  %.15f \n' % (absmaxrdiff,absmaxidiff))
          odbg.write('Values of max(rdiff), max(idiff): %.15f  %.15fi \n' % (maxrdiff,maxidiff))
          matdim = mat.shape[0]
          for i in range(matdim):
            for j in range(matdim):
               if (numpy.abs(rdiff[i,j]) == numpy.abs(maxrdiff)):
                 if (mat[i,j].real*mat_h[i,j].real>0.0):
                    odbg.write("same sign for element (%i,%i).real\n" % (i,j))
                    diffrel= (mat[i,j].real - mat_h[i,j].real)/mat[i,j].real
                    odbg.write("For element (%i,%i).real relative diff: %.15f\n" % (i,j,diffrel))
                 else: 
                    odbg.write("opposit sign for element (%i,%i).real\n" % (i,j))
               if (numpy.abs(idiff[i,j]) == numpy.abs(maxidiff)):
                 if (mat[i,j].imag*mat_h[i,j].imag>0.0):
                    odbg.write("same sign for element (%i,%i).imag\n" % (i,j))
                    idiffrel= (mat[i,j].imag - mat_h[i,j].imag)/mat[i,j].imag
                    odbg.write("For element (%i,%i).imag  relative diff: %.15f\n" % (i,j,idiffrel))
                    
                 else: 
                    odbg.write("opposit sign for element (%i,%i).imag\n" % (i,j))

    try: 
       w,v=numpy.linalg.eigh(mat)
    except numpy.linalg.LinAlgError:
        print("Error in numpy.linalg.eigh of inputted matrix")
        return None

    diag=numpy.exp(-1.j*w*dt)
    # build the diagonal matrix
    # use numpy.diagflat(w)
    #   dmat=numpy.zeros(mat.shape,dtype=float)
    #   for num in range(diag.shape[0]):
    #       dmat[num,num]=diag[num]
    
    dmat=numpy.diagflat(diag)
    
    # for a general matrix Diag = M^(-1) A M
    # M is v 
    #try:
    #   v_i=numpy.linalg.inv(v)
    #except numpy.linalg.LinAlgError:
    #   return None 
       
    # transform back
    # matmul introduced in numpy 1.10 is preferred with respect 
    # numpy.dot 
    #tmp = numpy.matmul(dmat,v_i)
    tmp = numpy.matmul(dmat,numpy.conjugate(v.T))
    
    #in an orthonrmal basis v_inv = v.H
    
    mat_exp = numpy.matmul(v,tmp)
    
    return mat_exp

#######################################################################

def kick (Fmax, w, t, t0=0.0, s=0.0):

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

def gauss_env (Fmax, w, t, t0=3.0, s=0.2):
    
    #s : sqrt(variance) of gaussian envelop
    #w : the frequency of the carrier
    #Fmax : the maximum amplitude of the field
    #t0 :center of Gaussian envelope (in au time)
    
    func=Fmax*numpy.exp(-(t-t0)**2.0/(2.0*s**2.0))*numpy.sin(w*t)
    
    return func

#######################################################################

def envelope (Fmax, w, t, t0=0.0, s=0.0):
   
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

   func = Amp*numpy.sin(w*t)

   return func

#######################################################################

def sin_oc (Fmax, w, t, t0=0.0, s=0.0):
   
   # 1-oscillation-cycle sinusoid
   if (t >= 0.0 and t<= 2.00*numpy.pi/w):
      Amp = Fmax
   else:
      Amp = 0.0

   func = Amp*numpy.sin(w*t)

   return func

#######################################################################

def cos_env(Fmax, w, t, t0=0.0, n=20.0):

   #define the period (time for an oscillation cycle)
   #n is the number of oscillation cycle in the 
   # envelope
   oc=2.00*numpy.pi/w
   s=oc*n/2.0
   if (abs(t-s)<= s):
      func=numpy.sin(w*t)*Fmax*(numpy.cos(numpy.pi/2.0/s*(s-t)))**2.0
   else:
      func=0.0

   return func
 
#######################################################################

def lin_cos(Fmax, w, t, t0=0.0, s=0.0):
    
    if ( t <= 2.0*numpy.pi/w and t>=0.0):
        func=Fmax*w*t/2.0/numpy.pi*numpy.cos(w*t)
    else:
        func=Fmax*numpy.cos(w*t)
     
    return func

#######################################################################

def analytic(Fmax, w, t, t0=0.0, s=0.0):

   func = 0.0
 
   return func

#######################################################################

def dipole_selection(dipole,nshift,nocc,occlist,virtlist,odbg=sys.stderr,debug=False):
    
    a= occlist.pop(0)
    if debug:
       odbg.write("Selected occ. Mo: %s \n"% str(occlist))
       odbg.write("Selected virt. Mo: %s \n"% str(virtlist))
    tmp = dipole[nshift:,nshift:]
    offdiag = numpy.zeros((nshift,nshift),dtype=numpy.complex128)
    #diag = numpy.diagonal(tmp)
    #diagonal = numpy.diagflat(diag)
    nvirt = tmp.shape[0]-nocc
    odbg.write("n. virtual orbitals : %i" % nvirt)
    if (a == -1):
      for b in range(nvirt):
        for j in occlist:
          offdiag[nocc+b,j-1] = tmp[nocc+b,j-1]
    else:
      for b in virtlist:
        for j in  occlist:
          offdiag[b-1,j-1] = tmp[b-1,j-1]
    offdiag=(offdiag+numpy.conjugate(offdiag.T))
    #offdiag+=diagonal
    dipole[nshift:,nshift:] = offdiag

    return dipole

#######################################################################

def dipoleanalysis(dipole,dmat,occlist,virtlist,nshift,odbg=sys.stderr,debug=False):
#def dipoleanalysis(dipole,dmat,nocc,occlist,virtlist,nshift,odbg=sys.stderr,debug=False):
    
    shortl=occlist[1:]
    tot = len(shortl)*len(virtlist)
    if debug:
       odbg.write("Selected occ. Mo: %s \n"% str(shortl))
       odbg.write("Selected virt. Mo: %s \n"% str(virtlist))
    #if (a == -1): #HOMO-LUMO 
    #  i = nshift+nocc
    #  a = nshift+nocc+1
    #  res = dipole[i-1,a-1]*dmat[a-1,i-1] + dipole[a-1,i-1]*dmat[i-1,a-1]
    #else:
    res = numpy.zeros(tot,dtype=numpy.complex128)
    count = 0
    for i in shortl:
      for j in virtlist:    
         res[count] = dipole[i+nshift-1,j+nshift-1]*dmat[j+nshift-1,i+nshift-1] + dipole[j+nshift-1,i+nshift-1]*dmat[i+nshift-1,j+nshift-1]
         count +=1
    return res

#######################################################################

funcswitcher = {
    "kick": kick,
    "gauss_env": gauss_env,
    "envelope": envelope,
    "sin_oc": sin_oc,
    "cos_env": cos_env,
    "lin_cos": lin_cos,
    "analytic": analytic
     }
   
#######################################################################

def mo_fock_mid_forwd_eval(bertha, D_ti, fock_mid_ti_backwd, i, delta_t,
    dipole_z, C, C_inv, S, ndim, debug=False, odbg=sys.stderr, 
    impulsefunc="kick", fmax=0.0001, w=0.0, t0=0.0, sigma=0.0, propthresh=1.0e-6): 

   func = funcswitcher.get(impulsefunc, lambda: kick)

   fock_inter = numpy.zeros((ndim,ndim),dtype=numpy.complex128)   
   
   # D_ti is in AO basis 
   # transform in the MO ref basis
   
   Dp_ti = numpy.matmul(C_inv,numpy.matmul(D_ti,numpy.conjugate(C_inv.T)))
   k = 1
   t_arg = numpy.float_(i) * numpy.float_ (delta_t)
   fockmtx = bertha.get_realtime_fock(D_ti.T)
   
   pulse = func(fmax, w, t_arg, t0, sigma)
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

        u = exp_opmat(fockp_guess,delta_t,debug,odbg)

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
        
        pulse = func (fmax, w, t_arg + delta_t, t0, sigma)
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
                odbg.write("Norm: %.10e Should be <= %.10e \n" %(norm_f, propthresh))
                odbg.flush()
            if norm_f < (propthresh):
                if debug:
                  odbg.write(" Converged after %i interpolations\n" % (k))
                  odbg.write("   i = %i" % i)
                  odbg.write("   norm of D_ti_dt_(%i)-D_ti_dt(%i) : %.8f\n" % (k,k-1,norm_f))
                tr_dt = numpy.trace(numpy.matmul(S,D_ti_dt))
                break

        dens_test = numpy.copy(D_ti_dt)
        k += 1

   return fock_inter
