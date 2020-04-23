import numpy as np
import psi4
from pkg_resources import parse_version

##################################################################

def set_input(fgeom):
  geomobj = str()
  with open(fgeom,"r") as f:
   next(f)
   next(f)
   for line in f:
    geomobj +=str(line)
  geomobj += "symmetry c1" +"\n" +"no_reorient" +"\n" +"no_com"
  print(geomobj)
  mol =psi4.geometry(geomobj)
  f.close()
  return geomobj, mol

##################################################################

def exp_opmat(mat,dt):
    #first find eigenvectors and eigenvalues of F mat
    try:
       w,v=np.linalg.eigh(mat)
    except np.linalg.LinAlgError:
        print("Error in numpy.linalg.eigh of inputted matrix")
        return None

    diag=np.exp(-1.j*w*dt)

    dmat=np.diagflat(diag)

    # for a general matrix Diag = M^(-1) A M
    # M is v
    #try:
    #   v_i=np.linalg.inv(v)
    #except np.linalg.LinAlgError:
    #   return None

    # transform back
    #tmp = np.matmul(dmat,v_i)
    tmp = np.matmul(dmat,np.conjugate(v.T))

    #in an orthonrmal basis v_inv = v.H

    mat_exp = np.matmul(v,tmp)

    return mat_exp

##################################################################

def get_Fock(D, Hcore, I, f_type, basisset):
    # Build J,K matrices
    J = np.einsum('pqrs,rs->pq', I, D)
    if (f_type=='hf'):
        K = np.einsum('prqs,rs->pq', I, D)
        F = Hcore + J*np.float_(2.0) - K
        Exc=0.0
        J_ene = 0.0
    else:
        #D must be a psi4.core.Matrix object not a numpy.narray
        restricted = True
        if parse_version(psi4.__version__) >= parse_version('1.3a1'):
                build_superfunctional = psi4.driver.dft.build_superfunctional
        else:
                build_superfunctional = psi4.driver.dft_funcs.build_superfunctional
        sup = build_superfunctional(f_type, restricted)[0]
        sup.set_deriv(2)
        sup.allocate()
        vname = "RV"
        if not restricted:
            vname = "UV"
        potential=psi4.core.VBase.build(basisset,sup,vname)
        Dm=psi4.core.Matrix.from_array(D.real)
        potential.initialize()
        potential.set_D([Dm])
        nbf=D.shape[0]
        V=psi4.core.Matrix(nbf,nbf)
        potential.compute_V([V])
        potential.finalize()
        F = Hcore + J*np.float_(2.0) +V
        Exc= potential.quadrature_values()["FUNCTIONAL"]
        if sup.is_x_hybrid():
          alpha = sup.x_alpha()
          K = np.einsum('prqs,rs->pq', I, D)
          F += -alpha*K
          Exc += -alpha*np.trace(np.matmul(D,K))
        J_ene=2.00*np.trace(np.matmul(D,J))
    return J_ene,Exc,F

##################################################################
def set_params(filename="input.inp"):

    my_dict = {}
    with open(filename) as fileobj:
      for line in fileobj:
        key, value = line.split(":")
        my_dict[key.strip()] = value.strip()
    fileobj.close()
    imp_params = {}
    imp_params['Fmax'] = float(my_dict['F_max'])
    imp_params['w'] = float(my_dict['freq_carrier'])
    imp_params['s'] = float(my_dict['sigma'])
    imp_params['t0'] = float(my_dict['t0']) 
    imp_params['imp_type'] = my_dict['imp_type']
    
    calc_params ={}    
    calc_params['time_int']=float(my_dict['time_int'])
    calc_params['delta_t']=float(my_dict['delta_t'])
    calc_params['func_type'] =my_dict['func_type'] 
    calc_params['method']=my_dict['method_type']
    return imp_params,calc_params
##################################################################

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

##################################################################

def gauss_env (Fmax, w, t, t0=3.0, s=0.2):

    #s : sqrt(variance) of gaussian envelop
    #w : the frequency of the carrier
    #Fmax : the maximum amplitude of the field
    #t0 :center of Gaussian envelope (in au time)

    func=Fmax*np.exp(-(t-t0)**2.0/(2.0*s**2.0))*np.sin(w*t)

    return func

##################################################################

def envelope (Fmax, w, t, t0=0.0, s=0.0):

   if (t >= 0.0 and t<= 2.00*np.pi/w):
      Amp =(w*t/(2.00*np.pi))*Fmax
   elif (t > 2.00*np.pi/w and t < 4.00*np.pi/w):
      Amp = Fmax
   elif ( t >= 4.00*np.pi/w and t <= 6.00*np.pi/w):
      Amp = (3.00 -w*t/(2.00*np.pi))*Fmax
   elif ( t > 6.00*np.pi/w):
      Amp = 0.0
   else :

      Amp = 0.0

   func = Amp*np.sin(w*t)

   return func

##################################################################

def sin_oc (Fmax, w, t, t0=0.0, s=0.0):

   # 1-oscillation-cycle sinusoid
   if (t >= 0.0 and t<= 2.00*np.pi/w):
      Amp = Fmax
   else:
      Amp = 0.0

   func = Amp*np.sin(w*t)

   return func

##################################################################

def cos_env(Fmax, w, t, t0=0.0, s=20.0):

   #define the period (time for an oscillation cycle)
   #n is the number of oscillation cycles in the
   # envelope
   oc=2.00*np.pi/w
   n=oc*s/2.0
   if (abs(t-s)<= s):
      func=np.sin(w*t)*Fmax*(np.cos(np.pi/2.0/n*(n-t)))**2.0
   else:
      func=0.0

   return func

##################################################################

def analytic(Fmax, w, t, t0=0.0, s=0.0):

   func = 0.0

   return func

##################################################################

funcswitcher = {
    "kick": kick,
    "gauss_env": gauss_env,
    "envelope": envelope,
    "sin_oc": sin_oc,
    "cos_env": cos_env,
    "analytic": analytic
     }
   
##################################################################

def mo_fock_mid_forwd_eval(D_ti,fock_mid_ti_backwd,i,delta_t,H,I,dipole,\
                               C,C_inv,S,nbf,imp_opts,f_type,fout,basisset):

    t_arg=np.float_(i)*np.float_(delta_t)
    
    func = funcswitcher.get(imp_opts['imp_type'], lambda: kick)
    
    pulse = func(imp_opts['Fmax'], imp_opts['w'], t_arg,\
                        imp_opts['t0'], imp_opts['s'])

    #D_ti is in AO basis
    #transform in the MO ref basis
    Dp_ti= np.matmul(C_inv,np.matmul(D_ti,np.conjugate(C_inv.T)))
    
    k=1
    
    J_i,Exc_i,fock_mtx=get_Fock(D_ti,H,I,f_type,basisset)
    #add -pulse*dipole
    fock_ti_ao = fock_mtx - (dipole*pulse)

    #if i==0:
    #    print('F(0) equal to F_ref: %s' % np.allclose(fock_ti_ao,fock_mid_ti_backwd))
    
    #initialize dens_test !useless
    dens_test=np.zeros(Dp_ti.shape)

    # set guess for initial fock matrix
    fock_guess = 2.00*fock_ti_ao - fock_mid_ti_backwd
    #if i==0:
    #   print('Fock_guess for i =0 is Fock_0: %s' % np.allclose(fock_guess,fock_ti_ao))
    #transform fock_guess in MO basis
    while True:
        fockp_guess=np.matmul(np.conjugate(C.T),np.matmul(fock_guess,C))
        u=exp_opmat(fockp_guess,delta_t)
        #u=scipy.linalg.expm(-1.j*fockp_guess*delta_t) ! alternative routine
        test=np.matmul(u,np.conjugate(u.T))
    #print('U is unitary? %s' % (np.allclose(test,np.eye(u.shape[0]))))
        if (not np.allclose(test,np.eye(u.shape[0]))):
            Id=np.eye(u.shape[0])
            diff_u=test-Id
            norm_diff=np.linalg.norm(diff_u,'fro')
            print('from fock_mid:U deviates from unitarity, |UU^-1 -I| %.8f' % norm_diff)
    #evolve Dp_ti using u and obtain Dp_ti_dt (i.e Dp(ti+dt)). u i s built from the guess fock
    #density in the orthonormal basis
        tmpd=np.matmul(Dp_ti,np.conjugate(u.T))
        Dp_ti_dt=np.matmul(u,tmpd)
    #backtrasform Dp_ti_dt
        D_ti_dt=np.matmul(C,np.matmul(Dp_ti_dt,np.conjugate(C.T)))
    #build the correspondig Fock : fock_ti+dt
        
        dum1,dum2,fock_mtx=get_Fock(D_ti_dt,H,I,f_type,basisset)
        #update t_arg+=delta_t
        pulse_dt = func(imp_opts['Fmax'], imp_opts['w'], t_arg+delta_t,\
                        imp_opts['t0'], imp_opts['s'])
        fock_ti_dt_ao=fock_mtx -(dipole*pulse_dt)
        fock_inter= 0.5*fock_ti_ao + 0.5*fock_ti_dt_ao
    #update fock_guess
        fock_guess=np.copy(fock_inter)
        if k >1:
        #test on the norm: compare the density at current step and previous step
        #calc frobenius of the difference D_ti_dt_mo_new-D_ti_dt_mo
            diff=D_ti_dt-dens_test
            norm_f=np.linalg.norm(diff,'fro')
            if norm_f<(1e-6):
                tr_dt=np.trace(np.matmul(S,D_ti_dt))
                fout.write('converged after %i interpolations\n' % (k-1))
                fout.write('i is: %d\n' % i)
                fout.write('norm is: %.12f\n' % norm_f)
                fout.write('Trace(D)(t+dt) : %.8f\n' % tr_dt.real)
                break
        dens_test=np.copy(D_ti_dt)
        k+=1
        if k > 20:
         raise Exception("Numember of iterations exceeded (k>20)")
    return J_i,Exc_i,pulse,fock_ti_ao,fock_inter

##################################################################
# analysis based on MO-weighted dipole

def dipoleanalysis(dipole,dmat,nocc,occlist,virtlist,debug=False,HL=False):
    #virtlist can also contain occupied orbitals !check
    #just HOMO-LUMO vertical transition
    tot = len(occlist)*len(virtlist)
    if HL:
      i = nocc
      a = nocc+1
      res = dipole[i-1,a-1]*dmat[a-1,i-1] + dipole[a-1,i-1]*dmat[i-1,a-1]
    else:
      res = np.zeros(tot,dtype=np.complex128)
      count = 0
      for i in occlist:
        for j in virtlist:
           res[count] = dipole[i-1,j-1]*dmat[j-1,i-1] + dipole[j-1,i-1]*dmat[i-1,j-1]
           count +=1
    return res
