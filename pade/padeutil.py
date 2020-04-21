import numpy as np

#######################################################################

def pade(data, diagnostic = False):

    M = len(data)
    N = int(np.floor(M/2))
    print('n : %i'% N)
    print('number of sampling point : %i'% M)
    G=np.zeros((N,N))
    d=-1.00*np.array(data[N:])
    #check dimension of d
    if (d.shape[0] !=N):
        print('wrong dim of d vec')
    for k  in range(0,N):
      for m in range(0,N):
        G[k,m]=data[N-1-m+k]

    if (diagnostic):
        #diagnostic
        #compute eigenvalue of G to check non singularity
        try:
            sigma = np.linalg.eigvals(G)
        except np.linalg.LinAlgError:
            print("Error in np.linalg.eigvals of G mat")
        #check G eigvals
        counter=0
        for l in sigma:
           if (l <= 1.0e-12):
              counter+=1
        print('Found %i eigvals near to zero' % counter)
        #check if G is (semi-)negative or (semi-)positive definite
          
        for l in sigma:
           if (np.abs(l)*np.sign(l) != l):
              print('Found a negative eigvalue')
              break
    

    try:
         b = np.linalg.solve(G,d)
    except np.linalg.LinAlgError:
        print("Error in np.linalg.solve")

    b=np.insert(b,0,1.0)
    print('b dim: %i' % b.shape[0])
    #a0=c0
    a=np.zeros(N+1)
    for k in range(1,N+1):
      for m in range(0,k+1):
       a[k]+=b[m]*data[k-m]
    a[0]=data[0]
    return a,b,N,G

#######################################################################

def nodamp(t=0.0, n1=0.0):
    
    func = 1.00    

    return func

#######################################################################

def expdmp(t, n1):
   
    func = np.exp(-t*n1)
    
    return func

#######################################################################

def gaussdmp(t, n1):
    
    func = np.exp(-(t**2.0)*n1)
    
    return func

#######################################################################
def hanndmp(t, n1):
  print("dt set to 0.1 ! check !")
  n = t/0.1
  N = n1 -1 #n1 has the same value of limit
  # l window lenght
  # N number of intervals
  func = 0.5*(1.0 - np.cos(2.0*np.pi*n/N))
  return func

#######################################################################

dampswitcher = {
    "nodamp": nodamp,
    "gauss": gaussdmp,
    "exp": expdmp,
    "hann": hanndmp
     }
   
#######################################################################
#Misc functions
def gauss_env (Fmax, w, t, t0=3.0, s=0.2):
    
    #s : the pulse width
    #w : the frequency of the carrier
    #Fmax : the maximum amplitude of the field
    #t0 :center of Gaussian envelope (in au time)
    
    func=np.exp(-(t-t0)**2.0/(2.0*s**2.0))*np.sin(w*t)
    
    return func

#######################################################################

def analytic(Fmax, w, t, t0=0.0, s=0.0):

   func = 0.0
 
   return func

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

funcswitcher = {
    "kick": kick,
    "gauss_env": gauss_env,
    "analytic": analytic
     }

#######################################################################

def fft_gauss_env(Fmax, w, k, t0, w0):
    func = -0.5*Fmax*1.0j*k*(np.exp(-0.5*(k**2.)*(w0 + w)**2. +1.0j*t0*(w0 + w)) - np.exp(-0.5*(k**2)*(w - w0)**2.0 +1.0j*t0*(w - w0)))
    return func 
#######################################################################
