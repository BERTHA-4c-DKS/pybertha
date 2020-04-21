#pade.py - Fourier Transform utility by means of Pade approximant.

import sys
import os.path
import argparse
import numpy as np
import padeutil

listdamp = ""
for key in padeutil.dampswitcher:
    listdamp += key
    listdamp += " "   

listpulse = ""
for key1 in padeutil.funcswitcher:
    listpulse += key1
    listpulse += " "  
 
parser = argparse.ArgumentParser()

parser.add_argument("-f","--inputfile", help="Specify the dipole file (default: dipole.txt)", required=False, 
        type=str, default="dipole.txt")
parser.add_argument("-o","--outfile", help="Specify the output file  (default: ftoutput.txt)", required=False, 
        type=str, default="ftoutput.txt")
parser.add_argument("--limit", help="Number (even) of sampling points to be included (default: 10000)", required=False,
                                 default=10000, type=np.int32)
parser.add_argument("--frequency", help="Maximum frequency in the spectrum (eV) (default : 10.0)", required=False,
                             default=10.0, type=np.float64)
parser.add_argument("-f0", help="Set the frequency analysis for a specific set of transitions (eV) (default : 10.0)", required=False,
                             default="10;0;0", type=str)
parser.add_argument("--damping", help="Specify the form of the damping factor [" + listdamp + "] (default : exp)", required=False,
                             default="exp", type=str)
parser.add_argument("--gamma", help="Specify the damping parameter (default for exponential damp : 1.0e-3)", required=False,
                                 default=1.0e-3, type=np.float64)
parser.add_argument("--pad", help="Padding (add n points)", required=False,
                                 default=0, type=int)
parser.add_argument("--diagn", help="Check properties of G matrix (dafault : False)", required=False, 
        default=False, action="store_true")
parser.add_argument("--pulseFmax", help="Specify the pulse Fmax value (default: 1.0)", 
        default=1.0, type=np.float64)

parser.add_argument("--dw", help="Specify the frequency spacing (in e.V)", 
        default=0.03, type=float)

parser.add_argument("-p","--preproc", help="Dipole won't be zero at t=0", required=False, 
        default=True, action="store_false")
args = parser.parse_args()
preproc_zero =args.preproc

freqarray = args.f0.split(";")
freqarray = [float(m) for m in freqarray]
freqarray = np.array(freqarray)/27.2114

fnameinput = args.inputfile
fnameout = args.outfile
if not os.path.isfile(fnameinput):
    print("File ", fnameinput, " does not exist")
    exit(1)
#data for input
raw=np.loadtxt(fnameinput)
#to be read from input
limit = args.limit
#get the time spacing
dt=raw[1,0]-raw[0,0]
print('dt: %.8f' % dt)
dipole=raw[0:limit,1]
if preproc_zero:
  #subtract static dipole (or from input)
  dipole = dipole - dipole[0]
if args.pad:
   print("Padding on")
   pad=np.zeros(args.pad)
   dipole = np.append(dipole,pad)
nbins = int(np.floor(args.frequency/args.dw))
ffreq = args.frequency/27.2114
###################################
#print options
print("limit : %i, frequency : %.4f, gamma %.3e, damping : %s, nbins: %i \
        Fmax: %.8f,preproc_zero: %s, inputfile: %s, outfile: %s\n" % (limit,ffreq, args.gamma,\
      args.damping,nbins,args.pulseFmax, preproc_zero,args.inputfile,args.outfile))
###################################

lendip=dipole.shape[0]
print("Sampling points: %i\n" % lendip)

n1 = args.gamma
#try a damping
tx=np.linspace(0.0,(lendip-1)*dt,lendip)
#check

if (tx.shape[0] != lendip):
   print('Check dipole lenght')

func = padeutil.dampswitcher.get(args.damping, lambda: expdamp)    

fdamp = func(tx,n1)

inp=dipole*fdamp

a,b,N,G=padeutil.pade(inp,args.diagn)

#define the frequency abscissa
freq=np.linspace(0.0,ffreq,nbins)

z=np.exp(1.0j*dt*freq)

print(z.shape[0])

fft=np.zeros(z.shape[0],dtype=np.complex128)
fftimp=np.zeros(z.shape[0],dtype=np.complex128)


for l in range(0,z.shape[0]):
    A=0.0
    B=0.0
    for k in range(0,N+1):
        A+=a[k]*z[l]**k
        B+=b[k]*z[l]**k
    fft[l]=A/B

fft_list = []
for k in freqarray:
    A = 0.0
    B = 0.0
    z0 = np.exp(1.0j*dt*k)
    for k in range(0,N+1):
        A+=a[k]*z0**k
        B+=b[k]*z0**k
    fft_list.append(A/B/args.pulseFmax)

fft = fft/args.pulseFmax


np.savetxt('damp_dip.txt', np.c_[tx,inp], fmt='%.12f')
w=1.0j/dt*np.log(z)
np.savetxt(fnameout, np.c_[freq,fft.real,fft.imag,np.abs(fft)], fmt='%.12f')
#np.savetxt(fnameout, np.c_[freq,fft], fmt='%.12f')
fft_list=np.array(fft_list).imag
print("Dipole strength funcion values for selected freqs\n")
for i in range(freqarray.shape[0]):

 print("%.12e %.12e\n" % (freqarray[i],fft_list[i]*freqarray[i]))
