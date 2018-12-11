import numpy 
import sys
import os.path
sys.path.insert(0, '../src/')
import berthamod
import cdautil
from scipy.linalg import eig
from scipy.linalg import eigh

ovapcmp = berthamod.read_ovapfile ("ovapab.out")

cmatab = berthamod.read_vctfile ("vctab.out")

cmata = berthamod.read_vctfile ("vcta.out")
cmatb = berthamod.read_vctfile ("vctb.out")

density = numpy.matmul(cmatab,numpy.conjugate(cmatab.T))
trace = numpy.trace(numpy.matmul(ovapcmp,density))
ndimab = cmatab.shape[0]
noccab = cmatab.shape[1]
print("Trace of DS: %.8f %.8fi\n" % (trace.real, trace.imag))
cmat_join = cdautil.join_cmat(cmata,cmatb,ndimab)
print("Enter Loewdin")
print("Compute O")
O = numpy.matmul(numpy.conjugate(cmat_join.T),numpy.matmul(ovapcmp,cmat_join))
print("Check Trace of O")
otrace = numpy.trace(O)
print(otrace)
print("Compute inverse of O : O^-1")
oinv = numpy.linalg.inv(O)

print("Check O inversion")
test = numpy.matmul(O,oinv)
print("O*oinv =  1 : %s\n" % numpy.allclose(test,numpy.eye(test.shape[0]),atol=1.0e-14))
#compute left eigenvectors of O-1
try:
    w,z = eig(oinv,left=True,right=False)
except LinAlgError:
     print "Error in scipy.linalg.eig of O^-1"
fo = open("nocv_info.txt", "w")
fo.write("eigs of O^-1\n")
for i in range(w.shape[0]):
 fo.write("%.10e %.10e i\n" % (w[i].real,w[i].imag))
fo.close()
#check left eigenvectors
test = numpy.matmul(numpy.conjugate(z[:,1].T),oinv)
res = numpy.allclose(numpy.conjugate(z[:,1].T)*w[1],test)
print("left eigs O^-1 success: %s\n" % res)
da = numpy.diagflat(numpy.sqrt(w.real))
#compute zinv
zinv = numpy.linalg.inv(z)
LoewdinMat = numpy.matmul(z,numpy.matmul(da,zinv))
vct0 = numpy.matmul(cmat_join,LoewdinMat)
#density of the promolecule
dmat0 = numpy.matmul(vct0,numpy.conjugate(vct0.T))
#density of abduct AB
dmat = numpy.matmul(cmatab,numpy.conjugate(cmatab.T))
#compute density difference
tmp = dmat -dmat0
#check the trace of dmat and dmat0
trdmat = numpy.trace(numpy.matmul(dmat,ovapcmp))
trdmat0 = numpy.trace(numpy.matmul(dmat0,ovapcmp))
print("Trace of DmatAB %.8f %.8fi\n" % (trdmat.real,trdmat.imag))
print("Trace of Dmat0 %.8f %.8fi\n" % (trdmat0.real,trdmat0.imag))
#compute vmat (V = SDS)
vmat = numpy.matmul(ovapcmp,numpy.matmul(tmp,ovapcmp))
#diagonalize vmat
try:
    eigenval, zmat = eigh(vmat,ovapcmp, eigvals_only=False)
except LinAlgError:
     print "Error in scipy.linalg.eigh of vmat"
fo = open("nocv_eigv.txt", "w")
i = 0
for j in eigenval:
 i += 1
 fo.write(" %i %.8f\n" % (i,j))

for i in range(0,eigenval.shape[0]/2):
  fo.write("pair (%i):%.8f (%i):%.8f\n" % (i+1,eigenval[i],-i-1 ,eigenval[eigenval.shape[0]-1-i]))

fo.close()
#check orthormality of zmat coeff
test=numpy.matmul(numpy.conjugate(zmat.T),numpy.matmul(ovapcmp,zmat))
print("NOCV orthonormal: %s\n" % (numpy.allclose(test,numpy.eye(zmat.shape[0]),atol=1.0e-12)))
npairs = 12 #to be read from input
if (npair > eigenval.shape[0]/2):
  print("Wrong n. of pairs\n")
#
#
