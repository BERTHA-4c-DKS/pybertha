import numpy 
import sys
import os.path
from scipy.linalg import eig
from scipy.linalg import eigh

sys.path.insert(0, '../src/')
import berthamod
import cdautil

ovapcmp = berthamod.read_ovapfile ("ovapab.out")
cmatab = berthamod.read_vctfile ("vctab.out")
cmata = berthamod.read_vctfile ("vcta.out")
cmatb = berthamod.read_vctfile ("vctb.out")

density = numpy.matmul(cmatab,numpy.conjugate(cmatab.T))
trace = numpy.trace(numpy.matmul(density,ovapcmp.T))
ndimab = cmatab.shape[0]
noccab = cmatab.shape[1]
print(("Trace of DS: %.8f %.8fi\n" % (trace.real, trace.imag)))
cmat_join = cdautil.join_cmat(cmata,cmatb,ndimab)
print("Enter Loewdin")
print("Compute O")
O = numpy.matmul(numpy.conjugate(cmat_join.T),numpy.matmul(ovapcmp.T,cmat_join))
print("Check Trace of O")
otrace = numpy.trace(O)
print(otrace)
print("Compute inverse of O : O^-1")
oinv = numpy.linalg.inv(O)

print("Compute trace of O^-1\n")
print(("Trace of O^-1 : %.14f, %.14f i\n" % (numpy.trace(oinv).real, numpy.trace(oinv).imag)))
#print("Check O inversion")
#test = numpy.matmul(O,oinv)
#print("O*oinv =  1 : %s\n" % numpy.allclose(test,numpy.eye(test.shape[0]),atol=1.0e-14))
#compute left eigenvectors of O-1
try:
    w,z = eig(oinv,left=True,right=False)
except LinAlgError:
     print("Error in scipy.linalg.eig of O^-1")

#compute zinv
zinv = numpy.linalg.inv(z)

#test eigenvector
print("Compute Z x O^-1 x Z^-1 to check eigenvector\n")
temp = numpy.matmul(z,numpy.matmul(oinv,zinv))
print(("trace of Z x O^-1 x Z^-1 : %.14f, %.14f i\n" % (numpy.trace(temp).real,numpy.trace(temp).imag)))
D = numpy.sqrt(w)

val = 0.0 + 0.0j

for i in D:
  val += i
print(("sum of D = sqrt(eigs of O^-1): %.14f %.14f\n" % (val.real,val.imag)))
da = numpy.diagflat(numpy.sqrt(w))

LoewdinMat = numpy.matmul(z,numpy.matmul(da,zinv))
vct0 = numpy.matmul(cmat_join,LoewdinMat)
#density of the promolecule
dmat0 = numpy.matmul(vct0,numpy.conjugate(vct0.T))
#density of abduct AB
dmat = numpy.matmul(cmatab,numpy.conjugate(cmatab.T))
#compute density difference
tmp = dmat -dmat0
#check the trace of dmat and dmat0
trdmat = numpy.trace(numpy.matmul(dmat,ovapcmp.T))
trdmat0 = numpy.trace(numpy.matmul(dmat0,ovapcmp.T))
print(("Trace of DmatAB %.8f %.8fi\n" % (trdmat.real,trdmat.imag)))
print(("Trace of Dmat0 %.8f %.8fi\n" % (trdmat0.real,trdmat0.imag)))
#compute vmat (V = SDS)
vmat = numpy.matmul(ovapcmp.T,numpy.matmul(tmp,ovapcmp.T))
#diagonalize vmat
try:
    eigenval, zmat = eigh(vmat,ovapcmp.T, eigvals_only=False)
except LinAlgError:
     print("Error in scipy.linalg.eigh of vmat")
fo = open("nocv_eigv.txt", "w")
i = 0
for j in eigenval:
 i += 1
 fo.write(" %i %.8f\n" % (i,j))

for i in range(0,eigenval.shape[0]/2):
  fo.write("pair (%i):%.8f (%i):%.8f\n" % (-i-1,eigenval[i],i+1 ,eigenval[eigenval.shape[0]-1-i]))

fo.close()
#check orthormality of zmat coeff
test=numpy.matmul(numpy.conjugate(zmat.T),numpy.matmul(ovapcmp.T,zmat))
print(("NOCV orthonormal: %s\n" % (numpy.allclose(test,numpy.eye(zmat.shape[0]),atol=1.0e-10))))
