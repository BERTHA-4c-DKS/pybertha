import argparse
import ctypes
import numpy
import sys
import re

import os.path

from numpy.linalg import eigvalsh
from scipy.linalg import eigh

import time

if __name__ == "__main__":
    
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
    parser.add_argument("-d", "--debug", help="Debug on, prints debug info to debug_info.txt", required=False, 
            default=False, action="store_true")
    parser.add_argument("-v", "--verbosity", help="Verbosity level 0 = minim, -1 = print iteration info, " + 
            "1 = maximum (defaul -1)", required=False, default=-1, type=int)
    parser.add_argument("--thresh", help="det threshold (default = 1.0e-12)", required=False, 
            type=numpy.float64, default=1.0e-12)
    parser.add_argument("--wrapperso", help="set wrapper SO (default = ../../lib/bertha_wrapper.so)", 
            required=False, type=str, default="../lib/bertha_wrapper.so")
    parser.add_argument("--berthamodpath", help="set berthamod path (default = ../src)", 
            required=False, type=str, default="../src")
    parser.add_argument("-np", "--npairs", help="Specify the numerber of nocv-pair density (default: 0)", required=False,
               default=0, type=int)
    parser.add_argument("--deltax", help="cube dx step (default = 0.2)", required=False, 
            type=numpy.float64, default=0.2)
    parser.add_argument("--deltay", help="cube dy step (default = 0.2)", required=False, 
            type=numpy.float64, default=0.2)
    parser.add_argument("--deltaz", help="cube dz step (default = 0.2)", required=False, 
            type=numpy.float64, default=0.2)
    parser.add_argument("--lmargin", help="cube margin parameter (default = 10.0)", required=False, 
            type=numpy.float64, default=10.0)
    args = parser.parse_args()

    sys.path.insert(0, args.berthamodpath)

    import berthamod
    import cdautil
    
    print("Options: ")
    print(args) 
    print("")
    print("")
    
    if not os.path.isfile(args.wrapperso):
        print("SO File ", args.wrapperso, " does not exist")
        exit(1)
    
    bertha = berthamod.pybertha(args.wrapperso)
    
    fittcoefffname = args.fitcoefffile
    vctfilename = args.vctfile
    ovapfilename = args.ovapfile
    
    fnameinput = args.inputfile
    if not os.path.isfile(fnameinput):
        print("File ", fnameinput, " does not exist")
        exit(1)
    
    fittfname = args.fittfile
    if not os.path.isfile(fittfname):
        print("File ", fittfname , " does not exist")
        exit(1)
    
    verbosity = args.verbosity
    dumpfiles = int(args.dumpfiles)
    
    bertha.set_fittcoefffname(fittcoefffname)
    bertha.set_ovapfilename(ovapfilename)
    bertha.set_vctfilename(vctfilename)
    bertha.set_fnameinput(fnameinput)
    bertha.set_fittfname(fittfname)
    bertha.set_thresh(args.thresh)
    
    bertha.set_verbosity(verbosity)
    bertha.set_dumpfiles(dumpfiles)
    
    bertha.set_densitydiff(1)
    
    bertha.init()
    
    ndim = bertha.get_ndim()
    nshift = bertha.get_nshift()
    nocc = bertha.get_nocc()
    sfact = bertha.get_sfact()
    nopen = bertha.get_nopen()
    
    print("Verbosity       : ", verbosity)
    print("Dumpfiles       : ", dumpfiles)
    print("")
    print("Matrix dimension: ", ndim)
    print("            nocc: ", nocc)
    print("          nshift: ", nshift)
    print("           nopen: ", nopen)
    print("     level shift: ", sfact)
    print("")
    
    start = time.time()
    cstart = time.process_time() 
    
    ovapm, eigem, fockm, eigen = bertha.run()
    
    end = time.time()
    cend = time.process_time() 
    
    print("Totaltime:    ", end - start, " (CPU time: " , cend - cstart, ") s ")
    print("MainRun Time: ", bertha.get_mainruntime() , \
            " (CPU time: " , bertha.get_mainrunctime(), ") s ")
    
    sys.stdout.flush()
    
    
    if (fockm is None) or (eigen is None) or (fockm is None) \
            or (eigen is None):
        print("Error in bertha run")
        exit(-1)
    
    print("")
    print("Final results ")
    for i in range(nocc+nopen):
        print("eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact))
        
    print("      lumo       %20.8f"%(eigen[i+nshift+1]))
    
    erep = bertha.get_erep()
    etotal = bertha.get_etotal()
    
    print("")
    print("total electronic energy  = %20.8f"%(etotal-(sfact*nocc)))
    print("nuclear repulsion energy = %20.8f"%(erep))
    print("total energy             = %20.8f"%(etotal+erep-(sfact*nocc)))
    
    #bertha.finalize()
    
    occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)
    
    iocc = 0
    for i in range(ndim):
        if i >= nshift and iocc < nocc:
            for j in range(ndim):
                occeigv[j, iocc] = eigem[j, i]
            iocc = iocc + 1
    
    print("")
    print("Compute density matrix ")
    density = numpy.matmul(occeigv, numpy.conjugate(occeigv.transpose()), out=None)
    density = numpy.matmul(density, ovapm)
    print("Trace  ")
    trace = density.trace()
    print("(%20.10f, %20.10fi)"%(trace.real, trace.imag))
    
    
    # NOCV analysis start here  ! check matrix from ovap file, transposition needed
    from scipy.linalg import eig

    #check vct and ovap abduct !REMOVE
    #if not os.path.isfile(vctfilename):
    #    print("File ", vctfilename, " does not exist")
    #    exit(1)
    
    #if not os.path.isfile(ovapfilename):
    #    print("File ", ovapfilename, " does not exist")
    #    exit(1)
    #check vct of frags

    if not os.path.isfile("vcta.out"):
        print("File vcta.out does not exist")
        exit(1)
    
    if not os.path.isfile("vctb.out"):
        print("File vctb.out does not exist")
        exit(1)
    
    npairs = args.npairs
    drx = args.deltax
    dry = args.deltay
    drz = args.deltaz
    margin = args.lmargin
    #ovapcmp = berthamod.read_ovapfile ("ovapab.out")
    
    #cmatab = berthamod.read_vctfile (vctfilename)
    cmatab = occeigv
    
    cmata = berthamod.read_vctfile ("vcta.out")
    cmatb = berthamod.read_vctfile ("vctb.out")
    
    density = numpy.matmul(cmatab,numpy.conjugate(cmatab.T))
    trace = numpy.trace(numpy.matmul(ovapm,density))
    ndimab = cmatab.shape[0]
    noccab = cmatab.shape[1]
    
    print(("Trace of DS: %.8f %.8fi\n" % (trace.real, trace.imag)))
    cmat_join = cdautil.join_cmat(cmata,cmatb,ndimab)

    if cmat_join is None:
        print("Error in joining cmat")
        exit(1)

    print("Enter Loewdin")
    print("Compute O")
    O = numpy.matmul(numpy.conjugate(cmat_join.T),numpy.matmul(ovapm,cmat_join))
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
    
    val = 0.0 + 0.0j
    
    for i in w:
      val += i
    print(("sum of eigs of O^-1: %.14f %.14f\n" % (val.real,val.imag)))
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
    #TODO add the D_anti contribution
    trdmat = numpy.trace(numpy.matmul(dmat,ovapm))
    trdmat0 = numpy.trace(numpy.matmul(dmat0,ovapm))
    print(("Trace of DmatAB %.8f %.8fi\n" % (trdmat.real,trdmat.imag)))
    print(("Trace of Dmat0 %.8f %.8fi\n" % (trdmat0.real,trdmat0.imag)))
    #compute vmat (V = SDS)
    vmat = numpy.matmul(ovapm,numpy.matmul(tmp,ovapm))
    #diagonalize vmat

    try:
        eigenval, zmat = eigh(vmat,ovapm, eigvals_only=False)
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
    test=numpy.matmul(numpy.conjugate(zmat.T),numpy.matmul(ovapm,zmat))
    print(("NOCV orthonormal: %s\n" % (numpy.allclose(test,numpy.eye(zmat.shape[0]),atol=1.0e-10))))
    #npairs = 2 #to be read from input
    if (npairs > eigenval.shape[0]/2):
        print("Wrong n. of pairs\n")
    for i in range(npairs):
        j = i + 1
        label = "pair"+str(j)
        tmp =  zmat[:,i]
        d1 = numpy.outer(tmp,numpy.conjugate(tmp))
        #check if d1 sum to 1
        trace = numpy.trace(numpy.matmul(d1,ovapm))
        print(("trace of nocv_-%i : %.8f\n" % (j,trace.real)))
        tmp =  zmat[:,-i-1]
        d2 = numpy.outer(tmp,numpy.conjugate(tmp))
        #check if d2 sum to 1
        trace = numpy.trace(numpy.matmul(d2,ovapm))
        print(("trace of nocv_+%i : %.8f\n" % (j,trace.real)))
        deltanocv = eigenval[i]*(d1 - d2)
    
        bertha.density_to_cube(d1.T, "nocv-"+str(j)+".cube", margin, drx, dry, drz )  
        bertha.density_to_cube(d2.T, "nocv+"+str(j)+".cube", margin, drx, dry, drz )  
        bertha.density_to_cube(deltanocv.T, label+".cube", margin, drx, dry, drz )  
    
    bertha.finalize()