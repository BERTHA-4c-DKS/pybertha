import argparse
import os.path
import ctypes
import numpy
import sys
import re
import os

import json
from json import encoder

import time

import scipy.linalg as scila
from numpy.linalg import eigvalsh
from scipy.linalg import eigh

berthamodpath = os.getenv("BERTHA_API_PATH", "../src")

sys.path.insert(0, berthamodpath)
import berthamod
import rtutil

##########################################################################################

def get_json_data(args, j, niter, ndim, ene_list, \
        dip_list, D_ti, fock_mid_backwd, dip_mat, C, \
        C_inv, ovapm, Dp_ti, weight_list):

    json_data = {}
    for arg in vars(args):
        json_data[arg] = getattr(args, arg)

    othervals = {
            'j': j,
            'niter': niter,
            'ndim' : ndim,
            'ene_list_REAL': numpy.real(ene_list).tolist(),
            'ene_list_IMAG': numpy.imag(ene_list).tolist(),
            'dip_list_REAL': numpy.real(dip_list).tolist(),
            'dip_list_IMAG': numpy.imag(dip_list).tolist(),
            'D_ti_REAL' : numpy.real(D_ti).tolist(),
            'D_ti_IMAG' : numpy.imag(D_ti).tolist(),
            'fock_mid_backwd_REAL' : numpy.real(fock_mid_backwd).tolist(),
            'fock_mid_backwd_IMAG' : numpy.imag(fock_mid_backwd).tolist(),
            'dip_mat_REAL': numpy.real(dip_mat).tolist(),
            'dip_mat_IMAG': numpy.imag(dip_mat).tolist(),
            'C_REAL': numpy.real(C).tolist(), 
            'C_IMAG': numpy.imag(C).tolist(),
            'C_inv_REAL': numpy.real(C_inv).tolist(),
            'C_inv_IMAG': numpy.imag(C_inv).tolist(),
            'ovapm_REAL': numpy.real(ovapm).tolist(),
            'ovapm_IMAG': numpy.imag(ovapm).tolist(),
            'Dp_ti_REAL': numpy.real(Dp_ti).tolist(),
            'Dp_ti_IMAG': numpy.imag(Dp_ti).tolist(),
            'weight_list_REAL': numpy.real(numpy.array(weight_list)).tolist(),
            'weight_list_IMAG': numpy.imag(numpy.array(weight_list)).tolist()
            }

    json_data.update(othervals)

    return json_data

##########################################################################################

def single_point (args, bertha):

    fittcoefffname = args.fitcoefffile
    vctfilename = args.vctfile
    ovapfilename = args.ovapfile
    
    fnameinput =""
    fittfname =""

    generatedinout = False

    if args.geometry != "":
        
        if not os.path.isfile (args.geometry):
            print("File ", args.geometry, " does not exist")
            exit(1)

        if args.fittset == "" or args.obs == "":
            print("Need to specify fittset and basis set for each elements")
            exit(1)

        # generate input 
        import pybgen

        pygenoption = pybgen.berthainputoption

        pygenoption.inputfile = args.geometry
        pygenoption.jsonbasisfile = args.jsonbasisfile
        pygenoption.fittset = args.fittset
        pygenoption.basisset = args.obs
        pygenoption.functxc = args.func
        pygenoption.convertlengthunit = args.convertlengthunit
        pygenoption.maxit = args.berthamaxit

        fnameinput = str(uuid.uuid4())
        fittfname = str(uuid.uuid4())  

        pygenoption.berthainfname = fnameinput
        pygenoption.berthafittfname = fittfname

        for filename in [fnameinput , fittfname]:

            if os.path.isfile(filename):
                print("File ", filename, " will be overwritten")
                try:
                    os.remove(filename)
                except OSError:
                    pass

        pybgen.generateinputfiles (pygenoption)  

        generatedinout = True
    else:
        fnameinput = args.inputfile
        if not os.path.isfile(fnameinput):
            print("File ", fnameinput, " does not exist")
            return None, 
    
        fittfname = args.fittfile
        if not os.path.isfile(fittfname):
            print("File ", fittfname , " does not exist")
            return None, 
    
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

    berthastart = time.time()
    berthacstart = time.process_time() 
    
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

    sys.stdout.flush()
    
    ovapm, eigem, fockm, eigen = bertha.run()
    if (fockm is None) or (eigen is None) or (fockm is None) \
            or (eigen is None):
        print("Error in bertha run")
        return None, 
    
    bertha.set_densitydiff(0)

    berthaend = time.time()
    berthacend = time.process_time() 
   
    print("Single-Point time: %15.5f"%(berthaend - berthastart), \
                    " (CPU time: %15.5f"%(berthacend - berthacstart), ") s ")
    
    sys.stdout.flush()

    if generatedinout:
        
        for filename in [fnameinput , fittfname]:
            if os.path.isfile(filename):
                print("File ", filename, " will be removed")
                try:
                    os.remove(filename)
                except OSError:
                    pass
    print("")
    print("Final results ")
    sum=0.0
    for i in range(nocc+nopen):
        print("eigenvalue %5d %20.8f"%(i+1, eigen[i+nshift]-sfact))
        sum+=eigen[i+nshift]-sfact
    print("      lumo       %20.8f"%(eigen[i+nshift+1]))
    print("      Sum of eigen : ", sum)
    erep = bertha.get_erep()
    etotal = bertha.get_etotal()
    
    print("")
    print("total electronic energy  = %20.8f"%(etotal-(sfact*nocc)))
    print("nuclear repulsion energy = %20.8f"%(erep))
    print("total energy             = %20.8f"%(etotal+erep-(sfact*nocc)))
    print(" ")
 
    return ovapm, eigem, fockm, eigen

##########################################################################################

def check_and_covert (mat_REAL, mat_IMAG, ndim):

    if ((mat_REAL.shape == mat_IMAG.shape) and 
        (mat_REAL.shape[0] == mat_REAL.shape[1]) and 
        (mat_REAL.shape[0] == ndim)):

        mat = numpy.zeros((ndim,ndim),dtype=numpy.complex128)

        for i in range(ndim):
            for j in range(ndim):
                mat[i, j] = numpy.complex128(complex(mat_REAL[i][j], 
                    mat_IMAG[i][j]))

        return mat
        
    return None


##########################################################################################

def main_loop (j, niter, bertha, pulse, pulseFmax, pulsew, propthresh, pulseS, t0,
        fo, D_ti, fock_mid_backwd, dt, dip_mat, C, C_inv, ovapm, 
        ndim, debug, Dp_ti, dip_list, ene_list, weight_list=None, 
        select="-2 ; 0 & 0"):

  
    fock_mid_tmp = rtutil.mo_fock_mid_forwd_eval(bertha, numpy.copy(D_ti), \
            fock_mid_backwd, j, numpy.float_(dt), dip_mat, C, C_inv, ovapm, \
            ndim, debug, fo, pulse, pulseFmax, pulsew, t0, pulseS, propthresh)
    
    if (fock_mid_tmp is None):
        print("Error accurs in mo_fock_mid_forwd_eval")
        return None, None, None
   
    if debug:
      fo.write('%.8f\n' % numpy.trace(numpy.matmul(ovapm,D_ti)).real)
      Ah = numpy.conjugate(fock_mid_tmp.T)
      fo.write('Fock_mid hermitian: %s\n' % numpy.allclose(fock_mid_tmp,Ah,atol=1.e-14))

    # transform fock_mid_tmp in MO basis
    fockp_mid_tmp = numpy.matmul(numpy.conjugate(C.T),numpy.matmul(fock_mid_tmp, C))
    u = rtutil.exp_opmat(numpy.copy(fockp_mid_tmp),numpy.float_(dt),debug,fo)
    # u=scila.expm(-1.0j*fockp_mid_tmp*dt)
    # check u is unitary
    test_u = numpy.matmul(u,numpy.conjugate(u.T))
    if (not numpy.allclose(numpy.eye(u.shape[0]),test_u,atol=1.e-14)):
        print("WARNING: U is not unitary\n")
   
    #evolve the density in orthonormal basis
    #check the trace of density to evolve
    if debug:
      fo.write('trace of density to evolve: %.8f\n' % numpy.trace(Dp_ti).real)
    temp=numpy.matmul(Dp_ti,numpy.conjugate(u.T))
    Dp_ti_dt = numpy.matmul(u,temp)
    #backtransform Dp_ti_dt
    D_ti_dt=numpy.matmul(C,numpy.matmul(Dp_ti_dt,numpy.conjugate(C.T)))
    if debug:
      fo.write('  Trace of D_ti_dt %.8f\n' % numpy.trace(Dp_ti_dt).real)
    #dipole expectation for D_ti_dt
    dip_list.append(numpy.trace(numpy.matmul(dip_mat,D_ti_dt)))

    if (pulse == "analytic"):
        molist = select.split("&")
        occlist = molist[0].split(";")
        occlist = [int(m) for m in occlist]
        virtlist = molist[1].split(";")
        virtlist = [int(m) for m in virtlist]
        if (occlist[0] != -2):
            #dipole analysis
            nshift = ndim/2
            dipz_mo = numpy.matmul(numpy.conjugate(C.T),numpy.matmul(dip_mat,C))
            res=rtutil.dipoleanalysis(dipz_mo,Dp_ti_dt,occlist,virtlist,int(nshift),fo,debug)
            if (weight_list != None):
                weight_list.append(res)
    if debug:
      fo.write("Dipole: %.12e\n"%(numpy.trace(numpy.matmul(dip_mat,D_ti_dt)).real))
   
    #Energy expectation value at t = t_i_dt 
    fockm_ti_dt = bertha.get_realtime_fock(D_ti_dt.T)
    
    """
    # to dump density og a grid 
    npoints = 50
    grid = numpy.zeros((npoints*npoints*npoints, 4))
    grid = numpy.ascontiguousarray(grid, dtype=numpy.double)
    
    dstep = 0.2
    w = 1.0
    idx = 0

    x = 8.0
    for ix in range(npoints):
        y = 4.0
        for iy in range(npoints):
            z = -4.0
            for iz in range(npoints):
                grid[idx,0] = x
                grid[idx,1] = y
                grid[idx,2] = z
                grid[idx,3] = w

                idx += 1

                z += dstep
            y += dstep
        x += dstep

    density = bertha.get_density_on_grid(grid)

    tosave = numpy.column_stack((grid, density))

    numpy.savetxt("step_%d.txt"%(j), tosave)
    """
    
    ene_list.append(numpy.trace(numpy.matmul(D_ti_dt,fockm_ti_dt)))
    
    # update D_ti and Dp_ti for the next step
    # message for debug
    # fo.write('here I update the matrices Dp_ti and D_ti\n')
    D_ti = numpy.copy(D_ti_dt)
    Dp_ti = numpy.copy(Dp_ti_dt)

    if debug:
      fo.write('  Trace of Dp_ti %.8f\n' % numpy.trace(Dp_ti).real)
      fo.write('  Trace of ovapm * D_ti  %.8f\n' % numpy.trace(numpy.matmul(ovapm,D_ti)).real)
      fo.flush()
   
    #update fock_mid_backwd for the next step
    fock_mid_backwd = numpy.copy(fock_mid_tmp)
   
    return fock_mid_backwd, D_ti, Dp_ti

##########################################################################################

def run_iterations_from_to (startiter, niter, bertha, args, fock_mid_backwd, dt, \
        dip_mat, C, C_inv, ovapm, ndim, debug, Dp_ti, dip_list, ene_list, weight_list, \
        fo, D_ti, occlist):

    dumpcounter = 0
    
    for j in range(startiter, niter):

        sys.stdout.flush()

        start = time.time()
        cstart = time.process_time() 
    
        fock_mid_backwd, D_ti, Dp_ti = main_loop(j, niter, bertha, 
                args.pulse, args.pulseFmax, args.pulsew, args.propthresh,
                args.pulseS, args.t0, fo, D_ti, 
                fock_mid_backwd, dt, dip_mat, C, C_inv, ovapm, 
                ndim, debug, Dp_ti, dip_list, ene_list, weight_list, args.select)

        if fock_mid_backwd is None:
            return False

        dumpcounter += 1

        if args.dumprestartnum > 0:
            if (dumpcounter == args.dumprestartnum) or \
                    (j == niter-1):
                
                encoder.FLOAT_REPR = lambda o: format(o, '.25E')

                json_data = get_json_data(args, j, niter, ndim, ene_list, 
                        dip_list, D_ti, fock_mid_backwd, dip_mat, C, 
                        C_inv, ovapm, Dp_ti, weight_list)

                with open(args.restartfile, 'w') as fp:
                    json.dump(json_data, fp, sort_keys=True, indent=4)

                dumpcounter = 0

        end = time.time()
        cend = time.process_time() 
   
        if (args.iterations):
            print("Iteration %10d"%j, " of %10d"%(niter-1), " ( %15.5f"%(end - start), \
                    " (CPU time: %15.5f"%(cend - cstart), ") s )")
        else:
            rtutil.progress_bar(j, niter-1)

        sys.stdout.flush()

    sys.stdout.flush()

    print("")
    print("")
    print("Dump density density.cube")
    bertha.density_to_cube(D_ti.T, "density.cube", margin = 5.0)
    
    print("Done")
    
    sys.stdout.flush()

    if debug:
        fo.close()
    
    t_point = numpy.linspace(0.0, niter*dt, niter+1)
    numpy.savetxt('dipole.txt',numpy.c_[t_point.real,numpy.array(dip_list).real], fmt='%.12e')
    numpy.savetxt('ene.txt',numpy.c_[t_point.real,numpy.array(ene_list).real], fmt='%.12e')
    if (args.pulse == "analytic"):
        if (occlist[0] != -2):
            numpy.savetxt('weighted_dip.txt', \
                    numpy.c_[t_point.real,numpy.array(weight_list).real], \
                    fmt='%.12e')
    #print(numpy.array(weight_list).shape)
    bertha.finalize()

    return True

##########################################################################################

def restart_run(args):
    
    fp = open(args.restartfile, 'r')
    json_data = json.load(fp)
    fp.close()

    jstart = int(json_data["j"])
    ndim = int(json_data["ndim"])
    niter = int(json_data["niter"])

    fock_mid_backwd_REAL = numpy.float_(json_data["fock_mid_backwd_REAL"])
    fock_mid_backwd_IMAG = numpy.float_(json_data["fock_mid_backwd_IMAG"])
    dip_mat_REAL = numpy.float_(json_data["dip_mat_REAL"])
    dip_mat_IMAG = numpy.float_(json_data["dip_mat_IMAG"])
    C_REAL = numpy.float_(json_data["C_REAL"])
    C_IMAG = numpy.float_(json_data["C_IMAG"])
    C_inv_REAL = numpy.float_(json_data["C_inv_REAL"])
    C_inv_IMAG = numpy.float_(json_data["C_inv_IMAG"])
    ovapm_REAL = numpy.float_(json_data["ovapm_REAL"])
    ovapm_IMAG = numpy.float_(json_data["ovapm_IMAG"])
    Dp_ti_REAL = numpy.float_(json_data["Dp_ti_REAL"])
    Dp_ti_IMAG = numpy.float_(json_data["Dp_ti_IMAG"])
    D_ti_REAL = numpy.float_(json_data["D_ti_REAL"])
    D_ti_IMAG = numpy.float_(json_data["D_ti_IMAG"])
 
    fock_mid_backwd = check_and_covert (fock_mid_backwd_REAL ,
            fock_mid_backwd_IMAG, ndim)
    dip_mat = check_and_covert (dip_mat_REAL, dip_mat_IMAG, ndim)
    C = check_and_covert (C_REAL, C_IMAG, ndim)
    C_inv = check_and_covert (C_inv_REAL, C_inv_IMAG, ndim)
    ovapm = check_and_covert (ovapm_REAL, ovapm_IMAG, ndim)
    Dp_ti = check_and_covert (Dp_ti_REAL, Dp_ti_IMAG, ndim)
    D_ti = check_and_covert (D_ti_REAL, D_ti_IMAG, ndim)

    ene_list_REAL = numpy.float_(json_data["ene_list_REAL"])
    ene_list_IMAG = numpy.float_(json_data["ene_list_IMAG"])
    dip_list_REAL = numpy.float_(json_data["dip_list_REAL"])
    dip_list_IMAG = numpy.float_(json_data["dip_list_IMAG"])

    if ((ene_list_REAL.shape != ene_list_IMAG.shape) or 
        (dip_list_REAL.shape != dip_list_IMAG.shape) or
        (ene_list_REAL.shape != dip_list_IMAG.shape)):
        return False

    ene_list = []
    dip_list = []

    for i in range(ene_list_REAL.shape[0]):
        ene_list.append(numpy.complex128(complex(ene_list_REAL[i],
            ene_list_IMAG[i])))
        dip_list.append(numpy.complex128(complex(dip_list_REAL[i],
            dip_list_IMAG[i])))

    weight_list_REAL = numpy.float_(json_data["weight_list_REAL"])
    weight_list_IMAG = numpy.float_(json_data["weight_list_IMAG"])

    weight_list = []

    if len(weight_list_REAL) != len(weight_list_IMAG):
        print("ERROR in weight_list size")
        exit(1)

    for i in range(len(weight_list_REAL)):
        if len(weight_list_REAL[i]) != len(weight_list_IMAG[i]):
            print("ERROR in weight_list size")
            exit(1)

        row = numpy.zeros(len(weight_list_REAL[i]), dtype=numpy.complex128)

        for j in range(len(weight_list_REAL[i])):
            row[j] = numpy.complex128(complex(weight_list_REAL[i][j],
                weight_list_IMAG[i][j]))

        weight_list.append(row)

    args.pulse = json_data['pulse']
    args.pulseFmax = json_data['pulseFmax']
    args.pulsew = json_data['pulsew'] 
    args.dt = json_data['dt']
    args.inputfile = json_data['inputfile']
    args.fittfile = json_data["fittfile"]
    args.fitcoefffile = json_data["fitcoefffile"]
    args.vctfile = json_data["vctfile"]
    args.ovapfile = json_data["ovapfile"]
    args.dumpfiles = json_data["dumpfiles"]
    args.direction = json_data["direction"]
    
    totaltime = json_data["totaltime"]
    if args.totaltime < totaltime:
        args.totaltime = json_data["totaltime"]
    else:
        niter = int(args.totaltime/args.dt)

    args.debug = json_data["debug"]
    args.verbosity = json_data["verbosity"]
    args.iterations = json_data["iterations"]
    args.select = json_data["select"]
    args.select_pert = json_data["select_pert"]
    args.propthresh = json_data["propthresh"]
    args.thresh = json_data["thresh"]
    args.wrapperso = json_data["wrapperso"]
    args.dumprestartnum = json_data["dumprestartnum"]
    args.pulseS = json_data["pulseS"]
    args.t0 = json_data["t0"]

    molist = args.select.split("&")
    occlist = molist[0].split(";")
    occlist = [int(m) for m in occlist]
    virtlist = molist[1].split(";")
    virtlist = [int(m) for m in virtlist]

    print("Options: ")
    print(args) 
    print("")
    print("")

    if not os.path.isfile(args.wrapperso):
        print("SO File ", args.wrapperso, " does not exist")
        return False
    
    bertha = berthamod.pybertha(args.wrapperso)

    # TODO to remove full run 
    ovapm, eigem, fockm, eigen = single_point (args, bertha)
    if ovapm is None:
        return False

    ndim = bertha.get_ndim()
    nshift = bertha.get_nshift()
    nocc = bertha.get_nocc()
 
    bertha.realtime_init()
    
    print("Start RT")
    
    debug = args.debug
    dt = args.dt
    t_int = args.totaltime
    
    print("Debug: ", debug)
    print("dt : ", dt)
    print("Total time  : ", t_int)
    print("Number of iterations: ", niter)
    
    sys.stdout.flush()
    
    fo = sys.stderr
    if debug:
        fo = open("debug_info.txt", "w")
    
    return run_iterations_from_to (jstart+1, niter, bertha, args, fock_mid_backwd, \
            dt, dip_mat, C, C_inv, ovapm, ndim, debug, Dp_ti, dip_list, ene_list, \
            weight_list, fo, D_ti, occlist)
 
##########################################################################################

def normal_run(args):

    print("Options: ")
    print(args) 
    print("")
    print("")

    if not os.path.isfile(args.wrapperso):
        print("SO File ", args.wrapperso, " does not exist")
        return False
    
    bertha = berthamod.pybertha(args.wrapperso)

    ovapm, eigem, fockm, eigen = single_point (args, bertha)
    if ovapm is None:
        return False

    ndim = bertha.get_ndim()
    nshift = bertha.get_nshift()
    nocc = bertha.get_nocc()
    
    bertha.realtime_init()
    
    print("Start RT")
    
    debug = args.debug
    dt = args.dt
    t_int = args.totaltime
    niter = int(t_int/dt)
    
    print("Debug: ", debug)
    print("dt : ", dt)
    print("Total time  : ", t_int)
    print("Number of iterations: ", niter)
    
    sys.stdout.flush()
    
    ene_list = []
    dip_list = []
    weight_list = []
    imp_list = []
    Enuc_list = []
    
    C = eigem
    D_0 = numpy.zeros((ndim,ndim), dtype=numpy.complex128)
    for num in range(nocc):
        D_0[num+nshift,num+nshift]=1.0+0.0j
    
    fo = sys.stderr
    if debug:
        fo = open("debug_info.txt", "w")
    
    #print type(eigem)
    #C_inv used to backtransform D(AO)
    
    try: 
        C_inv = numpy.linalg.inv(eigem)
    except LinAlgError:
        print("Error in numpy.linalg.inv of eigem") 
        return False
    
    if debug:
      test = numpy.matmul(C_inv, eigem)
      fo.write("Check if atol is 1.e-14 for inversion of C: %s\n"% \
            numpy.allclose(numpy.eye((ndim),dtype=numpy.complex128), \
            test,atol=1.e-14))
    
    if debug:
      diff = test - numpy.eye((ndim),dtype=numpy.complex128)
      mdiff = numpy.max(diff)
      fo.write("  maxdiff is: %.12e %.12e\n"%(mdiff.real, mdiff.imag))
    
    if debug:
      test = numpy.matmul(numpy.conjugate(C.T),numpy.matmul(ovapm,C))
      fo.write("Check orthonormal orbitals (atol = 1.e-14): %s\n"% \
            numpy.allclose(numpy.eye((ndim),dtype=numpy.complex128), \
            test, atol=1.e-14))
      diff = test - numpy.eye((ndim),dtype=numpy.complex128)
      mdiff = numpy.max(diff)
      fo.write("  maxdiff is: %.12e %.12e\n"%(mdiff.real, mdiff.imag))
    
    #build density in ao basis
    
    occeigv = numpy.zeros((ndim,nocc), dtype=numpy.complex128)
    iocc = 0
    
    for i in range(ndim):
        if i >= nshift and iocc < nocc:
            for j in range(ndim):
                occeigv[j, iocc] = eigem[j, i]
            iocc = iocc + 1
    
    Da = numpy.matmul(occeigv,numpy.conjugate(occeigv.transpose()))
    
    print("")
    print("Dump ground state density density0.cube")
    bertha.density_to_cube(Da.T, "density0.cube",margin=5.0)
    
    print("Done")

    if debug:
      #check trace(S Da)
      trace_ds = numpy.trace(numpy.matmul(Da,ovapm))
      trace_dsfock = numpy.trace(numpy.matmul(Da,fockm))
      fo.write("Density matrix trace at  t0: %.12e %.12e \n"%(trace_ds.real,trace_ds.imag))
      fo.write("Trace of fock*density at t0: %.12e %.12e \n"%(trace_dsfock.real, trace_dsfock.imag))
    
    normalise = 1
    dip_mat = None

    dipx_mat, dipy_mat, dipz_mat = \
            bertha.get_realtime_dipolematrix (0, normalise)
    
    if (args.direction == 4):
        dip_mat = dipx_mat
    elif (args.direction == 3):
        dip_mat = dipy_mat
    elif (args.direction == 2):
        dip_mat = dipz_mat
    
    if debug:
      fockmh = numpy.conjugate(fockm.T)
      diff_fockmh = fockm-fockmh
      mdiff = numpy.max(diff_fockmh)
      fo.write("Check max diff fockm-fockmh: %.12e %.12e\n"%\
            (mdiff.real, mdiff.imag))
      fo.write("Fockm (t=0) is hermitian: %s \n"%numpy.allclose(fockm,fockmh,atol=1.e-15))
    
    molist = args.select.split("&")
    occlist = molist[0].split(";")
    occlist = [int(m) for m in occlist]
    virtlist = molist[1].split(";")
    virtlist = [int(m) for m in virtlist]
    
    if (args.pulse == "analytic"):
        Amp=args.pulseFmax

        # to check 
        dipz_mo=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(dip_mat,C))

        if args.select_pert:
            dipz_mo=rtutil.dipole_selection(dipz_mo,nshift,nocc,occlist,virtlist,fo,debug)

        print(" Perturb with analytic kick ")
        u0=rtutil.exp_opmat(dipz_mo,numpy.float_(-Amp),debug,fo)
        Dp_init=numpy.matmul(u0,numpy.matmul(D_0,numpy.conjugate(u0.T)))
        #transform back Dp_int
        Da=numpy.matmul(C,numpy.matmul(Dp_init,numpy.conjugate(C.T)))
        D_0=Dp_init

    print("Start first mo_fock_mid_forwd_eval ")
    
    fock_mid_init = rtutil.mo_fock_mid_forwd_eval(bertha,Da,fockm,0,numpy.float_(dt),\
            dip_mat,C,C_inv,ovapm,ndim, debug, fo, args.pulse, args.pulseFmax, args.pulsew, args.t0, args.pulseS, 
            args.propthresh)
    
    if (fock_mid_init is None):
        print("Error accurs in mo_fock_mid_forwd_eval")
        return False
    
    if debug:
      fock_mid_h=numpy.conjugate(fock_mid_init.T)
      diff_fock_mid_h=fock_mid_init-fock_mid_h
      fo.write("Max diff fock_mid_init-fock_mid_h"%(numpy.max(diff_fock_mid_h)))
      fo.write('Fockm (t=1/2) is hermitian: %s\n'% 
            numpy.allclose(fock_mid_init,fock_mid_h,atol=1.e-14))
    
    fockp_mid_init=numpy.matmul(numpy.conjugate(C.T),numpy.matmul(fock_mid_init,C))
    u=rtutil.exp_opmat(fockp_mid_init,numpy.float_(dt),debug,fo)
    #u=rtutil.exp_opmat(fockp_mid_init,numpy.float_(dt),debug,fo)
    #u=scila.expm(-1.j*fockp_mid_init*dt)
    temp=numpy.matmul(D_0,numpy.conjugate(u.T))
    Dp_t1=numpy.matmul(u,temp)
    
    #check u if unitary
    if debug:
      test_u=numpy.matmul(u,numpy.conjugate(u.T))
      fo.write('U is unitary : %s' % 
            numpy.allclose(test_u,numpy.eye(u.shape[0]),atol=1.e-14))
    
    #backtrasform Dp_t1
    D_t1=numpy.matmul(C,numpy.matmul(Dp_t1,numpy.conjugate(C.T)))
    
    if debug:
      diff = D_t1 - Da 
      mdiff = numpy.max(diff)
      fo.write("Max diff density: %.12e %.12e \n"%(mdiff.real, mdiff.imag))
    
    dip_list.append(numpy.trace(numpy.matmul(Da,dip_mat)))
    dip_list.append(numpy.trace(numpy.matmul(D_t1,dip_mat)))
    if (args.pulse == "analytic"):
      if (occlist[0] != -2):
          #dipoleanalysis
          res=rtutil.dipoleanalysis(dipz_mo,D_0,occlist,virtlist,nshift,fo,debug)
          weight_list.append(res)            
          res=rtutil.dipoleanalysis(dipz_mo,Dp_t1,occlist,virtlist,nshift,fo,debug)
          weight_list.append(res)            
    if debug:                              
      fo.write("Dipole: %.12e\n"%(numpy.trace(numpy.matmul(Da,dip_mat))).real)
      fo.write("Dipole: %.12e\n"%(numpy.trace(numpy.matmul(D_t1,dip_mat))).real)
    
    if debug:
      tfock = numpy.trace(numpy.matmul(D_t1,fockm))
      fo.write("Trace fock*density t1: %.12e, %.12e\n"%(tfock.real, tfock.imag))
      trace_ds=numpy.trace(numpy.matmul(D_t1,ovapm))
      fo.write(" Traceds: %.12e %.12ei\n" % (trace_ds.real,trace_ds.imag))
    
    Ndip_z=0.0
    #estrarre le 3 componenti del dipolo nucleare
    
    D_ti=D_t1
    Dp_ti=Dp_t1
    #aggiungere repulsione nucleare
    #Enuc_list.append(-func_t0*Ndip_z+Nuc_rep) #just in case of non-zero nuclear dipole
    
    start = time.time()
    cstart = time.process_time() 

    fockm_ti=bertha.get_realtime_fock(D_ti.T)

    end = time.time()
    cend = time.process_time() 

    print("RealTime Fock Time:            %15.5f"%(end-start), " (CPU: %15.5f"%(cend-cstart), " ) s")  
    print("RealTime Fock Time Without Py: %15.5f"%(bertha.get_focktime()), \
            " (CPU: %15.5f"%(bertha.get_fockctime()), " ) s")
    print("")
    ene_list.append(numpy.trace(numpy.matmul(Da,fockm)))
    ene_list.append(numpy.trace(numpy.matmul(D_ti,fockm_ti)))
    
    print("Starting iterations ...")
    print("")
    
    fock_mid_backwd = numpy.copy(fock_mid_init)

    return run_iterations_from_to (1, niter, bertha, args, fock_mid_backwd, \
            dt, dip_mat, C, C_inv, ovapm, ndim, debug, Dp_ti, dip_list, ene_list, \
            weight_list, fo, D_ti, occlist)
    
##########################################################################################

def main():

   listpulses = ""
   for key in rtutil.funcswitcher:
       listpulses += key
       listpulses += " "
   
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
   parser.add_argument("-m", "--dt", help="Specify dt to be used (default: 0.1)", required=False,
           default=0.1, type=numpy.float64)
   parser.add_argument("-T", "--totaltime", help="Specify total time (default: 1.0)", required=False,
           default=1.0, type=numpy.float64)
   parser.add_argument("--pulse", help="Specify the pulse to use [" + listpulses + "] default kick", required=False, 
           type=str, default="kick")
   parser.add_argument("--pulseFmax", help="Specify the pulse Fmax value (default: 0.0001)", 
           default=0.0001, type=numpy.float64)
   parser.add_argument("--pulsew", help="Specify the pulse w value if needed (default: 0.0)", 
           default=0.0, type=numpy.float64)
   parser.add_argument("-d", "--debug", help="Debug on, prints debug info to debug_info.txt", required=False, 
           default=False, action="store_true")
   parser.add_argument("-v", "--verbosity", help="Verbosity level 0 = minim, -1 = print iteration info, " + 
           "1 = maximum (defaul -1)", required=False, default=-1, type=int)
   parser.add_argument("--iterations", help="Use iteration number instead of progressbar",
           required=False, default=False, action="store_true")
   parser.add_argument("--thresh", help="det threshold (default = 1.0e-12)", required=False, 
           type=numpy.float64, default=1.0e-12)
   parser.add_argument("--wrapperso", help="set wrapper SO (default = ../../lib/bertha_wrapper.so)", 
           required=False, type=str, default="../lib/bertha_wrapper.so")
   parser.add_argument("--select", help="Specify  occ. and virt. MO for selective perturbation (default: -2; 0 & 0)", 
           default="-2; 0 & 0", type=str)
   parser.add_argument("--select_pert", help="Turn on selective perturbation mode", required=False,
           default=False, action="store_true")
   parser.add_argument("--propthresh", help="threshold for midpoint iterative scheme (default = 1.0e-6)", required=False, 
           type=numpy.float64, default=1.0e-6)
   parser.add_argument("--t0", help="Specify the center of gaussian enveloped pulse (default: 0.0)", 
           default=0.0, type=numpy.float64)
   parser.add_argument("--pulseS", help="Specify the pulse s-parameter: for gaussian enveloped pulse, " + 
           "'s' is sqrt(variance) (default: 0.0)", default=0.0, type=numpy.float64)
   parser.add_argument("--direction", help="Specify the pulse direction (2 = z, 3 = y, 4 = x) ", 
           default=4, type=int)

   parser.add_argument("--restartfile", help="set a restart file (default: restart_pybertha.json)", 
           required=False, type=str, default="restart_pybertha.json")
   parser.add_argument("--dumprestartnum", help="dump restart file every N iterations (default: -1)",
           required=False, type=int, default=-1)
   parser.add_argument("--restart", help="restart run from file",
           required=False, default=False, action="store_true")

   parser.add_argument("-g","--geometry", help="Specify system (Angstrom) geometry", required=False, 
        type=str, default="")
   parser.add_argument("--obs", \
        help="Specify BERTHA (Active system) basisset \"atomname1:basisset1,atomname2:basisset2,...\"", \
        required=False, type=str, default="")
   parser.add_argument("--fittset", \
        help="Specify BERTHA (Active system) fitting set \"atomname1:fittset1,atomname2:fittset2,...\"", \
        required=False, type=str, default="")
   parser.add_argument("--func", 
	    help="Specify exchange–correlation energy functional for active system available: LDA,B88P86,HCTH93,BLYP (default=BLYP)", \
        type=str, default="BLYP")
   parser.add_argument("-j","--jsonbasisfile", \
        help="Specify BERTHA JSON file for fitting and basis (default: fullsets.json)", \
       required=False, type=str, default="fullsets.json")
   parser.add_argument("--convertlengthunit", help="Specify a length converter [default=1.0]", \
        type=float, default=1.0)
   parser.add_argument("--berthamaxit", help="set bertha maxiterations (default = %d)"%(MAXIT), 
        required=False, type=int, default=MAXIT)

   args = parser.parse_args()

   if (not args.restart):
       if args.totaltime < 0.0:
           args.totaltime = 1.0

       if (not normal_run (args)):
           exit(1)
   else:
      if (not restart_run (args)):
           exit(1)

##########################################################################################

if __name__ == "__main__":
    main()

