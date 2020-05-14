####################################
#Psi4RT (aka EZED): real-time td-hf/dft
#      based on Psi4NumPy
###################################
# Test for perturb-and-evolve approach
# parameters in input : set imp_type as "kick"
import os
import sys
import time
import scipy
import argparse
import numpy as np

import json
import pickle 
from json import encoder


##########################################################################################

def vctto_npcmplxarray (realp, imagp):

    list_REAL = np.float_(realp)
    list_IMAG = np.float_(imagp)

    if len(list_REAL) != len(list_IMAG):
        return None

    rlist = []
    for i in range(len(list_REAL)):
        rlist.append(np.complex128(complex(list_REAL[i], list_IMAG[i])))

    return rlist

##########################################################################################

def mtxto_npcmplxarray (realp, imagp):

    list_REAL = np.float_(realp)
    list_IMAG = np.float_(imagp)

    if len(list_REAL) != len(list_IMAG):
        return None

    rlist = []

    for i in range(len(list_REAL)):
        if len(list_REAL[i]) != len(list_IMAG[i]):
            return None
        
        row = np.zeros(len(list_REAL[i]), dtype=np.complex128)

        for j in range(len(list_REAL[i])):
            row[j] = np.complex128(complex(list_REAL[i][j],
                list_IMAG[i][j]))
                
        rlist.append(row)

    return rlist

##########################################################################################

def is_jsonable(x):
    
    try:
        json.dumps(x)
        return True
    except:
        return False

##########################################################################################

def get_json_data(args, D_ti, fock_mid_backwd, j, dt, H, I, dip_mat, \
        C, C_inv, S, nbf, imp_opts, func, Dp_ti, weighted_dip, dip_list, \
        ene_list, imp_list, dipmo_mat, ndocc, occlist, virtlist, debug, HL, \
        psi4options, geom, do_weighted, Enuc_list, imp_params, calc_params, \
        Ndip_dir, Nuc_rep):

    json_data = {}

    # in such a way some values is duobled to be improved  (see HL or debug)
    for arg in vars(args):
        json_data["args."+arg] = getattr(args, arg)

    othervalues = {
        "j" : j, 
        "Nuc_rep" : Nuc_rep,
        "Ndip_dir" : Ndip_dir,
        "imp_params" : imp_params, 
        "calc_params" : calc_params,
        "do_weighted" : do_weighted,
        "geom" : geom,
        "psi4options": psi4options,
        "dt" : dt, 
        "H" : H.tolist(), 
        "I" : I.tolist(), 
        "dip_mat" : dip_mat.tolist(), 
        "C" : C.tolist(), 
        "C_inv" : C_inv.tolist(), 
        "S" : S.tolist(), 
        "nbf" : nbf, 
        "func" : func,
        "imp_list" : imp_list, 
        "dipmo_mat" : dipmo_mat.tolist(), 
        "ndocc" : ndocc, 
        "occlist" : occlist, 
        "virtlist" : virtlist, 
        "debug" : debug, 
        "HL" : HL,
        "Enuc_list" : Enuc_list,
        "imp_opts" : imp_opts,
        'ene_list_REAL': np.real(ene_list).tolist(),
        'ene_list_IMAG': np.imag(ene_list).tolist(),
        'dip_list_REAL': np.real(dip_list).tolist(),
        'dip_list_IMAG': np.imag(dip_list).tolist(),
        'weighted_dip_REAL': np.real(weighted_dip).tolist(),
        'weighted_dip_IMAG': np.imag(weighted_dip).tolist(),
        'D_ti_REAL' : np.real(D_ti).tolist(),
        'D_ti_IMAG' : np.imag(D_ti).tolist(),
        'fock_mid_backwd_REAL' : np.real(fock_mid_backwd).tolist(),
        'fock_mid_backwd_IMAG' : np.imag(fock_mid_backwd).tolist(),
        'Dp_ti_REAL': np.real(Dp_ti).tolist(),
        'Dp_ti_IMAG': np.imag(Dp_ti).tolist()
        }

    json_data.update(othervalues)

    #check if it is seralizable
    for key in json_data:
        if not is_jsonable(json_data[key]):
            print("This ", key, " is not serializable ", type(json_data[key]), \
                    is_jsonable(json_data[key]))
            print (weighted_dip)

    return json_data
        
####################################################################################

def main_loop (D_ti, fock_mid_backwd, j, dt, H, I, dip_mat, C, C_inv, S, nbf, \
        imp_opts, func, fo, basisset, Dp_ti, weighted_dip, dip_list, ene_list, \
        imp_list, dipmo_mat, ndocc, occlist, virtlist, debug, HL, do_weighted, \
        Enuc_list, Ndip_dir, Nuc_rep):

    J_i,Exc_i,func_ti,F_ti,fock_mid_tmp=util.mo_fock_mid_forwd_eval(D_ti,\
                fock_mid_backwd,j,dt,H,I,dip_mat,C,C_inv,S,nbf,\
                imp_opts,func,fo,basisset)
        
    Ah=np.conjugate(fock_mid_tmp.T)
    fo.write('Fock_mid hermitian: %s\n' % np.allclose(fock_mid_tmp,Ah))
    #transform fock_mid_init in MO basis
    fockp_mid_tmp=np.matmul(np.conjugate(C.T),np.matmul(fock_mid_tmp,C))
    u=util.exp_opmat(np.copy(fockp_mid_tmp),dt)
    #u=scipy.linalg.expm(-1.0j*fockp_mid_tmp*dt)
    #check u is unitary
    test_u=np.matmul(u,np.conjugate(u.T))
    if (not np.allclose(np.eye(u.shape[0]),test_u)):
        print('U is not unitary\n')
    
    #check the trace of density to evolve
    fo.write('tr of density to evolve: %.8f\n' % np.trace(Dp_ti).real)
    
    #evolve the density in orthonormal basis
    temp=np.matmul(Dp_ti,np.conjugate(u.T))
    Dp_ti_dt=np.matmul(u,temp)
    
    #backtransform Dp_ti_dt
    D_ti_dt=np.matmul(C,np.matmul(Dp_ti_dt,np.conjugate(C.T)))
    fo.write('%.8f\n' % np.trace(Dp_ti_dt).real)
    #dipole expectation for D_ti
    dip_list.append(np.trace(np.matmul(dip_mat,D_ti)))
    
    if debug:
        fo.write('Dipole  %.8f %.15f\n' % (j*dt, 2.00*dip_list[j].real))
    
    if (do_weighted == -2):
        #weighted dipole 
        res = util.dipoleanalysis(dipmo_mat,Dp_ti,ndocc,occlist,virtlist,debug,HL)
        weighted_dip.append(res)
    #Energy expectation value at t = t_i 
    Enuc_list.append(-func_ti*Ndip_dir+Nuc_rep) #just in case of non-zero nuclear dipole
    if (func=='hf'):
        ene_list.append(np.trace(np.matmul(D_ti,(H+F_ti))))
    else:
        ene_list.append(2.00*np.trace(np.matmul(D_ti,H))+J_i+Exc_i-np.trace(np.matmul(D_ti,(func_ti*dip_mat))))
    imp_list.append(func_ti)
    
    #update D_ti and Dp_ti for the next step
    
    if debug :
        fo.write('here I update the matrices Dp_ti and D_ti\n')
    
    D_ti=np.copy(D_ti_dt)
    Dp_ti=np.copy(Dp_ti_dt)
    #update fock_mid_backwd for the next step
    fock_mid_backwd=np.copy(fock_mid_tmp)

    return fock_mid_backwd, D_ti, Dp_ti

####################################################################################

def normal_run_init (args):

    cfnames = args.cubefilenames.split(";")
    if len(cfnames) != 4:
        print("Error in --cube-filenames, you need to specify 4 filenames")
        exit(1)

    dbgfnames = args.dbgfnames.split(";")
    if len(dbgfnames) != 2:
        print("Error in --debug-filenames, you need to specify 2 filenames")
        exit(1)

    outfnames = args.outfilenames.split(";")
    if len(outfnames) != 4:
        print("Error in --out-filenames, you need to specify 4 filenames")
        exit(1)
   
    debug = args.debug
    fgeom = args.geom
    basis_set = args.obs
    direction = args.axis
    HL = args.principal
    
    if args.puream:
        use_am = "True"
    else:
        use_am = "False"
    
    molist = args.select.split("&")
    occlist = molist[0].split(";")
    occlist = [int(m) for m in occlist]
    do_weighted = occlist.pop(0)
    virtlist = molist[1].split(";")
    virtlist = [int(m) for m in virtlist]

    psi4options = {'basis': basis_set,
                   'puream': use_am,
                   #'DF_BASIS_SCF' : 'cc-pvqz-jkfit',
                   'dft_radial_scheme' : 'becke',
                   #'dft_radial_points': 49,
                   'dft_spherical_points' : 434,
                   #'dft_spherical_points': 146,     # Often needed
                   #'dft_radial_points': 55,         # Often needed
                   #'dft_radial_scheme': 'treutler',   # Rarely needed
                   #'dft_nuclear_scheme': 'treutler',  # Rarely needed
                   'scf_type': 'direct',
                   'DF_SCF_GUESS': 'False',
                   'cubeprop_tasks': ['density'],
                   'CUBIC_GRID_OVERAGE' : [7.0,7.0,7.0],
                   'CUBEPROP_ISOCONTOUR_THRESHOLD' : 1.0,
                   'CUBIC_GRID_SPACING' : [0.1,0.1,0.1],
                   'e_convergence': 1e-8,
                   'd_convergence': 1e-8}

    svwn5_func = {
        "name": "SVWN5",
        "x_functionals": {
            "LDA_X": {}
        },
        "c_functionals": {
            "LDA_C_VWN": {}
        }
    }
    
    if (do_weighted == -2):
        if debug:
            print("Selected transitions from %s to %s MOs"%(str(occlist), str(virtlist)))
    
    imp_opts, calc_params = util.set_params(args.inputfname)
    
    if imp_opts['imp_type'] == 'analytic' :
        analytic = True
    else:
        analytic = False

    #dt in a.u
    dt = np.float_(calc_params['delta_t'])
    #time_int in atomic unit
    time_int = calc_params['time_int']
    niter = int(time_int/dt)
    ######################################
    #memory
    psi4.set_memory(int(2e9))
    #output
    psi4.core.set_output_file(dbgfnames[1], False)
    #basis set options etc
    psi4.set_options(psi4options)
    #geometry set
    print("Reaading geometry from ", fgeom)
    geom, mol = util.set_input(fgeom)
    numpy_memory = 8
    #build dummy wfn
    mol_wfn = psi4.core.Wavefunction.build( \
            mol,psi4.core.get_global_option('basis'))
    mol_wfn.basisset().print_detail_out()
    #define ndocc and nbeta
    ndocc=mol_wfn.nalpha()
    nbeta=mol_wfn.nbeta()
    if (ndocc != nbeta):
        print('Not close-shell')
    
    #initialize mints object
    mints = psi4.core.MintsHelper(mol_wfn.basisset())
    S = np.array(mints.ao_overlap())
    #get T,V,dipole_z
    T = np.array(mints.ao_kinetic())
    V = np.array(mints.ao_potential())
    dipole=mints.ao_dipole()
    #dipole is  a list of psi4.core.Matrix objects
    # dipole[0] -> x, dipole[1] -> y, dipole[2] -> z
    dip_mat=np.copy(dipole[direction])
    H = T+V
    #internal defined functional 

    #compute ground state wfn (for additional quantities)
    if (calc_params['func_type']=='hf'):
        ene, wfn = psi4.energy('scf',return_wfn=True)
    elif (calc_params['func_type'] == 'svwn5'):
        ene, wfn = psi4.energy('scf',dft_functional=svwn5_func,return_wfn=True)
    else:
        ene, wfn = psi4.energy(calc_params['func_type'],return_wfn=True)
    print("Using ", calc_params['func_type']) 
    psi4.cubeprop(wfn)
    os.rename("Da.cube",cfnames[0])
    os.rename("Db.cube",cfnames[1])
    os.rename("Dt.cube",cfnames[2])
    os.rename("Ds.cube",cfnames[3])
    #C coefficients
    C=np.array(wfn.Ca())  
    
    dipmo_mat = np.matmul(np.conjugate(C.T),np.matmul(dip_mat,C))
    #get a scf alpha density
    Da = np.array(wfn.Da())
    
    ################################################################
    # Get nbf for closed shell molecules
    # ndocc already defined as ndocc=h2o_wfn.nalpha()
    nbf = S.shape[0]
    
    print('\nNumber of occupied orbitals: %d' % ndocc)
    print('Number of basis functions: %d' % nbf)
    
    # Run a quick check to make sure everything will fit into memory
    I_Size = (nbf**4) * 8.e-9
    print("\nSize of the ERI tensor will be %4.2f GB." % I_Size)
    
    # Estimate memory usage
    memory_footprint = I_Size * 1.5
    if I_Size > numpy_memory:
        psi4.core.clean()
        raise Exception("Estimated memory utilization (%4.2f GB) " +\
                "exceeds numpy_memory limit of %4.2f GB." % (memory_footprint, numpy_memory))
    #Get Eri (2-electron repulsion integrals)
    I = np.array(mints.ao_eri())
    
    #D_0 is the density in the reference molecular basis
    D_0=np.zeros((nbf,nbf))
    for num in range(int(ndocc)):
        D_0[num,num]=1.0
    #nuclear dipole for non-homonuclear molecules
    Ndip = mol.nuclear_dipole()
    Ndip_dir = Ndip[direction]
    #for the time being the nuclear dipole contribution to the dipole and energy
    # is not considered
    Ndip_dir = 0.0
    
    Nuc_rep = mol.nuclear_repulsion_energy()
    Enuc_list = []
    print("Number of iterations ", niter)
    Dp_0= D_0
    #set the functional type
    if (calc_params['func_type'] == 'svwn5'):
       func = svwn5_func
    else:
       func = calc_params['func_type']
    # the basisset object
    basisset = mol_wfn.basisset()
    #print("analytic : %i" % analytic)
    if (analytic):
        print('Perturb density with analytic delta')
        # set the perturbed density -> exp(-ikP)D_0exp(+ikP)
        k = imp_opts['Fmax']
       
        #dip_mat is transformed to the reference MO basis
        dip_mo=np.matmul(np.conjugate(C.T),np.matmul(dip_mat,C))
        u0 = util.exp_opmat(dip_mo,np.float_(-k))
        Dp_init= np.matmul(u0,np.matmul(Dp_0,np.conjugate(u0.T)))
        func_t0 = k
        #backtrasform Dp_init
        D_init = np.matmul(C,np.matmul(Dp_init,np.conjugate(C.T)))
        Da = D_init
        Dp_0 = Dp_init 
       
        #J0p,Exc0p,F_t0=util.get_Fock(D_ti,H,I,func,basisset)
        #if (func == 'hf'):                                  
        #    testene = np.trace(np.matmul(D_init,(H+F_t0)))  
        #else:                                               
        #    testene = 2.00*np.trace(np.matmul(D_init,H))+J0p+Exc0p
        #print('trace D(0+): %.8f' % np.trace(Dp_init).real)       
        #print(testene+Nuc_rep)                                    

    fo = open(dbgfnames[0], "w")
    
    #containers
    ene_list = []
    dip_list = []
    imp_list=[]
    weighted_dip = []
    #for molecules with permanent nuclear dipole add Enuc_list[k] to ene
    #note that : i.e in the case of linear symmetric molecule the nuclear dipole can be zero depending
    #on the choice of reference coordinates system
    
    #C_inv used to backtransform D(AO)
    
    C_inv=np.linalg.inv(C)
    print('Entering in the first step of propagation')
    J0,Exc0,func_t0,F_t0,fock_mid_init=util.mo_fock_mid_forwd_eval(Da,wfn.Fa(),\
            0,dt,H,I,dip_mat,C,C_inv,S,nbf,imp_opts,func,fo,basisset)
    
    #check the Fock
    if debug :
         print('F_t0 is equal to wfn.Fa() : %s' % np.allclose(wfn.Fa(),F_t0,atol=1.0e-12))
    
    #check hermicity of fock_mid_init
    
    Ah=np.conjugate(fock_mid_init.T)
    fo.write('Fock_mid hermitian: %s\n' % np.allclose(fock_mid_init,Ah))
    
    #propagate D_t0 -->D(t0+dt)
    #
    #fock_mid_init is transformed in the MO ref basis
    fockp_mid_init = np.matmul(np.conjugate(C.T),np.matmul(fock_mid_init,C))
    
    #u=scipy.linalg.expm(-1.j*fockp_mid_init*dt)
    u = util.exp_opmat(fockp_mid_init,dt)
    
    temp = np.matmul(Dp_0,np.conjugate(u.T))
    
    Dp_t1 = np.matmul(u,temp)
    
    #check u if unitary
    test_u = np.matmul(u,np.conjugate(u.T))
    fo.write('U is unitary :%s\n' % np.allclose(test_u,np.eye(u.shape[0])))
    
    fock_mid_backwd = np.copy(fock_mid_init)
    
    #backtrasform Dp_t1
    
    D_t1=np.matmul(C,np.matmul(Dp_t1,np.conjugate(C.T)))
    
    if (func == 'hf'):
        ene_list.append(np.trace(np.matmul(Da,(H+F_t0))))
    else:    
        ene_list.append(2.00*np.trace(np.matmul(Da,H))+J0+Exc0-np.trace(np.matmul(Da,(func_t0*dip_mat))))
    
    dip_list.append(np.trace(np.matmul(Da,dip_mat)))
    
    #weighted dipole
    if (do_weighted == -2):
        res = util.dipoleanalysis(dipmo_mat,Dp_0,ndocc,occlist,virtlist,debug,HL)
        weighted_dip.append(res)
    
    fock_mid_backwd=np.copy(fock_mid_init) #prepare the fock at the previous midpint
    D_ti=D_t1
    Dp_ti=Dp_t1
    Enuc_list.append(-func_t0*Ndip_dir+Nuc_rep) #just in case of non-zero nuclear dipole
    #
    imp_list.append(func_t0)
    if debug :  
        #trace of D_t1
        fo.write('%.8f\n' % np.trace(Dp_ti).real)
        fo.write('Trace of DS %.8f\n' % np.trace(np.matmul(S,D_ti)).real)
        fo.write('Trace of SD.real %.14f\n' % np.trace(np.matmul(S,D_ti.real)))
        fo.write('Trace of SD.imag %.14f\n' % np.trace(np.matmul(S,D_ti.imag)))
        fo.write('Dipole %.8f %.15f\n' % (0.000, 2.00*dip_list[0].real))
        
    return D_ti, fock_mid_backwd, dt, H, I, dip_mat, C, C_inv, S, nbf, \
            imp_opts, func, fo, basisset, Dp_ti, weighted_dip, dip_list, \
            ene_list, imp_list, dipmo_mat, ndocc, occlist, virtlist, debug, \
            HL, Ndip_dir, Nuc_rep, niter, do_weighted, Enuc_list, psi4options, \
            geom, outfnames, wfn, imp_params, calc_params

####################################################################################

def restart_init (args):

    fp = open(args.restartfile, 'r')
    json_data = json.load(fp)
    fp.close()

    # in such a way some values is duobled to be improved  (see HL or debug)
    #for arg in vars(args):
    #    if arg != "restart" and arg != "debug":
    #        toexe = "args." + arg + " = json_data[\"args."+arg+"\"]" 
    #        eval(toexe)

    args.debug = json_data["args.debug"]
    args.dbgfnames = json_data["args.dbgfnames"]
    args.principal = json_data["args.principal"]
    args.axis = json_data["args.axis"]
    args.select = json_data["args.select"]
    args.geom = json_data["args.geom"]
    args.obs = json_data["args.obs"]
    args.puream = json_data["args.puream"]
    args.psi4root = json_data["args.psi4root"]
    args.psi4basis = json_data["args.psi4basis"]
    args.inputfname = json_data["args.inputfname"]
    args.cubefilenames = json_data["args.cubefilenames"]
    args.outfilenames = json_data["args.outfilenames"]
    args.iterations = json_data["args.iterations"]
    args.restartfile = json_data["args.restartfile"]
    args.dumprestartnum = json_data["args.dumprestartnum"]
    #args.restart 

    fock_mid_backwd = None
    weighted_dip = None
    do_weighted = None
    calc_params = None
    imp_params = None  
    dipmo_mat = None
    Enuc_list = None
    outfnames = None
    Ndip_dir = None
    dip_list = None
    ene_list = None
    imp_list = None
    virtlist = None
    imp_opts = None
    Nuc_rep = None
    occlist = None
    dip_mat = None
    ndocc = None
    debug = None
    C_inv = None
    Dp_ti = None
    niter = None 
    func = None
    D_ti = None
    geom = None 
    nbf = None
    wfn = None 
    HL = None
    fo = None
    dt = None
    H = None
    I = None
    C = None
    S = None

    j = None 

    j            = json_data["j"]
    imp_params   = json_data["imp_params"] 
    calc_params  = json_data["calc_params"]
    do_weighted  = json_data["do_weighted"]
    geom         = json_data["geom"] 
    psi4options  = json_data["psi4options"]
    dt           = json_data["dt"]
    nbf          = json_data["nbf"] 
    func         = json_data["func"] 
    imp_list     = json_data["imp_list"] 
    ndocc        = json_data["ndocc"] 
    occlist      = json_data["occlist"] 
    virtlist     = json_data["virtlist"] 
    debug        = json_data["debug"]
    HL           = json_data["HL"] 
    Enuc_list    = json_data["Enuc_list"] 
    imp_opts     = json_data["imp_opts"] 
    Ndip_dir     = json_data["Ndip_dir"]
    Nuc_rep      = json_data["Nuc_rep"]

    imp_opts_new, calc_params_new = \
            util.set_params(args.inputfname)

    time_int = calc_params_new['time_int']
    niter = int(time_int/dt)

    dipmo_mat = np.asarray(json_data["dipmo_mat"])
    H = np.asarray(json_data["H"])
    I = np.asarray(json_data["I"])
    dip_mat = np.asarray(json_data["dip_mat"])
    C = np.asarray(json_data["C"])
    C_inv = np.asarray(json_data["C_inv"])
    S = np.asarray(json_data["S"])

    ene_list = vctto_npcmplxarray (json_data["ene_list_REAL"], \
            json_data["ene_list_IMAG"])
    dip_list = vctto_npcmplxarray (json_data["dip_list_REAL"], \
            json_data["dip_list_IMAG"])

    weighted_dip = mtxto_npcmplxarray (json_data["weighted_dip_REAL"], \
            json_data["weighted_dip_IMAG"])

    D_ti = np.array(mtxto_npcmplxarray (json_data["D_ti_REAL"], \
            json_data["D_ti_IMAG"]))
    Dp_ti = np.array(mtxto_npcmplxarray (json_data["Dp_ti_REAL"], \
            json_data["Dp_ti_IMAG"]))
    fock_mid_backwd = np.array(mtxto_npcmplxarray (json_data["fock_mid_backwd_REAL"], \
            json_data["fock_mid_backwd_IMAG"]))

    cfnames = args.cubefilenames.split(";")
    if len(cfnames) != 4:
        print("Error in --cube-filenames, you need to specify 4 filenames")
        exit(1)

    dbgfnames = args.dbgfnames.split(";")
    if len(dbgfnames) != 2:
        print("Error in --debug-filenames, you need to specify 2 filenames")
        exit(1)

    outfnames = args.outfilenames.split(";")
    if len(outfnames) != 4:
        print("Error in --out-filenames, you need to specify 4 filenames")
        exit(1)

    fo = open(dbgfnames[0], "w")

    return D_ti, fock_mid_backwd, dt, H, I, dip_mat, C, C_inv, S, nbf, imp_opts, func, fo, \
           Dp_ti, weighted_dip, dip_list, ene_list, imp_list, dipmo_mat, ndocc, occlist, \
           virtlist, debug, HL, Ndip_dir, Nuc_rep, niter, do_weighted, Enuc_list, psi4options, \
           geom, outfnames, wfn, imp_params, calc_params, j

####################################################################################

if __name__ == "__main__":

    ####################################
    # parse arguments from std input
    ####################################
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-d", "--debug", help="Debug on, prints extra debug info to the first --debug-filename", required=False,
            default=False, action="store_true")
    parser.add_argument("--debug-filenames", help="Debug filename [default=\"err.txt;test.dat\"]", required=False,
            default="err.txt;test.dat", type=str, dest="dbgfnames")
    parser.add_argument("-p", "--principal", help="Turn on homo-lumo weighted dipole moment", required=False,
            default=False, action="store_true")
    parser.add_argument("-a", "--axis", help="The axis of  electric field direction (x = 0, y = 1, z = 2, default 2)",
            default=2, type = int)
    parser.add_argument("--select", help="Specify the occ-virt MO weighted dipole moment. " + \
            "(-2; 0 & 0 to activate) (default: 0; 0 & 0)",
            default="0; 0 & 0", type=str)
    parser.add_argument("-g","--geom", help="Specify geometry file", required=False, 
            type=str, default="geom.xyz")
    parser.add_argument("-o","--obs", help="Specify the orbital basis set", required=False, 
            type=str, default="cc-pvdz")
    parser.add_argument("--puream", help="Use pure am basis set", required=False, 
            default=False, action="store_true" )
    parser.add_argument("--psi4root", help="Add psi4 rootdir if needed", required=False, 
            default="", type=str)
    parser.add_argument("--pybertharoot", help="Add pybertha rootdir if needed", required=False, 
            default="", type=str)
    parser.add_argument("--psi4basis", help="Add psi4 basis set rootdir", required=False, 
            default="/home/redo/anaconda3/pkgs/psi4-1.3.2+ecbda83-py37h31b3128_0/share/psi4/basis", type=str)
    parser.add_argument("--input-param-file", help="Add input parameters filename [default=\"input.inp\"]", 
            required=False, default="input.inp", type=str, dest='inputfname')
    parser.add_argument("--cube-filenames", help="Specify cube filenames " + \
            "[default=\"Da0.cube;Db0.cube;Dt0.cube;Ds0.cube\"]", \
            required=False, default="Da0.cube;Db0.cube;Dt0.cube;Ds0.cube", type=str, dest='cubefilenames')
    parser.add_argument("--out-filenames", help="Specify cout filenames " + \
            "[default=\"dipole.txt;imp.txt;ene.txt;weighteddip.txt\"]", \
            required=False, default="dipole.txt;imp.txt;ene.txt;weighteddip.txt", type=str, dest='outfilenames')
    parser.add_argument("--iterations", help="Use iteration number instead of progressbar",
            required=False, default=False, action="store_true")
    parser.add_argument("--restartfile", help="set a restart file (default: restart_psi4rt.json)", 
            required=False, type=str, default="restart_psi4rt.json")
    parser.add_argument("--dumprestartnum", help="dump restart file every N iterations (default: -1)",
            required=False, type=int, default=-1)
    parser.add_argument("--restart", help="restart run from file",
            required=False, default=False, action="store_true")
 
    args = parser.parse_args()

    fock_mid_backwd = None
    weighted_dip = None
    do_weighted = None
    calc_params = None
    imp_params = None  
    dipmo_mat = None
    Enuc_list = None
    outfnames = None
    Ndip_dir = None
    basisset = None
    dip_list = None
    ene_list = None
    imp_list = None
    virtlist = None
    imp_opts = None
    Nuc_rep = None
    occlist = None
    dip_mat = None
    ndocc = None
    debug = None
    C_inv = None
    Dp_ti = None
    niter = None 
    func = None
    D_ti = None
    geom = None 
    nbf = None
    wfn = None 
    HL = None
    fo = None
    dt = None
    H = None
    I = None
    C = None
    S = None

    if args.psi4root != "":
        sys.path.append(args.psi4root)
    elif "PSI4ROOT" in os.environ:
        sys.path.append(os.environ['PSI4ROOT'])

    if args.pybertharoot != "":
        sys.path.insert(0, args.pybertharoot+"/src")
    elif "PYBERTHAROOT" in os.environ:
        sys.path.append(os.environ['PYBERTHAROOT']+"/src")
        
    import rtutil
    import psi4
    import util
    
    os.environ['PSIPATH'] = args.psi4basis
    sys.path.append(os.environ['PSIPATH'])

    jstart = None 

    if args.restart:
        D_ti, fock_mid_backwd, dt, H, I, dip_mat, \
                C, C_inv, S, nbf, imp_opts, func, fo, \
                Dp_ti, weighted_dip, dip_list, \
                ene_list, imp_list, dipmo_mat, ndocc, occlist, \
                virtlist, debug, HL, Ndip_dir, Nuc_rep, niter, do_weighted, \
                Enuc_list, psi4options, geom, outfnames, wfn, imp_params, calc_params, j \
                = restart_init(args)

        print(j, niter)

        mol = psi4.geometry(geom)
        psi4.set_options(psi4options)
        mol_wfn = psi4.core.Wavefunction.build( \
                mol,psi4.core.get_global_option('basis'))
        basisset = mol_wfn.basisset()

        jstart = j
    else:
        D_ti, fock_mid_backwd, dt, H, I, dip_mat, \
                C, C_inv, S, nbf, imp_opts, func, fo, \
                basisset, Dp_ti, weighted_dip, dip_list, \
                ene_list, imp_list, dipmo_mat, ndocc, occlist, \
                virtlist, debug, HL, Ndip_dir, Nuc_rep, niter, do_weighted, \
                Enuc_list, psi4options, geom, outfnames, wfn, imp_params, calc_params \
                = normal_run_init (args)

        jstart = 0

    start = time.time()
    cstart = time.process_time()

    print("Start main iterations \n")
    dumpcounter = 0
    for j in range(jstart+1,niter+1):
        fock_mid_backwd, D_ti, Dp_ti = main_loop (D_ti, fock_mid_backwd, j, dt, \
                H, I, dip_mat, C, C_inv, S, nbf, imp_opts, func, fo, basisset, \
                Dp_ti, weighted_dip, dip_list, ene_list, imp_list, dipmo_mat, \
                ndocc, occlist, virtlist, debug, HL, do_weighted, Enuc_list, \
                Ndip_dir, Nuc_rep)

        dumpcounter += 1

        #if args.dumprestartnum > 0:
        if (dumpcounter == args.dumprestartnum) or \
                (j == niter):
            
            encoder.FLOAT_REPR = lambda o: format(o, '.25E')

            json_data = get_json_data(args, D_ti, \
                    fock_mid_backwd, j, dt, H, I, \
                    dip_mat, C, C_inv, S, nbf, imp_opts, \
                    func, Dp_ti, weighted_dip, dip_list, \
                    ene_list, imp_list, dipmo_mat, \
                    ndocc, occlist, virtlist, debug, HL, \
                    psi4options, geom, do_weighted, Enuc_list, \
                    imp_params, calc_params, Ndip_dir, Nuc_rep)

            with open(args.restartfile, 'w') as fp:
                json.dump(json_data, fp, sort_keys=True, indent=4)

            dumpcounter = 0

        if args.iterations:
            print ("Iter %10d od %10d"%(j,niter))
        else:
            rtutil.progress_bar(j, niter)

    print("")
    print("")
    fo.close()
    end = time.time()
    cend = time.process_time()
    print("Time for %10d time iterations : (%.5f s, %.5f s)\n" %(niter+1,end-start,cend-cstart))
    t_point=np.linspace(0.0,niter*dt,niter+1)
    dip_t=2.00*np.array(dip_list).real + Ndip_dir
    ene_t=np.array(ene_list).real + Nuc_rep
    imp_t=np.array(imp_list)

    print("Dumping output files")
    if (do_weighted == -2):
        wd_dip=2.00*np.array(weighted_dip).real
        np.savetxt(outfnames[3], np.c_[t_point,wd_dip], \
                fmt='%.12e')
    
    np.savetxt(outfnames[0], np.c_[t_point,dip_t], fmt='%.12e')
    np.savetxt(outfnames[1], np.c_[t_point,imp_t], fmt='%.12e')
    np.savetxt(outfnames[2], np.c_[t_point,ene_t], fmt='%.12e')

    if not args.restart:
        wfn.Da().copy(psi4.core.Matrix.from_array(D_ti.real))
        wfn.Db().copy(psi4.core.Matrix.from_array(D_ti.real))
        psi4.cubeprop(wfn)
