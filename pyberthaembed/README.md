conda activate p4env_oldversion
. /usr/local/adf2019.307/adfbashrc.sh 

python3 pyberthaembed.py -gA H2O.xyz -gB NH3.xyz  \
  --modpaths="/home/redo/BERTHA/pybertha/pyemb;/home/redo/BERTHA/xcfun/build/lib/python;/home/redo/BERTHA/pybertha/src;/home/redo/BERTHA/pyadf/src;/home/redo/BERTHA/berthaingen/pybgen;/home/redo/BERTHA/pybertha/psi4rt" \
  --act_fittset "H:aux4,N:aux4,O:aux4" --act_obs "H:aug-cc-pvdz,O:aug-cc-pvdz,N:aug-cc-pvdz" --convertlengthunit=1.8897259886

can use also PYBERTHA_MOD_PATH env variable to specify all PATHS:

export PYBERTHA_MOD_PATH="/home/redo/BERTHA/pybertha/pyemb;/home/redo/BERTHA/xcfun/build/lib/python;/home/redo/BERTHA/pybertha/src;/home/redo/BERTHA/pyadf/src;/home/redo/BERTHA/berthaingen/pybgen;/home/redo/BERTHA/pybertha/psi4rt"

python3 pyberthaembed.py -gA H2O.xyz -gB NH3.xyz --act_fittset "H:aux4,N:aux4,O:aux4" --act_obs "H:aug-cc-pvdz,O:aug-cc-pvdz,N:aug-cc-pvdz" --convertlengthunit=1.8897259886
