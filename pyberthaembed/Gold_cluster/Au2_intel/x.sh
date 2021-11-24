. /home/belp/INSTALL/adf2019.307/adfbashrc.sh


export PYBERTHA_MOD_PATH="/home/belp/EMBEDDING/pybertha/pyemb;/home/belp/EMBEDDING/xcfun/build/lib/python;/home/belp/EMBEDDING/pybertha/src;/home/belp/EMBEDDING/pyadf/src;/home/belp/EMBEDDING/berthaingen/pybgen;/home/belp/EMBEDDING/pybertha/psi4rt"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/belp/EMBEDDING/pybertha/lib

export outputfile=risudyallDZ_B16_thread3

# Run with embedding without field. Dipole for both free active system and embedded active system is also reported.
echo 'Run without external field'

# Set the number of threads to 1 to ensure serial
ulimit -s 4048576
export KMP_STACKSIZE=3400m

# Set the number of threads to 1 to ensure serial
export OMP_NUM_THREADS=3



python3 pyberthaembed_timing.py -gA au2.xyz -gB w8.xyz --act_fittset "Au:b16" --act_obs "Au:dyall_vdz" --convertlengthunit=1.8897259886 --act_func LXCPBE  --env_obs DZ --wrapperso /home/belp/EMBEDDING/pybertha/lib/bertha_wrapper.so -v 1  > $outputfile 2>&1

grep 'CPU'                         $outputfile >  $outputfile.timing
grep 'Time for onel_pot_from_grid' $outputfile >> $outputfile.timing
grep 'Memory usage'                $outputfile >> $outputfile.timing

