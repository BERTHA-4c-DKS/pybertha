Fourier transform based on Padé approximants
-f : input file containing dipole data points
--limits : total number of sampling points to be considered
--frequency : max frequency
python ./pade/pade_transform.py --limit 10000 -f dipole-zz.txt --gamma 3.0e-4 --frequency 20.0 --dw 0.002 -o zz_w.txt > pade.out &
