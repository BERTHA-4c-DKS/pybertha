python3 full_eda_nocv.py --berthamodpaths "../pybertha;../../berthaingen/pybgen"  --fragA ./SHE/AuCn+/Au+.xyz  \
   --fragB ./SHE/AuCn+/Cn.xyz  --molecule ./SHE/AuCn+/AuCn+.xyz \ 
   --basisset "Au:dyall_vtz,Cn:dyall_vtz" --fittset "Au:b20,Cn:cn" \
   --energyconverter 627.50961 --npairs=12 --convertlengthunit 1.8897259886 --cube \
   --lmargin 10.0 --deltax 0.15 --deltay 0.15 --deltaz 0.15 --externalproces | tee AuCn+out.txt 


python3 full_eda_nocv.py --berthamodpaths "../pybertha;../../berthaingen/pybgen"  --fragA ./SHE/AuCN/Au.xyz \
    --fragB ./SHE/AuCN/CN.xyz --molecule ./SHE/AuCN/AuCN.xyz --basisset "C:aug-cc-pVTZ-DK,N:aug-cc-pVTZ-DK,Au:dyall_vtz" \
    --fittset "Au:b20,C:fittset,N:fittset" --energyconverter 627.50961 --npairs=12 --convertlengthunit 1.8897259886 --cube \
    --lmargin 10.0 --deltax 0.15 --deltay 0.15 --deltaz 0.15 --externalproces | tee AuCNout.txt

pymol load_cube.py

python3 pycd.py -f diff_tot.cube
python3 pycd.py -f pair1.cube 
python3 pycd.py -f pair2.cube 
python3 pycd.py -f pair3.cube

gnuplot> plot "diff_tot.cube_cdz.txt" u 1:2 w l lw 4, 
           "pair1.cube_cdz.txt"  u 1:($2*2) w l lw 4, 
            "pair3.cube_cdz.txt" u 1:(2*$2) w l lw 4

