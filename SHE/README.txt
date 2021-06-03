python3 full_eda_nocv.py --berthamodpaths "../pybertha;../../berthaingen/pybgen"  --fragA ./SHE/AuCn+/Au+.xyz  \
   --fragB ./SHE/AuCn+/Cn.xyz  --molecule ./SHE/AuCn+/AuCn+.xyz \ 
   --basisset "Au:dyall_vtz,Cn:dyall_vtz" --fittset "Au:b20,Cn:cn" \
   --energyconverter 627.50961 --npairs=12 --convertlengthunit 1.8897259886 --cube | tee AuCn+out.txt 

