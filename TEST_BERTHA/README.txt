python3 pybgen.py -f fragA.xyz -b "H:basis1,O:basis" -t "H:fit,O:fittset" --convertlengthunit 1.8897259886
python3 pybertha.py  --eda_nocv_info --eda_nocv_frag_file info_eda_nocv_fragA.json

python3 pybgen.py -f fragB.xyz -b "H:basis1,O:basis" -t "H:fit,O:fittset" --convertlengthunit 1.8897259886
python3 pybertha.py  --eda_nocv_info --eda_nocv_frag_file info_eda_nocv_fragB.json

python3 pybgen.py -f full.xyz -b "H:basis1,O:basis" -t "H:fit,O:fittset" --convertlengthunit 1.8897259886
python3 py_eda_nocv.py --energyconverter 627.50961 -np 10

or in a single run:

$ python3 full_eda_nocv.py --berthamodpaths "../pybertha;../../berthaingen/pybgen" \
   --fragA ./TEST_BERTHA/fragA.xyz --fragB ./TEST_BERTHA/fragB.xyz --molecule ./TEST_BERTHA/full.xyz \
    -b "H:avdz,O:avdz" -t "H:fittset,O:fittset" --convertlengthunit 1.8897259886 \
    --energyconverter 627.50961 --npairs=5 | tee out_full 

