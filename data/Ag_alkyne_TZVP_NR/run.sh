
cd Ag+

echo 'Run pybertha in Water1'
python3 /home/belp/BERTHA/pybertha/pybertha/pybertha.py --wrapperso=/home/belp/BERTHA/pybertha/lib/bertha_wrapper.so --berthamodpath=/home/belp/BERTHA/pybertha/src --thresh=0.0000000001 --eda_nocv_info  

cd ..
cd ALKYNE

echo 'Run pybertha in Water2'
python3 /home/belp/BERTHA/pybertha/pybertha/pybertha.py --wrapperso=/home/belp/BERTHA/pybertha/lib/bertha_wrapper.so --berthamodpath=/home/belp/BERTHA/pybertha/src --thresh=0.0000000001 --eda_nocv_info
cd ..
cd NOCV

cp ../Ag+/info_eda_nocv_fragX.json info_eda_nocv_fragA.json
cp ../ALKYNE/info_eda_nocv_fragX.json info_eda_nocv_fragB.json

echo 'Run pyeda in NOCV'
python3  /home/belp/BERTHA/pybertha/pynocv/py_eda_nocv.py --wrapperso=/home/belp/BERTHA/pybertha/lib/bertha_wrapper.so -np 12 --thresh  0.0000000001 > eda.out

