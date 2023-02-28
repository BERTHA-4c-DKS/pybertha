
cd Water1

echo 'Run pybertha in Water1'
python3 /home/belp/BERTHA/pybertha/pybertha/pybertha.py --wrapperso=/home/belp/BERTHA/pybertha/lib/bertha_wrapper.so --berthamodpath=/home/belp/BERTHA/pybertha/src --thresh=0.00000000001 --eda_nocv_info  

cd ..
cd Water2

echo 'Run pybertha in Water2'
python3 /home/belp/BERTHA/pybertha/pybertha/pybertha.py --wrapperso=/home/belp/BERTHA/pybertha/lib/bertha_wrapper.so --berthamodpath=/home/belp/BERTHA/pybertha/src --thresh=0.00000000001 --eda_nocv_info
cd ..
cd NOCV

cp ../Water1/info_eda_nocv_fragX.json info_eda_nocv_fragA.json
cp ../Water2/info_eda_nocv_fragX.json info_eda_nocv_fragB.json

echo 'Run pyeda in NOCV'
python3  /home/belp/BERTHA/pybertha/pynocv/py_eda_nocv.py --wrapperso=/home/belp/BERTHA/pybertha/lib/bertha_wrapper.so -np 12 --thresh  0.00000000001 > eda.out

