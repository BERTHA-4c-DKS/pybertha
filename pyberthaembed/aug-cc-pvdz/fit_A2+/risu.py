import argparse
import os
import sys
import math
import numpy as np

filename = "risu_fde_nofield"
with open(filename,"r") as fp:
    Lines = fp.readlines()
count =0
for ln in Lines:
    count = count + 1
    if ln.startswith('MOLECULAR PROPERTIES of ACTIVE SYSTEM WITHOUT EMBEDDING'):
          linestart = count 
print(linestart)
dipx=float(Lines[linestart+3].split()[3])
dipy=float(Lines[linestart+4].split()[3])
dipz=float(Lines[linestart+5].split()[3])
print("Extracted from file: "+filename)
print('Dipole WITHOUT EMBEDDING (a.u.)')
print(dipx)
print(dipy)
print(dipz)
print(math.sqrt(dipx**2+dipx**2+dipx**2))

ln=''
count = 0
for ln in Lines:
    count = count + 1
    if ln.startswith('MOLECULAR PROPERTIES of ACTIVE SYSTEM IN EMBEDDING:...'):
          linestart = count
print(linestart)
dipx=float(Lines[linestart+3].split()[3])
dipy=float(Lines[linestart+4].split()[3])
dipz=float(Lines[linestart+5].split()[3])
print("Extracted from file: "+filename)
print('Dipole WITH EMBEDDING (a.u.)')
print(dipx)
print(dipy)
print(dipz)
print(math.sqrt(dipx**2+dipx**2+dipx**2))
ln=''
fp.close()


dipx=[]
dipy=[]
dipz=[]

for direction in [0,1,2]:

     for field in [0.001,-0.001]:

        filename = "risu_fde_dir_"+str(field)+"_field_"+str(direction)
        print(filename)
        with open(filename,"r") as fp:
           Lines = fp.readlines()
        
        count =0
        for ll in Lines:
            count = count + 1
            if ll.startswith('MOLECULAR PROPERTIES of ACTIVE SYSTEM IN EMBEDDING:...'):
                   linestart = count
          
        dipx.append(float(Lines[linestart+3].split()[3]))
        dipy.append(float(Lines[linestart+4].split()[3]))
        dipz.append(float(Lines[linestart+5].split()[3]))
        fp.close()

print((dipx[1]-dipx[0])/0.002)
print((dipy[3]-dipy[2])/0.002)
print((dipz[5]-dipz[4])/0.002)

