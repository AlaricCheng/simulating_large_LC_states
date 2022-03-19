
import json

import numpy as np

fdata=np.genfromtxt("../cache/big_circuit3_fidelity.txt", delimiter=",")[:-1]
print('4bit tomo: {:0.3f}'.format(np.mean(fdata[[*range(36),*range(48,48+36)]])))
print('3bit tomo: {:0.3f}'.format(np.mean(fdata[[*range(36,48),*range(48+36,48+48)]])))

with open('../cache/calibration.json') as fid:
    data = json.load(fid)

qs=["q3","q4","q5","q6"]
czs=["cz34","cz54","cz56"]

def pln(arr,num):
    print(' & '.join([('{:0.'+f'{num}'+'f}').format(ai) for ai in arr])+'  \\\\')
def pk(arr,k,num,scale=1):
    pln([data[ai][k]*scale for ai in arr],num)

print('s table 1:')

pk(qs,'f01',3)
pk(qs,'fah',0)
pk(qs,'t1',1)
pk(qs,'t2',2)
pk(qs,'f_00',1,scale=100)
pk(qs,'f_11',1,scale=100)
pk(qs,'rb_avg',2)
pk(qs,'rb_std',2)
pk(czs,'rb_avg',1)
pk(czs,'rb_std',1)
