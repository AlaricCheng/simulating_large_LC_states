import numpy as np
from local_lib import *


def circ1_vector(b,mask, circ1):
    V=np.zeros([1,6])
    for j in range(6):
        V[0,j]=exp4bit(2,j,b,mask,circ1)
    return V


def circ1_tensor(b,mask, circ1):
    coef=(1,1,1/2,-1/2,1/2,-1/2)
    T=np.zeros([6,6])
    for i in range(6):
        for j in range(6):
            T[i,j]=coef[i]*exp4bit(i,j,b,mask,circ1)
    return T


def circ2_vector(b,mask, circ2):
    coef=(1,1,1/2,-1/2,1/2,-1/2)
    V=np.zeros((6,1))
    for i in range(6):
        V[i]=coef[i]*exp3bit(i,b,mask,circ2)
    return V


def bitstr(n):
    return [(bin(i)[2:]).zfill(n)[::-1] for i in range(2**n)]


def compose_on_mask(mask_arr,base, V4, T4, V3): # L = 48 is the num of repititions
    n_qubit=len(mask_arr)
    n_circ=n_qubit//3
    base_index={"XZ":[0,1]*n_circ,"ZX":[1,0]*n_circ}[base]# b is not global
    
    M=[]
    M.append(V4[mask_arr[:3]+"1"][base_index[0]])
    
    for i in range(1,n_circ-1):
        M.append(T4[mask_arr[3*i:3*(i+1)]+"1"][base_index[i]])
    M.append(V3[mask_arr[3*(n_circ-1):3*n_circ]][base_index[n_circ-1]])

    
    return reduce(np.dot,M)



def compose_fidelity(n_qubit, V4, T4, V3, return_avg=True):
    global XZ_masks_cache
    global ZX_masks_cache
    
    try:
        XZ_mask_arr=XZ_masks_cache[f"n_qubit"]
        ZX_mask_arr=ZX_masks_cache[f"n_qubit"]
    except:
        XZ_mask_arr=load_mask_file("./cache/mask_cache/XZ_{0}_mask.csv".format(n_qubit),N=n_qubit)
        ZX_mask_arr=load_mask_file("./cache/mask_cache/ZX_{0}_mask.csv".format(n_qubit),N=n_qubit)
    
    XZ=np.array([compose_on_mask(m, "XZ", V4, T4, V3) for m in XZ_mask_arr])
    ZX=np.array([compose_on_mask(m, "ZX", V4, T4, V3) for m in ZX_mask_arr])
    if return_avg:
        return ZX.mean()+XZ.mean()-1
    else:
        return ZX.flatten(), XZ.flatten()