from qutip.qip.circuit import QubitCircuit, Gate
from qutip.qip.operations import gate_sequence_product
from itertools import product, takewhile
from functools import reduce
from qutip.operators import *
from qutip.states import *
from qutip.expect import expect
from qutip.tensor import *
import numpy as np

"""
This is a module based on QuTip to implement the LCC circuit composition .

"""


'''
specific function only apply on 12bit-LCC experimental data
'''

def exp4bit(a,n,b,mask,circ):
    '''
    a: 1-6 determine the input state |0x0| |1x1| |+x+| |-x-| |+ix+i| |-ix-i| 
    n: 1-6, map to 0,1,2 determine the output base |0x0| |1x1| X, X, Y, Y 
    b: 0/1 xz or zx
    '''
    c=(n//2+2)%3
    eigs=apply_mask(get_eigs(4,n),mask)
    return get_expectation(circ[6*a+3*b+c],eigs,reverse_order=True)

def exp3bit(a,b,mask,circ):
    '''
    a: 1-6, determine the input state 
    b: 0/1 xz or zx
    '''
    eigs=apply_mask(get_eigs(3,-1),mask)
    return get_expectation(circ[2*a+b],eigs,reverse_order=True)

def E_LCC_comb_3x(circuit1,circuit2,base="XZ",mask='111111111111'):
    '''
    cir1: 4bit 
    cir2: 3bit
    base: XZ bases
    '''
    n_qubit=len(mask)
    assert n_qubit%3==0
    N=n_qubit//3
    
    base_index={"XZ":[0,1]*N,"ZX":[1,0]*N}[base]# b is not global
    Coef=(1,1,1/2,-1/2,1/2,-1/2) # order Z, X, Y

    def indexed_expectation(index_tuple):
        E=reduce(lambda a,b: a*b, [Coef[i] for i in index_tuple])
        E*=exp4bit(2,index_tuple[0],base_index[0],mask[:3]+"1",circuit1)
        for i in range(1,N-1):
            E*=exp4bit(index_tuple[i-1],index_tuple[i],base_index[i],mask[3*i:3*(i+1)]+"1",circuit1)
        E*=exp3bit(index_tuple[N-2],base_index[N-1],mask[3*(N-1):3*N],circuit2) #last one is N-2
        return E

    
    E=0 
    for index_tuple in product(range(6),repeat=N-1):
        E+=indexed_expectation(index_tuple)
    return E


def LCC(n):
    '''
    Helper function to implement a n qubit linear coupled cluster circuit
    '''
    qc=QubitCircuit(N=n)
    for i in range(n):
        qc.add_gate("SNOT",targets=i)
    for j in range(n-1):
        qc.add_gate("CSIGN",controls=j,targets=j+1)
    return qc


def circuit_compose_expectation(circuits, bases):
    '''
    Given a list of circuits and corresponding measurement bases, 
    perform measurment on each circuits using provided bases and compose these result to a larger
    circuit.
    
    Args：
        circuits(list): list of input circuits, U1, U2, U3.. Qobj
        bases(list): list of measurement bases, eg. [sigmax(),sigmaz()]*12
    Return:
        E(real): expectation value of measurement of the composed circuit
    '''
    if type(circuits[0])==QubitCircuit:
        circuits=[gate_sequence_product(qc.propagators()) for qc in circuits]
    N=len(circuits)
    n_qb=[len(c.dims[1]) for c in circuits]
    accu_n_qb=[sum(n_qb[:i])-i for i in range(N)]
    assert sum(n_qb)-N+1==len(bases)

    Coef=(1,1,1/2,-1/2,1/2,-1/2)
    INPUT=[basis(2,0),basis(2,1),
           (basis(2,0)+basis(2,1)).unit(),
           (basis(2,0)-basis(2,1)).unit(),
           (basis(2,0)+1.0j*basis(2,1)).unit(),
          (basis(2,0)-1.0j*basis(2,1)).unit()] 
    OUTPUT=[basis(2,0)*basis(2,0).dag(),basis(2,1)*basis(2,1).dag(),sigmax(),sigmax(),sigmay(),sigmay()] 
    
    def indexed_expectation(index_tuple):
        IN=[tensor([basis(2,0)]*n_qb[0])]+\
        [tensor([INPUT[index_tuple[i-1]]]+[basis(2,0)]*(n_qb[i]-1)) for i in range(1,N)]
        
        OUT=[bases[accu_n_qb[i]:accu_n_qb[i+1]]+[OUTPUT[index_tuple[i]]] for i in range(N-1)]+\
        [bases[accu_n_qb[-1]:]]
        
        coef=reduce(lambda a,b: a*b, map(lambda index: Coef[index],index_tuple))
        exp=reduce(lambda a,b: a*b,[expect(tensor(OUT[i]),circuits[i]*IN[i]) for i in range(N)])
        return coef*exp
    
    E=0
    for index_tuple in product(range(6),repeat=N-1):
        E+=indexed_expectation(index_tuple)
    return E

def to_alphabet_base(base_str,mask_arr=[],return_qutip=True):
    '''
    Given a base string like "XZXZXZXZ" and mask string like "01011010", replace the measurement base from Pauli to 
    identity using the mask string
    
    Args:
        base_str(string):
        mask_arr(string):
    Return:
        res(string || list of Gates)
    '''
    if mask_arr:
        assert len(mask_arr)==len(base_str)
        res=''
        for i in range(len(mask_arr)):
            if mask_arr[i]=='0':
                res+='I'
            else:
                res+=base_str[i]
    else:
        res=base_str
        
    mapping={"X":sigmax(),"Y":sigmay(),"Z":sigmaz(),"I":identity(2)}
    if return_qutip:
        return list(map(lambda x: mapping[x], res))
    else:
        return res

def get_expectation(state_vector, eigs=[], reverse_order=False):
    '''
    Args:
        state_vector(array::float/complex): 2^n elements, quantum state vector in probability or amplitude 
        eigs(2darray::floats): (n,2) array to determine the corresponding eigen values for each qubit,
            default values set to be ((1,-1),(1,-1),...)
    Return:
        expectation(float): expectation of the state 
    '''
    if state_vector.dtype == 'complex':
        state_vector = np.vectorize(
            lambda x: np.linalg.norm(x)**2)(state_vector)

    n = int(np.log2(len(state_vector)))
    expectation = 0

    if eigs:
        if type(eigs) != np.ndarray:
            eigs = np.array(eigs)
        assert eigs.shape == (n, 2)
    else:
        eigs = np.array([(1, -1)]*n)

    for i in range(len(state_vector)):
        eig_val = 1
        # get the bitstring of corresponding vector element
        if reverse_order:
            bit_str = (bin(i)[2:]).zfill(n)[::-1]
        else:
            bit_str = (bin(i)[2:]).zfill(n)
        for j in range(n):
            # get eigen value of the bitstring
            eig_val *= eigs[j, int(bit_str[j])]
        # add probability times eigen value
        expectation += state_vector[i]*eig_val
    return expectation

def apply_mask(eig_arr, mask_arr):
    '''
    apply mask on the array of eigen values 
    
    >>> apply_mask([(1,-1),(1,0),(0,1),(1,-1)],"0101")
    [(1, 1), (1, 0), (1, 1), (1, -1)]
    '''
    assert len(eig_arr) == len(mask_arr)
    for i in range(len(mask_arr)):
        if mask_arr[i] == '0':
            eig_arr[i] = (1, 1)
    return eig_arr


def get_eigs(n, i):
    '''
    >>> get_eigs(4, 1)
    [(1, -1), (1, -1), (1, -1), (0, 1)]
    '''
    if i == 0:
        eigs = [(1, -1)]*(n-1)+[(1, 0)]
    elif i == 1:
        eigs = [(1, -1)]*(n-1)+[(0, 1)]
    else:
        eigs = [(1, -1)]*n
    return eigs

def reverse_binary_order(state):
    '''
    input a quantum state output the state with all it's elements reversed 
    '''
    m=len(state)
    n=int(np.log2(m))
    new=np.zeros(m)
    for i in range(m):
        j=int(bin(i)[2:].zfill(n)[::-1],2)
        new[i]=state[j]
    return new

def data_parser(data_path,return_predict=False):
    '''
    ONLY APPLY ON THE GIVEN TEXT DATE FORM IN THIS EXPERIMENT 
    
    Arg:
        data_path(str): directory of data file
    Return:
        data(array): m x n array with each row to be a state vector of measurement result
    '''
    data=[]
    predicted=[]
    
    with open(data_path, "r") as f:
        raw=f.readlines()
    for i in range(len(raw)):
        if "===" in raw[i]:
            data.append(list(map(float,raw[i+1].strip().split(" "))))
            predicted.append(list(map(float,raw[i+2].strip().split(" "))))
            
    data=np.array(data)
    predicted=np.array(predicted)
    
    if return_predict:
        return data, predicted
    return data

def load_mask_file(file_name,N):
    """
    load mask array from a file generate by the mathematica script
    """
    mask_set=[]
    with open(file_name,"r") as f:
        for line in f:
            mask=["1"]*N
#             print(str(line).strip()+",")
            try:
                m=eval(str(line).strip()+",")
                for loc in m:
                    mask[loc-1]='0'
            except:
                pass
            mask_set.append(reduce(lambda x,y: x+y, mask))
    return mask_set

if __name__ == "__main__":
    import doctest
    doctest.testmod()