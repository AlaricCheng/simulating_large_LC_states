from projectq import MainEngine
from projectq.ops import *
from projectq.meta import Dagger
import numpy as np

BaseX = [H]
BaseY = [H, S]
BaseZ = []

"""
ProjectQ library build for the experiment data analysis

"""


def io_circuit(eng, input_gate, output_gate):
    '''
    measure circuit1 on given output_gate, 
    
    Args:
        input_gate(iteratable::tuple::Gate): gate to apply to prepare the required state
        output_gate(iteratable::tuple::Gate): gate to apply to the prepared state before measure
    Return:
        prob(list::float): probability of the given circuit
    '''
    assert len(output_gate) == len(input_gate)
    n = len(output_gate)

    state = eng.allocate_qureg(n)
    for i in range(n):
        for gate in input_gate[i]:
            gate | state[i]
    for i in range(n-1):
        CZ | (state[i], state[i+1])
    with Dagger(eng):  # not sure wether this would work, apply dagger on qubits
        for i in range(n):
            #             output_gate.reverse()
            for gate in output_gate[i]:
                gate | state[i]
    return state

# depreciated, equals to eng.backend.cheat
def circuit_tomography(eng, init, base, reverse_order=False):
    '''
    measure circuit1 on given base, 
    
    Args:
        init(iteratable::Gate): gate to apply to prepare the required state
        base(iteratable::Gate): gate to apply to the prepared state before measure
    Return:
        prob(list::float): amplitude of the given circuit, with vector elements sorted by bitstring order 
    '''
    state = io_circuit(eng, init, base)
    n = len(state)
    eng.flush()

    if reverse_order:
        #         result=[eng.backend.get_amplitude(bin(i)[2:].replace('0','a').replace("1","0").replace("a","1").zfill(n),state) for i in range(2**n)]
        result = [eng.backend.get_amplitude(
            bin(i)[2:].zfill(n)[::-1], state) for i in range(2**n)]
    else:
        result = [eng.backend.get_amplitude(
            bin(i)[2:].zfill(n), state) for i in range(2**n)]
#     different vector order compared with the original order in HiQ
#     result=eng.backend.cheat()
    All(Measure) | state
    return np.array(result)
