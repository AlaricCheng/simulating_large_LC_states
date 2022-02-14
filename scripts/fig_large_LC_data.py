'''
Generate Fig. 4 and Fig. S1 (b)
'''
import sys
sys.path.append('./scripts/lib')

from local_lib import *
import matplotlib.pyplot as plt
from tensor_processing import *
from time import perf_counter
# from copy import deepcopy


##### Load data #####
circ1_arr, circ1_the_arr=data_parser("./cache/big_circuit3_check_temp_likelyrho_1.txt",return_predict=True)
circ2_arr, circ2_the_arr=data_parser("./cache/big_circuit3_check_temp_likelyrho_2.txt",return_predict=True)
n_rep = int(circ1_arr.shape[0]/36) # num of repetitions

n_scaling=[3*k for k in range(2,12)]
XZ_masks_cache={f"{n}": load_mask_file("./cache/mask_cache/XZ_{0}_mask.csv".format(n),N=n) for n in n_scaling}
ZX_masks_cache={f"{n}": load_mask_file("./cache/mask_cache/ZX_{0}_mask.csv".format(n),N=n) for n in n_scaling}


##### constructed from small circuits #####
mask4=bitstr(4)
mask3=bitstr(3)
data_XZ_ml = []
data_ZX_ml = []
for j in range(n_rep): 
    circ1,circ2=circ1_arr[j*36:(j+1)*36], circ2_arr[j*12:(j+1)*12]
    T4={m: (circ1_tensor(0,m, circ1),circ1_tensor(1,m, circ1)) for m in mask4}
    V4={m: (circ1_vector(0,m, circ1),circ1_vector(1,m, circ1)) for m in mask4}
    V3={m: (circ2_vector(0,m, circ2),circ2_vector(1,m, circ2)) for m in mask3}
    ZX_12, XZ_12 = compose_fidelity(12, V4, T4, V3, return_avg=False)
    data_XZ_ml.append(XZ_12)
    data_ZX_ml.append(ZX_12)

print(f"Data shape: {np.array(data_XZ_ml).shape}")
print(f"Fidelity from circuit-cutting scheme: {np.mean(data_XZ_ml) + np.mean(data_ZX_ml) - 1}")


##### 12-qubit data #####
from projectq_lib import *
eng=MainEngine()
LCC12_ZX=circuit_tomography(eng,[(H,)]*12,[BaseZ,BaseX]*6)
LCC12_XZ=circuit_tomography(eng,[(H,)]*12,[BaseX,BaseZ]*6)
dist_12 = np.genfromtxt("./cache/cluster12_likelyrho.csv", delimiter=",")[:,:-1]

data_12_ZX_ml = []
for mask in ZX_masks_cache["12"]:
    data_12_ZX_ml.append(get_expectation(dist_12[0], apply_mask(get_eigs(12,-1), mask)))

data_12_XZ_ml = []
for mask in XZ_masks_cache["12"]:
    data_12_XZ_ml.append(get_expectation(dist_12[1], apply_mask(get_eigs(12,-1), mask)))

print(f"Fidelity from the 12-qubit experiment: {np.mean(data_12_XZ_ml) + np.mean(data_12_ZX_ml) - 1}")

###############

def circuit_cutting_vs_12_qubit_XZ():
    '''Plot Fig. 4 (a)'''
    hfont = {"fontsize": 16}
    fig,ax=plt.subplots(figsize = (10, 4), dpi=120)

    bar1 = ax.bar(np.arange(len(XZ_masks_cache["12"])) - 0.2, np.mean(data_XZ_ml, axis = 0), yerr=np.std(data_XZ_ml, axis = 0), capsize = 2, width = 0.4, label = "Circuit-cutting")
    bar2 = ax.bar(np.arange(len(XZ_masks_cache["12"])) + 0.2, data_12_XZ_ml, width = 0.4, label = "12-qubit", color = "tab:orange")
    ax.set_ylim([0, 1.1])
    ax.set_ylabel("Expectation", fontsize = hfont["fontsize"])
    ax.set_title("Measurement Basis (XZ)", fontsize = hfont["fontsize"])
    ax.tick_params(axis='both', labelsize=hfont["fontsize"])

    fig.legend([bar1, bar2], ["Circuit cutting", "12-qubit"], fontsize = hfont["fontsize"], bbox_to_anchor = [1, 1.10], frameon = False)

    fig.tight_layout()
    fig.savefig("./output/circuit_cutting_vs_12_qubit_XZ.svg",bbox_inches = 'tight')


def circuit_cutting_vs_12_qubit_ZX():
    '''Plot Fig. S1. (b)'''
    hfont = {"fontsize": 16}
    fig,ax=plt.subplots(figsize = (10, 4), dpi=120)

    bar1 = ax.bar(np.arange(len(ZX_masks_cache["12"])) - 0.2, np.mean(data_ZX_ml, axis = 0), yerr=np.std(data_ZX_ml, axis = 0), capsize = 2, width = 0.4, label = "Circuit-cutting")
    bar2 = ax.bar(np.arange(len(ZX_masks_cache["12"])) + 0.2, data_12_ZX_ml, width = 0.4, label = "12-qubit", color = "tab:orange")
    ax.set_ylim([0, 1.1])
    ax.set_ylabel("Expectation", fontsize = hfont["fontsize"])
    ax.set_title("Measurement Basis (ZX)", fontsize = hfont["fontsize"])
    ax.tick_params(axis='both', labelsize=hfont["fontsize"])

    fig.legend([bar1, bar2], ["Circuit cutting", "12-qubit"], fontsize = hfont["fontsize"], bbox_to_anchor = [1, 1.10], frameon = False)


    fig.tight_layout()
    fig.savefig("./output/circuit_cutting_vs_12_qubit_ZX.svg",bbox_inches = 'tight')


def fidelity_decay():
    '''Plot Fig. 4 (b)'''
    ##### Process data with the circuit-cutting scheme #####
    data_ml=[]
    time_ml = []
    for j in range(n_rep): 
        circ1,circ2=circ1_arr[j*36:(j+1)*36],circ2_arr[j*12:(j+1)*12]
        T4={m: (circ1_tensor(0,m, circ1),circ1_tensor(1,m, circ1)) for m in mask4}
        V4={m: (circ1_vector(0,m, circ1),circ1_vector(1,m, circ1)) for m in mask4}
        V3={m: (circ2_vector(0,m, circ2),circ2_vector(1,m, circ2)) for m in mask3}
        tmp_data = []
        tmp_time = []
        for n in n_scaling:
            tick = perf_counter()
            tmp_data.append(compose_fidelity(n, V4, T4, V3))
            tmp_time.append(perf_counter() - tick)
        data_ml.append(tmp_data)
        time_ml.append(tmp_time)
    data_ml = np.array(data_ml)
    time_ml = np.array(time_ml)

    ##### Plot figure #####
    fig, ax=plt.subplots(figsize = (10,4), dpi=120)
    fontsize = 16

    line1 = ax.errorbar(np.array(n_scaling), data_ml.mean(axis=0), yerr=data_ml.std(axis=0), capsize = 3, color = "royalblue", label = "Fidelity bound", linestyle = "--")
    ax.set_xlabel("Effective number of qubits", fontsize = fontsize)
    ax.set_ylabel("Fidelity bound", fontsize = fontsize)

    ax.set_xticks(n_scaling, fontsize = fontsize)
    ax.tick_params(axis = "both", labelsize = fontsize)

    ax2 = ax.twinx()
    line2 = ax2.errorbar(np.array(n_scaling), time_ml.mean(axis = 0), yerr = time_ml.std(axis=0), capsize = 3, color = "lightcoral", label = "Processing time", linestyle = "--")
    ax2.set_ylabel("Processing time (s)", fontsize = fontsize)
    ax2.set_yscale("log")
    ax2.tick_params(axis = "y", labelsize = fontsize)
    # ax2.legend()

    ax.legend(handles = [line1, line2], loc = "lower center", fontsize = fontsize)
    fig.tight_layout()
    fig.savefig("./output/fidelity_decay.svg", bbox_inches = 'tight')



if __name__ == "__main__":
    circuit_cutting_vs_12_qubit_XZ()
    circuit_cutting_vs_12_qubit_ZX()
    fidelity_decay() 