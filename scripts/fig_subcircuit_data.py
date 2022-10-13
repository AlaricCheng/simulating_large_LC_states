'''
    Generate Fig. 3 and Fig. S1 (a)

    In `XZ_original_data` and `ZX_original_data`, the output probability distributions of 3-qubit and 4-qubit linear cluster states are plotted. The following parameters used in these two functions specify the distributions. 
    - a: 1-6 determine the input state |0x0| |1x1| |+x+| |-x-| |+ix+i| |-ix-i| 
    - n: 1-6, map to 0,1,2 determine the output base |0x0| |1x1| X, X, Y, Y 
    - b: 0/1 xz or zx
'''
import sys
sys.path.append('./scripts/lib')
import argparse

from local_lib import *
import matplotlib.pyplot as plt
import seaborn as sns
# from tensor_processing import *
# from copy import deepcopy


def to_observable(bitstring, basis = "XZ"):
    '''
    Convert a bitstring to an observable
    '''
    n = len(bitstring)
    if n % 2 == 0:
        basis = basis * int(n/2)
    else:
        basis = basis * int(n/2) + basis[0]
    tmp = [basis[i] if bitstring[i] == "1" else "I" for i in range(n)]
    return reduce(lambda x, y: x+y, tmp)    



##### Data that will be used globally #####
# load data
circ1_arr, circ1_the_arr=data_parser("./cache/big_circuit3_check_temp_likelyrho_1.txt",return_predict=True) # 25 groups of data. Rows for different circuits and columns for different measurement outcomes
circ2_arr, circ2_the_arr=data_parser("./cache/big_circuit3_check_temp_likelyrho_2.txt",return_predict=True)
n_rep = int(circ1_arr.shape[0]/36) # num of repetitions, which 25


# process data
circ1=[circ1_arr[j*36:(j+1)*36] for j in range(n_rep)] # only take the first and the last experiments
circ2=[circ2_arr[j*12:(j+1)*12] for j in range(n_rep)]
circ1_the=[circ1_the_arr[j*36:(j+1)*36] for j in range(n_rep)]
circ2_the=[circ2_the_arr[j*12:(j+1)*12] for j in range(n_rep)]


# average distribution
circ1_avg=sum(circ1)/n_rep
circ2_avg=sum(circ2)/n_rep
circ1_the_avg=sum(circ1_the)/n_rep
circ2_the_avg=sum(circ2_the)/n_rep





def XZ_original_data():
    '''Plot Fig. 3 (a)''' 
    fontsize = 16 # set fontsize of the figure
    
    ##### select the distributions corresponding to a linear cluster state with XZ measurement #####
    #XZXZ
    (a,n,b)=(2,0,0)
    c=(n//2+2)%3 # this specifies the distribution from sub-circuit 1
    bit4_XZ=circ1_avg[6*a+3*b+c].reshape((4,4))
    bit4_XZ_the=circ1_the_avg[6*a+3*b+c].reshape((4,4))

    #XZX
    (a,b)=(2,0) # this specifies the distribution from sub-circuit 2
    bit3_XZ=circ2_avg[2*a+b].reshape((2,4))
    bit3_XZ_the=circ2_the_avg[2*a+b].reshape((2,4))
    data_list_XZ = [bit4_XZ_the, bit4_XZ, bit3_XZ_the, bit3_XZ]
    min_val = np.min([data_list_XZ[i].min() for i in range(4)])
    max_val = np.max([data_list_XZ[i].max() for i in range(4)])

    ##### Plot figure #####
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (8, 6), dpi = 120)
    for i in range(4):
#         im = axes.flat[i].imshow(data_list_XZ[i] , vmin=min_val, vmax=max_val, cmap = "coolwarm")
        im = sns.heatmap(data_list_XZ[i] , annot=True, fmt='.2f', linewidths=.5, ax=axes.flat[i], vmin=min_val, vmax=max_val, cmap='viridis', cbar=False)
        axes.flat[i].set_xticks(range(4))
        axes.flat[i].tick_params(axis = "both", labelbottom = False, labelleft = False, bottom = False, left = False)
        
    axes.flat[0].set_xlabel("4-qubit theory", fontsize = fontsize)
    axes.flat[1].set_xlabel("4-qubit experiment", fontsize = fontsize)
    axes.flat[2].set_xlabel("3-qubit theory", fontsize = fontsize)
    axes.flat[3].set_xlabel("3-qubit experiment", fontsize = fontsize)

    cbar = fig.colorbar(im.collections[0], ax=axes.ravel().tolist())
    cbar.ax.tick_params(labelsize = fontsize)

    # fig.tight_layout()
    fig.savefig("./output/XZ_original_data.svg",bbox_inches = 'tight')
    print("XZ_original_data.svg")



def ZX_original_data():
    '''Plot Fig. S1. (a)'''
    fontsize = 16 # set fontsize of the figure

    ##### select the distributions corresponding to a linear cluster state with ZX measurement #####
    #ZXZX 
    (a,n,b)=(2,3,1)
    c=(n//2+2)%3
    bit4_ZX=circ1_avg[6*a+3*b+c].reshape((4,4))
    bit4_ZX_the=circ1_the_avg[6*a+3*b+c].reshape((4,4))

    #ZXZ 
    (a,b)=(2,1)
    bit3_ZX=circ2_avg[2*a+b].reshape((2,4))
    bit3_ZX_the=circ2_the_avg[2*a+b].reshape((2,4))

    data_list_ZX = [bit4_ZX_the, bit4_ZX, bit3_ZX_the, bit3_ZX]

    min_val = np.min([data_list_ZX[i].min() for i in range(4)])
    max_val = np.max([data_list_ZX[i].max() for i in range(4)])
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (8, 6), dpi = 120)
    for i in range(4):
#         im = axes.flat[i].imshow(data_list_ZX[i] , vmin=min_val, vmax=max_val, cmap = "coolwarm")
        im = sns.heatmap(data_list_ZX[i] , annot=True, fmt='.2f', linewidths=.5, ax=axes.flat[i], vmin=min_val, vmax=max_val, cmap='viridis', cbar=False)
        axes.flat[i].set_xticks(range(4))
        axes.flat[i].tick_params(axis = "both", labelbottom = False, labelleft = False, bottom = False, left = False)
        
    axes.flat[0].set_xlabel("4-qubit theory", fontsize = fontsize)
    axes.flat[1].set_xlabel("4-qubit experiment", fontsize = fontsize)
    axes.flat[2].set_xlabel("3-qubit theory", fontsize = fontsize)
    axes.flat[3].set_xlabel("3-qubit experiment", fontsize = fontsize)

    cbar = fig.colorbar(im.collections[0], ax=axes.ravel().tolist())
    cbar.ax.tick_params(labelsize = fontsize)

    # fig.tight_layout()
    fig.savefig("./output/ZX_original_data.svg",bbox_inches = 'tight')
    print("ZX_original_data.svg")




def subcircuit_expectations():
    '''Plot Fig. 3 (b)'''
    
    ##### 4-qubit linear cluster state #####
    n_qubit=4
    XZ_mask_arr_4=load_mask_file("./cache/mask_cache/XZ_{0}_mask.csv".format(n_qubit),N=n_qubit)
    ZX_mask_arr_4=load_mask_file("./cache/mask_cache/ZX_{0}_mask.csv".format(n_qubit),N=n_qubit)

    data_4bit = []
    # expectations for XZ basis
    avg_4_XZ=[]
    std_4_XZ=[]
    (a,n,b)=(2,0,0)
    for mask in XZ_mask_arr_4:
        f=list(map(lambda circ: exp4bit(a,n,b,mask,circ, final_Z = True),circ1))
        data_4bit.append(f)
        avg_4_XZ.append(np.mean(f))
        std_4_XZ.append(np.std(f))

    # expectations for ZX basis
    avg_4_ZX=[]
    std_4_ZX=[]
    (a,n,b)=(2,3,1)
    c=(n//2+2)%3
    for mask in ZX_mask_arr_4:
        f=list(map(lambda circ: exp4bit(a,n,b,mask,circ),circ1))
        data_4bit.append(f)
        avg_4_ZX.append(np.mean(f))
        std_4_ZX.append(np.std(f))

    data_4bit = np.array(data_4bit)
        
    # print(f"Fidelity of 4-qubit LC state: {np.mean(avg_4_XZ) + np.mean(avg_4_ZX) - 1}")
    fidelity_4 = [np.mean(data_4bit[:4, i]) + np.mean(data_4bit[4:, i]) - 1 for i in range(n_rep)]
    print(f"Fidelity of 4-qubit LC state: {np.mean(fidelity_4)}, Std: {np.std(fidelity_4)}")

    ##### 3-qubit linear cluster state #####
    n_qubit=3
    XZ_mask_arr_3 = load_mask_file("./cache/mask_cache/XZ_{0}_mask.csv".format(n_qubit),N=n_qubit)
    ZX_mask_arr_3 = load_mask_file("./cache/mask_cache/ZX_{0}_mask.csv".format(n_qubit),N=n_qubit)

    data_3bit = []
    # expectations for XZ basis
    avg_3_XZ=[]
    std_3_XZ=[]
    (a,b)=(2,0)
    for mask in XZ_mask_arr_3:
        f=list(map(lambda circ: exp3bit(a,b,mask,circ),circ2))
        data_3bit.append(f)
        avg_3_XZ.append(np.mean(f))
        std_3_XZ.append(np.std(f))

    # expectations for ZX basis
    avg_3_ZX=[]
    std_3_ZX=[]
    (a,b)=(2,1)
    for mask in ZX_mask_arr_3:
        f=list(map(lambda circ: exp3bit(a,b,mask,circ),circ2))
        data_3bit.append(f)
        avg_3_ZX.append(np.mean(f))
        std_3_ZX.append(np.std(f))
        
    data_3bit = np.array(data_3bit)
        
    # print(f"Fidelity of the 3-qubit LC state: {np.mean(avg_3_XZ) + np.mean(avg_3_ZX) - 1}")
    fidelity_3 = [np.mean(data_3bit[:4, i]) + np.mean(data_3bit[4:, i]) - 1 for i in range(n_rep)]
    print(f"Fidelity of 3-qubit LC state: {np.mean(fidelity_3)}, Std: {np.std(fidelity_3)}")

    ##### Plot figure #####
    hfont = {'fontname':'Consolas', "fontsize": 16}
    fig, axes = plt.subplots(ncols = 2, figsize = (8, 4), dpi=120)

    # 4-qubit
    xticks_label_4 = list(map(lambda mask: to_observable(mask), XZ_mask_arr_4)) + list(map(lambda mask: to_observable(mask, basis = "ZX"), ZX_mask_arr_4[:-1])) # group the labels
    axes[0].bar(range(len(xticks_label_4)), avg_4_XZ + avg_4_ZX[:-1], yerr=std_4_XZ + std_4_ZX[:-1], capsize = 4, width = 0.6)
    axes[0].set_ylim([0,1.1])
    axes[0].set_title("4-qubit linear-cluster state", fontsize = hfont["fontsize"])
    axes[0].set_ylabel("Expectation", fontsize = hfont["fontsize"])
    axes[0].set_xticks(range(len(xticks_label_4)), labels = xticks_label_4, rotation=30, **hfont)
    axes[0].tick_params(axis='y', labelsize=hfont["fontsize"])

    # 3-qubit
    xticks_label_3 = list(map(lambda mask: to_observable(mask), XZ_mask_arr_3)) + [to_observable(ZX_mask_arr_3[0], basis = "ZX")] # group the labels
    axes[1].bar(range(len(xticks_label_3)), avg_3_XZ+[avg_3_ZX[0]], yerr=std_3_XZ+[std_3_ZX[0]], capsize = 4, width = 0.4, color = "tab:red")
    axes[1].set_ylim([0,1.1])
    axes[1].set_title("3-qubit linear-cluster state", fontsize = hfont["fontsize"])
    axes[1].set_xticks(range(len(xticks_label_3)), labels = xticks_label_3, rotation=30, **hfont)
    axes[1].tick_params(axis='y', left = False, labelleft = False)

    fig.tight_layout()
    fig.savefig("./output/subcircuit_expectations.svg")
    print("subcircuit_expectations.svg")




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-p", choices=["3a", "S1a", "3b"], help = "parameters for different figures")

    args = parser.parse_args()
    if args.p == "3a":
        XZ_original_data()
    elif args.p == "S1a":
        ZX_original_data()
    elif args.p == "3b": 
        subcircuit_expectations()
