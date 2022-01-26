import sys
sys.path.append('lib')

from local_lib import *
import matplotlib.pyplot as plt
from tensor_processing import *
from copy import deepcopy

# https://matplotlib.org/stable/tutorials/colors/colormaps.html
fontsize = 16

data_list_XZ = [bit4_XZ_the, bit4_XZ, bit3_XZ_the, bit3_XZ]
min_val = np.min([data_list_XZ[i].min() for i in range(4)])
max_val = np.max([data_list_XZ[i].max() for i in range(4)])
fig, axes = plt.subplots(nrows=2, ncols=2, figsize = (8, 6), dpi = 120)
for i in range(4):
    im = axes.flat[i].imshow(data_list_XZ[i] , vmin=min_val, vmax=max_val, cmap = "coolwarm")
    axes.flat[i].set_xticks(range(4))
    axes.flat[i].tick_params(axis = "both", labelbottom = False, labelleft = False, bottom = False, left = False)
    
axes.flat[0].set_xlabel("4-qubit theory", fontsize = fontsize)
axes.flat[1].set_xlabel("4-qubit experiment", fontsize = fontsize)
axes.flat[2].set_xlabel("3-qubit theory", fontsize = fontsize)
axes.flat[3].set_xlabel("3-qubit experiment", fontsize = fontsize)

cbar = fig.colorbar(im, ax=axes.ravel().tolist())
cbar.ax.tick_params(labelsize = fontsize)

# fig.tight_layout()
fig.savefig("./fig/XZ_original_data.svg",bbox_inches = 'tight')

plt.show()