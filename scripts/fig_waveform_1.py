import os
from scipy.io import loadmat as load
import matplotlib.pyplot as plt
import numpy as np
basedir='../raw'
a=os.listdir(basedir)
b=load(basedir+'/waveform_filenames.mat')

es = []
for fi in range(48):
    c=b['big_circuit3_filenames'][0][0+fi][0].split('\\')[-1]
    d=load(basedir+'/waveform/'+c)
    e=d['sequenceSamples']
    es.append(e)
    # >>> e.shape
    # (12, 20626)

def hide0(arr,xyz):
    if xyz%3==2:
        return arr
    b=arr.copy()
    b[abs(b)<20]=np.nan
    return b

# default color cycle - 10 
# color table https://www.jianshu.com/p/d0b789ba5be9
# ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
color12 = ['pink', 'brown', 'goldenrod', 'blue', 'orange', 'gray', 'red', 'violet', 'olive', 'cyan', 'purple', 'green']
color9 = color12[3:]
DIS=60000
HDIS=DIS//2

# plt.figure()
# for fi in range(36):
#     e=es[fi]
#     # plt.figure()
#     for ii in range(e.shape[0]):
#         plt.plot(e[ii,:750]-(ii//3)*DIS)

# plt.figure()
# for fi in range(36,48):
#     e=es[fi]
#     # plt.figure()
#     for ii in range(e.shape[0]):
#         plt.plot(e[ii,:480]-(ii//3)*DIS)

def plotone(elist,maxx):
    plt.figure()
    for fi in elist:
        e=es[fi]
        # plt.figure()
        for ii in range(e.shape[0]):
            plt.plot(e[ii,:maxx]-(ii//3)*DIS)

case=6

if case==0:
    plotone(range(36),750)
    plotone(range(36,48),480)

if case==1:
    plotone([6],750)
    plotone([7],750)
    plotone([8],750)

if case==2:
    plotone([0],750)
    plotone([6],750)
    plotone([12],750)
    plotone([18],750)
    plotone([24],750)
    plotone([32],750)

if case==3:
    plotone([6],750)
    plotone([9],750)

if case==4:
    plotone([36],480)
    plotone([38],480)
    plotone([40],480)
    plotone([42],480)
    plotone([44],480)
    plotone([46],480)

if case==5:
    plotone([38],480)
    plotone([39],480)

if case==6:
    signlinewidth=[1.2,1.2,2]*4
    framelinewidth=1.5
    e=es[6].copy()

    e[0:3,346:484]=0
    e[0:9,612:750]=0
    e[9:12,612:750]=0

    plt.figure(figsize=[6.4*2.8/2,4.8*2.2/2])
    ax=plt.gca()

    # main
    for ii in range(e.shape[0]):
        plt.plot(hide0(e[ii,:750],ii)-(ii//3)*DIS,color=color12[ii],linewidth=signlinewidth[ii])

    rect=plt.Rectangle((346,-HDIS), 484-346, DIS, fill=False, edgecolor = 'gray',linewidth=framelinewidth)
    ax.add_patch(rect)
    rect=plt.Rectangle((612,-HDIS*7), 750-612, DIS-1000, fill=False, edgecolor = 'orange',linewidth=framelinewidth)
    ax.add_patch(rect)
    rect=plt.Rectangle((612,-HDIS*5), 750-612, DIS*3, fill=False, edgecolor = 'red',linewidth=framelinewidth)
    ax.add_patch(rect)

    # inital state
    IGAP1=14000
    IGAP2=50
    for jj,ei in enumerate([0,6,12,18,24,32]):
        e=es[ei][0:3,346:484].copy()
        if jj>1:
            e[:,422-484:]=0
        if jj==0:
            e[:,:]=0
        x=np.arange(346,484)-346+jj*(484-346+IGAP2)
        for ii in range(e.shape[0]):
            plt.plot(x,hide0(e[ii],ii)-(ii//3)*DIS+DIS+IGAP1,color=color12[ii],linewidth=signlinewidth[ii])
        rect=plt.Rectangle((x[0],HDIS+IGAP1), 484-346, DIS, fill=False, edgecolor = 'gray',linewidth=framelinewidth)
        ax.add_patch(rect)

    # oi 测量
    for jj,ei in enumerate([6,7,8]):
        e=es[ei][9:12,612-60:750-60].copy()
        if jj==0:
            e=es[ei][9:12,612:750].copy()
        x=np.arange(612,750)+(jj+1)*(750-612+IGAP2)
        for ii in range(e.shape[0]):
            plt.plot(x,hide0(e[ii],ii)-(ii//3)*DIS-3*DIS,color=color12[9:][ii],linewidth=signlinewidth[ii])
        rect=plt.Rectangle((x[0],-HDIS*7), 750-612, DIS-1000, fill=False, edgecolor = 'orange',linewidth=framelinewidth)
        ax.add_patch(rect)
    
    # xzx zxz
    for jj,ei in enumerate([6,9]):
        e=es[ei][:9,612:750].copy()
        x=np.arange(612,750)+(jj+1)*(750-612+IGAP2)
        for ii in range(e.shape[0]):
            plt.plot(x,hide0(e[ii],ii)-(ii//3)*DIS,color=color12[ii],linewidth=signlinewidth[ii])
        rect=plt.Rectangle((x[0],-HDIS*5), 750-612, DIS*3, fill=False, edgecolor = 'red',linewidth=framelinewidth)
        ax.add_patch(rect)

if case==7:
    signlinewidth=[1.2,1.2,2]*4
    framelinewidth=1.5
    e=es[38].copy()

    e[0:3,:138]=0
    e[0:9,412:480]=0

    plt.figure(figsize=[6.4*2.8/2,4.8*2.2/2])
    ax=plt.gca()

    # main
    for ii in range(e.shape[0]):
        plt.plot(hide0(e[ii,:480],ii)-(ii//3)*DIS,color=color9[ii],linewidth=signlinewidth[ii])

    rect=plt.Rectangle((0,-HDIS), 138, DIS, fill=False, edgecolor = 'gray',linewidth=framelinewidth)
    ax.add_patch(rect)
    rect=plt.Rectangle((412,-HDIS*5), 480-412, DIS*3, fill=False, edgecolor = 'red',linewidth=framelinewidth)
    ax.add_patch(rect)

    # inital state
    IGAP1=14000
    IGAP2=20
    for jj,ei in enumerate([36,38,40,42,44,46]):
        e=es[ei][0:3,:138].copy()
        if jj>1:
            e[:,78:]=0
        if jj==0:
            e[:,:]=0
        x=np.arange(0,138)-0+jj*(138-0+IGAP2)
        for ii in range(e.shape[0]):
            plt.plot(x,hide0(e[ii],ii)-(ii//3)*DIS+DIS+IGAP1,color=color9[ii],linewidth=signlinewidth[ii])
        rect=plt.Rectangle((x[0],HDIS+IGAP1), 138, DIS, fill=False, edgecolor = 'gray',linewidth=framelinewidth)
        ax.add_patch(rect)

    # xzx zxz
    for jj,ei in enumerate([38,39]):
        e=es[ei][:9,412:480].copy()
        x=np.arange(412,480)+(jj+1)*(480-412+IGAP2)+30
        for ii in range(e.shape[0]):
            plt.plot(x,hide0(e[ii],ii)-(ii//3)*DIS,color=color9[ii],linewidth=signlinewidth[ii])
        rect=plt.Rectangle((x[0],-HDIS*5), 480-412, DIS*3, fill=False, edgecolor = 'red',linewidth=framelinewidth)
        ax.add_patch(rect)

# 单比特(60)

# 0~35
# 创建12*750的0
# 开头以及cz门的部分用[6]的 0-346 484-612
# 6个初态,分别用 [0]空 [6]346-484 [12]346-422 [18]346-422 [24]346-422 [32]346-422
# oi测量, [6]612-690 [7]612-690 ~x+60 [8]612-690 ~x+60 
# xzx zxz [6]612-690 [9]612-690

# 36~47
# 创建9*500的0
# cz用[38]的 138-412
# 6个初态,分别用 [36]空 [38]0-138 [40]0-78 [42]0-78 [44]0-78 [46]0-78
# xzx zxz [38]412-480 [39]412-480

# plt.show()
plt.savefig('../output/plot_b_raw.pdf',dpi=600)