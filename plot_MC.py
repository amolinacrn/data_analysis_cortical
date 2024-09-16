from stat import SF_APPEND
from statistics import stdev
from turtle import fd
import numpy as np
import matplotlib
import os
import numpy 
import math
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import seaborn as sns
#from include.get_function_implementations import *
import networkx as nx
import networkx.algorithms.community as nx_comm
from community import community_louvain
import matplotlib as mpl
from matplotlib import cbook
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec
from numpy.linalg import norm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from statsmodels.graphics.gofplots import qqplot
matplotlib.rcParams['text.usetex'] = True 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def config_ejes(k,rango_en_y,rango_en_x,psup):
    '''
    rango_en_y: para los tiks en el eje y
    psup: grafique aqui

    '''
    
    minX0=0; minY0=0
    maxX0=len(rango_en_x)-1 
    maxY0=len(rango_en_y)-1
    
    psup.set_ylim(maxY0,minY0)
    psup.set_xlim(minX0,maxX0)

     #EN X FILA 0
    num_ejex =np.linspace(minX0,maxX0,8)
    lformat = r'${:.0f}$'
    psup.xaxis.set_major_locator(mticker.FixedLocator(num_ejex))
    psup.set_xticklabels([lformat.format(x) for x in num_ejex],fontsize=15)
    
    #EN Y FILA 0
    yfrmt = r'${:.0f}$'
    njey =np.linspace(minY0,maxY0,8)
    psup.yaxis.set_major_locator(mticker.FixedLocator(njey))
    psup.set_yticklabels([yfrmt.format(x) for x in njey],fontsize=15)

    psup.xaxis.set_tick_params(labelsize=15)
    psup.set_xlabel(r"neuron $j$",fontsize=20)

    if(k==0):
        psup.yaxis.set_tick_params(labelsize=15)
        psup.yaxis.set_tick_params(labelsize=15)
        psup.set_ylabel(r"neuron $i$",fontsize=20)

    else:
        psup.tick_params(axis="y", labelleft=False)
        psup.tick_params(axis="y", labelleft=False)


def matrices_tot(kk,mxc,psup,pinf,xcmap,trshl):

    psup.imshow(mxc,cmap=xcmap,aspect='auto')
    iG = nx.from_numpy_array(mxc)
    nx.draw_networkx(iG,node_size=15,node_color="#FFA500", with_labels = False,ax=pinf,width=1,label="jas")
    config_ejes(kk,mxc[:,0],mxc[:,0],psup)
    if(kk==0):
        pinf.set_xlabel(r'$H={:.0f}$'.format(trshl),fontsize=20)
    else:
        pinf.set_xlabel(r'$H={:.4f}$'.format(trshl),fontsize=20)
    #ax.set_title(f"p = {p:.3f}")
def def_color(xnum):
    
    cmap = ListedColormap(["black","#FFA500"])
    viridis =cm.get_cmap('coolwarm')
    newcolors = viridis(np.linspace(0, 1, 256))
    pink = np.array([0,0 ,0 ,1 ])
    if(xnum==0):
        xnum=0
        newcolors[:xnum, :] = pink
    else:
        xnum=1
        newcolors[:xnum, :] = pink
    
    newcmp = ListedColormap(newcolors)
    return newcmp

qzm=[]
fig = plt.figure(figsize=(10, 5))
gs = fig.add_gridspec(2,4)

namExp=["Mar07","Mar10","Mar04","Jan29","Jan21","Jan14","Dez20"]
namUmb=["Mar10","Mar10","Mar10","Mar10"]
for i in range(4):
    xsubdir="mcc_Exp_sua_2020"+namUmb[i]+".txt"
    xdir="mxc_tot/"
    fdat=np.loadtxt(xdir+str(xsubdir))#datos empiricos
    qzm+=fdat.flatten().tolist() 
    
# xnorm= matplotlib.colors.Normalize(vmin=0, vmax=max(qzm)) 


# zcax = plt.axes([0.20, 0.5, 0.65, 0.02])
# sm = matplotlib.cm.ScalarMappable(cmap=cmap,norm=xnorm) 
# xcbar=fig.colorbar(sm,extend='both',cax=zcax,ax=zcax,orientation='horizontal',label="")

treshold=[0,np.mean(qzm),np.mean(qzm)+stdev(qzm),
          np.mean(qzm)+2*stdev(qzm),np.mean(qzm),
          np.mean(qzm),np.mean(qzm)]
# Hb=[r"0",r"",r"",r"",r""]
for i in range(4):
    azsup = fig.add_subplot(gs[0,i])
    azinf = fig.add_subplot(gs[1,i])
    xsubdir="mcc_Exp_sua_2020"+namUmb[i]+".txt"
    xdir="mxc_tot/"
    
    fdx=abs(np.loadtxt(xdir+str(xsubdir)))#datos empiricos

    ik,il=np.where(fdx==0)
    fdx[ik,il]=0.026038427953472405
    xG = nx.from_numpy_array(fdx)  
    xG.remove_edges_from(list(nx.selfloop_edges(xG)))
    matrix=nx.to_scipy_sparse_array(xG)
    mx=np.array(matrix.todense())
    k,l=np.where(mx<treshold[i])
    mx[k,l]=0
    matrices_tot(i,mx,azsup,azinf,def_color(i),treshold[i])
plt.show()
