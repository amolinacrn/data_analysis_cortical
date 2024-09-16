# -*- coding: utf-8 -*-
from cProfile import label
from lib2to3.pytree import HUGE
from operator import imod
import re
from statistics import stdev
from turtle import fd, shape
import matplotlib

import numpy as np
import os
import numpy 
import math
import seaborn as sns
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
import networkx as nx
import matplotlib as mpl
#from clasificar import componente_gigante
matplotlib.rcParams['text.usetex'] = True    
def mean_degree_network(mxc):
    # #En esta seccio se calcula el grado medio de la red

    xGz = nx.from_numpy_array(mxc)
    strength = dict(xGz.degree(weight=1))
    normstrengthlist=list(strength.values())
    mean_degree=np.mean(normstrengthlist)
    std_degree=np.std(normstrengthlist)/np.sqrt(len(normstrengthlist))
    return normstrengthlist, mean_degree, std_degree

def mibarcax(cv,pargraf,cmap,label_format,ticks_loc,xcvd,fig):
    norm= matplotlib.colors.Normalize(vmin=min(cv), vmax=max(cv)) 
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm) 
    ax2_divider = make_axes_locatable(pargraf)
    cax0 = ax2_divider.append_axes("bottom", size="4%", pad="35%")
    cbar0=fig.colorbar(sm,extend='both' ,cax=cax0,orientation='horizontal')
    cbar0.ax.set_xlabel(xcvd,fontsize=15)   
    cbar0.ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar0.ax.set_xticklabels([label_format.format(x) for x in ticks_loc],fontsize=13)
    cax0.xaxis.set_ticks_position("bottom")

def fun_histogram(datos_histo,interval):
    fx, xb=np.histogram(datos_histo,bins=interval)
    sx=xb[:-1]; sy= fx/sum(fx)
    return sx, sy
   
def plot_metricas(maxd_lp,maxd_bins,pargraf,xm,interv_clust,
                    newcmp,i):
    
    Gtx = nx.from_numpy_array(xm)

    Gtc = {(e1, e2): abs(weight) for e1, e2, weight in Gtx.edges(data='weight')}

    #Luego agreguelos como atributos a los bordes del grafico.
    nx.set_edge_attributes(Gtx, Gtc, 'distance')
    vector_k_lover,nkover,err_kover= mean_degree_network(xm)
    lp=dict(nx.all_pairs_dijkstra_path_length(Gtx,weight='distance'))
    log_path=[]; matz_lp=[]

    for os in range(len(lp)):
        log_path += list(lp[os].values())
    log_path=np.array(log_path)[np.array(log_path)!=0]
    std_pl=np.std(log_path)/np.sqrt(len(log_path))
    path_lengh=np.mean(log_path)  
        
    clustering = list(nx.clustering(Gtx, weight='weight').values())
    xlp=np.array(log_path)

    '''
    
    #if(ft==0):
    sx,sy=fun_histogram(xlp,interv_clust)
    maxY0=len(interv_clust)-2
    # if(ft==1):
    #     sx,sy=fun_histogram(xlp,interv_lpath)
    #     maxY0=len(interv_lpath)-2

    #----------------------------------------------
    
    maxd_lp.append(sy) 
    maxd_bins.append(sx)   
    pargraf[1].plot(sx,sy,linewidth=0.3, color=newcmp(i))
    '''
    return path_lengh

def par_graf(k,maxd_bins,maxd_lp,pargraf,
            xnewcmp,altura_abajo,maxY0,zeje_name,eje_name,cv,
            lformat,yfrmt,lfoy,xfrmt,xeje_name,xmicol0,xmicol1,fig):
    
    mz_maxd_bins=np.array(maxd_bins).flatten()
   
    matriz_desity_lp=np.array(maxd_lp)
    
    xvc=matriz_desity_lp.flatten()

    
    # normx= matplotlib.colors.Normalize(vmin=min(xvc), vmax=max(xvc)) 
    # smx = matplotlib.cm.ScalarMappable(cmap=xcmap, norm=normx)

    pargraf[0].imshow(np.transpose(matriz_desity_lp),cmap=xnewcmp,aspect='auto')

    #-----------------------------------------------
    minX1=0;  minY1=0; minX0=0; minY0=0
    maxX1=max(mz_maxd_bins); maxY1=altura_abajo; maxX0=len(cv); 
    zmaxY1=max(xvc)

    pargraf[0].invert_yaxis()

    pargraf[1].set_ylim(minY1,maxY1)
    pargraf[1].set_xlim(minX1,maxX1)
    
    pargraf[0].set_ylim(minY0,maxY0)
    pargraf[0].set_xlim(minX0,maxX0)

    pargraf[0].xaxis.set_tick_params(labelsize=13)
    pargraf[0].set_xlabel(zeje_name,fontsize=15)
    pargraf[1].set_xlabel(eje_name,fontsize=20)
    
    pargraf[1].xaxis.set_tick_params(labelsize=13)
    pargraf[1].set_xlabel(eje_name,fontsize=20)
    pargraf[1].xaxis.set_tick_params(labelsize=13)

   
    num_ejex =np.linspace(minX0,maxX0,6)
    list_numers_x=np.linspace(min(cv),max(cv),6)

    pargraf[0].xaxis.set_major_locator(mticker.FixedLocator(num_ejex))
    pargraf[0].set_xticklabels([lformat.format(x) for x in list_numers_x],fontsize=13)
    
    njey =np.linspace(minY0,maxY0,8)
    list_numers_y=np.linspace(minX1,maxX1,8)
    pargraf[0].yaxis.set_major_locator(mticker.FixedLocator(njey))
    pargraf[0].set_yticklabels([yfrmt.format(x) for x in list_numers_y],fontsize=13)

    num_ejey =np.linspace(minY1,maxY1,6)
    pargraf[1].yaxis.set_major_locator(mticker.FixedLocator(num_ejey))
    pargraf[1].set_yticklabels([lfoy.format(x) for x in num_ejey],fontsize=13)

    njex =np.linspace(minX1,maxX1,6)
    pargraf[1].xaxis.set_major_locator(mticker.FixedLocator(njex))
    pargraf[1].set_xticklabels([xfrmt.format(x) for x in njex],fontsize=13)



    if(k==0):
        pargraf[0].yaxis.set_tick_params(labelsize=13)
        pargraf[0].set_ylabel(eje_name,fontsize=20)
        pargraf[1].yaxis.set_tick_params(labelsize=13)
        pargraf[1].set_ylabel(xeje_name,fontsize=20)

    else:
        pargraf[0].tick_params(axis="y", labelleft=False)
        pargraf[1].tick_params(axis="y", labelleft=False)

    # mibarcax(xvc,pargraf[0],xmicol1,r'${:.2f}$',np.linspace(minY1,zmaxY1,5),r"$\rho(k_i)$",fig)
    # mibarcax(cv,pargraf[1],xmicol0,r'${:.1f}$',np.linspace(min(cv),max(cv),6),r"$CV$",fig)




