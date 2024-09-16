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
from scipy.optimize import curve_fit
#from NAT2 import *
import statsmodels.api as sm
from statistics import stdev
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
from bct.algorithms import number_of_components

#adres_dir = os.listdir('ccExp_sua_2020Mar4w20')  
    
def  draw_comunidades(vxy,nam_exp,zname,mostrar=False):
    
    ffile_txt=nam_exp+".txt"
    
    ffdata = open(ffile_txt, 'w')  

    try:      
        ffdata.write("{")
        for wrfil in range(len(vxy)):
            ffdata.write("{")
            for wrcol in range(len(vxy)):
                ffdata.write(str(vxy[wrfil][wrcol]))
                if(wrcol==len(vxy)-1 and wrfil < len(vxy)-1):
                    ffdata.write("},")
                elif(wrfil == len(vxy)-1 and wrcol==len(vxy)-1):
                    ffdata.write("}")
                else:
                    ffdata.write(",")
        ffdata.write("}")
    finally:
        ffdata.close()

    # xintr=np.arange(0,2.5,0.01)
    # xfrec, xbins=np.histogram(vxy.flatten(),bins=xintr)
    # xx=xbins[:-1]; yy= xfrec/sum(xfrec)
    # #sns.kdeplot(data=D_kls.flatten(),multiple="stack")

    #plt.figure(figsize=(12,10))
    G_fb = nx.from_numpy_array(vxy) #ma
    #plt.figure(figsize=(12,10))
    pos_p = nx.spring_layout(G_fb)

    #prx=nx_comm.louvain_communities(G_fb)

    partitions = community_louvain.best_partition(G_fb)
    nodos_por_comundad=[]

    values = list(partitions.values())
    nNCom=list(set(values))
    
    try:  
        xo=0  
        file_txt=zname+".txt"
        fdata = open(file_txt, 'w')
        
        for cluster_id in nNCom:
            cluster = [node for node in G_fb.nodes() if partitions[node] == cluster_id]
            nodos_por_comundad.append(cluster)
            for wrfil in range(len(cluster)):
                fdata.write(str(cluster_id)+"\t"+str(cluster[wrfil])+"\n")#+str(int(sdt[xo]))+"\n")
                xo+=1
             
    finally:
        ffdata.close()
    
    qQ=nx.community.modularity(G_fb,nodos_por_comundad)
    print("Modularidad = ",qQ)
    
    #nuro de comunidades detectadas
    
    # cluster = [node for node in G_fb.nodes() if partitions[node] == cluster_id]
    # cluster = G_fb.subgraph(cluster)
    
    if(mostrar):
        plt.axis("off")

        nx.draw_networkx(G_fb, pos=pos_p, cmap=plt.get_cmap("jet"), 
                        node_color=values, node_size=10, with_labels = False,width=0.05)
        plt.show()
    
def componente_significativa(D_kl,tresh,cgig,mostrar=False):
    '''
    calculo de la componeten gigante en redes
    '''
    xD_kls=np.zeros((len(D_kl[0,:]),len(D_kl[0,:])))
   
    for i in range(len(D_kl[0,:])):
        for j in range(len(D_kl[0,:])):
            if(cgig==0):
                if(D_kl[i,j]<tresh and D_kl[i,j]!=0):
                    xD_kls[i,j]=D_kl[i,j]
                    xD_kls[j,i]=D_kl[i,j]
                else:
                    xD_kls[i,j]=0
                    xD_kls[j,i]=0
                    
            if(cgig==1):
                if(D_kl[i,j]>tresh and D_kl[i,j]!=0):
                    xD_kls[i,j]=1
                    xD_kls[j,i]=1
                else:           
                    xD_kls[i,j]=0
                    xD_kls[j,i]=0
    if(mostrar):            
        G_fb = nx.from_numpy_array(xD_kls)
        nx.draw_networkx(G_fb, node_size=10, with_labels = False,width=0.05)
        plt.show()
    return xD_kls
               
def componente_gigante(D,treshold,tresh=0.1,xnod=0,mostrar=False):
    '''
    calculo de la componeten gigante en redes
    '''
    #addres_dir = os.listdir('Correlation/ResMar07')
    gdata=[]
    pdata=[]
    xGr = nx.from_numpy_array(D)
    xD_kls=np.zeros((len(D),len(D)))
    xD_aux=np.zeros((len(D),len(D)))
    ox=0
    while treshold<tresh:
        pdata.append(treshold)
        for i in range(len(D)):
            for j in range(len(D)):
                if(D[i,j]<treshold and D[i,j]!=0):
                    xD_kls[i,j]=1
                    xD_kls[j,i]=1
                    xD_aux[i,j]= D[i,j]
                    xD_aux[j,i]= D[i,j]
                    
                else:
                    xD_kls[i,j]=0
                    xD_kls[j,i]=0
                    xD_aux[i,j]=0
                    xD_aux[j,i]=0
                    
        G_fb = nx.from_numpy_array(xD_kls)
        cpct = nx.connected_components(G_fb)
        
        gdata.append(max(len(container) for container in cpct) )
        treshold*=1.05
        if(gdata[ox]==len(xGr.nodes)-xnod):
            break
        ox+=1

    umbral=0

    for i in range(len(gdata)):
        if(gdata[i]==len(xGr.nodes)-xnod):
            print(gdata[i]," ",pdata[i])
            umbral=pdata[i]
            break

    if(mostrar):    
        plt.plot(pdata,gdata,"o-")
        plt.xscale("log")
        plt.show()
    # nx.draw_networkx(G_fb, node_size=10, with_labels = False,width=0.05)
    # plt.show()
    
    # umbral=[]
    # for sr in range(len(pdata)):
    #     if(gdata[sr]==len(D)):
    #         umbral.append(pdata[sr])
   
    xtreshl=round(round(np.min(umbral),5)+0.00002,5)
    print(xtreshl)
    #tcm_graficos(xD_kls)

    
    
    # fig, axs = plt.subplots(1,2 , layout='constrained')
    #axs[0].plot(xx,yy,"bo")
    # pc = axs[1].imshow(D_kls, cmap='plasma',
    #                           norm=mpl.colors.LogNorm(vmin=0.01, vmax=100))
    # fig.colorbar(pc, ax=axs[1], extend='both')
    # axs[1].set_title('imshow() with LogNorm()')
    
    vd=np.array(np.triu(D).flatten())
    frec, xbin=np.histogram(vd[vd>0],bins='auto')
    x=xbin[:-1]; y= frec/sum(frec)
    '''
    fig, main_ax = plt.subplots()
    main_ax.plot(x,y,"-")#,xx,yy,"o",markersize = 3,linewidth=1)
    main_ax.stackplot(x,y)
    # this is an inset axes over the main axes
    right_inset_ax = fig.add_axes([0.4, 0.3,0.5,0.5], facecolor='k')
    pc =right_inset_ax.imshow(xD_kls, cmap='plasma')
    fig.colorbar(pc, extend='both')
    plt.show()
    '''
    return xtreshl,xD_aux


def similitud_cv(cvs):
    
    interval=np.arange(0.01,0.3,0.001)
    c=np.zeros((len(cvs),len(cvs)))
    cv=np.zeros((len(cvs),len(cvs)))

    for i in range(len(cvs)-1):
        for j in range(i+1,len(cvs)):
            zs=abs(cvs[i]-cvs[j])
            c[i,j] = zs
            c[j,i] = zs

    # cx=c.copy()
    # mx=np.array(cx.flatten())
    # xmin=min(mx[mx>0]);xdif=(max(mx[mx>0])-min(mx[mx>0]))
    # print(xmin)
    # for i in range(len(cvs)):
    #     for j in range(len(cvs)):
    #         cv[i,j]=(c[i,j]-xmin)/xdif
    # fig, axs = plt.subplots(1,2 , layout='constrained')
    # axs[0].plot(xx,yy,"bo")
    # pc = axs[1].imshow(mxc, cmap='plasma',
    # norm=mpl.colors.LogNorm(vmin=0.01, vmax=100))
    # fig.colorbar(pc, ax=axs[1], extend='both')
    # axs[1].set_title('imshow() with LogNorm()')
            
    # frec, xbin=np.histogram(c.flatten(),bins='auto')
    # x=xbin[:-1]; y= frec/sum(frec)
    # plt.xscale("log")
    # #plt.yscale("log")
    # plt.plot(x,y,linewidth=0.5)

    plt.show()
    return c

    # return lcv
