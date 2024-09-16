from stat import SF_APPEND
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
from scipy.optimize import curve_fit
from NAT2 import *
import statsmodels.api as sm
from statistics import stdev
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
from bct.algorithms import number_of_components
from portrait_Djs.funct_PDjs import arbole_expansion_min
from libreriaMolCr import delete_nod
#adres_dir = os.listdir('ccExp_sua_2020Mar4w20')  
def fun_ajuste(x,mu,sg,c):
    fn=(c*(sg*math.sqrt(2*np.pi)))*np.exp(-((x-mu)**2/(2*sg**2)))
    return fn

def dist_grafo(nsubdir):
    mxD=(np.loadtxt("redesDjs/"+nsubdir+".txt"))#datos empiricos
    mxD=np.sqrt(np.array(mxD))
    mxc=mxD[mxD>0]

    frec, xbin=np.histogram(mxc,bins=100)
    xo=xbin[:-1]; yo= frec/sum(frec)

    viridis = mpl.colormaps['plasma'].resampled(len(xo))
    newcolors = viridis(np.linspace(0,1, len(xo)))
    
    newcmp = ListedColormap(newcolors)
    fig, ax0 = plt.subplots(figsize=(10,7))



    for row in range(len(xo)):
        ax0.bar(xo[row], yo[row], width=0.0085,ec='w',color=viridis.colors[row]) 
    
    label_format = r'${:.2f}$'

    #para los nodos
    label_format_nodos = r'${:.2f}$'
    xlistnod=np.linspace(min(xo),max(xo),5)
    ylistnod=np.linspace(min(yo),max(yo),5)

    ax0.xaxis.set_major_locator(mticker.FixedLocator(xlistnod))
    ax0.set_xticklabels([label_format_nodos.format(x) for x in xlistnod],fontsize=20)

    ax0.yaxis.set_major_locator(mticker.FixedLocator(ylistnod))
    ax0.set_yticklabels([label_format.format(x) for x in ylistnod],fontsize=20)

    ax0.set_xlabel(r"$D_{JS}$",fontsize=20)
    ax0.set_ylabel(r"$\rho(D_{JS})$",fontsize=20)

    xlist=np.linspace(0,260,8)
    ylist=np.linspace(0,260,8)
        
    
    
    right_inset_ax = fig.add_axes([0.4, 0.40, 0.5, 0.5], facecolor="None")
    divider = make_axes_locatable(right_inset_ax)
    ax_cb = divider.append_axes("right", size="5%", pad=0.05)

    pc =right_inset_ax.imshow(mxD, cmap=newcmp)
    cbar=fig.colorbar(pc,extend='both',cax=ax_cb)

    labelft = r'${:.0f}$'
    label_format = r'${:.2f}$'
    ticks_loc =np.linspace(0,max(mxc),6)#cbar.ax.get_yticks().tolist()

    cbar.ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar.ax.set_yticklabels([label_format.format(x) for x in ticks_loc],fontsize=17)
    cbar.ax.set_ylabel(r"$D_{JS}(i,j)$",fontsize=17)
   

    right_inset_ax.xaxis.set_major_locator(mticker.FixedLocator(xlist))
    right_inset_ax.set_xticklabels([labelft.format(x) for x in xlist],fontsize=17)

    right_inset_ax.yaxis.set_major_locator(mticker.FixedLocator(ylist))
    right_inset_ax.set_yticklabels([labelft.format(x) for x in ylist],fontsize=17)

    right_inset_ax.set_xlabel(r"$j$",fontsize=17)
    right_inset_ax.set_ylabel(r"$i$",fontsize=17)
    ax0.set_xlim(0,1)

    plt.show()

def plot_distribucion_continua(addres_dir,file_dir,interval,mdat,idcv):
    xlin=0.7
    #cmap = plt.colormaps["viridis"]#plasma
    viridis = mpl.colormaps['plasma'].resampled(len(idcv))
    newcolors = viridis(np.linspace(0,1, len(idcv)))
    newcmp = ListedColormap(newcolors)
    #fig, ax = plt.subplots()
    #fig.suptitle(" ", size=14)
    mxrango = slice(14,22)

    fig = plt.figure(figsize=(10,7))
    gs = gridspec.GridSpec(1, 3)
    # ax0 = fig.add_subplot(gs[0,0])
    # ax1 = fig.add_subplot(gs[0,1])
    # ax2 = fig.add_subplot(gs[0,2])
    ax3 = fig.add_subplot(gs[0,:])

    #fig, axx = plt.subplots(figsize=(10,5))
    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    thrho=np.mean(mdat[:,4])-2*np.std(mdat[:,4])# minimo de Desdidad
   
    #if(continuo):
    psMx=[]; mean_network=[]; dtao=[];zy=[]
    
    for i in range(len(idcv)):
        
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(int(idcv[i]))+'.txt')))
        xcv=nx.from_numpy_array(xmcc)
        dzx=nx.density(xcv)

        if(dzx > thrho):
            dtao.append(idcv[i])
            mxc=xmcc[xmcc!=0]
            frec, xbin=np.histogram(mxc,bins=interval)
            x=xbin[:-1]; y= frec/sum(frec)
            zy.append(max(y))
            if(max(y)<0.16):
                ax3.plot(x,y,linewidth=xlin, color=viridis.colors[i-1])
            mean_network.append(np.mean(xmcc))
        else:
            psMx.append(idcv[i])

    cof_var=[]

    for o in dtao:
        for j,ss in enumerate(mdat[:,0]):
            if(int(o)==int(ss)):
                cof_var.append(mdat[int(j),1])
    cof_var.sort()
    
    
    xpeso=np.linspace(0,max(zy),6)
    label_format_peso = r'${:.2f}$'
    ax3.yaxis.set_major_locator(mticker.FixedLocator(xpeso))
    ax3.set_yticklabels([label_format_peso.format(x) for x in xpeso],fontsize=15)

    ax3.xaxis.set_tick_params(labelsize=15)
    ax3.yaxis.set_tick_params(labelsize=15)
    ax3.set_xlabel(r"$\mathcal{W}$",fontsize=20)
    ax3.set_ylabel(r"$\rho(\mathcal{W})$",fontsize=20)


    viridis = mpl.colormaps['plasma'].resampled(len(idcv))
    newcolors = viridis(np.linspace(0,1,len(mdat[:,1])))
    newcmp = ListedColormap(newcolors)

    right_inset_ax = fig.add_axes([0.5, 0.50, 0.4, 0.45], facecolor="None")
    divider = make_axes_locatable(right_inset_ax)
    ax_cb = divider.append_axes("right", size="5%", pad=0.05)
    
  
    df = {'x': cof_var, 'y': mean_network}
    sns.scatterplot(data = df,x = "x", y = "y",marker = "o",hue=cof_var,palette=newcmp,legend=False,ax=right_inset_ax)
    
    cb=plt.cm.ScalarMappable(norm=Normalize(min(mdat[:,1]),max(mdat[:,1]) ), cmap=newcmp)
    
    cbar=fig.colorbar(cb, extend='both',cax=ax_cb,ax=ax3,orientation='vertical' ,label="")
    
    label_format = r'${:.2f}$'
    ticks_loc =np.linspace(min(mdat[:,1]),max(mdat[:,1]),6)#cbar.ax.get_yticks().tolist()

    cbar.ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar.ax.set_yticklabels([label_format.format(x) for x in ticks_loc],fontsize=15)
    cbar.ax.set_ylabel(r"$\langle CV \rangle$",fontsize=15)


    xlist=np.linspace(min(cof_var),max(cof_var),5)
    ylist=np.linspace(min(mean_network),max(mean_network),5)
    
    right_inset_ax.xaxis.set_major_locator(mticker.FixedLocator(xlist))
    right_inset_ax.set_xticklabels([label_format.format(x) for x in xlist],fontsize=15)

    right_inset_ax.yaxis.set_major_locator(mticker.FixedLocator(ylist))
    right_inset_ax.set_yticklabels([label_format.format(x) for x in ylist],fontsize=15)

    right_inset_ax.set_xlabel(r"$\langle CV \rangle$",fontsize=20)
    right_inset_ax.set_ylabel(r"$\langle\mathcal{W}\rangle$",fontsize=20)
    #ax3.set_xscale("log")
    #ax3.set_yscale("log")
    plt.tight_layout()
    plt.show()

def plot_weights_networks():
    '''
    fraccion de nodos en funcion del 
    tiempo color en el cv
    '''
    fig = plt.figure(figsize=(12, 4))
    gs = fig.add_gridspec(1, 3)
    mxrango = slice(14,22)

    custom_params = {"axes.spines.right": False, "axes.spines.top": False}       
    sns.set_theme(font_scale=1.2,style="ticks",rc=custom_params)
    
    min_cv=[3.6,3.7,3.8]

    for o,xn in enumerate([30,40,50]): 
        dataExp=[
                "Mar04/"+str(xn)+"sMar04","Mar10/"+str(xn)+"sMar10",
                "Jan14/"+str(xn)+"sJan14","Jan21/"+str(xn)+"sJan21",
                "Jan29/"+str(xn)+"sJan29","Dez20/"+str(xn)+"sDez20",
                "Mar07/"+str(xn)+"sMar07"
                ]
    
        axhx = fig.add_subplot(gs[0,o])
        
        for zexp in dataExp:
            nsubdir = zexp[slice(6,15)] 
            axdir=os.listdir('mcporcv/'+zexp)
            fdat=np.loadtxt("datMetx/"+str(nsubdir)+".txt")#datos empiricos
            fildir= axdir[0].split(str(int(axdir[0][mxrango].split('.')[0]))+'.txt')[0]
            thrho=0#np.mean(fdat[:,4])-2*np.std(fdat[:,4])# minimo de Desdidad
            fdiverg=delete_nod(np.loadtxt("retratostxt/"+str(nsubdir)+"retr.txt"))#datos empiricos
            
            psMx=[]; mean_network=[]; dtao=[]; divergen=[]

            for i in range(len(fdat[:,0])):

                zdir='mcporcv/'+zexp+"/"+fildir+str(int(fdat[i,0]))+'.txt'
                xmcc = abs(delete_nod(np.genfromtxt(zdir)))
                
                xcv=nx.from_numpy_array(xmcc)
                dzx=nx.density(xcv)

                if(dzx > thrho):
                    dtao.append(fdat[i,0])
                    mxc=xmcc[xmcc>0]
                    mean_network.append(np.mean(np.triu(mxc).flatten()))
                else:
                    psMx.append(fdat[i,0])

            cof_var=[]

            for zo in dtao:
                for j,ss in enumerate(fdat[:,0]):
                    if(int(zo)==int(ss)):
                        cof_var.append(fdat[int(j),1])

            for fil in range(len(fdat[:,0])):
                divergen.append(np.mean(fdiverg[fil:,]))

            cof_var.sort()

            df = {'x': cof_var, 'y': mean_network}
            sns.scatterplot(data = df,x = "x", y = "y",marker = "o",legend=False)#,palette=newcmp,legend=False,ax=right_inset_ax)
        
        label_format = r'${:.1f}$'
        label_format_ylist = r'${:.2f}$'

        xlist=np.linspace(0,min_cv[o],6)
        ylist=np.linspace(0,0.17,6)
        print(ylist)
        axhx.xaxis.set_major_locator(mticker.FixedLocator(xlist))
        axhx.set_xticklabels([label_format.format(x) for x in xlist],fontsize=14)

        axhx.yaxis.set_major_locator(mticker.FixedLocator(ylist))
        axhx.set_yticklabels([label_format_ylist.format(x) for x in ylist],fontsize=14)
        
        axhx.set_ylabel(r"Q",fontsize=15)
        axhx.yaxis.set_tick_params(labelsize=15)

        axhx.set_xlabel(r"$CV$",fontsize=15)
        axhx.xaxis.set_tick_params(labelsize=15)

        axhx.set_ylim(0,0.18)
        axhx.set_xlim(0,min_cv[o]+0.15)
            
        font = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 15,
                }
        plt.title(r'Network portrait divergence  vs. CV, $\Delta t=$ '+str(xn)+' s',fontdict=font)
    plt.tight_layout()
    plt.show()
 
def plot_distribuciones_clasificadas(addres_dir,file_dir,interval,cv_alto,cv_medio,cv_bajo,mdat):
    xlin=0.5
    interlog0=np.linspace(-3.2,-1.5,80)
    interlog1=np.linspace(-3.1,-0.9,80)
    interlog2=np.linspace(-3.1,-0.9,80)
 
    mxrango = slice(14,22)

    fig, axl = plt.subplots(figsize=(7,4))

    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    thrho=np.mean(mdat[:,3])-2*np.std(mdat[:,3])# minimo de Desdidad
  
    nNcbajo=0
    
    for i in range(len(cv_bajo)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_bajo[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1
    
    cv_frec_bajo=np.zeros((int(len(interval)-1),nNcbajo)); yfr_bajo=[];yfr_medio=[];yfr_alto=[]
    xcv_frec_bajo=np.zeros((int(len(interlog0)-1),nNcbajo)); xyfr_bajo=[];xyfr_medio=[];xyfr_alto=[]
    dst_alto=[];dst_bajo=[];dst_medio=[];ox=0

    for i in range(len(cv_bajo)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_bajo[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        
        if(rsz>thrho):
            frec, xbin=np.histogram(xmcc[xmcc>0],bins=interval)
            x=xbin[:-1]; y= frec/sum(frec)
            for yfre in range(len(y)):
                cv_frec_bajo[yfre,ox] = y[yfre] 
            axl.plot(x,y,linewidth=xlin,color="#ac0000")
            ox += 1
    
    for k in range(len(cv_frec_bajo[0,:])):
        yfr_bajo.append(np.mean(cv_frec_bajo[k,:]))
        
    for k in range(len(xcv_frec_bajo[0,:])):
        xyfr_bajo.append(np.mean(xcv_frec_bajo[k,:]))

    for k in range(len(xcv_frec_bajo[0,:])):
        dst_bajo.append(stdev(xcv_frec_bajo[k,:]))

    axl.plot(x,yfr_bajo,"--",linewidth=3,color="black",label="Average data")


    nNcbajo=0

    for i in range(len(cv_medio)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_medio[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1
    
    cv_frec_medio=np.zeros((int(len(interval)-1),nNcbajo))
    xcv_frec_medio=np.zeros((int(len(interlog1)-1),nNcbajo)); 
    
    ox=0
    
    for i in range(len(cv_medio)):
        
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_medio[i])+'.txt')))
        xft=nx.density(nx.from_numpy_array(xmcc))

        if(xft > thrho):
            mxc=xmcc[xmcc>0]
            frec, xbin=np.histogram(mxc,bins=interval)
            x1=xbin[:-1]; y= frec/sum(frec)
            for yfre in range(len(y)):
                cv_frec_medio[yfre,ox] = y[yfre] 
            axl.plot(x1,y,linewidth=xlin,color="#d1af00")
            ox+=1

    for i in range(len(xcv_frec_medio[0,:])):
        if(any(cv_frec_medio[:,i])!=True):
            np.delete(cv_frec_medio, i, axis=1)
            np.delete(xcv_frec_medio, i, axis=1)
        
    for k in range(len(cv_frec_medio[:,0])):
        yfr_medio.append(np.mean(cv_frec_medio[k,:]))
        
    for k in range(len(xcv_frec_medio[:,0])):
        xyfr_medio.append(np.mean(xcv_frec_medio[k,:]))

    for k in range(len(xcv_frec_medio[:,0])):
        dst_medio.append(stdev(xcv_frec_medio[k,:]))
    
    axl.plot(x1,yfr_medio,"--",linewidth=3,color="black",label="Average data")

    nNcbajo=0

    for i in range(len(cv_alto)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_alto[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1

    cv_frec_alto=np.zeros((int(len(interval)-1),nNcbajo))
    xcv_frec_alto=np.zeros((int(len(interlog2)-1),nNcbajo));  ox=0
    
    for i in range(len(cv_alto)):
        
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_alto[i])+'.txt')))
        dzx=nx.density(nx.from_numpy_array(xmcc))
    
        if(dzx > thrho):
            mxc=xmcc[xmcc>0]
            frec, xbin=np.histogram(mxc,bins=interval)
            x2=xbin[:-1]; y= frec/sum(frec);            
            for yfre in range(len(y)):
                cv_frec_alto[yfre,ox] = y[yfre]

            ox+=1

            axl.plot(x2, y, linewidth=xlin, color="#9c46b1")

    for i in range(len(xcv_frec_alto[0,:])):
        if(any(xcv_frec_alto[:,i])!=True):
            np.delete(cv_frec_alto, i, axis=1)
            np.delete(xcv_frec_alto, i, axis=1)

    for k in range(len(cv_frec_alto[:,0])):
        yfr_alto.append(np.mean(cv_frec_alto[k,:]))
        
    for k in range(len(xcv_frec_alto[:,0])):
        xyfr_alto.append(np.mean(xcv_frec_alto[k,:]))

    for k in range(len(xcv_frec_alto[:,0])):
        dst_alto.append(stdev(xcv_frec_alto[k,:]))
    
    axl.plot(x2,yfr_alto,"--",linewidth=3,color="black",label="Average data")
    plt.legend(loc = "upper left",fontsize='large')   
            
    axl.xaxis.set_tick_params(labelsize=15)
    axl.yaxis.set_tick_params(labelsize=15)
    axl.set_xlabel(r"$\mathcal{W}$",fontsize=20)
    axl.set_ylabel(r"$\rho(\mathcal{W})$",fontsize=20)
    
    plt.tight_layout()
    plt.show()
               
def plot_distribuciones_lognormal(addres_dir,nsubdir,file_dir,interval,cv_alto,cv_medio,cv_bajo,mdat):
    xlin=0.5
    interlog0=np.linspace(-3.2,-1.5,80)
    interlog1=np.linspace(-3.1,-0.9,80)
    interlog2=np.linspace(-3.1,-0.9,80)

    mxrango = slice(14,22)
    fig = plt.figure(figsize=(15,5))
    gs = gridspec.GridSpec(1, 3)

    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    thrho=np.mean(mdat[:,3])-2*np.std(mdat[:,3])# minimo de Desdidad
   
    nNcbajo=0
  
    for i in range(len(cv_bajo)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_bajo[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1
    
    cv_frec_bajo=np.zeros((int(len(interval)-1),nNcbajo)); yfr_bajo=[];yfr_medio=[];yfr_alto=[]
    xcv_frec_bajo=np.zeros((int(len(interlog0)-1),nNcbajo)); xyfr_bajo=[];xyfr_medio=[];xyfr_alto=[]
    dst_alto=[];dst_bajo=[];dst_medio=[]

    axy2 = fig.add_subplot(gs[2]); o=0;ox=0

    for i in range(len(cv_bajo)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_bajo[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        
        if(rsz>thrho):
          
            xrfrec, xrxbin=np.histogram(np.log(xmcc[xmcc>0]),bins=interlog0)
           
            zrx=xrxbin[:-1]; zry= xrfrec/sum(xrfrec)
            for yfr in range(len(zry)):
                xcv_frec_bajo[yfr,ox] = zry[yfr]

            axy2.plot(zrx,zry,linewidth=xlin,color="#bababa")
            ox += 1
    
            if(ox==nNcbajo):#esto es para leyenda
                axy2.plot(zrx[0],zry[0],linewidth=2,color="#bababa",label="Experimental data")
                plt.legend(loc = "upper right",fontsize='x-large')

    for k in range(len(cv_frec_bajo[:,0])):
        yfr_bajo.append(np.mean(cv_frec_bajo[k,:]))
        
    for k in range(len(xcv_frec_bajo[:,0])):
        xyfr_bajo.append(np.mean(xcv_frec_bajo[k,:]))

    for k in range(len(xcv_frec_bajo[:,0])):
        dst_bajo.append(stdev(xcv_frec_bajo[k,:]))

    axy2.plot(zrx,xyfr_bajo,"--",linewidth=3,color="#ac0000",label="Average data")
    plt.legend(loc = "upper right",fontsize='x-large')

    nNcbajo=0

    for i in range(len(cv_medio)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_medio[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1
    
    cv_frec_medio=np.zeros((int(len(interval)-1),nNcbajo))
    xcv_frec_medio=np.zeros((int(len(interlog1)-1),nNcbajo)); 
    
    axy1 = fig.add_subplot(gs[1]); ox=0
    
    for i in range(len(cv_medio)):
        
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_medio[i])+'.txt')))
        xft=nx.density(nx.from_numpy_array(xmcc))

        if(xft > thrho):
            mxc=xmcc[xmcc>0]
          
            xrfrec, xrxbin=np.histogram(np.log(mxc),bins=interlog1)
          
            zrx1=xrxbin[:-1]; zry1= xrfrec/sum(xrfrec) 
            
            for yfr in range(len(zry1)):
                xcv_frec_medio[yfr,ox] = zry1[yfr]
            
            axy1.plot(zrx1,zry1,linewidth=xlin,color="#bababa")
            ox+=1

            if(ox==nNcbajo):#esto es para leyenda
                axy1.plot(zrx1[0],zry1[0],linewidth=2,color="#bababa",label="Experimental data")
                plt.legend(loc = "upper right",fontsize='x-large')

    for i in range(len(xcv_frec_medio[0,:])):
        if(any(cv_frec_medio[:,i])!=True):
            np.delete(cv_frec_medio, i, axis=1)
            np.delete(xcv_frec_medio, i, axis=1)
        
    for k in range(len(cv_frec_medio[:,0])):
        yfr_medio.append(np.mean(cv_frec_medio[k,:]))
        
    for k in range(len(xcv_frec_medio[:,0])):
        xyfr_medio.append(np.mean(xcv_frec_medio[k,:]))

    for k in range(len(xcv_frec_medio[:,0])):
        dst_medio.append(stdev(xcv_frec_medio[k,:]))
    
    axy1.plot(zrx1,xyfr_medio,"--",linewidth=3,color="#d1af00",label="Average data") 
    plt.legend(loc = "upper right",fontsize='x-large')

    nNcbajo=0

    for i in range(len(cv_alto)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_alto[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1

    cv_frec_alto=np.zeros((int(len(interval)-1),nNcbajo))
    xcv_frec_alto=np.zeros((int(len(interlog2)-1),nNcbajo)); 
    
    axy0 = fig.add_subplot(gs[0]); ox=0; 
    
    for i in range(len(cv_alto)):
        
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_alto[i])+'.txt')))
        dzx=nx.density(nx.from_numpy_array(xmcc))
    
        if(dzx > thrho):
            mxc=xmcc[xmcc>0]
            xrfrec, xrxbin=np.histogram(np.log(mxc),bins=interlog2)
           
            zrx2=xrxbin[:-1]; zry2= xrfrec/sum(xrfrec)
                
            for yfr in range(len(zry2)):
                xcv_frec_alto[yfr,ox] = zry2[yfr]

            ox+=1

           
            axy0.plot(zrx2,zry2,linewidth=xlin,color="#bababa")
            
            if(ox==nNcbajo):#esto es para leyenda
                axy0.plot(zrx2[0],zry2[0],linewidth=2,color="#bababa",label="Experimental data")
                plt.legend(loc = "upper left",fontsize='x-large')

    for i in range(len(xcv_frec_alto[0,:])):
        if(any(xcv_frec_alto[:,i])!=True):
            np.delete(cv_frec_alto, i, axis=1)
            np.delete(xcv_frec_alto, i, axis=1)

    for k in range(len(cv_frec_alto[:,0])):
        yfr_alto.append(np.mean(cv_frec_alto[k,:]))
        
    for k in range(len(xcv_frec_alto[:,0])):
        xyfr_alto.append(np.mean(xcv_frec_alto[k,:]))

    for k in range(len(xcv_frec_alto[:,0])):
        dst_alto.append(stdev(xcv_frec_alto[k,:]))
    
    plt.legend(loc = "upper left",fontsize='x-large')  
    axy0.plot(zrx2,xyfr_alto,"--",linewidth=3,color="#9c46b1",label="Average data") 
    plt.legend(loc = "upper left",fontsize='x-large') 



    axy0.xaxis.set_tick_params(labelsize=15)
    axy0.yaxis.set_tick_params(labelsize=15)
    axy0.set_xlabel(r"$\Omega(\mathcal{W})$",fontsize=20)
    axy0.set_ylabel(r"$\rho(\Omega)$",fontsize=20)
    axy0.yaxis.set_data_interval(0, 0.071 ,True)

    axy1.xaxis.set_tick_params(labelsize=15)
    axy1.yaxis.set_tick_params(labelsize=15)
    axy1.set_xlabel(r"$\Omega(\mathcal{W})$",fontsize=20)
    axy1.set_ylabel(r"$\rho(\Omega)$",fontsize=20)
    axy1.yaxis.set_data_interval(0, 0.071 ,True)

    axy2.xaxis.set_tick_params(labelsize=15)
    axy2.yaxis.set_tick_params(labelsize=15)
    axy2.set_xlabel(r"$\Omega(\mathcal{W})$",fontsize=20)
    axy2.set_ylabel(r"$\rho(\Omega)$",fontsize=20)
    axy2.yaxis.set_data_interval(0, 0.071 ,True)

    plt.tight_layout()
    
    plt.show()
    
    bfile_txt=nsubdir+"_bajo.txt"
    bfdata = open(bfile_txt, 'w')
    
    try:
        for i in range(len(xyfr_bajo)):
            bfdata.write(str(xyfr_bajo[i])+"\t"+str(zrx[i])+"\t"+str(dst_bajo[i])+"\t"+str(0)+"\n")      
    finally:
        bfdata.close()


    file_txt=nsubdir+"_medio.txt"
    fdata = open(file_txt, 'w')
    
    try:
        for i in range(len(xyfr_medio)):
            fdata.write(str(xyfr_medio[i])+"\t"+str(zrx1[i])+"\t"+str(dst_medio[i])+"\t"+str(0)+"\n")      
    finally:
        fdata.close()

    fxile_txt=nsubdir+"_alto.txt"
    fxdata = open(fxile_txt, 'w')
    
    try:
        for i in range(len(xyfr_alto)):
            fxdata.write(str(xyfr_alto[i])+"\t"+str(zrx2[i])+"\t"+str(dst_alto[i])+"\t"+str(0)+"\n")      
    finally:
        fxdata.close()

def similitud_cv(cvs):
    
    interval=np.arange(0.01,0.3,0.001)
    c=np.zeros((len(cvs),len(cvs)))
    cv=np.zeros((len(cvs),len(cvs)))

    for i in range(len(cvs)-1):
        for j in range(i+1,len(cvs)):
            zs=abs(cvs[i]-cvs[j])
            c[i,j] = zs
            c[j,i] = zs
    
    return c


