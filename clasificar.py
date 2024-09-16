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
def fun_ajuste(x,mu,sg,c):
    fn=(c*(sg*math.sqrt(2*np.pi)))*np.exp(-((x-mu)**2/(2*sg**2)))
    return fn
    
def delete_nod(M):

    if number_of_components(np.abs(M)) > 1:

        #print("Matriz desconectada \t")

        for j in range(len(M)):
                
            for i in range(len(M)):
                if(any(M[:,int(i)])!=True):
                    A=np.delete(M, i, axis=1)
                    break  
                
            for i in range(len(M)):    
                if(any(A[int(i),:])!=True):
                    B=np.delete(A, i, axis=0) 
                    break
        
            M = np.resize(M,(len(B),len(B)))
            M=B; contador=-1
            
            for a_i in range(len(M)):   
                if(any(M[int(a_i),:])==True):
                    contador+=1
            
            if(a_i==contador):
                break
        
   
        return M
    else:
        
        #print("Matriz conectada \t") 
        
        return M

def dist_psCV(addres_dir,nsubdir,file_dir,interval,cv_alto,cv_medio,cv_bajo,mdat,idcv,cv_otras,comunidad_tres=False,comunidad_otras=False):
    xlin=0.5
    interlog0=np.linspace(-3.5,-0.8,80)
    interlog1=np.linspace(-3.4,-1.8,80)
    interlog2=np.linspace(-3.5,-1,80)

    #cmap = plt.colormaps["viridis"]#plasma
    viridis = mpl.colormaps['plasma'].resampled(len(idcv))
    newcolors = viridis(np.linspace(0,1, len(idcv)))
    newcmp = ListedColormap(newcolors)
    #fig, ax = plt.subplots()
    #fig.suptitle(" ", size=14)
    mxrango = slice(14,22)
    fig = plt.figure(figsize=(13,8))
    gs = gridspec.GridSpec(2, 3)
    #axc = figx.add_subplot()
    axx = fig.add_subplot(gs[0,:2])
    axl = fig.add_subplot(gs[0, 2])
    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    print(fildir)
    thrho=np.mean(mdat[:,4])-2*np.std(mdat[:,4])# minimo de Desdidad
   
    #if(continuo):
    psMx=[]

    for i in range(len(idcv)):
        
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(int(idcv[i]))+'.txt')))
        xcv=nx.from_numpy_array(xmcc)
        dzx=nx.density(xcv)


        if(dzx > thrho):
            mxc=xmcc[xmcc>0]
            frec, xbin=np.histogram(mxc,bins=interval)
            x=xbin[:-1]; y= frec/sum(frec)
            axx.plot(x,y,linewidth=xlin, color=viridis.colors[i-1])
        else:
            psMx.append(idcv[i])
    
    # xDaux[psMx,:]=0
    # xDaux[:,psMx]=0
    #xDux=abs(delete_nod(np.array(xDaux)))
    
    axx.xaxis.set_tick_params(labelsize=15)
    axx.yaxis.set_tick_params(labelsize=15)
    #cbar.ax.set_ylabel(r"$CV$",fontsize=15)
    axx.set_xlabel(r"$\mathcal{W}$",fontsize=20)
    axx.set_ylabel(r"$\rho(\mathcal{W})$",fontsize=20)
    
    frec, xbin=np.histogram(mdat[:,4][mdat[:,4]>thrho],bins=20)
    xo=xbin[:-1]; yo= frec/sum(frec)
    viridis = mpl.colormaps['plasma'].resampled(len(yo))
    newcolors = viridis(np.linspace(0,1,len(mdat[:,1])))
    newcmp = ListedColormap(newcolors)

    right_inset_ax = fig.add_axes([0.447, 0.76, 0.15, 0.22], facecolor="None")
  
    for row in range(len(xo)):
        right_inset_ax.bar(xo[row], yo[row], width=.01,ec='w',color=viridis.colors[row]) 

    cb=plt.cm.ScalarMappable(norm=Normalize(min(mdat[:,1]),max(mdat[:,1]) ), cmap=newcmp)
    cax = plt.axes([0.6, 0.76, 0.01, 0.22])
    cbar=fig.colorbar(cb, cax=cax,ax=cax, label="")
    label_format = r'${:.2f}$'
    ticks_loc =cbar.ax.get_yticks().tolist()

    xlist=np.linspace(min(xo),max(xo),5)
    ylist=np.linspace(min(yo),max(yo),5)

    cbar.ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar.ax.set_yticklabels([label_format.format(x) for x in ticks_loc],fontsize=14)
    cbar.ax.set_ylabel(r"$CV$",fontsize=15)

    right_inset_ax.xaxis.set_major_locator(mticker.FixedLocator(xlist))
    right_inset_ax.set_xticklabels([label_format.format(x) for x in xlist],fontsize=14)

    right_inset_ax.yaxis.set_major_locator(mticker.FixedLocator(ylist))
    right_inset_ax.set_yticklabels([label_format.format(x) for x in ylist],fontsize=14)

    right_inset_ax.set_xlabel(r"$\sigma_{\textup{network}}$",fontsize=20)
    right_inset_ax.set_ylabel(r"$\rho(\sigma)$",fontsize=20)

    '''
    right_inset_ax = fig.add_axes([0.42, 0.73, 0.24, 0.24], facecolor='k')
    pc =right_inset_ax.imshow(xDux, cmap=newcmp,vmin=0, vmax=cxmax)
    divider = make_axes_locatable(right_inset_ax)
    cax = divider.append_axes("right" ,size="7%", pad=0.04)
    cbar=fig.colorbar(pc, cax=cax)
    
    label_format = r'${:.1f}$'
    ticks_loc =cbar.ax.get_yticks().tolist()
    
    cbar.ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar.ax.set_yticklabels([label_format.format(x) for x in ticks_loc],fontsize=13)
    
    xlist=np.linspace(0,262,7)
    rlt=[int(ix) for ix in xlist]
    
    right_inset_ax.xaxis.set_tick_params(labelsize=13) 
    right_inset_ax.xaxis.set_major_locator(mticker.FixedLocator(rlt)) 
        
    right_inset_ax.yaxis.set_tick_params(labelsize=13)
    right_inset_ax.yaxis.set_major_locator(mticker.FixedLocator(rlt))
    
    cbar.ax.set_ylabel(r"$\mathcal{D}_{euc}(CV_i,CV_j)$",fontsize=13)
    right_inset_ax.set_xlabel(r"$CV_j$",fontsize=13)
    right_inset_ax.set_ylabel(r"$CV_i$",fontsize=13)
    '''
    #else:  

    nNcbajo=0
    distrib_dat=np.zeros((len(interlog2)*len(cv_alto),2))
    
    for i in range(len(cv_bajo)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_bajo[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1
    
    cv_frec_bajo=np.zeros((int(len(interval)-1),nNcbajo)); yfr_bajo=[];yfr_medio=[];yfr_alto=[]
    xcv_frec_bajo=np.zeros((int(len(interlog0)-1),nNcbajo)); xyfr_bajo=[];xyfr_medio=[];xyfr_alto=[]
    dst_alto=[];dst_bajo=[];dst_medio=[]

    
    axy2 = fig.add_subplot(gs[1, 2]); o=0;ox=0

    for i in range(len(cv_bajo)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_bajo[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        
        if(rsz>thrho):
            frec, xbin=np.histogram(xmcc[xmcc>0],bins=interval)
            xrfrec, xrxbin=np.histogram(np.log(xmcc[xmcc>0]),bins=interlog0)
            x=xbin[:-1]; y= frec/sum(frec)
            zrx=xrxbin[:-1]; zry= xrfrec/sum(xrfrec)
            
            for yfre in range(len(y)):
                cv_frec_bajo[yfre,ox] = y[yfre] 

            for yfr in range(len(zry)):
                xcv_frec_bajo[yfr,ox] = zry[yfr]
                            
            axl.plot(x,y,linewidth=xlin,color="#9c46b1")
            axy2.plot(zrx,zry,linewidth=xlin,color="#bababa")
            ox += 1
    
            if(ox==nNcbajo):#esto es para leyenda
                axl.plot(x[0],y[0],linewidth=2,color="#9c46b1",label="Experimental data")
                plt.legend(loc = "upper right")
                axy2.plot(zrx[0],zry[0],linewidth=2,color="#bababa",label="Experimental data")
                plt.legend(loc = "upper right",ncols=2,fontsize='large')

    for k in range(len(cv_frec_bajo[:,0])):
        yfr_bajo.append(np.mean(cv_frec_bajo[k,:]))
        
    for k in range(len(xcv_frec_bajo[:,0])):
        xyfr_bajo.append(np.mean(xcv_frec_bajo[k,:]))

    for k in range(len(xcv_frec_bajo[:,0])):
        dst_bajo.append(stdev(xcv_frec_bajo[k,:]))

    axl.plot(x,yfr_bajo,"--",linewidth=3,color="black",label="Average data")
    plt.legend(loc = "upper right",fontsize='large')
    axy2.plot(zrx,xyfr_bajo,"--",linewidth=3,color="#9c46b1",label="Average data")
    plt.legend(loc = "upper right",ncols=2,fontsize='large')

    nNcbajo=0

    for i in range(len(cv_medio)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_medio[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1
    
    cv_frec_medio=np.zeros((int(len(interval)-1),nNcbajo))
    xcv_frec_medio=np.zeros((int(len(interlog1)-1),nNcbajo)); 
    
    axy1 = fig.add_subplot(gs[1, 1]); ox=0
    
    for i in range(len(cv_medio)):
        
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_medio[i])+'.txt')))
        xft=nx.density(nx.from_numpy_array(xmcc))

        if(xft > thrho):
            mxc=xmcc[xmcc>0]
            frec, xbin=np.histogram(mxc,bins=interval)
            xrfrec, xrxbin=np.histogram(np.log(mxc),bins=interlog1)
            x1=xbin[:-1]; y= frec/sum(frec)
            zrx1=xrxbin[:-1]; zry1= xrfrec/sum(xrfrec)
            
            for yfre in range(len(y)):
                cv_frec_medio[yfre,ox] = y[yfre] 
            
            for yfr in range(len(zry1)):
                xcv_frec_medio[yfr,ox] = zry1[yfr]
            
            axl.plot(x1,y,linewidth=xlin,color="#d1af00")
            axy1.plot(zrx1,zry1,linewidth=xlin,color="#bababa")
            ox+=1

            if(ox==nNcbajo):#esto es para leyenda
                axy1.plot(zrx1[0],zry1[0],linewidth=2,color="#bababa",label="Experimental data")
                plt.legend(loc = "upper right",ncols=2,fontsize='large')

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
    axy1.plot(zrx1,xyfr_medio,"--",linewidth=3,color="#d1af00",label="Average data") 
    plt.legend(loc = "upper right",ncols=2,fontsize='large')

    nNcbajo=0

    axy0 = fig.add_subplot(gs[1, 0]); ox=0; 

    if(comunidad_tres):
        
        for i in range(len(cv_alto)):
            xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_alto[i])+'.txt')))
            rsz=nx.density(nx.from_numpy_array(xmcc))
            if(rsz>thrho):
                nNcbajo+=1

        cv_frec_alto=np.zeros((int(len(interval)-1),nNcbajo))
        xcv_frec_alto=np.zeros((int(len(interlog2)-1),nNcbajo)); 
        
        for i in range(len(cv_alto)):
            
            xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_alto[i])+'.txt')))
            dzx=nx.density(nx.from_numpy_array(xmcc))
        
            if(dzx > thrho):
                mxc=xmcc[xmcc>0]
                frec, xbin=np.histogram(mxc,bins=interval)
                xrfrec, xrxbin=np.histogram(np.log(mxc),bins=interlog2)
                x2=xbin[:-1]; y= frec/sum(frec);
                zrx2=xrxbin[:-1]; zry2= xrfrec/sum(xrfrec)
                
                for yfre in range(len(y)):
                    cv_frec_alto[yfre,ox] = y[yfre]
                    
                for yfr in range(len(zry2)):
                    xcv_frec_alto[yfr,ox] = zry2[yfr]
                    distrib_dat[o,0]=zrx2[yfr]
                    distrib_dat[o,1]=zry2[yfr]
                    o+=1
                ox+=1
    
                axl.plot(x2, y, linewidth=xlin, color="#ac0000")
                axy0.plot(zrx2,zry2,linewidth=xlin,color="#bababa")
                
                if(ox==nNcbajo):#esto es para leyenda
                    axy0.plot(zrx2[0],zry2[0],linewidth=2,color="#bababa",label="Experimental data")
                    plt.legend(loc = "upper left",ncols=2,fontsize='large')

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
        plt.legend(loc = "upper left",ncols=2,fontsize='large')  
        axy0.plot(zrx2,xyfr_alto,"--",linewidth=3,color="#ac0000",label="Average data") 
        plt.legend(loc = "upper left",ncols=2,fontsize='large')   

    if(comunidad_otras):
        for i in range(len(cv_otras)):
            
            xmcc = abs(delete_nod(np.genfromtxt(fildir+str(cv_otras[i])+'.txt')))
            dzx=nx.density(nx.from_numpy_array(xmcc))
        
            if(dzx > thrho):
                mxc=xmcc[xmcc>0]
                frec, xbin=np.histogram(mxc,bins=interval)
                x=xbin[:-1]; y= frec/sum(frec)
                axx.plot(x2, y, linewidth=0.5, color="black")
             
    axl.xaxis.set_tick_params(labelsize=15)
    axl.yaxis.set_tick_params(labelsize=15)
    axl.set_xlabel(r"$\mathcal{W}$",fontsize=20)
    axl.set_ylabel(r"$\rho(\mathcal{W})$",fontsize=20)
    axl.xaxis.set_data_interval(0.025, 0.35 ,True)
    axl.set_xlim(0.025,0.35)

    l_fr = r'${:.2f}$'
    xl_ist=np.linspace(0.025,0.35,6)
    axl.xaxis.set_major_locator(mticker.FixedLocator(xl_ist))
    axl.set_xticklabels([l_fr.format(x) for x in xl_ist],fontsize=14)

    axy0.xaxis.set_tick_params(labelsize=15)
    axy0.yaxis.set_tick_params(labelsize=15)
    axy0.set_xlabel(r"$\Omega(\mathcal{W})$",fontsize=20)
    axy0.set_ylabel(r"$\rho(\Omega)$",fontsize=20)
    axy0.yaxis.set_data_interval(0, 0.065 ,True)
   
    axy1.xaxis.set_tick_params(labelsize=15)
    axy1.yaxis.set_tick_params(labelsize=15)
    axy1.set_xlabel(r"$\Omega(\mathcal{W})$",fontsize=20)
    axy1.set_ylabel(r"$\rho(\Omega)$",fontsize=20)
    axy1.yaxis.set_data_interval(0, 0.065 ,True)
    

    axy2.xaxis.set_tick_params(labelsize=15)
    axy2.yaxis.set_tick_params(labelsize=15)
    axy2.set_xlabel(r"$\Omega(\mathcal{W})$",fontsize=20)
    axy2.set_ylabel(r"$\rho(\Omega)$",fontsize=20)
    axy2.yaxis.set_data_interval(0, 0.065 ,True)
   

    # xvz=np.array(xDaux)
    # norm = matplotlib.colors.Normalize(vmin=np.min(xDaux), vmax=np.max(xDaux)) 
    # sxm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm) 
    # divider = make_axes_locatable(axx)
    # cax = divider.append_axes("right" ,size="4%", pad=0.04)
    #cbar=fig.colorbar(sxm, cax=cax)
    plt.tight_layout()
    
    plt.show()
    
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
        
def dist_pesos(addres_dir,file_dir,min_nNodos,interval):
    
    #fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))  
    mxrango = slice(14,22)
    sub_cadena=[]
    
    nNum_inter=[]
    
    nNfiles=len(addres_dir)

    
    for i in range(nNfiles):
        sub_cadena.append(int(addres_dir[i][mxrango].split('.')[0]))

    kz=0
    for i in range(nNfiles):
        xm=np.array(delete_nod(np.genfromtxt(file_dir+"_"+str(sub_cadena[i])+'.txt')))
        if(len(xm)>=min_nNodos):
            kz+=1
    
    mxz=np.zeros((kz,len(interval)))
    
    contador=0
    
    sub_cadena.sort()
    
    for i in range(nNfiles):

        xmcc = np.genfromtxt(file_dir+"_"+str(sub_cadena[i])+'.txt')
        
        xmcnect =np.array(delete_nod(xmcc))

        if(len(xmcnect)>=min_nNodos):
            
            #print(file_dir+"_"+str(sub_cadena[i])+'.txt')
            
            nNum_inter.append(int(sub_cadena[i]))
            
            mxCC = np.triu(xmcnect).flatten()
            
            mxc=[]   
            
            for k in range(len(mxCC)):
                if(mxCC[k]!=0):
                    mxc.append(abs(mxCC[k]))
            
            
            frec, xbin=np.histogram(mxc,bins=interval)
            x=xbin[:-1]; y= frec/sum(frec);
            #print(len(x),' ',len(y))
            for ox in range(len(y)):
                mxz[contador,ox]=y[ox]+0.00001;
            contador+=1
            
            plt.plot(x,y,linewidth=0.5)
    
    plt.show()
    return nNum_inter
    
def  draw_comunidades(vxy,nam_exp,zname,mostrar=False):
    
    ffile_txt="comu_correlacion/"+nam_exp+".txt"
    
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
        file_txt=zname+nam_exp+"ok.txt"
        fdata = open(file_txt, 'w')
        
        for cluster_id in nNCom:
            cluster = [node for node in G_fb.nodes() if partitions[node] == cluster_id]
            nodos_por_comundad.append(cluster)
            for wrfil in range(len(cluster)):
                fdata.write(str(cluster_id)+"\t"+str(cluster[wrfil])+"\n")#+str(int(sdt[xo]))+"\n")
                xo+=1
             
    finally:
        fdata.close()
    
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
                    xD_kls[i,j]=1#D_kl[i,j]
                    xD_kls[j,i]=1#D_kl[i,j]
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
        treshold*=1.08
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
        plt.plot(pdata,gdata,"-")
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
