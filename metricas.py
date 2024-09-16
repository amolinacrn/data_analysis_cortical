# -*- coding: utf-8 -*-
from cProfile import label
from lib2to3.pytree import HUGE
from operator import imod
from statistics import stdev
from turtle import fd, shape
import matplotlib

import numpy as np
import os
import numpy 
import math
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import seaborn as sns
from include.get_function_implementations import *
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
from NAT2 import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from include.swp.small_world_propensity import *
from include.swp.modelos_nulos_SWP import *
import matplotlib.ticker as mticker
from portrait_Djs.funct_PDjs import arbole_expansion_min
#from clasificar import componente_gigante

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

def ajuste_distribucion_metricas(dt,cv,file_dir,addres_dir,nam_experiment):

    interv_lpath=np.linspace(0,81,81)
    
    mxrango = slice(14,22)

    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    
    maxd_lp=[]; maxd_bins=[]#np.zeros((len(dt),len(interv_lpath)-1))

    for i in range(len(dt)):
        
        xm=abs(delete_nod(np.genfromtxt(fildir+str(int(dt[i]))+'.txt')))
        
        
        # io,jo=np.where(xm!=0)
        # xm[io,jo]=1
        Gx = nx.from_numpy_array(xm)
              
        Gtc = {(e1, e2): 1 / abs(weight) for e1, e2, weight in Gx.edges(data='weight')}
        nx.set_edge_attributes(Gx, Gtc, 'distance')
        
        vector_k_lover,nkover,err_kover= mean_degree_network(xm)
        
        lp=dict(nx.all_pairs_dijkstra_path_length(Gx,weight='distance'))
        
        log_path=[]; matz_lp=[]
        
        for os in range(len(lp)):
            matz_lp.append(list(lp[os].values()))
            log_path += list(lp[os].values())  
        
        clustering = list(nx.clustering(Gx, weight='weight').values())
        xlpz=np.array(log_path)
        xlp=list(xlpz[xlpz!=0])

        sx,sy=fun_histogram(xlp,interv_lpath)
       
        maxd_lp.append(sy) 
   
    matriz_desity_lp=np.transpose(np.array(maxd_lp))
    ko,lo=np.where(matriz_desity_lp<0.07)
    matriz_desity_lp[ko,lo]=0

    file_txt="metrica_"+nam_experiment+".txt"
    fdata = open(file_txt, 'w')
    print(np.shape(matriz_desity_lp))
    try:    
        for io in range(len(matriz_desity_lp[:,0])):
            for jo in range(len(matriz_desity_lp[0,:])):
                fdata.write(str(cv[jo])+"\t"+str(matriz_desity_lp[io,jo])+"\n")      
    finally:
        fdata.close()

    

def mibarcax(cv,pargraf,cmap,label_format,ticks_loc,xcvd):
    norm= matplotlib.colors.Normalize(vmin=min(cv), vmax=max(cv)) 
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm) 
    ax2_divider = make_axes_locatable(pargraf)
    cax0 = ax2_divider.append_axes("bottom", size="4%", pad="35%")
    cbar0=fig.colorbar(sm,extend='both' ,cax=cax0,orientation='horizontal')
    cbar0.ax.set_xlabel(xcvd,fontsize=20)   
    cbar0.ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar0.ax.set_xticklabels([label_format.format(x) for x in ticks_loc],fontsize=13)
    cax0.xaxis.set_ticks_position("bottom")

def fun_histogram(datos_histo,interval):
    fx, xb=np.histogram(datos_histo,bins=interval)
    sx=xb[:-1]; sy= fx/sum(fx)
    return sx, sy
   
def plot_metricas(k,pargraf,dt,cv,file_dir,addres_dir,nam_experiment):
    xmicol0="plasma"
    xmicol1="afmhot"#hot,afmhot
    altura_abajo=0.12
    #ft = 0, 1, 2. cl, lp ,k 
    ft=1; eje_name=r"$\epsilon_{i}$"; 
    xeje_name=r"$\rho(\epsilon_{i})$"
    zeje_name=r"$\langle CV \rangle$"

    #EN Y FILA 0
    yfrmt = r'${:.2f}$'

    # EN X FILA 0
    lformat = r'${:.2f}$'

    #EN X FILA 1
    xfrmt = r'${:.2f}$'

    #EN Y FILA 1
    lfoy = r'${:.2f}$'
    
    # cmap = plt.colormaps[xmicol0]
    #xcmap = plt.colormaps[xmicol1]
    
    cmap = mpl.colormaps[xmicol0].resampled(len(dt))
    newcolors = cmap(np.linspace(0,1,len(dt)))
    newcmp = ListedColormap(newcolors)

    xcmap = mpl.colormaps[xmicol1].resampled(len(dt))
    xnewcolors = xcmap(np.linspace(0,1,len(dt)))
    xnewcmp = ListedColormap(xnewcolors)

    interv_lpath=np.linspace(0,0.1,70)
    interv_clust=np.linspace(0,0.6,50)
    interv_grade=np.linspace(0,15,50)
    interv_witges=np.linspace(0,.25,50)
    
    mxrango = slice(14,22)

    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    
    maxd_lp=[]; maxd_bins=[]#np.zeros((len(dt),len(interv_lpath)-1))

    for i in range(len(dt)):
        
        mspikes=abs(delete_nod(np.genfromtxt(fildir+str(int(dt[i]))+'.txt')))
        xm=np.zeros((len(mspikes[:,0]),len(mspikes[0,:])))
        io,jo=np.where(mspikes)
        xm[io,jo]=mspikes[io,jo]
        #print(file_dir+"_"+str(int(dt[i]))+'.txt')
        # io,jo=np.where(xm!=0)
        # xm[io,jo]=1
        Gx = nx.from_numpy_array(xm)
               
        Gtc = {(e1, e2): abs(weight) for e1, e2, weight in Gx.edges(data='weight')}
        nx.set_edge_attributes(Gx, Gtc, 'distance')
        
        vector_k_lover,nkover,err_kover= mean_degree_network(xm)
        
        lp=dict(nx.all_pairs_dijkstra_path_length(Gx,weight='distance'))
        
        log_path=[]; matz_lp=[]
        
        for os in range(len(lp)):
            matz_lp.append(list(lp[os].values()))
            log_path += list(lp[os].values())  

        clustering = list(nx.clustering(Gx, weight='weight').values())
        xlpz=np.array(log_path)
        xlp=list(xlpz[xlpz!=0])

        if(ft==0):
            sx,sy=fun_histogram(clustering,interv_clust)
            maxY0=len(interv_clust)-2
        if(ft==1):
            sx,sy=fun_histogram(xlp,interv_lpath)
            maxY0=len(interv_lpath)-2
        if(ft==2):
            sx,sy=fun_histogram(vector_k_lover,interv_grade)
            maxY0=len(interv_grade)-2

        if(ft==3):
            sx,sy=fun_histogram(xm[xm!=0],interv_witges)
            maxY0=len(interv_witges)-2

        #----------------------------------------------1
        
        maxd_lp.append(sy) 
        maxd_bins.append(sx)   
        pargraf[1].plot(sx,sy,linewidth=0.3, color=newcmp(i))
       
        # ay0.imshow(matz_path_lengt,cmap=cmap,
        #     vmin=np.min(matz_path_lengt[matz_path_lengt>0]),
        #     vmax=np.max(matz_path_lengt),
        #     aspect='auto')#plot(rx,ry,linewidth=0.7, color=cmap(i))
    mz_maxd_bins=np.array(maxd_bins).flatten()
   
    matriz_desity_lp=np.array(maxd_lp)

    xvc=matriz_desity_lp.flatten()

    
    # normx= matplotlib.colors.Normalize(vmin=min(xvc), vmax=max(xvc)) 
    # smx = matplotlib.cm.ScalarMappable(cmap=xcmap, norm=normx)

    pargraf[0].imshow(np.transpose(matriz_desity_lp),cmap=xnewcmp,aspect='auto')
    pargraf[0].set_title(nam_experiment,fontsize="x-large")


    #-----------------------------------------------
    minX1=0;  minY1=0; minX0=0; minY0=0
    maxX1=max(mz_maxd_bins); maxY1=altura_abajo; maxX0=len(cv)-1; 
    zmaxY1=max(xvc)

    pargraf[0].invert_yaxis()

    pargraf[1].set_ylim(minY1,maxY1)
    pargraf[1].set_xlim(minX1,maxX1)
    
    pargraf[0].set_ylim(minY0,maxY0)
    pargraf[0].set_xlim(minX0,maxX0)

    pargraf[0].xaxis.set_tick_params(labelsize=13)
    pargraf[0].set_xlabel(zeje_name,fontsize=20)
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
        pargraf[0].set_ylabel(eje_name,fontsize=30)
        pargraf[1].yaxis.set_tick_params(labelsize=13)
        pargraf[1].set_ylabel(xeje_name,fontsize=20)

    else:
        pargraf[0].tick_params(axis="y", labelleft=False)
        pargraf[1].tick_params(axis="y", labelleft=False)

    mibarcax(xvc,pargraf[0],xmicol1,r'${:.2f}$',np.linspace(minY1,zmaxY1,5),r"$\rho(\epsilon_{i})$")
    mibarcax(cv,pargraf[1],xmicol0,r'${:.1f}$',np.linspace(min(cv),max(cv),6),r"$CV$")
   

def promedio_error_metricas(addres_dir,filexdir,min_nNodos=10):
    '''
    retorna:numero de nodos, 
            intervalo de tiempo
            densidad de red
            grado medio de red
            propension mundo pequeño
            longitud de camino
            clustering
    '''
    #fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))  
    mxrango = slice(14,22)
    sub_cadena=[]
   
    nNfiles=len(addres_dir); 
    for i in range(nNfiles):
        sub_cadena.append(int(addres_dir[i][mxrango].split('.')[0]))
    
    file_dir= filexdir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]

    sub_cadena.sort()

    for i in range(1):#nNfiles):
        datx=np.array(abs(delete_nod(np.genfromtxt(file_dir+str(sub_cadena[i])+'.txt'))))
        xGz = nx.from_numpy_array(datx)
        if((len(xGz.nodes)>min_nNodos)):
            nN=len(xGz.nodes)
            Gtc = {(e1, e2): 1 / abs(weight) for e1, e2, weight in xGz.edges(data='weight')}
            nx.set_edge_attributes(xGz, Gtc, 'distance')

            lp=dict(nx.all_pairs_dijkstra_path_length(xGz,weight='distance'))
            vector_k_lover, nkover= mean_degree_network(datx)
            clus=list(nx.clustering(xGz, weight='weight').values())
            #pl, cl = metricas_L_C(datx); log_path=[]
            
            # for os in range(len(clus)):
            #     log_path += list(clus[os].values())

            # # log_path=np.array(log_path)[np.array(log_path)!=0]
            # # stdev_pat_lengt=np.mean(log_path)

            # clustering = nx.average_clustering(xGz, weight='weight')
            # mean_clust = np.sum(clus)/len(clus)
            strength = dict(xGz.degree(weight='weight'))
            xgrad=list(strength.values())
            print(xgrad)

def calcular_metricas_red_total(addres_dir,filexdir,thsrho,min_nNodos,xnz,nNodos):
    
    '''
    retorna:numero de nodos, 
            intervalo de tiempo
            densidad de red
            grado medio de red
            propension mundo pequeño
            longitud de camino
            clustering
    '''
    #fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))  
    mxrango = slice(14,22)
    sub_cadena=[]; dtime=[]; mxnkover=[]; nNum_total=[]; densidad_de_red=[];modul=[];wpesos=[];dst_mean=[]
    num_conex_exist=[]; num_conex_posibles=[]
    nNfiles=len(addres_dir); mxswp=[]; mxlongpath=[]; mxclustering=[]; nNum_desconect=[]
    error_xclust=[];error_xpl=[];error_grado_medio=[];comp_gig=[]; xlp=[];std_lp=[]
    for i in range(nNfiles):
        sub_cadena.append(int(addres_dir[i][mxrango].split('.')[0]))
    
    file_dir= filexdir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    
    sub_cadena.sort()
    for i in range(nNfiles):

        xmc=np.array(abs((np.genfromtxt(file_dir+str(sub_cadena[i])+'.txt'))))
        #zzzGr = nx.from_numpy_array(xmc)
        dtx=np.array(abs(delete_nod(np.genfromtxt(file_dir+str(sub_cadena[i])+'.txt'))))
        xGr = nx.from_numpy_array(dtx)
        if((len(xGr.nodes)>min_nNodos)):# and (dzx >= thsrho)):
            nN=len(xGr.nodes)
            #nNd=len(zzzGr.nodes)
            #ncp=nN*(nN-1)/2
            cg=nN/nNodos
            comp_gig.append(cg)
            #print(file_dir,sub_cadena[i],"\t",cg,"\t",nN,"\t",nNodos)
            nNum_total.append(nN)
            nNum_desconect.append(nNodos-nN)
            print(nNodos)
            # vector_k_lover= mean_degree_network_tot(dtx)
            # mean_degree=sum(vector_k_lover)/nNodos
            # std_degree=np.std(vector_k_lover)/np.sqrt(len(vector_k_lover))
            # error_grado_medio.append(std_degree)
            # mxnkover.append(mean_degree)
            # dtime.append(sub_cadena[i])
            #long_path, error_lp = path_lenght(xmc)
        

    return dtime,comp_gig,mxnkover,error_grado_medio

def calcular_metricas(addres_dir,filexdir,thsrho,min_nNodos,xnz):
    '''
    retorna:numero de nodos, 
            intervalo de tiempo
            densidad de red
            grado medio de red
            propension mundo pequeño
            longitud de camino
            clustering
    '''
    #fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))  
    mxrango = slice(14,22)
    sub_cadena=[]; dtime=[]; mxnkover=[]; nNum_total=[]; densidad_de_red=[];modul=[];wpesos=[];dst_mean=[]
    num_conex_exist=[]; num_conex_posibles=[]
    nNfiles=len(addres_dir); mxswp=[]; mxlongpath=[]; mxclustering=[]; nNum_desconect=[]
    error_xclust=[];error_xpl=[];error_grado_medio=[];fraccion_nod=[]; xlp=[];std_lp=[]
    for i in range(nNfiles):
        sub_cadena.append(int(addres_dir[i][mxrango].split('.')[0]))
    
    file_dir= filexdir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]

    sub_cadena.sort()
    for i in range(nNfiles):

        xmc=np.array(abs((np.genfromtxt(file_dir+str(sub_cadena[i])+'.txt'))))
        zzzGr = nx.from_numpy_array(xmc)
        dtx=np.array(abs(delete_nod(np.genfromtxt(file_dir+str(sub_cadena[i])+'.txt'))))
        # datx=np.zeros((len(dtx[:,0]),len(dtx[:,0])))
        # io,jo=np.where(dtx!=0)
        # datx[io,jo]=1
        xGr = nx.from_numpy_array(dtx)
        #fxGr = nx.from_numpy_array(fxmc)
        dzx=nx.density(xGr)
        #if(len(fxGr.nodes)>min_nNodos):
        if((len(xGr.nodes)>min_nNodos)):# and (dzx >= thsrho)):
            #xmc=arbole_expansion_min(datx,thsrho,xnz)
            nN=len(xGr.nodes)
            nNd=len(zzzGr.nodes)
            ncp=nN*(nN-1)/2
            nodefrac=nN/nNd
            fraccion_nod.append(nodefrac)
            print(file_dir,sub_cadena[i])
            num_conex_exist.append(np.count_nonzero(np.triu(dtx).flatten()))
            nNum_total.append(nN)
            nNum_desconect.append(nNd)
            num_conex_posibles.append(ncp)
            densidad_de_red.append(nx.density(nx.from_numpy_array(dtx)))
            vector_k_lover, nkover, error_degree= mean_degree_network(dtx)
            error_grado_medio.append(error_degree)
            mxnkover.append(nkover)
            swp=netSWP(dtx)
            long_path, error_lp, mean_cl, error_mean_cl = metricas_L_C(dtx)
            lg_path,std_path_leght=path_lenght(dtx)
            mxswp.append(swp);mxlongpath.append(long_path);mxclustering.append(mean_cl)
            error_xpl.append(error_lp);error_xclust.append(error_mean_cl)
            #mxswp=0;mxlongpath=0;mxclustering=0
            dtime.append(sub_cadena[i])
            modul.append(modularidad(dtx))
            wpesos.append(np.mean(dtx))
            dst_mean.append(np.std(dtx)/np.sqrt(len(dtx.flatten())))
            xlp.append(lg_path)
            std_lp.append(std_path_leght)
            
            #arbol_exp_min.append(arbole_expansion_min(xm))
    return dtime,nNum_total,nNum_desconect,fraccion_nod, densidad_de_red,mxnkover,error_grado_medio,num_conex_exist, num_conex_posibles,mxlongpath,error_xpl,mxclustering,error_xclust,mxswp,modul,wpesos,dst_mean,xlp,std_lp

def normalizar_metricas(addres_dir,filexdir,xthsrho,min_nNodos,xnz=1):
    
    nNfiles=len(addres_dir); mxrango = slice(14,22)
    sub_cadena=[];

    sub_cadena=[]; dtime=[]; nNum=[]; densidad_de_red=[];modul=[]
    num_conex_exist=[]; num_conex_posibles=[]; mxnkover=[]; nod_tot=[]

    for i in range(nNfiles):
        sub_cadena.append(int(addres_dir[i][mxrango].split('.')[0]))
    
    file_dir= filexdir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]

    sub_cadena.sort()

    LP_norm=[]; CL_norm=[];densidad=[]

    dtime=[]
    for i in range(nNfiles):
        datx=np.array(abs(delete_nod(np.genfromtxt(file_dir+str(sub_cadena[i])+'.txt'))))
        xzm=np.array(abs(np.genfromtxt(file_dir+str(sub_cadena[i])+'.txt')))
        xGr = nx.from_numpy_array(datx)
        zG = nx.from_numpy_array(xzm)
        num_nodos=len(xGr.nodes)
        mod_nul_LP=[]; mod_nul_CL=[]
        dzx=nx.density(xGr)

        if(num_nodos>min_nNodos and (dzx >= xthsrho)):
            xmc=arbole_expansion_min(datx,xthsrho,xnz)
            xG = nx.from_numpy_array(xmc)
            densidad.append(nx.density(xG))
            # print(rt , dzx , xthsrho)
            vectk_over,k_over_2= mean_degree_network(xmc)
            print(file_dir,sub_cadena[i])
            xlp,xcl=metricas_L_C(xmc)
            zrx=len(zG.nodes)
            nNum.append(num_nodos)
            mxnkover.append(k_over_2)
            nod_tot.append(zrx)
            nun_conec=len(xG.edges)
            dtime.append(sub_cadena[i])
            # lp_ramd,cl_ramd=metricas_L_C(null_ramdomize)
            # lp_latt,cl_latt=metricas_L_C(null_lattice)
            # xLp=(xlp-lp_ramd)/(lp_latt-lp_ramd)
            # xCl=(cl-cl_ramd)/(cl_latt-cl_ramd)
            # if(xLp<1 and xCl<1):
            #     if(xLp>0 and xCl>0)
            #null_lattice,null_ramdomize=random_lattice(xmc)

            for k in range(5):
                
                nulo_bin=generate_small_world_networks(num_nodos,int(k_over_2),0.025,int(nun_conec))              
                modelo_nulo=construir_enlaces(xmc,nulo_bin)
                xLp,xCl=metricas_L_C(modelo_nulo)
                mod_nul_LP.append(xLp)
                mod_nul_CL.append(xCl)

            media_LP=np.mean(mod_nul_LP); fpl=xlp/media_LP
            media_CL=np.mean(mod_nul_CL); fcl=xcl/media_CL
            LP_norm.append(fpl)
            CL_norm.append(fcl)  

    return dtime, LP_norm, CL_norm, mxnkover,nNum,nod_tot,densidad

def xplot_metricas(dataExp,xsubC,numCol,ancho_bin_hitograma):
    '''
    fraccion de nodos en funcion del 
    tiempo color en el cv
    '''
    
    fig = plt.figure(figsize=(7, 5))
    gs = fig.add_gridspec(3,1)

    delt=[30,40,50]
    
    for i,xsubdir in enumerate(dataExp):
        
        #fdat=np.loadtxt("datMetxUnc/"+str(xsubdir)+"_densid.txt")#datos empiricos
        fdat=np.loadtxt("datMetx/"+str(xsubdir)+".txt")#datos empiricos
        
        frac_nod=(np.array(fdat[:,numCol])-np.array(fdat[:,numCol+1]))/np.array(fdat[:,numCol])
        frec, xbin=np.histogram(frac_nod,bins=30)
        xo=xbin[:-1]; yo= frec/sum(frec)

        viridis = mpl.colormaps['plasma'].resampled(len(xo))
        newcolors = viridis(np.linspace(0,1,len(xo)))
        newcmp = ListedColormap(newcolors)


        custom_params = {"axes.spines.right": False, "axes.spines.top": False}       
        sns.set_theme(font_scale=1.2,style="white",rc=custom_params)
        
        axhx = fig.add_subplot(gs[i,:4])


        miscorlores=["k","b","r"]
        dfx = {'xx': delt[i]*fdat[:,0], 'yy': fdat[:,1]}
        df = {'x':delt[i]*fdat[:,0], 'y': frac_nod}
        
        #sns.scatterplot(x = "x", y = "y", data = df,marker = "o",color=miscorlores[i],ax=axhx,label=r"Experimental data "+str(xsubC)+", $\Delta{t}=$"+str(delt[i])+" s")
        sns.lineplot(x = "x", y = "y", data = df,marker = "o",color="black",ax=axhx,label=r"CV ")
        sns.lineplot(x = "xx", y = "yy", data = dfx,marker = "o",color="blue",ax=axhx,label=r"node Fraction ")
        axhx.legend(title="Exp. "+str(xsubC)+", $\Delta{t}=$"+str(delt[i])+" s")
        label_format = r'${:.2f}$'
        label_format_ylist = r'${:.1f}$'
        #ticks_loc =cbar.ax.get_yticks().tolist()

        ticks_loc=np.linspace(min(fdat[:,1]),max(fdat[:,1]),5)
        xlist=np.linspace(min(yo),max(yo),4)
        ylist=np.linspace(0,max(fdat[:,1]),5)

        axhx.yaxis.set_major_locator(mticker.FixedLocator(ylist))
        axhx.set_yticklabels([label_format_ylist.format(x) for x in ylist],fontsize=14)
        
        axhx.set_ylabel(r"",fontsize=15)
        axhx.yaxis.set_tick_params(labelsize=15)

        axhx.set_xlabel(r"Time(s)",fontsize=15)
        axhx.xaxis.set_tick_params(labelsize=15)
    
    plt.tight_layout()
    plt.show()
    
def plot_metricas_red(dataExp,xsubC,numCol,ancho_bin_hitograma):
    '''
    fraccion de nodos en funcion del 
    tiempo color en el cv
    '''
    
    fig = plt.figure(figsize=(12, 8))
    #gs = fig.add_gridspec(3, 5)
    gs = fig.add_gridspec(2, 3)
    # axhx0 = fig.add_subplot(gs[0,:4])
    # axhy0 = fig.add_subplot(gs[0,4])

    mibar=[[0.90, 0.704, 0.01, 0.228],
           [0.90, 0.401, 0.01, 0.228],
           [0.90, 0.1, 0.01, 0.228]]

    delt=[30,40,50]
    for i,xsubdir in enumerate(dataExp):
        
        #fdat=np.loadtxt("datMetxUnc/"+str(xsubdir)+"_densid.txt")#datos empiricos
        fdat=np.loadtxt("restW10T50/"+str(xsubdir)+".txt")#datos empiricos
        
        frac_nod=(np.array(fdat[:,numCol]))#-np.array(fdat[:,numCol+1]))/np.array(fdat[:,numCol])
        frec, xbin=np.histogram(frac_nod,bins=30)
        xo=xbin[:-1]; yo= frec/sum(frec)

        viridis = mpl.colormaps['plasma'].resampled(len(xo))
        newcolors = viridis(np.linspace(0,1,len(xo)))
        newcmp = ListedColormap(newcolors)

        # axhx1 = fig.add_subplot(gs[1,:4])
        # axhy1 = fig.add_subplot(gs[1, 4])
        # axhx2 = fig.add_subplot(gs[2,:4])
        # axhy2 = fig.add_subplot(gs[2, 4])

        custom_params = {"axes.spines.right": False, "axes.spines.top": False}       
        sns.set_theme(font_scale=1.2,style="ticks",rc=custom_params)
        
        ax0 = fig.add_subplot(gs[0,0])
        ax1 = fig.add_subplot(gs[0,1])
        ax2 = fig.add_subplot(gs[0,2])
        ax3 = fig.add_subplot(gs[1,0])
        ax4 = fig.add_subplot(gs[1,1])
        ax5 = fig.add_subplot(gs[1,2])
        #axhy = fig.add_subplot(gs[i,4])

        miscorlores=["k","b","r"]
        df = {'t': 30*fdat[:,0], 'CV': fdat[:,1]}

        x0df = {'t': 30*fdat[:,0], 'Grade': fdat[:,5]}
        x1df = {'t': 30*fdat[:,0], 'Path length': fdat[:,8]}
        x2df = {'t': 30*fdat[:,0], 'clustering': fdat[:,9]}
        x3df = {'t': 30*fdat[:,0], 'Nodes': fdat[:,2]}
        x4df = {'t': 30*fdat[:,0], 'Density': fdat[:,4]}

        #sns.scatterplot(x = "x", y = "y", data = df,marker = "o",color=miscorlores[i],ax=axhx,label=r"Experimental data "+str(xsubC)+", $\Delta{t}=$"+str(delt[i])+" s")
        sns.lineplot(x = "t", y = "CV", data = df,marker = "o",color="red",ax=ax0)#r"Experimental data "+str(xsubC)+", $\Delta{t}=$"+str(delt[i])+" s")
        sns.lineplot(x = "t", y = "Grade", data = x0df,marker = "o",color=miscorlores[i],ax=ax1)#r"Experimental data "+str(xsubC)+", $\Delta{t}=$"+str(delt[i])+" s")
        sns.lineplot(x = "t", y = "Path length", data = x1df,marker = "o",color=miscorlores[i],ax=ax2)#r"Experimental data "+str(xsubC)+", $\Delta{t}=$"+str(delt[i])+" s")
        sns.lineplot(x = "t", y = "clustering", data = x2df,marker = "o",color=miscorlores[i],ax=ax3)#r"Experimental data "+str(xsubC)+", $\Delta{t}=$"+str(delt[i])+" s")
        sns.lineplot(x = "t", y = "Nodes", data = x3df,marker = "o",color=miscorlores[i],ax=ax4)
        sns.lineplot(x = "t", y = "Density", data = x4df,marker = "o",color=miscorlores[i],ax=ax5)
        '''
        #for row in range(len(xo)):
        #    axhy.barh(xo[row],yo[row],height=ancho_bin_hitograma[i],ec='w',color=viridis.colors[len(xo)-row-1])
        axhy.barh(xo,yo,height=ancho_bin_hitograma[i],ec='w')
        axhy.tick_params(axis="y", labelleft=False)
        # axhy.invert_yaxis()
        # if(i<2):
        #     axhx.tick_params(axis='x', labelbottom=False) 
        #     axhy.tick_params(axis='x', labelbottom=False)  

    
        # cb=plt.cm.ScalarMappable(norm=Normalize(min(fdat[:,1]),max(fdat[:,1]) ), cmap=newcmp)
        # cax = plt.axes(mibar[i])
        # cbar=fig.colorbar(cb, cax=cax,ax=axhy, label="")
        label_format = r'${:.2f}$'
        label_format_ylist = r'${:.1f}$'
        #ticks_loc =cbar.ax.get_yticks().tolist()

        ticks_loc=np.linspace(min(fdat[:,1]),max(fdat[:,1]),5)
        xlist=np.linspace(min(yo),max(yo),4)
        ylist=np.linspace(min(frac_nod),max(frac_nod),5)

        # cbar.ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        # cbar.ax.set_yticklabels([label_format.format(x) for x in ticks_loc],fontsize=14)
        # cbar.ax.set_ylabel(r"$CV$",fontsize=15)
        
        axhx.yaxis.set_major_locator(mticker.FixedLocator(ylist))
        axhx.set_yticklabels([label_format_ylist.format(x) for x in ylist],fontsize=14)
        
        axhx.set_ylabel(r"Average grade",fontsize=15)
        axhx.yaxis.set_tick_params(labelsize=15)
        
        axhy.xaxis.set_major_locator(mticker.FixedLocator(xlist))
        axhy.set_xticklabels([label_format.format(x) for x in xlist],fontsize=14)     

        axhy.xaxis.set_tick_params(labelsize=15)
        axhx.set_xlabel(r"Time(s)",fontsize=15)
        axhx.xaxis.set_tick_params(labelsize=15)
        axhy.xaxis.set_tick_params(labelsize=15)
        axhy.set_xlabel(r"Density",fontsize=15)
        '''
    plt.tight_layout()
    plt.show()

def xcolors(szM):
    xo=np.loadtxt("datMetx/"+str(szM)+".txt")#datos empiricos
    viridis = mpl.colormaps['plasma'].resampled(len(xo[:,1]))
    newcolors = viridis(np.linspace(0,1,len(xo[:,1])))
    newcmp = ListedColormap(newcolors)
    return viridis

def plot_metricas_compare(xsubC,numCol):
    '''
    fraccion de nodos en funcion del 
    tiempo color en el cv
    '''
    xfig = plt.figure()
    ax = xfig.add_subplot(projection='3d')
    #["Mar10","Mar04","Jan29","Jan21","Jan14","Dez20"]
    for namExp in ["Mar07","Mar10","Mar04","Jan29","Jan21","Jan14"]:
        fig = plt.figure(figsize=(15, 5.5))
        gs = fig.add_gridspec(1, 3)
  

        xlim=[4,4,4]; max_cv=[];  max_net=[]
        
        custom_params = {"axes.spines.right": False, "axes.spines.top": False}       
        sns.set_theme(font_scale=1.2,style="ticks",rc=custom_params)


        fig.suptitle(r'Variation of network metrics on experimental data '+namExp+'.',size="x-large")
        #fig.suptitle(r"Variation of network metrics on experimental data",size="x-large")

        for o in range(3): 
            xn=250
            dataExp=[
                    str(xn)+"s"+namExp+"",
                    # str(xn)+"sMar04",
                    # str(xn)+"sMar10",
                    # str(xn)+"sJan29",
                    # str(xn)+"sMar07",
                    
                
                    
                    #str(xn)+"sJan14",
                    #str(xn)+"sJan21",
                    #str(xn)+"sDez20"
                    ]

            mibar=[[0.90, 0.704, 0.01, 0.228],
                [0.90, 0.401, 0.01, 0.228],
                [0.90, 0.1, 0.01, 0.228]]

            delt=7*[30]
            
            axhx = fig.add_subplot(gs[0,o])
            
            for i,xsubdir in enumerate(dataExp):
                #viridis=xcolors(xsubdir)
                #fdat=np.loadtxt("Exp_wG_densy/"+str(xsubdir)+".txt")#datos empiricos
                #fdat=np.loadtxt("metric_norm/"+str(xsubdir)+".txt")#datos empiricos
                fdat=np.loadtxt("restW10T250/"+str(xsubdir)+".txt")#datos empiricos
                #fcv=np.loadtxt("retratostxt/"+str(xsubdir)+"retr.txt")#datos empiricos
            
                # divergr=[]; dcv=[];divJS=[]

                # for fil in range(len(fdat[:,1])):
                #     divergr.append(np.mean(fcv[fil,:]))#-np.array(fdat[:,numCol+1]))/np.array(fdat[:,numCol])
                    
                # for xfil in range(len(fdat[:,1])):
                #     if(divergr[xfil]!=0):
                #         dcv.append(fdat[xfil,1])
                #         divJS.append(divergr[xfil])
                # frec, xbin=np.histogram(fdat[:,10],bins=30)
                # xo=xbin[:-1]; yo= frec/sum(frec)
                # plt.bar(xo,yo, width=0.02,alpha=1,ec='w')

                zn=np.array([[1,14],[1,14],[1,14]])
                
                xdf = {'x': (fdat[:,zn[o,0]]), 'y': (fdat[:,zn[o,1]])}
                xdrf = {'x': (fdat[:,zn[o,0]]), 'y': fdat[:,zn[o,1]]}

                # # sns.displot(df, x="x", y="y",kind="kde")

                # #sns.kdeplot(data=df, x="y", fill=True, common_norm=False, alpha=.5, linewidth=0)
                # # plt.xscale("log")
                micolo=["blue","red","black","gray","green","purple","orange"]

                #sns.lineplot(data = xdf,x = "x", y = "y",legend=False)#,label=r"Experimental data "+str(xsubC[i]))#+", $\Delta{t}=$"+str(delt[i])+" s")
                miscorlores=["m","b","r"]
                sns.lineplot(data = xdf,x = "x", y = "y",marker = "o",color=miscorlores[o],ax=axhx,legend=False,markersize=9,linewidth = 3)#,label=r"Experimental data "+str(xsubC[i]))#+", $\Delta{t}=$"+str(delt[i])+" s")
                #sns.lineplot(data = xdrf,x = "x", y = "y",marker = "o",color="black",ax=axhx,legend=False,markersize=9,linewidth = 3)#,label=r"Experimental data "+str(xsubC[i]))#+", $\Delta{t}=$"+str(delt[i])+" s")

                #for fil in range(len(fdat[:,1])):
                #axhx.plot((fdat[:,zn[o,0]]), (fdat[:,zn[o,1]]),"o",markersize=3,color=micolo[i])#,color=viridis.colors[fil])
                # axhx.plot((fdat[:,1]), (fdat[:,5]),"o",markersize=3,color=micolo[i+1])#,color=viridis.colors[fil])
                #axhx.plot((xfdat[:,5]), (xfdat[:,8]),"o",markersize=3,color=micolo[i+1])#,color=viridis.colors[fil])            
                # plt.yscale("log")
                #ax.scatter(fdat[:,1], fdat[:,4], fdat[:,8])


                # for fil in range(len(fdat[:,1])):
                    #ax.plot((fdat[fil,1]), (fdat[fil,1]), fdat[fil,2],"o",markersize=2,color=viridis.colors[fil])
                #ax.plot(fdat[:,1], fdat[:,5], fdat[:,9],"o",markersize=2,color=miscorlores[o])

                #axhx.plot(dcv,divJS,"o", markersize=2,alpha=0.8)
                #axhx.plot(fdat[:,1],s_w_p)
                # max_net.append(max(divJS))
                # max_cv.append(max(dcv))
            # ax.set_xlabel(r"Nodes",fontsize=20)
            # ax.set_ylabel(r"Grade",fontsize=20)
            # ax.set_zlabel(r"Average clustering",fontsize=20)
            
            # ax.yaxis.set_tick_params(labelsize=17)
            # ax.xaxis.set_tick_params(labelsize=17)
            # ax.zaxis.set_tick_params(labelsize=17)

            label_format = r'${:.1f}$'
            label_format_ylist = r'${:.2f}$'
            if(o==1):
                label_format_ylist = r'${:.1f}$'

            # xlist=np.linspace(min(fdat[:,zn[o,0]]),max(fdat[:,zn[o,0]]),7)
            # ylist=np.linspace(min(fdat[:,zn[o,1]]),max(fdat[:,zn[o,1]]),6)
            
            # axhx.xaxis.set_major_locator(mticker.FixedLocator(xlist))
            # axhx.set_xticklabels([label_format.format(x) for x in xlist],fontsize=20)

            # axhx.yaxis.set_major_locator(mticker.FixedLocator(ylist))
            # axhx.set_yticklabels([label_format_ylist.format(x) for x in ylist],fontsize=20)
            
            #rxm=np.array([[r"$\langle L \rangle$",r"$\langle k \rangle$"],[r"$\langle k \rangle$",r"$\langle CV \rangle$"]])
            rxm=np.array([[r"$k$",r"$\langle CV \rangle$"],
                        [r"$\langle L \rangle$",r"$CV$"],
                        [r"$\langle C \rangle$",r"$CV$"]])
            
            axhx.set_ylabel(r""+str(rxm[o,0])+"",fontsize=20)
            axhx.yaxis.set_tick_params(labelsize=20)

            axhx.set_xlabel(r""+str(rxm[o,1])+"",fontsize=20)
            axhx.xaxis.set_tick_params(labelsize=20)

            epsilon=np.array([[0.1,0.1],[0.1,0.1],[0.1,0.1]])
            axhx.set_xlim(min(fdat[:,zn[o,0]])-epsilon[o,0],max(fdat[:,zn[o,0]])+epsilon[o,0])
            axhx.set_ylim(min(fdat[:,zn[o,1]])-epsilon[o,1],max(fdat[:,zn[o,1]])+epsilon[o,1])

            zmx=["Experimental data Mar. 04","Experimental data Mar. 04"]
            font = {'family': 'serif',
                    'color':  'black',
                    'weight': 'normal',
                    'size': 15,
                    }
            # plt.title(r""+str(zmx[o])+"",fontdict=font)
        plt.tight_layout()
        plt.show()

###########################################################################
###########################################################################
###########################################################################

matplotlib.rcParams['text.usetex'] = True    

#fig.suptitle(" ", size=14)
#fig.subplots_adjust(wspace=0.5)
interv=np.arange(10,16,1)

adrespr =["Mar07","Mar10","Mar04","Jan29","Jan21","Jan14"]#os.listdir('mcporcv')
nome_exper=["SUAM07","SUAM10","SUAM04","SUAJ29","SUAJ21","SUAJ14"]
bin_cv=0.050; subxc="Mar. 07"
cv_time=10

num_Col=1
bin_hist=[0.1,0.1,0.1]
dataExp=[152,144,73,76,88,76]
# for o,xn in enumerate([50]): 
#     dataExp.append(str(xn)+"sMar07")
#plot_metricas_red(dataExp,subxc,num_Col,bin_hist)
#xplot_metricas(dataExp,subxc,num_Col,bin_hist)
addir = os.listdir('datMetx')

#plot_metricas_compare(subxc,num_Col)
nNods=[]
fig, zx = plt.subplots(2, 6, figsize=(17,10 ))
difile='datW10T250s/'
rest_file="restW10T250/"
janela=250
def metrik():
    x_nz=4
    xcv=[]
    for l,ndir in enumerate(adrespr):
        #xdr=os.listdir("datMetx/"+ndir)
        data=np.loadtxt('datos_cv/datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
        for xsubdir in [str(janela)+'s'+str(ndir)]:        
            addres_dir=os.listdir(difile+adrespr[l]+"/"+str(xsubdir))        
            dir_file=difile+adrespr[l]+"/"+str(xsubdir)
            fdat=np.loadtxt(rest_file+str(xsubdir)+".txt")#datos empiricos
            print(ndir)
            a_x = zx[0,l]
            a_y = zx[1,l]
            pgraf=[a_x,a_y]
            
           
            thrho=0#min(fdat[:,4])#-2*np.std(fdat[:,4])
            # print(cv_time," ",xsubdir)
            #dT,nNodos_tot,nNodos_desconet,frac_nod, densidad_de_red,nkover,err_nkover,nNoconex_exist,nNum_posib,long_path,error_long_path,xcluster,error_cluster,swpx,modul,pesos,xdst_mean,phtL,xdst_Phtl = calcular_metricas(addres_dir,dir_file,thrho,10,x_nz)        
            dT,frac_nod, nkover,err_nkover = calcular_metricas_red_total(addres_dir,dir_file,thrho,10,x_nz,dataExp[l])        
           
            #dT,pl_norm,cl_norm,xgrad,knods,no_desct,densid=normalizar_metricas(addres_dir,dir_file,thrho,10,x_nz)
'''  
            #int(xsubdir.split('s')[0])
            #calcular_metricas(addres_dir,dir_file,thrho,10,x_nz)
            #promedio_error_metricas(addres_dir,dir_file)
            #power_law(addres_dir,dir_file,10)
            #fdat=np.loadtxt(rest_file+str(xsubdir)+".txt")#datos empiricos

            #plot_metricas(l,pgraf,fdat[:,0],fdat[:,1],dir_file,addres_dir,nome_exper[l])
            #ajuste_distribucion_metricas(fdat[:,0],fdat[:,1],dir_file,addres_dir,nome_exper[l])
            #print(np.mean(fdat[:,3])," ",np.std(fdat[:,3])," ",np.mean(fdat[:,3])-2*np.std(fdat[:,3]))
            
            # frec, xbin=np.histogram(fdat[:,3],bins=50)
            # xo=xbin[:-1]; yo= frec/sum(frec)
            # plt.bar(xo,yo, width=1,alpha=1,ec='w')

            #plt.tight_layout()
            #fig.savefig('fhc.jpg')
            #plt.show()
            

            if(ndir=="Mar04" or ndir=="Mar07"):
                tcb2=NAT2(data[0,:])
            else:
                tcb2=NAT2(data[:,0])       
                    
            CV=tcb2.CV_calc(cv_time,bin_cv, show=False, save=False, ADFA=False, 
                                    showbin=False, bin_size=None)
            

            cv_dat=[]; xmx=tcb2.CV; 
            
            for tcv in dT: 
                cv_dat.append(xmx[tcv])

            xmx_sort=np.sort(xmx)
            cv_sort=xmx_sort[xmx_sort<3]
            average_CV=[]
            for cvi in range(len(dT)):
                a=25*(cvi)
                b=25*(cvi+1)
                avercv=[]
                for tk in range(a,b):
                    avercv.append(cv_sort[tk])
                average_CV.append(np.mean(avercv))
            
            file_txt=rest_file+xsubdir+"Tot.txt"
            fdata = open(file_txt, 'w')
        
            try:
                for i in range(len(dT)):

                    fdata.write(str(dT[i])
                                +"\t"+str(average_CV[i])
                                +"\t"+str(frac_nod[i])
                                +"\t"+str(nkover[i])
                                +"\t"+str(err_nkover[i])
                                +"\n")                        
            finally:
                fdata.close()

'''
'''            
            file_txt=rest_file+xsubdir+".txt"
            #file_txt="datMetx/"+xsubdir+"_fracion.txt"
            fdata = open(file_txt, 'w')
            
            #swp=np.array(swpx)
        
            try:
                #for k in range(len(cv_dat)):
                for i in range(len(dT)):#len(cv_dat)):   
                    #if(xmx_sort[k]==cv_dat[i] and (cv_dat[i]<3)):
                    fdata.write(str(dT[i])
                                +"\t"+str(average_CV[i])
                                +"\t"+str(nNodos_tot[i])
                                +"\t"+str(nNodos_desconet[i])
                                +"\t"+str(frac_nod[i])
                                +"\t"+str(densidad_de_red[i])
                                +"\t"+str(nkover[i])
                                +"\t"+str(err_nkover[i])
                                +"\t"+str(nNoconex_exist[i])
                                +"\t"+str(nNum_posib[i])
                                +"\t"+str(long_path[i])
                                +"\t"+str(error_long_path[i])
                                +"\t"+str(xcluster[i])
                                +"\t"+str(error_cluster[i])
                                +"\t"+str(swp[i][0])
                                +"\t"+str(modul[i])
                                +"\t"+str(pesos[i])
                                +"\t"+str(xdst_mean[i])
                                +"\t"+str(phtL[i])
                                +"\t"+str(xdst_Phtl[i])
                                +"\n")                        
                        #break
            finally:
                fdata.close()
'''
'''
            
            # file_txt="metric_norm/"+xsubdir+".txt"
            # fdata = open(file_txt, 'w')

            # try:
            #     #for k in range(len(cv_dat)):
            #         for i in range(len(dT)):   
            #             #if(xmx_sort[k]==cv_dat[i] and (cv_dat[i]<5)):
            #             fdata.write(str(dT[i])+"\t"+str(average_CV[i])+"\t"+str(pl_norm[i])
            #                         +"\t"+str(cl_norm[i])+"\t"+str(xgrad[i])+"\t"+str(knods[i])
            #                         +"\t"+str(no_desct[i])+"\t"+str(densid[i])+"\n") 
            #             #break
            # finally:
            #    fdata.close()


for m in addir:#["30sMar07.txt"]:#addir:
    print(m)
    fig = plt.figure(tight_layout=True)
    gs = gridspec.GridSpec(2, 3)
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[0,1])
    ax2 = fig.add_subplot(gs[0,2])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    ax5 = fig.add_subplot(gs[1,2])
    
    f = open("datMetx/"+str(m), 'r')
    contenido = f.read()
    
    if(contenido!=''):
        dt,cv,nodos,densid,gradmedio,swp,lp,clus,arb_exp_min,modulard =np.loadtxt("datMetx/"+str(m),skiprows=0, usecols=[0,1,2,3,4,5,6,7,8,9],unpack=True)#datos empiricos

        tips=np.loadtxt("datMetx/"+str(m))
        #tips=np.loadtxt("datMetx/10sMar07Ret.txt")
        
        # norm = matplotlib.colors.Normalize(vmin=min(cv), vmax=max(cv)) 
        
        # sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
        yxc=tips[:,6]
        sns.scatterplot(data=tips, x=tips[:,1], y=yxc, hue=tips[:,1],palette=cmap,legend=False,ax=ax0)
        sns.scatterplot(data=tips, x=tips[:,1], y=tips[:,7], hue=tips[:,1],palette=cmap,legend=False,ax=ax1)
        sns.scatterplot(data=tips, x=tips[:,1], y=tips[:,4], hue=tips[:,1],palette=cmap,legend=False,ax=ax2)
        sns.scatterplot(data=tips, x=tips[:,1], y=tips[:,5], hue=tips[:,1],palette=cmap,legend=False,ax=ax3)
        sns.scatterplot(data=tips, x=tips[:,1], y=tips[:,1], hue=tips[:,1],palette=cmap,legend=False,ax=ax4)
        sns.scatterplot(data=tips, x=tips[:,1], y=tips[:,9], hue=tips[:,1],palette=cmap,legend=False,ax=ax5)


        ax0.xaxis.set_tick_params(labelsize=13)
        ax0.yaxis.set_tick_params(labelsize=13)

        ax0.set_xlabel(r"$CV$",fontsize=15)
        ax0.set_ylabel(r"Average path lengt",fontsize=15)

        ax0.xaxis.set_tick_params(labelsize=13)
        ax0.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax0.set_xlabel(r"$CV$",fontsize=15)
        ax0.set_ylabel(r"Average path lengt",fontsize=15)

        ax1.xaxis.set_tick_params(labelsize=13)
        ax1.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax1.set_xlabel(r"$tiempo$",fontsize=15)
        ax1.set_ylabel(r"densidad",fontsize=15)

        ax2.xaxis.set_tick_params(labelsize=13)
        ax2.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax2.set_xlabel(r"$CV$",fontsize=15)
        ax2.set_ylabel(r"Average path lengt",fontsize=15)


        ax3.xaxis.set_tick_params(labelsize=13)
        ax3.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax3.set_xlabel(r"$CV$",fontsize=15)
        ax3.set_ylabel(r"Average path lengt",fontsize=15)


        ax4.xaxis.set_tick_params(labelsize=13)
        ax4.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax4.set_xlabel(r"$CV$",fontsize=15)
        ax4.set_ylabel(r"Average path lengt",fontsize=15)


        ax5.xaxis.set_tick_params(labelsize=13)
        ax5.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax5.set_xlabel(r"$CV$",fontsize=15)
        ax5.set_ylabel(r"Average path lengt",fontsize=15)
        fig.tight_layout()
        plt.show()
     

for i in range(len(cv[cv<5])):
    if(swp[i]>0):
        ax0.plot(cv[i],clus[i],"o",linewidth=1,color=cmap(i))      
        ax1.plot(cv[i],lp[i],"o",linewidth=1,color=cmap(i))
        ax2.plot(cv[i],densid[i],"o",linewidth=1,color=cmap(i))
        ax3.plot(cv[i],gradmedio[i],"o",linewidth=1,color=cmap(i))
        ax4.plot(cv[i],nodos[i],"o",linewidth=1,color=cmap(i))
        ax5.plot(cv[i],swp[i],"o",linewidth=1,color=cmap(i))

        #axx = plt.axes(projection='3d')
        #for i in range(len(cv)):
        #axx.scatter3D(nodos[i],gradmedio[i], lp[i],color=cmap(i))#, c=zdata, cmap='Greens');
        #ax.plot(x,y,linewidth=1)#,color=cmap(i))

        # axx.set_xlabel('Nodos')
        # axx.set_ylabel('Gradio medio')
        # axx.set_zlabel('clustering')

        #ax.set_yscale('log')
        #ax.set_xscale('log')

        ax0.xaxis.set_tick_params(labelsize=13)
        ax0.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax0.set_xlabel(r"$CV$",fontsize=15)
        ax0.set_ylabel(r"Average path lengt",fontsize=15)

        ax1.xaxis.set_tick_params(labelsize=13)
        ax1.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax1.set_xlabel(r"$CV$",fontsize=15)
        ax1.set_ylabel(r"Average path lengt",fontsize=15)

        ax2.xaxis.set_tick_params(labelsize=13)
        ax2.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax2.set_xlabel(r"$CV$",fontsize=15)
        ax2.set_ylabel(r"Average path lengt",fontsize=15)


        ax3.xaxis.set_tick_params(labelsize=13)
        ax3.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax3.set_xlabel(r"$CV$",fontsize=15)
        ax3.set_ylabel(r"Average path lengt",fontsize=15)


        ax4.xaxis.set_tick_params(labelsize=13)
        ax4.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax4.set_xlabel(r"$CV$",fontsize=15)
        ax4.set_ylabel(r"Average path lengt",fontsize=15)


        ax5.xaxis.set_tick_params(labelsize=13)
        ax5.yaxis.set_tick_params(labelsize=13)
        #cbar.ax.set_ylabel(r"$CV$",fontsize=15)

        ax5.set_xlabel(r"$CV$",fontsize=15)
        ax5.set_ylabel(r"Average path lengt",fontsize=15)
        
divider = make_axes_locatable(ax0)

cax = divider.append_axes("right" ,size="4%", pad=0.04)

cbar=fig.colorbar(sm, cax=cax)

label_format = r'${:.1f}$'
ticks_loc = cbar.ax.get_yticks().tolist()
cbar.ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
cbar.ax.set_yticklabels([label_format.format(x) for x in ticks_loc],fontsize=13)
'''
        
metrik()