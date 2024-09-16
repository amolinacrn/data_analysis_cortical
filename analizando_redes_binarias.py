from cProfile import label
import re
from tkinter import Listbox
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True 
from email import contentmanager
import scipy.io
import os
import h5py
from smallworld import get_smallworld_graph
from bct.algorithms import number_of_components,randmio_und_connected
from bct.utils import BCTParamError,get_rng
import networkx as nx # version 2.4
from libreriaMolCr import *
from scipy.signal import savgol_filter
import matplotlib as mpl
import matplotlib.ticker as mticker
import random

def smallworld_graph(nNod, grad, pb):
    xG_ran = get_smallworld_graph(nNod, grad, pb)
    mx=nx.to_scipy_sparse_array(xG_ran)
    mod=mx.todense()
    return mod

def model_weight(io,jo,psf):
    miGgraf = nx.Graph()
    for i in range(len(io)):
        vf = np.random.uniform(min(psf),max(psf))# random.choice(psf)
        miGgraf.add_edge(io[i], jo[i], weight=vf)
    mxx=nx.to_scipy_sparse_array(miGgraf)
    xmod=mxx.todense()
    return xmod

def modelo_network(psf,nNod,grad,pb):
    
    mod = smallworld_graph(nNod, grad, pb)
    io,jo=np.where(mod!=0)
    xmod=model_weight(io,jo,psf)
    
    grado_var=np.arrange(grad,(nNod-1)/2,2)

    return xmod


def weight_network(psf,nNod,grad,pb):
    miGgraf = nx.Graph()
    xG_ran = get_smallworld_graph(nNod, grad, pb)
    mx=nx.to_scipy_sparse_array(xG_ran)
    mod=mx.todense()
    io,jo=np.where(mod!=0)

    for i in range(len(io)):
        vf = np.random.choice(psf)
        miGgraf.add_edge(io[i], jo[i], weight=vf)

    mxx=nx.to_scipy_sparse_array(miGgraf)
    xmod=mxx.todense()
    return xmod


def  graficar_comunidades(vxy,nam_exp,mostrar=False):
    
    ffile_txt=nam_exp+"w250.txt"
    
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

    G_fb = nx.from_numpy_array(vxy)

    pos_p = nx.spring_layout(G_fb)

    partitions = community_louvain.best_partition(G_fb)
    nodos_por_comundad=[]

    values = list(partitions.values())
    nNCom=list(set(values))

        
    for cluster_id in nNCom:
        cluster = [node for node in G_fb.nodes() if partitions[node] == cluster_id]
        nodos_por_comundad.append(cluster)

    qQ=nx.community.modularity(G_fb,nodos_por_comundad)
    print("Modularidad = ",qQ)
    
    if(mostrar):
        plt.axis("off")

        nx.draw_networkx(G_fb, pos=pos_p, cmap=plt.get_cmap("jet"), 
                        node_color=values, node_size=15, with_labels = False,width=0.1)
        plt.show()

def comunidades_redes(vxy):
    G_fb = nx.from_numpy_array(vxy)    
    partitions = community_louvain.best_partition(G_fb)
    nodos_por_comundad=[]

    values = list(partitions.values())
    nNCom=list(set(values))
    
    for cluster_id in nNCom:
        cluster = [node for node in G_fb.nodes() if partitions[node] == cluster_id]
        nodos_por_comundad.append(cluster)
    
    qQ=nx.community.modularity(G_fb,nodos_por_comundad)
    #print("Modularidad = ",qQ)
    return qQ
  

def randomize_matrix(A):
    if not np.allclose(A, A.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(A) > 1:
        raise BCTParamError("Input is not connected")
    
    num_nodes=len(A)
    A_rand=np.zeros((num_nodes,num_nodes))
    mask=np.triu(np.ones(num_nodes),1)
    grab_indices=np.extract(mask>0, A)
    num_edges=len(grab_indices)
    
    rand_index=np.random.permutation(num_edges)
    
    randomized_edges=np.zeros(len(grab_indices))
    
    for i in range(len(grab_indices)):
        randomized_edges[rand_index[i]]=grab_indices[i]
 
    #print(randomized_edges," ",grab_indices )
    
    edge=0
    for i in range(num_nodes-1):
        for j in range(i+1,num_nodes,1):
            A_rand[i][j]=randomized_edges[edge]
            A_rand[j][i]=randomized_edges[edge]
            edge+=1
        
    return A_rand 

def metricas_cl(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    
    G_tcm = nx.from_numpy_array(mxc)
    G_tcm_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in G_tcm.edges(data='weight')}
    nx.set_edge_attributes(G_tcm, G_tcm_distance_dict, 'distance')

    cl=nx.average_clustering(G_tcm,weight='weight')

    return cl

def metricas_pl(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    
    G_tcm = nx.from_numpy_array(mxc)
    G_tcm_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in G_tcm.edges(data='weight')}
    nx.set_edge_attributes(G_tcm, G_tcm_distance_dict, 'distance')    
    pl=nx.average_shortest_path_length(G_tcm,weight='distance')
    return pl


def fun_histogram(datos_histo,interval):
    fx, xb=np.histogram(datos_histo,bins=interval)
    sx=xb[:-1]; sy= fx/sum(fx)
    return sx, sy

def treshold(adrespr,result,ndir):

    mxrango = slice(0,9)
    kpath="{}/{}".format(adrespr,ndir) 
    conjunto_datos=os.listdir(kpath)
    trh_univ=[]
    for dw,xsubdir in enumerate(conjunto_datos):
        kfile="{}/{}/{}".format(adrespr,ndir,xsubdir[mxrango])
        reswt="{}/{}.txt".format(result,xsubdir[mxrango])
        print(reswt); print(kfile)
        dt_cv=np.loadtxt(reswt)
        mean_pl_model=[];cl_normalizada=[]; mean_network=[]; std_pl_model=[]
        mc_files=os.listdir(kfile)
        #name_spike=mc_files[0].split('_')[0]+"_"+mc_files[0].split('_')[1]
        name_spike=mc_files[0].split('_')[0]+"_Exp"+"_"+mc_files[0].split('_')[2]

        for _,u in enumerate(dt_cv[:,0]):
            load_file="{}/{}/{}/{}_{}.txt".format(adrespr,ndir,xsubdir[mxrango],name_spike,int(u))
            xpk=np.array(abs(delete_nod(np.genfromtxt(load_file))))
            trh_univ+=list(xpk[xpk!=0].flatten())
    return trh_univ

def Funct_BCTParamError(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    
 

def mean_degree(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    # #En esta seccio se calcula el grado medio de la red

    xGz = nx.from_numpy_array(mxc)
    strength = dict(xGz.degree(weight='weight'))
    normstrengthlist=list(strength.values())
    mean_degree=np.mean(normstrengthlist)
    std_degree=np.std(normstrengthlist)/np.sqrt(len(normstrengthlist))
    return mean_degree
    
    #np.round(mean_degree/2.0,2)#devuelve k/2 = k_over
def func_tiks_axes(axx,par,sd,ntks):
    axx.set_xlabel(r'$\langle CV \rangle$')#,fontsize=25)
    axx.set_ylabel('clustering')
    
    label_ylist = r'${:.1f}$'
    label_xlist = r'${:.1f}$' 
    axx.set_xlim(min(par)-sd,max(par)+sd)

    xlist=np.linspace(min(par)-sd,max(par)+sd,ntks)
   
    axx.xaxis.set_major_locator(mticker.FixedLocator(xlist))
    axx.set_xticklabels([label_xlist.format(x) for x in xlist])#,fontsize=20)
    axx.legend(fontsize='x-large')
    
    #axx.legend(loc = "lower right",fontsize='x-large')
    

def fungraf(cv,metr_real,metr_model_ran,metr_model_lat,std_model_lat,std_model_ramd,ayx):
    rxc = np.array(metr_real)-np.array(metr_model_lat)
    dC=rxc/(np.array(metr_model_lat)-np.array(metr_model_ran))
    for i in range(len(dC)):
        if(dC[i]>1):
            dC[i]=1
        elif(dC[i]<0):
            dC[i]=0
        print(dC[i])
    ayx.plot(cv,dC,"o-",color="#663399",label=r'$\langle C \rangle $ ')
    # ayx.plot(cv,metr_model_lat,"o-",color="#ff9415",label=r'${\langle \bar{C}_{lattice}\rangle}$')
    # ayx.plot(cv,metr_model_ran,"o-",color="#ff9415",label=r'${\langle \bar{C}_{ramdon}\rangle}$')
    #ayx.errorbar(cv, metr_model, std_model, fmt="none", linewidth=0.5, capsize=1,color="red")#,barsabovebool=False)
    # ayx.fill_between(cv, metr_model_ran-np.array(std_model_ramd), metr_model_ran+np.array(std_model_ramd),color="gray",alpha=0.2,hatch="") 
    # ayx.fill_between(cv, metr_model_lat-np.array(std_model_lat), metr_model_lat+np.array(std_model_lat),color="gray",alpha=0.2,hatch="") 
    #ayx.hist2d(metr_model,metr_real, bins=(10,10), cmap=plt.cm.plasma)


def fdc(C_obs,C_latt,C_ran):
    rxc=np.mean(C_latt)-C_obs
    fx=np.mean(C_latt)-np.mean(C_ran)
    if(fx ==0):
        fx=1;rxc=0
    dc=(rxc)/(fx)
    if(dc>1):
        dc=1
    if(dc<0):
        dc=0

    return dc

def fdl(L_obs,L_latt,L_ran):
    rxc=L_obs-np.mean(L_ran)
    fx=np.mean(L_latt)-np.mean(L_ran)
    if(fx==0):
        fx=1;rxc=0
    dl=(rxc)/(fx)
    if(dl>1):
        dl=1
    if(dl<0):
        dl=0

    return dl

def dldc(num_nodes,k_over_2,dL,dC,xGr,cv):
    lp_ran=[]; c_ran=[]; lp_latt=[]; c_latt=[]
    for o in range(5):
        xG_ran = get_smallworld_graph(num_nodes, k_over_2, 1)
        xG_latt = get_smallworld_graph(num_nodes, k_over_2, 0)
        matrix_ran=nx.to_scipy_sparse_array(xG_ran)
        matrix_latt=nx.to_scipy_sparse_array(xG_latt)
        modelo_ran=matrix_ran.todense()
        modelo_latt=matrix_latt.todense()
        #num_EDGES_mofrl=len(nx.from_numpy_array(modelo).edges())#//average_shortest_path_length
        c_ran.append(nx.average_clustering(nx.from_numpy_array(modelo_ran)))
        c_latt.append(nx.average_clustering(nx.from_numpy_array(modelo_latt)))
        lp_ran.append(nx.average_shortest_path_length(nx.from_numpy_array(modelo_ran)))
        lp_latt.append(nx.average_shortest_path_length(nx.from_numpy_array(modelo_latt)))

    #mean_network.append(nx.average_clustering(xGr))
    dC.append(fdc(nx.average_clustering(xGr),c_latt,c_ran))
    dL.append(fdl(nx.average_shortest_path_length(xGr),lp_latt,lp_ran))

def calculo_swp(adrespr,result,ndir):
    mxrango = slice(0,9)
    kpath="{}/{}".format(adrespr,ndir) 
    conjunto_datos=os.listdir(kpath)
    
    for dw,xsubdir in enumerate(conjunto_datos):
        kfile="{}/{}/{}".format(adrespr,ndir,xsubdir[mxrango])
        reswt="{}/{}.txt".format(result,xsubdir[mxrango])
        print(kfile)
        dt_cv=np.loadtxt(reswt)
        mean_model_ramd=[]; mean_model_latt=[]; pl_mean_network=[]; 
        rxan_model_ramd=[]; rxan_model_latt=[]; cl_mean_network=[]; 
        std_model_ramd=[]; std_model_latt=[]
        r_std_model_ramd=[]; r_std_model_latt=[]
        mc_files=os.listdir(kfile)
        #name_spike=mc_files[0].split('_')[0]+"_"+mc_files[0].split('_')[1]
        name_spike=mc_files[0].split('_')[0]+"_Exp"+"_"+mc_files[0].split('_')[2]
        
        for _,u in enumerate(dt_cv[:,0]):
            load_file="{}/{}/{}/{}_{}.txt".format(adrespr,ndir,xsubdir[mxrango],name_spike,int(u))
            mspikes=np.array(abs(delete_nod(np.genfromtxt(load_file))))
            xpk=np.zeros((len(mspikes[:,0]),len(mspikes[0,:])))
            io,jo=np.where(mspikes>0)
            xpk[io,jo]=1#mspikes[io,jo]
            xGr = nx.from_numpy_array(mspikes)
            k_over_2=int(np.round(mean_degree(xpk)/2))
            xGr = nx.from_numpy_array(mspikes)
            num_nodes=len(xGr.nodes)
            mxmod_ramd=[]; mxmod_latt=[]
            rxmod_ramd=[]; rxmod_latt=[]
            #Funct_BCTParamError(xpk)
            for o in range(50):       
                modelo_Latt=weight_network(mspikes[mspikes>0],num_nodes,k_over_2,0)               
                modelo_ramd=weight_network(mspikes[mspikes>0],num_nodes,k_over_2,1)

                mxmod_latt.append(metricas_cl(modelo_Latt))
                mxmod_ramd.append(metricas_cl(modelo_ramd))

                rxmod_latt.append(metricas_pl(modelo_Latt))
                rxmod_ramd.append(metricas_pl(modelo_ramd))
            
            mean_model_ramd.append(np.mean(mxmod_ramd))
            mean_model_latt.append(np.mean(mxmod_latt))

            rxan_model_ramd.append(np.mean(rxmod_ramd))
            rxan_model_latt.append(np.mean(rxmod_latt))

            std_model_latt.append(np.std(mxmod_latt))
            std_model_ramd.append(np.std(mxmod_ramd))
            
            r_std_model_latt.append(np.std(rxmod_latt))
            r_std_model_ramd.append(np.std(rxmod_ramd))

            cl_mean_network.append(metricas_cl(mspikes))
            pl_mean_network.append(metricas_pl(mspikes))

        rxc = np.array(mean_model_latt)-np.array(cl_mean_network)
        dC=rxc/(np.array(mean_model_latt)-np.array(mean_model_ramd))

        for i in range(len(dC)):
            if(dC[i]>1):
                dC[i]=1
            elif(dC[i]<0):
                dC[i]=0
        dC=np.array(dC)

        rxp = np.array(pl_mean_network)-np.array(rxan_model_ramd)
        dL=rxp/(np.array(rxan_model_latt)-np.array(rxan_model_ramd))

        for i in range(len(dL)):
            if(dL[i]>1):
                dL[i]=1
            elif(dL[i]<0):
                dL[i]=0
        dL=np.array(dL)

        phi=1-np.sqrt((dL*dL+dC*dC)/2)

        file_txt="SWP/"+ndir+".txt"
            
        fdata = open(file_txt, 'w')
    
        try:
            for i in range(len(phi)):
                fdata.write(str(phi[i])+"\t"+str(dC[i])+"\t"+str(dL[i])+"\t"+str(r_std_model_latt[i])+"\t"+str(r_std_model_ramd[i])+"\n")        
        finally:
            fdata.close()

        


def analisis_redes_pesadas(adrespr,result,ndir,pltf):#,foldprincp,foldresult,interval_bins):

    mxrango = slice(0,9)
    kpath="{}/{}".format(adrespr,ndir) 
    conjunto_datos=os.listdir(kpath)
    
    for dw,xsubdir in enumerate(conjunto_datos):
        kfile="{}/{}/{}".format(adrespr,ndir,xsubdir[mxrango])
        reswt="{}/{}.txt".format(result,xsubdir[mxrango])
        print(kfile)
        dt_cv=np.loadtxt(reswt)
        mean_model_ramd=[]; mean_model_latt=[]; mean_network=[]; std_model_latt=[]
        std_model_ramd=[]
        mc_files=os.listdir(kfile)
        #name_spike=mc_files[0].split('_')[0]+"_"+mc_files[0].split('_')[1]
        name_spike=mc_files[0].split('_')[0]+"_Exp"+"_"+mc_files[0].split('_')[2]
        
        for _,u in enumerate(dt_cv[:,0]):
            load_file="{}/{}/{}/{}_{}.txt".format(adrespr,ndir,xsubdir[mxrango],name_spike,int(u))
            mspikes=np.array(abs(delete_nod(np.genfromtxt(load_file))))
            xpk=np.zeros((len(mspikes[:,0]),len(mspikes[0,:])))
            io,jo=np.where(mspikes>0)
            xpk[io,jo]=1#mspikes[io,jo]
            xGr = nx.from_numpy_array(mspikes)
            k_over_2=int(np.round(mean_degree(xpk)/2))
            xGr = nx.from_numpy_array(mspikes)
            num_nodes=len(xGr.nodes)
            mxmod_ramd=[]; mxmod_latt=[]
            #Funct_BCTParamError(xpk)
            for o in range(5):       
                modelo_Latt=weight_network(mspikes[mspikes>0],num_nodes,k_over_2,0)               
                modelo_ramd=weight_network(mspikes[mspikes>0],num_nodes,k_over_2,1)
                mxmod_latt.append(metricas_pl(modelo_Latt))
                mxmod_ramd.append(metricas_pl(modelo_ramd))
            mean_model_ramd.append(np.mean(mxmod_ramd))
            mean_model_latt.append(np.mean(mxmod_latt))
            std_model_latt.append(np.std(mxmod_latt))
            std_model_ramd.append(mxmod_ramd)
            mean_network.append(metricas_pl(mspikes))
            
        cv=dt_cv[:,1]

        fungraf(cv,mean_network,mean_model_ramd,mean_model_latt,std_model_latt,std_model_ramd,pltf)
        return cv


def cargar_datos(adrespr,result,ndir,pltf):#,foldprincp,foldresult,interval_bins):

    mxrango = slice(0,9)
    kpath="{}/{}".format(adrespr,ndir) 
    conjunto_datos=os.listdir(kpath)
    
    for dw,xsubdir in enumerate(conjunto_datos):
        kfile="{}/{}/{}".format(adrespr,ndir,xsubdir[mxrango])
        reswt="{}/{}.txt".format(result,xsubdir[mxrango])
        print(kfile)
        dt_cv=np.loadtxt(reswt)
        mean_model_pl=[];cl_normalizada=[] 
        real_network_pl=[]; std_model_pl=[]
        mean_model_cl=[]; std_model_cl=[]
        dC=[];dL=[]; modular_red=[]; real_network_cl=[];
        mc_files=os.listdir(kfile)
        #name_spike=mc_files[0].split('_')[0]+"_"+mc_files[0].split('_')[1]
        name_spike=mc_files[0].split('_')[0]+"_Exp"+"_"+mc_files[0].split('_')[2]
        modular_red_modelo=[]
        for _,u in enumerate(dt_cv[:,0]):
            load_file="{}/{}/{}/{}_{}.txt".format(adrespr,ndir,xsubdir[mxrango],name_spike,int(u))
            mspikes=np.array(abs(delete_nod(np.genfromtxt(load_file))))
            xpk=np.zeros((len(mspikes[:,0]),len(mspikes[0,:])))
            io,jo=np.where(mspikes>0.01)
            xpk[io,jo]=1#mspikes[io,jo]
            k_over_2=int(np.round(mean_degree(xpk)/2))
            xGr = nx.from_numpy_array(xpk)
            num_nodes=len(xGr.nodes)
            num_EDGES_exp=len(xGr.edges())
            #Funct_BCTParamError(xpk)
            mxmod_pl=[]; mxmod_cl=[]
            #modular_red.append(nx.average_shortest_path_length(xpk))
            #graficar_comunidades(xpk,ndir,True)
            modular_modelo=[]
            for o in range(50):
                xG_ran = get_smallworld_graph(num_nodes, k_over_2, 0)
                matrix=nx.to_scipy_sparse_array(xG_ran)
                modelo_ran=matrix.todense()
                #modular_modelo.append(nx.average_shortest_path_length(modelo_ran))
                #num_EDGES_mofrl=len(nx.from_numpy_array(modelo).edges())#//average_shortest_path_length
                mxmod_pl.append(nx.average_shortest_path_length(nx.from_numpy_array(modelo_ran)))
                mxmod_cl.append(nx.average_clustering(nx.from_numpy_array(modelo_ran)))

            mean_model_pl.append(np.mean(mxmod_pl))
            std_model_pl.append(2*np.std(mxmod_pl))
            mean_model_cl.append(np.mean(mxmod_cl))
            std_model_cl.append(2*np.std(mxmod_cl))
            real_network_pl.append(nx.average_shortest_path_length(xGr))
            real_network_cl.append(nx.average_clustering(xGr))
            #modular_red_modelo.append(np.mean(mxmod))
        Pth_Legth=np.array(real_network_pl)/np.array(mean_model_pl)
        Clst=np.array(real_network_cl)/np.array(mean_model_cl)
        cv=dt_cv[:,1]
        fungraf(cv,Pth_Legth,Clst,std_model_pl,pltf)
        
        
        return cv
        
        
dir_datos="datW10T250s"
resultwt="restW10T250"

fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(2,3)
nNo=[0,1,2,0,1,2]; r=0

params = {'xtick.labelsize': 20, 'ytick.labelsize': 20}
mpl.rcParams.update(params)

experimentos=["Mar07","Mar10","Mar04","Jan29","Jan21","Jan14"]
nome_exper=["SUAM07","SUAM10","SUAM04","SUAJ29","SUAJ21","SUAJ14"]
totales=[]; treshol_universal=[]
ct=[0.03,0];limit_Y=[0,.4]
for l,xndr in enumerate(experimentos):
    print(xndr)
    #treshol_universal+=list(treshold(dir_datos,resultwt,xndr))
    

    if(l>2):
        r=1
    o=nNo[l]

    azh = fig.add_subplot(gs[r,o])
    #limit_X=analisis_redes_pesadas(dir_datos,resultwt,xndr,azh)#,foldprincp,foldresult,interval_bins)
    calculo_swp(dir_datos,resultwt,xndr)
    #configure_axes_canvas(o,azh,limit_Y,limit_X,nome_exper[l],ct)

#     totales+=list(lp)

# sx,sy=funct_histograma(totales,50)
# cmean=np.mean(totales)
# dev_mas=cmean+2*np.std(totales)
# dev_men=cmean-2*np.std(totales)
# print(dev_mas,dev_men)
# plt.tight_layout()
# plt.show()


