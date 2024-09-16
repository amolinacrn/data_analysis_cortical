#!/usr/bin/env python
# -*- coding: utf-8 -*

##//////////////////////////////////////////////////////////////////##
##     Small-World Propensity and Weighted Brain Networks           ##
##   Sarah Feldt Muldoon, Eric W. Bridgeford & Danielle S. Bassett  ##
##   Traducido de matlab Por Miguel Alejandro Molina  03/11/2022    ##
##   referencia, https://www.nature.com/articles/srep22057#Sec15     ##
##//////////////////////////////////////////////////////////////////##

#------------------------------------------------------------------------------------------


import numpy as np
import math
from bct.utils import BCTParamError
from bct.algorithms import number_of_components
import random
import sys
# Import required Python packages 
import networkx as nx # version 2.4 
import community # version 0.13 (python-louvain) 
#import gudhi # version 3.3.0 
import scipy.io # version 1.4.1 from sklearn 
import preprocessing # version 0.23.1 
import itertools 
import seaborn as sns # version 0.11.0 
import matplotlib.pyplot as plt # version 3.3.2 
import pandas as pd
import plotly.graph_objects as go
import os
from smallworld.draw import draw_network
from smallworld import get_smallworld_graph
import matplotlib.pyplot as pl
import include.swp.randmio as func_randmio
from community import community_louvain
import networkx.algorithms.community as nx_comm

def density_network(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")
    
    G = nx.from_numpy_array(mxc)
    densid=nx.density(G)
    return densid

def centralidad_intermedia(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    xG = nx.from_numpy_array(mxc)
    
    Gdistance = {(e1, e2): 1 / abs(weight) for e1, e2, weight in xG.edges(data='weight')}
    nx.set_edge_attributes(xG, Gdistance, 'distance')
    
    resp=nx.betweenness_centrality(xG,weight='distance',normalized=True)
    return resp 

def centralidad_proximidad(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    xG = nx.from_numpy_array(mxc)

    Gdistance = {(e1, e2): 1 / abs(weight) for e1, e2, weight in xG.edges(data='weight')}
    nx.set_edge_attributes(xG, Gdistance, 'distance')

    resp=nx.closeness_centrality(xG, distance='distance')
    return resp
    
def arbol_expansion_minimo(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    
    xG = nx.from_numpy_array(mxc)

    Gdistance = {(e1, e2): 1 / abs(weight) for e1, e2, weight in xG.edges(data='weight')}
    nx.set_edge_attributes(xG, Gdistance, 'distance')

    resp=nx.minimum_spanning_tree(xG, weight="distance")
    matrix=nx.to_scipy_sparse_array(resp)
    mx=matrix.todense()
    ne=int(np.count_nonzero(np.triu(mx)))
    return ne

def modularidad(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    
    xG = nx.from_numpy_array(mxc)
    
    part = community_louvain.best_partition(xG, weight="weight")
    resp=len(set(part.values()).union())

    return resp

def mean_degree_network_tot(mxc):
    # #En esta seccio se calcula el grado medio de la red

    xGz = nx.from_numpy_array(mxc)
    strength = dict(xGz.degree(weight='weight'))
    normstrengthlist=list(strength.values())
    return normstrengthlist

def mean_degree_network(mxc):
    # #En esta seccio se calcula el grado medio de la red

    xGz = nx.from_numpy_array(mxc)
    strength = dict(xGz.degree(weight='weight'))
    normstrengthlist=list(strength.values())
    mean_degree=np.mean(normstrengthlist)
    std_degree=np.std(normstrengthlist)/np.sqrt(len(normstrengthlist))
    return normstrengthlist, mean_degree, std_degree
    
    #np.round(mean_degree/2.0,2)#devuelve k/2 = k_over


def metricas_L_C(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    
    G_tcm = nx.from_numpy_array(mxc)
    G_tcm_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in G_tcm.edges(data='weight')}

    nx.set_edge_attributes(G_tcm, G_tcm_distance_dict, 'distance')


    lp=dict(nx.all_pairs_dijkstra_path_length(G_tcm,weight='distance'))
    log_path=[]
    for os in range(len(lp)):
        log_path += list(lp[os].values())
    log_path=np.array(log_path)[np.array(log_path)!=0]
    std_pl=np.std(log_path)/np.sqrt(len(log_path))
    path_lengh=np.mean(log_path)

    #path_lengh=nx.average_shortest_path_length(G_tcm, weight='distance') # L-obs
    #clustering = nx.average_clustering(G_tcm, weight='weight')
    clust=list(nx.clustering(G_tcm, weight='weight').values())
    mean_clust = np.mean(clust)
    std_clust = np.std(clust)/np.sqrt(len(clust))

    return  path_lengh, std_pl, mean_clust, std_clust

def path_lenght(mxc):
    # if not np.allclose(mxc, mxc.T):
    #     raise BCTParamError("Input must be undirected")

    # if number_of_components(mxc) > 1:
    #     raise BCTParamError("Input is not connected")
    
    G_tcm = nx.from_numpy_array(mxc)
    G_tcm_distance_dict = {(e1, e2): abs(weight) for e1, e2, weight in G_tcm.edges(data='weight')}

    nx.set_edge_attributes(G_tcm, G_tcm_distance_dict, 'distance')

    lp=dict(nx.all_pairs_dijkstra_path_length(G_tcm,weight='distance'))
    
    log_path=[]
    for os in range(len(lp)):
        log_path += list(lp[os].values());print(list(lp[os].values()))
    log_path=np.array(log_path)[np.array(log_path)!=0]
    std_pl=np.std(log_path)/np.sqrt(len(log_path))
    path_lengh=np.mean(log_path)

    return  path_lengh, std_pl

def regular_matrix_generator(G,r):

    '''
    Inputs:
        G    the adjacency matrix for the given network; must be symmmeterized
        r    the approximate radius of the regular network 

    Outputs:
        M    the regular matrix for the given network, where all 
            weights are sorted such that the inner radius has the
            highest weights randomly distributed across the nodes, 
            and so on
    '''
    n = len(G)
    G = np.triu(G)
    #reshape the matrix G into an array, B
    B = np.reshape(G,(len(G)**2,1)) 
    B=sorted(B,reverse=True)  
    num_els =np.ceil(np.size(G)/(2*n))
    num_zeros = 2*n*num_els - np.size(G)
    B = np.concatenate((B,np.zeros((math.trunc(num_zeros),1))))
    B = np.reshape(B,(-1,n)).T 
        
    M = np.zeros((n,n))
    
    if(r<=0):
        print("Error!, parámetro 'r' no válido")
        sys.exit()
    
    
    for i in range(n):
        
        for z in range(int(r)):
            
            a=round(random.uniform(0,n-1))
                        
            while((B[a][z]==0 and z != r-1) or (B[a][z] == 0 and z==r-1 and any(B[:,int(r-1)]))):
                a=round(random.uniform(0,n-1))    
            
            y_coor_1 = ((i+z+1)%(len(G)))
           
            M[i][y_coor_1] = B[a][z]
            M[y_coor_1][i] = B[a][z]
            B[a][z] = 0;  
    
    # file_txt="lattice_mc.txt"
            
    # fdata = open(file_txt, 'w')
    
    # try:
    #     for zz_i in range(n):
    #         for zz_j in range(n):
    #             fdata.write( str(M[zz_i][zz_j])+"\t")
    #         fdata.write("\n")
                    
    # finally:
    #     fdata.close()


    return M 

       
def random_lattice(A):
    '''
    Inputs:
        A           the connectivity matrix, weighted or binary
    '''
    
    #the matrix must be symmetrical and undirected
           
    if not np.allclose(A, A.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(A) > 1:
        raise BCTParamError("Input is not connected")
                        
    W=A; n = len(W)

    #compute the average degree of the unweighted network, to give
    #the approximate radius
    
    nn=len(W)
    numb_connections=np.count_nonzero(W)
    avg_deg_unw = numb_connections/nn
    avg_rad_unw = avg_deg_unw/2.0
    avg_rad_eff = np.ceil(avg_rad_unw)
    
    #se llama a las funciones
    lattice=regular_matrix_generator(W, avg_rad_eff) #latizar una red
    ramdomize=randomize_matrix(W)#randomizar una red
    
    return lattice,ramdomize
    

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
    
    # file_txt="randomize_mc.txt"
            
    # fdata = open(file_txt, 'w')
   
    # try:
    #     for zz_i in range(num_nodes):
    #         for zz_j in range(num_nodes):
    #             fdata.write( str(A_rand[zz_i][zz_j])+"\t")
    #         fdata.write("\n")
                    
    # finally:
    #     fdata.close()   
        
    return A_rand 

def SWP_phi(mxc,regular_mc,ramdom_mc):
    
    #//////////en esta funcion se calcula la propension de mundo pequeño de la red ///////////
    
    G_tcm = nx.from_numpy_array(mxc)
    G_tcm_lat = nx.from_numpy_array(regular_mc)
    G_tcm_rdm = nx.from_numpy_array(ramdom_mc)

    G_tcm_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in G_tcm.edges(data='weight')}
    G_tcm_lat_distance_dict = {(ee1, ee2): 1 / abs(weight) for ee1, ee2, weight in G_tcm_lat.edges(data='weight')}
    G_tcm_rdm_distance_dict = {(eee1, eee2): 1 / abs(weight) for eee1, eee2, weight in G_tcm_rdm.edges(data='weight')}

    #Luego agréguelos como atributos a los bordes del gráfico.

    nx.set_edge_attributes(G_tcm, G_tcm_distance_dict, 'distance')

    nx.set_edge_attributes(G_tcm_lat, G_tcm_lat_distance_dict, 'distance')

    nx.set_edge_attributes(G_tcm_rdm, G_tcm_rdm_distance_dict, 'distance')

    L_tcm=nx.average_shortest_path_length(G_tcm, weight='distance') # L-obs
    L_tcm_lat=nx.average_shortest_path_length(G_tcm_lat, weight='distance') # L-lattice
    L_tcm_rdm=nx.average_shortest_path_length(G_tcm_rdm, weight='distance') # L-ramdom

    rxp = L_tcm-L_tcm_rdm

    if(rxp < 0):
        rxp = 0
    
    dL=rxp/(L_tcm_lat-L_tcm_rdm)
    
    if(L_tcm == math.inf or L_tcm_rdm == math.inf or L_tcm_lat ==math.inf):
        dL=1
    
    if(dL>1):
        dL=1
    
    if(dL<0):
        dL=0

    clus_tcm = nx.average_clustering(G_tcm, weight='weight')
    clus_tcm_lat = nx.average_clustering(G_tcm_lat, weight='weight')
    clus_tcm_rdm = nx.average_clustering(G_tcm_rdm, weight='weight')

    rxc = clus_tcm_lat-clus_tcm
    
    if(rxc<0):
        rxc=0
    
    dC=rxc/(clus_tcm_lat-clus_tcm_rdm)
    
    if(math.isnan(clus_tcm_lat) or math.isnan(clus_tcm_rdm) or math.isnan(clus_tcm)):
        dC=0
        
    if(dC>1):
        dC=1

    if(dC<0):
        dC=0
        
    pfi=1-np.sqrt((dC*dC+dL*dL)/2)

    return pfi, L_tcm, clus_tcm
       
def small_World_Propensity(A):
    '''
    Inputs:
        A           the connectivity matrix, weighted or binary
    '''
    
    #the matrix must be symmetrical and undirected
           
    if not np.allclose(A, A.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(A) > 1:
        raise BCTParamError("Input is not connected")
                        
    W=A; n = len(W)

    #compute the average degree of the unweighted network, to give
    #the approximate radius
    
    nn=len(W)
    numb_connections=np.count_nonzero(W)
    avg_deg_unw = numb_connections/nn
    avg_rad_unw = avg_deg_unw/2.0
    avg_rad_eff = np.ceil(avg_rad_unw)
    
    #se llama a las funciones
    lattice=regular_matrix_generator(W, avg_rad_eff) #latizar una red
    ramdomize=randomize_matrix(W)#randomizar una red
    
    return SWP_phi(A,lattice,ramdomize)
   

def metrics_set_null_networks():#(matrz conectividad, cantidad de archivos nulos)
    
    '''
    Esta funcion obtiene las metricas de las matrices nulas
    y retorna un vector con las metricas deseadas
    '''
    
    adres_dir = os.listdir('include/swp/datrandom')
    n=len(adres_dir)

    #file_txt="cfswp.txt"
            
    #fdata = open(file_txt, 'w')
    
    coeficientes=[]
    mx_long_path =[]
    mx_clustering =[]
    
    for i in range(n):
        mc_conectividad=np.genfromtxt("include/swp/datrandom/rdm"+str(i)+".txt")#"MC_pesada.txt")  
        mtric_swp, mtric_Lp,mtric_cl = small_World_Propensity(mc_conectividad)
        coeficientes.append(mtric_swp)
        mx_long_path.append(mtric_Lp)
        mx_clustering.append(mtric_cl)
        
    '''        
    try:
        for zz_i in range(len(coeficientes)):
                fdata.write( str(coeficientes[zz_i])+"\n")                        
    finally:
        fdata.close()   
    '''
    return coeficientes,mx_long_path,mx_clustering

def network_metrics(mc_conectividad):#(matrz conectividad, cantidad de archivos nulos)
    '''
    calcular metrica de una matriz de conectividad deseada
    '''
    
    s_w_p=[]
    mx_long_path =[]
    mx_clustering =[]

    mtric_swp, mtric_Lp,mtric_cl = small_World_Propensity(mc_conectividad)
    
    s_w_p.append(mtric_swp)
    mx_long_path.append(mtric_Lp)
    mx_clustering.append(mtric_cl)
    
    return s_w_p, mx_long_path, mx_clustering


def netSWP(mc_conectividad):#(matrz conectividad, cantidad de archivos nulos)
    '''
    calcular metrica de una matriz de conectividad deseada
    '''
    
    s_w_p=[]
    mx_long_path =[]
    mx_clustering =[]

    mtric_swp, mtric_Lp,mtric_cl = small_World_Propensity(mc_conectividad)
    
    s_w_p.append(mtric_swp)
    
    return s_w_p