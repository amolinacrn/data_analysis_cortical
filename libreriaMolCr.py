import numpy as np
from bct.algorithms import number_of_components
from bct.utils import BCTParamError,get_rng
import math
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
from community import community_louvain
import networkx.algorithms.community as nx_comm
import matplotlib.ticker as mticker

def delete_nod(M):

    if number_of_components(np.abs(M)) > 1:

        #print("Matriz desconectada \t")

        for j in range(len(M[:,0])):
                
            for i in range(len(M[:,0])):
                if(any(M[:,int(i)])!=True):
                    A=np.delete(M, i, axis=1)
                    break  
                
            for i in range(len(M[:,0])):    
                if(any(A[int(i),:])!=True):
                    B=np.delete(A, i, axis=0) 
                    break
        
            M = np.resize(M,(len(B),len(B)))
            M=B; contador=-1
            
            for a_i in range(len(M[:,0])):   
                if((any(M[int(a_i),:])==True) and (any(M[:,int(a_i)]))==True):
                    contador+=1
            
            if(a_i==contador):
                break
        
   
        return M
    else:
        
        #print("Matriz conectada \t") 
        
        return M

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

def mean_degre(mxc):
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
    #return normstrengthlist, mean_degree, std_degree
    
    #np.round(mean_degree/2.0,2)#devuelve k/2 = k_over
    return mean_degree


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

    return  clust,log_path,path_lengh, std_pl, mean_clust, std_clust


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
    
    return M 

       
def red_reticulada(A):
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

    
    return lattice

def funct_histograma(datos_histo,interval):
    fx, xb=np.histogram(datos_histo,bins=interval)
    sx=xb[:-1]; sy= fx/sum(fx)
    return sx,sy

def configure_axes_canvas(ko,axhx,dy0,zrtx,nam_experiment,cx):
    
    axhx.set_title(nam_experiment,fontsize="x-large")
    
    axhx.legend(loc = "upper right",fontsize='x-large')
    
    label_ylist = r'${:.0f}$'
    label_xlist = r'${:.0f}$'
    
    epsilon=np.array([cx,cx,cx])
    ylist=np.linspace(min(dy0)-epsilon[ko,1],max(dy0)+epsilon[ko,1],7)
    xlist=np.linspace(min(zrtx)-epsilon[ko,0],max(zrtx)+epsilon[ko,0],7)
    
    rxm=np.array([[r"Variation coefficient",r"$\langle CV \rangle$"],
            [r"Variation coefficient",r"$\langle CV \rangle$"],
            [r"Variation coefficient",r"$\langle CV \rangle$"]])

    axhx.xaxis.set_major_locator(mticker.FixedLocator(xlist))
    axhx.set_xticklabels([label_xlist.format(x) for x in xlist],fontsize=14)

    axhx.yaxis.set_major_locator(mticker.FixedLocator(ylist))
    axhx.set_yticklabels([label_ylist.format(x) for x in ylist],fontsize=14)

    axhx.set_ylabel(r""+str(rxm[ko,0])+"",fontsize=20)
    axhx.yaxis.set_tick_params(labelsize=20)

    axhx.set_xlabel(r""+str(rxm[ko,1])+"",fontsize=20)
    axhx.xaxis.set_tick_params(labelsize=20)

   
    axhx.set_xlim(min(zrtx)-epsilon[ko,0],max(zrtx)+epsilon[ko,0])
    axhx.set_ylim(min(dy0)-epsilon[ko,1],max(dy0)+epsilon[ko,1])
