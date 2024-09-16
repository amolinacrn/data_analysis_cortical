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
    log_path=[]; logitud_camino=[]
    for os in range(len(lp)):
        log_path += list(lp[os].values())
        logitud_camino.append(list(lp[os].values()))
    log_path=np.array(log_path)[np.array(log_path)!=0]
    std_pl=np.std(log_path)/np.sqrt(len(log_path))
    path_lengh=np.mean(log_path)

    #path_lengh=nx.average_shortest_path_length(G_tcm, weight='distance') # L-obs
    #clustering = nx.average_clustering(G_tcm, weight='weight')
    clust=list(nx.clustering(G_tcm, weight='weight').values())
    mean_clust = np.mean(clust)
    std_clust = np.std(clust)/np.sqrt(len(clust))

    return  logitud_camino,path_lengh, std_pl, mean_clust, std_clust

def path_length(mxc):
    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")
    
    G_tcm = nx.from_numpy_array(mxc)
    G_tcm_distance_dict = {(e1, e2): 1 / abs(weight) for e1, e2, weight in G_tcm.edges(data='weight')}

    nx.set_edge_attributes(G_tcm, G_tcm_distance_dict, 'distance')

    lp=dict(nx.all_pairs_dijkstra_path_length(G_tcm,weight='distance'))
    
    Nnl=len(list(lp.keys()))
    logitud_camino=np.zeros((Nnl,Nnl))

    for os in range(len(lp)):
        for key in sorted(lp[os].keys()) :
            logitud_camino[os,key] =(lp[os][key])
            
    return  logitud_camino
