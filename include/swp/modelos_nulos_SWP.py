#!/usr/bin/env python
# -*- coding: utf-8 -*

##/////////////////////////////////////////////////////##
##     Por Miguel Alejandro Molina  03/11/2022         ##
##/////////////////////////////////////////////////////##

#from smallworld.draw import draw_network
from smallworld import get_smallworld_graph
import networkx as nx
import matplotlib.pyplot as pl
import numpy as np
import matplotlib.pyplot as plt # version 3.3.2
from matplotlib.colors import ListedColormap
from matplotlib import cm
import  time
from numpy import random
import os

def delet_zeros(original):
    nueva = []
    for dato in original:
        if dato != 0:
            nueva.append(dato)
    return nueva

def construir_enlaces(mxTCM,MxNull):
    
    n=len(mxTCM)
    mxNull=np.zeros((len(MxNull),len(MxNull)))
    for i in range(len(MxNull)):
        for j in range(len(MxNull)):
            mxNull[i][j]=MxNull[i][j]

    for i in range(n*n):
        if(np.triu(mxTCM).flatten()[i]!=0):
            for k in range(n):
                for l in range(n):
                    salir_for=False
                    if(mxNull[k][l]==1.0 ):                       
                        mxNull[k][l] = np.triu(mxTCM).flatten()[i]
                        mxNull[l][k] = np.triu(mxTCM).flatten()[i]
                        salir_for=True
                        break
                if(salir_for):    
                    break

    return mxNull
    # file_txt=str(rtFilNull)+str(kl)+".txt"      
    # fdata = open(file_txt, 'w')

    # try:
    #     for zz_i in range(len(mxNull)):
    #         for zz_j in range(len(mxNull)):
    #             fdata.write( str(mxNull[zz_i][zz_j])+"\t")
    #         fdata.write("\n")
                            
    # finally:
    #     fdata.close()    

                   
def generate_small_world_networks(num_nodes,k_over_2,beta,nun_conec,grafico=False):
   
    '''
    Esta funcion genera modelos de red de mundo pequeño, 
    
    num_nodes: numero de nodos de la red a generar
    k_over_2: grado medio de la red
    beta, probabildad de conexion de largo alcence
    nN: cantidad de modelo nulos
    nun_conec: numero de conexiones de la red real
    grafico=False: si desea un grafico escriba True
    '''
    
    #for i in range(nN):
    #time.sleep()
    xG = get_smallworld_graph(num_nodes, k_over_2, beta)
    #draw_network(G,k_over_2,focal_node=0)
    #pl.subplots_adjust(wspace=0.3)
            
    matrix=nx.to_scipy_sparse_array(xG)
    mx=matrix.todense()

    ne=int(np.count_nonzero(np.triu(mx)))
    
    while(ne < nun_conec):
        xri = random.randint(0,len(mx))
        xrj = random.randint(0,len(mx))
        
        while(xri==xrj):
            xri = random.randint(0,len(mx))
            xrj = random.randint(0,len(mx))
                        
        mx[xri,xrj]=1
        mx[xrj,xri]=1   
    
        ne=int(np.count_nonzero(np.triu(mx)))   
    
    while(ne > nun_conec):
        ri = random.randint(0,len(mx))
        rj = random.randint(0,len(mx))
        
        mx[ri,rj]=0
        mx[rj,ri]=0  
        
        ne=int(np.count_nonzero(np.triu(mx)))   

    
    return mx

    # file_txt=str(addr)+str(i)+".txt"
    # print(file_txt)          
    # fdata = open(file_txt, 'w')

    # try:
    #     for zz_i in range(len(mx)):
    #         for zz_j in range(len(mx)):
    #             fdata.write( str(mx[zz_i][zz_j])+"\t")
    #         fdata.write("\n")
                            
    # finally:
    #     fdata.close() 
    
    # G = nx.from_numpy_matrix(mx)                
    # #draw_network(G,k_over_2,focal_node=0)
    # pl.subplots_adjust(wspace=0.3)
    
    # if(grafico):    
    #     pl.show()
        

