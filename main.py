#!/usr/bin/env python
# -*- coding: utf-8 -*

##/////////////////////////////////////////////////////##
##     Por Miguel Alejandro Molina  03/11/2022         ##
##/////////////////////////////////////////////////////##

#------------------------------------------------------------------------------------------
# importar achivos necesarios para el analisis
from pathlib import Path
from include.get_function_implementations import *

from include.insignificant_clusters import delete_insignificant_clusters 

from include.swp.small_world_propensity import *

#from include.swp.rootgraf import graf_root
from include.swp.modelos_nulos_SWP import *
from clasificar import *
import include.swp.randmio as func_randmio

from extract_experimental_sets import extract_cluster_id, extract_experimental_set

import numpy as np

import os

import statistics

import csv

from shutil import rmtree

import time

import numpy as np

import matplotlib.pyplot as plt


import seaborn as sb
#----------------------------------------------------------------------------------------------------------------
#  se llama las funciones.
#----------------------------------------------------------------------------------------------------------------
'''
#tmc_exit= np.loadtxt("ids/mc_Exit.txt",dtype=float)
#tmc_inhib = np.loadtxt(,dtype=float)

#matrix_conec = np.loadtxt("resultados/mc_pesada.txt",dtype=float)

#correlation_Excitator_Inhibitor_Matrix(id_clust_inhib,id_SUA,adress_mc,sub_cadena,solo_exit, solo_inhi)

#connectivity_matrix_linkage(tmc_inhib,tmc_exit,treshold_hard)

#resultado = numero_conexiones_tcm(matrix_conec,tmc_inhib,tmc_exit,treshold_hard,solo_exit, solo_inhi)

#print(resultado)

#----------------------------------------------------------------------------------


#metrics_set_null_networks()
#network_metrics(adress_mc)
#func_randmio(matrix_conec,10,100)



#-----------------------------------------------------------------------------------

#delete_insignificant_clusters(adress_mc,10)

#graficos()
'''



def build_connectivity_matrix():
    adres_dir = os.listdir('../principal')
    n=len(adres_dir);j=0#
    w_corr=30
    mxrango = slice(7,22)
    nxrango = slice(0,6)
    xmd = slice(11,22)
    w_shuf=20
    info_dir="ncc_0.1-1-20"#ncc_(frecuencia muestreo)-(bin)-(ventana correlacion)
    num_surr=8
    solo_exit=1
    solo_inhi=1

    nivel = ["hight","inter","lower"]

    for i in range(n):
        
        if(adres_dir[i][nxrango]=='id_sua' and j < 3):  
 
                sub_cadena = adres_dir[i][mxrango]; j+=1
                
                #print(sub_cadena)
                #if(sub_cadena == "2020Mar07_lower"):
                   
                id_SUA= np.loadtxt("id_sua_"+str(sub_cadena)+".txt") # id clusters de un experimento 
            
                id_clust_tot = np.loadtxt("Correlation/Experiments/Exp_sua_"+str(sub_cadena)+"/id_"+str(sub_cadena)+".txt")# todos los id clusters para obtenr una mc

                treshold_hard = 0 # valores entre 0 y 1. En 1 no el treshold para mc es 0
                adres_datos = "Correlation/Experiments/Exp_sua_"+str(sub_cadena)+"/cc"+(sub_cadena)+"w"+str(w_shuf) #directorio donde estan los datos de correlacion
                
                id_clust_inhib = np.loadtxt("Correlation/Experiments/Exp_sua_"+str(sub_cadena)+"/id_"+str(sub_cadena)+"_inhib.txt") # id clusters inhbitorios
                
                num_clusters =len(id_SUA)
                num_cl=len(id_clust_tot)
                
                dir_surr = "Correlation/Experiments/Exp_sua_"+str(sub_cadena)+"/sr"+str(sub_cadena)+"w"+str(w_shuf)
            
                #graf_inhib(adres_datos,138,74)
                
                #if(corrx):

                #correlation_Matrix(id_clust_tot, adres_datos,id_clust_inhib,num_clusters,num_cl,num_surr,dir_surr,w_corr,sub_cadena,True)
                correlation_Matrix_grafics(id_clust_tot,adres_datos,id_clust_inhib,num_clusters,num_cl,dir_surr,True)
                
                #return sub_cadena        

#build_connectivity_matrix()




def graphics_and_metrics(tcmx,adres):
    mxrango = slice(12,27)
    sub_cadena = adres[mxrango]
    #sub_cadena=build_connectivity_matrix(False)    
    #adress_mc = np.genfromtxt("ResMC/tmc_"+str(sub_cadena)+"_conect.txt")#cargo matriz de conectiviidad 
    #print(np.shape(tcmx))
    swp, long_path, clusterin = network_metrics(tcmx)
    #print(np.round(swp,2),np.round(long_path,2),np.round(clusterin,2))
    print(np.round(swp,2),",",sub_cadena)
    
        
    #tcm_graficos(tcmx)
    
#adress_mc = np.genfromtxt("ResMC/conectadas/tmc_2020Jan21_inter.txt")

#graphics_and_metrics(adress_mc)

adres_dir = os.listdir("Correlation/Experiments")

nsurr=100; xbin = 1; wind = 30

mxrango = slice(7,23)

for i, xdir in enumerate(adres_dir):
    
    if("mcc_2020Jan14_hight.txt" == xdir):
        
        sub_cadena = adres_dir[i][mxrango]
 
        
        dir_matriz_empirica = "Correlation/Experiments/"+str(xdir)+"/"+str(xdir)+"/"+"ncc_0.1-"+str(xbin)+"-"+str(wind)+"msec/"+str(xdir)+"_CC_Symmetric.txt"
        
        dir_datos_surrogate = "Correlation/Experiments/"+str(xdir)+"/Exp_surr/"+str(xdir)
         
        xmcc = np.genfromtxt(dir_matriz_empirica)

       
#thresholded_conetivity_matrix(xmcc,dir_datos_surrogate,nsurr,dir_resultados,xbin,wind)
        
#graphics_and_metrics(mxc)

#build_connectivity_matrix()

'''



def graphics_metrics():
    interval=np.arange(0,0.3,0.01)
    adres_dir = os.listdir("desconect")
    o=0
    adres_dir.sort()
    
    cjt= slice(14,19)
    
    #cjt= slice(16,21)

    mxrango = slice(4,19)
    
    #mxrango = slice(6,21)

    mx_cjt =["hight","inter","lower"]

    #for sk,xjt in enumerate(mx_cjt): 
    nNUmbrail=[]  
    
    for j, xdir in enumerate(adres_dir):
    
        sub_cadena = adres_dir[j][mxrango]
            
        mxdir="include/swp/modnulos/"+str(sub_cadena)


        #adres="connect/mcc_"+str(sub_cadena)+".txt"
        
        adres="desconect/mcc_"+str(sub_cadena)+".txt"
        print(adres)
        mcc = np.genfromtxt(adres)
        mxc=[]
        dir_resultados="connect/mcc_"+str(sub_cadena)+".txt"
        xmcc = np.genfromtxt(dir_resultados)
        cmpxmcc=xmcc.copy()
        
        umbral=componente_gigante(xmcc,0.1)
        
        vxy=componente_significativa(xmcc,umbral,1)
        
        xxmcc = np.triu(vxy).flatten()
        for k in range(len(xxmcc)):
            if(xxmcc[k]!=0):
                mxc.append(abs(xxmcc[k]))
        
        colores=["red","blue","green"]
        frec, xbin=np.histogram(mxc,bins=interval)
        x=xbin[:-1]; y= frec/sum(frec)
        
        
        plt.plot(x,y,linewidth=1,label=str(sub_cadena))
        plt.legend()
        if(o==2):
            #plt.savefig('fig'+sub_cadena+'_.jpg')
            plt.show()  
            o=-1   
        o += 1
      
graphics_metrics()
                
'''                  
def modelos_nulos_ponderados():
    
    adres_dir = os.listdir("connect")

 
    adres_dir.sort()
    cjt= slice(14,19)
    #cjt= slice(10,15)

    mxrango = slice(4,19)
    #mxrango = slice(0,19)

        
    for j in range(len(adres_dir)):

        sub_cadena = adres_dir[j][mxrango] 
        
        mxdir="include/swp/rtFilNull/"+str(sub_cadena)
        
        rmtree(mxdir)
        
        path = Path(mxdir)
        
        path.mkdir(parents=True)

        mxfiles="include/swp/modnulos/"+str(sub_cadena)+"/"+str(sub_cadena)
        print(mxfiles)
        rtFilNull="include/swp/rtFilNull/"+str(sub_cadena)

        mcc = np.genfromtxt("connect/mcc_"+str(sub_cadena)+".txt")
        
        for kl in range(200):
         
            mxfiles="include/swp/modnulos/"+str(sub_cadena)+"/"+str(sub_cadena)+str(kl)+".txt"
        
            mxNull = np.genfromtxt(mxfiles) 
            
            construir_enlaces(kl,rtFilNull+"/"+str(str(sub_cadena)),mcc,mxNull)
'''
def  rel_freq (x): 
    freqs = [(value, x.count (value) / len (x)) for value in set (x)] 
    return freqs
             
def modelos_nulos(nNul):
    
    adres_dir = os.listdir("connect")
    #fig, ax = plt.subplots(3, 3, figsize=(8, 8))
 
    adres_dir.sort()
    cjt= slice(14,19)
    #cjt= slice(10,15)

    mxrango = slice(4,19)
    #mxrango = slice(0,19)
    #fig, ax = plt.subplots(3, 3, figsize=(8, 8))
 
    o=0 ; r=0   
    for j in range(len(adres_dir)):
        
        sub_cadena = adres_dir[j][mxrango] 

        milistavacia=[]     
                    
        dir_resultados="connect/mcc_"+str(sub_cadena)+".txt"
        
        mcc = np.genfromtxt(dir_resultados)
        
        cpmx=mcc.copy()
        x_x=[]
        vector_k_lover,k_over,otra_opcions = mean_degree_network(mcc)
        densidad_de_red = density_network(cpmx)
        for k in range(len(vector_k_lover)):
            x_x.append(vector_k_lover[k])
            
        #print(densidad_de_red," ",sub_cadena)
        #print(k_over,"\t ",otra_opcions,sub_cadena)
        #print(k_over," ",len(cpmx)," ",sub_cadena)
        
        N=len(vector_k_lover)
        vector_k_lover.sort()
        
        x=vector_k_lover
        
        h,x1 = np.histogram( x_x,density=True)
        x=[];y=[]
        for i in range(len(rel_freq(x_x))):
            x.append(rel_freq(x_x)[i][0])
            y.append(rel_freq(x_x)[i][1])
            
        #plt.plot(x,y)    
        
        f1 = np.cumsum(y)
        # #method 2
        x2 = np.sort(vector_k_lover)
        f2 = np.array(range(N))/float(N)
        
        #plt.plot(x, f1)
        plt.plot(np.log(x2), np.log(1-f2),label=sub_cadena,linewidth=3)
        plt.ylabel('Cumulative distribution')
        plt.xlabel('Grade k')
        plt.legend()
        #plt.yscale("log")
        #plt.xscale("log")

        
        
        mxleng_path=[]
        mxcrusterig=[]
        small_world_propensity=[]
        mcMx = np.zeros((len(cpmx),len(cpmx)))
        
        mx_result=mcMx.copy()
        
        for i in range(nNul):
            
            mcMx=generate_small_world_networks(len(cpmx),k_over,0.025,np.count_nonzero(np.triu(cpmx)))
            
            mx_result=construir_enlaces(cpmx,mcMx)    

            swp, long_path, clusterin = network_metrics(cpmx)
        
            print(np.round(swp,2),",",np.round(long_path,2),",",np.round(clusterin,2),",",sub_cadena)
            small_world_propensity.append(swp)
            mxleng_path.append(long_path)
            mxcrusterig.append(clusterin)
            
            
            #print(np.round(swp,2),",",sub_cadena)
            
        
        media_mxleng_path=np.round(np.mean(mxleng_path),2) 
        media_mxcrusterig=np.round(np.mean(mxcrusterig),2) 
        
        #print(media_mxleng_path,",",media_mxcrusterig,",",sub_cadena)  
            
        
        if(o==2):
            plt.savefig('fig'+sub_cadena+'_.jpg')
            plt.show()  
            o=-1   
            r+=1
        o += 1
        
    #plt.show()              
#modelos_nulos(1)


#seccion para eliminar clusters sobrantes
 
  
#eliminate_excess_clusters("Correlation/Experiments")

#extract_experimental_set()
#extract_cluster_id()     
  
#esta seccion es para obtener enlaces inhibitorios despues del treshold'''

# a_id=np.genfromtxt("Correlation/verificando/Exp_sua_2020Mar07/Exp_sua_2020Mar07_inter/id_2020Mar07_inter_inhib.txt")
# xih= np.genfromtxt("Correlation/verificando/Exp_sua_2020Mar07/Exp_sua_2020Mar07_inter/Exp_sua_2020Mar07_inter/ncc_0.1-1-30msec/Exp_sua_2020Mar07_inter_Links_inhib.txt")
# c_tcm=np.genfromtxt("Correlation/Resultados/mcc_2020Mar07_inter.txt")

# real_inhibitory_links(c_tcm,xih,a_id)       

eliminate_excess_clusters("Experiments")


#evaluar tendencia en los datos prueba MK

adres_dir = os.listdir("tmp")

infdat= slice(0,4)
inf_dat= slice(0,38)

ri=0

# for r_k in range(4):#len(adres_dir)):
#     if('time'==adres_dir[r_k][infdat]):
#         mcc = np.genfromtxt("tmp/"+str(adres_dir[r_k][inf_dat])+".txt")
#         mxc=np.zeros(len(mcc[:,0]),dtype='uint8')
#         for i_m in range(len(mcc[:,0])):
#             mxc[i_m]=math.trunc(1000*mcc[:,0][i_m])
#         print(mxc) 

#mxc = np.genfromtxt("time_cluster_shank_sua_2020Mar4.txt")   
#Nd=mxc[0,:][-1]-mxc[0,:][0]+1
#ns_detect(mxc[0,:],Nd)  

mxc = np.genfromtxt("matriz_red.txt")   
#tcm_graficos(mxc)

#tcm_graficos(mxc)
dir_resultados="matriz_red_conectada.txt"
