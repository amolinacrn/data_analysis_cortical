from cProfile import label
from enum import unique
from statistics import mean
import tempfile
from traceback import print_tb
from turtle import color
import numpy
import math
import re
import pandas as pd
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from libreriaMolCr import delete_nod,modularidad,metricas_L_C,path_length
import seaborn as sns
import networkx as nx
from pathlib import Path
from community import community_louvain
import matplotlib as mpl
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def mean_degree_network(mxc):
    # #En esta seccio se calcula el grado medio de la red

    xGz = nx.from_numpy_array(mxc)
    strength = dict(xGz.degree(weight='weight'))
    normstrengthlist=list(strength.values())
    mean_degree=np.mean(normstrengthlist)
    std_degree=np.std(normstrengthlist)/np.sqrt(len(normstrengthlist))
    return normstrengthlist, mean_degree, std_degree

def fun_histo(datos_histo,interval):
    fx, xb=np.histogram(datos_histo,bins=interval)
    sx=xb[:-1]; sy= fx/sum(fx)
    return sx,sy

def dist_ucl(X,Y):
    distancia_Euclidiana=np.zeros((len(X),len(X)))
    for k in range(len(X)-1):
        for l in range(k+1,len(X)):
            d=math.sqrt((X[k]-X[l])*(X[k]-X[l]) + (Y[k]-Y[l])*(Y[k]-Y[l]))
            distancia_Euclidiana[k,l]=d
            distancia_Euclidiana[l,k]=d
    return distancia_Euclidiana

def  rel_freq (x): 
    freqs = [(value, x.count (value) / len (x)) for value in set (x)] 
    return freqs


def mi_colors(cv):
    viridis = mpl.colormaps['plasma'].resampled(len(cv))
    newcolors = viridis(np.linspace(0,1, len(cv)))
    return viridis

def simil_cos(xx,yy):
    cos_sim=np.zeros((len(xx),len(yy)))
    for i in range(len(yy)-1):
        va=[xx[i], yy[i]]
        norm_a=np.sqrt(np.dot(va,va))
        for j in range(i+1,len(yy)):
            vb=[xx[j], yy[j]]
            norm_b=np.sqrt(np.dot(vb,vb))
            pp=norm_a*norm_b
            rest=(np.dot(va,vb)/pp)
            cos_sim[i,j]=1-(rest) 
            cos_sim[j,i]=1-(rest)
    return cos_sim

def graf_dist_fisica(impDat,importar_datos,dats):#otrdist,impDat,importar_datos,dats):
        #for ks in range(1,2):
        meanMC=[]; meanEUCL=[]; todwij=[]; xtwij=[]
        cv=[];modul=[];nNf=[];theta=[];ejey=[];ejex=[]
        log_camin=[];error_kover=[];fxc=[]
        newcolor=mi_colors(impDat[:,0])
        for l,i in enumerate(impDat[:,0]):#range(len(impDat[:,0])):
            # fig, zx = plt.subplots(1, 2, figsize=(5,5 )
            mcnsx=[]; mcnsy=[]
            kpath="datW10T250s/{}/250s{}/mcc_Exp_{}_{}.txt".format(dats,dats,dats,int(i)) 
            #kpath="datW10T50s/{}/50s{}/mcc_2020{}_{}.txt".format(dats,dats,dats,int(i)) 
            matMC=abs(delete_nod(np.loadtxt(kpath)))
            print(kpath)


            #matMC=abs(delete_nod(np.loadtxt("50s"+str(dats)+"/mcc_2020"+str(dats)+"_"+str(int(impDat[i,0]))+".txt")))
            #print("30sMar04/mcc__2020Mar4_"+str(int(impDat[i,0]))+".txt")
            # xGr = nx.from_numpy_array(matMC)
            # nN=len(xGr.nodes)
            # if(nN>10):
            id_no_silen=np.loadtxt("id_notsilen_W10T250/no_silenciosoExp_"+str(dats)+"_"+str(int(i))+".txt")
            
            #id_no_silen=np.loadtxt("id_notsilen_W10T50/no_silencioso2020"+str(dats)+"_"+str(int(impDat[i,0]))+".txt")
            VmatMC=np.zeros((len(matMC[:,0]),len(matMC[0,:])))
            matrz_dist=np.zeros((len(matMC[:,0]),len(matMC[0,:])))

            for j,datj in enumerate(id_no_silen):
                for k,datk in enumerate(importar_datos[:,3]):
                    if(datj==datk):
                        mcnsx.append(importar_datos[k,0])
                        mcnsy.append(importar_datos[k,1])
                        break

            
            frac_mdcl=dist_ucl(list(importar_datos[:,0]),list(importar_datos[:,1]))
            mdcl=dist_ucl(mcnsx,mcnsy)

            #io,jo = np.where(matMC)
            ko,lo =np.where(mdcl!=0)
            
            vf=matMC[ko,lo]
            vp=mdcl[ko,lo]
            #matrz_dist[ko,lo]=1
            # print("id_notsilen_W10T50/no_silencioso2020"+str(dats)+"_"+str(int(i))+".txt")
            # print(len(matMC),"\t",len(mdcl))

            #matrz_dist[io,jo]=mdcl[io,jo]
            # meanMC += list(matMC[io,jo])
            # meanEUCL += list(mdcl[io,jo])
            
            # meanMC.append(np.mean((matMC)))
            #meanEUCL.append(mdcl)
        
            # w_ij=[]; #VmatMC=matMC
                        
            # vector_k_lover, nkover,_= mean_degree_network(VmatMC)
            
            # for kn in range(len(VmatMC[:,0])):
            #     w_ij.append(VmatMC[kn,:]/vector_k_lover[kn])
                
            #logitud_camino=path_length(matMC)
            # W_ij=np.array(w_ij)
           
            # vector_k_lover, xnkover,_= mean_degree_network(W_ij)
            # print(xnkover)
            #W_ij=np.array(logitud_camino)
            #print(vector_k_lover)
            # angulo=matMC[io,jo]/mdcl[io,jo]
            # xtheta=(list(np.arctan(angulo)))
            
            # ejey+=list((matMC[io,jo]))
            # ejex.append(0)
            
            #xtwij+=list(np.mean(matMC[matMC>0]))
            
            #todwij.append(nkover)#np.mean(W_ij))
            # print(todwij)
            # cv.append(impDat[l,1])

            # vect_dist_k, mean_dist_k,_= mean_degree_network(matrz_dist)
            
            # meanEUCL.append(np.mean(1/mean_dist_k))
            #sxz,syz=fun_histo(xtheta,100)
            #tx=np.concatenate(theta)
            #plt.hist2d()
        # for i in range(len(cv)):
        #     plt.plot(todwij[i],meanEUCL[i],"d",color=newcolor.colors[i])
        # plt.yscale("log")
        # plt.xscale("log")
        # plt.title(dats,fontsize="x-large")
        # plt.show()

        # vx,vy=fun_histo(theta,3000)
        # plt.plot(vx,vy,"-")
        
        #return cv,ejex,ejey, error_kover,theta#,log_camin,(todwij)

def smil_dist_fisica(impDat,importar_datos,dats):#otrdist,impDat,importar_datos,dats):
    msim=[]
    path="restW10T250/250s{}.txt".format(dats) 
    cv=np.loadtxt(path)[:,1]
    for l,i in enumerate(impDat[:,0]):#range(len(impDat[:,0])):
        # fig, zx = plt.subplots(1, 2, figsize=(5,5 )
        mcnsx=[]; mcnsy=[]
        kpath="datW10T250s/{}/250s{}/mcc_Exp_{}_{}.txt".format(dats,dats,dats,int(i)) 
        #kpath="datW10T50s/{}/50s{}/mcc_2020{}_{}.txt".format(dats,dats,dats,int(i)) 
        matMC=abs(delete_nod(np.loadtxt(kpath)))
        #matMC=abs(delete_nod(np.loadtxt(kpath)))
        #print(kpath)
        id_no_silen=np.loadtxt("id_notsilen_W10T250/no_silenciosoExp_"+str(dats)+"_"+str(int(i))+".txt")
        
        #id_no_silen=np.loadtxt("id_notsilen_W10T50/no_silencioso2020"+str(dats)+"_"+str(int(impDat[i,0]))+".txt")
        VmatMC=np.zeros((len(matMC[:,0]),len(matMC[0,:])))
        matrz_dist=np.zeros((len(matMC[:,0]),len(matMC[0,:])))

        for j,datj in enumerate(id_no_silen):
            for k,datk in enumerate(importar_datos[:,3]):
                if(datj==datk):
                    mcnsx.append(importar_datos[k,0])
                    mcnsy.append(importar_datos[k,1])
                    break

        mdcl=dist_ucl(mcnsx,mcnsy)

        io,jo =np.where(matMC!=0)
        msim.append(np.mean(mdcl[io,jo]))

        #f,b=np.histogram(xtheta,100)


        plt.plot(np.log(matMC[io,jo]),np.log(mdcl[io,jo]),"o")
    # print(m)
    #plt.imshow(m)
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.show()
        
