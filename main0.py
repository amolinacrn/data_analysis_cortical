#!/usr/bin/env python
# -*- coding: utf-8 -*

##/////////////////////////////////////////////////////##
##     Por Miguel Alejandro Molina  03/11/2022         ##
##/////////////////////////////////////////////////////##

#------------------------------------------------------------------------------------------
# importar achivos necesarios para el analisis

#from ..include.tcm_analisis import correlation_Matrix
#from ..include.swp.small_world_propensity import network_metrics, metrics_set_null_networks
#from include.swp.rootgraf import graf_root
#from ..include.swp.modelos_nulos_SWP import generate_small_world_networks
#import ..include.swp.randmio as func_randmio
from cProfile import label
from enum import Flag
from extract_experimental_sets import extract_cluster_id, extract_experimental_set
import numpy as np
import os
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 10:52:38 2022

@author: CSN Admin
"""
from turtle import pos
import scipy.io
import os
from os import PRIO_PGRP, listdir
from os.path import isfile, join
import h5py
import numpy as np

import scipy.io
import os
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import numpy.matlib
import scipy.optimize as opt
import random
import pandas as pd
import matplotlib.gridspec as gridspec
#import seaborn as sns
from scipy import stats
import glob
from glob import glob
from pathlib import Path
from scipy.stats import entropy
from deteccion_comunitaria import *
from LibreriaMolCr import *
from collections import Counter
from mpl_toolkits.mplot3d import Axes3D
import re

class XNAT2:
    def __init__(self, data_original):
        self.data = data_original[0,:] #leitura do arquivo
        self.iD = data_original[1,:]
        #print('File loaded successfully!')


    def coef_variacion(self,cv_time,bin_cv, XADFA=False):
        #Retornar CV
        self.cv_time = cv_time
        #janela para o cálculo do CV, normalmente 10s
        self.bin_cv = bin_cv


        self.n_cv = int(np.round((self.data[-1]-self.data[0])/cv_time))

        self.CV = []#np.zeros(self.n_cv)
        self.SIZES = []
        self.DURATIONS = []
        self.bin_sizes = []
        matriz_spikes=[]
        matriz_id=[]
        if XADFA:
            self.ADFA = True
            self.DFA = np.zeros(self.n_cv)
        else:
            self.ADFA = False

        for t in range(self.n_cv-25):
           # print('janela: ', t)
            a = t*cv_time + self.data[0]
            b = a + cv_time
            spk = self.data[(self.data>=a) & (self.data<=b)]
            id_spk = self.iD[(self.data>=a) & (self.data<=b)]

            if (len(spk)==0 or len(spk)==1 or len(spk)==2):
                print('Janela vacia: ', t)
            #     self.CV.append(100)
            #     matriz_spikes.append(-100)
            #     matriz_id.append(-100)
            #    # print(spk)
                continue

            matriz_spikes.append(spk.tolist())
            matriz_id.append(id_spk.tolist())

            nbins = np.arange(a,b,bin_cv)#numpy.arange([start, ]stop, [step, ]dtype=None)
            count, edges=np.histogram(spk, bins = nbins)

            self.CV.append(np.std(count)/np.mean(count))

        self.CV = np.array(self.CV)

        return np.array(range(len(self.CV)))*cv_time,self.CV

    def CV_calc(self,cv_time,bin_cv, show=False, save=False, ADFA=False, showbin=False, bin_size=None):
        #Retornar CV
        self.cv_time = cv_time
        #janela para o cálculo do CV, normalmente 10s
        self.bin_cv = bin_cv
        #binarização para gerar a serie temporal para o cálculo de CV,
        #em geral 0.05 = 50 ms

        #Configurações das figuras
        #plt.rcParams["figure.figsize"] = [10,6]
        #plt.rcParams.update({'font.size': 18})


        self.n_cv = int(np.round((self.data[-1]-self.data[0])/cv_time))
        # O número de cv's é equivalente ao tempo do último spike menos o primeiro spike dividido pelo bin de 10 segundos
        #print("Total CV's = ", self.n_cv)

        self.CV = []#np.zeros(self.n_cv)
        self.SIZES = []
        self.DURATIONS = []
        self.bin_sizes = []
        matriz_spikes=[]
        matriz_id=[]
        if ADFA:
            self.ADFA = True
            self.DFA = np.zeros(self.n_cv)
        else:
            self.ADFA = False

        for t in range(0,self.n_cv):
           # print('janela: ', t)
            a = t*cv_time + self.data[0]
            b = a + cv_time
            spk = self.data[(self.data>=a) & (self.data<=b)]
            id_spk = self.iD[(self.data>=a) & (self.data<=b)]

            if (len(spk)==0 or len(spk)==1 or len(spk)==2):
                print('Janela vacia: ', t)
                self.CV.append(100)
                matriz_spikes.append(-100)
                matriz_id.append(-100)
               # print(spk)
                continue

            matriz_spikes.append(spk.tolist())
            matriz_id.append(id_spk.tolist())

            nbins = np.arange(a,b,bin_cv)#numpy.arange([start, ]stop, [step, ]dtype=None)
            count, edges=np.histogram(spk, bins = nbins)

            self.CV.append(np.std(count)/np.mean(count))

        #print('<CV>: {:.2f} ± {:.2f}'.format(np.mean(self.CV), np.std(self.CV)))
        #print('<ISI>: {:.2f} ± {:.2f} (ms)'.format(np.mean(self.bin_sizes)*1000, np.std(self.bin_sizes)*1000))
        self.CV = np.array(self.CV)
        self.I = np.argsort(self.CV)
        self.n_cv = len(self.CV)
        if show:
            plt.plot(np.array(range(len(self.CV)))*cv_time, self.CV)
            plt.ylabel('CV')

            plt.xlabel('t (s)')

            plt.tight_layout()
            if save:
                plt.savefig('CV_t.png', dpi = 150)
            plt.show()
            n, bins, patches = plt.hist(self.CV, bins = 50, density=1, facecolor='green', alpha=1)
            plt.title('CV histogram')
            plt.xlabel("CV")
            plt.ylabel("Density")
            plt.tight_layout()
            if save:
                plt.savefig('hist_CV_t.png', dpi = 150)
            plt.show()

        return(self.CV,matriz_spikes,matriz_id)


    def similitud_coseno(self,cv_time,bin_cv,nexp,u,ad,ADFA=False):
        #Retornar CV
        self.cv_time = cv_time
        #janela para o cálculo do CV, normalmente 10s
        self.bin_cv = bin_cv
        #binarização para gerar a serie temporal para o cálculo de CV,
        #em geral 0.05 = 50 ms

        #Configurações das figuras
        #plt.rcParams["figure.figsize"] = [10,6]
        #plt.rcParams.update({'font.size': 18})


        self.n_cv = int(np.round((self.data[-1]-self.data[0])/cv_time))
        # O número de cv's é equivalente ao tempo do último spike menos o primeiro spike dividido pelo bin de 10 segundos
        #print("Total CV's = ", self.n_cv)

        self.CV = []#np.zeros(self.n_cv)
        self.SIZES = []
        self.DURATIONS = []
        self.bin_sizes = []
        matriz_spikes=[]
        matriz_id=[]
        if ADFA:
            self.ADFA = True
            self.DFA = np.zeros(self.n_cv)
        else:
            self.ADFA = False
        distrib=[];
        Long_Path=[]; meansimil=[]; Mean_Degree=[]; Modul=[]; comunid=[];xcv=[]
        count_frec=[]; numC=[]; count_frec_extend=[];  nx_count_frec=[]
        #fig, zx = plt.subplots(1, 2, figsize=(15,5 ))

        for t in range(self.n_cv-25):
           # print('janela: ', t)
            ar=self.data[0]; br=ar + cv_time
            binfijo=np.arange(ar,br,bin_cv)
            nN=len(binfijo)
            a = t*cv_time + self.data[0]
            b = a + cv_time
            spk = self.data[(self.data>=a) & (self.data<=b)]
            id_spk = self.iD[(self.data>=a) & (self.data<=b)]

            if (len(spk)==0 or len(spk)==1 or len(spk)==2):
                print('Janela vacia: ', t)
                # self.CV.append(100)
                matriz_spikes.append(-100)
                matriz_id.append(-100)
               # print(spk)
                continue

            matriz_spikes.append(spk.tolist())
            matriz_id.append(id_spk.tolist())

            nbins = np.arange(a,b,bin_cv)
            count, nbin=np.histogram(spk, bins = nbins)

            self.CV.append(np.std(count)/np.mean(count))

            vect=[]; nx=[]
            for n,c in Counter(count).items():
                #if(n!=0):
                vect.append(c)
                nx.append(n)
            count_frec.append(vect)
            count_frec_extend.extend(nx)
            nx_count_frec.append(nx)

        valr_max=max(count_frec_extend)

        num_dat_frec=np.zeros((len(nx_count_frec),valr_max+1))

        for r in range(len(nx_count_frec)):
            xm=np.array(nx_count_frec[r])
            rtx=np.array(count_frec[r])
            for g in range(valr_max+1):
                io=np.where(xm==g)[0]
                if(len(xm[io])!=0):
                    num_dat_frec[r,g]=rtx[io][0]
                else:
                    num_dat_frec[r,g]=0
        return num_dat_frec
        #plot_cv(num_dat_frec,self.CV,nexp,ad)

        # cos_sim=np.zeros((len(num_dat_frec[:,0]),len(num_dat_frec[:,0])))
        # xcosim=np.zeros((len(num_dat_frec[:,0]),len(num_dat_frec[:,0])))

        # for i in range(len(num_dat_frec[:,0])-1):
        #     for j in range(i+1,len(num_dat_frec[:,0])):
        #         va=num_dat_frec[i,:]/sum(num_dat_frec[i,:])
        #         vb=num_dat_frec[j,:]/sum(num_dat_frec[j,:])
        #         rest=D_js(va,vb)#np.dot(va,vb)/(norm(va)*norm(vb))
        #         cos_sim[i,j]=(rest)
        #         cos_sim[j,i]=(rest)

        #xcosim=np.zeros((len(self.CV),len(self.CV)))


        #cv=similitud_cvs(self.CV)

        # plt.imshow(cv)
        # plt.show()

        #Hu=humbral(cos_sim)
        #tH=[0.06059,0.02707,0.0242,0.03015,0.02738,0.04303]
        #tHcv=[0.07264,0.09648,0.03461,0.11621,0.39662,0.43975]
        #xo,yo=np.where((cos_sim<Hu) & (cos_sim>0))
        #xcosim[xo,yo]=1#cos_sim[xo,yo]
        #xcosim=delete_nod(xcosim)
        #graficar_matrices_similtud(xcosim,nexp,ad)
        #graf_red(xcosim,nexp,True,False)
        #plot_comunidades(xcosim,nexp)

def similitud_cvs(cvs):

    c=np.zeros((len(cvs),len(cvs)))

    for i in range(len(cvs)-1):
        for j in range(i+1,len(cvs)):
            zs=abs(cvs[i]-cvs[j])
            c[i,j] = zs
            c[j,i] = zs
    return c

def graficar_matrices_similtud(mx,naExp,ax):
    ax.set_title(naExp,fontsize="x-large")
    #ax.imshow(mx,cmap="afmhot")
    ax.pcolormesh(mx, cmap="afmhot")
    label_ylist = r'${:.0f}$'
    label_xlist = r'${:.0f}$'
    xlist=np.linspace(0,len(mx[0,:]),7)
    ylist=np.linspace(0,len(mx[0,:]),7)
    titulos_ejes=[r"$i$","$j$ "]
    lab_format=[label_xlist,label_ylist]
    lab_list_1=[xlist,ylist]

    xmicol1="afmhot"
    xcmap = mpl.colormaps[xmicol1].resampled(250)
    xnewcolors = xcmap(np.linspace(0,1,250))
    xnewcmp = ListedColormap(xnewcolors)

    num_ejex =np.linspace(0,len(mx[0,:]),7)
    list_numers_x=np.linspace(0,len(mx[0,:]),7)
    num_ejey =np.linspace(len(mx[0,:]),0,7)
    list_numers_y=np.linspace(len(mx[0,:]),0,7)

    ax.set_xlim(0,len(mx[0,:]))
    ax.set_ylim(len(mx[0,:]),0)
    ax.xaxis.set_major_locator(mticker.FixedLocator(num_ejex))
    ax.set_xticklabels([label_xlist.format(x) for x in list_numers_x],fontsize=20)

    ax.yaxis.set_major_locator(mticker.FixedLocator(num_ejey))
    ax.set_yticklabels([r'${:.0f}$'.format(x) for x in list_numers_y],fontsize=20)

    camvax(ax,lab_format,lab_list_1,titulos_ejes)
    mibarcax(mx.flatten(),ax,"afmhot",r'${:.2f}$',np.linspace(0,np.max(mx),6),r"$D_{JS}(i,j)$")



def  plot_comunidades(mxc,nam_exp):

    if not np.allclose(mxc, mxc.T):
        raise BCTParamError("Input must be undirected")

    if number_of_components(mxc) > 1:
        raise BCTParamError("Input is not connected")

    ffile_txt=nam_exp+"_QSil.txt"

    ffdata = open(ffile_txt, 'w')

    try:
        ffdata.write("{")
        for wrfil in range(len(mxc)):
            ffdata.write("{")
            for wrcol in range(len(mxc)):
                ffdata.write(str(mxc[wrfil][wrcol]))
                if(wrcol==len(mxc)-1 and wrfil < len(mxc)-1):
                    ffdata.write("},")
                elif(wrfil == len(mxc)-1 and wrcol==len(mxc)-1):
                    ffdata.write("}")
                else:
                    ffdata.write(",")
        ffdata.write("}")
    finally:
        ffdata.close()

def plot_cv(num_dat_frec,cv,nam_experiment,ax):
    ax.set_title(nam_experiment,fontsize="x-large")
    '''
    d=[]
    zx=np.delete(num_dat_frec,0, axis=1)
    for k in range(len(zx[:,0])):
        d.append(np.mean(zx[k,:]))
    fcv=(np.array(d)-np.min(np.array(d)))/(np.max(np.array(d))-np.min(np.array(d)))
    rf=(np.array(num_dat_frec[:,0])-np.min(np.array(num_dat_frec[:,0])))/(np.max(num_dat_frec[:,0])-np.min(np.array(num_dat_frec[:,0])))
    ax[1].plot(cv,fcv,"o",color="#9467bd")
    ax[1].plot(cv,rf,"o",color="#8c564b")
    '''
    xcv = np.array(cv)
    cv=sorted(cv)
    twf=[]

    for o,ix in enumerate(cv):
        zo=np.where(xcv==ix)[0]
        zr=num_dat_frec[zo,:][0]/sum(num_dat_frec[zo,:][0])
        twf.append(zr)
    xft=np.array(twf)

    x=np.linspace(0,len(xft[0,:]),len(xft[0,:]))
    for i in range(len(cv)):
        ax.plot(x,xft[i,:])

    label_ylist = r'${:.2f}$'
    label_xlist = r'${:.0f}$'
    xlist=np.linspace(0,len(xft[0,:]),8)
    ylist=np.linspace(0,np.max(xft),10)
    titulos_ejes=[r"$\rho(\eta)$","$\eta$ "]
    lab_format=[label_xlist,label_ylist]
    lab_list_1=[xlist,ylist]

    '''
    ax.imshow(np.transpose(xft),cmap="afmhot",aspect="auto")
    ax.invert_yaxis()

    xmicol1="afmhot"
    xcmap = mpl.colormaps[xmicol1].resampled(250)
    xnewcolors = xcmap(np.linspace(0,1,250))
    xnewcmp = ListedColormap(xnewcolors)

    num_ejex =np.linspace(0,len(cv),8)
    list_numers_x=np.linspace(min(cv),max(cv),8)
    num_ejey =np.linspace(0,len(xft[0,:])-2,8)
    list_numers_y=np.linspace(0,len(xft[0,:]),8)
    print(list_numers_y)
    ax.set_xlim(0,len(cv))
    ax.set_ylim(-0.5,len(xft[0,:])-2)
    ax.xaxis.set_major_locator(mticker.FixedLocator(num_ejex))
    ax.set_xticklabels([label_xlist.format(x) for x in list_numers_x],fontsize=20)

    ax.yaxis.set_major_locator(mticker.FixedLocator(num_ejey))
    ax.set_yticklabels([r'${:.0f}$'.format(x) for x in list_numers_y],fontsize=20)

    ax.set_ylabel(r"titles[0]",fontsize=20)
    ax.yaxis.set_tick_params(labelsize=20)
    ax.set_xlabel(r"CV",fontsize=20)
    '''
    camvax(ax,lab_format,lab_list_1,titulos_ejes)
    #mibarcax(xft.flatten(),ax,"afmhot",r'${:.0f}$',np.linspace(0,np.max(xft),8),r"$CV$")


def mibarcax(cv,pargraf,cmap,label_format,ticks_loc,xcvd):
    norm= matplotlib.colors.Normalize(vmin=min(cv), vmax=max(cv))
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    ax2_divider = make_axes_locatable(pargraf)
    cax0 = ax2_divider.append_axes("bottom", size="4%", pad="35%")
    cbar0=fig.colorbar(sm,extend='both' ,cax=cax0,orientation='horizontal')
    cbar0.ax.set_xlabel(xcvd,fontsize=20)
    cbar0.ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar0.ax.set_xticklabels([label_format.format(x) for x in ticks_loc],fontsize=20)
    cax0.xaxis.set_ticks_position("bottom")

def graf_red(mx,nexp,xwrite,graf):

    xG = nx.from_numpy_array(mx)

    pos_p = nx.spring_layout(xG)

    part = community_louvain.best_partition(xG, weight="weight")
    resp=len(set(part.values()).union())

    nodos_por_comundad=[]
    values = list(part.values())
    nNCom=list(set(values))


    for cluster_id in nNCom:
        cluster = [node for node in xG.nodes() if part[node] == cluster_id]
        nodos_por_comundad.append(cluster)

    qQ=nx.community.modularity(xG,nodos_por_comundad)
    print("Modularidad = ",qQ)

    if(xwrite):

        try:
            xo=0
            file_txt=nexp+"SilPy.txt"
            print(file_txt)
            fdata = open(file_txt, 'w')

            for cluster_id in nNCom:
                cluster = [node for node in xG.nodes() if part[node] == cluster_id]
                nodos_por_comundad.append(cluster)
                for wrfil in range(len(cluster)):
                    fdata.write(str(cluster_id)+"\t"+str(cluster[wrfil])+"\n")#+str(int(sdt[xo]))+"\n")
                    xo+=1

        finally:
            fdata.close()

    if(graf):
        nx.draw_networkx(xG, pos=pos_p, cmap=plt.get_cmap("jet"),
                        node_color=values, node_size=10, with_labels = False,width=0.05)
        plt.show()


def humbral(xdm):
    mxD=np.array(xdm)
    vd=np.array(np.triu(mxD).flatten())
    umbral_min=np.min(vd[vd>0])
    par_red=np.max(vd[vd>0])
    umbral,xDaux=componente_gigante(mxD,umbral_min,par_red,0,False)
    return umbral

def D_js(P, Q):
    xM = 0.5*(P+Q)
    KLDpm = entropy(P, xM)
    KLDqm = entropy(Q, xM)
    JSDpq = 0.5*(KLDpm + KLDqm)
    return JSDpq

def data_Div(total_data,cv_time):
    data=total_data[:,0]
    id=total_data[:,1]

    # for i in range(50):#len(total_data)):
    #     print(data[i],"\t",id[i])

    n_cv = int(np.round((data[-1]-data[0])/cv_time))

    for t in range(0,n_cv):
        # print('janela: ', t)
        a = t*cv_time + data[0]
        b = a + cv_time

        fdata = open('Mar10/time_cluster_shank_sua_2020Mar10_'+str(t)+'.txt', 'w')
        try:
            for j in range(len(data)):
                if(data[j]>=a and data[j]<=b):
                    fdata.write(str(data[j])+"\t"+str(id[j])+'\n')
        finally:
            fdata.close()


def delete_cols(M):
    print(np.shape(M))
    for j in range(len(M[0,:])):
        for i in range(len(M[0,:])):
            if(any(M[:,int(i)])!=True):
                A=np.delete(M, i, axis=1)
                break

        M = np.resize(M,(len(A[:,0]),len(A[0,:])))
        M=A
    print(np.shape(M))
    return M


def data_cv_dt(data,cv_time,bin_cv,xdir,wt):

    tcb2=XNAT2(data)

    CV,mspike,mxid=tcb2.CV_calc(cv_time,bin_cv, show=False, save=False, ADFA=False,
                            showbin=False, bin_size=None)
    px=np.where(np.array(mspike,dtype=object))[0]

    xspk=np.array(mspike,dtype=object)
    idspk=np.array(mxid,dtype=object)

    mx=tcb2.CV

    cv_sr=np.array(sorted(mx.copy()))

    cv_sort=cv_sr[cv_sr<3]
    cop_cv_sort=np.unique(cv_sort.copy())

    zpy=np.where(np.array(mx,dtype=object))[0]

    nNwt=int(len(cv_sort)/wt)

    dtcv_sort=[]; iDcv_sort=[]

    for i in range(len(cop_cv_sort)):
        for j in range(len(mx)):
            if(cop_cv_sort[i] == mx[j]):
                dtcv_sort.append(zpy[j])
                iDcv_sort.append(zpy[j])

    for cvi in range(nNwt):
        a=wt*(cvi)
        b=wt*(cvi+1)
        cv_tiepo_wt=[];  cv_id_wt=[]
        for tk in range(a,b):
            cv_tiepo_wt += xspk[dtcv_sort[tk]]
            cv_id_wt += idspk[iDcv_sort[tk]]

    mx_sort=np.sort(mx.copy())
    cv_sort=mx_sort[mx_sort<3]
    average_CV=[]; cvs_agrupados=[]

    for cvi in range(nNwt):
        a=wt*(cvi)
        b=wt*(cvi+1)
        avercv=[]
        for tk in range(a,b):
            avercv.append(cv_sort[tk])

        average_CV.append(np.mean(avercv))
        cvs_agrupados.append(avercv)

        '''
        ordTime=sorted(np.unique(cv_tiepo_wt))

        try:
            fdata = open(xdir+str(cvi)+'.txt', 'w')

            for orx in ordTime:
                for ory in range(len(cv_tiepo_wt)):
                    if(orx==cv_tiepo_wt[ory]):
                        fdata.write(str(cv_tiepo_wt[ory])+"\t"+str(cv_id_wt[ory])+'\n')
        finally:
            fdata.close()
        '''

    return average_CV, cvs_agrupados

def polygon_under_graph(x, y):
    y




def gafplot(x,y,zx):
    zx.plot(x,y,"o",color="#7D119B")
    return x,y

def config_axes_canvas(ko,axhx,dy0,zrtx,nam_experiment,cx):

    axhx.set_title(nam_experiment,fontsize="x-large")

    #axhx.legend(loc = "upper right",fontsize='x-large')

    label_ylist = r'${:.2f}$'
    label_xlist = r'${:.0f}$'

    epsilon=np.array([cx,cx,cx])
    ylist=np.linspace(min(dy0)-epsilon[ko,1],max(dy0)+epsilon[ko,1],7)
    xlist=np.linspace(min(zrtx)-epsilon[ko,0],max(zrtx)+epsilon[ko,0],7)

    rxm=np.array([[r"Variation coefficient",r"Time $(s)$"],
            [r"Variation coefficient",r"Time $(s)$"],
            [r"Variation coefficient",r"Time $(s)$"]])

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

def camvax(axhx,label_format,label_list,titles,axes_limit=np.array([[0,0],[0,0]]),limit=False):

    axhx.xaxis.set_major_locator(mticker.FixedLocator(label_list[0]))
    axhx.set_xticklabels([label_format[0].format(x) for x in label_list[0]],fontsize=15)
    axhx.yaxis.set_major_locator(mticker.FixedLocator(label_list[1]))
    axhx.set_yticklabels([label_format[1].format(x) for x in label_list[1]],fontsize=15)

    axhx.set_ylabel(titles[0],fontsize=20)
    axhx.yaxis.set_tick_params(labelsize=20)

    axhx.set_xlabel(titles[1],fontsize=20)
    axhx.xaxis.set_tick_params(labelsize=20)

    if(limit):
        axhx.set_xlim(axes_limit[0,0],axes_limit[0,1])
        axhx.set_ylim(axes_limit[1,0],axes_limit[1,1])

def ant():

    dirExp=["Mar07","Mar10","Mar04","Jan29","Jan21"]#,"Jan14"]
    for i in dirExp:
        xdm=np.genfromtxt("datos_spikes/"+i+".txt")
        plt.plot(xdm[:,0],xdm[:,3],"o")
        plt.show()

def mbarcax(cv,pargraf,cmap,label_format,ticks_loc,xcvd):
    norm= matplotlib.colors.Normalize(vmin=min(cv), vmax=max(cv))
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    ax2_divider = make_axes_locatable(pargraf)
    cax0 = ax2_divider.append_axes("bottom", size="4%", pad="20%")
    cbar0=fig.colorbar(sm,extend='both' ,cax=cax0,orientation='horizontal')
    cbar0.ax.set_xlabel(xcvd,fontsize=20)
    cbar0.ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar0.ax.set_xticklabels([label_format.format(x) for x in ticks_loc],fontsize=13)
    cax0.xaxis.set_ticks_position("bottom")

def mbar(ccv,cmap,label_format,ticks_loc,xcvd):
    #dimensions [left, bottom, width,height]
    norm= matplotlib.colors.Normalize(vmin=min(ccv), vmax=max(ccv))
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    cax0 = plt.axes([0.25, 0.07, 0.5, 0.015])
    cbar0=fig.colorbar(sm,extend='both' ,cax=cax0,orientation='horizontal')
    cbar0.ax.set_xlabel(xcvd,fontsize=20)
    cbar0.ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    cbar0.ax.set_xticklabels([label_format.format(x) for x in ticks_loc],fontsize=15)
    cax0.xaxis.set_ticks_position("bottom")

'''
fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(2,5)
nNo=[0,1,2,0,1,2];r=0
ct=[0,0.2]
#for expr in ["Exp_Mar07","Exp_Mar10","Exp_Mar04","Exp_Jan29","Exp_Jan21","Exp_Jan14"]:
def pp():

    for i,expr in enumerate(["Exp_Mar07","Exp_Mar10","Exp_Mar04","Exp_Jan29","Exp_Jan21","Exp_Jan14"]):

        dirExp=["Mar07","Mar10","Mar04","Jan29","Jan21","Jan14"]
        #dirExp=["Mar04"]

        nome_exper=["SUAM07","SUAM10","SUAM04","SUAJ29","SUAJ21","SUAJ14"]
        #nome_exper=["SUAM04"]

        rdir=expr+"/time_cluster_shank_sua_2020"+expr+"_"

        if(os.path.isdir(expr)==False):
            path = Path(expr)
            path.mkdir(parents=True)

        if(dirExp[i]=="Mar04" or dirExp[i]=="Mar07"):
            data=np.loadtxt("datas/time_cluster_shank_sua_2020"+dirExp[i]+".txt")
        else:
            data=np.transpose(np.loadtxt("datas/time_cluster_shank_sua_2020"+dirExp[i]+".txt"))#np.transpose()

        xwt=5; xbin_cv=0.05; xcv_time=10

        txb=XNAT2(data)

        if(i>2):
            r=1
        o=nNo[i]

        ax0 = fig.add_subplot(gs[:,0])
        ax1 = fig.add_subplot(gs[r,o])

        azh=[ax0,ax1]

        #xx,yy=

        txb.similitud_coseno(xcv_time,xbin_cv,nome_exper[i],i,azh)#,azh)

        #limit_X,limit_Y=gafplot(xx[:-25],yy[:-25],azh)

        #config_axes_canvas(o,azh,limit_Y,limit_X,nome_exper[i],ct)
        #configure_canvas()
        #data_cv_dt(data,xcv_time,xbin_cv,rdir,xwt)

        # dir_coplet='../datos_cv/Exp_30s/'+expr
        # #extract_cluster_id(dir_coplet,expr)
        # extract_experimental_set(dir_coplet,expr)
    # plt.yscale("log")
    # plt.xscale("log")
    #plt.tight_layout()
    #plt.show()



    # ant()
'''


def dat_grup(dat_csv):
    with open(dat_csv) as fichero:
        df = re.split("\n",fichero.read())# ' '.join(map(str, df)))
    xdat=[]; tsr=[]
    for u in df:
        fdat=[]
        for at in re.split("\t",u):
            if at != "":
                fdat.append(int(at))

        xdat.append(fdat)
        tsr.extend(fdat)

    fDatx=np.zeros((len(tsr),2)); k=0

    for i,xlist in enumerate(xdat):
        for j in xlist:
            fDatx[k,0]=i
            fDatx[k,1]=j-1
            k+=1
    return fDatx

# def graf3D(x,y,vx):

#     fig = plt.figure()
#     ax = fig.add_subplot(projection='3d')

#     for i,j in enumerate(vx):
#         ax.plot(x, y[i,:], zs=j, zdir='x', label='curve in (x,y)')

#     plt.show()
'''
def qq():

    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2,3)
    nNo=[0,1,2,0,1,2];r=0
    ct=[0,0.2]
    expr="Exp_Mar04"

    dirExp=["Mar04"]

    rdir=expr+"/time_cluster_shank_sua_2020"+expr+"_"

    if(os.path.isdir(expr)==False):
        path = Path(expr)
        path.mkdir(parents=True)

    if(dirExp[0]=="Mar04" or dirExp[0]=="Mar07"):
        data=np.loadtxt("datas/time_cluster_shank_sua_2020"+dirExp[0]+".txt")
    else:
        data=np.transpose(np.loadtxt("datas/time_cluster_shank_sua_2020"+dirExp[0]+".txt"))#np.transpose()

    xwt=25; xbin_cv=0.05; xcv_time=10

    txb=XNAT2(data)

    # params = {'xtick.labelsize': 20, 'ytick.labelsize': 15}
    # mpl.rcParams.update(params)

    t,cv=txb.coef_variacion(xcv_time,xbin_cv,dirExp[0])

    label_ylist = r'${:.2f}$'
    label_xlist = r'${:.0f}$'
    ylist=np.linspace(min(cv),max(cv),10)
    xlist=np.linspace(min(t),max(t),10)
    titulos_ejes=["Coefficient of variation ","Time "+r"$(s)$"]
    lab_format=[label_xlist,label_ylist]
    lab_list=[xlist,ylist]

    xmicol1="plasma"
    xcmap = mpl.colormaps[xmicol1].resampled(250)
    xnewcolors = xcmap(np.linspace(0,1,250))
    xnewcmp = ListedColormap(xnewcolors)

    colors=["#0000FF","#008000","#BFBF00","#BF00BF","#FF0000"]

    ax0 = fig.add_subplot(gs[0,:])
    for j,nit in enumerate([0,10,19,28,43]):
        mc=abs(delete_nod(np.loadtxt("mc/mcc_Exp_Mar04_"+str(nit)+".txt")))
        mcf=np.zeros((len(mc[:,0]),len(mc[:,0])))
        io,jo=np.where(mc<0.35)
        mcf[io,jo]=mc[io,jo]
        # ax1 = fig.add_subplot(gs[1,j])
        # ax1.imshow(mcf,cmap=xnewcmp)
        # ax1.set_title(r"$\bullet$",color=colors[j],fontweight ="bold",fontsize="xx-large")

        # camvax(ax1,[r'${:.0f}$',r'${:.0f}$'],
        #             [np.linspace(0,len(mcf[:,0]),9),np.linspace(0,len(mcf[:,0]),9)],
        #                 ["neuron "+r"$i$","neuron "+r"$j$"],np.array([[0,len(mcf[:,0])-0.5],
        #                                                             [0,len(mcf[:,0])-0.5]]),False)

    #mbar([0,np.max(mcf)],xnewcmp,r'${:.2f}$',np.linspace(0,np.max(mcf),15),"functional weight")
    #t,cv=txb.coef_variacion(xcv_time,xbin_cv,dirExp[0])
    #ax0.plot(t,cv,"o",color="black")

    #fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    #ax.plot(t,cv,"o",color="black")

    dat_csv="SUAM04.txt"
    mx= dat_grup(dat_csv)
    xo=mx[:,0]; yo=mx[:,1]

    #xo,yo=np.loadtxt(dat_csv, usecols=(0,1), unpack=True)

    yo = np.array([int(q) for q in yo])

    for o in np.unique(xo):
        io=np.where(xo==int(o))
        cv[yo[io]]
        ax0.plot(t[yo[io]],cv[yo[io]],"o",color=colors[int(o)])

    # for m,u in enumerate([0,10,19,28,43]):
    #     dtim=[]
    #     for i,p in enumerate(cV[u]):
    #         io=np.where(cv==p)[0]
    #         dtim.extend(t[io].tolist())
    #     ax0.plot(dtim,cV[u],"o",color=colors[m])



    # mean_cv,cV=data_cv_dt(data,xcv_time,xbin_cv,rdir,xwt)
    # for m,u in enumerate([0,10,19,28,43]):
    #     dtim=[]
    #     for i,p in enumerate(cV[u]):
    #         io=np.where(cv==p)[0]
    #         dtim.extend(t[io].tolist())
    #     ax0.plot(dtim,cV[u],"o",color=colors[m])

    camvax(ax0,lab_format,lab_list,titulos_ejes)

    #config_axes_canvas(o,azh,limit_Y,limit_X,nome_exper[i],ct)
    #configure_canvas()

    # # dir_coplet='../datos_cv/Exp_30s/'+expr
    # # #extract_cluster_id(dir_coplet,expr)
    # # extract_experimental_set(dir_coplet,expr)
    # # plt.yscale("log")
    # # plt.xscale("log")
    plt.tight_layout()
    plt.show()
    # # ant()
'''

def  qgrafh_plot(i,u,gs):

    expr=np.array([["Exp_Mar07","Exp_Mar10","Exp_Mar04"],
                  ["Exp_Jan29","Exp_Jan21","Exp_Jan14"]])

    dirExp=np.array([["Mar07","Mar10","Mar04"],
                    ["Jan29","Jan21","Jan14"]])

    if(os.path.isdir(expr[i,u])==False):
        path = Path(expr[i,u])
        path.mkdir(parents=True)

    if(dirExp[i,u]=="Mar04" or dirExp[i,u]=="Mar07"):
        data=np.loadtxt("datas/time_cluster_shank_sua_2020"+dirExp[i,u]+".txt")
    else:
        data=np.transpose(np.loadtxt("datas/time_cluster_shank_sua_2020"+dirExp[i,u]+".txt"))

    xwt=25; xbin_cv=0.05; xcv_time=10

    txb=XNAT2(data)

    dat_csv=np.array([["SUAM07","SUAM10","SUAM04"],
                    ["SUAJ29","SUAJ21","SUAJ14"]])


    t,cv=txb.coef_variacion(xcv_time,xbin_cv,dirExp[i,u])

    # label_ylist = r'${:.2f}$'
    # label_xlist = r'${:.0f}$'
    # ylist=np.linspace(min(cv),max(cv),8)
    # xlist=np.linspace(min(t),max(t),6)
    # titulos_ejes=["Coefficient of variation ","Time "+r"$(s)$"]
    # lab_format=[label_xlist,label_ylist]
    # lab_list=[xlist,ylist]
    ax0 = fig.add_subplot(gs[i,u])
    
    rx=txb.similitud_coseno(xcv_time,xbin_cv,dat_csv[i],i,ax0)
    x=np.linspace(0,len(rx[0,:]),len(rx[0,:]))

    label_ylist = r'${:.0f}$'
    label_xlist = r'${:.0f}$'
    ylist=np.linspace(0,np.max(rx),8)
    xlist=np.linspace(0,max(x),6)
    titulos_ejes=[r"$\rho(\eta)$",r"$\eta$"]
    lab_format=[label_xlist,label_ylist]
    lab_list=[xlist,ylist]


    xmicol1="plasma"
    xcmap = mpl.colormaps[xmicol1].resampled(250)
    xnewcolors = xcmap(np.linspace(0,1,250))
    xnewcmp = ListedColormap(xnewcolors)

    colors=["#BA55D3","#800000","#FF4500","#4169E1","#FF0000","#FF8C00","#BC8F8F","#F4A460","#4682B4"]
    #colors=["#a3b4cc","#e6c100","#cc0000","#b451cc"]
    
    mx= dat_grup(dat_csv[i,u]+"Qsil.txt")

    xo=mx[:,0]; yo=mx[:,1]

    print(dat_csv[i,u]+".txt")
    #xo,yo=np.loadtxt(dat_csv[i,u]+"Py.txt", usecols=(0,1), unpack=True)

    yo = np.array([int(q) for q in yo])

    for o in np.unique(xo):
        io=np.where(xo==int(o))
        cv[yo[io]];p=[]
        for d in rx[yo[io],:]:
            d
            #ax0.plot(x,d,color=colors[int(o)])
        ax0.plot(t[yo[io]],cv[yo[io]],"o",color=colors[int(o)])


    # ax0.set_yscale("log")

    # ax0.set_yscale("log")
    # ax0.set_xscale("log")
    ax0.set_title(dat_csv[i,u],fontsize="x-large")
    
    camvax(ax0,lab_format,lab_list,titulos_ejes)

    #ax0 = fig.add_subplot(gs[i,u])
    #txb.similitud_coseno(xcv_time,xbin_cv,dat_csv[i,u],i,ax0)#,azh)

fig = plt.figure(figsize=(17, 10))
gs = fig.add_gridspec(2,3)

for i in [1]:#range(2):
    for u in [2]:#range(3):
        qgrafh_plot(i,u,gs)

plt.tight_layout()
fig.savefig('fhc.jpg')

plt.show()



# qq()



