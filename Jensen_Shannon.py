#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# example_use.py
# Jim Bagrow
# Last Modified: 2018-04-22

import sys, os
import itertools
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from bct.algorithms import number_of_components
from bct.utils import BCTParamError,get_rng
from scipy.stats import entropy
from portrait_Djs.funct_PDjs import arbole_expansion_min
from libreriaMolCr import delete_nod

def D_js(XG, XH,xbins):
    zP = XG[XG>0]
    zQ = XH[XH>0]
   
    pfr, pbin=np.histogram(zP,bins=xbins)
    P= pfr/sum(pfr)
    qfr, qbin=np.histogram(zQ,bins=xbins)
    Q= qfr/sum(qfr)

    M = 0.5*(P+Q)
    KLDpm = entropy(P, M)
    KLDqm = entropy(Q, M)
    JSDpq = 0.5*(KLDpm + KLDqm)
    
    return JSDpq

def jensen_shannon(dt,file_dir,addres_dir,xname,xbins,xfile,thsrho=0,xnz=0,denvar=False):

    mxrango = slice(14,22)

    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    MxRet=np.zeros((len(dt),len(dt)))
    
    if(denvar):
        for r in range(len(dt)-1):
            datx=abs(delete_nod(np.genfromtxt(fildir+str(int(dt[r]))+'.txt')))
            xzx=nx.from_numpy_array(datx)
            dzx=nx.density(xzx)
            if(dzx > thsrho):
                #print(dt[r])
                rhoxG=arbole_expansion_min(datx,thsrho,xnz)              
                print(fildir+str(int(dt[r]))+'.txt')
                for s in range(r+1,len(dt)): 
                    daty=abs(delete_nod(np.genfromtxt(fildir+str(int(dt[s]))+'.txt')))
                    yzj=nx.from_numpy_array(daty)
                    drx=nx.density(yzj)
                    if(drx > thsrho):
                        rhoyG=arbole_expansion_min(daty,thsrho,xnz)
                        xrest=D_js(rhoxG,rhoyG,xbins)
                        MxRet[r,s]=xrest; MxRet[s,r]=xrest
                        print(int(dt[s]))
        if(xfile):
            fxile_txt="mxsimilt_w10n50/"+xname+"w10n50.txt"
            fxdata = open(fxile_txt, 'w')
            
            try:
                for r in range(len(dt)):
                    for s in range(len(dt)):
                        fxdata.write(str(MxRet[r,s])+"\t")    
                    fxdata.write("\n")    
            finally:
                fxdata.close()  
    else:
        for r in range(len(dt)-1):
            datx=abs(delete_nod(np.genfromtxt(fildir+str(int(dt[r]))+'.txt')))
            #xzx=nx.from_numpy_array(datx)
            #dzx=nx.density(xzx)
            #if(dzx > thsrho):
            # fx, xb=np.histogram(datx[datx>0],bins=xbins)
            # sx=xb[:-1]; sy= fx/sum(fx)
            # plt.plot(sx,sy)
            
            print(fildir+str(int(dt[r]))+'.txt')
            
            for s in range(r+1,len(dt)): 
                daty=abs(delete_nod(np.genfromtxt(fildir+str(int(dt[s]))+'.txt')))
                #yzj=nx.from_numpy_array(daty)
                #drx=nx.density(yzj)
                #if(drx > thsrho):
                xrest=D_js(datx,daty,xbins)
                MxRet[r,s]=xrest; MxRet[s,r]=xrest

        # plt.xscale("log")
        # plt.yscale("log")
        # plt.show()
                    
        if(xfile):
            fxile_txt="mxsimilt_w10n50/"+xname+"w10n50.txt"
            fxdata = open(fxile_txt, 'w')
            
            try:
                for r in range(len(dt)):
                    for s in range(len(dt)):
                        fxdata.write(str(MxRet[r,s])+"\t")    
                    fxdata.write("\n")    
            finally:
                fxdata.close()  
              
   

                    
