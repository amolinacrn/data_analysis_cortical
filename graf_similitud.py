#!/usr/bin/env python
# -*- coding: utf-8 -*

##/////////////////////////////////////////////////////##
##     Por Miguel Alejandro Molina  03/11/2022         ##
##/////////////////////////////////////////////////////##

#------------------------------------------------------------------------------------------
# importar achivos necesarios para el analisis

import os
import numpy as np
import matplotlib.pyplot as plt

def gafplot(x,y,colr):
    
    plt.plot(x,y,"o",color=colr)
dcv=10; col=3
colors=["#0000FF","#008000","#BFBF00","#BF00BF","#FF0000","#00BFBF"]
dirExp=["Jan14","Mar07","Mar10","Mar04",'Jan29',"Jan21"]

for i,fil in enumerate(dirExp):
    dat=np.loadtxt('datos_spikes/'+fil+'D.txt')
    xcv=np.array(sorted(dat[:,0]))
    datos=[]
    
    for l in range(len(xcv)):
        for k in range(len(xcv)):
            if(xcv[l]==dat[:,0][k]):
                datos.append(dat[k,:])
    datos=np.array(datos) 
    
    io=np.where(np.array(sorted(dat[:,0])))[0]
    
    nN=int(len(dat[:,0])/dcv)
    cv=[]; var_int=[]
    xdat=datos[:,col]
    
    for t in range(nN):
        a = t*dcv
        b = a + dcv
        jo= io[(io>=a) & (io<b)]
        cv.append(np.mean(xcv[jo]))
        var_int.append(np.mean(xdat[jo]))
                
    gafplot(cv,var_int,colors[i])
plt.show()
