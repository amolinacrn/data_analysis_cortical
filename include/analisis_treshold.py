#!/usr/bin/env python
# -*- coding: utf-8 -*

##/////////////////////////////////////////////////////##
##     Por Miguel Alejandro Molina  03/11/2022         ##
##/////////////////////////////////////////////////////##

#------------------------------------------------------------------------------------------


import numpy as np
import os
import numpy
import math
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import seaborn as sns
import plotly.graph_objects as go
import networkx as nx
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors
import seaborn as sns
#import ROOT # paquete root/cern debe instalarse por separado
import statistics
from scipy.signal import savgol_filter
from scipy.interpolate import make_interp_spline
from scipy import stats


def thresholded_conetivity_matrix(xdir,exp,xmc,xbin,xw,file_txt):
    '''
    xdir: directorio principal de datos sustitutos
    exp: tipo de experimento
    xmc: matriz de conectividad empirica
    xbin: tamaño de bin 
    xw: venta de correlacion
    file_txt: nombre del archivo de la matriz de conectividad   
    '''
    
    mxc_surr= []
     
    for k in range(len(xmc)):
        for m in range(len(xmc)):  
            mxfile=str(xdir)+"/"+str(exp)+"/surr_"+str(1)+"/ncc_0.1-"+str(xbin)+"-"+str(xw)+"msec/surr_"+str(1)+"_CC_Symmetric.txt" 
            print(mxfile)
            
    '''        
            
            mcc_surr=np.genfromtxt(mxfile)     
            mxc_surr.append([mcc_surr[k,m]]) 

    treshold = statistics.mean(mxc_surr)+2*statistics.stdev(mxc_surr)
    
    for i in range(len(xmc)):

        for j in range(len(xmc)):
        
            if (xmc[i,j] > treshold):    
                xmc[i,j]=xmc[i,j]
            else:  
                xmc[i,j]=0
        
    try:
            
        fdata = open(file_txt, 'w')
        
        for wrfil in range(len(xmc)):
            for wrcol in range(len(xmc)):
                fdata.write(str(xmc[wrfil][wrcol]))
                fdata.write("\t")
            fdata.write("\n")  
    
    finally:
        fdata.close()    
    '''