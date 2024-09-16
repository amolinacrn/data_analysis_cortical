#!/usr/bin/env python
# -*- coding: utf-8 -*

from cProfile import label
from math import log
from operator import truediv
from turtle import color, fd
import xxlimited
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
matplotlib.rcParams['text.usetex'] = True 
import matplotlib.ticker as mticker
import seaborn as sns
from scipy.stats import norm
from dist_eucl import graf_dist_fisica,fun_histo,smil_dist_fisica
from pathlib import Path
import os
from MAIN import miplot_main

def gen_datos_fit():
    params = {'xtick.labelsize': 20, 'ytick.labelsize': 20}
    mpl.rcParams.update(params)
    sdats=["Jan29","Mar04","Mar10","Jan21","Jan14"]

    #mdat=miplot_main()
  
    bfile_txt="datos_power_law.txt"
    bfdata = open(bfile_txt, 'w')


    for idats in sdats:

        simportar_datos=np.loadtxt("coordenadas_clusters/crsua"+str(idats)+".txt")
        simpDat=np.loadtxt("restW10T250/250s"+str(idats)+".txt")
        
        #cv,vx,vy,erry,ang=
        smil_dist_fisica(simpDat,simportar_datos,idats)
        #plt.hist2d(simpDat[:,12],simpDat[:,10],bins = [50,50])

        #plt.plot(cv,vx,"o")
        # plt.yscale("log")

        #plt.title(idats,fontsize="x-large")
        plt.show()

gen_datos_fit()