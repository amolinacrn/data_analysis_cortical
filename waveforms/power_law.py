#!/usr/bin/env python
# -*- coding: utf-8 -*

from cProfile import label
from math import log
from operator import truediv
from turtle import color, fd
import xxlimited
import ROOT
from ROOT import gPad
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib
matplotlib.rcParams['text.usetex'] = True 
import matplotlib.ticker as mticker
import seaborn as sns
from scipy.stats import norm
from dist_eucl import graf_dist_fisica,fun_histo
from pathlib import Path
import os
from MAIN import miplot_main

def fun_print_parametros(zxprt,datafolder):

    mi_path="{}".format(datafolder)
    print(mi_path)
    
    if(os.path.isdir(mi_path)==False):
        path = Path(mi_path)
        path.mkdir(parents=True)#crear capetas
    
    fdata = open(mi_path+"/parametros.txt", 'w')

    try:
        #fdata.write("$b$\t" +"&\t" +"$m$\t"+ "&\t"+ "$& \Delta b$\t"+"& \t"+"& \Delta m"+"& \chi^2 "+"& $\ "+"nu$"+"\n")     
        for i in range(len(zxprt[:,0])):
            for j in range(len(zxprt[0,:])):  
                fdata.write(str(np.format_float_scientific(zxprt[i,j], precision = 1, exp_digits=2)))
                
                if(j<len(zxprt[0,:])-1):
                    fdata.write("\t"+"&\t")

            fdata.write("\\"+"\\"+"\t"+"\hline")     
                          
            fdata.write("\n")       
    finally:
        fdata.close()


def funcion_normal(x,par):
    a=par[0]
    b=par[1]
    # c=par[2]
    
    z=x[0]
    f = a*pow(z,b)
    
    return f

def graf_root(x,y,eyy,aa,bb,cc,x_min,x_max,num_par,name_ajs):
    c=ROOT.TCanvas("cV5","migrafico",300,500,550,470)
    xpar=[]; chisqrt=[]; nDf=[]
    mx = np.array(x,dtype = float)
    my = np.array(y,dtype = float)
    ex = np.array(len(x)*[0],dtype = float)
    ey = np.array(len(x)*[0],dtype = float)

    gr = ROOT.TGraphErrors(len(x),mx,my,ex,ey)
    gr.SetMinimum(1E-5)
    gr.SetMaximum(2E-1)
    gr.GetMarkerStyle(23)
    c.SetLogy()
    c.SetLogx()
    
    f3 = ROOT.TF1(name_ajs,funcion_normal,x_min,x_max,num_par)
    f3.SetParameters(aa,bb)
    gr.Fit(f3,"","",x_min,x_max)
    gr.Draw("");input("Prompt: ")
    a=f3.GetParameters()#c.Update()
    erp0=f3.GetParError(0)
    erp1=f3.GetParError(1)
    chisq=f3.GetChisquare()
    ndf=f3.GetNDF()

    for p in [0,1]:
        xpar.append(a[p])
    return xpar, erp0, erp1, chisq, ndf,  chisq/ndf


def funplots(ko,axhx,dx0,dy0,ery,xx,xfpar):#,zrtx,zrty):

    f = xfpar[0]*pow(xx,xfpar[1])
    #xdf = {'x': np.log(dx0), 'y': np.log(dy0)}

    ery=np.array(ery)#/np.array(dy0)
    dx0=(np.array(dx0))
    dy0=(np.array(dy0))
    
    miscorlores=["#7D119B","#7D119B","#7D119B","#7D119B"]
    
    #axhx.errorbar(dx0, dy0, ery, fmt="none", linewidth=0.5, capsize=1,color=miscorlores[ko],barsabovebool=False)
    axhx.plot(dx0,dy0,"v",color=miscorlores[ko])
    axhx.plot(xx,f,"-",color="black",label="fitting model")
    
    # label_ylist = r'${:.2f}$'
    # label_xlist = r'${:.1f}$'

    # epsilon=np.array([[0.0,1],[0.01,0.01],[0.1,0.1]])
    # ylist=np.linspace(min(dy0)-epsilon[ko,1],max(dy0)+epsilon[ko,1],7)

    rxm=np.array([[r"$\langle C \rangle$",r"$CV$"],
            [r"$\langle C \rangle$",r"$CV$"],
            [r"$\langle C \rangle$",r"$CV$"],
            [r"$\langle C \rangle$",r"$CV$"]])
    
    # axhx.yaxis.set_major_locator(mticker.FixedLocator(ylist))
    # axhx.set_yticklabels([label_ylist.format(x) for x in ylist],fontsize=14)

    axhx.set_ylabel(r""+str(rxm[ko,0])+"",fontsize=20)
    axhx.yaxis.set_tick_params(labelsize=20)

    axhx.set_xlabel(r""+str(rxm[ko,1])+"",fontsize=20)
    axhx.xaxis.set_tick_params(labelsize=20)
    # axhx.set_ylim(min(dy0)-epsilon[ko,1],max(dy0)+epsilon[ko,1])
    # axhx.set_xlim(min(dx0)-epsilon[ko,0],max(dx0)+epsilon[ko,0])
    axhx.set_yscale("log")
    axhx.set_xscale("log")

def gen_datos_fit():
    params = {'xtick.labelsize': 20, 'ytick.labelsize': 20}
    mpl.rcParams.update(params)
    sdats=["Jan29"]#,"Mar04","Mar10","Jan21","Jan14"]

    mdat=miplot_main()
  
    bfile_txt="datos_power_law.txt"
    bfdata = open(bfile_txt, 'w')

    try:
        for idats in sdats:

            simportar_datos=np.loadtxt("coordenadas_clusters/crsua"+str(idats)+".txt")
            simpDat=np.loadtxt("restW10T50/50s"+str(idats)+".txt")
            
            vx,vy,erry,ang=graf_dist_fisica(mdat,simpDat,simportar_datos,idats)
            
            #vx,vy=fun_histo(ang,3000)
            plt.plot(vx,vy,"-")
            
            for k in range(len(vx)):
                bfdata.write(str(vx[k])+"\t"+str(vy[k])+"\n") 
        
        x_max=4E-2; x_min=6E-4
        x = np.linspace(x_min,x_max,100)
        a=8.66E-07; b=-1.450
        fx=(a+0.000005)*pow(x,b)
        #plt.plot(x,fx,"--",linewidth=2,color="black")
        plt.yscale("log")
        plt.xscale("log")
        #plt.text(0.05, 0.04, r'$\rho\left(\Theta_{ij}\right)\propto\Theta_{ij}^{-\alpha}$', fontsize=20, color='blue')
        # plt.text(0.05, 0.2, r'$\Theta_{ij}=\tan^{-1}\left(\frac{w_{ij}}{d_{ij}}\right)$', fontsize=20, color='blue')
        # plt.text(0.05, 0.005, r'$\alpha \approx 1.40$', fontsize=20, color='blue')

        plt.xlabel(r'$\Theta_{ij}$',fontsize=35)
        plt.ylabel(r'$\rho(\Theta_{ij})$',fontsize=35)
        # plt.tick_params(labelsize = 20)
        plt.tight_layout()
        plt.show()

    finally:
        bfdata.close()

def graf_power_law():

    params = {'xtick.labelsize': 20, 'ytick.labelsize': 20}
    mpl.rcParams.update(params)
    datos=np.loadtxt("datos_power_law.txt")
    vx=datos[:,0]; vy=datos[:,1]
    erry=len(vx)*[0]
    x_max=4E-2; x_min=6E-4
    x = np.linspace(x_min,x_max,100)
    a=8.66E-07; b=-1.50
    fx=(a+0.000005)*pow(x,b)

    plt.plot(x,fx,"--",linewidth=2)
    plt.plot(vx,vy,"-")
    plt.yscale("log")
    plt.xscale("log")
    plt.text(0.05, 0.04, r'$\rho\left(\Theta\right)\propto\Theta^{-\alpha}$', fontsize=20, color='blue')
    plt.text(0.05, 0.2, r'$\Theta=\tan^{-1}\left(\frac{w_{ij}}{d_{ij}}\right)$', fontsize=20, color='blue')
    plt.text(0.05, 0.005, r'$\alpha \approx 1.40$', fontsize=20, color='blue')

    plt.xlabel(r'$\Theta_{ij}$',fontsize=35)
    plt.ylabel(r'$\rho(\Theta_{ij})$',fontsize=35)
    # plt.tick_params(labelsize = 20)
    plt.tight_layout()
    plt.show()

    #xmin=min(vx); xmax=max(vx)
    # x_max=3.5E-3
    # x_min=2.5E-4

    # xpr0, xpr1, xpr2, xpr3, xpr4, xpr5 = _=graf_root(vx,vy,erry,aax,bbx,ccx,x_min,x_max,num_par,"ajuste")
    # tx.append([xpr0[0],xpr0[1],xpr1 ,xpr2, xpr3, xpr4, xpr5])
    # zx = np.array(np.linspace(1E-4,x_max,100))   
    # funplots(oo,azh,vx,vy,erry,zx,xpr0)

    # # parats=np.array(tx)

    # # fun_print_parametros(parats,"datos_ajustes")

    # plt.tight_layout()
    # plt.show()

    # num_par=2

    # aax=1.5E-6; bbx=-1.2; ccx=1 ; xn=50

    # fig = plt.figure(figsize=(10, 5))
    # gs = fig.add_gridspec(1,4)

    # Xfig, rzx = plt.subplots(1, 1, figsize=(5,5 ))

gen_datos_fit()