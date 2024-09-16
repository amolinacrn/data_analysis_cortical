
from cProfile import label
from math import log
from operator import truediv
from turtle import color, fd
import xxlimited
import ROOT
from ROOT import gPad
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams['text.usetex'] = True 
import matplotlib.ticker as mticker
import seaborn as sns
from scipy.stats import norm
from dist_eucl import graf_dist_fisica,fun_histo
from pathlib import Path
import os
from MAIN import miplot_main


def funcion_normal(x,par):
    a=par[0]
    b=par[1]
    # c=par[2]
    
    z=x[0]
    f = a*pow(z,b)
    
    return f

def graf_root(x,y,eyy,aa,bb,x_min,x_max,num_par,name_ajs):
    c=ROOT.TCanvas("cV5","migrafico",300,500,550,470)
    xpar=[]; chisqrt=[]; nDf=[]
    mx = np.array(x,dtype = float)
    my = np.array(y,dtype = float)
    ex = np.array(len(x)*[0],dtype = float)
    ey = np.array(len(x)*[0],dtype = float)

    gr = ROOT.TGraphErrors(len(x),mx,my,ex,ey)
    gr.SetMarkerStyle(23)
    gr.SetMinimum(1E-5)
    gr.SetMaximum(1E-0)
    c.SetLogy()
    c.SetLogx()
    
    f3 = ROOT.TF1(name_ajs,funcion_normal,x_min,x_max,num_par)
    f3.SetParameters(aa,bb)
    gr.Fit(f3,"","",x_min,x_max)
    gr.Draw("AP");input("Prompt: ")
    a=f3.GetParameters()#c.Update()
    erp0=f3.GetParError(0)
    erp1=f3.GetParError(1)
    chisq=f3.GetChisquare()
    ndf=f3.GetNDF()

    for p in [0,1]:
        xpar.append(a[p])
    return xpar, erp0, erp1, chisq, ndf,  chisq/ndf

datos=np.loadtxt("datos_power_law.txt")
vx=datos[:,0]; vy=datos[:,1]
erry=len(vx)*[0]

num_par=2; aax=1.5E-6; bbx=-1.2
x_max=1E-2; x_min=1.5E-4

xpr0, xpr1, xpr2, xpr3, xpr4, xpr5 =graf_root(vx,vy,erry,aax,bbx,x_min,x_max,num_par,"ajuste")
