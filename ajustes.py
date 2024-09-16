#!/usr/bin/env python
# -*- coding: utf-8 -*

from cProfile import label
from turtle import color
import ROOT
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams['text.usetex'] = True 
import matplotlib.ticker as mticker
from MAIN import *


par_x,par_y, par_z,par_v,par_w=miplot_main()

from scipy.stats import norm

def funcion_normal(x,par):
    const=par[0]
    sg=par[1]
    mu=par[2]
    z=x[0]
    f = const/(sg*np.sqrt(2*np.pi))*np.exp(-(z-mu)**2/(2*sg**2))
    return f

def graf_root(xdat,const,dvstr,media,x_min,x_max,num_par,name_ajs):
    #c=ROOT.TCanvas()
    xpar=[]; erpar=[]
    mx = np.array(xdat[:,1],dtype = float)
    my = np.array(xdat[:,0],dtype = float)
    ex = np.array(xdat[:,3],dtype = float)
    ey = np.array(xdat[:,2],dtype = float)

    gr = ROOT.TGraphErrors(len(xdat),mx,my,ex,ey)

    f3 = ROOT.TF1(name_ajs,funcion_normal,x_min,x_max,num_par)
    f3.SetParameters(const,dvstr,media)
    gr.Fit(f3,"","",x_min,x_max)
    #gr.Draw("AP*");input("Prompt: ")
    a=f3.GetParameters();#c.Update()
    err=f3.GetParError()
    for p in [0,1,2]:
        xpar.append(a[p])
        erpar.append(err[p])
    return xpar,erpar
    
               
def funct_lognormal(addres_dir,nsubdir,file_dir,mdat,x_data,p_min,p_max,axy2,x_color,xcomu,y_max):
    xbins=80
    xlin=0.5
    nbins=np.linspace(p_min,p_max,xbins+1)

    mxrango = slice(14,22)

    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    thrho=np.mean(mdat[:,4])-2*np.std(mdat[:,4])# minimo de Desdidad
   
    nNcbajo=0
  
    for i in range(len(x_data)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(x_data[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        if(rsz>thrho):
            nNcbajo+=1
    
    xcv_frec_bajo=np.zeros((int(xbins),nNcbajo)); xyfr_bajo=[]
    dst_bajo=[]

    ox=0

    for i in range(len(x_data)):
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(x_data[i])+'.txt')))
        rsz=nx.density(nx.from_numpy_array(xmcc))
        
        if(rsz>thrho):
          
            xrfrec, xrxbin=np.histogram(np.log(xmcc[xmcc>0]),bins=nbins)
           
            zrx=xrxbin[:-1]; zry= xrfrec/sum(xrfrec)
            for yfr in range(len(zry)):
                xcv_frec_bajo[yfr,ox] = zry[yfr]

            axy2.plot(zrx,zry,linewidth=xlin,color="#bababa")
            ox += 1

    for k in range(len(xcv_frec_bajo[:,0])):
        xyfr_bajo.append(np.mean(xcv_frec_bajo[k,:]))

    for k in range(len(xcv_frec_bajo[:,0])):
        dst_bajo.append(stdev(xcv_frec_bajo[k,:]))

    axy2.plot(zrx[0],zry[0],linewidth=3,color="#bababa",label="Experimental data")
    axy2.plot(zrx,xyfr_bajo,"--",linewidth=3,color=x_color,label="Average data")
    axy2.legend(loc = "upper right",ncol=2,fontsize='large')

    label_ylist = r'${:.2f}$'
    label_xlist = r'${:.1f}$'

    xlist=np.linspace(p_min,p_max,5)
    ylist=np.linspace(0,y_max,5)

    axy2.xaxis.set_major_locator(mticker.FixedLocator(xlist))
    axy2.set_xticklabels([label_xlist.format(x) for x in xlist],fontsize=14)

    axy2.yaxis.set_major_locator(mticker.FixedLocator(ylist))
    axy2.set_yticklabels([label_ylist.format(x) for x in ylist],fontsize=14)

    axy2.set_ylim(0,y_max+0.005)

    axy2.set_xlim(p_min-0.1,p_max+0.1)

    axy2.xaxis.set_tick_params(labelsize=15)
    axy2.yaxis.set_tick_params(labelsize=15)
    axy2.set_xlabel(r"$\Omega(\mathcal{W})$",fontsize=20)
    axy2.set_ylabel(r"$\rho(\Omega)$",fontsize=20)

    #axy2.legend(loc = "upper right",ncols=2,fontsize='large')
    fig.tight_layout()


    bfile_txt=nsubdir+xcomu+".txt"
    bfdata = open(bfile_txt, 'w')
    
    try:
        for i in range(len(xyfr_bajo)):
            bfdata.write(str(xyfr_bajo[i])+"\t"+str(zrx[i])+"\t"+str(dst_bajo[i])+"\t"+str(0)+"\n")      
    finally:
        bfdata.close()


def funplots(xdat,xx,fpar,ax0,xcol,x_min,x_max,yx_max):
    f = fpar[0]*norm.pdf(xx, fpar[2], fpar[1])
    mean, var, skew, kurt = norm.stats(moments='mvsk')
    rx = np.array(xdat[:,1],dtype = float)
    ry = np.array(xdat[:,0],dtype = float)
    erx = np.array(xdat[:,3],dtype = float)
    ery = np.array(xdat[:,2],dtype = float)
     
    plt.plot(rx,ry,"--",linewidth=3,color=str(xcol),label="average data") 
    plt.plot(xx,f,"-",color="black",label="fitting model")
    plt.fill_between(rx, ry-ery, ry+ery,color="gray",alpha=0.2)
    
    label_format = r'${:.2f}$'
    label_format_xlist = r'${:.1f}$'
    xlist=np.linspace(x_min,x_max,5)
    ylist=np.linspace(0,yx_max,5)

    ax0.xaxis.set_major_locator(mticker.FixedLocator(xlist))
    ax0.set_xticklabels([label_format_xlist.format(x) for x in xlist],fontsize=14)

    ax0.yaxis.set_major_locator(mticker.FixedLocator(ylist))
    ax0.set_yticklabels([label_format.format(x) for x in ylist],fontsize=14)

    ax0.set_ylim(0,yx_max+0.005)
    ax0.set_xlim(x_min-0.1,x_max+0.1)
    ax0.set_ylabel(r"$\rho(\Omega)$",fontsize=20)
    ax0.yaxis.set_tick_params(labelsize=17)

    ax0.set_xlabel(r"$\Omega(\mathcal{W})$",fontsize=20)
    ax0.xaxis.set_tick_params(labelsize=17)
    ax0.legend(loc = "upper right",ncol=2,fontsize='large')
    fig.tight_layout()
    
#intervalo de busqueda ajuste

x_min=[-3.2,-3.2,-3.2]; x_max=[-1.6,-0.90,-0.90]; num_par=3 
#parametros de inicializacion
const=0.019; dvstr=0.3; media=-2.5
#datos
fig = plt.figure(figsize=(12, 4))
gs = fig.add_gridspec(2, 3)
xcolr=["#9c46b1","#d1af00","#ac0000"]
y_max=0.06;xy_max=0.05

for i,xcom in enumerate(["bajo","medio","alto"]):
    ax0 = fig.add_subplot(gs[0,i])
    ax1 = fig.add_subplot(gs[1,i])

    funct_lognormal(par_x,par_y, par_z,par_v,par_w[i],x_min[i],x_max[i],ax0,xcolr[i],xcom,y_max)

    fdat = np.genfromtxt("50sMar07"+str(xcom)+".txt")
    xpr=graf_root(fdat,const,dvstr,media,x_min[i],x_max[i],num_par,"ajuste")
    zx = np.linspace(x_min[i],x_max[i],100)   
    funplots(fdat,zx,xpr,ax1,xcolr[i],x_min[i],x_max[i],xy_max)

plt.show()
