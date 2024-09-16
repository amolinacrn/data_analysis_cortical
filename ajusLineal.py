#!/usr/bin/env python
# -*- coding: utf-8 -*

from cProfile import label
from operator import truediv
from os import PRIO_PGRP
from turtle import color, fd
import ROOT
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams['text.usetex'] = True 
import matplotlib.ticker as mticker
import seaborn as sns
from scipy.stats import norm


def fun_print_parametros(zxprt,name_exp):
    
    file_txt="par_ajust/parametros.txt"
    fdata = open(file_txt, 'w')

    try:
        #fdata.write(str(np.format_float_scientific(zxprt[i,j], precision = 1, exp_digits=2)))
        #fdata.write("$b$\t" +"&\t" +"$m$\t"+ "&\t"+ "$& \Delta b$\t"+"& \t"+"& \Delta m"+"& \chi^2 "+"& $\ "+"nu$"+"\n")     
        for i in range(len(zxprt[:,0])):
            fdata.write(name_exp[i]+"\t"+"&"+"\t")
            for j in range(len(zxprt[0,:])):  
                fdata.write(str(np.round(zxprt[i,j]+0.0003,3)))
                #fdata.write(str(np.format_float_scientific(zxprt[i,j], precision = 1, exp_digits=2)))
                
                if(j<len(zxprt[0,:])-1 and j>2):
                    fdata.write("\t"+"&\t")

                elif(j==0):
                    fdata.write("\t"+"$\pm$ \t")

                elif(j==1 ):
                    fdata.write("\t"+"&\t")

                
                elif(j==2 ):
                    fdata.write("\t"+"$\pm$\t")

            fdata.write("\\"+"\\"+"\t"+"\hline")     
                          
            fdata.write("\n")       
    finally:
        fdata.close()


def funcion_normal(x,par):
    a=par[0]
    b=par[1]
    z=x[0]
    f =a+b*z#a*np.exp(b*z)
    #a+b*z
    return f

def graf_root(x,y,eyy,aa,bb,x_min,x_max,num_par,name_ajs):
    xpar=[]; chisqrt=[]; nDf=[]
    mx = np.array(x,dtype = float)
    my = np.array(y,dtype = float)
    ex = np.array(len(y)*[0],dtype = float)
    ey = np.array(eyy,dtype = float)

    gr = ROOT.TGraphErrors(len(x),mx,my,ex,ey)

    f3 = ROOT.TF1(name_ajs,funcion_normal,x_min,x_max,num_par)
    f3.SetParameters(aa,bb)
    gr.Fit(f3,"","",x_min,x_max)
    #gr.Draw("AP*");input("Prompt: ")
    a=f3.GetParameters()#c.Update()
    erp0=f3.GetParError(0)
    erp1=f3.GetParError(1)
    chisq=f3.GetChisquare()
    ndf=f3.GetNDF()

    for p in [0,1]:
        xpar.append(a[p])

    return xpar, erp0, erp1, chisq, ndf,  chisq/ndf
def intercepto(xp):
    po=xp[0,0]; p1=xp[1,0]; dpo=xp[0,2]; dp1=xp[1,2]
    mo=xp[0,1]; m1=xp[1,1]; dmo=xp[0,3]; dm1=xp[1,3]
    dp=np.sqrt(dpo*dpo+dp1*dp1)
    dm=np.sqrt(dmo*dmo+dm1*dm1)
    xo=(po-p1)/(m1-mo)
    er_xo=xo*np.sqrt(pow(dp/(po-p1),2)+pow(dm/(m1-mo),2))
    return xo,er_xo

def plot_datas(ko,axhx,dx0,dy0,ery,cx,cy):
    miscorlores=["#7D119B","#7D119B","#7D119B"]
    axhx.plot(cx,cy,"o",color="gray")
    axhx.plot(dx0,dy0,"o",color=miscorlores[ko],label="Model fitted data.")
    #axhx.errorbar(dx0, dy0, ery, fmt="none", linewidth=0.5, capsize=1,color=miscorlores[ko])#,barsabovebool=False)
    #dz=dy0.copy()
    #i=np.where(dz.copy())[0]
    #dz[i]=0.6
    #axhx.plot(dx0,dz,"0",color="black",label="Experimental data")
    #sns.scatterplot(data = xdf,x = "x", y = "y",sizes=(200,200),marker = "v",color=miscorlores[ko],ax=axhx,legend=False)

def funct_model_completa(xp,x):
    po=xp[0]
    mo=xp[1]
    f=po+mo*(x)

    return x,f

def funct_modelo(xp,xo,x1):
    po=xp[0,0]; p1=xp[1,0]
    mo=xp[0,1]; m1=xp[1,1]

    fo=po+mo*(xo)
    f1=p1+m1*(x1)

    return xo,fo,x1,f1

def plot_model(axhx,dafunct):
    axhx.plot(dafunct[0],dafunct[1],"-",color="black",label=r"\textup{Fitting model}")


def plot_model_trozos(axhx,dafunct):
    axhx.plot(dafunct[0],dafunct[1],"-",color="black",label=r"\textup{Fitting model}")
    axhx.plot(dafunct[2],dafunct[3],"-",color="black")

def configure_axes_canvas(so,ko,axhx,dx0,dy0,ery,zrtx,nam_experiment,parametrs,cx,inpto=0,err_pto=0):
    
    axhx.set_title(nam_experiment,fontsize="x-large")

    # axhx.fill_between(dx0, dy0-ery, dy0+ery,color="gray",alpha=0.2,hatch="+") 
    # axhx.legend(loc = "upper left",fontsize='x-large')
    text4=0.05; text5=0.2
    # axhx.text(text4, text5, r'$\langle CV \rangle_{0}=$'+str(np.round(inpto,1))+r"$\pm$"+str(np.round(err_pto,1)), fontsize=15, color='red',transform = axhx.transAxes)
   

    if(so==3 or so==5 or so==4):
        axhx.legend(loc = "upper left",fontsize='x-large')
        text0=0.63; text2=0.05; text4=0.05
        text1=0.2; text3=0.7;  text5=0.6

        text0=0.63; text2=0.05; text4=0.05
        text1=0.2; text3=0.60;  text5=0.50
    
    # elif(so==5):
    #     axhx.legend(loc = "upper right",fontsize='x-large')
    #     text0=0.7; text2=0.03; text4=0.03
    #     text1=0.35; text3=0.8;  text5=0.7

    # elif(so==4):
    #     axhx.legend(loc = "lower right",fontsize='x-large')
    #     text0=0.7; text2=0.03; text4=0.03
    #     text1=0.35; text3=.9;  text5=0.8

    else:
        axhx.legend(loc = "upper left",fontsize='x-large')
        text0=0.1; text2=0.55; text4=0.55
        text1=0.6; text3=0.2;  text5=0.1
    
    axhx.text(text0, text1, r'$\chi^2/\nu\approx$'+str(np.round(parametrs[3],2)),fontsize=15, color='red', transform = axhx.transAxes)
    axhx.text(text2, text3, r'$A=$'+str(np.round(parametrs[0][0],3))+r"$\pm$"+str(np.round(parametrs[1],3)),fontsize=15,color='black' ,transform = axhx.transAxes)
    axhx.text(text4, text5, r'$B=$'+str(np.round(parametrs[0][1],3))+r"$\pm$"+str(np.round(parametrs[2],3)), fontsize=15, color='black',transform = axhx.transAxes)
    
    label_ylist = r'${:.2f}$'
    label_xlist = r'${:.2f}$'
    
    epsilon=np.array([cx,cx,cx])
    
    if(so==4 or so==5):
        ylist=np.linspace(0,max(dy0)+epsilon[ko,1],7)
    else:
        ylist=np.linspace(min(dy0)-epsilon[ko,1],max(dy0)+epsilon[ko,1],7)
    xlist=np.linspace(min(zrtx)-epsilon[ko,0],max(zrtx)+epsilon[ko,0],7)
    
    # rxm=np.array([[r"$\ln(\langle L \rangle)$",r"$\ln(\langle CV \rangle)$"],
    #             [r"$\ln(\langle L \rangle)$",r"$\ln(\langle CV \rangle)$"],
    #             [r"$\ln(\langle L \rangle)$",r"$\ln(\langle k \rangle)$"]])

    rxm=np.array([[r"$\langle E \rangle$",r"$\langle CV \rangle$"],
            [r"$\langle E \rangle$",r"$\langle CV \rangle$"],
            [r"$\langle E \rangle$",r"$\langle CV \rangle$"]])

    # rxm=np.array([[r"$\ln(\langle C \rangle)$",r"$k$"],
    #         [r"$\ln(\langle C \rangle)$",r"$k$"],
    #         [r"$\ln(\langle C \rangle)$",r"$k$"]])

    # rxm=np.array([[r"$\langle C \rangle$",r"$CV$"],
    #         [r"$\langle C \rangle$",r"$CV$"],
    #         [r"$\langle C \rangle$",r"$CV$"]])
    
    axhx.xaxis.set_major_locator(mticker.FixedLocator(xlist))
    axhx.set_xticklabels([label_xlist.format(x) for x in xlist],fontsize=14)

    axhx.yaxis.set_major_locator(mticker.FixedLocator(ylist))
    axhx.set_yticklabels([label_ylist.format(x) for x in ylist],fontsize=14)

    axhx.set_ylabel(r""+str(rxm[ko,0])+"",fontsize=20)
    axhx.yaxis.set_tick_params(labelsize=20)

    axhx.set_xlabel(r""+str(rxm[ko,1])+"",fontsize=20)
    axhx.xaxis.set_tick_params(labelsize=20)

    if(so == 4 or so == 5):
        axhx.set_ylim(0,max(dy0)+epsilon[ko,1])
    else:
        axhx.set_ylim(min(dy0)-epsilon[ko,1],max(dy0)+epsilon[ko,1])
    axhx.set_xlim(min(zrtx)-epsilon[ko,0],max(zrtx)+epsilon[ko,0])
    
    # plt.yscale("log")
    # plt.xscale("log")

def grafico_funcion_trozos(tx,io,oo,azhx,ovx,ovy,ozvx,oerry,oaax,obbx,onum_par,xnome_exper):

    plot_datas(oo,azhx,ovx,ovy,oerry) 
    
    mpart=[]    
    
    for ot in [0,1]:
        
        if(ot==0):
            print("cv<1.35")
            ox_min_ajuste=min(ovx)-0.03; ox_max_ajuste=max(vx)+0.03
            xpr0, xpr1, xpr2, xpr3, xpr4, xpr5 = graf_root(ovx,ovy,oerry,oaax,obbx,ox_min_ajuste,ox_max_ajuste,onum_par,"ajuste")

        # if(ot==0):
        #     print("cv>1.35")
        #     ox_min_ajuste=1.25; ox_max_ajuste=max(ovx)+0.03
        #     xpr0, xpr1, xpr2, xpr3, xpr4, xpr5 = graf_root(ovx,ovy,oerry,oaax,obbx,ox_min_ajuste,ox_max_ajuste,onum_par,"ajuste")
        mpart.append([xpr0[0],xpr0[1],xpr1,xpr2])
        parts=[xpr0,xpr1,xpr2,xpr5]
        
    Xo,erXo=intercepto(np.array(mpart))
    
    x_min=min(ovx)-0.03; x_max=max(ovx)+0.03
        
    plot_model_trozos(azhx,funct_modelo(np.array(mpart),np.linspace(x_min,Xo,100),np.linspace(Xo,x_max,100)))
    
    tx.append([xpr0[0],xpr1,xpr0[1] ,xpr2, xpr3, xpr4, xpr5])
    configure_axes_canvas(io,oo,azhx,ovx,ovy,oerry,ozvx,xnome_exper[io],Xo,erXo,parts) 
    return tx


def grafico_funcion(tx,io,oo,azhy,ovx,ovy,ozvx,ozvy,oerry,oaax,obbx,onum_par,xnome_exper):
    if(io==4 or io == 5):
        ct=[.1,.03]
    else:
        ct=[.1,.01]
        
    plot_datas(oo,azhy,ovx,ovy,oerry,ozvx,ozvy) 
    
    mpart=[]    
    

    # ox_min_ajuste=min(ovx)-0.03; ox_max_ajuste=max(ovx)+0.03
    # xpr0, xpr1, xpr2, xpr3, xpr4, xpr5 = graf_root(ovx,ovy,oerry,oaax,obbx,ox_min_ajuste,ox_max_ajuste,onum_par,"ajuste")
    # mpart.append([xpr0[0],xpr0[1],xpr1,xpr2])
    # parts=[xpr0,xpr1,xpr2,xpr5]
            
    # xo,fo=funct_model_completa(xpr0,np.linspace(ox_min_ajuste,ox_max_ajuste,100))
    # plot_model(azhy,[xo,fo])
    # tx.append([xpr0[0],xpr1,xpr0[1] ,xpr2, xpr3, xpr4, xpr5])
    # configure_axes_canvas(io,oo,azhy,ovx,ovy,oerry,ozvx,xnome_exper,parts,ct) 
    # return tx

num_par=2


fig = plt.figure(figsize=(17,9))
gs = fig.add_gridspec(2,3)
xcolr=["#9c46b1","#d1af00","#ac0000"]
yy_max=3;xx_max=2

#dataExp=["50sMar07","50sMar10","50sMar04","50sJan29","50sJan21","50sJan14","50sDez20"]
namExp=["Mar07","Mar10","Mar04","Jan29","Jan21","Jan14"]
nom_exper=["SUAM07","SUAM10","SUAM04","SUAJ29","SUAJ21","SUAJ14"]
nNo=[0,1,2,0,1,2]

lisdat=[3,2]; col_err=[1,4]; k=5

zn=np.array([lisdat,lisdat,lisdat])
erzn=np.array([col_err,col_err,col_err])

r=0; tx=[]
aax=0.02; bbx=0.05; xn=250

for i in range(6):
    xsubdir=str(xn)+"s"+namExp[i]

    if(i>2):
        r=1
    oo=nNo[i]      

    azh = fig.add_subplot(gs[r,oo])

    xdir="restW10T250/"
    fdat=np.loadtxt(xdir+str(xsubdir)+"Tot.txt")#datos empiricos
    #fdat=np.loadtxt("metrica_"+nom_exper[i]+".txt")#datos empiricos
    #frd=ajuste_distribucion_metricas()
    print(xdir+str(xsubdir)+".txt")
    zvx=(fdat[:,zn[oo,0]])
    zvy=(fdat[:,zn[oo,1]])
    err=np.array(fdat[:,erzn[oo,1]])

    dat_vx=(fdat[:,zn[oo,0]]) #eje X   
    dat_vy=((fdat[:,zn[oo,1]]))#eje Y
    
    # if(i==0):
    #     k=30
    # if(i==1):
    #     k=5.3
    # if(i==2):
    #     k=5
    # if(i==3):
    #     k=10
    # if(i==4):
    #     k=10
    # if(i==5):
    #     k=10

    dat_erry=k*np.array(fdat[:,erzn[oo,1]])

    

    mxtz=np.array(dat_vy)
    xcv=np.array(dat_vx)
    vx=dat_vx#list(xcv[(mxtz>5) & (mxtz<26)])
    vy=dat_vy#list(np.log(mxtz[(mxtz>5) & (mxtz<26)]))
    erry=k*dat_erry#list(dat_erry[(mxtz>5) & (mxtz<26)])


    # if(i==1):
    #     vx=(fdat[:,zn[oo,0]])
    #     vy=(fdat[:,zn[oo,1]])
    #     erry=k*(fdat[:,erzn[oo,1]])

    # if(i==4):
    #     erry=k*(fdat[:,erzn[oo,1]][vy<0.05])
    #     vx=(fdat[:,zn[oo,0]][vy<0.05])
    #     vy=(fdat[:,zn[oo,1]][vy<0.05])
        
    

    # if(i==0):
    #     erry=k*(fdat[:,erzn[oo,1]][vx>0.59])
    #     vy=(fdat[:,zn[oo,1]][vx>0.59])
    #     vx=(fdat[:,zn[oo,0]][vx>0.59])

    # if(i==3):
    #     erry=k*(fdat[:,erzn[oo,1]][vx>0.9])
    #     vy=(fdat[:,zn[oo,1]][vx>0.9])
    #     vx=(fdat[:,zn[oo,0]][vx>0.9])
        
    # if(i==4):
    #     erry=k*(fdat[:,erzn[oo,1]][vx<1.28])
    #     vy=(fdat[:,zn[oo,1]][vx<1.28])
    #     vx=(fdat[:,zn[oo,0]][vx<1.28])
        
    # if(i==5):
    #     erry=k*(fdat[:,erzn[oo,1]][vx>1])
    #     vy=(fdat[:,zn[oo,1]][vx>1])
    #     vx=(fdat[:,zn[oo,0]][vx>1])
    
# print(len(tx))
    #px=
    grafico_funcion(tx,i,oo,azh,vx,vy,zvx,zvy,erry,aax,bbx,num_par,nom_exper[i])


# parats=np.array(px)
# fun_print_parametros(parats,nom_exper)


plt.tight_layout()
fig.savefig('fhc.jpg')
plt.show()


