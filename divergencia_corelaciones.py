from email.mime import base
from math import pi
import numpy as np
import os
from scipy.stats import entropy
from libreriaMolCr import funct_histograma, delete_nod
import matplotlib.pyplot as plt
import math
from clasificar import *
from collections import Counter

def D_js(P, Q):
    xM = 0.5*(P+Q)
    KLDpm = entropy(P, xM,base=2)
    KLDqm = entropy(Q, xM,base=2)
    JSDpq = 0.5*(KLDpm + KLDqm)
    return JSDpq

def cargar_comunidades():
    
    adrespr=os.listdir('ccMar4')
    xm=np.genfromtxt("ccMar4/"+adrespr[0])[:,0]
    ixd=np.genfromtxt("comunidades.txt")  
    fig, zx = plt.subplots(1,6 , figsize=(15,5 ))
    comuz=[]; comux=[]; comuy=[]; xCVs=[]; idCVs=[];comuw=[]; 
    fdat=[]
    for k,r in enumerate([4]):#,158,381,540,884]):
        io,jo=np.where(ixd==r)
        for i in io:
           
            cardat=np.genfromtxt("ccMar4/rsurr_"+str(int(ixd[i,1]))+"_"+str(int(ixd[i,2]))+".txt")[:,1]
            fdat.append(cardat/np.sum(cardat))
            #print(str(int(ixd[i,1]))+"\t"+str(int(ixd[i,2])))
            cardat=cardat-min(cardat)
            zx[k].plot(xm,cardat,"d")
    plt.show()
    fun_corr=np.array(fdat)
    return xm,fun_corr

def corr_ncc():

    adrespr=os.listdir('ccMar4')

    xm=np.genfromtxt("ccMar4/"+adrespr[0])[:,0]
    ixd=np.genfromtxt("id_sua_2020Mar04.txt")    
    fdat=[]
    for i in range(len(ixd)-1):
        for j in range(i+1,len(ixd)):
            cardat=np.genfromtxt("ccMar4/rsurr_"+str(int(ixd[i]))+"_"+str(int(ixd[j]))+".txt")[:,1]
            if(any(cardat)==True):
                fdat.append(cardat/np.sum(cardat))
                #print(str(int(ixd[i]))+"\t"+str(int(ixd[j])))
    fun_corr=np.array(fdat)
    return xm,fun_corr

def distancias_corr():

    _,fdat=cargar_comunidades()

    fun_corr=(np.array(fdat))

    MxRet=np.zeros((len(fun_corr),len(fun_corr)))
    
    for i in range(len(fun_corr)-1):
        print(i)
        for j in range(i+1,len(fun_corr)):
            xrest=D_js((fun_corr[i]),abs(fun_corr[j]))
            if(math.isnan(xrest)==True):
                print("fila nula")
            MxRet[i,j]=xrest; MxRet[j,i]=xrest
        
    fxile_txt="CorrMar4_tct.txt"
    fxdata = open(fxile_txt, 'w')

    try:
        for r in range(len(MxRet[:,0])):
            for s in range(len(MxRet[:,0])):
                fxdata.write(str(MxRet[r,s])+"\t")    
            fxdata.write("\n")    
    finally:
        fxdata.close()  
#distancias_corr()


xdm=(np.genfromtxt("DJSspikes.txt"))
x,y=funct_histograma(xdm,500)
plt.plot(x,y)
plt.show()
#x,y=cargar_comunidades()
#plt.matshow(corr(),cmap="plasma",aspect='auto')
# mxD=np.array(xdm)
# # vd=np.array(np.triu(mxD).flatten())
# # umbral_min=np.min(vd[vd>0])
# # par_red=np.max(vd[vd>0])
# # #umbral,xDaux=componente_gigante(mxD,umbral_min,par_red,0,False)
# xmz=np.array(componente_significativa(mxD,0.0075,0,False))
# draw_comunidades(xmz,"CorrMar4.txt","w10n50",True)
