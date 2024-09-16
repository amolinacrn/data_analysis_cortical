from email import contentmanager
import numpy as np
import scipy.io
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def id_clust(kfile):
    id_clust=[]
    for time_file in os.listdir(kfile):
        id=(time_file.split('_')[4]).split('.')[0]
        id_clust.append(id)
    id_clust.sort()
    return id_clust

def cargar_datos(adrespr,ndir,rbin):#,foldprincp,foldresult,interval_bins):

    mxrango = slice(0,22)
    kpath="{}/{}".format(adrespr,ndir) 
    conjunto_datos=os.listdir(kpath)
    
    for dw,xsubdir in enumerate(conjunto_datos):
        kfile="{}/{}/{}{}/{}{}".format(adrespr,ndir,xsubdir[mxrango],dw,xsubdir[mxrango],dw)
        conteos_spk=[];print(kfile)
        cluster_activos=0; 
        longitud_matriz=[]
        for id in id_clust(kfile):
            #kfil_spk="{}/{}/{}/{}/sua_2020{}_{}_{}.txt".format(adrespr,ndir,xsubdir,xsubdir,ndir,dw,id)
            kfil_spk="{}/{}/{}{}/{}{}/sua_2020{}_{}_{}.txt".format(adrespr,ndir,xsubdir[mxrango],dw,xsubdir[mxrango],dw,ndir,dw,id)
            xpk=np.loadtxt(kfil_spk)
            mispk=xpk[2:len(xpk[:-1])]
            nbins=np.arange(mispk[0],mispk[-1],rbin)               
            count, edges=np.histogram(mispk, bins = nbins)
            
            if(len(count)!=0):
                conteos_spk.append(count)
                longitud_matriz.append(len(count))
                cluster_activos+=1
        # print(cluster_activos)    
        # mx_convariance=np.zeros((cluster_activos,max(longitud_matriz)))

    
dir_datos="datos_cv/Experimts_250s"
experimentos=os.listdir(dir_datos)

for l,xndr in enumerate(["Exp_Mar04"]):#experimentos):
    cargar_datos(dir_datos,xndr,0.05)#,foldprincp,foldresult,interval_bins)