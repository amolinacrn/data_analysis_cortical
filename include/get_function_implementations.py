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
from bct.algorithms import number_of_components
import re
from os import remove


def delet_zeros(original):
    nueva = []
    for dato in original:
        if dato != 0:
            nueva.append(dato)
    return nueva

def numero_de_componentes(M,dirfile):
    print(dirfile)
    if number_of_components(np.abs(M)) > 1:
        remove(dirfile)
    else:
        print("la matriz esta conectada \t")

def ADinf(AN):
    if AN<2:
        val=1/math.sqrt(AN)*math.exp(-1.2337141/AN)*(2.00012+(0.247105-(0.0649821-(0.0347962-(0.0116720-0.00168691*AN)*AN)*AN)*AN)*AN)
    else:
        val=math.exp(-math.exp(1.0776-(2.30695-(0.43424-(0.082433-(0.008056-0.0003146*AN)*AN)*AN)*AN)*AN))
    
    return val


def ns_detect(A,T):
    tn=np.zeros(len(A))
    
    for i in range(len(A)):
        tn[i]=A[i]
    
    Zn=tn/T
 
    if Zn[-1]>=1:
        print('Error: T debe ser estrictamente superior a la hora de llegada del último pico');
 
    N=len(Zn)
    
    AN=-N-1/N*np.dot([(2*i-1) for i in range(1,N+1)],np.log(Zn*(1-Zn[::-1])))

    pvalue=1-ADinf(AN)
   
    if pvalue<math.pow(10,-5):
        print('pvalue<10^-5,limit of accuracy exceeded, pvalue and NS values may not be precise'); 
    #NS=-math.log(pvalue)
    #print("NS: ",NS," pvalue: ",pvalue)

def real_inhibitory_links(cmt,idmx,dfil):
    contar=0 
    
    for j in range(len(idmx)):
        for k in range(len(idmx)):
            if(cmt[j,k]==0):
                idmx[j,k]=0   
                        
    for i in range(len(dfil[:,0])):
        for j in range(len(idmx)):
                for k in range(len(idmx)):
                    if(str(math.trunc(dfil[i,0]))+"."+str(math.trunc(dfil[i,1])))==str(idmx[j,k]):
                        contar+=1
                        print(idmx[j,k])
    print(contar)
                    
def conver_datos(mx_address,dirfile):
    mcc = np.genfromtxt(mx_address)
    
    nNf=len(mcc[:,1])
   
    try:
            
        fdata = open(dirfile+"/"+mx_address+".txt", 'w')
        print(dirfile)
        
        for i in range(nNf):
            for j in range(1,4):
                fdata.write(str(mcc[i][j]))
                fdata.write("\t")
            fdata.write("\n")  
    
    finally:
        fdata.close()
  
def eliminate_excess_clusters(path_dir_files):
        
    path_principal = os.listdir(path_dir_files)
    
    for num_exp in path_principal: 
           
        adres_path=path_dir_files+"/"+num_exp   
           
        adres_dir = os.listdir(adres_path)

        mxrango = slice(20,27)

        nfiles1=[]; nfiles2=[]; nfiles3=[]

        addres=[]

        for i, xdir in enumerate(adres_dir):
        
            dir_matriz_empirica = str(adres_path)+"/"+str(xdir)+"/"+str(xdir)
            #print(dir_matriz_empirica)
            xadres_dir = os.listdir(dir_matriz_empirica)
            
            addres.append(dir_matriz_empirica)
        
            if (i==0):
                for j in range(len(xadres_dir)):
                    nfiles1.append(xadres_dir[j])
            if (i==1):
                for j in range(len(xadres_dir)):
                    nfiles2.append(xadres_dir[j])
            if (i==2):
                for j in range(len(xadres_dir)):
                    nfiles3.append(xadres_dir[j])            
        
        mxc=[nfiles1,nfiles2,nfiles3]; id1=[]

        mxsize = [len(mxc[0]),len(mxc[1]),len(mxc[2])]; mxsize.sort()
        
        for mxfil in mxc:
            if(np.min(mxsize)==len(mxfil)):            
                for i in range(len(mxfil)):
                    id1.append(int(re.sub(r'.txt','', mxfil[i][mxrango])))
                break
        
        id1.sort()               
        
        #print("Archivos eliminados") 
        
        for a_k in range(len(mxc)):
            files_iguales=[]
            for k in range(len(id1)):
                for a_i,mxfil1 in enumerate(mxc[a_k]):
                    if(id1[k]==int(re.sub(r'.txt','', mxc[a_k][a_i][mxrango]))):
                        files_iguales.append(mxfil1)
  
            for k in range(len(files_iguales)):
                for a_i in range(len(mxc[a_k])):
                    if(files_iguales[k]==mxc[a_k][a_i]):
                        mxc[a_k][a_i]='o'          
            #print(mxc[a_k])            
            for a_i in range(len(mxc[a_k])):
                if(mxc[a_k][a_i]!='o'):
                    xddir =str(addres[a_k])+"/"+str(mxc[a_k][a_i])
                    remove(xddir) 
                    #print(mxc[a_k][a_i])
        print(mxsize)                        
                   
            
def graf_inhib(address,mclus1,mclus2):
       
    id_cluster_1=mclus1
    id_cluster_2=mclus2
    
    mx_address=str(address+'/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+ '.txt')
    print(mx_address)
    xs, ys = np.loadtxt(mx_address,skiprows=0, usecols=[0,1],unpack=True)#datos empiricos
     
    plt.plot(xs,ys) 
        
    plt.show()            



def delete_unconnected_nodes(M,dirfile):

    if number_of_components(np.abs(M)) > 1:

        print("Matriz desconectada \t",dirfile)

        for j in range(len(M)):
                
            for i in range(len(M)):
                if(any(M[:,int(i)])!=True):
                    A=np.delete(M, i, axis=1)
                    break  
                
            for i in range(len(M)):    
                if(any(A[int(i),:])!=True):
                    B=np.delete(A, i, axis=0) 
                    break
        
            M = np.resize(M,(len(B),len(B)))
            M=B; contador=-1
            
            for a_i in range(len(M)):   
                if(any(M[int(a_i),:])==True):
                    contador+=1
            
            if(a_i==contador):
                break
        
        try:
                
            fdata = open(dirfile, 'w')
            #print(dirfile)
            
            for wrfil in range(len(M)):
                for wrcol in range(len(M)):
                    fdata.write(str(M[wrfil][wrcol]))
                    fdata.write("\t")
                fdata.write("\n")  
        
        finally:
            fdata.close()     
    
    else:
        
        print("Matriz conectada \t",dirfile) 
           
        try:
            fdata = open(dirfile, 'w')
            #print(dirfile)
            
            for wrfil in range(len(M)):
                for wrcol in range(len(M)):
                    fdata.write(str(M[wrfil][wrcol]))
                    fdata.write("\t")
                fdata.write("\n")  
        
        finally:
            fdata.close() 

def thresholded_conetivity_matrix(mcc,xdir,nsurr, dir_file,xbin,xw):
    '''
    mcc: matriz de conectividad empirica
    xdir: directorio principal de datos sustitutos
    nsurr: numero de experimentos sustitutos
    xbin: archivo de matrix de conectividad  
    '''
    mxc_cuadrado = np.zeros((len(mcc),len(mcc)))
    mxc_surr = np.zeros((len(mcc),len(mcc)))
    media_mcc = np.zeros((len(mcc),len(mcc)))
    devstd_mcc = np.zeros((len(mcc),len(mcc)))
    treshold = np.zeros((len(mcc),len(mcc)))
    
    for k in range(nsurr):
        
        dir_mcc_surr = xdir+"/surr_"+str(k+1)+"/"+"ncc_0.1-"+str(xbin)+"-"+str(xw)+"msec/surr_"+str(k+1)+"_CC_Symmetric.txt"
        print(dir_mcc_surr)
        for i in range(len(mcc)):         
            for j in range(len(mcc)):
                print(i,j)
                mx=np.genfromtxt(dir_mcc_surr)[i,j]
                mxc_surr[i,j] += mx 
                mxc_cuadrado[i,j] += mx*mx
                
    for i in range(len(mcc)):
        for j in range(len(mcc)): 
            
            media_mcc[i,j] = (1/nsurr)*mxc_surr[i,j]
            devstd_mcc[i,j] =(1/nsurr)*mxc_cuadrado[i,j]            
            treshold[i,j] = media_mcc[i,j]+2*devstd_mcc[i,j]

            if (mcc[i,j] > treshold[i,j]):    
                mcc[i,j]
            else:  
                mcc[i,j]=0

    try:
            
        fdata = open(dir_file+".txt", 'w')
        print(fdata)
        
        for wrfil in range(len(mcc)):
            for wrcol in range(len(mcc)):
                fdata.write(str(mcc[wrfil][wrcol]))
                fdata.write("\t")
            fdata.write("\n")  
    
    finally:
        fdata.close()    
    
def treshold_shufling( nfilas,num_surr,dir_surr,id_cluster_1, id_cluster_2, enlacex):
    
    #dir_contar_filas= str(dir_surr+"/rsurr_"+str(id_cluster_1)+"_"+str(id_cluster_2)+".txt")
    
    #nfilas = np.loadtxt(dir_contar_filas,skiprows=0, usecols=[1],unpack=True) #datos surrogate
    
    mxsurr = np.zeros((nfilas,num_surr))
    
    sigf=np.zeros(nfilas)
    
    vmedia=np.zeros(nfilas)
        
    for columnas in range(num_surr):
        addirsurr = dir_surr+"/surr_"+str(columnas+1)
       
        mx_adir=str(addirsurr+'/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+'.txt')
              
        dxsr, dysr = np.loadtxt(mx_adir,skiprows=0, usecols=[0,1],unpack=True)   
   
        
        if(enlacex):                        
            for filas in range(len(dysr)):
                baxar=abs(dysr[filas]-np.max(dysr))
                mxsurr[filas][columnas] = baxar+np.max(dysr) #inhibitoria
               
        else:                       
            for filas in range(len(dysr)):
                mxsurr[filas][columnas] = dysr[filas]#(dysr[filas]-np.min(dysr)) #exitatoria
        
                            
    for contfil in range(len(dysr)):         
        vmedia[contfil]=statistics.mean(mxsurr[contfil])
        dstdar=statistics.stdev(mxsurr[contfil])
        sigf[contfil]=vmedia[contfil]+2*dstdar # significancia 5% 
    
    '''
    k=0
    
    for i in range(len(sigf)):
        if(k<10 or k>40):
            sigf[i]=np.min(sigf) 
        k+=1 
    '''    
    return sigf,vmedia,mxsurr


def correlation_Matrix_grafics(id_tot,address,id_inhib,num_clusters,num_cl,dir_surr,versurr = False):
    
    #fig, (axx1) = plt.subplots(1, 1, figsize=(5, 5))
    
    #------------------------  inicializacion ---------------------------------
       
    id_exitator = id_tot.copy()
    id_inhibitor = id_tot.copy()
    
    r=0
    for iz in range(len(id_inhib)):
        for im in range(len(id_tot)):
           
            if id_inhib[iz][0] == id_tot[im][0] and id_inhib[iz][1] == id_tot[im][1]:
                id_exitator[im][0]=0 
                id_exitator[im][1]=0 
                
    for iz in range(len(id_exitator)):
        for im in range(len(id_tot)):
            if id_exitator[iz][0] == id_tot[im][0] and id_exitator[iz][1] == id_tot[im][1]:
                id_inhibitor[im][0]=0 
                id_inhibitor[im][1]=0
            
    for i in range(num_cl-1):
        
        id_cluster_1=math.trunc(id_tot[i][0])
        id_cluster_2=math.trunc(id_tot[i][1])

        for j in range(num_cl):
            
            id_cluster_exit_1=math.trunc(id_exitator[j][0])
            id_cluster_exit_2=math.trunc(id_exitator[j][1])
            
            id_cluster_inhib_1=math.trunc(id_inhibitor[j][0])
            id_cluster_inhib_2=math.trunc(id_inhibitor[j][1])
                            
            if (id_cluster_1 == id_cluster_exit_1!=0 and id_cluster_2 == id_cluster_exit_2!=0):
               
                mx_address=str(address+'/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+ '.txt')
                xs, ys = np.loadtxt(mx_address,skiprows=0, usecols=[0,1],unpack=True)#datos empiricos
            
                n_dat=len(ys)
                y_s = n_dat*[0]
            
                for s in range(n_dat):     
                    y_s[s] = ys[s] #- np.min(ys) # np.min para exitatorios
                
                tresholdh = 0.05
               
                #--------------------------------------------------------------
                print(id_cluster_1," ",id_cluster_2," ",mx_address," --> ", "Exitatorio ") 
                            
                if (versurr):
                    
                    fig, ax1 = plt.subplots(1, 1, figsize=(8, 7)) 
                    
                    if(np.max(y_s) > tresholdh ): 
                        ax1.plot(xs,y_s,label="Empirical data,exitatoria") 
                        varimg = "graficas_exit_aceptadas"
                    else:
                        ax1.plot(xs,y_s,label="Empirical data,exitatoria") 
                        varimg = "graficas_exit_rechazadas"
                        
                    fig.savefig("include/"+str(varimg)+"/file"+str(id_cluster_1)+"_"+str(id_cluster_2)+"_.jpg")
                                        
            
            elif id_cluster_1 == id_cluster_inhib_1!=0 and id_cluster_2 == id_cluster_inhib_2!=0:
                    
                mx_address=str(address+'/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+ '.txt')              
                vxs, ys = np.loadtxt(mx_address,skiprows=0, usecols=[0,1],unpack=True)
            
                n_dat=len(ys)
                y_ss = n_dat*[0]
            
                for s in range(n_dat): 
                    mbaxar=abs(ys[s] - np.max(ys))    
                    y_ss[s] = mbaxar + np.max(ys)# np.max para inibitorios
                
                tresholdh = 0.05

                print(id_cluster_1," ",id_cluster_2," ",mx_address,"-->  ", "Inhibitorio")
                
                if (versurr):
                
                    fig, ax2 = plt.subplots(1, 1, figsize=(8, 7)) 
            
                    if(np.max(y_ss) > tresholdh ): 
                        ax2.plot(vxs,y_ss,label="Empirical data,inhibitoria") 
                        varimg = "graficas_inhib_aceptadas"
                    else:
                        ax2.plot(vxs,y_ss,label="Empirical data,inhibitoria") 
                        varimg = "graficas_inhib_rechazadas"
                        
                    fig.savefig("include/"+str(varimg)+"/file"+str(id_cluster_1)+"_"+str(id_cluster_2)+"_.jpg")

                                    
            else:
                j
            matplotlib.pyplot.close()    




def correlation_Matrix(id_tot,address,id_inhib,num_clusters,num_cl,num_surr,dir_surr,Wc,experit,versurr = False):
    
    #fig, (axx1) = plt.subplots(1, 1, figsize=(5, 5))
    
    #------------------------  inicializacion ---------------------------------
    
    file_txt="ResMC/tmc_"+str(experit)+".txt" # matrix de correlaciones no umbralizada 
    
    fdata = open(file_txt, 'w')
    
    #------------------------------------------------------------------------
       
    id_exitator = id_tot.copy()
    id_inhibitor = id_tot.copy()
    
    r=0
    for iz in range(len(id_inhib)):
        for im in range(len(id_tot)):
           
            if id_inhib[iz][0] == id_tot[im][0] and id_inhib[iz][1] == id_tot[im][1]:
                id_exitator[im][0]=0 
                id_exitator[im][1]=0 
                
    for iz in range(len(id_exitator)):
        for im in range(len(id_tot)):
            if id_exitator[iz][0] == id_tot[im][0] and id_exitator[iz][1] == id_tot[im][1]:
                id_inhibitor[im][0]=0 
                id_inhibitor[im][1]=0
            
    TCM_trsupr = np.zeros(( num_clusters, num_clusters ))
    TCM_trinfr  = np.zeros((  num_clusters, num_clusters ))
    TCM_treshold  = np.zeros((  num_clusters, num_clusters ))

    
    l, k=0, 0
    
    try:
        for i in range(num_cl-1):
            
            id_cluster_1=math.trunc(id_tot[i][0])
            id_cluster_2=math.trunc(id_tot[i][1])
 
            for j in range(num_cl):
                
                id_cluster_exit_1=math.trunc(id_exitator[j][0])
                id_cluster_exit_2=math.trunc(id_exitator[j][1])
                
                id_cluster_inhib_1=math.trunc(id_inhibitor[j][0])
                id_cluster_inhib_2=math.trunc(id_inhibitor[j][1])
                              
                if (id_cluster_1 == id_cluster_exit_1!=0 and id_cluster_2 == id_cluster_exit_2!=0):
                    
                    tipo_enlace = False #para exitatoria
                    
                    mx_address_surr=str(dir_surr+'/surr_1/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+ '.txt')
                    xsr, ysr = np.loadtxt(mx_address_surr,skiprows=0, usecols=[0,1],unpack=True) #datos surrogate
                    
                    mx_address=str(address+'/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+ '.txt')
                    xs, ys = np.loadtxt(mx_address,skiprows=0, usecols=[0,1],unpack=True)#datos empiricos
                
                    n_dat=len(ys)
                    y_s = n_dat*[0]
                    y_sr = n_dat*[0]
                
                    for s in range(n_dat):     
                        y_s[s] = ys[s] #- np.min(ys) # np.min para exitatorios
                        y_sr[s] =ysr[s]# - np.min(ysr) # np.min para exitatorios
                   
                    #y_s=savgol_filter(y_s,20,3, mode="nearest") #  surrogate
                    
                    fstd,fmed,mxsurr_hipt=treshold_shufling(2*Wc+1,num_surr, dir_surr, id_cluster_1, id_cluster_2, tipo_enlace)  
                    
                    #fstd=savgol_filter(fstd,20,3, mode="nearest") #  surrogate

                    ko=0; cop_y_s = y_s.copy(); ffstd_copi = fstd.copy()
                    '''
                    for io in range(len(y_s)):
                        if(ko<25 or ko>175):
                            cop_y_s[io]=0
                            ffstd_copi[io]=0
                        ko+=1 
                    '''
                    #tresholdh = np.max(ffstd_copi)
                    tresholdh = np.max(fstd)
                    
                    muestra_hopotesis=[]
                    for ky in range(len(mxsurr_hipt[0,:])):
                        muestra_hopotesis.append(np.max(mxsurr_hipt[:,ky]))
                    
                    
                    #tresholdh_hipt=statistics.mean(muestra_hopotesis)+2*statistics.stdev(muestra_hopotesis)
                    
                    zr = (statistics.mean(muestra_hopotesis)-np.max(y_s))/(statistics.stdev(muestra_hopotesis)/math.sqrt(len(mxsurr_hipt[0,:])))  
                    
                    if( zr >=0):
                        pr = stats.norm.sf(zr) 
                    else:
                        pr = stats.norm.cdf(zr)
                    
                    k+=1
                    #-------------------------------------------------------------       
                    if(np.max(y_s) > tresholdh and pr < 0.05):
                        TCM_trsupr[l][k] =  np.max(y_s)
                    else:
                        TCM_trsupr[l][k]=0
                    
                    #--------------------------------------------------------------
                    print(id_cluster_1," ",id_cluster_2," ",mx_address," --> ", "Exitatorio ") 
                    
                    #if (id_tot[i][0]!=id_tot[i+1][0]):
                        #l+=1; k=l
                    
                    if (versurr):
                        
                        fig, ax1 = plt.subplots(1, 1, figsize=(8, 7)) 
                        
                        if(np.max(y_s) > tresholdh and pr < 0.05): 
                            ax1.plot(xs,y_s,label="Empirical data,exitatoria") 
                            varimg = "graficas_exit_aceptadas"
                        else:
                            ax1.plot(xs,y_s,label="Empirical data,exitatoria") 
                            varimg = "graficas_exit_rechazadas"
                            
                        fig.savefig("include/"+str(varimg)+"/file"+str(id_cluster_1)+"_"+str(id_cluster_2)+"_.jpg")
                        
           
                    
                    if(versurr):
                        fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))    
                        #ax1.plot(mxc,myc,"^")
                        ax1.plot(xs,y_s,label="Empirical data,exitatoria")    
                        #ax1.plot(xs,yy_s,label="Empirical data,exitatoria")                   
                        ax1.plot(xs,fstd,"-",label="$\mu siginificancia data",color="blue")    
                        #ax1.plot(xsr,y_sr,label="surrogate data w=20 ms")
                        ax1.plot(xsr,fmed,label="$\mu$ surrogate data ")
                        ax1.legend()
                        fig.savefig('include/graficas_surr/file'+str(id_cluster_1)+'_'+str(id_cluster_2)+'_.jpg')
                                            
                
                elif id_cluster_1 == id_cluster_inhib_1!=0 and id_cluster_2 == id_cluster_inhib_2!=0:
                    
                    tipo_enlace = True #para inhibitoria
                    
                    mx_address_surr=str(dir_surr+'/surr_1/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+ '.txt')
                    xsr, ysrr = np.loadtxt(mx_address_surr,skiprows=0, usecols=[0,1],unpack=True) #datos surrogate
                     
                    mx_address=str(address+'/rsurr_' + str(id_cluster_1)+ '_' + str(id_cluster_2)+ '.txt')              
                    vxs, ys = np.loadtxt(mx_address,skiprows=0, usecols=[0,1],unpack=True)
                
                    n_dat=len(ys)
                    y_ss = n_dat*[0]
                    y_ssr = n_dat*[0]
                    
                    k+=1
                
                    for s in range(n_dat): 
                        mbaxar=abs(ys[s] - np.max(ys))    
                        y_ss[s] = mbaxar + np.max(ys)# np.max para inibitorios
                        mbaxarx = abs(ysrr[s] - np.max(ysrr))
                        y_ssr[s] = mbaxarx + np.max(ysrr) # np.max para inibitorios
                    
                    #y_ss=savgol_filter(y_ss,20,3) #  surrogate  

                    fstdd,fmedd ,mxsurr_hipt=treshold_shufling(2*Wc+1,num_surr, dir_surr, id_cluster_1, id_cluster_2, tipo_enlace)

                    koo=0; cop_y_ss = y_ss.copy(); fstd_copi = fstdd.copy()
                    
                    #fstdd=savgol_filter(fstdd,20,3, mode="nearest") #  surrogate
                    
                    '''
                    for ioo in range(len(y_ss)):
                        if(koo<25 or koo>175):
                            cop_y_ss[ioo]=0
                            fstd_copi[ioo]=0    
                        koo+=1 
                    '''
                    #tresholdh = np.max(fstd_copi)
                    tresholdh = np.max(fstdd)
                    
                    muestra_hopotesis=[]
                    for ky in range(len(mxsurr_hipt[0,:])):
                        muestra_hopotesis.append(np.max(mxsurr_hipt[:,ky]))
                    
                    #tresholdh_hipt=statistics.mean(muestra_hopotesis)+2*statistics.stdev(muestra_hopotesis) 
                    
                    zr = (statistics.mean(muestra_hopotesis)-np.max(y_s))/(statistics.stdev(muestra_hopotesis)/math.sqrt(len(mxsurr_hipt[0,:])))  
                    
                    if( zr >=0):
                        pr = stats.norm.sf(zr) 
                    else:
                        pr = stats.norm.cdf(zr) 
                    
                    #print(tresholdh," ",np.max(y_ss))
                    #-------------------------------------------------------------       
                    if(np.max(y_ss) > tresholdh and pr < 0.05):
                        TCM_trsupr[l][k] =  np.max(y_ss)
                    else:
                        TCM_trsupr[l][k]=0
                    #--------------------------------------------------------------
                    print(id_cluster_1," ",id_cluster_2," ",mx_address,"-->  ", "Inhibitorio")
                    
                    #varimg="o"
                    
                    if (versurr):
                    
                        fig, ax2 = plt.subplots(1, 1, figsize=(8, 7)) 
                
                        if(np.max(y_ss) > tresholdh and pr < 0.05): 
                            ax2.plot(vxs,y_ss,label="Empirical data,inhibitoria") 
                            varimg = "graficas_inhib_aceptadas"
                        else:
                            ax2.plot(vxs,y_ss,label="Empirical data,inhibitoria") 
                            varimg = "graficas_inhib_rechazadas"
                            
                        fig.savefig("include/"+str(varimg)+"/file"+str(id_cluster_1)+"_"+str(id_cluster_2)+"_.jpg")
                    
                    
                    if (versurr):
                        
                        fig, ax2 = plt.subplots(1, 1, figsize=(8, 8))
                        #ax2.plot(xmxc,xmyc,"^")
                        ax2.plot(vxs,y_ss,label="Empirical data,inhibitoria")                        
                        ax2.plot(vxs,fstdd,"-",label="$\mu + 2\sigma$ surrogate data",color="blue")    
                        #ax2.plot(xsr,y_ssr,label="surrogate data w=20 ms")
                        ax2.plot(xsr,fmedd,label="$\mu$ surrogate data ")
                        ax2.legend()
                        fig.savefig('include/graficas_surr/file'+str(id_cluster_1)+'_'+str(id_cluster_2)+'_.jpg')
                                     
                else:
                    j
                matplotlib.pyplot.close()    
                
            if (id_tot[i][0]!=id_tot[i+1][0]):
                l+=1; k=l        
                    
        for xfil in range(num_clusters):
            for xcol in range(num_clusters):
                TCM_trinfr[xcol][xfil] = TCM_trsupr[xfil][xcol]    
                
        for yfil in range(num_clusters):
            for ycol in range(num_clusters):
                TCM_treshold[yfil][ycol] = TCM_trsupr[yfil][ycol]+ TCM_trinfr[yfil][ycol] 
                            
        for wrfil in range(num_clusters):
            for wrcol in range(num_clusters):
                fdata.write(str(TCM_treshold[wrfil][wrcol]))
                fdata.write("\t")
            fdata.write("\n")    
            
        #plt.show()
        
    finally:
        fdata.close()
    
    return TCM_treshold
        
def correlation_Excitator_Inhibitor_Matrix(idd_clusters_inhibitorios,id_clusters,loadmatrizmc,exprto,unic_exit,unic_inh):
    #fddata = open(fille_txt, 'w')
    
    varm= len(idd_clusters_inhibitorios)-1

    id_clusters_inhibitorios =np.zeros((varm,2))
    
    rcont_j = 0    
            
    for i_m in range(len(idd_clusters_inhibitorios)):            
        if idd_clusters_inhibitorios[i_m][0] !=0 and idd_clusters_inhibitorios[i_m][1]!=0:
            id_clusters_inhibitorios[rcont_j][0]=idd_clusters_inhibitorios[i_m][0]
            id_clusters_inhibitorios[rcont_j][1]=idd_clusters_inhibitorios[i_m][1]
            rcont_j +=1        
    
    n_mc=len(loadmatrizmc)
    n_c=len(id_clusters_inhibitorios)

    mx_conectivity=loadmatrizmc.copy()#loadmatrizmc[im][jm]
    matriz_conectividad=loadmatrizmc.copy()#loadmatrizmc[im][jm]
    copia_matrix_conectividad=loadmatrizmc.copy()#loadmatrizmc[im][jm]
          
    id_clusters_inhibitorios_col_1=n_c*[0]
    id_clusters_inhibitorios_col_2=n_c*[0]
     
    for i in range(n_c):
        id_clusters_inhibitorios_col_1[i]=math.trunc(id_clusters_inhibitorios[i][0])
        id_clusters_inhibitorios_col_2[i]=math.trunc(id_clusters_inhibitorios[i][1])

    for a_i in range(n_c):
        a_k=0
        for b_i in id_clusters:
            a_m=0
            for b_j in id_clusters:         
                if(id_clusters_inhibitorios_col_1[a_i] == b_i and id_clusters_inhibitorios_col_2[a_i] == b_j):
                    matriz_conectividad[a_k][a_m]=0  
                a_m +=1
            a_k +=1    
        
    # para matriz inhibitoria     
    if(unic_inh==1 and unic_exit ==0):  
        ffile_txt="mc_Wght_"+str(exprto)+"ihhib.txt"
        ffdata = open(ffile_txt, 'w')
        
        try:
    
            for z_i in range(n_mc):
                for z_j in range(n_mc):
                    matriz_conectividad[z_i][z_j]=copia_matrix_conectividad[z_i][z_j]-matriz_conectividad[z_i][z_j]     
            
            transpose_matriz_conectividad= np.zeros((n_mc,n_mc))
            
            for z_i in range(n_mc):
                for z_j in range(n_mc):
                    transpose_matriz_conectividad[z_j][z_i]= matriz_conectividad[z_i][z_j]      
                    
            for z_i in range(n_mc):
                for z_j in range(n_mc):
                    matriz_conectividad[z_i][z_j]=matriz_conectividad[z_i][z_j]+transpose_matriz_conectividad[z_i][z_j]          
        
            for wrfil in range(n_mc):
                for wrcol in range(n_mc):
                    ffdata.write(str(matriz_conectividad[wrfil][wrcol]))
                    ffdata.write("\t")
                ffdata.write("\n") 
                
            #mct_inhib= matriz_conectividad.copy()
                  
        finally:
            ffdata.close() 
        
    #para matriz exitatoria        
    elif(unic_inh==0 and unic_exit ==1):
        file_txt="mc_Wght_"+str(exprto)+"Exit.txt"
        fdata = open(file_txt, 'w')
        
        #matriz_conectividad = mx_conectivity  
        
        for z_i in range(n_mc):
            for z_j in range(n_mc):
                if(matriz_conectividad[z_i][z_j]!=0):
                    mx_conectivity[z_i][z_j]= 0     
        #--------------------------------------
    
        try:
            for wrfil in range(n_mc):
                for wrcol in range(n_mc):
                    fdata.write(str(mx_conectivity[wrfil][wrcol]))
                    fdata.write("\t")
                fdata.write("\n")
                
            #mct_exit = matriz_conectividad.copy()    
        finally:
            fdata.close() 
        
        
    else:
        print("la opcion que eligio no es correcta")
        print("intente:")
        print("     Exit = 1, inhib = 0")
        print("     Exit = 0, inhib = 1")
       
    
def connectivity_matrix_linkage (mctt_inhib,mctt_exit,tresshol_hard,stirp):

    fille_txt="mc_Wght_"+str(stirp)

    fddata = open(fille_txt, 'w')
    
    n_mc=len(mctt_exit)
      
    mx_exit= np.zeros((n_mc,n_mc))
    mx_inib= np.zeros((n_mc,n_mc))
    loadmatrizmc= np.zeros((n_mc,n_mc))
    
    for im in range(n_mc):
        for jm in range(n_mc):
            mx_inib[im][jm]=mctt_inhib[im][jm]#(mctt_inhib[im][jm]-np.min(mctt_inhib))/(np.max(mctt_inhib)-np.min(mctt_inhib))
            mx_exit[im][jm]=mctt_exit[im][jm]#(mctt_exit[im][jm]-np.min(mctt_exit))/(np.max(mctt_exit)-np.min(mctt_exit))
    
    rg = n_mc*n_mc*[0]     
    kj=0   
    for im in range(n_mc):
        for jm in range(n_mc):
            rg[kj]=mx_inib[im][jm]
            kj+=1
    x=0
    for i in range(len(delet_zeros(rg))):  
        if delet_zeros(rg)[i]<0.1:    
            x+=1  
            #print(x,delet_zeros(rg)[i])
            
    matrix_conectividade = np.zeros((n_mc,n_mc))
    matrix_addres = np.zeros((n_mc,n_mc))
    
    try:
        
        for xfil in range(n_mc):
            for xcol in range(n_mc):
                matrix_conectividade[xfil][xcol]=mx_exit[xfil][xcol] + mx_inib[xfil][xcol]
                matrix_addres[xfil][xcol] = matrix_conectividade[xfil][xcol]
                if matrix_conectividade[xfil][xcol]<tresshol_hard:
                    matrix_conectividade[xfil][xcol]=0
                fddata.write(str(matrix_conectividade[xfil][xcol]))
                fddata.write("\t")
            fddata.write("\n")
            
    finally:
        fddata.close()
    
    mxcontar=((len(matrix_addres))**2)*[0]
    
    contador=0
    for i in range(len(matrix_addres)):
        for j in range(len(matrix_addres)): 
            mxcontar[contador] = matrix_addres[i][j]
            contador+=1 

 
def tcm_graficos(matriz_conectividad):
    cmap = ListedColormap(["white","black"])#, "gold", "lawngreen", "lightseagreen","blue","magenta","cyan","green","red"])
    #cmap = ListedColormap([ "jet"])

    viridis =cm.get_cmap('plasma')
    newcolors = viridis(np.linspace(0, 1, 256))
    pink = np.array([1])#}[248/256, 24/256, 148/256, 1])
    newcolors[:0, :] = pink
    newcmp = ListedColormap(newcolors)

    '''c=ROOT.TCanvas("cV5","migrafico",300,500,550,470)'''
    
    n_mc=len(matriz_conectividad)

    mx_conectividad = np.zeros((n_mc,n_mc))

    mxc = n_mc*n_mc*[0] 

    contador_mxc=0

    for i in range(n_mc):
        for j in range(n_mc):
            mxc[contador_mxc] =  matriz_conectividad[i][j] 
            contador_mxc+=1

    ccont_mxc =0 

    for o in range(contador_mxc):
        if mxc[o]!=0:
            ccont_mxc +=1

    tcm = ccont_mxc*[0]
    contador_tcm=0

    for io in range(contador_mxc):
        if mxc[io] !=0:
            tcm[contador_tcm]=mxc[io]                   
            contador_tcm +=1   
    conx=0

    for i in range(n_mc):
        for j in range(n_mc):
            mx_conectividad[i][j] =  matriz_conectividad[i][j]  
                    
    for i in range(n_mc):
        for j in range(n_mc):
            if(mx_conectividad[i][j]==0):
                mx_conectividad[i][j]
            else:
                mx_conectividad[i][j]
    
    plt.title("Connectivity matrix with 20% of the most significant connections")        
    plt.imshow(mx_conectividad,cmap =newcmp)##si quiero en colores variado cambio cmap="jet"
    #plt.subplot(111)
    #plt.subplots_adjust(bottom=0.1,left=-0.065 ,right=1, top=0.989)
    #plt.xlabel("")
    #plt.ylabel("")
    cax = plt.axes([0.81, 0.1, 0.03, 0.82])
    
    plt.colorbar(cax=cax)
    #plt.tight_layout()
    plt.show()        
    
    ################################################################
    ############### grafico de red topologica  #####################
    ################################################################

    G = nx.from_numpy_matrix(mx_conectividad) #matriz_conectividad)
    pos = nx.spring_layout(G)
    for n, p in pos.items():
        G.nodes[n]['pos'] = p

    edge_x = []
    edge_y = []

    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.3, color='black'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='plasma',
            #reversescale=True,
            color=[],
            size=9,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'),
            line_width=0.4))
            
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))
        
    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text


    title = "Network Graph Demonstration"
    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(height=800, width=800,
                    paper_bgcolor='white',
                    plot_bgcolor='white',
                    #title=title,
                    xaxis=dict(showgrid=False, zeroline=False,
                            showticklabels=False, mirror=True),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, mirror=True)
                    ))

    fig.show()
    #fig.savefig('file.png')
    
    '''    
    H_x = np.array(np.log(tcm),dtype = float)
    H_y = np.ones(contador_tcm)

    h1 = ROOT.TH1D("h1","h1 title",100,-7,-0.02);

    h1.FillN(contador_tcm, H_x, H_y)
    h1.SetFillColor(47)
    h1.Draw()

    c.Update()
    '''

    input("Prompt: ") 
    
  
def numero_conexiones_tcm(addrres,tmc_inhib,tmc_exit,treshold_hard,solo_exit, solo_inhi):
     
    print("\n Tamaño de la matriz = ",len(addrres),"x", len(addrres))
    mxcontar=((len(addrres))**2)*[0]
    contador=0
    
    for i in range(len(addrres)):
        for j in range(len(addrres)): 
            mxcontar[contador] = addrres[i][j]
            contador+=1 
    num_contivity=len(delet_zeros(mxcontar))             
    n_conexxiones =(int)(num_contivity/2)
    n_itereaciones_latezar = n_conexxiones*100
    salida= "\n número de conexioes con treshold hard = " + str(n_conexxiones) +" -->  recablear " + str(n_itereaciones_latezar)+ "  veces"
       
    print(salida)
       
    return "\n"

     
  
