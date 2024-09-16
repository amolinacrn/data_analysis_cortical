import numpy as np
import os
import math
from pathlib import Path 
from shutil import rmtree
import random

def extract_clusters(xddir,dir,zdirect):
    mx = np.genfromtxt(xddir)
    print(mx[:,1])
    result = []
    
    for item in mx[:,1]:
    #for item in mmx:
        if item not in result:
            result.append(item)
    m=sorted(result)       
    for i in range(len(m)):
        print((math.trunc(m[i])))
    
    
        
    file_txt=zdirect+"/id_sua_"+str(dir)+".txt"
            
    fdata = open(file_txt, 'w')

    try:

        for zz_i in range(len(m)):
            fdata.write( str(math.trunc(m[zz_i]))+"\n")
                            
    finally:
        fdata.close()
    

def extract_data_cluster(tot_datos,id_clusters,a_i):
   
    n = len(tot_datos[:,0])

    mx_tiempo_id = []

    n_temporal_x_cluster=0
        
    for a_j in range(n):
        if (tot_datos[:,1][a_j]==id_clusters[a_i]):
            mx_tiempo_id.append(tot_datos[:,0][a_j])
            n_temporal_x_cluster +=1
    return mx_tiempo_id

def write_times_cluster (id_clusters,tot_datos,dir_adr,dirpath):
    
    print(dirpath)
    
    if(os.path.exists(dirpath)):
        rmtree(dirpath)
    path = Path(dirpath)
    path.mkdir(parents=True)#crear capetas 
    print("The directory ",dirpath," was created successfully")
    
    #t_maximo_local=[]
    
    #for x_i in range(len(id_clusters)):
    #t_maximo_local.append(np.max(extract_data_cluster(tot_datos,id_clusters,x_i)))
    t_maximo=np.max(tot_datos[:,0])
    
    for a_i in range(len(id_clusters)):
        
        fdata = open(str(dir_adr+str(math.trunc(id_clusters[a_i]))+'.txt'), 'w')
        print(str(dir_adr+str(math.trunc(id_clusters[a_i]))+'.txt'))
        try:

            mx_tiempo_id=extract_data_cluster(tot_datos,id_clusters,a_i)                      
                        
            fdata.write(str(math.trunc(float((t_maximo*1E+3)+25))) +'\n')
            fdata.write(str(len(mx_tiempo_id))+'\n')
                                
            for c_i in range(len(mx_tiempo_id)):
                fdata.write(str(math.trunc(float((mx_tiempo_id[c_i]*1E+3))))+'\n')
            fdata.write("0"+'\n')
                
        finally:
            fdata.close()

def extract_experimental_set(xdirectio_comple,Expr):            
    adres_dir = os.listdir(xdirectio_comple)
    n=len(adres_dir)
    
    rango = slice(7,25)
    xrango = slice(0,6)

    for i in range(n):
        if(adres_dir[i][xrango]=='id_sua'):     
            sub_cadena = adres_dir[i][rango]
            directorio = Expr+"/Exp_sua_"+str(sub_cadena.split(".")[0])+"/Exp_sua_"+str(sub_cadena.split(".")[0])
            addres = str(directorio)+"/sua_"+str(sub_cadena.split(".")[0])+"_"
            id_clusters= np.genfromtxt(xdirectio_comple+"/id_sua_"+str(sub_cadena.split(".")[0])+".txt")        
            tot_datos = (np.genfromtxt(xdirectio_comple+"/time_cluster_shank_sua_"+str(sub_cadena.split(".")[0])+".txt"))
            print(id_clusters)
            write_times_cluster (id_clusters,tot_datos,addres,directorio)

def extract_cluster_id(xdirectio_comple,experimento):
    xadres_dir = os.listdir(xdirectio_comple)
    nx=len(xadres_dir)
    mxrango = slice(23,40)
    nxrango = slice(0,22)
    for i in range(nx):
        if(xadres_dir[i][nxrango]=='time_cluster_shank_sua'):     
            zxaddres=xdirectio_comple+"/time_cluster_shank_sua_"+str(xadres_dir[i][mxrango].split(".")[0])+".txt"
            extract_clusters(zxaddres,xadres_dir[i][mxrango].split(".")[0],xdirectio_comple)
            
            
            

