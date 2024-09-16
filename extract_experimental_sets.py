import numpy as np
import os
import math
from pathlib import Path 
from shutil import rmtree

def extract_clusters(dir):
    mx = np.genfromtxt("time_cluster_shank_sua_"+str(dir)+".txt")
    print(mx[:,1])
    #print(mx[1])
    result = []
    
    for item in mx[:,1]:
    #for item in mmx:
        if item not in result:
            result.append(item)
    result.sort()       
    for i in range(len(result)):
        print((math.trunc(result[i])))
    
    
        
    file_txt="id_sua_"+str(dir)+".txt"
            
    fdata = open(file_txt, 'w')

    try:

        for zz_i in range(len(result)):
            fdata.write( str(math.trunc(result[zz_i]))+"\n")
                            
    finally:
        fdata.close()
    

def extract_data_cluster(tot_datos,id_clusters,a_i):
   
    n = len(tot_datos[:,0])

    mx_tiempo_id = []
       
    for a_j in range(n):
        if (tot_datos[a_j][1]==id_clusters[a_i]):
            mx_tiempo_id.append(tot_datos[a_j][0])
    print(len(mx_tiempo_id)) 
    return mx_tiempo_id

def write_times_cluster(id_clusters,tot_datos,dir_adr,dirpath):
    
    if(os.path.exists(dirpath)):
        rmtree(dirpath)
    path = Path(dirpath)
    path.mkdir(parents=True)#crear capetas 
    print("The directory ",dirpath," was created successfully")

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

def extract_experimental_set():            
    adres_dir = os.listdir('../probando_estraer_datos')
    n=len(adres_dir)

    rango = slice(7,22)
    xrango = slice(0,6)

    for i in range(n):
        if(adres_dir[i][xrango]=='id_sua'):     
            sub_cadena = adres_dir[i][rango]
            directorio = "Experiments/Exp_sua_"+str(sub_cadena)+"/Exp_sua_"+str(sub_cadena)
            addres = str(directorio)+"/sua_"+str(sub_cadena)+"_"
            id_clusters= np.genfromtxt("id_sua_"+str(sub_cadena)+".txt")
            print(id_clusters)        
            tot_datos = np.genfromtxt("time_cluster_shank_sua_"+str(sub_cadena)+".txt")
            write_times_cluster (id_clusters,tot_datos,addres,directorio)
           

def extract_cluster_id():            
    xadres_dir = os.listdir('../probando_estraer_datos')
    nx=len(xadres_dir)

    mxrango = slice(23,38)
    nxrango = slice(0,22)
    for i in range(nx):
        if(xadres_dir[i][nxrango]=='time_cluster_shank_sua'):    
            print("time_cluster_shank_sua_"+str(xadres_dir[i][mxrango])+".txt")     
            extract_clusters(xadres_dir[i][mxrango])
            

