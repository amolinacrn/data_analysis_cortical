import numpy as np
import math
import re
from pylab import *
import pandas as pd
import csv
from pathlib import Path
import os
import errno
from shutil import rmtree
import pathlib


def cluster_Info_Csv(id_cluster_info_csv,i):
    '''
    devuelve todos los clusters good
    '''
    contador=0
    n_id_sua=0
        
    dat_csv=str(id_cluster_info_csv)+str(i)+"_cluster_info.csv"
    with open(dat_csv) as fichero, open(dat_csv) as fichero_id:
        lectura_fichero_csv = csv.DictReader(fichero)
        for linea in lectura_fichero_csv:
            grupo =linea['group']
            if(grupo=="good"):
                contador+=1
                
        id_clusters_sua=contador*[0]
        lectura_fichero_id = csv.DictReader(fichero_id)
        
        for clinea in lectura_fichero_id:
            idclusters =clinea['id']
            grupo =clinea['group']
            if(grupo=="good"):
                #print(grupo)              
                id_clusters_sua[n_id_sua]=int(idclusters)
                n_id_sua+=1
        
    return id_clusters_sua       

#def id_extract_Waveforms_ALL(file_waveforms,file_ALL_txt):
def id_extract_Waveforms_ALL(mes_experimentos,n_shanks):
    '''
    retorna todos los datos correspondientes a las formas de onda, limpios
    '''
   
    for g_o in range(n_shanks):
  
        file_waveforms="datoswaveforms/waveforms_mean_ALL/waveforms_mean_ALL_2020"+str(mes_experimentos)+"_g"+str(g_o)+".txt"#datos crudos
        file_ALL_txt="datoswaveforms/waveforms_ALL/wvf_ALL_"+str(mes_experimentos)+"_g"+str(g_o)+".txt"
        #resultado_all=id_extract_Waveforms_ALL(fdata_waveforms,file_all_txt) #CARGAR FORMAS DE ONDA
    
    
        f_fdata = open(file_ALL_txt, 'w')
    
        try:
            #print("ch_0, ch_1, ch_2, ch_3, ch_4, ch_5, ch_6, ch_7, ch_8, ch_9"+"\n")
            with open(file_waveforms,'r', encoding='utf-8') as f:
                my_str = f.read()
        
                #for i_d in  range(0,len(id_clusters_sua)-69,1):
                #m_idetif_t=str(re.findall(str(i_d)+":|array\(", my_str))
        
                m_del_t= re.sub("\{|\}","", my_str)
                m_del_t = re.sub("\[\[|\]\],","", m_del_t)
                m_del_t= re.sub("\n","", m_del_t)
                #m_del_t= re.sub(":",",#,#,#,#,#,#,#,#,#,#\n", m_del_t)
                m_del_t= re.sub(":",",0,0,0,0,0,0,0,0,0,0\n", m_del_t)
                m_del_t= re.sub("\],","\n", m_del_t)
                m_del_t= re.sub("\[","0,", m_del_t)
                m_del_t = re.sub("array\(","0,",m_del_t)
                m_del_t= re.sub("dtype=float32\),","\n", m_del_t)
                m_del_t= re.sub("dtype=float32\)","\n", m_del_t)
                m_del_t= re.sub("        ","", m_del_t)
                m_del_t= re.sub(" |       ","", m_del_t)
                m_del_t= re.sub(",","  ", m_del_t)
                #print("id_cluster,channel_0,channel_1,channel_2,channel_3,channel_4,channel_5,channel_6,channel_7,channel_8,channel_9")
                print(m_del_t)
                #m_del_t = re.sub("\{|array\(\[|\[|\],|\]\],|dtype=float32\),|\n|\}","", my_str)
                #m_del_t = re.sub("\{|array\(\[|\[|\],|dtype=float32\),|\n|\}","", my_str)
                #m_del_t = re.sub("\{|array\(|dtype=float32\)|\}","", my_str)
        
                f_fdata.write(str(m_del_t))
        
        finally:
            f_fdata.close()


def id_Cluster_Waveforms_SUA(el_mes):
    
    '''
    retorna datos correspondientes a las formas de onda,  clusters good
    '''
    for g_i in range(6):
    
        dir_clusters_sua="datoswaveforms/waveforms_SUA/wf2020"+el_mes+"/wvf_SUA_2020"+str(el_mes)+"_g"+str(g_i)

        ffile_all="datoswaveforms/waveforms_ALL/wvf_ALL_"+str(el_mes)+"_g"+str(g_i)+".txt"
        
        # rmtree(dir_clusters_sua)#eliminar carpetas
        # path = Path(dir_clusters_sua)
        # path.mkdir(parents=True)#crear capetas
        
        if(os.path.isdir(dir_clusters_sua)==False):
            path = Path(dir_clusters_sua)
            path.mkdir(parents=True)

        directorio_cluster_info="datos_info_csv/"+str(el_mes)+"/g0"
        # for i in range(1,7,1):
      
        id_cluster_sua=cluster_Info_Csv(directorio_cluster_info,g_i+1)
        
        
        id_datos = np.loadtxt(ffile_all,dtype=float)

        n_dc=len(id_datos)

        wvfj,wvf0,wvf1,wvf2,wvf3,wvf4,wvf5,wvf6,wvf7,wvf8,wvf9,wvf10=n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0],n_dc*[0]

        #for i in range(1,12,1):
        for j in range(len(id_datos)):
            wvfj[j]=j
            wvf0[j]=id_datos[j][0]
            wvf1[j]=id_datos[j][1]
            wvf2[j]=id_datos[j][2]
            wvf3[j]=id_datos[j][3]
            wvf4[j]=id_datos[j][4]
            wvf5[j]=id_datos[j][5]
            wvf6[j]=id_datos[j][6]
            wvf7[j]=id_datos[j][7]
            wvf8[j]=id_datos[j][8]
            wvf9[j]=id_datos[j][9]
            wvf10[j]=id_datos[j][10]

        o=0
        
        for a_l in id_cluster_sua:
        
            fil_SUA_txt="datoswaveforms/waveforms_SUA/wf2020"+str(el_mes)+"/wvf_SUA_2020"+str(el_mes)+"_g"+str(g_i)+"/wvf_sua_"+str(el_mes)+"_g"+str(g_i)+"_"+str(a_l)+".txt"
        
            ffdata_sua = open(fil_SUA_txt, 'w')
            
            try:
                
                if(a_l == wvf0[o]):
                
                    for o in range(o+1,o+33,1):
                                                    
                        ffdata_sua.write(#str(wvfj[o])+"\t"+
                                         str(wvf1[o])+'\t'+
                                         str(wvf2[o])+'\t'+
                                         str(wvf3[o])+'\t'+
                                         str(wvf4[o])+'\t'+
                                         str(wvf5[o])+'\t'+
                                         str(wvf6[o])+'\t'+
                                         str(wvf7[o])+'\t'+
                                         str(wvf8[o])+'\t'+
                                         str(wvf9[o])+'\t'+
                                         str(wvf10[o])+'\t'+"\n")
                        o+=1
                    
                
            finally:
                ffdata_sua.close()
        


#aqui se cambia el mes del experimento
mes_experimentos='Jan29'
n_shanks=6
#creando directorio y archivos de datos "semicrudos"
dir_clus_t_all="datoswaveforms/waveforms_ALL/"
rmtree(dir_clus_t_all)#eliminar carpetas
path = Path(dir_clus_t_all)
path.mkdir(parents=True)#crear capetas  

#cargar y arreglar datos formas de onda.
resultado_all=id_extract_Waveforms_ALL(mes_experimentos,n_shanks) #CARGAR FORMAS DE ONDA


            
#separamos formas de onda por cluster y las guardamos en un fichero .txt
resultado_sua= id_Cluster_Waveforms_SUA(mes_experimentos)


        



