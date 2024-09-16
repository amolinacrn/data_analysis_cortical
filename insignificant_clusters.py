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

def delete_insignificant_clusters(mxct,n_x):

    contar_elementos=[]
    
    elements_salida=[]
    
    for i in range(len(mxct)):
        
        elements=0
        
        for j in range(len(mxct)):
            if mxct[i][j]!=0:
               elements+=1
               
        contar_elementos.append(elements)   
        elements_salida.append("Fila: "+str(i+1)+": "+str(elements))
        
    copi_contar_elementos=contar_elementos.copy()  

    elementos_ordenados= sorted(contar_elementos) 
          
    mxr=[]
    
    print(elementos_ordenados,"\n")
    
    for i in range(n_x):
        for j in range(len(elementos_ordenados)):
            if(copi_contar_elementos[j]==elementos_ordenados[i]):
                mxr.append(elements_salida[j])
    
    result = []
    
    for item in mxr:
        if item not in result:
            result.append(item)
    for r in range(len(result)):
        print(result[r])

#def delete_unconnected_nodes(tmc,delet_fil,delet_col):
    
    
    