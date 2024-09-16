import numpy as np
from bct.utils import BCTParamError, get_rng
import bct.algorithms as reference

def func_randmio(mx,n,iter):
    '''
    esta funcion crea matrices nulas mediante un proceso de aleatorizacion, 
    manteniendo el mismo numero de conexionesy el grado, no 
    mantiene la distruibucion de pesos
    esta funcion llama a la funcion "randmio_und_connected"
    
    mx: es la matrix de concetividad para la cual deseamos crear los mdelos nulos
    n : es el numero de matrices nulas que deseamos
    iter: numero de iteraciones por matriz nula, como minimo iter = 100  
    '''
    #mx = np.loadtxt("mwc.txt",dtype=float)

    n_mc=len(mx)

    for i in range(n):

        mxramdom = reference.randmio_und_connected(mx,iter)
       
        tmc_ramdom = mxramdom[0]

        ffile_txt="include/swp/datrandom/rdm"+str(i)+".txt"
        
        print(ffile_txt)        
        
        ffdata = open(ffile_txt, 'w')
        
        try:           
            for z_i in range(n_mc):
                for z_j in range(n_mc):
                    ffdata.write( str(tmc_ramdom[z_i][z_j])+"\t")
                ffdata.write("\n")
        finally:
            ffdata.close()
        
        