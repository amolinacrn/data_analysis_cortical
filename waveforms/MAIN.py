from math import pi
import numpy as np
import os
import networkx as nx
from libreriaMolCr import delete_nod
#PARA CALCULAR EL RETRATO DE UNA RED Y LA DIVERGENCIA DE RETRATO


def mean_degree_network(mxc):
    # #En esta seccio se calcula el grado medio de la red

    xGz = nx.from_numpy_array(mxc)
    strength = dict(xGz.degree(weight='weight'))
    normstrengthlist=list(strength.values())
    mean_degree=np.mean(normstrengthlist)
    std_degree=np.std(normstrengthlist)/np.sqrt(len(normstrengthlist))
    return normstrengthlist, mean_degree, std_degree

def dist_critico(addres_dir,nsubdir,file_dir,interval,cv_alto,cv_medio,cv_bajo,mdat,idcv,cv_otras,comunidad_tres=False,comunidad_otras=False):
    mxrango = slice(14,22)
    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    xtwij=[]
    for i in range(len(idcv)):
        w_ij=[]
        xmcc = abs(delete_nod(np.genfromtxt(fildir+str(int(idcv[i]))+'.txt')))
        vector_k_lover, nkover, error_degree= mean_degree_network(xmcc)
        for kn in range(len(xmcc[:,0])):
            w_ij.append(xmcc[kn,:]/vector_k_lover[kn])
        W_ij=np.array(w_ij)
        xtwij+=list(W_ij.flatten())

    # frec, xbin=np.histogram(xtwij,bins=interval)
    # x=xbin[:-1]; y= frec/sum(frec)
    return xtwij

##////////////////////////////////////////////////////////
##  para crear toda clase de plots     
##////////////////////////////////////////////////////////
xbins=np.linspace(0,0.35,50)

def miplot_main():

    adrespr =["Mar04"]#,"Mar04","Mar10","Jan14","Jan29","Jan21","Dez20"]#os.listdir('expUnico')

    bin_cv=0.030
    #plot_weights_networks()

    folder_principal="datW10T50s/"
    folder_resultados="restW10T50/"

    for l,ndir in enumerate(adrespr):
        xdr=os.listdir(folder_principal+str(ndir))
        #data=np.loadtxt('datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
        #print('datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
        for xi,nsubdir in enumerate(["50sMar04"]):#enumerate(xdr):#["30sMar07","40sMar07","50sMar07"]):#,"20sMar07","30sMar07","40sMar07","50sMar07"]:#xdr:
            axdir=os.listdir(folder_principal+adrespr[l]+"/"+str(nsubdir))
            dir_file=folder_principal+adrespr[l]+"/"+str(nsubdir)
            xCV=(np.loadtxt(folder_resultados+str(nsubdir)+".txt"))#datos empiricos
            fdat=np.loadtxt(folder_resultados+str(nsubdir)+".txt")#datos empiricos
            thrho=np.mean(fdat[:,4])-2*np.std(fdat[:,4])
            #detao=filtrar(fdat[:,0],dir_file,axdir,thrho)        
            
            comutxt=np.genfromtxt("mxsimilt_w10n50/"+nsubdir+".txt")
            #print(nsubdir)
            comuz=[]; comux=[]; comuy=[]; xCVs=[]; idCVs=[];comuw=[]; 

            for i , j in enumerate(comutxt[:,0]):
                if(j==0):
                    comux.append(int(comutxt[:,1][i]))
                if(j==2):
                    comuy.append(int(comutxt[:,1][i]))
                if(j==1):
                    comuz.append(int(comutxt[:,1][i]))
                if(j==3):
                    comuw.append(int(comutxt[:,1][i]))

            tx=np.concatenate((comuz,comux))
            empqdata=[comuy,comuz,comux]      
            # dx=[];dy=[]           
            # for lista in empqdata:
            #     x=dist_critico(axdir,nsubdir,dir_file,xbins,comuz,comuy,comux,xCV,lista,comuw,False,False)
            #     #plt.plot(x[x>0.075],y[x>0.075],"^")
            #     dx.append(x)
            
            #dy.append(y)
            return empqdata


