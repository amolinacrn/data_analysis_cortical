import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import entropy
from libreriaMolCr import *
from MAIN import componente_gigante, componente_significativa,draw_comunidades

def D_js(P, Q):
    M = 0.5*(P+Q)
    KLDpm = entropy(P, M,base=2)
    KLDqm = entropy(Q, M,base=2)
    JSDpq = 0.5*(KLDpm + KLDqm)
    
    return JSDpq

def entropia_shanon(dats,dt,file_dir,addres_dir,xbins):

    mxrango = slice(14,22)

    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]
    MxRet=np.zeros((len(dt),len(dt)))
    restrop=[]; Mb=len(xbins[:-1]);du=Mb*[1/Mb];Nm=2*(Mb-1)/Mb
    for r in range(len(dt)):
        zdatx=(delete_nod(np.genfromtxt(fildir+str(int(dt[r]))+'.txt')))
        dtx=np.zeros((len(zdatx[:,0]),len(zdatx[:,0])))
        io,jo=np.where(zdatx>0); 
        dtx[io,jo]=zdatx[io,jo]#(zdatx[io,jo]-np.min(zdatx))/(np.max(zdatx)-np.min(zdatx))
        


        vector_k_lover, nkover, error_degree= mean_degree_network(dtx)
        
        # frec=[]      
        # for go in range(len(dtx[:,0])):
        #     xg=np.count_nonzero(dtx[go,:])
        #     frec.append(xg)
        
        frec, xbin=np.histogram(dtx,bins=xbins)
        x=xbin[:-1]; y= frec/sum(frec)

        frec, xbin=np.histogram(dtx,bins=np.linspace(0,1,30))
        zx=xbin[:-1]; zy= frec/sum(frec)
        
        #cls,lkp,long_path,error_lp, mean_cl, error_mean_cl = metricas_L_C(zdatx)
        



        # desimilid=1-D_js(y,du)
        #plt.plot(x,y,"-")
        
        restrop.append(entropy(y,base=2))
    return restrop

def define_entropia_shanon(adrespr,foldprincp,foldresult,interval_bins):
    for l,ndir in enumerate(adrespr):
        xdr=os.listdir(foldprincp+str(ndir))
        #print('datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
        for xsubdir in xdr:
            addres_dir=os.listdir(foldprincp+adrespr[l]+"/"+str(xsubdir))
            dir_file=foldprincp+adrespr[l]+"/"+str(xsubdir)       
            fdat=np.loadtxt(foldresult+str(xsubdir)+".txt")#datos empiricos
            thrho=0#np.mean(fdat[:,3])-2*np.std(fdat[:,3])
            print(xsubdir)#,"desidad: ",thrho);
            
            par_y=entropia_shanon(fdat,fdat[:,0],dir_file,addres_dir,interval_bins)     
            plt.plot(fdat[:,1],par_y,"^")
            
zinterval_bins=np.linspace(0,1,30)
for zadrespr in ["Mar07","Mar10","Mar04","Jan29","Jan21","Jan14"]:
    zfoldprincp="datW10T250s/"
    zfoldresult="restW10T250/"

    define_entropia_shanon([zadrespr],zfoldprincp,zfoldresult,zinterval_bins)
    # plt.yscale("log")
    # plt.xscale("log")
    plt.xlabel(r'$\langle CV \rangle$')
    plt.ylabel('Entropy')
    
    plt.show()