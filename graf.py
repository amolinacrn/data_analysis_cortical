import numpy as np
import os
import numpy 
import math
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import seaborn as sns
from include.get_function_implementations import *
import networkx as nx
import networkx.algorithms.community as nx_comm
from community import community_louvain
import matplotlib as mpl
from matplotlib import cbook
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
#adres_dir = os.listdir('ccExp_sua_2020Mar4w20')  

addres_dir = os.listdir('Correlation/Resultados_07')
#fig, ax1 = plt.subplots(1, 1, figsize=(8, 8))  

nNfiles=len(addres_dir)

daots=['alto','medio','bajo']
o=0 ; r=0; 

valores_max=[]

interval=np.arange(0,0.3,0.01)

mxz=np.zeros((nNfiles,nNfiles))

mxrango = slice(14,22)

   
#x2 = np.sort(mxc)
#f2 = np.array(range(len(mxc)))/float(len(mxc))
#plt.plot(x2, 1-f2,label=addres_dir[i],linewidth=2)
#sns.kdeplot(data=(mxc))#,label=addres_dir[i]) 

#plt.yscale("log")
#plt.xscale("log")
#plt.ylabel('cumulative weight distribution')
#plt.xlabel('weights')
#plt.legend()

# if(o==2):
#     plt.savefig('fig'+addres_dir[i]+'_.jpg')
#     plt.show()  
#     o=-1   
#     r+=1
# o += 1


mxLinks = np.triu(np.genfromtxt("Exp_sua_2020Mar4_CC_Symmetric.txt"))    
mxCC = np.triu(np.genfromtxt("mcc_Exp_sua_2020Mar4.txt"))
fmxCC = np.triu(np.genfromtxt("Links_inhib_Exp_sua_2020Mar4.txt"))

mxc=[]    
    
for k in range(len(mxCC)):
    for m in range(len(mxCC)):
        if(mxCC[k,m]!=0):
            mxLinks[k,m]=0
      
premxc=[]    
    
for k in range(len(mxLinks)):
    for m in range(len(mxLinks)):
        if(mxLinks[k,m]!=0):
            premxc.append(str(fmxCC[k,m]))

for i in range(len(premxc)):  
    
    for j in range(len(adres_dir)):   
        
        if ("rsurr_"+str(int(premxc[i].split('.')[0]))+"_"+str(int(premxc[i].split('.')[1]))+".txt"==str(adres_dir[j])):
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))  
            xs, ys = np.loadtxt("ccExp_sua_2020Mar4w20/"+adres_dir[j],skiprows=0, usecols=[0,1],unpack=True)#datos empiricos
            for l in range(1,2): 
                #x_s, y_s = np.loadtxt("surr_"+str(l)+"/"+addres_dir[j],skiprows=0, usecols=[0,1],unpack=True)#datos empiricos
            
                #ax1.plot(mxc,myc,"^")
                #ax1.plot(xs,yy_s,label="Empirical data,exitatoria")    
                ax.plot(xs,ys)#,"-",x_s,y_s,label=str(adres_dir[j]))
                ax.legend()
                fig.savefig('graf/fig'+str(int(premxc[i].split('.')[0]))+"_"+str(int(premxc[i].split('.')[1]))+'_.jpg')
            
            #fig.savefig('graf/fig'+str((premxc[i].split('.')[0]))+'_'+str((premxc[i].split('.')[1]))+'_.jpg')
            
# #plt.show()
    frec, xbin=np.histogram(mxc,bins=interval)
    x=xbin[:-1]; y= frec/sum(frec)
    #print(len(x),' ',len(y))
    for ox in range(len(y)):
        mxz[i,ox]=y[ox]+0.001

    plt.plot(x,y,linewidth=0.5)