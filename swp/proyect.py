import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import cm
from pylab import *
import os

def func_simul():

    adres_dir = os.listdir('datrandom')
    n=len(adres_dir)

    cmap = ListedColormap(["white","black"])
    viridis =cm.get_cmap('plasma')
    newcolors = viridis(np.linspace(0, 1, 256))
    pink = np.array([0])#[248/256, 24/256, 148/256, 1])
    newcolors[:0, :] = pink
    newcmp = ListedColormap(newcolors)

    address_mc="datrandom/rdm0.txt"
    address_mc=np.loadtxt(address_mc,dtype=float) 
    grafmx=imshow(address_mc,cmap=viridis)

    for i in range(n):
        fig, (ax) = plt.subplots(1, 1, figsize=(5, 5))
        address_mc="datrandom/rdm"+str(i)+".txt"
        address_mc=np.loadtxt(address_mc,dtype=float)
        plt.imshow(address_mc,cmap =newcmp)##si quiero en colores variado cambio cmap="jet"
        fig.savefig('graf/file'+str(i)+'.jpg')
        '''
        fig, (ax) = plt.subplots(1, 1, figsize=(5, 5))
        address_mc="datrandom/rdm"+str(i)+".txt"
        address_mc=np.loadtxt(address_mc,dtype=float) 
        grafmx.set_data(address_mc)
        fig.savefig('graf/file'+str(i)+'.jpg')
        plt.pause(1e-30)
        draw()
        '''

def graf_mc_conectividad(mx_conectividad):
    fig, (ax) = plt.subplots(1, 1, figsize=(5, 5))
    viridis =cm.get_cmap('plasma')
    newcolors = viridis(np.linspace(0, 1, 256))
    pink = np.array([1])#}[248/256, 24/256, 148/256, 1])
    newcolors[:1, :] = pink
    newcmp = ListedColormap(newcolors)

    plt.title("Connectivity matrix with 20% of the most significant connections")        
    plt.imshow(mx_conectividad,cmap =newcmp)##si quiero en colores variado cambio cmap="jet"
    #plt.subplot(111)
    #plt.subplots_adjust(bottom=0.1,left=-0.065 ,right=1, top=0.989)
    #plt.xlabel("")
    #plt.ylabel("")
    fig.savefig('graf/file.jpg')
    cax = plt.axes([0.81, 0.1, 0.03, 0.82])
    
    plt.colorbar(cax=cax)
    #plt.tight_layout()
    plt.show()


address_mc=np.loadtxt("datrandom/rdm0.txt",dtype=float)
adr_mc=np.loadtxt("mwc.txt",dtype=float)
#graf_mc_conectividad(adr_mc)
#graf_mc_conectividad(address_mc)
func_simul()
   
     
  
