from pty import spawn
import numpy as np
import os
import numpy 
import math
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from LibreriaMolCr import *
import matplotlib as mpl
from clasificar import componente_gigante

def simil_cos(xx,yy):
    
    cos_sim=np.zeros((len(xx),len(yy)))
    for i in range(len(yy)-1):
        va=[xx[i], yy[i]]
        norm_a=np.sqrt(np.dot(va,va))
        for j in range(i+1,len(yy)):
            vb=[xx[j], yy[j]]
            norm_b=np.sqrt(np.dot(vb,vb))
            pp=norm_a*norm_b
            rest=(np.dot(va,vb)/pp)
            cos_sim[i,j]=1-rest 
            cos_sim[j,i]=1-rest
    return cos_sim

def humbral(xdm):
    mxD=np.array(xdm)
    vd=np.array(np.triu(mxD).flatten())
    umbral_min=0.005+np.min(vd[vd>0])
    par_red=np.max(vd[vd>0])
    umbral,xDaux=componente_gigante(mxD,umbral_min,par_red,0,False)

    return umbral
def humbral_corr(cos_sim):
    xcosim=np.zeros((len(cos_sim[:,0]),len(cos_sim[0,:])))
    Hu=humbral(cos_sim)
    io,jo=np.where(cos_sim<Hu)
    xcosim[io,jo]=cos_sim[io,jo]

def simil_corr(adrespr,result,ndir,pltf):

    mxrango = slice(0,9)
    kpath="{}/{}".format(adrespr,ndir) 
    conjunto_datos=os.listdir(kpath)
    
    for dw,xsubdir in enumerate(conjunto_datos):
        kfile="{}/{}/{}".format(adrespr,ndir,xsubdir[mxrango])
        reswt="{}/{}.txt".format(result,xsubdir[mxrango])
        print(kfile)
        dt_cv=np.loadtxt(reswt)
        mc_files=os.listdir(kfile)
        #name_spike=mc_files[0].split('_')[0]+"_"+mc_files[0].split('_')[1]
        name_spike=mc_files[0].split('_')[0]+"_Exp"+"_"+mc_files[0].split('_')[2]
        xgrade=[]
        for _,u in enumerate(dt_cv[:,0]):
            load_file="{}/{}/{}/{}_{}.txt".format(adrespr,ndir,xsubdir[mxrango],name_spike,int(u))
            mspikes=np.array(abs(delete_nod(np.genfromtxt(load_file))))
            
            dv=np.linspace(np.min(mspikes[mspikes>0]),np.max(mspikes[mspikes>0]),101)
            xm=mspikes[mspikes>0]
            xfx, xxb=np.histogram(xm,bins=dv)
            sx=xxb[:-1]; sy= xfx/sum(xfx)
            
            mx=np.zeros((len(sx),len(sx))); xp=sx*sy
            for i in range(len(sx)-1):
                for j in range(i+1,len(sx)):
                    z=abs(xp[i]-xp[j])
                    mx[i,j]=z; mx[j,i]=z

            xgrade.append(path_lenght(mx))

        pltf.plot(dt_cv[:,1],xgrade,"o")
        # plt.imshow(mx)

            
dir_datos="datW10T250s"
resultwt="restW10T250"

fig = plt.figure(figsize=(13, 9))
gs = fig.add_gridspec(2,3)
nNo=[0,1,2,0,1,2]; r=0

params = {'xtick.labelsize': 20, 'ytick.labelsize': 20}
mpl.rcParams.update(params)

experimentos=["Mar07","Mar10","Mar04","Jan29","Jan21","Jan14"]
nome_exper=["SUAM07","SUAM10","SUAM04","SUAJ29","SUAJ21","SUAJ14"]
totales=[]; treshol_universal=[]
ct=[0.03,0];limit_Y=[0,.4]
for l,xndr in enumerate(experimentos):
    if(l>2):
        r=1
    o=nNo[l]

    azh = fig.add_subplot(gs[r,o])
    limit_X=simil_corr(dir_datos,resultwt,xndr,azh)

    # #configure_axes_canvas(o,azh,limit_Y,limit_X,nome_exper[l],ct)

plt.tight_layout()
plt.show()