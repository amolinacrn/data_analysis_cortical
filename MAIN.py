from math import pi
import numpy as np
import os
from portrait_Djs.funct_PDjs import * 
from clasificar import *
from plotsGraf import *
from Jensen_Shannon import *

#PARA CALCULAR EL RETRATO DE UNA RED Y LA DIVERGENCIA DE RETRATO

adrespr =["Mar07"]#os.listdir('../principal/mcporcv')

bin_cv=0.030; pdensidad=1.1

cfile=True; fsNone=None

addir = os.listdir('datMetx')

nterv=np.linspace(1,36,100)

def filtrar(dt,file_dir,addres_dir,thsrho):

    mxrango = slice(14,22)

    fildir= file_dir+"/"+addres_dir[0].split(str(int(addres_dir[0][mxrango].split('.')[0]))+'.txt')[0]

    dtao=[]
    for r in range(len(dt)):
        datx=abs(delete_nod(np.genfromtxt(fildir+str(int(dt[r]))+'.txt')))
        xzx=nx.from_numpy_array(datx)
        dzx=nx.density(xzx)
        if(dzx > thsrho):
            dtao.append(dt[r])
    return dtao
'''

for l,ndir in enumerate(adrespr):
    xdr=os.listdir('mcporcv/'+str(ndir))
    data=np.loadtxt('datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
    #print('datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
    for xsubdir in ["50sMar07"]:#xdr:
        addres_dir=os.listdir('mcporcv/'+adrespr[l]+"/"+str(xsubdir))
        dir_file='mcporcv/'+adrespr[l]+"/"+str(xsubdir)       
        fdat=np.loadtxt("datMetx/"+str(xsubdir)+".txt")#datos empiricos
        thrho=np.mean(fdat[:,3])-2*np.std(fdat[:,3])
        #print(xsubdir,"desidad: ",thrho);
        #retrato_red(fdat,nterv,dir_file,addres_dir,thrho,pdensidad)    
        networkPDjs(fdat[:,0],nterv,dir_file,addres_dir,xsubdir,thrho,pdensidad,cfile)            
                 
#adrespr =[]#os.listdir('mcporcv')
bin_cv=0.030;klz=1
addir = os.listdir('datMetx')


interval=np.linspace(0,0.25,20)
adrespr=["Mar04"]
foldprincp="datW10T50s/"
foldresult="restW10T50/"
datos_Exper='datas/time_cluster_shank_sua_2020'

for l,ndir in enumerate(adrespr):
    xdr=os.listdir(foldprincp+str(ndir))
    data=np.loadtxt(datos_Exper+str(ndir)+'.txt')
    #print('datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
    for xsubdir in xdr:
        addres_dir=os.listdir(foldprincp+adrespr[l]+"/"+str(xsubdir))
        dir_file=foldprincp+adrespr[l]+"/"+str(xsubdir)       
        fdat=np.loadtxt(foldresult+str(xsubdir)+".txt")#datos empiricos
        thrho=0#np.mean(fdat[:,3])-2*np.std(fdat[:,3])
        print(xsubdir)#,"desidad: ",thrho);
        jensen_shannon(fdat[:,0],dir_file,addres_dir,xsubdir,interval,True)     
'''
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

def miplot_main(adrespr,folder_principal,folder_resultados,trdartos):
    # cmap = plt.colormaps["plasma"]
    # #fig, srx = plt.subplots()#1, 1, figsize=(15, 6))
        
    # umbral_min=0.01
    # #interval=np.arange(0.05,0.35,0.02)
    # par_red=1
    # xbinr=np.arange(0.0,1,0.02)
    #nNinter = dist_pesos(addres_dir,dir_file,10,interval,par_red,umbral_min,False,False)
  
    for l,ndir in enumerate(adrespr):
        xdr=os.listdir(folder_principal+str(ndir))
        #data=np.loadtxt('datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
        #print('datas/time_cluster_shank_sua_2020'+str(ndir)+'.txt')
        for xi,nsubdir in enumerate(["100sMar07"]):#enumerate(xdr):#["30sMar07","40sMar07","50sMar07"]):#,"20sMar07","30sMar07","40sMar07","50sMar07"]:#xdr:
            axdir=os.listdir(folder_principal+adrespr[l]+"/"+str(nsubdir))
            dir_file=folder_principal+adrespr[l]+"/"+str(nsubdir)
            xCV=(np.loadtxt(folder_resultados+str(nsubdir)+".txt"))#datos empiricos
            fdat=np.loadtxt(folder_resultados+str(nsubdir)+".txt")#datos empiricos
            thrho=np.mean(fdat[:,4])-2*np.std(fdat[:,4])

            #detao=filtrar(fdat[:,0],dir_file,axdir,thrho)        
            #nNodos,sub_cadena, densidad_de_red,nkover,swpx, long_path, clustering,arb_exp_min,moduld = calcular_metricas(addres_dir,dir_file,10)
            
            # densy=[]
            # for o in detao:
            #     for j,ss in enumerate(xCV[:,0]):
            #         if(int(o)==int(ss)):
            #             densy.append(xCV[int(j),4])
            # densy.sort()
                
            # print(nsubdir)

            #mxret=np.genfromtxt("retratostxt/"+str(nsubdir)+"_retrato.txt")#+nsubdir+".txt")
            #fdat=np.loadtxt("datMetx/"+str(nsubdir)+".txt")#datos empiricos
            
            #mxD=similitud_cv(xCV[:,1])
            
            # tips=np.zeros((len(mxret[0,:])+2,len(mxret[0,:])))

            # tips[0,:]=np.transpose(np.array(xCV[:,0])); 
            # tips[1,:]=np.transpose(np.array(xCV[:,1]))
            
            # for q in range(len(mxret[0,:])):
            #     tips[q+2,:]=mxret[q,:]

            # frec, xbin=np.histogram(xCV[:,5],bins=70)
            # xo=xbin[:-1]; yo= frec/sum(frec)
            # plt.bar(xo,yo, width=0.001,alpha=1,ec='w')
            # plt.show()
            #n, bins, patches = plt.hist(xCV[:,3][xCV[:,3]>0.76], bins = 30, density=1, facecolor='green', alpha=1)
            # frec, xbin=np.histogram(np.triu(xret).flatten()[np.triu(xret).flatten()>0],bins="auto")
            # x=xbin[:-1]; y= frec/sum(frec)
            #sns.displot(data=xCV[:,3],bins=20,palette=True,norm=Normalize(0, 1))

            # plt.xlabel('Smarts')
            # plt.ylabel('Probability')
            # plt.title('Histogram of IQ')
            # plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
            # plt.axis([40, 120, 0, 0.3])
            # # plt.grid(True)
            # plt.show()

            #tips=np.transpose(tips)
        
            # for oz,xcv in enumerate(xCV[:,1]):#range(len(fdat[:,0])):   
            #     nNcolor=(xcv-min(xCV[:,1]))/(max(xCV[:,1])-min(xCV[:,1]))
            #     kys=savgol_filter(mxret[oz,:],5,3)
            #     #fig, rs = plt.subplots(1,1 , figsize=(8, 8)) 
            #     xz[0].plot(mxD[oz,:],mxret[oz,:],"o", markersize=2,color=cmap(nNcolor),alpha=0.8)
                #sns.scatterplot(data=tips, x=tips[:,0], y=tips[:,tv+2], hue=tips[:,1],palette=cmap,legend=False)
                
            # norm = matplotlib.colors.Normalize(vmin=min(xCV[:,1]), vmax=max(xCV[:,1])) 
            # sxm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm) 
            # divider = make_axes_locatable(xz[0])
            # cax = divider.append_axes("right" ,size="4%", pad=0.04)
            # cbar=fig.colorbar(sxm, cax=cax)

            #fig.savefig('grar_retr/fig'+str(tv)+'.jpg')
            # xz[1].plot(mxret[tv,:],xCV[:,1],"o",linewidth=0.5,color=cmap(tv))
        
            #fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap=cmap),ax=ax, label="Normalized Thrust [a.u.]")

            # sns.histplot(data=comutxt.flatten(), bins=150,kde=True)
            #plt.plot(x,y,linewidth=0.5)
            #xz[0].plot(xo,yo,linewidth=0.5)
            #plt.show()

            #cmap = ListedColormap(["white","black"])#, "gold", "lawngreen", "lightseagreen","blue","magenta","cyan","green","red"])
            #cmap = ListedColormap([ "jet"])
            
            # viridis =cm.get_cmap('plasma')
            # newcolors = viridis(np.linspace(0, 1, 256))
            # pink = np.array([1])#}[248/256, 24/256, 148/256, 1])
            # newcolors[:0, :] = pink
            # newcmp = ListedColormap(newcolors)

            # plt.title("Connectivity matrix with 20% of the most significant connections")        
            #plt.imshow(mxret,cmap =newcmp)##si quiero en colores variado cambio cmap="jet"
            #xz[1].imshow(mxD,cmap =newcmp)##si quiero en colores variado cambio cmap="jet"
            #plt.subplot(111)
            # plt.subplots_adjust(bottom=0.1,left=-0.065 ,right=1, top=0.989)
            # plt.xlabel("")
            # plt.ylabel("")
            # cax0 = plt.axes([0.81, 0.1, 0.03, 0.82])
            # cax1 = plt.axes([0.81, 0.1, 0.03, 0.82])
            
            # plt.colorbar(cax=cax0)
            # plt.colorbar(cax=cax1)
            # plt.tight_layout()
            # plt.show() 
                    
            
            #////////////////////////////////////////////////////////////////#
            #//////////// SPACION PARA CALCULAR COMUNIDADES /////////////////#
            #////////////////////////////////////////////////////////////////#
            
            # MxRet=np.genfromtxt(trdartos+str(nsubdir)+".txt")
            # mxD=np.array(delete_nod(MxRet))
            # vd=np.array(np.triu(mxD).flatten())
            # umbral_min=np.min(vd[vd>0])
            # par_red=np.max(vd[vd>0])
            # umbral,xDaux=componente_gigante(mxD,umbral_min,par_red,0,False)
            # xmz=np.array(componente_significativa(mxD,umbral,0,False))
            # draw_comunidades(xmz,nsubdir,trdartos,True)

            
            # comutxt=np.genfromtxt("mxsimilt_w10n50/"+nsubdir+".txt")
            # print(nsubdir)
            # comuz=[]; comux=[]; comuy=[]; xCVs=[]; idCVs=[];comuw=[]; 

            # for i , j in enumerate(comutxt[:,0]):
            #     if(j==2):
            #         comux.append(int(comutxt[:,1][i]))
            #     if(j==1):
            #         comuy.append(int(comutxt[:,1][i]))
            #     if(j==0):
            #         comuz.append(int(comutxt[:,1][i]))
            #     if(j==3):
            #         comuw.append(int(comutxt[:,1][i]))

            # tx=np.concatenate((comuy,comuz,comux))
            # empqdata=[comuy,comuz,comux]      
            # dx=[];dy=[]       
            
            # for lista in empqdata:
            #     x,y=dist_critico(axdir,nsubdir,dir_file,xbins,comuz,comuy,comux,xCV,lista,comuw,False,False)
            #     dx.append(x)
            #     dy.append(y)
           
            #plt.xscale("log")
            # plt.yscale("log")
            # plt.show()
            
            #dist_psCV(axdir,nsubdir,dir_file,xbins,comuz,comuy,comux,xCV,tx,comuw,False,False)
            plot_distribucion_continua(axdir,dir_file,xbins,xCV,xCV[:,0])
            #plot_distribuciones_lognormal(axdir,nsubdir,dir_file,xbins,comux,comuy,comuz,fdat)
            #plot_distribuciones_clasificadas(axdir,dir_file,xbins,comuy,comuz,comux,xCV)

            # for tim in xCV[:,1]:
            #     print(tim)

    # plt.xscale("log")
    # plt.yscale("log")
    # plt.show()

    # fig, axx = plt.subplots(figsize=(10,5))
    # xCV=np.loadtxt("datMetx/30sMar07.txt")#datos empiricos
    # MxRet=np.genfromtxt("redesDjs/30sMar07rDjs.txt")
    # #sns.displot(data=MxRet[MxRet>0], kind="kde")
    # frec, xbin=np.histogram(MxRet[MxRet>0],bins=100)
    # xo=xbin[:-1]; yo= frec/sum(frec)
    # #plt.plot(xo,yo)
    # viridis = mpl.colormaps['viridis'].resampled(len(yo))
    # newcolors = viridis(np.linspace(0,1,len(xCV[:,1])))
    # newcmp = ListedColormap(newcolors)

    # for row in range(len(xo)):
    #     axx.bar(xo[row], yo[row], width=.007,ec='k',color=viridis.colors[row]) 
    # plt.tight_layout()


adrespr =["Mar07"]#,"Mar04","Mar10","Jan14","Jan29","Jan21","Dez20"]#os.listdir('expUnico')
folder_principal="datW10T100s/"
folder_resultados="restW10T100/"
trdartos="mxsimilt_w10n50/"
miplot_main(adrespr,folder_principal,folder_resultados,trdartos)
