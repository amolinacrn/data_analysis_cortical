import numpy as np
import os
import numpy
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors
# Import required Python packages 
import statistics
import networkx as nx # version 2.4 
#import community # version 0.13 (python-louvain) 
#import gudhi # version 3.3.0 
import scipy.io # version 1.4.1 from sklearn 
#import preprocessing # version 0.23.1 
import itertools 
import seaborn as sns # version 0.11.0 
import pandas as pd
from scipy.interpolate import make_interp_spline
import wvf_datos


def wvf_maxima(dat_exp):

    #fig, [[ax1,ax2,ax3],[ax4,ax5,ax6]] = plt.subplots(3, 3, figsize=(12, 10))


    d_ch=200

    #canal = 0     1      2     3      4      5     6      7      8      9
               

    cx0 = [0.0,  32.0,  2.5,  29.0,  5.0,  27.0,  7.5,   24.5,  10.0,  22.0 ]

    chy   = [0.0,  20.0,  40.0, 60.0,  80.0, 100.0, 120.0, 140.0, 160.0, 180.0 ]

    cx1 = [0.0+d_ch,  32.0+d_ch,  2.5+d_ch,  29.0+d_ch,  5.0+d_ch,  27.0+d_ch,  7.5+d_ch, 24.5+d_ch, 10.0+d_ch, 22.0+d_ch]

 
    cx2 = [ 0.0+2*d_ch,32.0+2*d_ch,2.5+2*d_ch,29.0+2*d_ch,5.0+2*d_ch,27.0+2*d_ch,7.5+2*d_ch,24.5+2*d_ch,10.0+2*d_ch, 22.0+2*d_ch]


    cx3 = [ 0.0+3*d_ch,32.0+3*d_ch, 2.5+3*d_ch,29.0+3*d_ch, 5.0+3*d_ch,27.0+3*d_ch,7.5+3*d_ch,24.5+3*d_ch,10.0+3*d_ch,22.0+3*d_ch]


    cx4 = [0.0+4*d_ch,32.0+4*d_ch, 2.5+4*d_ch,29.0+4*d_ch,5.0+4*d_ch,  27.0+4*d_ch, 7.5+4*d_ch, 24.5+4*d_ch, 10.0+4*d_ch, 22.0+4*d_ch ]


    cx5 = [0.0+5*d_ch,32.0+5*d_ch,2.5+5*d_ch, 29.0+5*d_ch, 5.0+5*d_ch, 27.0+5*d_ch, 7.5+5*d_ch, 24.5+5*d_ch, 10.0+5*d_ch,22.0+5*d_ch ]


    X=[0]

    Y=[0]

    contador=0
    
    for k in range(len(X)-1):
        for l in range(k+1,len(X),1):
            distancia_Euclidiana=math.sqrt(  (X[k]-X[l])*(X[k]-X[l])   +    (Y[k]-Y[l])*(Y[k]-Y[l])  )
            #contador +=1
            #print(contador)
        
    max_wvf_sua=10*[0]

    file_ALL_txt="coordenadas_clusters/crsua"+dat_exp+".txt"
    
    f_fdata = open(file_ALL_txt, 'w')
    
    #f_fdata.write(str(crda_X)+'\t'+str(crda_Y)+'\t'+'g'+str()+'\n')

    try:

        for g_i in range(6):
        
            directorio_cluster_info="datos_info_csv/"+dat_exp+"/g0"
        
            id_cluster_sua=wvf_datos.cluster_Info_Csv(directorio_cluster_info,g_i+1)
            
            for i_d in id_cluster_sua:
                
                d_dir="datoswaveforms/waveforms_SUA/wf2020"+dat_exp+"/wvf_SUA_2020"+dat_exp+"_g"+str(g_i)+"/wvf_sua_"+dat_exp+"_g"+str(g_i)+"_"+str(i_d)+".txt"
                
                wvf=np.loadtxt(d_dir,dtype=float)
    
                xs=len(wvf)*[0]; ys=len(wvf)*[0]; xss=len(wvf)*[0]; am=len(wvf)*[0]; amm=10*[0]; xx=10*[0]; ax=6*[0]; var=0
    
                sum_wx=0
                sum_wy=0
    
                for i in range(10):
                    for j in range(32):         
                        xs[j]=wvf[j][i]
                        xss[j]=wvf[j][i]
                        ys[j]=0.5*j
                
                    wvf_min=abs(numpy.min(xs)-xs[0])
                    wvf_max=abs(numpy.max(xs)-xs[0])    
                
                    if wvf_max > wvf_min:
                        max_wvf_sua[i]=wvf_max
                    else:
                        max_wvf_sua[i]=wvf_min
                        
                    if(g_i==0):     
                        sum_wx+=max_wvf_sua[i]*cx0[i]
                        sum_wy+=max_wvf_sua[i]*chy[i]

                    if(g_i==1):     
                        sum_wx+=max_wvf_sua[i]*cx1[i]
                        sum_wy+=max_wvf_sua[i]*chy[i]

                    if(g_i==2):     
                        sum_wx+=max_wvf_sua[i]*cx2[i]
                        sum_wy+=max_wvf_sua[i]*chy[i]

                    if(g_i==3):     
                        sum_wx+=max_wvf_sua[i]*cx3[i]
                        sum_wy+=max_wvf_sua[i]*chy[i]

                    if(g_i==4):     
                        sum_wx+=max_wvf_sua[i]*cx4[i]
                        sum_wy+=max_wvf_sua[i]*chy[i]

                    if(g_i==5):     
                        sum_wx+=max_wvf_sua[i]*cx5[i]
                        sum_wy+=max_wvf_sua[i]*chy[i]

    
                Rcm0x_g0=sum_wx/sum(max_wvf_sua)
                Rcm0y_g0=sum_wy/sum(max_wvf_sua)
        
                f_fdata.write(str(Rcm0x_g0)+'\t'+str(Rcm0y_g0)+'\t'+str(g_i)+str(i_d)+'\n')
        
    finally:
        f_fdata.close()

for zdatp in ["Jan29"]: 
    wvf_maxima(zdatp)
'''

    d_ch=200

    #canal = 0     1      2     3      4      5     6      7      8      9
    fig, [[ax1,ax2,ax3],[ax4,ax5,ax6]] = plt.subplots(2, 3, figsize=(12, 10))

    cx0 = [0.0,  32.0,  2.5,  29.0,  5.0,  27.0,  7.5,   24.5,  10.0,  22.0 ]

    chy   = [0.0,  20.0,  40.0, 60.0,  80.0, 100.0, 120.0, 140.0, 160.0, 180.0 ]

    cx1 = [0.0+d_ch,  32.0+d_ch,  2.5+d_ch,  29.0+d_ch,  5.0+d_ch,  27.0+d_ch,  7.5+d_ch, 24.5+d_ch, 10.0+d_ch, 22.0+d_ch]


    cx2 = [ 0.0+2*d_ch,32.0+2*d_ch,2.5+2*d_ch,29.0+2*d_ch,5.0+2*d_ch,27.0+2*d_ch,7.5+2*d_ch,24.5+2*d_ch,10.0+2*d_ch, 22.0+2*d_ch]


    cx3 = [ 0.0+3*d_ch,32.0+3*d_ch, 2.5+3*d_ch,29.0+3*d_ch, 5.0+3*d_ch,27.0+3*d_ch,7.5+3*d_ch,24.5+3*d_ch,10.0+3*d_ch,22.0+3*d_ch]


    cx4 = [0.0+4*d_ch,32.0+4*d_ch, 2.5+4*d_ch,29.0+4*d_ch,5.0+4*d_ch,  27.0+4*d_ch, 7.5+4*d_ch, 24.5+4*d_ch, 10.0+4*d_ch, 22.0+4*d_ch ]


    cx5 = [0.0+5*d_ch,32.0+5*d_ch,2.5+5*d_ch, 29.0+5*d_ch, 5.0+5*d_ch, 27.0+5*d_ch, 7.5+5*d_ch, 24.5+5*d_ch, 10.0+5*d_ch,22.0+5*d_ch ]


    for ti in range(6):
        #mxcm=np.loadtxt("coordenadas_clusters/crsuag20M4.txt", dtype=float)
        mxcm=np.loadtxt("coordenadas_clusters/crsuag"+str(ti)+zdatp+".txt", dtype=float)

        m_xx=len(mxcm)*[0]
        m_yy=len(mxcm)*[0]
        cxx0=10*[0];cxx1=10*[0];cxx2=10*[0];cxx3=10*[0];cxx4=10*[0];cxx5=10*[0];chyy=10*[0]

        for a_i in range(len(mxcm)):
            m_xx[a_i]=mxcm[a_i][0]
            m_yy[a_i]=mxcm[a_i][1]

        if(ti==0):
            ax1.scatter(cx0,chy,color="blue",marker="o")
            ax1.scatter(m_xx,m_yy,color="red",marker="v")
        if(ti==1):
            ax2.scatter(cx1,chy,color="blue",marker="o")
            ax2.scatter(m_xx,m_yy,color="red",marker="v")
        if(ti==2):
            ax3.scatter(cx2,chy,color="blue",marker="o")
            ax3.scatter(m_xx,m_yy,color="red",marker="v")
        if(ti==3):
            ax4.scatter(cx3,chy,color="blue",marker="o")
            ax4.scatter(m_xx,m_yy,color="red",marker="v")
        if(ti==4):
            ax5.scatter(cx4,chy,color="blue",marker="o")
            ax5.scatter(m_xx,m_yy,color="red",marker="v")
        if(ti==5):
            ax6.scatter(cx5,chy,color="blue",marker="o")
            ax6.scatter(m_xx,m_yy,color="red",marker="v")
        
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Scatter Plot")
    plt.show()
    

    print(amp1,amp2)
    #x, y = np.loadtxt(xs,ys,unpack=True)
    a=ys[0]#abs(numpy.min(ys)-wvf[0][i])
    for m in range(32):
        am[m]=a

    amm[i]=a
    xx[i]=i
    model=make_interp_spline(xs, ys)
    model1=make_interp_spline(xs, yss)
    
    x=np.linspace(0,15,1000)
    y=model(x)
    xx=np.linspace(0,15,1000)
    yy=model(xx)

    
    
   
    if(i==0):
        ax=ax1
        ax.plot(x,y)
        ax.plot(xs,am)
    if(i==1):
        ax=ax2
        ax.plot(x,y)
        ax.plot(xs,am)
    if(i==2):
        ax=ax3
        ax.plot(x,y)
        ax.plot(xs,am)
    if(i==3):
        ax=ax4
        ax.plot(x,y)
        ax.plot(xs,am)

    if(i==4):
        ax=ax5
        ax.plot(x,y)
        ax.plot(xs,am)
    if(i==5):
        ax=ax6
        ax.plot(x,y)
        ax.plot(xs,am)
    if(i==6):
        ax=ax7
        ax.plot(x,y)
        ax.plot(xs,am)
    if(i==7):
        ax=ax8
        ax.plot(x,y)
        ax.plot(xs,am)
'''



    
