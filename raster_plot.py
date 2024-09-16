import numpy as np

#@title binarize time stamps
def binary_array_from_stamps(x, N, nSteps, binSize):
    """
    Creates a binary array out of timestamps for a chosen bin size

    Args:
        x(numpy array of integers)            : timeStamps in (t,neuron_j) form
        N(scalar)                             : number of variables
        nSteps(scalar)                        : time series length
        binSize(scalar)                       : amount of time in each array slot (ms)

    Returns:
        binary_array(numpy array of integers) : (N,nbins) array with 1's on spike times t for neuron N
    """
    nbins = int(nSteps/binSize)
    binary_array = np.zeros((N,nbins),dtype=int)
    for i,j in x:
        timeSlot = int(i/binSize)
        if timeSlot<len(binary_array[0,:]):
            binary_array[j,timeSlot] = 1

    return binary_array

#@title rasterplots by CV

def get_CV(x):
    """
    Computes the Coefficient of Variation (CV) of the time series' summed activity
    Args:
        x (numpy array)                                 : data (N,t)
    Returns:
        CV (scalar)                                     : coefficient of variation
    """
    total_activity = np.sum(x,axis=0)
    CV = np.std(total_activity)/np.mean(total_activity)

    return CV


#@title raster plots (fig 3)
def plot_rasterCV(dataFile, binSize, trialLength):
    sliceSize = int(trialLength/binSize) #trial duration in number of array slots
    nTrials = int(len(dataFile[0,:])/sliceSize)
    dt = binSize*0.001

    dats = [[] for i in range(nTrials)]
    CVs = np.zeros((nTrials))
    for i in range(nTrials):
        next_slice = slice((i*sliceSize),(i+1)*sliceSize)
        dats[i] = dataFile[:,next_slice]
        CVs[i] = get_CV(dats[i])

    idx= [np.argmin(CVs), 180, np.nanargmax(CVs)]
    print(idx)

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(131)
    color = plt.cm.magma(np.linspace(0, 1, 12))
    aux = np.ones((len(dats[idx[0]][0,1:])))
    t = np.arange(0.05,trialLength*0.001,dt)
    print(np.shape(t))
    for i in range(len(dats[idx[0]][:,0])):
        x = dats[idx[0]][i,1:]*t
        l1, = ax.plot(x, i*aux, linestyle = 'none', marker = ".", color = color[0], markersize=2)
    plotname = 'CV = %.2f' %(CVs[idx[0]])
    plt.title(plotname,fontsize=18)
    ax.set_xticks([])
    ax.axis('off')

    ax = fig.add_subplot(132)
    for i in range(len(dats[idx[1]][:,0])):
        x = dats[idx[1]][i,1:]*t
        l1, = ax.plot(x, i*aux, linestyle = 'none', marker = ".", color = color[4], markersize=2)
    plotname = 'CV = %.2f' %(CVs[idx[1]])
    plt.title(plotname,fontsize=18)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.axis('off')
    fig.subplots_adjust(wspace=0, hspace=0)

    ax = fig.add_subplot(133)
    for i in range(len(dats[idx[2]][:,0])):
        x = dats[idx[2]][i,1:]*t
        l1, = ax.plot(x, i*aux, linestyle = 'none', marker = ".", color = color[8], markersize=2)
    plotname = 'CV = %.2f' %(CVs[idx[2]])
    plt.title(plotname,fontsize=18)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.axis('off')
    fig.subplots_adjust(wspace=0, hspace=0)

trialLength = 30*1000   #ms
binSize = 50
plot_rasterCV(binary_array, binSize, trialLength)


def raster(spk,clu):
  clu_id=np.unique(clu)
  spk_total=[]
  
  for index,k in enumerate(clu_id):
    index_clu=np.argwhere(clu==k)
    spk_total.append(spk[index_clu].flatten())
  return spk_total,clu_id

def matriz_spk(fildat,bin_option):
  spk=fildat[:,0]
  clu=fildat[:,1]
  
  spk_total,_=raster(spk,clu)# matriz em que cada linha é um cluster e cada coluna é um spike. Mas cada linha tem um tamanho diferente.
  ####Definir o tamanho que a matriz binarizada precisa ter para abarcar todos os spikes.############################
  spk_concatenate=[]
  for i in range(len(spk_total)):
    spk_concatenate.extend(spk_total[i])
  
  spk_concatenate=np.sort(spk_concatenate)
  interval=int((spk_concatenate[-1]-spk_concatenate[0])/bin_option)
  return spk_total,interval,spk_concatenate

#### Binarizar datps ##############################
def fun_binarizar(xfil,xbin,f):
  spk_total,interval,spk_concatenate=matriz_spk(xfil,xbin)
  spk_total_binary=np.zeros((len(spk_total),(interval+1)))
  for i in range(len(spk_total)):
    for j in range(interval+1):
      count_neuro=np.argwhere((spk_total[i]>=j*xbin+spk_concatenate[0]) & (spk_total[i]<=((j+1)*xbin)+spk_concatenate[0])) 
      if len(count_neuro)>=1:
        spk_total_binary[i,j]=1  
        print(spk_total_binary[i,j],i,j)

  np.save(f,spk_total_binary)
              
# def plot_raster_sort(raster,sort_per_minute=True):

#   if sort_per_minute:
#     pop_FR=raster.sum(axis=0)
#     desc_FR_neu=np.argsort(raster.sum(axis=1))
#     raster=raster[desc_FR_neu,:]

#   f=np.argwhere(raster!=0)

#   neuro=f[:,0]
#   spike=f[:,1]
#   return neuro,spike


nbins=0.05
f= "raster_Jan21_12_{}ms".format(int(nbins*5000))
xfil="datos_cv/Exp_250s/Exp_Jan21/time_cluster_shank_sua_2020Exp_Jan21_12.txt"
datostxt=(np.loadtxt(xfil))
#fun_binarizar(datostxt,nbins,f)


import matplotlib.pyplot as plt

filebin=np.load('raster_Jan21_12_50ms.npy')
for i in range(len(filebin)):
  for j in range(len(filebin)):
    print(np.shape(filebin))