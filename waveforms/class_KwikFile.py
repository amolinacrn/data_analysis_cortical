#@title Importar pacotes
import matplotlib.pyplot as plt
import numpy as np
#%matplotlib inline

import math
import numpy.matlib
import scipy.optimize as opt
import random

from klusta.kwik import KwikModel
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import os
from pathlib import Path
from shutil import rmtree
#%matplotlib inline



plt.rcParams["figure.figsize"] = [10,6]
plt.rcParams.update({'font.size': 18})

#@title class KwikFile
class KwikFile:
    """!  @brief Model for Kwik file, strongly based on KwikModel from phy project

    The main purpose of this class is provide an abstraction for kwik files provided by phy project. The current version contains a basic set of fundamental methods used in kwik file      
    @author: Nivaldo A P de Vasconcelos
    @date: 2018.Feb.02
    """
    #get_path
    def __init__(self,kpath=None,name=None):  
      self.kwik_model=None
      self.name = name
      self.kpath=None
      if (kpath is not None):
          self.kwik_model=KwikModel(kpath)
          self.kpath=kpath
          if (name is None):
              self.name = self.kwik_model.name
          #print ("Created class on = %s !" % kpath)
      else:
          print ("It still with no path:(")
    
    def get_name(self):
        """! @brief Returns the found in name field in kwik file.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        return (self.name)
    def set_kwik_file (self,kpath):
        """! @brief Defines the corresponding kwik file

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        self.kwik_model=KwikModel(kpath)
        self.name = self.kwik_model.name
        self.kpath=kpath
    def sampling_rate (self):
        """! @brief Returns the sampling rate used during the recordings 

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        return (self.kwik_model.sample_rate)
    def shank (self):
        """! @brief Returns the shank/population's id used to group the recordings.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        return (self.kwik_model.name)
    def get_spike_samples (self):
        """! @brief Returns the spike's samples on the recordings.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """

        return (self.kwik_model.spike_samples)
    def get_spike_clusters (self):
        """! @brief Returns the corresponding spike's clusters on the recordings.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        return (self.kwik_model.spike_clusters)
    def describe(self):
        """! @brief Describes the kwik file

        It calls the describe method in KwikModel

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        self.kwik_model.describe()
    def close (self):
        """! @brief Closes the corresponding kwik model

        It calls the close method in KwikModel

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        self.kwik_model.close()
    def list_of_groups (self):
        """! @brief Returns the list of groups found in kwik file

        The result has a list's form.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """

        lgroups = list(self.groups().values())
        lgroups = list(set(lgroups))

        return (lgroups)
    def list_of_non_noisy_groups (self):

        """! @brief Returns the list of groups found in kwik file which are not called noise

        The result has a list's form.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """

        lgroups = list(self.groups().values())
        lgroups = list(set(lgroups)-set(['noise',]))
        return (lgroups)
    def all_clusters (self):
        """! @brief Returns the list of all clusters in kwik file

        The result has a list's form.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """

        llabels = list(self.groups().keys())
        llabels = list(set(llabels))

        return (llabels)
    def groups(self):
        """! @brief Returns a dict with cluster label and its respective group

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        if not(isinstance(self.kwik_model,KwikModel)):
            raise ValueError("There is no KwikModel assigned for this object.")
        return (self.kwik_model.cluster_groups)
    

    def data_info_clusters(self,group_name=None):
        """! 
        esta funcion me devuelve el nombre del cluster, si es sua o mua
        y el id de cada cluster.
        """
        sgrj=[]
        if (group_name is None):
            return (self.all_clusters())
 
        if not(group_name in self.list_of_groups()):
            raise ValueError("\nThis group was not found in kwik file: %s\n" % group_name)
        group=self.groups()
        return group

    def clusters (self,group_name=None):
        """! @brief Returns the list of clusters on kwik file

        It can be used to get the list of clusters for a given group by pproviding
        this information the group_name.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02
        """
        if (group_name is None):
            return (self.all_clusters())
 
        if not(group_name in self.list_of_groups()):
            raise ValueError("\nThis group was not found in kwik file: %s\n" % group_name)
        group=self.groups()
        
        clusters=[]
        for c in self.all_clusters():
            if (group[c]==group_name):
                clusters.append(c)
        clusters.sort()
        return (clusters)
    print
    def get_cluster_waveforms (self,cluster_id):


        try:
            if (not(type(self.kwik_model) is KwikModel)):
                raise ValueError       
        except ValueError:
                print ("Exception: internal Kwik Model is not valid !!")
                return
        
        clusters = self.kwik_model.spike_clusters
        try:
            if ((not(cluster_id in clusters))):
                raise ValueError       
        except ValueError:
                print ("Exception: cluster_id (%d) not found !! " % cluster_id)
                return

        w=self.kwik_model.all_waveforms[clusters==cluster_id]
        return w
    def plot_cluster_waveforms (self,cluster_id,nspikes, save=False,save_path=None, wave_color="gray"):
    
    
        w = self.get_cluster_waveforms (cluster_id)
        model=self.kwik_model
        y_scale = .7
        x_scale = 4
        num_channels = w.shape[2]
        waveform_size = w.shape[1]
        np.random.seed(int(time.time()))
        
        fig=plt.figure(num=None, figsize=(6, 18), dpi=50, facecolor='w', edgecolor='k')
        plt.clf()
        spike_id = np.arange(w.shape[0])
        np.random.shuffle(spike_id)
        spike_id = spike_id[0:nspikes]
        posx=np.flipud (model.channel_positions [:,0])
        posy=np.flipud (model.channel_positions [:,1])
        for ch in range (0,num_channels):
            x_offset = model.channel_positions [ch,0]
            y_offset =posy[ch]*y_scale
            #y_offset = model.channel_positions [ch,1]*y_scale
            mu_spikes = np.mean(w[:,:,ch],0)
            for i in spike_id:
                spike = w[i,:,ch]
                x=x_scale*x_offset+range(0,waveform_size)
                plt.plot (x,0.05*spike+y_offset,color=wave_color,alpha=0.5)
            plt.plot (x,0.05*mu_spikes+y_offset,"--",color="black",linewidth=3,alpha=0.3)
        plt.tight_layout()
        plt.axis("off")
        plt.show()

        if (save):
            if (save_path):
                filename = "%s/%s_waveform_%02d.pdf" % (save_path,self.kwik_model.name,cluster_id)
            else:
                filename = "waveform_%02d.pdf" % cluster_id
            fig.savefig (filename)
    def all_spikes_on_groups (self,group_names):
        """! @bri\ef Returns the all spike samples within a list of groups

        Usually the clusters are organized in groups. Ex: noise, mua, sua,
        unsorted This method returns, in a single list of spike samples, all
        spikes found in a lists of groups (group_names). 

        Parameters:
        group_names: list of group names, where the spikes will be searched. 

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02 
        """

        spikes=[]
        all_spikes=self.get_spike_samples()
        all_labels=self.get_spike_clusters()
        
        if not(isinstance(group_names,list)):
           raise ValueError("\nThe argument must be a list.") 
        
        for group_name in group_names:
            if not(group_name in self.list_of_groups()):
                raise ValueError("\nThis group was not found in kwik file: %s\n" % group_name)
            for c in self.clusters(group_name=group_name):
                spikes=spikes+list(all_spikes[all_labels==c])
        spikes.sort()
        return (spikes)
    def group_firing_rate (self,group_names=None,a=None,b=None): 
        """! @brief Returns firing rate in a given set of groups found in kwik file.

        Usually, the clusters are organized in groups. Ex: noise, mua, sua,
        unsorted. This method returns, in a doubled dictionary, the firing rate
        for each cluster, organized by groups.

        Parameters: 
        group_names: list of group names, where the spikes will be
        searched. When this input is 'None' all groups are taken. The resulting
        dictionary has the first keys as groups, and the second keys as the
        respective cluster id's, whereas the value, is the corresponding firing
        rate within [a,b].


        Please refer to the method cluster_firing_rate in order to get more 
        details about the firing calculation.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02 
        """
        if not(isinstance(group_names,list)) and not(group_names is None):
           raise ValueError("\nThe argument must be a list or a None.") 
        spk=dict()
        if group_names is None:
            group_names=self.list_of_non_noisy_groups()
        for group_name in group_names:
            if not(group_name in self.list_of_groups()):
                raise ValueError("\nThis group was not found in kwik file: %s\n" % group_name)
            spk[group_name]=dict()
            for c in self.clusters(group_name=group_name):
                spk[group_name][c]=self.cluster_firing_rate(c,a=a,b=b)
        return (spk)
    def cluster_firing_rate (self,cluster_id,a=None,b=None):
        """! @brief Returns firing rate in a given cluster_id found in the kwik file

        In the kwik file, a cluster stores the spike times sorted for a given neuronal
        unit. The firing rate here is calculated by dividing the number of spike times
        by the number of seconds of the time period defined by [a,b]. 
        If a is 'None' a is assingned to zero; if b is 'None', it is assigned to the time
        of the last spike within the cluster.

        Parameters:
        cluster_id: id which identifies the cluster.
        a,b: limits of the time period where the firing rate must be calculated.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02 
        """
        sr=self.sampling_rate()
        spikes=np.array(self.spikes_on_cluster (cluster_id))/sr
        if a is None:
            a=0
        if b is None:
            b=spikes[-1]
        if (a==b):
            raise ValueError ("\nThe limits of the time interval are equal\n")
        piece=spikes[(spikes>=a)]
        piece=piece[piece<=b]
        return (len(piece)/(b-a))
    def spikes_on_cluster (self,cluster_id):
        """! @brief Returns the all spike samples within a single cluster

        Parameters:
        cluster_id: id used to indentify the cluster.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02 
        """

        if not(cluster_id in self.all_clusters()):
            raise ValueError("\nThis cluster was not found in kwik file: %s\n" % cluster_id)
        all_spikes=self.get_spike_samples()
        all_labels=self.get_spike_clusters()
        spikes=list(all_spikes[all_labels==cluster_id])
        spikes.sort()
        return (spikes)
    def group_firing_rate_to_dataframe (self,group_names=None,a=None,b=None):
        """! @brief Exports the group's firing rate into a pandas dataframe


        Usually, the clusters are organized in groups. Ex: noise, mua, sua,
        unsorted. This method returns, in a pandas dataframe, which contains the
        following information for each unit: 'shank', 'group', 'label', and 'fr';

        
        Parameters:
        group_names: list of group names, where the spikes will be
        searched. When this input is 'None' all groups are taken. The resulting
        dictionary has the first keys as groups, and the second keys as the
        respective cluster id's, whereas the value, is the corresponding firing
        rate within [a,b].

        Please refer to the method cluster_firing_rate in order to get more 
        details about the firing calculation.

        Author: Nivaldo A P de Vasconcelos
        Date: 2018.Feb.02 
        """
        d=self.group_firing_rate (group_names=group_names,a=a,b=b)

        shank_id = self.name
        group_names = d.keys()
        data=[]
        for group_name in group_names:
            for label in d[group_name].keys():        
                fr=d[group_name][label]
                data.append ({"shank_id":shank_id, "group":group_name,"label":label,"fr":fr})
        return (pd.DataFrame(data))  
    def crop(self, cluster_id, a, b):
        spikes_j=np.array(self.spikes_on_cluster(cluster_id=cluster_id))/self.sampling_rate()
        spikes_j=spikes_j[(spikes_j>=a) & (spikes_j<=b)]
        return spikes_j  
    def corr(self, lista, a, b, sigma=0.05, fr_th = 0.01, bit_size = 1e-03):
        N = len(lista)

        correlacoes = []
        t = np.arange(-4*sigma, 4*sigma, bit_size)
        gaussian_kernel = ((1.0/np.sqrt(2*np.pi))/sigma)*np.exp(-(t**2.)/(2.*sigma**2.0))

        for i in range(0,N):
            spikes_i=self.crop(lista[i],a,b)
            fr_i = len(spikes_i)/(b-a)

            if (fr_i > fr_th):
                binnep_pop_i,bins=np.histogram(spikes_i,bins=np.arange(a,b,bit_size))
                aux_i = np.convolve(binnep_pop_i, gaussian_kernel, 'same')
                for j in range(i+1,N):
                    spikes_j=self.crop(lista[j],a,b)
                    fr_j = len(spikes_j)/(b-a)
                    if (fr_j > fr_th):
                        binnep_pop_j,bins=np.histogram(spikes_j,bins=np.arange(a,b,bit_size))
                        aux_j = np.convolve(binnep_pop_j, gaussian_kernel, 'same')
                        correlacoes.append(np.corrcoef(aux_i,aux_j)[1,0])
        return(np.array(correlacoes)) 
    def corr_mexican(self, lista, a, b, sigma=0.05, fr_th = 0.01, bit_size = 1e-03):
        N = len(lista)

        correlacoes = []
        t = np.arange(-4*sigma, 4*sigma, bit_size)
        mexican_kernel = ((2.0/np.sqrt(3*sigma))/np.pi**0.25)*(1 - (t/sigma)**2)*np.exp(-(t**2.)/(2.*sigma**2.0))

        for i in range(0,N):
            spikes_i=self.crop(lista[i],a,b)
            fr_i = len(spikes_i)/(b-a)

            if (fr_i > fr_th):
                binnep_pop_i,bins=np.histogram(spikes_i,bins=np.arange(a,b,bit_size))
                aux_i = np.convolve(binnep_pop_i, mexican_kernel, 'same')
                for j in range(i+1,N):
                    spikes_j=self.crop(lista[j],a,b)
                    fr_j = len(spikes_j)/(b-a)
                    if (fr_j > fr_th):
                        binnep_pop_j,bins=np.histogram(spikes_j,bins=np.arange(a,b,bit_size))
                        aux_j = np.convolve(binnep_pop_j, mexican_kernel, 'same')
                        correlacoes.append(np.corrcoef(aux_i,aux_j)[1,0])
        return(np.array(correlacoes))
        





#title Função para extrair informações do kwik para um csv


def export_kwikTOdf(datafolder,date,num_shanks,only_SUA,chan_info=True,save=True):
    
    spikes_per_shank=[]
    clu_id=[]
    shank=[]
    chan_max=[]
    for i in range(num_shanks):
        #carregar o dado
        kpath="{}/{}/g0{}/g0{}.kwik".format(datafolder,date,i+1,i+1) 
        print(kpath)
      
        #chamar a classe
        km=KwikFile(kpath)
    
        
        list_of_non_noisy_groups=km.list_of_non_noisy_groups()
        
        if len(list_of_non_noisy_groups)==1: # Exemplo: ['good', 'mua', 'unsorted']
            print("The Manual Spike Sorting Has not been made")
            break
        if only_SUA:
            type_clu="SUA"
            if "usorted" in list_of_non_noisy_groups: 
                id_clu=km.clusters(group_name="unsorted")+km.clusters(group_name="good")
            else:
                id_clu=km.clusters(group_name="good")
        else:
            type_clu="ALL"
            if "usorted" in list_of_non_noisy_groups: 
                id_clu=km.clusters(group_name="unsorted")+km.clusters(group_name="good")+km.clusters(group_name="mua")
                
            else:
                id_clu=km.clusters(group_name="good")+km.clusters(group_name="mua")
                
        dict_waveMean={}
      
        for k in id_clu:
            print(k)
            w=km.get_cluster_waveforms(k) #(n_spikes, n_samples, n_channels)`
            print(w)
            w_mean=w.mean(axis=0)
            dict_waveMean[k]=w_mean


            w=[] # uma forma de apagar o w porque ele é muito grande
            amplitude_channels=w_mean[0,:]-w_mean.min(axis=0) # Pode ser alterado a forma de calcular a amplitude aqui
            channel_max=amplitude_channels.argmax()


            spikes_per_shank.extend(km.spikes_on_cluster(k))
            shank.extend(i*np.ones(len(km.spikes_on_cluster(k))))
            clu_id.extend(k*np.ones(len(km.spikes_on_cluster(k))))
            chan_max.extend(channel_max*np.ones(len(km.spikes_on_cluster(k))))

        f = open("{}/{}/waveforms_mean_{}_{}_g{}.txt".format(datafolder,date,type_clu,date,i),'w')
        f.write(str(dict_waveMean))
        f.close()
  
    if len(list_of_non_noisy_groups)==1: # Exemplo: ['good', 'mua', 'unsorted']
        print("The Manual Spike Sorting Has not been made")
    else:
        # Sortear os spikes    
        sampling_rate=30000
        spikes_per_shank = [i/sampling_rate for i in spikes_per_shank] 
        if chan_info:
            spikes_per_shank , clu_id, shank,chan_max = zip(*sorted(zip(spikes_per_shank , clu_id, shank,chan_max)))
            df=pd.DataFrame({"time(s)":spikes_per_shank,"cluID":clu_id,"shankID":shank,"chanMAX":chan_max})
        else:
            spikes_per_shank , clu_id, shank = zip(*sorted(zip(spikes_per_shank , clu_id, shank)))
            df=pd.DataFrame({"time(s)":spikes_per_shank,"cluID":clu_id,"shankID":shank})      
        if save:
            filename="{}/{}/time_cluster_shank_{}_{}.csv".format(datafolder,date,type_clu,date)
            df.to_csv(filename)

  
    return df       
        

def info_tsv(datafolder,date,num_shanks,only_SUA,chan_info=True,save=True):  

    for i in range(num_shanks):
        #carregar o dado
        kpath="{}/{}/g0{}.kwik".format(datafolder,date,i+1,i+1) 
        print(kpath)
        #chamar a classe
        km=KwikFile(kpath)
        xinfo=km.data_info_clusters(group_name="good")
        spikes_clu_id=list(xinfo.keys())
        spikes_grupo=list(xinfo.values())

        mi_path="{}/info_{}".format(datafolder,date) 
        #z_path="{}/info_{}/g0{}_cluster_info.csv".format(datafolder,date,i+1)
        
        if(os.path.isdir(mi_path)==False):
            path = Path(mi_path)
            path.mkdir(parents=True)#crear capetas

        df=pd.DataFrame({"id":spikes_clu_id,"group":spikes_grupo}) 
        filename="{}/info_{}/g0{}_cluster_info.csv".format(datafolder,date,i+1)
        df.to_csv(filename)

    return xinfo

datafolder="Files_kwik"
only_SUA=False # Se for True, os clusters usados serão só SUA e se for False, todos os clusters serão usados exceto os noise
chan_info=True # Se for True, o canal de máxima amplitude será informado, False essa informação não será salva
#date_opt=["2019Dez20","2020Jan29"]
date_opt=["Jan29"]#,"2020Mar10","2020Jan14","2020Jan21"]#,"2019Dez20","2020Jan29"]

num_shanks=6
for date in date_opt:
    kpath="{}/{}".format(datafolder,date)  
     
    ykpath=info_tsv(datafolder,date,num_shanks,only_SUA)
    
    file_txt="id_forma_ondas_date_opt.txt"
    #file_txt="datMetx/"+xsubdir+"_fracion.txt"
    fdata = open(file_txt, 'w')
    
    # try:
    #     for i in range(len(dT)):
    #             fdata.write(str(dT[i])+"\t"+str(average_CV[i]))    
    # finally:
    #     fdata.close()


    