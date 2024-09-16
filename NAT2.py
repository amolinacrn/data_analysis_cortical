# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 10:52:38 2022

@author: CSN Admin
"""
from turtle import pos
import scipy.io
import os
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np

import scipy.io
import os
from os import listdir
from os.path import isfile, join
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import numpy.matlib
import scipy.optimize as opt
import random
import pandas as pd
import matplotlib.gridspec as gridspec
#import seaborn as sns
from scipy import stats
import glob
import seaborn as sns
from glob import glob
#import Packages
#@title NAT2 

# plt.rcParams["figure.figsize"] = [10,6]
# plt.rcParams.update({'font.size': 18})
#title NAT2 

# plt.rcParams["figure.figsize"] = [10,6]
# plt.rcParams.update({'font.size': 18})
class NAT2:
    def __init__(self, data_original):
        self.data = data_original #leitura do arquivo
        #print('File loaded successfully!')
        
        
    def CV_calc(self,cv_time,bin_cv, show=False, save=False, ADFA=False, showbin=False, bin_size=None): 
        #Retornar CV
        self.cv_time = cv_time 
        #janela para o cálculo do CV, normalmente 10s 
        self.bin_cv = bin_cv 
        #binarização para gerar a serie temporal para o cálculo de CV,
        #em geral 0.05 = 50 ms
        
        #Configurações das figuras
        #plt.rcParams["figure.figsize"] = [10,6]
        #plt.rcParams.update({'font.size': 18})
        
        
        self.n_cv = int(np.round((self.data[-1]-self.data[0])/cv_time))
        # O número de cv's é equivalente ao tempo do último spike menos o primeiro spike dividido pelo bin de 10 segundos
        #print("Total CV's = ", self.n_cv)
        
        self.CV = []#np.zeros(self.n_cv)
        self.SIZES = []
        self.DURATIONS = []
        self.bin_sizes = []
        matriz_spikes=[]
        if ADFA:
            self.ADFA = True
            self.DFA = np.zeros(self.n_cv)
        else:
            self.ADFA = False
     
        for t in range(0,self.n_cv):
           # print('janela: ', t)
            a = t*cv_time + self.data[0]
            b = a + cv_time
            spk = self.data[(self.data>=a) & (self.data<=b)]
            
            
            if (len(spk)==0 or len(spk)==1 or len(spk)==2):
                print('Janela vacia: ', t)
                self.CV.append(100)
                matriz_spikes.append(-100)
               # print(spk)
                continue

            matriz_spikes.append(spk.tolist())
            
            nbins = np.arange(a,b,bin_cv)#numpy.arange([start, ]stop, [step, ]dtype=None)
            count, edges=np.histogram(spk, bins = nbins)

            self.CV.append(np.std(count)/np.mean(count))
            
            if ADFA:
                aux = self.avalanches(spk, ADFA=True, show=showbin, bin_size=bin_size)
            else:
                aux = self.avalanches(spk,show=showbin, bin_size=bin_size)
            
           # print(len(aux[0]), len#,label=r"Experimental data "+str(xsubC)+", $\Delta{t}=$"+str(delt[i])+" s")
            
            self.SIZES.append(aux[0])
            self.DURATIONS.append(aux[1])
            self.bin_sizes.append(aux[2])
            if ADFA:
                self.DFA[t]=aux[3]
        
        #print('<CV>: {:.2f} ± {:.2f}'.format(np.mean(self.CV), np.std(self.CV)))
        #print('<ISI>: {:.2f} ± {:.2f} (ms)'.format(np.mean(self.bin_sizes)*1000, np.std(self.bin_sizes)*1000))
        self.CV = np.array(self.CV)
        self.I = np.argsort(self.CV)
        self.n_cv = len(self.CV)
        if show:
            df = {'x':np.array(range(len(self.CV)))*cv_time, 'y': self.CV}
            
            sns.lineplot(x = "x", y = "y", data = df,marker = "o",color="blue")
            # plt.plot(np.array(range(len(self.CV)))*cv_time, self.CV,"o")
            plt.ylabel('CV')
            
            plt.xlabel('t (s)')
            
            plt.tight_layout()
            if save:
                plt.savefig('CV_t.png', dpi = 150)
            plt.show()
            n, bins, patches = plt.hist(self.CV, bins = 50, density=1, facecolor='green', alpha=1)
            plt.title('CV histogram')
            plt.xlabel("CV")
            plt.ylabel("Density")
            plt.tight_layout()
            if save:
                plt.savefig('hist_CV_t.png', dpi = 150)
            plt.show()    
            
        return(self.CV,matriz_spikes)
    
    
    def raster(self, nCV, graph = True):
        t = self.I[nCV]
        a = t*self.cv_time
        b = a + self.cv_time
        spk = self.data[(self.data>=a) & (self.data<=b)]
        nbins = np.arange(a,b,self.bin_cv)
        count, edges=np.histogram(spk, bins = nbins)
        #print('CV = {:.2f}'.format(np.std(count)/np.mean(count)))
        samples = self.data[(self.data>=a) & (self.data<=b)]
        mfr=len(samples)/(samples[-1]-samples[0])
        bin_size=1.0/mfr
        self.count=count
        #print('<ISI> = {:.2f} ms'.format(bin_size*1000))
        if graph:
            #plt.title('CV = {}'.format(nCV))
            plt.ylabel('#spikes')
            plt.xlabel('Time (s)')
            plt.plot(np.arange(a,b,self.bin_cv)[:-1], self.count)
            plt.tight_layout()
            plt.savefig('raster_{}.png'.format(nCV), dpi = 150)
            plt.show()
        return self.count
        
    def concatenation(self, nCV1, nCV2, xmin=None, tmin=None, xmax=None, tmax=None, show=False):
        PS = np.zeros(0)
        PD = np.zeros(0)
        
        for i in range(nCV1, nCV2):
            PS = np.concatenate((PS, self.SIZES[self.I[i]]), axis=None)
            PD = np.concatenate((PD, self.DURATIONS[self.I[i]]), axis=None)
           
        self.PD=PD#testes
        self.PS=PS
        if (xmin==None):
            xmin = np.min(PS)
        if (xmax==None):
            xmax = np.max(PS)

        if (tmin==None):
            tmin = np.min(PD)

        if (tmax==None):
            tmax = np.max(PD)

        aux_s = self.mle(PS,xmin,xmax)
        aux_d = self.mle(PD,tmin,tmax)
        tau_s = aux_s[0]
        tau_d = aux_d[0]

        #bootstrap
        
        exp_erro_s = np.zeros(25)
        exp_erro_d = np.zeros(25)
        
        for i in range(25):
            test_PS = np.random.choice(PS, size=int(len(PS)))
            test_PD = np.random.choice(PD, size=int(len(PD)))
            aux_s = self.mle(test_PS,xmin,xmax)
            aux_d = self.mle(test_PD,tmin,tmax)
            exp_erro_s[i] = aux_s[0]
            exp_erro_d[i] = aux_d[0]
        
        erro_s = np.std(exp_erro_s)
        erro_d = np.std(exp_erro_d)
        
        AIC_s = (2*2 - 2*self.lnormal1(PS, xmin, xmax)) - (2*1 - 2*aux_s[1])
        AIC_d = (2*2 - 2*self.lnormal1(PD, tmin, tmax)) - (2*1 - 2*aux_d[1])
        beta_calc = (tau_d - 1.0)/(tau_s - 1.0)
        
        erro_relation = beta_calc*np.sqrt((erro_s/tau_s)**2 + (erro_d/tau_d)**2)

        #print(i,len(PS),len(PD))
        aux = self.size_avg(PD, PS)
        beta_data = self.fitt(aux[0],aux[1],tmin,tmax) 
        
        if show:############################################################
            plt.rcParams["figure.figsize"] = [14,10]
            plt.rcParams.update({'font.size': 14})
            f, axarr = plt.subplots(2, 3)
            f.suptitle("xmin: {}, xmax = {}, tmin = {}, tmax = {}, #avalanches = {}".format(xmin,xmax,tmin,tmax,len(PS)))
        # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
            if (max(PS)>10000):
                max_PS = 10000
            else:
                max_PS = int(max(PS))
            
            n, bins = np.histogram(PS, max_PS, density=True)

            def coeff_s(x, a):
                return (a*x**(-tau_s))

            fit = opt.curve_fit(coeff_s, bins[:-1][(bins[:-1]>=xmin) & (bins[:-1]<=xmax)], n[(bins[:-1]>=xmin) & (bins[:-1]<=xmax)])
            a_fit = fit[0]

            
            #n, bins = np.histogram(PS, max_PS, normed=1)
            x = np.arange(xmin,xmax+1,1)
            y = a_fit[0]*x**(-tau_s)
            
            axarr[0, 0].plot(bins[:-1], n, 'o')
            axarr[0, 0].plot(x,y,'--', linewidth=3, label = r"$\tau$ = {:.2f} ± {:.2f}".format(tau_s, erro_s))
            #plt.plot(x,y2,'--')
            axarr[0, 0].set_xlabel('Size')
            axarr[0, 0].set_ylabel('Probability')
            axarr[0, 0].set_xscale('log')
            axarr[0, 0].set_yscale('log')
            axarr[0, 0].legend()
            #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
            #plt.axis([40, 160, 0, 0.03])
            axarr[0, 0].grid(True)
            
            log_min_size = np.log10(min(PS))
            log_max_size = np.log10(max(PS))
            number_of_bins = 35 #35 pela estetica
            #number_of_bins = np.ceil((log_max_size-log_min_size)*10) #comando Antonio
            bins_ton=np.unique(np.floor(np.logspace(log_min_size, log_max_size, num=number_of_bins)))
            prob, edges = np.histogram(PS, bins_ton, density=True)
            bin_centers = 0.5 * (bins_ton[:-1] + bins_ton[1:])

            #x = np.arange(xmin,xmax,1)
            
            fit = opt.curve_fit(coeff_s, bin_centers[(bin_centers>=xmin) & (bin_centers<=xmax)], prob[(bin_centers>=xmin) & (bin_centers<=xmax)])
            a_fit = fit[0]

            y = a_fit[0]*x**(-tau_s)            
            
            axarr[1, 0].plot(bin_centers, prob, 'o')
            axarr[1, 0].plot(x,y,'--', linewidth = 3, label = r"$\tau$ = {:.2f} ± {:.2f}".format(tau_s, erro_s))
            # 100 = bins, alpha, tonalidade da cor de saida
            #plt.plot(x,y2,'--')
            axarr[1, 0].set_xlabel('Size')
            axarr[1, 0].set_ylabel('Probability')
            axarr[1, 0].set_xscale('log')
            axarr[1, 0].set_yscale('log')
            axarr[1, 0].grid(True)

            #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
            #plt.axis([40, 160, 0, 0.03])
            axarr[1, 0].legend()
            #plt.savefig('ps_350.png')
            plt.grid(True)
            
            n, bins = np.histogram(PD, int(max(PD)), density=True)
            
            x = np.arange(tmin,tmax+1,1)
            
            def coeff_d(x, a):
                return (a*x**(-tau_d))

            fit = opt.curve_fit(coeff_d, bins[:-1][(bins[:-1]>=tmin) & (bins[:-1]<=tmax)], n[(bins[:-1]>=tmin) & (bins[:-1]<=tmax)])
            a_fit = fit[0]

            y = a_fit[0]*x**(-tau_d)            
            
            
            axarr[0, 1].plot(bins[:-1], n, 'o')
            axarr[0, 1].plot(x,y,'--', linewidth=3, label = r"$\tau_t$ = {:.2f} ± {:.2f}".format(tau_d, erro_d))
            #plt.plot(x,y2,'--')
            axarr[0, 1].set_xlabel('Duration')
            axarr[0, 1].set_ylabel('Probability')
            axarr[0, 1].set_xscale('log')
            axarr[0, 1].set_yscale('log')
            axarr[0, 1].legend()
            #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
            #plt.axis([40, 160, 0, 0.03])
            axarr[0, 1].grid(True)
            
            log_min_size = np.log10(min(PD))
            log_max_size = np.log10(max(PD))
            number_of_bins = 35 #35 pela estetica
            #number_of_bins = np.ceil((log_max_size-log_min_size)*10) #comando Antonio
            bins_ton=np.unique(np.floor(np.logspace(log_min_size, log_max_size, num=number_of_bins)))
            prob, edges = np.histogram(PD, bins_ton, density=True)
            bin_centers = 0.5 * (bins_ton[:-1] + bins_ton[1:])

            fit = opt.curve_fit(coeff_d, bin_centers[(bin_centers>=tmin) & (bin_centers<=tmax)], prob[(bin_centers>=tmin) & (bin_centers<=tmax)])
            a_fit = fit[0]

            y = a_fit[0]*x**(-tau_d)
            
            
            axarr[1, 1].plot(bin_centers, prob, 'o')
            axarr[1, 1].plot(x,y,'--', linewidth = 3, label = r"$\tau_t$ = {:.2f} ± {:.2f}".format(tau_d, erro_d))
            # 100 = bins, alpha, tonalidade da cor de saida

            #plt.plot(x,y2,'--')
            axarr[1, 1].set_xlabel('Duration')
            axarr[1, 1].set_ylabel('Probability')
            axarr[1, 1].set_xscale('log')
            axarr[1, 1].set_yscale('log')

            #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
            #plt.axis([40, 160, 0, 0.03])
            axarr[1, 1].legend()
            #plt.savefig('ps_350.png')
            axarr[1, 1].grid(True)
            
            def coeff_beta(x, a):
                return (a*x**(beta_data[0]))

            fit = opt.curve_fit(coeff_beta, aux[0][(aux[0]>=tmin) & (aux[0]<=tmax)], aux[1][(aux[0]>=tmin) & (aux[0]<=tmax)])
            a_fit = fit[0]

            y = a_fit[0]*x**(beta_data[0])
            
            #y = 1.3*x**(beta_data[0])
            #g = 2.00
            #y2 = 6*x**(g)
            #plt.scatter(a, b)
            axarr[0, 2].plot(aux[0],aux[1], 'o', markersize=7)
            axarr[0, 2].plot(x,y,'--', label = r'$\frac{1}{\sigma \nu z}$ = ' + '{0:.2f}'.format(beta_data[0]), linewidth=3)
            #plt.plot(x,y2,'g--', label = r'$\frac{1}{\sigma \nu z}$ = ' + '{0:.2f}'.format(g), linewidth=3)

            axarr[0, 2].set_xlabel('Duration')
            axarr[0, 2].set_ylabel('Size')
            axarr[0, 2].set_xscale('log')
            axarr[0, 2].set_yscale('log')
            axarr[0, 2].legend()
            axarr[0, 2].grid(True)
            
            
            #axarr[1,2].plot([0],[beta_calc], 'o', markersize=7, label=r'$\frac{\tau_t - 1}{\tau -1}$')
            axarr[1,2].errorbar([0],[beta_calc],yerr=[erro_relation], markersize=7, fmt="o", label=r'$\frac{\tau_t - 1}{\tau -1}$')
            axarr[1,2].errorbar([0],[beta_data[0]],yerr=[beta_data[1]], markersize=7, fmt="o", label=r'$\frac{1}{\sigma\nu z}$')
            axarr[1,2].set_ylabel('Scaling relation')
            axarr[1,2].set_xlabel('CV range: {0:.2f}/{1:.2f}'.format(self.CV[self.I[nCV1]],self.CV[self.I[nCV2]]))
            axarr[1, 2].legend()
            axarr[1, 2].grid(True)
            
            plt.tight_layout()
            f.tight_layout()
            f.subplots_adjust(top=0.94)
            plt.savefig('concatenation.png', dpi=150)
            plt.show()
            plt.rcParams["figure.figsize"] = [10,6]
            plt.rcParams.update({'font.size': 18})
            return()
        return(tau_s, tau_d, beta_calc, beta_data, AIC_s, AIC_d, PS, PD, erro_s, erro_d, erro_relation)
    
    def concatenation2(self, nCV1, nCV2, xmin=None, tmin=None, xmax=None, tmax=None, show=False):
        '''Vários dados num plot só'''
        PS = np.zeros(0)
        PD = np.zeros(0)
        
        for i in range(nCV1, nCV2):
            PS = np.concatenate((PS, self.SIZES[self.I[i]]), axis=None)
            PD = np.concatenate((PD, self.DURATIONS[self.I[i]]), axis=None)
            
        #self.PD=4*PD#testes

        if (xmin==None):
            xmin = np.min(PS)
        if (xmax==None):
            xmax = np.max(PS)

        if (tmin==None):
            tmin = np.min(PD)
        if (tmax==None):
            tmax = np.max(PD)

        aux_s = self.mle(PS,xmin,xmax)
        aux_d = self.mle(PD,tmin,tmax)
        tau_s = aux_s[0]
        tau_d = aux_d[0]

        #bootstrap
        
        exp_erro_s = np.zeros(25)
        exp_erro_d = np.zeros(25)
        
        for i in range(25):
            test_PS = np.random.choice(PS, size=int(len(PS)))
            test_PD = np.random.choice(PD, size=int(len(PD)))
            aux_s = self.mle(test_PS,xmin,xmax)
            aux_d = self.mle(test_PD,tmin,tmax)
            exp_erro_s[i] = aux_s[0]
            exp_erro_d[i] = aux_d[0]
        
        erro_s = np.std(exp_erro_s)
        erro_d = np.std(exp_erro_d)
        
        AIC_s = (2*2 - 2*self.lnormal1(PS, xmin, xmax)) - (2*1 - 2*aux_s[1])
        AIC_d = (2*2 - 2*self.lnormal1(PD, tmin, tmax)) - (2*1 - 2*aux_d[1])
        beta_calc = (tau_d - 1.0)/(tau_s - 1.0)
        
        erro_relation = beta_calc*np.sqrt((erro_s/tau_s)**2 + (erro_d/tau_d)**2)

        #print(i,len(PS),len(PD))
        aux = self.size_avg(PD, PS)
        beta_data = self.fitt(aux[0],aux[1],tmin,tmax) 
        
        if show:############################################################
            plt.rcParams["figure.figsize"] = [14,10]
            plt.rcParams.update({'font.size': 14})
            f, axarr = plt.subplots(2, 3)
            f.suptitle("xmin: {}, xmax = {}, tmin = {}, tmax = {}, #avalanches = {}".format(xmin,xmax,tmin,tmax,len(PS)))
        # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
            if (max(PS)>10000):
                max_PS = 10000
            else:
                max_PS = int(max(PS))
            
            n, bins = np.histogram(PS, max_PS, density=True)

            def coeff_s(x, a):
                return (a*x**(-tau_s))

            fit = opt.curve_fit(coeff_s, bins[:-1][(bins[:-1]>=xmin) & (bins[:-1]<=xmax)], n[(bins[:-1]>=xmin) & (bins[:-1]<=xmax)])
            a_fit = fit[0]

            
            #n, bins = np.histogram(PS, max_PS, normed=1)
            x = np.arange(xmin,xmax+1,1)
            y = a_fit[0]*x**(-tau_s)
            
            axarr[0, 0].plot(bins[:-1], n, 'o')
            axarr[0, 0].plot(x,y,'--', linewidth=3, label = r"$\tau$ = {:.2f} ± {:.2f}".format(tau_s, erro_s))
            #plt.plot(x,y2,'--')
            axarr[0, 0].set_xlabel('Size')
            axarr[0, 0].set_ylabel('Probability')
            axarr[0, 0].set_xscale('log')
            axarr[0, 0].set_yscale('log')
            axarr[0, 0].legend()
            #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
            #plt.axis([40, 160, 0, 0.03])
            axarr[0, 0].grid(True)
            
            log_min_size = np.log10(min(PS))
            log_max_size = np.log10(max(PS))
            number_of_bins = 35 #35 pela estetica
            #number_of_bins = np.ceil((log_max_size-log_min_size)*10) #comando Antonio
            bins_ton=np.unique(np.floor(np.logspace(log_min_size, log_max_size, num=number_of_bins)))
            prob, edges = np.histogram(PS, bins_ton, density=True)
            bin_centers = 0.5 * (bins_ton[:-1] + bins_ton[1:])

            #x = np.arange(xmin,xmax,1)
            
            fit = opt.curve_fit(coeff_s, bin_centers[(bin_centers>=xmin) & (bin_centers<=xmax)], prob[(bin_centers>=xmin) & (bin_centers<=xmax)])
            a_fit = fit[0]

            y = a_fit[0]*x**(-tau_s)            
            
            axarr[1, 0].plot(bin_centers, prob, 'o')
            axarr[1, 0].plot(x,y,'--', linewidth = 3, label = r"$\tau$ = {:.2f} ± {:.2f}".format(tau_s, erro_s))
            # 100 = bins, alpha, tonalidade da cor de saida
            #plt.plot(x,y2,'--')
            axarr[1, 0].set_xlabel('Size')
            axarr[1, 0].set_ylabel('Probability')
            axarr[1, 0].set_xscale('log')
            axarr[1, 0].set_yscale('log')
            axarr[1, 0].grid(True)

            #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
            #plt.axis([40, 160, 0, 0.03])
            axarr[1, 0].legend()
            #plt.savefig('ps_350.png')
            plt.grid(True)
            
            n, bins = np.histogram(PD, int(max(PD)), density=True)
            
            x = np.arange(tmin,tmax+1,1)
            
            def coeff_d(x, a):
                return (a*x**(-tau_d))

            fit = opt.curve_fit(coeff_d, bins[:-1][(bins[:-1]>=tmin) & (bins[:-1]<=tmax)], n[(bins[:-1]>=tmin) & (bins[:-1]<=tmax)])
            a_fit = fit[0]

            y = a_fit[0]*x**(-tau_d)            
            
            
            axarr[0, 1].plot(bins[:-1], n, 'o')
            axarr[0, 1].plot(x,y,'--', linewidth=3, label = r"$\tau_t$ = {:.2f} ± {:.2f}".format(tau_d, erro_d))
            #plt.plot(x,y2,'--')
            axarr[0, 1].set_xlabel('Duration')
            axarr[0, 1].set_ylabel('Probability')
            axarr[0, 1].set_xscale('log')
            axarr[0, 1].set_yscale('log')
            axarr[0, 1].legend()
            #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
            #plt.axis([40, 160, 0, 0.03])
            axarr[0, 1].grid(True)
            
            log_min_size = np.log10(min(PD))
            log_max_size = np.log10(max(PD))
            number_of_bins = 35 #35 pela estetica
            #number_of_bins = np.ceil((log_max_size-log_min_size)*10) #comando Antonio
            bins_ton=np.unique(np.floor(np.logspace(log_min_size, log_max_size, num=number_of_bins)))
            prob, edges = np.histogram(PD, bins_ton, density=True)
            bin_centers = 0.5 * (bins_ton[:-1] + bins_ton[1:])

            fit = opt.curve_fit(coeff_d, bin_centers[(bin_centers>=tmin) & (bin_centers<=tmax)], prob[(bin_centers>=tmin) & (bin_centers<=tmax)])
            a_fit = fit[0]

            y = a_fit[0]*x**(-tau_d)
            
            
            axarr[1, 1].plot(bin_centers, prob, 'o')
            axarr[1, 1].plot(x,y,'--', linewidth = 3, label = r"$\tau_t$ = {:.2f} ± {:.2f}".format(tau_d, erro_d))
            # 100 = bins, alpha, tonalidade da cor de saida

            #plt.plot(x,y2,'--')
            axarr[1, 1].set_xlabel('Duration')
            axarr[1, 1].set_ylabel('Probability')
            axarr[1, 1].set_xscale('log')
            axarr[1, 1].set_yscale('log')

            #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
            #plt.axis([40, 160, 0, 0.03])
            axarr[1, 1].legend()
            #plt.savefig('ps_350.png')
            axarr[1, 1].grid(True)
            
            def coeff_beta(x, a):
                return (a*x**(beta_data[0]))

            fit = opt.curve_fit(coeff_beta, aux[0][(aux[0]>=tmin) & (aux[0]<=tmax)], aux[1][(aux[0]>=tmin) & (aux[0]<=tmax)])
            a_fit = fit[0]

            y = a_fit[0]*x**(beta_data[0])
            
            #y = 1.3*x**(beta_data[0])
            #g = 2.00
            #y2 = 6*x**(g)
            #plt.scatter(a, b)
            axarr[0, 2].plot(aux[0],aux[1], 'o', markersize=7)
            axarr[0, 2].plot(x,y,'--', label = r'$\frac{1}{\sigma \nu z}$ = ' + '{0:.2f}'.format(beta_data[0]), linewidth=3)
            #plt.plot(x,y2,'g--', label = r'$\frac{1}{\sigma \nu z}$ = ' + '{0:.2f}'.format(g), linewidth=3)

            axarr[0, 2].set_xlabel('Duration')
            axarr[0, 2].set_ylabel('Size')
            axarr[0, 2].set_xscale('log')
            axarr[0, 2].set_yscale('log')
            axarr[0, 2].legend()
            axarr[0, 2].grid(True)
            
            
            #axarr[1,2].plot([0],[beta_calc], 'o', markersize=7, label=r'$\frac{\tau_t - 1}{\tau -1}$')
            axarr[1,2].errorbar([0],[beta_calc],yerr=[erro_relation], markersize=7, fmt="o", label=r'$\frac{\tau_t - 1}{\tau -1}$')
            axarr[1,2].errorbar([0],[beta_data[0]],yerr=[beta_data[1]], markersize=7, fmt="o", label=r'$\frac{1}{\sigma\nu z}$')
            axarr[1,2].set_ylabel('Scaling relation')
            axarr[1,2].set_xlabel('CV range: {0:.2f}/{1:.2f}'.format(self.CV[self.I[nCV1]],self.CV[self.I[nCV2]]))
            axarr[1, 2].legend()
            axarr[1, 2].grid(True)
            
            plt.tight_layout()
            f.tight_layout()
            f.subplots_adjust(top=0.94)
            plt.savefig('concatenation.png', dpi=150)
            plt.show()
            plt.rcParams["figure.figsize"] = [10,6]
            plt.rcParams.update({'font.size': 18})
            return()
        return(tau_s, tau_d, beta_calc, beta_data, AIC_s, AIC_d, PS, PD, erro_s, erro_d, erro_relation)
    
    
    def concatenation3(self, nCV1, nCV2, xmin=2, tmin=2, xmax=100, tmax=30, show=False):
        PS = np.zeros(0)
        PD = np.zeros(0)
        aux = 0
        for j in range(nCV1,nCV2):
                aux+=self.CV[self.I[j]]
        
        CV_medio = aux/(nCV2 - nCV1)
        
        for i in range(nCV1, nCV2):
            PS = np.concatenate((PS, self.SIZES[self.I[i]]), axis=None)
            PD = np.concatenate((PD, self.DURATIONS[self.I[i]]), axis=None)
        self.PD=PD#testes
        self.PS=PS 
        
        #aux_s = self.mle(PS,xmin,xmax)
        #aux_d = self.mle(PD,tmin,tmax)
        #tau_s = aux_s[0]
        #tau_d = aux_d[0]
        
        #AIC_s = (2*2 - 2*self.lnormal1(PS, xmin, xmax)) - (2*1 - 2*aux_s[1])
       # AIC_d = (2*2 - 2*self.lnormal1(PD, tmin, tmax)) - (2*1 - 2*aux_d[1])
        #beta_calc = (tau_d - 1.0)/(tau_s - 1.0)
        #print(i,len(PS),len(PD))
       # aux = self.size_avg(PD, PS)
       # beta_data = self.fitt(aux[0],aux[1],tmin,tmax) 
       
        return(CV_medio, PS, PD)
    

    def concatenation4(self, nCV1, nCV2):#, xmin=2, tmin=2, xmax=100, tmax=30, show=False):
      """Para fazer a figura fora da classe"""
      PS = np.zeros(0)
      PD = np.zeros(0)
      
      for i in range(nCV1, nCV2):
          PS = np.concatenate((PS, self.SIZES[self.I[i]]), axis=None)
          PD = np.concatenate((PD, self.DURATIONS[self.I[i]]), axis=None)
          
  

      self.PS_plot=PS#,bins]#,x,y]
    
      self.PD_plot=PD#n,bins]#,x,y]
      
      return 

  
    def scaling(self, juntos_cv, xmim=None, tmim=None, xmax0=100, tmax0=30, save=True): #juntos_cv = quantos CV's serão agrupados
        
        if (xmim==None):
            xmim = np.min(self.PS)
        if (xmax0==None):
            xmax0 = np.max(self.PS)

        if (tmim==None):
            tmim = np.min(self.PD)

        if (tmax0==None):
            tmax0 = np.max(self.PD)
        
        
        
        total = int(np.floor(self.n_cv/juntos_cv))#quantidade de pontos
        self.c_calc = np.zeros(total)
        self.c_data = np.zeros(total)
        self.c_data_erro = np.zeros(total)
        self.c_CV = np.zeros(total)
        self.c_CV_erro = np.zeros(total)
        self.c_AICs = np.zeros(total)
        self.c_AICd = np.zeros(total)
        self.c_tau = np.zeros(total)
        self.c_tau_t = np.zeros(total)
        self.c_erro_s = np.zeros(total)
        self.c_erro_d = np.zeros(total)
        self.c_erro_relation = np.zeros(total)
        
        if self.ADFA:
            self.c_DFA = np.zeros(total)
            for i in range(total):
                aux = 0
                for j in range(i*juntos_cv,(i+1)*juntos_cv):
                    aux+=self.DFA[self.I[j]]
                self.c_DFA[i] = aux/juntos_cv
            
        for i in range(total):
            c = self.concatenation(i*juntos_cv, (i+1)*juntos_cv, xmin=xmim, tmin=tmim, xmax=xmax0, tmax=tmax0)
            aux = 0
            for j in range(i*juntos_cv,(i+1)*juntos_cv):
                aux+=self.CV[self.I[j]]
                
            self.c_tau[i] = c[0]
            self.c_tau_t[i] = c[1]
            self.c_calc[i] = c[2]
            self.c_data[i] = c[3][0]
            self.c_data_erro[i] = c[3][1]
            self.c_AICs[i] = c[4]
            self.c_AICd[i] = c[5]
            self.c_CV[i] = aux/juntos_cv
            self.c_CV_erro[i] = np.std(self.CV[self.I[i*juntos_cv:(i+1)*juntos_cv]])
            self.c_erro_s[i] = c[8]
            self.c_erro_d[i] = c[9]
            self.c_erro_relation[i] = c[10] 
        
        # c=tau_s, tau_d, beta_calc, beta_data, AIC_s, AIC_d, PS, PD, erro_s, erro_d, erro_relation)

        plt.xlabel('<CV>')
        plt.ylabel('Scaling relation')
        #plt.plot(self.c_CV, self.c_data, 'o',markersize=8, label=r'$\frac{1}{\sigma\nu z}$')
        #plt.plot(self.c_CV, self.c_calc, 'o', markersize=8, label=r'$\frac{\tau_t - 1}{\tau -1}$')
        
        #plt.errorbar(self.c_CV, self.c_calc, xerr= self.c_CV_erro, fmt="o", label=r'$\frac{\tau_t - 1}{\tau -1}$') 
        plt.errorbar(self.c_CV, self.c_calc, xerr= self.c_CV_erro, yerr= self.c_erro_relation, fmt="o", markersize=8,label=r'$\frac{\tau_t - 1}{\tau -1}$') 
        plt.errorbar(self.c_CV, self.c_data, xerr= self.c_CV_erro, yerr= self.c_data_erro, fmt="o", markersize=8, label=r'$\frac{1}{\sigma\nu z}$') 

        plt.legend()
        plt.tight_layout()
        if save:
            plt.savefig('scaling_relation.png', dpi = 150)
        plt.show()
        
        plt.xlabel(r'$\tau$')
        plt.ylabel(r'$\tau_t$')
        #plt.plot(self.c_tau, self.c_tau_t, 'o',markersize=8)
        plt.errorbar(self.c_tau, self.c_tau_t, xerr= self.c_erro_s, yerr= self.c_erro_d, fmt="o", markersize=8)
        plt.plot([1,2.5],[1,2.92],"--",lw=3, label=r"$\frac{\tau_t - 1}{\tau -1} = $ 1.28")
        plt.legend()
        plt.tight_layout()
        if save:
            plt.savefig('tau.png', dpi = 150)
        plt.show()
        
        plt.title('AIC')
        plt.xlabel('<CV>')
        plt.ylabel('<Δ>')
        plt.plot(self.c_CV, self.c_AICs, '-o', markersize=8, label='Size')
        plt.plot(self.c_CV, self.c_AICd, '-o', markersize=8, label='Duration')
        plt.legend()
        plt.tight_layout()
        if save:
            plt.savefig('AIC.png', dpi = 150)
        plt.show()
        
        if self.ADFA:
            plt.xlabel('<CV>')
            plt.ylabel(r'$\alpha$')
            plt.plot(self.c_CV, self.c_DFA, 'o', markersize=8)
            plt.tight_layout()
            plt.savefig('DFA.png', dpi = 150)
            plt.show()
        
    def avalanches(self, samples, bin_size=None, ADFA=False, show=False):
        if (bin_size is None):
            #interspike interval
            #bin a definição de Beggs e Plenz bin=<ISI> 
            bin_size=(samples[-1]-samples[0])/len(samples)

        if show:
            print(bin_size)
        step=np.arange(samples[0],samples[-1],bin_size)
        #print(samples[0],samples[-1])
        spk_count, edges=np.histogram(samples,bins=step)
        I=np.argwhere(spk_count==0).T[0]
        self.duration=np.diff(I)   
        self.duration=self.duration-1;
        self.duration=self.duration[self.duration!=0]
        
        #duration=duration*bin_size
        
        somas=[]
        for i in range(0,len(I)-1):
            somaX=np.sum(spk_count[I[i]:I[i+1]])
            somas.append(somaX);
            #somas=np.concatenate((somas,somaX))
        somas=np.array(somas, dtype=float)
        self.size=somas[somas!=0]
        if ADFA:
            alphas = self.dfa(spk_count, [4,np.log2(len(spk_count))], 0.25)
            return(self.size, self.duration, bin_size, alphas[2])
        else:
            return(self.size, self.duration, bin_size)
        #return(self.size, self.duration, alphas[2])
    
    

    
    def duration(self,samples,bin_size=None): #Retorna as durações de cada avalanche 
        if (bin_size is None):
            #interspike interval
            #bin a definição de Beggs e Plenz bin=<ISI> 
            mfr=len(samples)/(samples[-1]-samples[0])
            bin_size=1.0/mfr

        step=np.arange(samples[0],samples[-1],bin_size)
        spk_count, edges=np.histogram(samples,bins=step)
        I=np.argwhere(spk_count==0).T[0]
        self.duration=np.diff(I)   
        self.duration=self.duration[self.duration!=0]
        #duration=duration*bin_size
        print("Duração de avalanches calculados com sucesso!")
        print('bin usado para duração = ', bin_size)
        return(self.duration)    

    def size (self, samples, bin_size=None):
        if (bin_size is None):
            #interspike interval
            #bin a definição de Beggs e Plenz bin=<ISI> 
            mfr=len(samples)/(samples[-1]-samples[0])
            bin_size=1.0/mfr  
        step=np.arange(samples[0],samples[-1],bin_size)
        spk_count, edges=np.histogram(samples,bins=step)
        I=np.argwhere(spk_count==0).T[0]
        somas=[]
        for i in range(0,len(I)-1):
            somaX=np.sum(spk_count[I[i]:I[i+1]])
            somas.append(somaX);
            #somas=np.concatenate((somas,somaX))
        somas=np.array(somas, dtype=float)
        self.size=somas[somas!=0]
        print("Tamanho de avalanches calculados com sucesso!")
        print('bin usado para tamanho = ', bin_size)
        return(self.size)
    
    
    def calc_rms(self, x, scale):
        """
        windowed Root Mean Square (RMS) with linear detrending.

        Args:
        -----
          *x* : numpy.array
            one dimensional data vector
          *scale* : int
            length of the window in which RMS will be calculaed
        Returns:
        --------
          *rms* : numpy.array
            RMS data in each window with length len(x)//scale
        """
        # making an array with data divided in windows
        shape = (x.shape[0]//scale, scale)
        X = np.lib.stride_tricks.as_strided(x,shape=shape)
        # vector of x-axis points to regression
        scale_ax = np.arange(scale)
        rms = np.zeros(X.shape[0])
        for e, xcut in enumerate(X):
            coeff = np.polyfit(scale_ax, xcut, 1)
            xfit = np.polyval(coeff, scale_ax)
            # detrending and computing RMS of each window
            rms[e] = np.sqrt(np.mean((xcut-xfit)**2))
        return rms

    def dfa(self, x,scale_lim, scale_dens, show=False):
        """
        Detrended Fluctuation Analysis - algorithm with measures power law
        scaling of the given signal *x*.
        More details about algorithm can be found e.g. here:
        Hardstone, R. et al. Detrended fluctuation analysis: A scale-free
        view on neuronal oscillations, (2012).

        Args:
        -----
        *x* : numpy.array
            one dimensional data vector
        *scale_lim* = [5,9] : list of lenght 2
            boundaries of the scale where scale means windows in which RMS
            is calculated. Numbers from list are indexes of 2 to the power
            of range.
        *scale_dens* = 0.25 : float
            density of scale divisions
        *show* = False
            if True it shows matplotlib picture
        Returns:
        --------
        *scales* : numpy.array
            vector of scales
        *fluct* : numpy.array
            fluctuation function
        *alpha* : float
            DFA exponent
        """
        # cumulative sum of data with substracted offset

        #print (scale_lim)

        scale_lim = [2, 7]
        y = np.cumsum(x - np.mean(x))
        scales = (2**np.arange(scale_lim[0], scale_lim[1], scale_dens)).astype(np.int)

        fluct = np.zeros(len(scales))
        # computing RMS for each window
        for e, sc in enumerate(scales):
            fluct[e] = np.mean(np.sqrt(self.calc_rms(y, sc)**2))
        # fitting a line to rms data
        #print(scales)
        #I = np.argwhere(scales<=300).T[0]
        coeff = np.polyfit(np.log2(scales), np.log2(fluct), 1)
            #coeff = np.polyfit(np.log2(x[I]), np.log2(y_avg[I]), 1) 	 
            #fluctfit = 2**np.polyval(coeff,np.log2(x[I]))

        if show:
            fluctfit = 2**np.polyval(coeff,np.log2(scales))
            plt.loglog(scales, fluct, 'bo')
            plt.loglog(scales, fluctfit, 'r', label=r'$\alpha$ = %0.2f'%coeff[0])
            plt.title('DFA')
            plt.xlabel(r'$\log_{10}$(time window)')
            plt.ylabel(r'$\log_{10}$<F(t)>')
            plt.legend()
            plt.show()
        return (scales, fluct, coeff[0],coeff)

    def mle(self, x, xmin=None, xmax=None):
        if (xmin==None):
            xmin = np.min(x)
        if (xmax==None):
            xmax = np.max(x)
        tauRange=np.array([1,12])
        precision = 10**(-3)
        # Error check the precision
        if math.log10(precision) != round(math.log10(precision)):
            print('The precision must be a power of ten.')

        x = np.reshape(x, len(x))

      #Determine data type
        if np.count_nonzero(np.absolute(x - np.round(x)) > 3*(np.finfo(float).eps)) > 0:
            dataType = 'CONT'
        else:
            dataType = 'INTS'
            x = np.round(x)

       # print(dataType)
        #Truncate
        z = x[(x>=xmin) & (x<=xmax)]
        unqZ = np.unique(z)
        nZ = len(z)
        nUnqZ = len(unqZ)
        allZ = np.arange(xmin,xmax+1)
        nallZ = len(allZ)

        #MLE calculation

        r = xmin / xmax
        nIterations = int(-math.log10(precision))

        for iIteration in range(1, nIterations+1):

            spacing = 10**(-iIteration)

            if iIteration == 1:
                taus = np.arange(tauRange[0], tauRange[1]+spacing, spacing)

            else: 
                if tauIdx == 0:
                    taus = np.arange(taus[0], taus[1]+spacing, spacing)
                    #return (taus,0,0,0)
                elif tauIdx == len(taus):    
                    taus = np.arange(taus[-2], taus[-1]+spacing, spacing)#####
                else:
                    taus = np.arange(taus[tauIdx-1], taus[tauIdx+1]+spacing, spacing)

            #return(dataType)        
            nTaus = len(taus)

            if dataType=='INTS':
                #replicate arrays to equal size
                allZMat = np.matlib.repmat(np.reshape(allZ,(nallZ,1)),1,nTaus)
                tauMat = np.matlib.repmat(taus,nallZ,1)

                #compute the log-likelihood function
                #L = - np.log(np.sum(np.power(allZMat,-tauMat),axis=0)) - (taus/nZ) * np.sum(np.log(z))
                L = - nZ*np.log(np.sum(np.power(allZMat,-tauMat),axis=0)) - (taus) * np.sum(np.log(z))


            elif dataType=='CONT':
                #return (taus,r, nZ,z)
                L = np.log( (taus - 1) / (1 - r**(taus - 1)) )- taus * (1/nZ) * np.sum(np.log(z)) - (1 - taus) * np.log(xmin)

                if numpy.in1d(1,taus):
                    L[taus == 1] = -np.log(np.log(1/r)) - (1/nZ) * np.sum(np.log(z))
            tauIdx=np.argmax(L)

        tau = taus[tauIdx]
        Ln = L[tauIdx]
        #return (taus,L,tau
        return (tau, Ln)
    
    def lnormal1(self, x, xmin, xmax):
        x = x[(x>=xmin) & (x<=xmax)]
        nX = len(x)
        μ = (sum(np.log(x)))/nX
        sig2 = (sum((np.log(x)-μ)**2))/nX
        σ=np.sqrt(sig2)
        A=σ*np.sqrt(2*math.pi)
        B=2*(σ**2)
        C = np.sum(np.log(x*A)+((np.log(x)-μ)**2)/B)
        l = -C
        return(l) 
    
    def size_avg(self, duracao, tamanho_ss):

        x=np.unique(duracao)
        y_avg=np.zeros_like(x)
        y_err=np.zeros_like(x)
        for j in range (0,len(x)):
            y_avg[j] = np.mean(tamanho_ss[duracao==x[j]])
            y_err[j] = np.std(tamanho_ss[duracao==x[j]])
        return(x,y_avg) 

    def fitt(self, time,size_avg,xmin,xmax):    
        logx = np.log10(time[(time>=xmin) & (time<=xmax)])
        logy = np.log10(size_avg[(time>=xmin) & (time<=xmax)])
        p, V = np.polyfit(logx, logy, 1, cov=True)
        error = np.sqrt(V[0][0])
        index = p[0]
        return(index,error)
    
    
#////////////////////////////////////////////////////////////////////////////////////
######################       Espacio de trabajo          ############################
#####################################################################################    

# data=np.transpose(np.loadtxt("datos_cv/datas/time_cluster_shank_sua_2020Mar07.txt"))
# cv_time=10
# wt=50
# bin_cv=0.030
# tcb2=NAT2(data[0,:])

# CV,mspike=tcb2.CV_calc(cv_time,bin_cv, show=False, save=False, ADFA=False, 
#                         showbin=False, bin_size=None)
# px=np.where(np.array(mspike,dtype=object))[0]

# xspk=np.array(mspike,dtype=object)

# mx=tcb2.CV
# cv_sort=np.array(sorted(mx.copy()))[np.array(sorted(mx.copy()))<3]
# zpy=np.where(np.array(mx,dtype=object))[0]

# nNwt=int(np.round(len(cv_sort)/wt))

# matriz_cv=np.zeros((len(mx),2))
# cw=0; dtcv_sort=[]

# for i in range(len(cv_sort)):
#     for j in range(len(cv_sort)):    
#         if(cv_sort[i] == mx[j]):
#             dtcv_sort.append(zpy[j])
#             break
        
# average_CV=[]
# for cvi in range(nNwt):
#     a=wt*(cvi)
#     b=wt*(cvi+1)
#     cv_tiepo_wt=[]; avercv=[]
#     for tk in range(a,b):
#         avercv.append(cv_sort[tk])
#         xspk[dtcv_sort[tk]]
#         cv_tiepo_wt += xspk[dtcv_sort[tk]]
#     average_CV.append(np.mean(avercv))

        
#     cv_tiepo_wt.sort()    
#     for o in cv_tiepo_wt:
#         o
#         #print(o)
# for i in average_CV:
#     print(i)
       
