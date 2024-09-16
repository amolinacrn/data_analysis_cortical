import numpy as np
import plotly.graph_objects as go
import networkx as nx
import os
import numpy
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors
import seaborn as sns
import test
import bct.algorithms as reference
from bct.utils import BCTParamError, binarize, get_rng


address_mc="MC_normalizada_pesada.txt"


recablear=55300

tcm = np.loadtxt(address_mc,dtype=float)


mxramdom = reference.randmio_und_connected(tcm,recablear)

#mxlatice = reference.latmio_und(tcm,recablear)

mxlatice = reference.latmio_und_connected(tcm,recablear)

#tmc_latice = mxramdom[0]
 
tmc_latice = mxlatice[0]
tmc_ramdom = mxramdom[0]

n_mc=len(tmc_latice)

############ escriba la matriz binarizada ################

file_txt="tcm_lattice.txt"
ffile_txt="tcm_ramdom.txt"
        
fdata = open(file_txt, 'w')
ffdata = open(ffile_txt, 'w')

try:
    for zz_i in range(n_mc):
        for zz_j in range(n_mc):
            fdata.write( str(tmc_latice[zz_i][zz_j])+"\t")
        fdata.write("\n")
                   
finally:
    fdata.close()
    
try:           
    for z_i in range(n_mc):
        for z_j in range(n_mc):
            ffdata.write( str(tmc_ramdom[z_i][z_j])+"\t")
        ffdata.write("\n")
finally:
    ffdata.close()

######################################################

#print(tmc_latice[0][1])

