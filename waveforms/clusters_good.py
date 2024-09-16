##codigo para importar datos
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import pandas as pd
import math
import cmath
import csv
import numpy
import re

#convertir datos a csv


for i in range(1,7):
    datas_dat="datas.txt"
   
    datas_dat="2020Jan29/g0"+str(i)+"_cluster_info.tsv"
    dat_csv="datos_info_csv/g0"+str(i)+"_cluster_info.csv"
    #dat_csv="waveforms_ALL_g0.csv"
    #datas_dat="cluster_info.csv"
    with open(datas_dat, 'r') as myfile: 
        with open(dat_csv, 'w') as csv_file:
            for line in myfile:
                fileContent = re.sub("\t", ",", line)
                csv_file.write(fileContent)
    print("Successfully made csv file")



       




