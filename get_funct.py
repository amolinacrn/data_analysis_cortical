#!/usr/bin/env python
# -*- coding: utf-8 -*

##/////////////////////////////////////////////////////##
##     Por Miguel Alejandro Molina  03/11/2022         ##
##/////////////////////////////////////////////////////##

#------------------------------------------------------------------------------------------


import numpy as np
import os
import numpy
import math
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import seaborn as sns
import plotly.graph_objects as go
import networkx as nx
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors
import seaborn as sns
#import ROOT # paquete root/cern debe instalarse por separado
import statistics
from scipy.signal import savgol_filter
from scipy.interpolate import make_interp_spline
from scipy import stats
from bct.algorithms import number_of_components
import re
from os import remove


def grafico_swp(mx):

    G = nx.from_numpy_matrix(mx) #matriz_conectividad)
    pos = nx.kamada_kawai_layout(G)
    for n, p in pos.items():
        G.nodes[n]['pos'] = p

    edge_x = []
    edge_y = []

    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.3, color='black'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='viridis',
            #reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'),
            line_width=0.4))
            
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))
        
    node_trace.marker.color = "blue"#node_adjacencies
    node_trace.text = node_text

    title = "Network Graph Demonstration"
    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(height=800, width=800,
                    paper_bgcolor='white',
                    plot_bgcolor='white',
                    #title=title,
                    xaxis=dict(showgrid=False, zeroline=False,
                            showticklabels=False, mirror=True),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, mirror=True)
                    ))

    fig.show()

adm = np.genfromtxt("SWP.txt")
  

grafico_swp(adm)