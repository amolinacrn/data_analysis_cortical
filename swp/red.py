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

def redtop(mx_conectividad):

    ################################################################
    ############### grafico de red topologica  #####################
    ################################################################

    G = nx.from_numpy_matrix(mx_conectividad) #matriz_conectividad)
    pos = nx.spring_layout(G)
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
            colorscale='plasma',
            #reversescale=True,
            color=[],
            size=9,
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
        
    node_trace.marker.color = node_adjacencies
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
    #fig.savefig('file.png')
def graf_mc_conectividad(mx_conectividad):
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
    cax = plt.axes([0.81, 0.1, 0.03, 0.82])
    
    plt.colorbar(cax=cax)
    #plt.tight_layout()
    plt.show()

adr_mc=np.loadtxt("mwc.txt",dtype=float)

redtop(adr_mc)
graf_mc_conectividad(adr_mc)