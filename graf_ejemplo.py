import numpy as np
from igraph import *
import igraph as ig
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
g = ig.Graph([(0,1),(1,1),(1,2)])
visual_style = {}
g.get_adjacency()
layout = g.layout("kk")
g.vs["label"] = range(3)
visual_style["vertex_size"] = 10
color_dict = {"m": "blue", "f": "pink"}
visual_style["layout"] = layout
visual_style["bbox"] = (300, 300)
visual_style["margin"] = 20


ig.plot(g, layout=layout, target=ax)
plt.show()
