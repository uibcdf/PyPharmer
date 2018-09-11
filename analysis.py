import numpy as np
import simplejson as json
from pprint import pprint
import networkx as nx
import itertools as it

def get_networkx (self, algorithm='networkx'):
    self.graph    = nx.Graph (pharmacophore='descriptors')
    for point  in self.points:
        node=(point[1],point[2],point[3],point[4])
        self.graph.add_node (node_for_adding=node, name=point[0])
        result=list(it.combinations(self.graph.nodes,2))
        for j in list(result):
            self.graph.add_edge (j[0],j[1])
