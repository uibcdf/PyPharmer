import numpy as np
import simplejson as json
from pprint import pprint
import networkx as nx
import itertools as it
from .viewer import  scatter_plot as _scatter_plot


# Descriptors list:

descriptors_list={ 'Hydrophobic':        1.0,\
             'HydrogenAcceptor':    2.0,\
             'HydrogenDonor':       3.0,\
             'Aromatic':            4.0,\
             'NegativeIon':         5.0,\
             'PositiveIon':         6.0,\
             'InclusionSphere':     7.0,\
             'Other':               8.0,\
             'PhenylalanineAnalog': 9.0,\
             'LeuValAnalog':        10.0 \
             }
color_code={ 1.0:        [0.100, 1.000, 0.000],\
             2.0:   [1.000, 0.84, 0.000],\
             3.0:      'gray',\
             4.0:           [0.627, 0.1254, 0.941],\
             5.0:        [1.00, 0.00, 0.00],\
             6.0:        [0.00, 0.00, 1.00],\
             7.0:    [0.00, 1.00, 1.00],\
             8.0:              [0.74, 0.74, 0.74],\
             9.0:[1.0, 1.0, 0.0],\
             10.0:       [1.0, 1.0, 0.0] \
             }

## This funtion was added to generate clique graph.
def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a
class Descriptors:

    def __init__(self,input_file='',json_object=None, select='enabled'):

        self.points         = None
        self.input_file     = input_file
        self.json_object    = None
        #select            = should be 'enabled','disabled' or 'all'
        self._json          = json.loads(open(input_file).read())
        tmp_json            = self._json
        tmp_json_keys       = tmp_json.keys()

        if select is 'enabled':

            if 'points' in tmp_json_keys:

                descriptors=[]

                for json_object in tmp_json.get ('points'):

                    if json_object.get ('enabled')==True:
                        name=descriptors_list[json_object.get ('name')]
                        vector=json_object.get ('vector')
                        vector_on=json_object.get ('vector_on')

                        if 'vector' and 'vector_on' in json_object.keys ():
                            if vector_on ==1:
                                for element in vector:
                                    descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'), element['x'], element['y'],element['z']])


                            if vector_on ==0:
                                vector=[0]
                                descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),vector[0],vector[0],vector[0]])

                            descriptors.append (descriptor)

                        if 'vector' and not 'vector_on' in json_object.keys ():
                            if vector!=None:
                                for element in vector:
                                    descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'), element['x'], element['y'],element['z']])

                            if vector==None:
                                vector=[0]
                                descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'), vector[0],vector[0],vector[0]])

                            descriptors.append (descriptor)

                self.points=np.array (descriptors)

        if select is 'disabled':

            if 'points' in tmp_json_keys:

                descriptors=[]

                for json_object in tmp_json.get ('points'):

                    if json_object.get ('enabled')==False:
                        name=descriptors_list[json_object.get ('name')]
                        vector=json_object.get ('vector')
                        vector_on=json_object.get ('vector_on')

                        if 'vector' and 'vector_on' in json_object.keys ():
                            if vector_on ==1:
                                for element in vector:
                                    descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),element['x'], element['y'],element['z']])

                            if vector_on ==0:
                                vector=[0]
                                descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),vector[0],vector[0],vector[0]])

                            descriptors.append (descriptor)

                        if 'vector' and not 'vector_on' in json_object.keys ():
                            if vector!=None:
                                for element in vector:
                                    descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),element['x'], element['y'],element['z']])

                            if vector==None:
                                vector=[0]
                                descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),vector[0],vector[0],vector[0]])

                            descriptors.append (descriptor)

                self.points=np.array (descriptors)

        if select is 'all':

            if 'points' in tmp_json_keys:

                descriptors=[]

                for json_object in tmp_json.get ('points'):
                    name=descriptors_list[json_object.get ('name')]
                    vector=json_object.get ('vector')
                    vector_on=json_object.get ('vector_on')

                    if 'vector' and 'vector_on' in json_object.keys ():
                        if vector_on ==1:
                            for element in vector:
                                descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'), element['x'], element['y'],element['z']])

                        if vector_on ==0:
                            vector=[0]
                            descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'), vector[0],vector[0],vector[0]])

                        descriptors.append (descriptor)

                    if 'vector' and not 'vector_on' in json_object.keys ():
                        if vector!=None:
                            for element in vector:
                                descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'), element['x'], element['y'],element['z']])

                        if vector==None:
                            vector=[0]
                            descriptor=([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'), vector[0],vector[0],vector[0]])

                        descriptors.append (descriptor)

                self.points=np.array (descriptors)

        pass

    def plot_points ():
        return _scatter_plot ()
    pass

    def get_cliques ():
        G= nx.Graph (pharmacophore='descriptors')
        for i in range (len(self.points)):
            node=totuple(self.points[i])
            G.add_node (node_for_adding=node)
            color=color_code[node[0]]
            result=list(it.combinations(G.nodes,2))
            for j in list(result):
                G.add_edge (j[0],j[1])
        return nx.draw (G,pos=nx.spectral_layout(G),node_color=color)
        print ('File:',file,' ','number or nodes:',G.number_of_nodes(),' ','Number of edges:',G.number_of_edges ())
    pass
