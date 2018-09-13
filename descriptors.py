import numpy as np
import simplejson as json
from pprint import pprint
import networkx as nx
import itertools as it
#from .viewer import  scatter_plot as _scatter_plot
from .analysis import get_networkx as _get_networkx

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
                        vector=json_object.get ('vector')
                        vector_on=json_object.get ('vector_on')

                        if 'vector' and 'vector_on' in json_object.keys ():
                            if vector_on ==1:
                                for element in vector:
                                    descriptor=([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                    element['x'], element['y'],element['z']])

                            if vector_on ==0:
                                vector=[0]
                                descriptor=([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                vector[0],vector[0],vector[0]])
                            descriptors.append (descriptor)

                        if 'vector' and not 'vector_on' in json_object.keys ():
                            if vector!=None:
                                for element in vector:
                                    descriptor=n([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                    element['x'], element['y'],element['z']])

                            if vector==None:
                                vector=[0]
                                descriptor=([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                vector[0],vector[0],vector[0]])
                            descriptors.append (descriptor)

                self.points=np.array (descriptors,dtype=object)

        if select is 'disabled':

            if 'points' in tmp_json_keys:

                descriptors=[]

                for json_object in tmp_json.get ('points'):

                    if json_object.get ('enabled')==False:
                        vector=json_object.get ('vector')
                        vector_on=json_object.get ('vector_on')

                        if 'vector' and 'vector_on' in json_object.keys ():
                            if vector_on ==1:
                                for element in vector:
                                    descriptor=np.array([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                     element['x'], element['y'],element['z']],dtype=object)

                            if vector_on ==0:
                                vector=[0]
                                descriptor=np.array([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                vector[0],vector[0],vector[0]],dtype=object)

                            descriptors.append (descriptor)

                        if 'vector' and not 'vector_on' in json_object.keys ():
                            if vector!=None:
                                for element in vector:
                                    descriptor=np.array([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                     element['x'], element['y'],element['z']],dtype=object)

                            if vector==None:
                                vector=[0]
                                descriptor=np.array([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                vector[0],vector[0],vector[0]],dtype=object)

                            descriptors.append (descriptor)

                self.points=np.array (descriptors)

        if select is 'all':

            if 'points' in tmp_json_keys:

                descriptors=[]

                for json_object in tmp_json.get ('points'):
                    vector=json_object.get ('vector')
                    vector_on=json_object.get ('vector_on')

                    if 'vector' and 'vector_on' in json_object.keys ():
                        if vector_on ==1:
                            for element in vector:
                                descriptor=np.array([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                 element['x'], element['y'],element['z']],dtype=object)

                        if vector_on ==0:
                            vector=[0]
                            descriptor=np.array([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                            vector[0],vector[0],vector[0]],dtype=object)

                        descriptors.append (descriptor)

                    if 'vector' and not 'vector_on' in json_object.keys ():
                        if vector!=None:
                            for element in vector:
                                descriptor=np.array([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                                 element['x'], element['y'],element['z']],dtype=object)


                        if vector==None:
                            vector=[0]
                            descriptor=np.array([json_object.get('name'),json_object.get('x'),json_object.get('y'),json_object.get('z'),json_object.get ('radius'),
                            vector[0],vector[0],vector[0]],dtype=object)

                        descriptors.append (descriptor)

                self.points=np.array (descriptors)

    def get_cliques (self,algorithm='networkx'):
        if algorithm == 'networkx':
            return _get_networkx(self,algorithm='networkx')
        pass

    def plot ():
        return _plot_cliques ()
    pass

    def plot_points ():
        return _scatter_plot ()
    pass
