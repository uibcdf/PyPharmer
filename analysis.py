import numpy as np
import simplejson as json


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

class Descriptors:
    input_file  = None
    enabled     = True
    all_points  = False
    _json       = None

    def __init__(self,input_file='',json_object=None, enabled=True, all_points=False):

        self.input_file    = input_file
        self._json         = json.loads(open(input_file).read())
        tmp_json           = self._json
        tmp_json_keys      = tmp_json.keys()

        if enabled is True:

             if 'points' in tmp_json_keys:
                for json_object in tmp_json.get ('points'):

                    if json_object.get ('enabled')==True:
                        name=descriptors_list [json_object.get ('name')]
                        vector=json_object.get ('vector')
                        vector_on=json_object.get ('vector_on')
            
                        if 'vector' and ('vector_on') in json_object.keys ():
                            if vector_on ==1:
                                descriptor = ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'), vector])
                                print (descriptor)
        
                        if 'vector' and 'vector_on' in json_object.keys ():
                            if vector_on ==0:
                                vector=[0]
                                descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),vector[0],vector[0],vector[0]])
                                print (descriptor)
            
                        if 'vector' and not 'vector_on' in json_object.keys ():
                            if vector!=None:
                                descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'), vector])
                                print (descriptor)
                            if vector==None:
                                vector=[0]
                                descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),vector[0],vector[0],vector[0]])
                                print (descriptor)
        if enabled is False:

            if 'points' in tmp_json_keys:
                for json_object in tmp_json.get ('points'):

                    if json_object.get ('enabled')==False:
                        name=descriptors_list [json_object.get('name')]
                        vector=json_object.get ('vector')
                        vector_on=json_object.get ('vector_on')
            
                        if 'vector' and ('vector_on') in json_object.keys ():
                            if vector_on ==1:
                                descriptor = ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'), vector])
                                print (descriptor)
        
                        if 'vector' and 'vector_on' in json_object.keys ():
                            if vector_on ==0:
                                vector=[0]
                                descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),vector[0],vector[0],vector[0]])
                                print (descriptor)
            
                        if 'vector' and not 'vector_on' in json_object.keys ():
                            if vector!=None:
                                descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'), vector])
                                print (descriptor)
                            if vector==None:
                                vector=[0]
                                descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),vector[0],vector[0],vector[0]])
                                print (descriptor)

        if all_points is True:

            if 'points' in tmp_json_keys:
                for json_object in tmp_json.get ('points'):

                    name=descriptors_list[json_object.get('name')]
                    vector=json_object.get ('vector')
                    vector_on=json_object.get ('vector_on')
            
                    if 'vector' and ('vector_on') in json_object.keys ():
                        if vector_on ==1:
                            descriptors = ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'), vector])
                            print (descriptor)
        
                    if 'vector' and 'vector_on' in json_object.keys ():
                        if vector_on ==0:
                            vector=[0]
                            descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),vector[0],vector[0],vector[0]])
                            print (descriptor)
            
                    if 'vector' and not 'vector_on' in json_object.keys ():
                        if vector!=None:
                            descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'), vector])
                            print (descriptor)
                        if vector==None:
                            vector=[0]
                            descriptor= ([name,json_object.get('x'),json_object.get('y'),json_object.get('z'),vector[0],vector[0],vector[0]])
                            print (descriptor)

        if all_points is False:

            pass
        

        pass






class Analyze:
    pass



