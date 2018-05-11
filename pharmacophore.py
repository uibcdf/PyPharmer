import numpy as np
import simplejson as json
from .viewer import get_nglview as _get_nglview

class Receptor:
    pass

class Ligand:
    pass

class Anchor:

    def __init__(self,json_object=None):

        self.name     = None
        self.wiggle   = None
        self.enabled  = True
        self.triangle = None

        if json_object is not None:

            tmp_json_keys = json_object.keys()

            self.name     = json_object.get('name')

            if 'wiggle' in tmp_json_keys:
                self.wiggle   = json_object.get('wiggle')

            if 'enabled' in tmp_json_keys:
                self.enabled  = json_object.get('enabled')

            x1 = np.array([json_object.get('x1'),json_object.get('y1'),json_object.get('z1')])
            x2 = np.array([json_object.get('x2'),json_object.get('y2'),json_object.get('z2')])
            x3 = np.array([json_object.get('x3'),json_object.get('y3'),json_object.get('z3')])

            self.triangle = np.v_stack([x1,x2,x3])

        pass

class Point:

    def __init__(self,json_object=None):

        self.name       = None
        self.position   = None
        self.radius     = None
        self.svector    = None
        self.vector     = None
        self.enabled    = True

        if json_object is not None:

            tmp_json_keys = json_object.keys()

            self.name     = json_object.get('name')
            self.position = np.array([json_object.get('x'),json_object.get('y'),json_object.get('z')])
            self.radius   = json_object.get('radius')

            if 'enabled' in tmp_json_keys:
                self.enabled = json_object.get('enabled')

            if 'svector' in tmp_json_keys:
                tmp_vector   = json_object.get('svector')
                if tmp_vector is not None:
                    self.svector = np.array([tmp_vector.get('x'),tmp_vector.get('y'),tmp_vector.get('z')])

            if 'vector' in tmp_json_keys:
                tmp_vector   = json_object.get('vector')[0]
                self.vector = np.array([tmp_vector.get('x'),tmp_vector.get('y'),tmp_vector.get('z')])


        pass

class Pharmacophore:

    input_file     = None

    anchors        = []
    points         = []

    with_receptor  = False
    with_ligand    = False
    _json          = None

    def __init__(self, input_file=''):

        self.input_file    = input_file

        self._json          = json.loads(open(input_file).read())

        tmp_json           = self._json
        tmp_json_keys      = tmp_json.keys()

        # If there is a receptor to be loaded (and we want it).
        if 'receptor' in tmp_json_keys:
            self.with_receptor = True
            self.receptor      = tmp_json.get('receptor')


        # If there is a ligand to be loaded (and we want it).
        if 'sdf' in tmp_json_keys:
                self.with_ligand   = True
                self.ligand        = tmp_json.get('sdf')
                #sdf = urllib.unquote(sdf) 
                ###write it out to a temp file
                #tfout = tempfile.NamedTemporaryFile(suffix = '.pdb', delete=False)
                #tmpname = tfout.name 
                #tfout.write(sdf)
                #tfout.close()

        if 'ligand' in tmp_json_keys:
                self.with_ligand   = True
                self.ligand        = tmp_json.get('ligand')

        # If there are anchors
        if 'anchors' in tmp_json_keys:
            for json_object in tmp_json.get('anchors'):
                self.anchors.append(Anchor(json_object))

        # List of points
        if 'points' in tmp_json_keys:
            for json_object in tmp_json.get('points'):
                self.points.append(Point(json_object))

        pass

    def get_view(self,viewer='nglview',arrow_norm=2.0,arrow_radius=0.2):

        if viewer == 'nglview':
            return _get_nglview(self,arrow_norm,arrow_radius)

        pass
