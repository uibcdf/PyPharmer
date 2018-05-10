import nglview
import mdtraj
import numpy as np

color_code={ 'Hydrophobic':        [0.100, 1.000, 0.000],\
             'HydrogenAcceptor':   [1.000, 0.84, 0.000],\
             'HydrogenDonor':      [1.000, 1.000, 1.000],\
             'Aromatic':           [0.627, 0.1254, 0.941],\
             'NegativeIon':        [1.00, 0.00, 0.00],\
             'PositiveIon':        [0.00, 0.00, 1.00],\
             'Other':              [0.74, 0.74, 0.74],\
             'PhenylalanineAnalog':[1.0, 1.0, 0.0],\
             'LeuValAnalog':       [1.0, 1.0, 0.0] \
             }

def get_nglview(query,arrow_norm=2.0,arrow_radius=0.2):

    tmp_view=nglview.NGLWidget()
    arrow_norm=np.sqrt(arrow_norm)

    for point in query.points:
        color=color_code[point.name]
        tmp_view.shape.add_sphere(point.position.tolist(),color,point.radius,point.name)
        if point.svector is not None:
            source_array=point.position
            end_array=(point.position+arrow_norm*point.svector)
            tmp_view.shape.add_arrow(source_array.tolist(),end_array.tolist(),color,arrow_radius)

    if receptor in query.receptor:
        tmp_view.add_cartoon (color='blue')
        tmp_view.add_surface (color='gray', opacity=0.5)

    return tmp_view
