import nglview
import mdtraj
import numpy as np
import tempfile
import os

color_code={ 'Hydrophobic':        [0.100, 1.000, 0.000],\
             'HydrogenAcceptor':   [1.000, 0.84, 0.000],\
             'HydrogenDonor':      [1.000, 1.000, 1.000],\
             'Aromatic':           [0.627, 0.1254, 0.941],\
             'NegativeIon':        [1.00, 0.00, 0.00],\
             'PositiveIon':        [0.00, 0.00, 1.00],\
             'InclusionSphere':    [0.00, 1.00, 1.00],\
             'Other':              [0.74, 0.74, 0.74],\
             'PhenylalanineAnalog':[1.0, 1.0, 0.0],\
             'LeuValAnalog':       [1.0, 1.0, 0.0] \
             }
# The color for InclusionSphere was given at random



def get_nglview(pharmacophore,receptor=True,ligand=True,arrow_norm=2.0,arrow_radius=0.2):

    view_with_ligand=False
    view_with_receptor=False

    tmp_view=nglview.NGLWidget()
    arrow_norm=np.sqrt(arrow_norm)

    if pharmacophore.with_receptor and receptor:
        view_with_receptor=True
        tmp_receptor_file = tempfile.NamedTemporaryFile(mode='w+',delete=False,suffix='.pdb')
        tmp_receptor_file.write(pharmacophore.receptor)
        tmp_receptor_file.close()
        rec=tmp_view.add_component(tmp_receptor_file.name)
        #rec.add_surface ('protein', opacity=0.4)
        rec.add_line ('protein')

    if pharmacophore.with_ligand and ligand:
        view_with_ligand=True
        tmp_ligand_file = tempfile.NamedTemporaryFile(mode='w+',delete=False,suffix='.mol2')
        tmp_ligand_file.write(pharmacophore.ligand)
        tmp_ligand_file.close()
        tmp_view.add_component(tmp_ligand_file.name)

    for point in pharmacophore.points:
        color=color_code[point.name]
        tmp_view.shape.add_sphere(point.position.tolist(),color,point.radius,point.name)
        if point.svector is not None:
            source_array=point.position
            end_array=(point.position+arrow_norm*point.svector)
            tmp_view.shape.add_arrow(source_array.tolist(),end_array.tolist(),color,arrow_radius)

    if view_with_receptor:
        tmp_receptor_file.close()
        os.unlink(tmp_receptor_file.name)

    if view_with_ligand:
        tmp_ligand_file.close()
        os.remove(tmp_ligand_file.name)


    return tmp_view
