import nglview
import mdtraj
import matplotlib as mt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
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
        rec.add_line ('protein')

    if pharmacophore.with_ligand and ligand:
        view_with_ligand=True
        tmp_ligand_file = tempfile.NamedTemporaryFile(mode='w+',delete=False,suffix='.mol2')
        tmp_ligand_file.write(pharmacophore.ligand)
        tmp_ligand_file.close()
        tmp_view.add_component(tmp_ligand_file.name)

    for point in pharmacophore.points:
        color=color_code[point.name]
        tmp_sphere=tmp_view.shape.add_sphere (point.position.tolist(),color,point.radius,point.name)
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

def scatter_plot ():

    Hydrophobic         = [[],[],[]]
    HydrogenAcceptor    = [[],[],[]]
    HydrogenDonor       = [[],[],[]]
    Aromatic            = [[],[],[]]
    NegativeIon         = [[],[],[]]
    PositiveIon         = [[],[],[]]
    InclusionSphere     = [[],[],[]]
    Other               = [[],[],[]]
    PhenylalanineAnalog = [[],[],[]]
    LeuValAnalog        = [[],[],[]]

    for element in self.points:

        if element[0]==1.0:
            Hydrophobic[0].append (element[1])
            Hydrophobic[1].append (element[2])
            Hydrophobic[2].append (element[3])
        if element[0]==2.0:
            HydrogenAcceptor[0].append (element[1])
            HydrogenAcceptor[1].append (element[2])
            HydrogenAcceptor[2].append (element[3])
        if element[0]==3.0:
            HydrogenDonor[0].append (element[1])
            HydrogenDonor[1].append (element[2])
            HydrogenDonor[2].append (element[3])
        if element[0]==4.0:
            Aromatic[0].append (element[1])
            Aromatic[1].append (element[2])
            Aromatic[2].append (element[3])
        if element[0]==4.0:
            NegativeIon[0].append (element[1])
            NegativeIon[1].append (element[2])
            NegativeIon[2].append (element[3])
        if element[0]==5.0:
            PositiveIon[0].append (element[1])
            PositiveIon[1].append (element[2])
            PositiveIon[2].append (element[3])
        if element[0]==6.0:
            Other[0].append (element[1])
            Other[1].append (element[2])
            Other[2].append (element[3])
        if element[0]==7.0:
            PhenylalanineAnalog[0].append (element[1])
            PhenylalanineAnalog[1].append (element[2])
            PhenylalanineAnalog[2].append (element[3])
        if element[0]==8.0:
            LeuValAnalog[0].append (element[1])
            LeuValAnalog[1].append (element[2])
            LeuValAnalog[2].append (element[3])

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    ax.scatter(Hydrophobic[0],Hydrophobic[1],Hydrophobic[2], color=[0.100, 1.000, 0.000])
    ax.scatter(HydrogenAcceptor[0],HydrogenAcceptor[1],HydrogenAcceptor[2], color=[1.000, 0.84, 0.000])
    ax.scatter(HydrogenDonor[0],HydrogenDonor[1],HydrogenDonor[2], color='black') #[1.000, 1.000, 1.000]
    ax.scatter(Aromatic[0],Aromatic[1],Aromatic[2], color=[0.627, 0.1254, 0.941])
    ax.scatter(NegativeIon[0],NegativeIon[1],NegativeIon[2], color=[1.00, 0.00, 0.00])
    ax.scatter(PositiveIon[0],PositiveIon[1],PositiveIon[2], color=[0.00, 0.00, 1.00])
    ax.scatter(Other[0],Other[1],Other[2], color=[0.74, 0.74, 0.74])
    ax.scatter(PhenylalanineAnalog[0],PhenylalanineAnalog[1],PhenylalanineAnalog[2], color=[1.0, 1.0, 0.0])
    ax.scatter(LeuValAnalog[0],LeuValAnalog[1],LeuValAnalog[2], color=[1.0, 1.0, 0.0])
    ax.set_xlabel('X Coordinates')
    ax.set_ylabel('Y Coordinates')
    ax.set_zlabel('Z Coordinates')
    return plt.show ()


def cliques ():
    plt.subplot(121)
    nx.draw_networkx (GRAPH, pos=nx.random_layout(GRAPH), with_labels=False)
    return plt.show()
