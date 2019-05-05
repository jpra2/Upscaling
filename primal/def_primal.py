import numpy as np
import yaml
from pymoab import core, types, rng, topo_util
import sys
import os
from utils import pymoab_utils as utpy

parent_dir = os.path.dirname(os.path.abspath(__file__))
parent_parent_dir = os.path.dirname(parent_dir)
input_dir = os.path.join(parent_parent_dir, 'input')
flying_dir = os.path.join(parent_parent_dir, 'flying')
# utils_dir = os.path.join(parent_parent_dir, 'utils')
# mono_dir = os.path.join(flying_dir, 'monofasico')
# bif_dir = os.path.join(flying_dir, 'bifasico')

class DefPrimal:
    def __init__(self):
        os.chdir(input_dir)
        with open("input_mesh_generator.yaml", 'r') as stream:
            data_loaded_mesh = yaml.load(stream)
            # data_loaded = yaml.load(stream, Loader=yaml.FullLoader)
            # data_loaded = yaml.full_load(stream)

        with open("input_wells.yaml", 'r') as stream:
            data_loaded_wells = yaml.load(stream)
            # data_loaded = yaml.load(stream, Loader=yaml.FullLoader)
            # data_loaded = yaml.full_load(stream)

        self.data_loaded = data_loaded_mesh
        file_name = data_loaded_mesh['file_name']
        self.file_name = file_name
        # self.Lz = float(data_loaded_mesh['mesh_size'][2])
        name_mesh = file_name + '_wells.h5m'
        self.mb = core.Core()
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        mesh_file = name_mesh
        os.chdir(flying_dir)
        self.mb.load_file(mesh_file)
        names_tags = np.load('tags_wells.npy')
        self.tags = utpy.get_all_tags_2(self.mb, names_tags)
        self.all_nodes, self.all_edges, self.all_faces, self.all_volumes = utpy.get_all_entities(self.mb)

        import pdb; pdb.set_trace()
