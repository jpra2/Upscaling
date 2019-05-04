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

class def_wells:
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

        self.data_loaded = data_loaded_wells
        file_name = data_loaded_mesh['file_name']
        self.file_name = file_name
        self.Lz = float(data_loaded_mesh['mesh_size'][2])
        name_mesh = file_name + '_mesh1.h5m'
        self.mb = core.Core()
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        mesh_file = name_mesh
        os.chdir(flying_dir)
        self.mb.load_file(mesh_file)

        names_tags = np.load('tags_mesh_generator.npy')
        self.tags = utpy.get_all_tags_2(self.mb, names_tags)
        self.gama = data_loaded_wells['dados_monofasico']['gama']
        self.all_nodes, self.all_edges, self.all_faces, self.all_volumes = utpy.get_all_entities(self.mb)

    def create_tags(self):
        self.tags['P'] = self.mb.tag_get_handle('P', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.tags['Q'] = self.mb.tag_get_handle('Q', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.tags['WELLS_INJECTOR'] = self.mb.tag_get_handle('WELLS_INJECTOR', 1, types.MB_TYPE_HANDLE, types.MB_TAG_MESH, True)
        self.tags['WELLS_PRODUCER'] = self.mb.tag_get_handle('WELLS_PRODUCER', 1, types.MB_TYPE_HANDLE, types.MB_TAG_MESH, True)

    def create_wells(self):
        self.wells_injector = self.mb.create_meshset()
        self.wells_producer = self.mb.create_meshset()
        self.mb.tag_set_data(self.tags['WELLS_PRODUCER'], 0, self.wells_producer)
        self.mb.tag_set_data(self.tags['WELLS_INJECTOR'], 0, self.wells_injector)

        all_centroids = self.mb.tag_get_data(self.tags['CENT'], self.all_volumes)
        Lz = self.Lz

        data_loaded = self.data_loaded
        gravity = data_loaded['gravity']
        all_volumes = self.all_volumes
        gama = self.gama

        bvd = []
        bvn = []
        mb = self.mb
        volumes_d = []
        volumes_n = []

        inds_wells = []
        for well in data_loaded['Wells_structured']:
            w = data_loaded['Wells_structured'][well]
            if w['type_region'] == 'box':
                box_volumes = np.array([np.array(w['region1']), np.array(w['region2'])])
                inds0 = np.where(all_centroids[:,0] > box_volumes[0,0])[0]
                inds1 = np.where(all_centroids[:,1] > box_volumes[0,1])[0]
                inds2 = np.where(all_centroids[:,2] > box_volumes[0,2])[0]
                c1 = set(inds0) & set(inds1) & set(inds2)
                inds0 = np.where(all_centroids[:,0] < box_volumes[1,0])[0]
                inds1 = np.where(all_centroids[:,1] < box_volumes[1,1])[0]
                inds2 = np.where(all_centroids[:,2] < box_volumes[1,2])[0]
                c2 = set(inds0) & set(inds1) & set(inds2)
                inds_vols = list(c1 & c2)
                inds_wells += inds_vols
                volumes = rng.Range(np.array(all_volumes)[inds_vols])
            else:
                raise NameError("Defina o tipo de regiao em type_region: 'box'")

            value = float(w['value'])

            if w['type_prescription'] == 'dirichlet':
                bvd.append(box_volumes)
                volumes_d += list(volumes)
                if gravity == False:
                    pressao = np.repeat(value, len(volumes))

                else:
                    z_elems_d = -1*mb.tag_get_data(self.tags['CENT'], volumes)[:,2]
                    delta_z = z_elems_d + Lz
                    pressao = gama*(delta_z) + value

                mb.tag_set_data(self.tags['P'], volumes, pressao)

            elif  w['type_prescription'] == 'neumann':
                volumes_n += list(volumes)
                bvn.append(box_volumes)
                value = value/len(volumes)
                if w['type_well'] == 'producer':
                    value = -value
                mb.tag_set_data(self.tags['Q'], volumes, np.repeat(value, len(volumes)))

            else:
                raise NameError("type_prescription == 'neumann' or 'dirichlet'")

            if w['type_well'] == 'injector':
                mb.add_entities(self.wells_injector, volumes)
            else:
                mb.add_entities(self.wells_producer, volumes)

    def finish(self):
        os.chdir(flying_dir)
        names_tags = np.array(list(self.tags.keys()))
        np.save('tags_wells', names_tags)
        ext_h5m = self.file_name + '_wells.h5m'
        self.mb.write_file(ext_h5m)
