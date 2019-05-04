import numpy as np
import yaml
from pymoab import core, types, rng, topo_util
import sys


# parent_dir = os.path.dirname(os.path.abspath(__file__))
# parent_parent_dir = os.path.dirname(parent_dir)
# input_dir = os.path.join(parent_parent_dir, 'input')
# flying_dir = os.path.join(parent_parent_dir, 'flying')
# utils_dir = os.path.join(parent_parent_dir, 'utils')
# mono_dir = os.path.join(flying_dir, 'monofasico')
# bif_dir = os.path.join(flying_dir, 'bifasico')

class MeshGenerator:

    def __init__(self):
        self.__verif = True
        with open("input_mesh_generator.yaml", 'r') as stream:
            data_loaded = yaml.load(stream)
            # data_loaded = yaml.load(stream, Loader=yaml.FullLoader)
            # data_loaded = yaml.full_load(stream)

        self.data_loaded = data_loaded
        self.tags = dict()
        self.mb = core.Core()
        self.mtu = topo_util.MeshTopoUtil(self.mb)
        self.nblocks = data_loaded['nblocks']
        self.mesh_size = data_loaded['mesh_size']
        self.blocksize = self.mesh_size/self.nblocks
        hx = self.blocksize[0]
        hy = self.blocksize[1]
        hz = self.blocksize[2]
        self.Area = np.array([hy*hz, hx*hz, hy*hx])

    def create_fine_vertices(self):

        x = np.linspace(0, self.mesh_size[0], self.nblocks[0]+1)
        y = np.linspace(0, self.mesh_size[1], self.nblocks[1]+1)
        z = np.linspace(0, self.mesh_size[2], self.nblocks[2]+1)
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        xx = xx.flatten()
        yy = yy.flatten()
        zz = zz.flatten()
        coords_verts = np.array([xx, yy, zz]).T
        # for i in coords_verts:
        #     print(i)
        #     import pdb; pdb.set_trace()
        self.verts = self.mb.create_vertices(coords_verts.flatten())
        self.coords_verts = coords_verts
        self.ind_verts = np.arange(1,len(self.verts)+1).astype(np.uint64)

    def create_fine_blocks(self):
        # hexas = []
        # for i in range(self.nblocks[0]):
        #     for j in range(self.nblocks[1]):
        #         for k in range(self.nblocks[2]):
        #             hexa = self.create_hexa(i, j, k)
        #             hexas.append(hexa)

        hexas = np.array([self.create_hexa(i, j, k) for i in range(self.nblocks[0])
                 for j in range(self.nblocks[1])
                 for k in range(self.nblocks[2])])
        [self.mb.create_element(types.MBHEX, hexa) for hexa in hexas]

    def create_hexa(self, i, j, k):
        indices =np.array([(k) + (j)*(self.nblocks[2]+1) + (i)*(self.nblocks[2]+1)*(self.nblocks[1]+1),
                           (k+1) + (j)*(self.nblocks[2]+1) + (i)*(self.nblocks[2]+1)*(self.nblocks[1]+1),
                           (k+1) + (j+1)*(self.nblocks[2]+1) + (i)*(self.nblocks[2]+1)*(self.nblocks[1]+1),
                           (k) + (j+1)*(self.nblocks[2]+1) + (i)*(self.nblocks[2]+1)*(self.nblocks[1]+1),
                           (k) + (j)*(self.nblocks[2]+1) + (i+1)*(self.nblocks[2]+1)*(self.nblocks[1]+1),
                           (k+1) + (j)*(self.nblocks[2]+1) + (i+1)*(self.nblocks[2]+1)*(self.nblocks[1]+1),
                           (k+1) + (j+1)*(self.nblocks[2]+1) + (i+1)*(self.nblocks[2]+1)*(self.nblocks[1]+1),
                           (k) + (j+1)*(self.nblocks[2]+1) + (i+1)*(self.nblocks[2]+1)*(self.nblocks[1]+1)
                           ]).astype(np.uint64)

        hexa = self.ind_verts[indices]
        return hexa

    def create_tags(self):
        self.tags['GIDV'] = self.mb.tag_get_handle('GIDV', 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
        self.tags['GIDF'] = self.mb.tag_get_handle('GIDF', 1, types.MB_TYPE_INTEGER, types.MB_TAG_SPARSE, True)
        self.tags['CENT'] = self.mb.tag_get_handle('CENT', 3, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.tags['AREA'] = self.mb.tag_get_handle('AREA', 1, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.tags['NORMAL'] = self.mb.tag_get_handle('NORMAL', 3, types.MB_TYPE_DOUBLE, types.MB_TAG_SPARSE, True)
        self.tags['BOUND_FACES'] = self.mb.tag_get_handle('BOUND_FACES', 1, types.MB_TYPE_HANDLE, types.MB_TAG_MESH, True)

    def get_all_entities(self):
        self.all_volumes = self.mb.get_entities_by_dimension(0, 3)
        self.mtu.construct_aentities(self.verts)
        self.all_faces = self.mb.get_entities_by_dimension(0, 2)
        self.all_edges = self.mb.get_entities_by_dimension(0, 1)

    def finish(self):
        self.__verif = False

    def set_area(self):
        gids_volumes = self.mb.tag_get_data(self.tags['GIDV'], self.all_volumes, flat=True)
        map_volumes = dict(zip(self.all_volumes, gids_volumes))
        areas = np.zeros(len(self.all_faces))
        normais = np.zeros((len(self.all_faces), 3))
        Adjs = [self.mb.get_adjacencies(f, 3) for f in self.all_faces]
        vertss = [self.mtu.get_bridge_adjacencies(f, 2, 0) for f in self.all_faces]
        lim = 1e-9
        bound_faces = self.mb.create_meshset()
        bfs = []
        for i, f in enumerate(self.all_faces):
            elems = Adjs[i]
            if len(elems) < 2:
                bfs.append(f)
            verts = vertss[i]
            coords = self.mb.get_coords(verts).reshape([len(verts), 3])
            mins = coords.min(axis=0)
            maxs = coords.max(axis=0)
            dx, dy, dz = np.absolute(maxs - mins)
            if dx < lim:
                dx = 1.0
            if dy < lim:
                dy = 1.0
            if dz < lim:
                dz = 1.0
            area = dx*dy*dz
            areas[i] = area
            if len(elems) > 1:

                id0 = map_volumes[elems[0]]
                id1 = map_volumes[elems[1]]
                cent0 = self.all_centroids[id0]
                cent1 = self.all_centroids[id1]
                direction = cent1 - cent0
                norma = np.linalg.norm(direction)
                uni = np.absolute(direction/norma)
                normais[i] = uni
            else:
                p0 = coords[0]
                p1 = coords[1]
                p2 = coords[2]
                direction = np.cross(p0-p1, p0-p2)
                norma = np.linalg.norm(direction)
                uni = np.absolute(direction/norma)
                normais[i] = uni

        self.mb.tag_set_data(self.tags['AREA'], self.all_faces, areas)
        self.mb.tag_set_data(self.tags['NORMAL'], self.all_faces, normais)
        bfs = rng.Range(bfs)
        self.mb.add_entities(bound_faces, bfs)
        self.mb.tag_set_data(self.tags['BOUND_FACES'], 0, bound_faces)
        self.bound_faces = bfs

    def set_centroids(self):
        all_centroids = np.array([self.mtu.get_average_position([v]) for v in self.all_volumes])
        for cent, v in zip(all_centroids, self.all_volumes):
            self.mb.tag_set_data(self.tags['CENT'], v, cent)
        self.all_centroids = all_centroids

    def set_gids(self):
        self.mb.tag_set_data(self.tags['GIDV'], self.all_volumes, np.arange(0, len(self.all_volumes)))
        self.mb.tag_set_data(self.tags['GIDF'], self.all_faces, np.arange(0, len(self.all_faces)))

    def nblocks():
        doc = "The nblocks property."
        def fget(self):
            return self._nblocks
        def fset(self, value):
            if self.__verif:
                try:
                    self._nblocks = np.array([int(v) for v in value])
                except:
                    print('Dados invalidos para nblocks. Deve ser uma lista de inteiros')
                    sys.exit(0)
            else:
                print('nao eh possivel alterar nblocks')

        return locals()
    nblocks = property(**nblocks())

    def mesh_size():
        doc = "The meshsize property."
        def fget(self):
            return self._mesh_size
        def fset(self, value):
            if self.__verif:
                try:
                    self._mesh_size = np.array([float(v) for v in value])
                except:
                    print('Dados invalidos para meshsize. Deve ser uma lista de floats')
                    sys.exit(0)
            else:
                print('nao eh possivel alterar meshsize')

        return locals()
    meshsize = property(**mesh_size())

    def blocksize():
        doc = "The blocksize property."
        def fget(self):
            return self._blocksize
        def fset(self, value):
            if self.__verif:
                self._blocksize = value
            else:
                print('Nao pode alterar blocksize')
        def fdel(self):
            del self._blocksize
        return locals()
    blocksize = property(**blocksize())
