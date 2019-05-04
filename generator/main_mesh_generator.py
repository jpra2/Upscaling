from generator.mesh_generator import MeshGenerator as mm
import os
import numpy as np

parent_dir = os.path.dirname(os.path.abspath(__file__))
parent_parent_dir = os.path.dirname(parent_dir)
input_dir = os.path.join(parent_parent_dir, 'input')
flying_dir = os.path.join(parent_parent_dir, 'flying')

class Mesh:
    def __init__(self, mb, mtu, all_volumes, all_faces, all_edges, all_nodes, tags, data_loaded):
        self.mb = mb
        self.mtu = mtu
        self.all_volumes = all_volumes
        self.all_faces = all_faces
        self.all_edges = all_edges
        self.all_nodes = all_nodes
        self.tags = tags
        self.data_loaded = data_loaded

def generate_mesh():
    os.chdir(input_dir)
    mesh = mm()
    mesh.create_fine_vertices()
    mesh.create_tags()
    mesh.create_fine_blocks()
    mesh.get_all_entities()
    mesh.set_gids()
    mesh.set_centroids()
    mesh.set_area()
    mesh.finish()
    msh = Mesh(mesh.mb, mesh.mtu, mesh.all_volumes, mesh.all_faces, mesh.all_edges, mesh.verts, mesh.tags, mesh.data_loaded)
    os.chdir(flying_dir)
    with open('__init__.py', 'w') as ff:
        pass
    file_name = mesh.data_loaded['file_name']
    ext_h5m = file_name + '_mesh1.h5m'
    tags_mesh_generator = np.array(list(mesh.tags.keys()))
    np.save('tags_mesh_generator', tags_mesh_generator)
    mesh.mb.write_file(ext_h5m)
