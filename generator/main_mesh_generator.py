from generator.mesh_generator import MeshGenerator as mm

class Mesh:
    def __init__(self, mb, mtu, all_volumes, all_faces, all_edges, all_nodes, tags):
        self.mb = mb
        self.mtu = mtu
        self.all_volumes = all_volumes
        self.all_faces = all_faces
        self.all_edges = all_edges
        self.all_nodes = all_nodes
        self.tags = tags

def generate_mesh():
    mesh = mm()
    mesh.create_fine_vertices()
    mesh.create_tags()
    mesh.create_fine_blocks()
    mesh.get_all_entities()
    mesh.set_gids()
    mesh.set_centroids()
    mesh.set_area()
    mesh.finish()
    msh = Mesh(mesh.mb, mesh.mtu, mesh.all_volumes, mesh.all_faces, mesh.all_edges, mesh.verts, mesh.tags)
    return msh
