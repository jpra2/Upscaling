from generator.main_mesh_generator import generate_mesh as gm
from definir_pocos.main_def1 import run as run_wells
from primal.main_primal import gen_primal

gm()
run_wells()
gen_primal()
