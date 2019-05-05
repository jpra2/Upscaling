from definir_pocos.def1 import def_wells

class Wells:
    def __init__(self):
        pass


def run():
    wells = def_wells()
    wells.create_tags()
    wells.create_wells()
    wells.finish()


    print('saiu def_wells')
