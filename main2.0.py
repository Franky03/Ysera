from ysera2 import Nodes, Edges
#import pymol
import time
import os


def run_software(name_=False, file=None, hydrogenate=False):
    start = time.time()

    if file is None:
        raise Exception("Load File")
    """
    if not hydrogenate:
        pymol.cmd.load(file, 'myprotein')
        pymol.cmd.h_add()
        pymol.cmd.save('./temp/input_file.pdb')
        time.sleep(2)
    """
    edges = Edges(name_, './temp/input_file.pdb')
    nodes= Nodes(name_, './temp/input_file.pdb')
    edges.to_file()
    nodes.to_file()

    finish = (time.time() - start)

    print(f"---{finish} seconds ---")

    with open("time_log.txt", "a") as file:
        file.write(f"---{finish} seconds ---\n")

    os.remove('./temp/input_file.pdb')


filename = 'file_30.pdb'
path = f'./temp/{filename}'
run_software('file_30', path, hydrogenate=False)
