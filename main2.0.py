from ysera2 import Nodes, Edges
import pymol
import time

def run_software(name_= False, file= None, hydrogenate = False):
    start= time.time()

    if file is None:
        raise Exception("Load File")
    
    if not hydrogenate:
        pymol.cmd.load(file, 'myprotein')
        pymol.cmd.h_add()
        pymol.cmd.save('./temp/input_file.pdb')
        time.sleep(2)

    edges= Edges(name_, './temp/input_file.pdb', multiple= True)
    edges.to_file()
    

    print(f"---{(time.time() - start)} seconds ---")


run_software('3og7', './temp/3og7.pdb', hydrogenate= True)

