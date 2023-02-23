from ysera2 import Nodes, Edges
import pymol
import time
import os

class MakeFiles:
    def __init__(self,name=False, file=None, hydrogenate=False):

        if file is None:
            raise Exception("Load File")
        
        if not hydrogenate:
            pymol.cmd.load(file, 'myprotein')
            pymol.cmd.h_add()
            pymol.cmd.save('./temp/input_file.pdb')
            
        
        self.edges = Edges(name, './temp/input_file.pdb')
        self.nodes = Nodes(name, './temp/input_file.pdb')

    def get_node_degrees(self):
        # getting the number of ligands in a residue 
        for node in self.nodes.nodes_id:
            degree = 0
            degree += self.edges.nodes_id1.count(node)
            degree += self.edges.nodes_id2.count(node)

            self.nodes.degrees.append(degree)

    def run_software(self):

        start = time.time()
        self.edges.to_file()
        self.nodes.search_nodes()
        self.get_node_degrees()
        self.nodes.to_file()

        finish = (time.time() - start)

        print(f"---{finish} seconds ---")


filename = '3og7.pdb'
path = f'./temp/{filename}'

run = MakeFiles('3og7', path, hydrogenate=False)
run.run_software()
