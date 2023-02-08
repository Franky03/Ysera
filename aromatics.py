from Bio.PDB import PDBParser, NeighborSearch


class FormatAromatics:
    def __init__(self):

        self.aromatic_pos = []
        self.aromatic_points = []
        self.aromatic_array = {}
        self.aromatic_normals = {}
        self.invalids = []

    def calcula_array(self):
        lig_032 =[]
        for residue in self.structure.get_residues():
            if str(residue.resname) == "032":
                lig_032.append(residue)

        for residue in lig_032:
            for atom in residue:
                for neighbor_pair in self.ns.search(atom.coord, 6.5, level='R'):
                    print(neighbor_pair)


if __name__ == '__main__':
    inst = FormatAromatics()
    inst.calcula_array()






