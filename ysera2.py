from Bio.PDB import PDBParser, Selection, NeighborSearch
from Bio.PDB.vectors import calc_angle
import numpy as np
import pandas as pd
import time
import mdtraj as md
from aromatics import AromaticsFormat


class Nodes:
    def __init__(self, name_=None, file_=False):
        self.name = name_
        self.parser = PDBParser(PERMISSIVE=1)
        self.structure = self.parser.get_structure(name_, file_)
        self.ns = NeighborSearch(list(self.structure.get_atoms()))
        # mdtraj
        self.pdb = md.load_pdb(file_)
        self.dssp_md = md.compute_dssp(self.pdb, simplified=False)

        self.model = self.structure[0]
        # self.dssp = DSSP(self.model, './Codes/3og7.pdb', file_type='PDB', dssp='dssp')
        self.all_dssps = []
        self.nodes_id, self.chains, self.positions, self.residues = [], [], [], []
        # Degrees não está com uma precisão boa, é preciso saber como pode ser feito igual o do RINGs
        self.degrees = []
        self.cut_dist = 8.0  # definir um limite de distância de corte (escolhido com base na literatura)
        # B-Factor, coords and filenames
        self.bfactors, self.coords, self.pdb_filename = [], [], []
        self.rapdfs = []

    def dssp_test(self):
        self.search_nodes()
        for residue in self.dssp:
            print(residue)

    def search_nodes(self):
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if str(residue.resname) != 'HOH':  # ignore solvent
                        # Node_ID, Chain, Position and Residue
                        self.nodes_id.append(f"{chain.id}:{str(residue.id[1])}:_:{str(residue.resname)}")
                        self.chains.append(str(chain.id))
                        self.positions.append(residue.id[1])
                        if str(residue.resname) == '032':
                            # Se o resíduo for o 032 então ele não terá um Bfactor-CA nem coordenadas
                            self.bfactors.append('NaN')
                            self.coords.append(np.array(['NaN', 'NaN', 'NaN']))
                        self.residues.append(str(residue.resname))

                        # Degree -> número de residuos diretamente conectados

                        degree = 0

                        if 'CA' in residue:
                            # R ->  list of residues list of modules
                            for residue_2 in Selection.unfold_entities(model, 'R'):
                                if (residue_2.get_id()[1] != residue.get_id()[1]):
                                    # CA -> Carbono Alfa (mais preciso que outros átomos)
                                    # Calcular a distancia euclidiana das coordenadas desse carbono alfa
                                    if 'CA' in residue_2:
                                        distance = np.linalg.norm(residue["CA"].coord - residue_2["CA"].coord)
                                        if (distance < self.cut_dist):
                                            degree += 1

                        self.degrees.append(degree)

                        # Bfactor_CA
                        b_factor = 0
                        count = 0
                        for atom in residue:
                            if (atom.get_name() == 'CA'):
                                b_factor += atom.get_bfactor()
                                count += 1

                                coords = atom.get_coord()
                                self.coords.append(coords)

                        if (count != 0):
                            bf_average = b_factor / count
                            self.bfactors.append(f"{bf_average:.3f}")

                        # pdb filenames

                        self.pdb_filename.append(f"input_file.cif#{str(residue.id[1])}.{str(chain.id)}")

                        # Rapdf
                        sum_of_dis = 0
                        rapdf_count = 0
                        for residue_2 in Selection.unfold_entities(model, 'R'):
                            if 'CA' in residue and 'CA' in residue_2 and residue_2.id[1] != residue.id[1]:
                                sum_of_dis += np.linalg.norm(residue["CA"].coord - residue_2["CA"].coord)
                                rapdf_count += 1
                        if rapdf_count != 0:
                            rapdf = sum_of_dis / rapdf_count
                            self.rapdfs.append(rapdf)
                        else:
                            self.rapdfs.append(0.0)

        # Pegando o DSSP
        for i in range(len(self.dssp_md[0])):
            if i > 0:
                if self.dssp_md[0][i] == 'NA':
                    pass
                else:
                    self.all_dssps.append(self.dssp_md[0][i] if self.dssp_md[0][i] != 'C' else '\t')
            else:

                self.all_dssps.append(self.dssp_md[0][i] if self.dssp_md[0][i] not in ['C', 'NA'] else '\t')

    def print_output(self):
        self.search_nodes()
        for n in range(len(self.nodes_id)):
            try:
                print(
                    f"{self.nodes_id[n]}\t{self.chains[n]}\t\t{self.positions[n]}\t\t{self.residues[n]}\t{self.all_dssps[n]}\t" +
                    f"{self.degrees[n]}\t{self.bfactors[n]:.3f}\t{self.coords[n][0]:.3f}\t{self.coords[n][1]:.3f}\t{self.coords[n][2]:.3f}\t" +
                    f"{self.pdb_filename[n]}\t{self.rapdfs[n]}"
                    )
            except Exception as e:
                print(
                    f"{self.nodes_id[n]}\t{self.chains[n]}\t\t{self.positions[n]}\t\t{self.residues[n]}\t{self.all_dssps[n]}\t" +
                    f"{self.degrees[n]}\t{self.bfactors[n]}\t{self.coords[n][0]}\t{self.coords[n][1]}\t{self.coords[n][2]}\t" +
                    f"{self.pdb_filename[n]}\t{self.rapdfs[n]}"
                    )

    def to_file(self):
        self.search_nodes()

        colunas = ["NodeId", "Chain", "Position", "Residue", "Dssp", "Degree", "Bfactor_CA", "x", "y", "z",
                   "pdbFileName", "Model"]

        x, y, z = [], [], []
        for coord in self.coords:
            if coord[0] != 'NaN':
                x.append(f"{coord[0]:.3f}")
                y.append(f"{coord[1]:.3f}")
                z.append(f"{coord[2]:.3f}")
            else:
                x.append(coord[0])
                y.append(coord[1])
                z.append(coord[2])

        data = pd.DataFrame(list(zip(self.nodes_id, self.chains, self.positions, self.residues,
                                     self.all_dssps, self.degrees, self.bfactors, x, y, z, self.pdb_filename,
                                     self.models
                                     )), columns=colunas)

        data.to_csv(f'./{self.name}_nodes.csv', sep='\t', index=False)


class Edges(Nodes):
    def __init__(self, name, file_pdb, multiple=True):
        af = AromaticsFormat('3og7.pdb')
        self.aromatic_array, self.aromatic_normals, self.invalids = af.get_data()
        Nodes.__init__(self, name_=name, file_=file_pdb)
        self.lighbdonor = {'ARG': ['NE', 'NH1', 'NH2'],
                           'ASN': ['ND2'],
                           'HIS': ['NE2', 'ND1'],
                           'SER': ['OG'],
                           'TYR': ['OH'],
                           'CYS': ['SG'],
                           'THR': ['OG1'],
                           'GLN': ['NE2'],
                           'LYS': ['NZ'],
                           'TRP': ['NE1']
                           }
        self.lighbac = {'ASN': ['OD1'],
                        'GLN': ['OE1'],
                        'MET': ['SD'],
                        'ASP': ['OD1', 'OD2'],
                        'GLU': ['OE1', 'OE2'],
                        'SER': ['OG'],
                        'THR': ['OG1'],
                        'HIS': ['ND1'],
                        'TYR': ['OH']
                        }
        self.ligvdw = ['C', 'CB', 'CG1', 'CG2', 'CD1', 'CD2', 'CE']
        self.ligpipi = ['HIS', 'TYR', 'TRP', 'PHE']
        self.nodes_id1, self.nodes_id2, self.bonds = [], [], []
        self.distances, self.donors, self.angles = [], [], []
        self.atom1, self.atom2 = [], []
        self.bonds_check, self.energies, self.orientation = [], [], []
        self.analyzed_pairs = set()
        self.multiple = multiple
        self.ligands = {'hb': 0, 'vdw': 0, 'ionic': 0, 'sbond': 0, 'pi_stacking': 0}
        self.exclusions = []

    def Iac(self):

        # Eles não calculam mais essas ligações na nova versão 

        lig_032 = []
        for residue in self.structure.get_residues():
            if str(residue.resname) == "032":
                lig_032.append(residue)

        for residue in lig_032:
            for atom in residue:
                for neighbor_pair in self.ns.search(atom.coord, 6.5, level='R'):
                    for atom2 in neighbor_pair:
                        if atom2.get_name() == 'CA':
                            distance = np.linalg.norm(atom.coord - atom2.coord)
                            # verificando se o átomo vizinho é de outro resíduo
                            if neighbor_pair != residue:
                                print(residue.resname, neighbor_pair.resname, neighbor_pair.id[1], distance)

    def add_bond(self, config, ligand):
        """
        config's params: chain, residue, neighbor residue,
                        bond, distance, angle, energie, 
                        atom1, atom2, donor.
        ligand: hb, vdw, ionic... 
        """

        self.nodes_id1.append(f"{config[0].id}:{str(config[1].id[1])}:_:{str(config[1].resname)}")
        self.nodes_id2.append(f"{config[2].get_parent().id}:{str(config[2].id[1])}:_:{str(config[2].resname)}")
        self.bonds.append(config[3])
        self.distances.append(config[4])
        self.angles.append(config[5])
        if ligand == 'hb':
            if float(config[4]) <= 1.5:
                self.energies.append(config[6][0])
            elif float(config[4]) >= 2.2:
                self.energies.append(config[6][1])
            else:
                self.energies.append(config[6][1])
        else:
            self.energies.append(config[6])

        self.atom1.append(config[7])
        self.atom2.append(config[8])
        self.donors.append(config[9])
        self.orientation.append(config[10])
        self.ligands[ligand] += 1

    def _hydrogen_bond(self, chain, residue, atom):
        chain1 = ''
        chain2 = ''

        global n_or_o_donor
        global h_donor

        atom_name = atom.get_name()

        if atom.fullname[1] in ['N', 'O'] or (atom_name == 'SG' and residue.resname == 'CYS'):
            neighbors = self.ns.search(atom.coord, 5.5)
            for neighbor in neighbors:

                neig_name = neighbor.get_name()
                neig_res = neighbor.get_parent()
                if neig_res.resname in ['HOH', '032']:
                    continue
                if neig_res.id[1] == residue.id[1] or neig_name[0] == atom_name[0]:
                    continue

                pair = (int(residue.id[1]), int(neig_res.id[1]))

                if pair in self.analyzed_pairs:
                    continue
                else:
                    self.analyzed_pairs.add((int(neig_res.id[1]), int(residue.id[1])))

                if neighbor.fullname[1] in ['N', 'O'] or (neighbor.get_name() == 'SG' and neig_res.resname == 'CYS'):
                    distance = np.linalg.norm(atom.coord - neighbor.coord)

                    # Verificando quem é doador
                    if (atom_name[0] == 'N' or (atom_name in ['OG', 'OH', 'OG1', 'SG'] and residue.resname in list(
                            self.lighbdonor.keys()))) and (neig_name[0] == 'O' or (
                            neig_name in ['SD', 'ND1'] and neig_res.resname in list(self.lighbac.keys()))):
                        # Aqui o doador vai ser o atomo principal
                        n_or_o_donor = atom
                        try:
                            h_list = [a for a in residue if a.element == 'H']
                            h_distances = {}
                            for h_atom in h_list:
                                h_dist = np.linalg.norm(atom.coord - h_atom.coord)
                                h_distances[h_dist] = h_atom
                            min_h = min(list(h_distances.keys()))
                            h_donor = h_distances[min_h]
                        except:
                            raise Exception("Hydrogens not found, hydrogenate the pdb file first!")


                    elif (neig_name[0] == 'N' or (neig_name in ['OG', 'OH', 'OG1', 'SG'] and neig_res.resname in list(
                            self.lighbdonor.keys()))) and (atom_name[0] == 'O' or (
                            atom_name in ['SD', 'ND1'] and residue.resname in list(self.lighbac.keys()))):
                        # Aqui o doador vai ser o atomo vizinho
                        n_or_o_donor = neighbor
                        try:

                            h_list = [a for a in neig_res if a.element == 'H']

                            h_distances = {}
                            for h_atom in h_list:
                                h_dist = np.linalg.norm(neighbor.coord - h_atom.coord)
                                h_distances[h_dist] = h_atom
                            min_h = min(list(h_distances.keys()))
                            h_donor = h_distances[min_h]
                        except:
                            raise Exception("Hydrogens not found, hydrogenate the pdb file first!")

                    terceiro_vetor = h_donor.get_vector()
                    neighbor_vector = neighbor.get_vector()
                    a_vector = atom.get_vector()

                    angle = 0.0
                    if n_or_o_donor == atom:
                        angle = np.degrees(calc_angle(terceiro_vetor, a_vector, neighbor_vector))
                    else:
                        angle = np.degrees(calc_angle(terceiro_vetor, neighbor_vector, a_vector))

                    if 2.5 < distance <= 3.5 and angle <= 63.0:

                        # Verificando se é a cadeia principal ou lateral
                        if atom.name in ["N", "O"]:
                            chain1 = 'MC'
                        else:
                            chain1 = 'SC'

                        if neighbor.name in ['N', 'O']:
                            chain2 = 'MC'
                        else:
                            chain2 = 'SC'

                        if self.multiple:
                            self.bonds_check.append((f"{chain.id}:{str(residue.id[1])}:_:{str(residue.resname)}",
                                                     f"{neig_res.get_parent().id}:{str(neig_res.id[1])}:_:{str(neig_res.resname)}"))

                        self.add_bond([
                            chain, residue, neig_res,
                            f"HBOND:{chain1}_{chain2}",
                            f"{distance:.3f}",
                            f"{angle:.3f}",
                            (f"{115.000:.3f}", f"{17.000:.3f}", f"{40.000:.3f}"),
                            atom_name,
                            neig_name,
                            f"{chain.id}:{str(n_or_o_donor.get_parent().id[1])}:_:{str(n_or_o_donor.get_parent().resname)}",
                            "NaN"
                        ], 'hb')

    def _vanderwaals(self, chain, residue, atom):
        chain1 = ''
        chain2 = ''

        vdw_radii = {'C': 1.77, 'S': 1.89, 'N': 1.8, 'O': 1.4}
        is_vdw = False

        atom_name = atom.get_name()
        if atom.fullname[1] in ['C', 'S', 'O', 'N']:
            neighbors = self.ns.search(atom.coord, 3.9)
            for neighbor in neighbors:
                is_vdw = False

                neig_name = neighbor.get_name()
                neig_res = neighbor.get_parent()
                distance = np.linalg.norm(atom.coord - neighbor.coord)
                if neig_res.id[1] == residue.id[1] or neig_name in ["CA", "CH2"] or atom_name in ["CA", "CH2"] or (
                        atom_name == 'C' and neig_name == 'C'):
                    continue

                if neig_res.resname in ['HOH', '032']:
                    continue

                pair = (int(residue.id[1]), int(neig_res.id[1]))

                if pair in self.analyzed_pairs:
                    continue
                else:
                    self.analyzed_pairs.add((int(neig_res.id[1]), int(residue.id[1])))

                if neighbor.fullname[1] in ['C', 'S', 'O', 'N']:

                    if atom.name in ["C", "S"]:
                        chain1 = 'MC'
                    else:
                        chain1 = 'SC'

                    if neighbor.name in ['C', 'S']:
                        chain2 = 'MC'
                    else:
                        chain2 = 'SC'

                    if (atom.fullname[1] == "C" and neighbor.fullname[1] == "C") or (
                            atom.fullname[1] == "C" and neighbor.fullname[1] == "S") or (
                            atom.fullname[1] == "S" and neighbor.fullname[1] == "C"):
                        is_vdw = True

                    elif (atom_name[0] == "N" or atom_name[0] == "O") and neig_name[0] == "C":
                        if (residue.resname == 'GLN' and (atom_name == "OE1" or atom_name == "NE2")) or (
                                residue.resname == 'ASN' and (atom_name == "OD1" or atom_name == "ND2")):
                            is_vdw = True

                    elif (neig_name[0] == "N" or neig_name[0] == "O") and atom_name[0] == "C":
                        if (neig_res.resname == 'GLN' and (neighbor.name == "OE1" or neighbor.name == "NE2")) or (
                                neig_res.resname == 'ASN' and (neighbor.name == "OD1" or neighbor.name == "ND2")):
                            is_vdw = True

                if is_vdw:
                    check_dist = distance - vdw_radii[atom.name[0]] - vdw_radii[neighbor.name[0]]

                    if check_dist <= 0.5:
                        if self.multiple:
                            self.bonds_check.append((f"{chain.id}:{str(residue.id[1])}:_:{str(residue.resname)}",
                                                     f"{neig_res.get_parent().id}:{str(neig_res.id[1])}:_:{str(neig_res.resname)}"))

                        self.add_bond([
                            chain, residue, neig_res,
                            f"VDW:{chain1}_{chain2}",
                            f"{distance:.3f}",
                            "NaN",
                            f"{6.000:.3f}",
                            atom_name,
                            neig_name,
                            "NaN",
                            "NaN"
                        ], 'vdw')

    def _dissulfide_bond(self, chain, residue, atom):
        # identificar as cadeias na ligação dissufide
        chain1 = ''
        chain2 = ''

        atom_name = atom.get_name()

        if atom_name[0] == 'S':
            neighbors = self.ns.search(atom.coord, 3.5)
            for neighbor in neighbors:
                neig_res = neighbor.get_parent()

                pair = (int(residue.id[1]), int(neig_res.id[1]))

                if neig_res.id[1] == residue.id[1]:
                    continue
                if pair in self.analyzed_pairs:
                    continue
                else:
                    self.analyzed_pairs.add((int(neig_res.id[1]), int(residue.id[1])))

                neig_name = neighbor.get_name()
                neig_res = neighbor.get_parent()
                distance = np.linalg.norm(atom.coord - neighbor.coord)
                if neig_name[0] == 'S' and distance <= 2.5:

                    if self.multiple:
                        self.bonds_check.append((f"{chain.id}:{str(residue.id[1])}:_:{str(residue.resname)}",
                                                 f"{neig_res.get_parent().id}:{str(neig_res.id[1])}:_:{str(neig_res.resname)}"))

                    self.add_bond([
                        chain, residue, neig_res,
                        f"SBOND:{chain1}_{chain2}",
                        f"{distance:.3f}",
                        "NaN",
                        f"{167.000:.3f}",
                        atom_name,
                        neig_name,
                        "NaN",
                        "NaN"
                    ], 'sbond')

    def _salt_bridge(self, chain, residue, atom):
        global ionic_donor
        chain1 = ''
        chain2 = ''

        atom_name = atom.get_name()
        if residue.resname in ['ARG', 'LYS', 'HIS', 'ASP', 'GLU']:
            analyzed_ionic = set()
            neighbors = self.ns.search(atom.coord, 8)
            for neighbor in neighbors:
                neig_res = neighbor.get_parent()
                neig_name = neighbor.get_name()
                if neig_res.id[1] == residue.id[1]:
                    continue
                if atom_name in ['CZ', 'NZ'] or neig_res in ['CZ', 'NZ']:
                    pair = (int(residue.id[1]), int(neig_res.id[1]))

                    if pair in self.analyzed_pairs or pair in analyzed_ionic:
                        continue
                    else:
                        self.analyzed_pairs.add((int(neig_res.id[1]), int(residue.id[1])))
                        analyzed_ionic.add(pair)

                    if neig_res.resname in ['ARG', 'LYS', 'HIS', 'ASP', 'GLU']:

                        if residue.resname in ['ARG', 'LYS', 'HIS'] and neig_res.resname in ['ASP', 'GLU']:
                            # aqui o átomo vai ser o doador 
                            ionic_donor = atom

                        elif neig_res.resname in ['ARG', 'LYS', 'HIS'] and residue.resname in ['ASP', 'GLU']:
                            # aqui o neighbor vai ser o doador
                            ionic_donor = neighbor

                        chain1 = 'MC' if len(atom_name) == 1 else 'SC'
                        chain2 = 'MC' if len(neig_name) == 1 else 'SC'

                        distance = np.linalg.norm(atom.coord - neighbor.coord)
                        # angle = np.degrees(calc_angle())

                        if distance <= 4.0:
                            if self.multiple:
                                self.bonds_check.append((f"{chain.id}:{str(residue.id[1])}:_:{str(residue.resname)}",
                                                         f"{neig_res.get_parent().id}:{str(neig_res.id[1])}:_:{str(neig_res.resname)}"))
                                # ionic_check.add(((residue.id[1], neig_res.id[1]),atom, neighbor, distance))

                            if atom_name in ['CZ', 'NZ']:

                                self.add_bond([
                                    chain, residue, neig_res,
                                    f"IONIC:{chain1}_{chain2}",
                                    f"{distance:.3f}",
                                    "NaN",
                                    f"{20.000:.3f}",
                                    atom_name,
                                    f"{neighbor.get_coord()[0]:.3f},{neighbor.get_coord()[1]:.3f},{neighbor.get_coord()[2]:.3f}",
                                    f"{chain.id}:{str(ionic_donor.get_parent().id[1])}:_:{str(ionic_donor.get_parent().resname)}",
                                    "NaN"
                                ], 'ionic')

                            elif atom_name not in ['CZ', 'NZ'] and neig_name in ['CZ', 'NZ']:

                                self.add_bond([
                                    chain, residue, neig_res,
                                    f"IONIC:{chain1}_{chain2}",
                                    f"{distance:.3f}",
                                    "NaN",
                                    f"{20.000:.3f}",
                                    f"{atom.get_coord()[0]:.3f},{atom.get_coord()[1]:.3f},{atom.get_coord()[2]:.3f}",
                                    neig_name,
                                    f"{chain.id}:{str(ionic_donor.get_parent().id[1])}:_:{str(ionic_donor.get_parent().resname)}",
                                    "NaN"
                                ], 'ionic')

    def _pi_stacking(self, chain, residue, atom):
        neighbors = self.ns.search(atom.coord, 7.2)
        amin = f'{chain.id} {residue.id[1]}'
        orient_type = ''
        for neighbor in neighbors:

            neig_res = neighbor.get_parent()
            neig_chain = neig_res.get_parent().id
            neig_amin = f'{neig_chain} {neig_res.id[1]}'

            if residue.get_resname() in ['TYR', 'PHE', 'TRP'] and neig_res.get_resname() in ['TYR', 'PHE', 'TRP']:
                if (amin not in self.invalids and neig_amin not in self.invalids) & \
                        ([amin, neig_amin] not in self.exclusions and [neig_amin, amin] not in self.exclusions):

                    coord_1 = np.array(self.aromatic_array[amin])
                    coord_2 = np.array(self.aromatic_array[neig_amin])
                    aromatic_distance = np.linalg.norm(coord_1 - coord_2)
                    if aromatic_distance < 5.5 and amin != neig_amin:
                        pair = (int(residue.id[1]), int(neig_res.id[1]))

                        if pair in self.analyzed_pairs:
                            continue
                        else:
                            self.analyzed_pairs.add((int(neig_res.id[1]), int(residue.id[1])))

                        normal_1 = self.aromatic_normals[amin] / np.linalg.norm(self.aromatic_normals[amin])
                        normal_2 = self.aromatic_normals[neig_amin] / np.linalg.norm(self.aromatic_normals[neig_amin])
                        angle = np.arccos(np.clip(np.dot(normal_1, normal_2), -1.0, 1.0))
                        if angle > 50:
                            # Tshaped
                            orient_type = 'T'
                        elif 30 < angle < 50:
                            # Inter (stacked no parallel)
                            orient_type = 'I'
                        elif angle < 30:
                            # Parallel
                            orient_type = 'P'
                        chain1 = 'MC' if len(atom.get_name()) == 1 else 'SC'
                        chain2 = 'MC' if len(neighbor.get_name()) == 1 else 'SC'
                        coord_1 = f'{coord_1[0]:.3f},{coord_1[1]:.3f},{coord_1[2]:.3f}'
                        coord_2 = f'{coord_2[0]:.3f},{coord_2[1]:.3f},{coord_2[2]:.3f}'
                        if self.multiple:
                            self.bonds_check.append((f"{chain.id}:{str(residue.id[1])}:_:{str(residue.resname)}",
                                                     f"{neig_res.get_parent().id}:{str(neig_res.id[1])}:_:{str(neig_res.resname)}"))
                            self.add_bond([
                                chain, residue, neig_res,
                                f"PIPISTACK:{chain1}_{chain2}",
                                f"{aromatic_distance:.3f}",
                                f"{angle:.3f}",
                                f"{9.4:.3f}",
                                coord_1,
                                coord_2,
                                "NaN",
                                f"{orient_type}"
                            ], 'pi_stacking')
                        self.exclusions.append([amin, neig_amin])

    def Bonds(self):

        for chain in self.structure.get_chains():
            for residue in chain:

                if residue.resname in ['032', 'HOH']:
                    continue
                for atom in residue:
                    atom_name = atom.get_name()
                    is_vdw = False

                    # Looking for HBOND
                    self._hydrogen_bond(chain, residue, atom)
                    # Looking for VDW
                    self._vanderwaals(chain, residue, atom)
                    # Looking for SBOND
                    self._dissulfide_bond(chain, residue, atom)
                    # Salt Bridges
                    self._salt_bridge(chain, residue, atom)
                    # Pi Stacking
                    self._pi_stacking(chain, residue, atom)

    def analyse(self, bond, lig):
        for pair in self.bonds_check:
            pair_dist, pair_idx = [], []

            for line in range(len(self.nodes_id1)):

                if (pair == (self.nodes_id1[line], self.nodes_id2[line]) and (bond in self.bonds[line])):
                    pair_dist.append(self.distances[line])
                    pair_idx.append(line)

            if len(pair_dist) > 1:
                min_idx = np.argmin(pair_dist)
                min_pair = pair_idx[min_idx]
                for i in pair_idx:
                    if i != min_pair:
                        self.nodes_id1.pop(i)
                        self.nodes_id2.pop(i)
                        self.donors.pop(i)
                        self.angles.pop(i)
                        self.energies.pop(i)
                        self.bonds.pop(i)
                        self.distances.pop(i)
                        self.atom1.pop(i)
                        self.atom2.pop(i)
                        self.orientation.pop(i)
                        self.ligands[lig] -= 1

    def multiple_mode(self):
        bonds = [("HBOND", "hb"), ("VDW", "vdw"), ("SBOND", "sbond"), ("IONIC", "ionic"),
                 ("PIPISTACK", "pi_stacking")]
        for b in bonds:
            self.analyse(b[0], b[1])

    def to_file(self):
        self.Bonds()
        if self.multiple:
            self.multiple_mode()
        colunas = ["NodeId1", "Interaction", "NodeId2", "Distance", "Angle", "Energy", "Atom1", "Atom2", "Donor"]

        data = pd.DataFrame(list(zip(self.nodes_id1, self.bonds, self.nodes_id2, self.distances,
                                     self.angles, self.energies, self.atom1, self.atom2, self.donors)), columns=colunas)

        data.to_csv(f'./{self.name}_edges.csv', sep='\t', index=False)

    def print_output(self):
        self.Bonds()

        if self.multiple:
            self.multiple_mode()

        print(len(self.nodes_id1), len(self.donors))
        time.sleep(2)
        for n in range(len(self.nodes_id1)):
            try:
                print(
                    f"{self.nodes_id1[n]}\t{self.bonds[n]}\t{self.nodes_id2[n]}\t{self.distances[n]}\t{self.angles[n]}\t\t{self.energies[n]}\t\t{self.atom1[n]}\t{self.atom2[n]}\t{self.donors[n]}")
            except Exception as e:
                print(e)
                print(
                    f"{self.nodes_id1[n]}\t{self.bonds[n]}\t{self.nodes_id2[n]}\t{self.distances[n]}\t{self.angles[n]}\t\t{self.energies[n]}\t\t{self.atom1[n]}\t{self.atom2[n]}\t{self.donors[n]}")
        print(self.ligands)
