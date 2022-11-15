import numpy as np
import pandas as pd


class Thread:
    def __init__(self, aromatic_array, aromatic_normals, invalids, total, total_dist):
        self.aromatic_array = aromatic_array
        self.aromatic_normals = aromatic_normals
        self.invalids = invalids
        self.total = total
        self.total_dist = total_dist
        self.line = pd.Series(dtype='float64')
        self.text = ''
        self.dist_old = []
        self.exclusions = []
        self.params = {
            'Hydrogen_Bond': [3.1, 0],
            'Salt_Bridge': [4.0, 0],
            'Dissulfide_Bond': [2.2, 0],
            'Van_der_Waals': [3.2, 0],
            'Pi_Stacking': [7.2, 0],
            'Sulfur_Aryl': [7.2, 0],
            'Cation_Aryl': [6, 0],
            'Anion_Aryl': [3.5, 0],
            "tshaped": 0,
            "inter": 0,
            "paralel": 0,
            'lighbacep': ['OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1'],
            'lighbdono': ['HE1', 'H', 'HE', 'HH11', 'HH12', 'HH21',
                          'HH22', 'HD1', 'HE2', 'HG', 'HG1', 'HG21', 'HG22',
                          'HG23', 'HH', 'HD21', 'HD22', 'HE21', 'HE22'],
            'ligsb1': ['OD2', 'OE2', 'OH'],
            'ligsb2': ['NZ', 'NH2', 'NH1'],
            'aasbneg': ['ASP', 'GLU'],
            'aasbpos': ['LYS', 'ARG'],
            'ligvdw': ['CB', 'CG1', 'CG2', 'CD1', 'CD2', 'CE'],
            'aavdw': ['VAL', 'TRE', 'MET', 'LEU', 'ILE'],
            'aa_arom': ['TYR', 'PHE', 'TRP'],
            'ligctn': ['MG', 'CU', 'K', 'FE2', 'FE', 'NI', 'NA', 'MO1', 'MO3', 'MO4',
                       'MO5', 'MO6', 'MO7', 'MO8', 'MO9', 'NZ', 'NH2', 'NH1'],
            'ligctn2': ['CG', 'CE2', 'CG'],
            'aactn_beg': 3.4,
            'aactn_end': 4.0,
            'ligan1': ['CL', 'BR', 'I', 'OD2', 'OE2', 'OH'],
            'ligspi1': ['SG'],
            'ligspi2': ['CG', 'CE2', 'CG'],
            'aaan_beg': 2.0,
            'aaan_end': 3.0,
            'aaspi': 5.3
            }

    def _format_line(self, n_line):
        """Retorna uma linha contendo apenas as distâncias possíveis de serem ligações"""
        local_line = self.total_dist.iloc[n_line]
        colunas_total = list(self.total.columns)
        colunas_liga = np.asarray(np.where(local_line.loc[0: (len(self.total_dist) - 1)] <= 7.2))
        index = np.where(colunas_liga == n_line)
        colunas_liga = np.delete(colunas_liga, index[1])
        for i in colunas_liga:
            colunas_total.append(i)
        local_line = local_line.filter(colunas_total)
        return local_line

    def _atom_info(self, param):
        """Recebe a distância máxima da ligação a ser checada, retorna apenas
        distâncias daquela ligação com os e as inforomações da linha analisada
        """
        dist_line = []
        dist_value = []
        atom_1 = self.line['atom_name']
        aa_1 = self.line['residue_name']
        amin_1 = self.line['amin']
        index = np.where(self.line.iloc[7: len(self.line)] <= self.params[param][0])
        if len(index[0]) != 0:
            for i in index[0]:
                dist_value.append(self.line.iloc[i+7])
                dist_line.append(self.line.index[i+7])

            return atom_1, aa_1, amin_1, dist_line, dist_value
        else:
            return 0, 0, 0, 0

    def _format_text(self, bond, atom_1, aa_1, amin_1, atom_2, aa_2, amin_2, dist):
        """Formata o texto que será colocado no arquivo"""
        if dist == 0:
            self.text = self.text + f'{bond}\t\t{atom_1}\t\t{aa_1}\t\t{amin_1}\t\t ' \
                                    f'{atom_2}\t\t {aa_2}\t\t{amin_2}\t\t{dist}\n'
            self.params['Pi_Stacking'][1] += 1
        elif dist not in self.dist_old:
            self.dist_old.append(dist)
            self.params[bond][1] += 1
            self.text = self.text + f'{bond}\t\t{atom_1}\t\t{aa_1}\t\t{amin_1}\t\t ' \
                                    f'{atom_2}\t\t {aa_2}\t\t{amin_2}\t\t{dist}\n'
    # Deteção das ligações:

    def _hydrogen_bond(self):
        bond = 'Hydrogen_Bond'
        atom_1, aa_1, amin_1, dist_line, dist_value = self._atom_info(bond)
        if dist_value != 0:
            for i, j in zip(dist_line, dist_value):
                atom_2 = self.total.iloc[i]['atom_name']
                aa_2 = self.total.iloc[i]['residue_name']
                amin_2 = self.total.iloc[i]['amin']
                if (atom_1 in self.params['lighbacep']) & (atom_2 in self.params['lighbdono']):
                    self._format_text(bond, atom_1, aa_1,
                                      amin_1, atom_2, aa_2, amin_2, j)
                elif (atom_2 in self.params['lighbacep']) & (atom_1 in self.params['lighbdono']):
                    self._format_text(bond, atom_1, aa_1,
                                      amin_1, atom_2, aa_2, amin_2, j)

    def _salt_bridge(self):
        bond = 'Salt_Bridge'
        atom_1, aa_1, amin_1, dist_line, dist_value = self._atom_info(bond)
        if dist_value != 0:
            for i, j in zip(dist_line, dist_value):
                atom_2 = self.total.iloc[i]['atom_name']
                aa_2 = self.total.iloc[i]['residue_name']
                amin_2 = self.total.iloc[i]['amin']

                if (aa_1 in self.params['aasbpos'] and atom_1 in
                        self.params['ligsb2']) & (aa_2 in self.params['aasbneg'] and atom_2 in self.params['ligsb1']):
                    self._format_text(bond, atom_1, aa_1, amin_1,
                                      atom_2, aa_2, amin_2, j)

                elif (aa_1 in self.params['aasbneg'] and atom_1 in
                        self.params['ligsb1']) & (aa_2 in self.params['aasbpos'] and atom_2 in self.params['ligsb2']):
                    self._format_text(bond, atom_1, aa_1, amin_1,
                                      atom_2, aa_2, amin_2, j)

    def _dissulfide_bond(self):
        bond = 'Dissulfide_Bond'
        atom_1, aa_1, amin_1, dist_line, dist_value = self._atom_info(bond)
        if dist_value != 0:
            if atom_1 == 'SG':
                for i, j in zip(dist_line, dist_value):

                    if self.total.iloc[i]['atom_name'] == 'SG':
                        atom_2 = self.total.iloc[i]['atom_name']
                        aa_2 = self.total.iloc[i]['residue_name']
                        amin_2 = self.total.iloc[i]['amin']
                        self._format_text(bond, atom_1, aa_1,
                                          amin_1, atom_2, aa_2, amin_2, j)

    def _vanderwaals(self):
        bond = 'Van_der_Waals'
        atom_1, aa_1, amin_1, dist_line, dist_value = self._atom_info(bond)
        if dist_value != 0:
            for i, j in zip(dist_line, dist_value):
                atom_2 = self.total.iloc[i]['atom_name']
                aa_2 = self.total.iloc[i]['residue_name']
                amin_2 = self.total.iloc[i]['amin']
                if (aa_1 in self.params['aavdw'] and atom_1 in
                        self.params['ligvdw']) & (aa_2 in self.params['aavdw'] and atom_2 in self.params['ligvdw']):
                    self._format_text(bond, atom_1, aa_1,
                                      amin_1, atom_2, aa_2, amin_2, j)

    def _pi_stacking(self):
        bond = 'Pi_Stacking'
        atom_1, aa_1, amin_1, dist_line, dist_value = self._atom_info(bond)
        if dist_value != 0:
            for i, j in zip(dist_line, dist_value):
                atom_2 = self.total.iloc[i]['atom_name']
                aa_2 = self.total.iloc[i]['residue_name']
                amin_2 = self.total.iloc[i]['amin']
                if aa_1 in self.params['aa_arom'] and aa_2 in self.params['aa_arom']:
                    if (amin_1 not in self.invalids and amin_2 not in self.invalids) &\
                            ([amin_1, amin_2] not in self.exclusions and [amin_2, amin_1] not in self.exclusions):
                        coord_1 = np.array(self.aromatic_array[amin_1])
                        coord_2 = np.array(self.aromatic_array[amin_2])
                        aromatic_distance = np.linalg.norm(coord_1 - coord_2)
                        if aromatic_distance < self.params['aaspi']:
                            self._format_text(bond, atom_1, aa_1,
                                              amin_1, atom_2, aa_2, amin_2, aromatic_distance)
                            normal_1 = self.aromatic_normals[amin_1] / np.linalg.norm(self.aromatic_normals[amin_1])
                            normal_2 = self.aromatic_normals[amin_2] / np.linalg.norm(self.aromatic_normals[amin_2])
                            angle = np.arccos(np.clip(np.dot(normal_1, normal_2), -1.0,
                                                      1.0))
                            if angle > 50:
                                self.params['tshaped'] += 1
                            elif 30 < angle < 50:
                                self.params['inter'] += 1
                            elif angle < 30:
                                self.params['paralel'] += 1
                            self.exclusions.append([amin_1, amin_2])

    def _cation_aryl(self):
        bond = 'Cation_Aryl'
        atom_1, aa_1, amin_1, dist_line, dist_value = self._atom_info(bond)
        if dist_value != 0:
            for i, j in zip(dist_line, dist_value):
                atom_2 = self.total.iloc[i]['atom_name']
                aa_2 = self.total.iloc[i]['residue_name']
                amin_2 = self.total.iloc[i]['amin']
                if aa_1 in self.params['aa_arom'] and atom_2 in self.params['ligctn']:
                    if (amin_1 not in self.invalids and amin_2 not in self.invalids) & \
                            ([amin_1, amin_2] not in self.exclusions and [amin_2, amin_1] not in self.exclusions):
                        coord_1 = np.array(self.aromatic_array[amin_1])
                        coord_2 = np.array([self.total.iloc[i]['x_coord'],
                                            self.total.iloc[i]['y_coord'], self.total.iloc[i]['z_coord']])
                        aromatic_distance = np.linalg.norm(coord_1 - coord_2)
                        if self.params['aactn_beg'] < aromatic_distance < self.params['aactn_end']:
                            self._format_text(bond, atom_1, aa_1,
                                              amin_1, atom_2, aa_2, amin_2, aromatic_distance)
                            self.exclusions.append([amin_1, amin_2])
                elif atom_1 in self.params['ligctn'] and aa_2 in self.params['aa_arom']:
                    if (amin_1 not in self.invalids and amin_2 not in self.invalids) & \
                            ([amin_1, amin_2] not in self.exclusions and [amin_2, amin_1] not in self.exclusions):
                        coord_1 = np.array([self.line['x_coord'], self.line['y_coord'], self.line['z_coord']])
                        coord_2 = np.array(self.aromatic_array[amin_2])
                        aromatic_distance = np.linalg.norm(coord_1 - coord_2)
                        if self.params['aactn_beg'] < aromatic_distance < self.params['aactn_end']:
                            self._format_text(bond, atom_1, aa_1,
                                              amin_1, atom_2, aa_2, amin_2, aromatic_distance)
                            self.exclusions.append([amin_1, amin_2])

    def _sulfur_aryl(self):
        bond = 'Sulfur_Aryl'
        atom_1, aa_1, amin_1, dist_line, dist_value = self._atom_info(bond)
        if dist_value != 0:
            for i, j in zip(dist_line, dist_value):
                atom_2 = self.total.iloc[i]['atom_name']
                aa_2 = self.total.iloc[i]['residue_name']
                amin_2 = self.total.iloc[i]['amin']
                if aa_1 in self.params['aa_arom'] and atom_2 in self.params['ligspi1']:
                    if (amin_1 not in self.invalids and amin_2 not in self.invalids) & \
                            ([amin_1, amin_2] not in self.exclusions and [amin_2, amin_1] not in self.exclusions):
                        coord_1 = np.array(self.aromatic_array[amin_1])
                        coord_2 = np.array([self.total.iloc[i]['x_coord'],
                                            self.total.iloc[i]['y_coord'], self.total.iloc[i]['z_coord']])
                        aromatic_distance = np.linalg.norm(coord_1 - coord_2)
                        if 0 < aromatic_distance < self.params['aaspi']:
                            self._format_text(bond, atom_1, aa_1,
                                              amin_1, atom_2, aa_2, amin_2, aromatic_distance)
                            self.exclusions.append([amin_1, amin_2])
                        elif atom_1 in self.params['ligspi1'] and aa_2 in self.params['aaspi']:
                            coord_1 = np.array([self.line['x_coord'], self.line['y_coord'], self.line['z_coord']])
                            coord_2 = np.array(self.aromatic_array[amin_2])
                            aromatic_distance = np.linalg.norm(coord_1 - coord_2)
                            if 0 < aromatic_distance < self.params['aaspi']:
                                self._format_text(bond, atom_1, aa_1,
                                                  amin_1, atom_2, aa_2, amin_2, aromatic_distance)
                                self.exclusions.append([amin_1, amin_2])

    def _anion_aryl(self):
        bond = 'Anion_Aryl'
        atom_1, aa_1, amin_1, dist_line, dist_value = self._atom_info(bond)
        if dist_value != 0:
            for i, j in zip(dist_line, dist_value):
                atom_2 = self.total.iloc[i]['atom_name']
                aa_2 = self.total.iloc[i]['residue_name']
                amin_2 = self.total.iloc[i]['amin']
                if aa_1 in self.params['aa_arom'] and atom_2 in self.params['ligan1']:
                    if (amin_1 not in self.invalids and amin_2 not in self.invalids) & \
                            ([amin_1, amin_2] not in self.exclusions and [amin_2, amin_1] not in self.exclusions):
                        coord_1 = np.array(self.aromatic_array[amin_1])
                        coord_2 = np.array([self.total.iloc[i]['x_coord'],
                                            self.total.iloc[i]['y_coord'], self.total.iloc[i]['z_coord']])
                        aromatic_distance = np.linalg.norm(coord_1 - coord_2)
                        if self.params['aaan_beg'] < aromatic_distance < self.params['aaan_end']:
                            self._format_text(bond, atom_1, aa_1,
                                              amin_1, atom_2, aa_2, amin_2, aromatic_distance)
                            self.exclusions.append([amin_1, amin_2])
                elif atom_1 in self.params['ligan1'] and aa_2 in self.params['aa_arom']:
                    if (amin_1 not in self.invalids and amin_2 not in self.invalids) & \
                            ([amin_1, amin_2] not in self.exclusions and [amin_2, amin_1] not in self.exclusions):
                        coord_1 = np.array([self.line['x_coord'], self.line['y_coord'], self.line['z_coord']])
                        coord_2 = np.array(self.aromatic_array[amin_2])
                        aromatic_distance = np.linalg.norm(coord_1 - coord_2)
                        if self.params['aaan_beg'] < aromatic_distance < self.params['aaan_end']:
                            self._format_text(bond, atom_1, aa_1,
                                              amin_1, atom_2, aa_2, amin_2, aromatic_distance)
                            self.exclusions.append([amin_1, amin_2])

    def run(self, beg, end, name):
        for i in range(beg, end):
            self.line = self._format_line(i)
            self._hydrogen_bond()
            self._salt_bridge()
            self._dissulfide_bond()
            self._vanderwaals()
            self._pi_stacking()
            self._cation_aryl()
            # self._sulfur_aryl()
            # self._anion_aryl()
        f = open('output2.0/' + name + ".txt", "w")
        f.write(self.text)
        f.close()
        print(f'Hydrogen_Bond: {self.params["Hydrogen_Bond"][1]}; '
              f'Salt_Bridge: {self.params["Salt_Bridge"][1]}; Dissulfide_Bond: {self.params["Dissulfide_Bond"][1]}; '
              f'Van_der_Waals: {self.params["Van_der_Waals"][1]}; Pi_Stacking: {self.params["Pi_Stacking"][1]}; '
              f'Cation_Aryl: {self.params["Cation_Aryl"][1]}; '
              f'Sulfur_Aryl: {self.params["Sulfur_Aryl"][1]}; Anion_Aryl: {self.params["Anion_Aryl"][1]}')
