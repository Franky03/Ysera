import pandas as pd
import numpy as np
import os
from sklearn.metrics.pairwise import euclidean_distances
from biopandas.pdb import PandasPdb


class AromaticsFormat:
    """Gera uma matriz com as coordenadas e todas as distâncias de par em par dos átomos,
    também gera o aromatic arrays e o aromatic normals"""
    def __init__(self, filename):
        self.filename = filename
        self.project_home = os.path.dirname(os.path.realpath(__file__))
        self.path = self.project_home + '/temp/' + self.filename
        self.aromatic_pos = []
        self.aromatic_points = []
        self.invalids = []
        self.aromatic_array = {}
        self.aromatic_normals = {}
        self.df_total = pd.DataFrame()
        self.aminos = pd.DataFrame()
        self.arom_phe_tyr = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
        self.arom_trp = ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']

    def _formata_arquivo(self, path):
        """Formata o dataframe inicial usando o biopandas, cria um dataframe só com os aminoácidos
        e os átomos necessários"""
        file = open(self.path, 'r')
        new_name = f'{self.filename}new.pdb'
        with open(new_name, 'w') as f:
            for line in file:
                if "ENDMDL" in line:
                    break
                else:
                    f.write(line)
        ppdb = PandasPdb()
        ppdb.read_pdb(path)
        os.remove(self.project_home + '/' + new_name)
        atom = ppdb.df['ATOM']
        hetatm = ppdb.df['HETATM']
        # Cria um dataframe apenas com ATOM E HETATM
        self.df_total = pd.concat([atom, hetatm], sort=False)
        self.df_total = self.df_total.reset_index()
        self.df_total = self.df_total.drop(['index'], axis=1)
        # Gera uma nova coluna 'amin' que vai ter o 'id' da cadeia e o número do resíduo
        self.df_total['amin'] = self.df_total['chain_id'] + ' ' + self.df_total['residue_number'].astype(str)
        self.df_total = self.df_total.drop(['atom_number', 'b_factor', 'alt_loc', 'line_idx',
                                            'occupancy', 'element_symbol', 'charge', 'insertion',
                                            'segment_id', 'blank_1', 'blank_2', 'blank_3', 'blank_4',
                                            'chain_id', 'residue_number'], axis=1)
        self.aminos = self.df_total[self.df_total['residue_name'].isin(['TYR', 'PHE', 'TRP'])]
        self.aminos = self.aminos.loc[
            (self.aminos['residue_name'].isin(['TYR'])) & (self.aminos['atom_name'].isin(self.arom_phe_tyr)) |
            (self.aminos['residue_name'].isin(['PHE'])) & (self.aminos['atom_name'].isin(self.arom_phe_tyr)) |
            (self.aminos['residue_name'].isin(['TRP'])) & (self.aminos['atom_name'].isin(self.arom_trp))]

    def _calcula_array(self, amin):
        """Calcula o aromaticpos e o aromaticpoints de um determinado Amin"""
        df = self.aminos[self.aminos['amin'].isin([amin])]
        self.aromatic_pos = [df.iloc[0]['x_coord'], df.iloc[0]['y_coord'], df.iloc[0]['z_coord']]
        for index, linha in df.iterrows():
            coordenada = self._gera_coord(linha)
            self.aromatic_pos = [(x + y) / 2 for x, y in zip(self.aromatic_pos, coordenada)]
            self.aromatic_array[amin] = self.aromatic_pos

        if len(df) < 3:
            self.invalids.append(amin)
        else:
            self.aromatic_points = self._gera_coord(df.iloc[0:3])
            veca = np.subtract(self.aromatic_points[1], self.aromatic_points[0])
            vecb = np.subtract(self.aromatic_points[2], self.aromatic_points[0])
            self.aromatic_normals[amin] = np.cross(veca, vecb)

    def _gera_coord(self, linhas):
        """Retorna as coordenadas de uma determinada linha ou listas de linhas"""
        coord = []
        if isinstance(linhas, pd.DataFrame):
            for index, linha in linhas.iterrows():
                coord.append([linha['x_coord'], linha['y_coord'], linha['z_coord']])
        else:
            coord = [linhas['x_coord'], linhas['y_coord'], linhas['z_coord']]
        return coord

    def _calcula_dist(self):
        """Calcula a distância euclidiana de todos os átomos de par em par"""
        dist = euclidean_distances(
            np.float32(self.df_total[["x_coord", "y_coord", "z_coord"]].to_numpy()),
            np.float32(self.df_total[["x_coord", "y_coord", "z_coord"]].to_numpy()))
        df_dist = pd.DataFrame(data=dist)
        df_dist = pd.concat([self.df_total, df_dist], axis=1, sort=False)
        return df_dist

    def get_data(self):
        """Roda os métodos da classe e retorna o dataframe final com todas as distâncias
        além do aromatic array e aromatic normals"""
        self._formata_arquivo(self.path)
        amin_list = list(dict.fromkeys(self.aminos['amin'].values))
        for i in amin_list:
            self._calcula_array(i)
        df_total_dist = self._calcula_dist()
        return self.aromatic_array, self.aromatic_normals, self.invalids, self.df_total, df_total_dist


if __name__ == '__main__':
    af = AromaticsFormat('file_30.pdb')
    array, normals, invalids,  total, total_dist = af.get_data()
