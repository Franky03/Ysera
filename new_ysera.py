import pandas as pd
import numpy as np
import math
import os
import asyncio
from os import listdir
from os.path import isfile, join
from sklearn.metrics.pairwise import euclidean_distances
import time
import threading
import _thread
import gzip
import json
import concurrent.futures
from multiprocessing import set_start_method, Pool, Manager

PROJECT_HOME= os.path.dirname(os.path.realpath(__file__)) #Pega a pasta em que o arquivo atual está

HB_DEFAULT= 3.1
SB_DEFAULT = 4.0
DB_DEFAULT = 2.2
VDW_DEFAULT = 5.5
PS_DEFAULT = 7.2 #Pi Stacking, variável definida como lpi
AAAN_BEG_DEFAULT = 2.0
AAAN_END_DEFAULT = 3.0
AASPI_DEFAULT = 5.3
AACTN_BEG_DEFAULT = 3.4
AACTN_END_DEFAULT = 4.0

hb=0 #Hydrogen Bond

# Nas listas são os nomes de alguns átomos presentes no arquivo pdb
lighbacep = ['OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1'] 
lighbdono = ['HE1', 'H', 'HE', 'HH11', 'HH12', 'HH21', 
            'HH22', 'HD1', 'HE2', 'HG', 'HG1', 'HG21', 'HG22',
             'HG23', 'HH', 'HD21', 'HD22', 'HE21', 'HE22']

sb=0 #Salt Bridge

ligsb1 = ['OD2', 'OE2', 'OH']
ligsb2 = ['NZ', 'NH2', 'NH1']
# Nas listas são tipos de aminoácidos
aasbneg = ['ASP', 'GLU']
aasbpos = ['LYS', 'ARG']


db= 0 #Dissulfide Bond

ligdb = ['SG']

vdw = 0 #Van der Waals
ligvdw = ['CB', 'CG1', 'CG2', 'CD1', 'CD2', 'CE']
aavdw = ['VAL', 'TRE', 'MET', 'LEU', 'ILE']

lpi = 0 #Pi_Stacking
#------------------
#Variações dos angulos da ligação pi: 

tshaped = 0
inter = 0
paralel = 0

#------------------
aapi = ['TYR', 'PHE', 'TRP']

ctn = 0 # Cation Aryl
ligctn = ['MG', 'CU', 'K', 'FE2', 'FE', 'NI', 'NA', 'MO1', 'MO3', 'MO4', 'MO5', 'MO6', 'MO7', 'MO8', 'MO9',
          'NZ', 'NH2', 'NH1']
ligctn2 = ['CG', 'CE2', 'CG']
aactn = ['TYR', 'PHE', 'TRP']

an = 0 # Anion Aryl
ligan1 = ['CL', 'BR', 'I', 'OD2', 'OE2', 'OH']
ligan2 = ['CG', 'CE2', 'CG']
aaan = ['TYR', 'PHE', 'TRP']

spi = 0 # Sulfur_Aryl
ligspi1 = ['SG']
ligspi2 = ['CG', 'CE2', 'CG']
aaspi = ['TYR', 'PHE', 'TRP']

Invalids = [] #List of incompletes Aromatic Structures
#Dicts to Get the Aromatics
AromaticArray = {}
AromaticNormals = {}
#--------------------------
Exclusions = []
exitFlag = 0

#Estava acontecendo um erro no Manager

#É preciso colocar essa condição para não dar erro
if __name__== '__main__':
    manager = Manager()
    #manager.start()
    d = manager.dict()
    n = manager.dict()

def myfunction(filename, params):
    string1= ""
    string2= ""
    path= PROJECT_HOME + '/temp/' + filename
    pathoutput = PROJECT_HOME + '/output/' + filename

    #Aromático dos aminoácidos:
    AROMTRP = ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'] # Triptofano
    AROMPHE = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'] # Fenilalanina
    AROMTYR = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'] # Tirosina

    datapdb= pd.DataFrame(columns=['Type', 'Atom ID', 'Atom Type', 'aa', 'Chain ID', 'X', 'Y', 'Z', 'occupancy', 'temperature factor',
                 'element symbol']) 
    
    i=0

    print("Starting to load data")
    start_time = time.time()
    Amin = '' 
    Aromaticpos = [] #Posição do Aromático 
    AromaticPoints = [] #Pontos do Aromatico
    
    file= open(path, 'r')
    print(file)

    #Aqui começaremos a ler o arquivo pdb
    
    for line in file:
        if ("ENDMDL" in line):
            pass
        if (line[16:20].strip() in "HOH"):
            continue
        if ("ATOM" in line[:6] or "HETATM" in line[:6]):
            if (line[17:20].strip() in ['TYR', 'PHE', 'TRP']):
                #Amin é aminoácido ?
                if (Amin == ''):
                    Amin = line[20:27].strip()
                elif (Amin != line[20:27].strip()):
                    AromaticArray[Amin] = np.asarray(Aromaticpos)
                    if (len(AromaticPoints) == 3):
                        veca= np.subtract(AromaticPoints[1], AromaticPoints[0])
                        vecb = np.subtract(AromaticPoints[2], AromaticPoints[0])
                        AromaticNormals[Amin] = np.cross(veca, vecb)
                    else:
                        print("Incomplete Aromatic Structure at {}".format(Amin))
                        Invalids.append(Amin)
                    Amin= line[20:27].strip()
                    Aromaticpos = []
                    AromaticPoints = []
                
                elif ((line[17:20].strip() == "TYR" and line[11:16].strip() in AROMTYR) or (line[17:20].strip() == 'PHE' and line[11:16].strip() in AROMPHE) or (line[17:20].strip() == 'TRP' and line[11:16].strip() in AROMTRP)):
                    if (Aromaticpos == []):
                        Aromaticpos = [float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())] #posição do aromático dos aminoácidos que entraram na condição
                    else:
                        Aromaticpos = [(x + y) / 2 for x, y in zip(Aromaticpos, [float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])]
                    
                    if (len(AromaticPoints) < 3):
                        AromaticPoints.append([float(line[27:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
            datapdb.loc[i] = [line[:6].strip(), line[6:11].strip(), line[11:16].strip(), line[17:20].strip(),line[20:27].strip(), line[27:38].strip(), line[38:46].strip(), line[46:54].strip(), line[54:60].strip(), line[60:66].strip(), line[66:80].strip()]

            # datapdb.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

            i += 1
    
    AromaticArray[Amin] = np.asarray(Aromaticpos)

    if (len(AromaticPoints) == 3):
        veca = np.subtract(AromaticPoints[1], AromaticPoints[0])
        vecb = np.subtract(AromaticPoints[2], AromaticPoints[0])
        AromaticNormals[Amin] = np.cross(veca, vecb)
    elif (0 < len(AromaticPoints) < 3):
        print("Incomplete Aromatic Structure at {}".format(Amin))
        Invalids.append(Amin)
    print("DataLoaded")
    print("---%s seconds ---" % (time.time() - start_time))
    if 'Atom ID' in datapdb.columns:
        datapdb = datapdb.drop(['Atom ID', 'occupancy', 'temperature factor', 'element symbol'], axis=1)
    
    #Tirei o parametro downcast pois é de outra versão do numpy

    datapdb["X"] = np.float32(datapdb["X"])
    datapdb["Y"] = np.float32(datapdb["Y"])
    datapdb["Z"] = np.float32(datapdb["Z"])

    dist= euclidean_distances(np.float32(datapdb[["X", "Y", "Z"]].to_numpy()),
                              np.float32(datapdb[["X", "Y", "Z"]].to_numpy()))

    Dist= pd.DataFrame(data= dist, index=None)

    New = pd.concat([datapdb, Dist], axis=1, sort= False)
    print(New)
    print(len(New))

    return New #Retornará um dataset que irá para a função mythread


def mythread(New, params, i, a, filename, string2):
    d= AromaticArray.copy()
    n= AromaticNormals.copy()
    f= open('output/' + filename + '.txt', 'w') #arquivo de saída
    while i< a:
        for j in range(i + 1, len(New)):
            distance= New[j].iloc[i]
            if (distance > 8 or distance == 0):
                continue
            global hb
            global sb
            global db
            global lpi
            global tshaped
            global inter
            global paralel
            global vdw
            global ctn
            global an
            global spi
            # print(tuple(hb,sb,db,lpi,tshaped,inter,paralel,vdw,ctn,an,spi))

            atom1= New['Atom Type'].iloc[i]
            atom2= New['Atom Type'].iloc[j]
            
            aa1= New['aa'].iloc[i] #aminoácido do atomo 1
            aa2= New['aa'].iloc[j] #aminoácido do atomo 2
            
            chaincode1= New['Chain ID'].iloc[i] #cadeia do atomo 1
            chaincode2= New['Chain ID'].iloc[j] #cadeia do atomo 2

            distance= New[j].iloc[i]
            # print(f"down distance: {distance}\n")
            #Looking for Hydrogen Bond
            
            if atom1 in lighbacep[:] and atom2 in lighbdono[:] and 0.0 < distance < params['hb']:
                hb += 1
                string2 = string2 + ('Hydrogen_Bond' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(distance) + '\n')
            
            elif atom1 in lighbdono[:] and atom2 in lighbacep[:] and 0.0< distance < params['hb'] :
                hb += 1
                string2 = string2 + ('Hydrogen_Bond' + '\t\t' + atom1 + '\t\t'+ aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(distance) + '\n')
            #Looking for Salt Bridge
            
            if aa1 in aasbpos[:] and atom1 in ligsb2[:] and aa2 in aasbneg[:] and atom2 in ligsb1[:] and 0.0< distance< params['sb']:
                sb+=1
                string2 = string2 + ('Salt_Bridge' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 +'\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(distance)+ '\n')
            
            elif aa1 in aasbneg[:] and atom1 in ligsb1[:] and aa2 in aasbpos[:] and atom2 in ligsb2[:] and 0.0 < distance < params['sb']:
                sb+=1
                string2 = string2 + (
                        'Salt_Bridge' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(
                    distance) + '\n')
            #Looking for Dissulfide_Bond

            if atom1 in ligdb[:] and atom2 in ligdb[:] and 0.0< distance < 2.2:
                db += 1
                string2= string2 + ('Dissulfide_bond' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(distance) + '\n')

            #Looking for Van Der Walls

            if aa1 in aavdw[:] and atom1 in ligvdw[:] and aa2 in aavdw[:] and atom2 in ligvdw[:] and 0.0< distance < params['vdw']:
                vdw += 1
                string2 = string2 + ('Van_der_Waals' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(distance) + '\n')
            
            # Pi Stacking

            if aa1 in aapi[:] and aa2 in aapi[:]:
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1, chaincode2] not in Exclusions and [chaincode2, chaincode1] not in Exclusions):
                    coordinates1= d[chaincode1]
                    coordinates2= d[chaincode2]
                    aromaticdistance=  np.linalg.norm(coordinates1 - coordinates2) # calculate one of the eight different matrix norms or one of the vector norms.
                    if(aromaticdistance < params['aaspi']):
                        string2= string2 + ('Pi_stacking  ' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance) + '\n')
                        NormalVector1= n[chaincode1] / np.linalg.norm(n[chaincode1])
                        NormalVector2= n[chaincode2] / np.linalg.norm(n[chaincode2])
                        Angle= np.arccos(np.clip(np.dot(NormalVector1, NormalVector2), -1.0, 1.0)) #clip tem a fução de retornar um valor entre um intervalo, nesse caso entre -1 e 1

                        if (Angle> 50):
                            tshaped += 1
                        elif (30< Angle < 50):
                            inter += 1
                        elif (Angle<30):
                            paralel += 1
                        
                        lpi += 1
                        Exclusions.append([chaincode1, chaincode2])

            # Looking for Cation Aryl 
            if aa1 in aactn[:] and atom2 in ligctn[:]:
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1, chaincode2] not in Exclusions and [chaincode2, chaincode1] not in Exclusions):
                    coordinates1= d[chaincode1]
                    coordinates2= np.array([New['X'].iloc[j], New['Y'].iloc[j], New['Z'].iloc[j]]) #Coordenadas 2 colocadas em um array 
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)

                    #Se a distancia aromática estiver entre aa começo e o aa final
                    if (params['aactn_beg'] < aromaticdistance < params['aactn_end']):
                        ctn+=1
                        string2 = string2 + ('Cation_Aryl' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])

            elif (atom1 in ligctn[:] and aa2 in aactn[:]):
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1, chaincode2] not in Exclusions and [chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = np.array([New['X'].iloc[i], New['Y'].iloc[i], New['Z'].iloc[i]])
                    coordinates2 = d[chaincode2]
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (params['aactn_beg']< aromaticdistance < params['aactn_end']):
                        ctn += 1
                        string2= string2 + ('Cation_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            #Looking for Sulfur Aryl
            if (aa1 in aaspi[:] and atom2 in ligspi1[:]):
                if(chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1, chaincode2] not in Exclusions and [chaincode2, chaincode1] not in Exclusions):
                    coordinates1= d[chaincode1]
                    coordinates2= np.array([New['X'].iloc[j], New['Y'].iloc[j], New['Z'].iloc[j]])
                    aromaticdistance= np.linalg.norm(coordinates1 - coordinates2)
                    if (0 < aromaticdistance < params['aaspi']):
                        spi += 1
                        string2= string2 + ('Sulfur_Aryl  ' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            elif (atom1 in ligspi1[:] and aa2 in aaspi[:]):
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1, chaincode2] not in Exclusions and [chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = np.array([New['X'].iloc[i], New['Y'].iloc[i], New['Z'].iloc[i]])
                    coordinates2 = d[chaincode2]
                    aromaticdistance= np.linalg.norm(coordinates1 - coordinates2)
                    if(0 < aromaticdistance < params['aaspi']): 
                        spi += 1
                        string2 = string2 + ('Sulfur_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            #Looking for Anion Aryl
            if aa1 in aaan[:] and atom2 in ligan1[:]:
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,chaincode2] not in Exclusions and [chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = d[chaincode1]
                    coordinates2 = np.array([New['X'].iloc[j], New['Y'].iloc[j], New['Z'].iloc[j]])
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (params['aaan_beg'] < aromaticdistance < params['aaan_end']):
                        an += 1
                        string2 = string2 + ('Anion_Aryl' + '\t\t' + 'centroid' + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + atom2 + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
            elif (atom1 in ligan1[:] and aa2 in aaan[:]):
                if (chaincode1 not in Invalids and chaincode2 not in Invalids and [chaincode1,chaincode2] not in Exclusions and [chaincode2, chaincode1] not in Exclusions):
                    coordinates1 = np.array([New['X'].iloc[i], New['Y'].iloc[i], New['Z'].iloc[i]])
                    coordinates2 = d[chaincode2]
                    aromaticdistance = np.linalg.norm(coordinates1 - coordinates2)
                    if (params['aaan_beg']< aromaticdistance < params['aaan_end']):
                        an += 1
                        string2 = string2 + ('Anion_Aryl' + '\t\t' + atom1 + '\t\t' + aa1 + '\t\t' + chaincode1 + '\t\t' + 'centroid' + '\t\t' + aa2 + '\t\t' + chaincode2 + '\t\t' + str(aromaticdistance) + '\n')
                        Exclusions.append([chaincode1, chaincode2])
        i += 1
    print("Teste Aqui !")

    f.write(string2)
    f.close()
    #Formatando a string para a saída desejada
    
    try:
        final= pd.read_table(('output/' + filename + '.txt'), delimiter="\t\t", header= None, encoding='utf-8', engine='python')
        colunas = ["Interaction", "Atom1", "AA1", "Chaincode1", "Atom2", "AA2", "Chaincode2", "Distance"]
        for i in range(len(colunas)):
            final.rename(columns={i: colunas[i]}, inplace=True)
        final.to_csv(('output/' + filename + '.txt'), sep='\t', index=False)
    except Exception as e:
        print(e)
        
    string1= {
        "filename": filename,
        "hb": hb,
        "sb": sb,
        "db": db,
        "lpi": lpi,
        "tshaped": tshaped,
        "inter": inter,
        "paralel": paralel,
        "vdw": vdw,
        "ctn": ctn,
        "an": an,
        "spi": spi
    }
    print(string1) #printando o total de cada tipo de ligação
    

def ysera(filename, params):
    #Corrigindo alguns valores das ligações nos parâmetros casos eles não estejam presentes

    if(not ("hb" in params)):
        params['hb']= HB_DEFAULT
    if (not ("sb" in params)):
        params['sb'] = SB_DEFAULT
    if (not ("db" in params)):
        params['db'] = DB_DEFAULT
    if (not ("vdw" in params)):
        params['vdw'] = VDW_DEFAULT
    if (not ("ps" in params)):
        params['ps'] = PS_DEFAULT
    if (not ("aaan_beg" in params)):
        params['aaan_beg'] = AAAN_BEG_DEFAULT
    if (not ("aaan_end" in params)):
        params['aaan_end'] = AAAN_END_DEFAULT
    if (not ("aaspi" in params)):
        params['aaspi'] = AASPI_DEFAULT
    if (not ("aactn_beg" in params)):
        params['aactn_beg'] = AACTN_BEG_DEFAULT
    if (not ("aactn_end" in params)):
        params['aactn_end'] = AACTN_END_DEFAULT

    print(params)

    res= myfunction(filename, params) #O que faz essa função?
    return res


def lerDados(filename, params):

    if (not ("hb" in params)):
        params['hb'] = HB_DEFAULT #Hydrogen Bond Default
    if (not ("sb" in params)):
        params['sb'] = SB_DEFAULT #Salt Brige Default 
    if (not ("db" in params)):
        params['db'] = DB_DEFAULT #Dissulfide Bond
    if (not ("vdw" in params)):
        params['vdw'] = VDW_DEFAULT #Van de Waals
    if (not ("ps" in params)):
        params['ps'] = PS_DEFAULT 
    if (not ("aaan_beg" in params)):
        params['aaan_beg'] = AAAN_BEG_DEFAULT #Anioan Aril aminoacido 
    if (not ("aaan_end" in params)):
        params['aaan_end'] = AAAN_END_DEFAULT
    if (not ("aaspi" in params)):
        params['aaspi'] = AASPI_DEFAULT
    if (not ("aactn_beg" in params)):
        params['aactn_beg'] = AACTN_BEG_DEFAULT
    if (not ("aactn_end" in params)):
        params['aactn_end'] = AACTN_END_DEFAULT

    new = myfunction(filename, params)

    return new, params