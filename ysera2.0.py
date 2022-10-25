from turtle import distance
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
PS_DEFAULT = 7.2
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

def mythread(New, params, i, a, filename, string2):
    d= AromaticArray.copy()
    n= AromaticNormals.copy()
    f= open('output/' + filename + '.txt', 'w+')
    while i< a:
        for j in range(i + 1, len(New)):
            distance= New[j].iloc[i]
            print(f"top distance: {distance}\n")
            if distance > 8 or distance==0:
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
            print(tuple(hb,sb,db,lpi,tshaped,inter,paralel,vdw,ctn,an,spi))

            atom1= New['Atom Type'].iloc[i]
            atom2= New['Atom Type'].iloc[j]
            #O que é aa1 e aa2?
            aa1= New['aa'].iloc[i]
            aa2= New['aa'].iloc[j]
            #O que é chaincode?
            chaincode1= New['Chain ID'].iloc[i]
            chaincode2= New['Chain ID'].iloc[j]

            distance= New[j].iloc[i]
            print(f"down distance: {distance}\n")
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
                pass