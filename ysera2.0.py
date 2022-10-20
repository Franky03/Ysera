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
    manager.start()
    d = manager.dict()
    n = manager.dict()

def mythread(New, params, i, a, filename, string2):
    d= AromaticArray.copy()
    n= AromaticArray.copy()
    f= open('output/' + filename + '.txt', 'w+')
    