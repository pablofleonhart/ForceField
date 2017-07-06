import copy
from math import *
import random

def calculaLennarJones(a1, a2, dist):
    #a1 e a2: informacoes dos atomos ['atom', sigma, epsilon] --> bib: ffnonbonded
    #dist: distancia entre os atomos, valor retornado da funcao 
    sigma = (a1[1] + a2[1]) / 2 #media aritmetica
    epsilon = sqrt(a1[2]*a2[2]) #media geometrica
    r = dist * 0.1 ##transformar em nm
    C12 = (sigma/r)**12
    C6 = (sigma/r)**6
    LJ = 4*epsilon * (C12 - C6)
    return LJ ##kJ / mol

def calculaCoulomb (c1, c2, dist):
    #c1 e c2: informacoes dos atomos (precisa apenas da carga) --> aminoacids.rtp
    #dist: distancia entre os atomos, valor retornado da funcao
    ke = 138.935485 #kJ mol-1 nm e-2
    r = dist * 0.1 ##transformar em nm
    VC = (ke * c1 * c2) / r
    return VC ##kJ/mol
    
def pesquisaFfnonbonded(atomo, ffnonbonded):
    
    for i in range(len(ffnonbonded)):
        if ffnonbonded[i][0] == atomo:
            return copy.deepcopy(ffnonbonded[i])            
    
    return None

def pesquisaCargaAA(atm, res, aminoacids):
    #print aminoacids
    for i in range(len(aminoacids[res][0])):
        if aminoacids[res][0][i][0] == atm:
            return copy.deepcopy(aminoacids[res][0][i][2]), copy.deepcopy(aminoacids[res][0][i][1])
            
    for i in range(len(aminoacids[res][0])):
        if aminoacids[res][0][i][0] == atm[0]:
            return copy.deepcopy(aminoacids[res][0][i][2]), copy.deepcopy(aminoacids[res][0][i][1])
            
    return None
