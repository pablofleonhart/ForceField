import copy
from math import *
import random
import calculos as c

def distancia(atm1, atm2):
    #a1 e a2: informacoes dos atomos ['atom', 'aa', pos, x, y, z]
    x = (atm1[3] - atm2[3])**2
    y = (atm1[4] - atm2[4])**2
    z = (atm1[5] - atm2[5])**2 
    dist = sqrt(x+y+z)
    
    return dist ## angstrom 

def procuraDihedral(a1, a2, a3, a4, ffbonded):
	# Procura a ligacao a1 - a2 - a3 - a4
	
	for i in range(len(ffbonded['dihedraltypes'])):
		if ffbonded['dihedraltypes'][i][0] == a1 and ffbonded['dihedraltypes'][i][1] == a2 and ffbonded['dihedraltypes'][i][2] == a3 and ffbonded['dihedraltypes'][i][3] == a4:
			return copy.deepcopy(ffbonded['dihedraltypes'][i])

def procuraAngle(a1, a2, a3, ffbonded):
	# Procura a ligacao a1 - a2 - a3
	
	for i in range(len(ffbonded['angletypes'])):
		if ffbonded['angletypes'][i][0] == a1 and ffbonded['angletypes'][i][1] == a2 and ffbonded['angletypes'][i][2] == a3:
			return copy.deepcopy(ffbonded['angletypes'][i])

def procuraLig(a1, a2, ffbonded):
	# Procura a ligacao a1 - a2
	
	for i in range(len(ffbonded['bondtypes'])):
		if (ffbonded['bondtypes'][i][0] == a1 and ffbonded['bondtypes'][i][1] == a2) or (ffbonded['bondtypes'][i][0] == a2 and ffbonded['bondtypes'][i][1] == a1):
			return copy.deepcopy(ffbonded['bondtypes'][i])			

def dihedralAngles(at1, at2, at3, at4, b):
	#lig a1 - a2 - a3 - a4
	#b: informacao da ligacao [a1, a2, a3, a4, phase, kd, n]
	
	a1 = [at1[3], at1[4], at1[5]]
	a2 = [at2[3], at2[4], at2[5]]
	a3 = [at3[3], at3[4], at3[5]]
	a4 = [at4[3], at4[4], at4[5]]
	angle = c.calcPhi(a1, a2, a3, a4) # Calcula angulo dihedral para a1 - a2 - a3 - a4
	phase = radians(b[4])
	
	Vda = (b[5] / 2) * (1 + cos((b[6] * angle) - phase))
	
	return Vda ## kJ / mol 

def bondAngles(a1, a2, a3, b):
	#lig a1 - a2 - a3
	#b: informacao da ligacao [atom1, atom2, atom3, theta, kt]
	
	A = distancia(a1, a2)
	B = distancia(a2, a3)
	C = distancia(a1, a3)
	
	aux = (A**2 + B**2 - C**2) / (2*A*B)
	theta = acos(aux) # Calcula angulo de ligacao 
	theta = degrees(theta)
	
	Va = b[4] * (theta - b[3])**2
	return Va ## kJ / mol

def covalentBond(b, dist):
	#b: informacao da ligacao ['atom1', 'atom2', b0, kb]
	
	r = dist * 0.1
	Vb = b[3] * (r - b[2])**2
	return Vb ## kJ / mol
