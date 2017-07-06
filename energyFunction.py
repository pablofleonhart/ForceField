import copy
from math import *
import random
import bonded as b
import nonbonded as nb

ffbonded = {}
ffnonbonded = []
aminoacids = {}
proteina = []
cutoff = 8.0

##########################
# Codigo do atomo pdb --> arquivos

def conversor(atom):
	global aminoacids
	
	key = atom[1]
	for i in range(len(aminoacids[key][0])):
		if aminoacids[key][0][i][0] == atom[0]:
			return aminoacids[key][0][i][1]

##########################

def estaLigado(atm1, atm2):
    global aminoacids
    
    ligado = False
    key = atm1[1]
    for i in range(len(aminoacids[key][1])):
        if (aminoacids[key][1][i][0] == atm1[0] and aminoacids[key][1][i][1] == atm2[0]) or (aminoacids[key][1][i][0] == atm2[0] and aminoacids[key][1][i][1] == atm1[0]):
            ligado = True
            break
    return ligado

##########################
# Retorna o potencial dihedral para 4 atomos

def estaDihedral(i, j, k, m):
	global proteina
	global ffbonded
	
	soma = 0.0
	
	while 1:
		if estaLigado(proteina[k], proteina[m]) == True: # Testa se os atomos k-m estao ligados (lig i - j - k existe) 
			a1 = conversor(proteina[i])
			a2 = conversor(proteina[j])
			a3 = conversor(proteina[k])
			a4 = conversor(proteina[m])
			die = b.procuraDihedral(a1, a2, a3, a4, ffbonded) # Procura a ligacao no arq ffbonded
			if die != None:
				Vda = b.dihedralAngles(proteina[i], proteina[j], proteina[k], proteina[m], die) # Calculo do Vda
				soma += Vda
			else:
				die = b.procuraDihedral('X', a2, a3, a4, ffbonded) # Procura lig X - j - k - m
				if die != None:
					Vda = b.dihedralAngles(proteina[i], proteina[j], proteina[k], proteina[m], die)
					soma += Vda
				else: 
					die = b.procuraDihedral('X', 'X', a3, a4, ffbonded) # Procura lig X - X - k - m
					if die != None:
						Vda = b.dihedralAngles(proteina[i], proteina[j], proteina[k], proteina[m], die)
						soma += Vda
					else:
						die = b.procuraDihedral('X', a2, a3, 'X', ffbonded) # Procura lig X - j - k - X
						if die != None:
							Vda = b.dihedralAngles(proteina[i], proteina[j], proteina[k], proteina[m], die)
							soma += Vda
		m += 1
		if m == (len(proteina)-1) or proteina[k][2] == proteina[m+1][2]:
			break
	return soma
	

##########################

def renameAtom( atom ):
	if atom == "HB3":
		return "HB1"

	elif atom == "2H":
		return "H2"

	return atom

def funEnergia(proteina):
	global ffbonded
	global ffnonbonded
	global aminoacids
		
	somatorio = 0.0
	
	for i in range(len(proteina)-1):
		for j in range(i+1, len(proteina)):
			if not(proteina[i][0] == "C" and proteina[j][0] == "N" and proteina[i][2] - proteina[j][2] == -1): # Nao lig.peptidica 
				if (((proteina[i][2] == proteina[j][2] and estaLigado(proteina[i], proteina[j]) == False)\
				 or (proteina[i][2] != proteina[j][2])) and b.distancia(proteina[i], proteina[j]) <= cutoff): # atomos no mesmo residuo nao ligados ou atomos em residuo diferente com d < cutoff
					dist = b.distancia(proteina[i], proteina[j])
					print dist
					print proteina[i][0], proteina[i][1]
					print proteina[j][0], proteina[j][1]
					atomi = renameAtom( proteina[i][0] )
					atomj = renameAtom( proteina[j][0] )

					c1, a1 = nb.pesquisaCargaAA( atomi, proteina[i][1], aminoacids) 
					c2, a2 = nb.pesquisaCargaAA( atomj, proteina[j][1], aminoacids)					
					atm1 = nb.pesquisaFfnonbonded(a1, ffnonbonded)
					atm2 = nb.pesquisaFfnonbonded(a2, ffnonbonded)
					
					Vlj = nb.calculaLennarJones(atm1, atm2, dist)
					Vc = nb.calculaCoulomb(c1, c2, dist) 
					somatorio += (Vlj + Vc)
				if (proteina[i][2] == proteina[j][2] and estaLigado(proteina[i], proteina[j]) == True): # atomos ligados no mesmo residuo
					dist = b.distancia(proteina[i], proteina[j])
					a1 = conversor(proteina[i])
					a2 = conversor(proteina[j])
					bond = b.procuraLig(a1, a2, ffbonded) # procura ligacao no ffbonded
					if bond != None:
						Vb = b.covalentBond(bond, dist) # calculo Vb
						somatorio += Vb
					k = j
					while 1:
						if estaLigado(proteina[j], proteina[k]) == True: # Testa ligacao i - j - k
							a3 = conversor(proteina[k])
							ang = b.procuraAngle(a1, a2, a3, ffbonded) # Procura ligacao i - j - k
							if ang != None:
								Va = b.bondAngles(proteina[i], proteina[j], proteina[k], ang) # Calcula Va
								somatorio += Va
							Vda = estaDihedral(i, j, k, k+1) # Chama funcao para calcular Vda
							somatorio += Vda
						k += 1
						#print k, len( proteina )
						if k >= (len(proteina)-1) or proteina[j][2] == proteina[k][2]:
							break
			else: # ligacao peptidica
				dist = b.distancia(proteina[i], proteina[j])
				bond = b.procuraLig(proteina[i][0], proteina[j][0], ffbonded) # procura ligacao no ffbonded
				if bond != None:
					Vb = b.covalentBond(bond, dist) # calculo Vb
					somatorio += Vb
				for l in range(0,10):
					if proteina[l][0] == 'CA' or proteina[l][0] == 'H':
						if estaLigado(proteina[j], proteina[l]) == True: # Testa ligacao i - j - k
							a3 = conversor(proteina[l])
							ang = b.procuraAngle(proteina[i][0], proteina[j][0], a3, ffbonded) # Procura ligacao i - j - k
							if ang != None:
								Va = b.bondAngles(proteina[i], proteina[j], proteina[l], ang) # Calcula Va
								somatorio += Va
						Vda = Vda = estaDihedral(i, j, l, l+1) # Chama funcao para calcular Vda
						somatorio += Vda
	return somatorio
	
##########################
# Leitura dos arquivos: reference.pdb, ffbonded.itp, ffnonbonded.itp, aminoacids.rtp

def files():
	global proteina
	global ffbonded
	global ffnonbonded
	global aminoacids
#-------------------------

	arq = open('reference.pdb', 'r')
	texto = arq.readlines()
	arq.close()
    
	linhas = []
	for i in range(len(texto)):
		linhas.append(texto[i].split())

	for i in range(len(linhas)):
		if linhas[i][0] == 'ATOM':
			if linhas[i][5] == '1':
				proteina.append([copy.deepcopy(linhas[i][2]), "N"+copy.deepcopy(linhas[i][3]), int(linhas[i][4]),float(linhas[i][5]), float(linhas[i][6]), float(linhas[i][7])])
			else:
				if linhas[i][5] == linhas[len(linhas)-2][5] :
					proteina.append([copy.deepcopy(linhas[i][2]), "C"+copy.deepcopy(linhas[i][3]), int(linhas[i][4]),float(linhas[i][5]), float(linhas[i][6]), float(linhas[i][7])])
				else:
					proteina.append([copy.deepcopy(linhas[i][2]), copy.deepcopy(linhas[i][3]), int(linhas[i][4]),float(linhas[i][5]), float(linhas[i][6]), float(linhas[i][7])])
	#print proteina               
#------------------------

	arq = open('AMBER99/ffbonded.itp','r')
	texto = arq.readlines()
	arq.close()
	
	linhas = []
	for i in range(len(texto)):
		texto[i] = texto[i].split()
		if len(texto[i]) != 0:
			linhas.append(texto[i])	
	
	for i in range(len(linhas)):
		if linhas[i][0] == '[':
			key = linhas[i][1]
			if key == 'bondtypes':
				ffbonded[key] = [] # [atom1, atom2, b0, kb]
				while 1:
					i += 1
					if linhas[i][0] != ';': 
						if linhas[i][0] == '[': 
							break
						ffbonded[key].append([linhas[i][0],linhas[i][1], float(linhas[i][3]), float(linhas[i][4])])
			if key == 'angletypes':
				ffbonded[key] = [] # [atom1, atom2, atom3, theta, kt]
				while 1:
					i += 1
					if linhas[i][0] != ';': 
						if linhas[i][0] == '[': 
							break
						ffbonded[key].append([linhas[i][0],linhas[i][1],linhas[i][2], float(linhas[i][4]), float(linhas[i][5])])
			if key == 'dihedraltypes':# [a1, a2, a3, a4, phase, kd, n]
				if ffbonded.get( key ) is None:
					ffbonded[key] = []
				while 1:
					i += 1
					if i != len(linhas)-1:
						if linhas[i][0] != ';i': 
							if linhas[i][0] == '[': 
								break
							ffbonded[key].append([linhas[i][0],linhas[i][1],linhas[i][2],linhas[i][3], float(linhas[i][5]), float(linhas[i][6]), float(linhas[i][7])])
					else:
						ffbonded[key].append([linhas[i][0],linhas[i][1],linhas[i][2],linhas[i][3], float(linhas[i][5]), float(linhas[i][6]), float(linhas[i][7])])
						break

	#print ffbonded
	
#-------------------------

	arq = open('AMBER99/ffnonbonded.itp','r')
	texto = arq.readlines()
	arq.close()
	
	linhas = []
	for i in range(len(texto)):
		linhas.append(texto[i].split())
	
	for i in range(len(linhas)):
		if linhas[i][0][0] != '[' and linhas[i][0][0] != ';' :
			ffnonbonded.append([copy.deepcopy(linhas[i][0]), float(linhas[i][5]), float(linhas[i][6])]) #([Nome, Sigma, Epsilon])

	#print ffnonbonded

#-------------------------
	
	arq = open('AMBER99/aminoacids.rtp','r')
	texto = arq.readlines()
	arq.close()
	
	indice = 0
	while True:
		if texto[indice] == "; Next are non-terminal AA's\n":
			break
		else:
			indice += 1
	indice += 2
	
	atom = False
	bond = False
	
	key = ""
	
	for i in range(indice, len(texto)):
		if texto[i].split() != []:
			if texto[i][0] == '[':
				key = texto[i].split()[1]
				aminoacids[key] = [[],[]] # [atoms, bonds]
			elif texto[i].split()[1] == "atoms":
				atom = True
				bond = False
			elif texto[i].split()[1] == "bonds":
				atom = False
				bond = True
			elif texto[i].split()[1] == "impropers":
				atom = False
				bond = False
			elif atom == True:
				aminoacids[key][0].append([texto[i].split()[0], texto[i].split()[1], float(texto[i].split()[2])]) # atoms
			elif bond == True:
				aminoacids[key][1].append([texto[i].split()[0], texto[i].split()[1]]) # bonds

	#print aminoacids
##########################

def main():
	
	files()
	print '\n'
	print '-----------------------'
	print '\n'
	print 'Funcao de Energia: ' + str(funEnergia(proteina))		
	print '\n'
	print '-----------------------'
	print '\n'	
#-------------------------

if __name__ == '__main__':
    main()