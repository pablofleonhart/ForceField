import copy
import math

class EnergyFunction( object ):
	aminoFile 	  = "amber99/aminoacids.rtp"
	bondedFile 	  = "amber99/ffbonded.itp"
	nonbondedFile = "amber99/ffnonbonded.itp"

	cutoff = 8.0

	def __init__( self, pdb ):
		self.protein = []
		self.aminoacids = {}
		self.ffbonded = {}
		self.ffnonbonded = []
		self.pdb = pdb

		self.readProtein()
		self.readAmberFiles()

	def areBonded( self, aminoA, atomA, atomB ):
		bonded = False
		key = aminoA
		for i in range( len( self.aminoacids[key][1] ) ):
		    if ( self.aminoacids[key][1][i][0] == atomA and self.aminoacids[key][1][i][1] == atomB )\
		    	or ( self.aminoacids[key][1][i][0] == atomB and self.aminoacids[key][1][i][1] == atomA ):
		        bonded = True
		        break

		return bonded

	def distance( self, atomA, atomB ):
	    x = ( atomA[0] - atomB[0] )**2
	    y = ( atomA[1] - atomB[1] )**2
	    z = ( atomA[2] - atomB[2] )**2 
	    distance = math.sqrt( x + y + z )
	    
	    return distance

	def findCharge( self, atom, aminoAcid ):
		#print self.aminoacids[aminoAcid]
		
		for i in range( len( self.aminoacids[aminoAcid][0] ) ):
			if self.aminoacids[aminoAcid][0][i][0] == atom:
				return copy.deepcopy( self.aminoacids[aminoAcid][0][i][2] ), copy.deepcopy( self.aminoacids[aminoAcid][0][i][1] )
	            
		for i in range( len( self.aminoacids[aminoAcid][0] ) ):
			if self.aminoacids[aminoAcid][0][i][0] == atom[0]:
				return copy.deepcopy( self.aminoacids[aminoAcid][0][i][2] ), copy.deepcopy( self.aminoacids[aminoAcid][0][i][1] )
	            
		return None

	def renameAtom( self, atom ):
		if atom == "HB3":
			return "HB1"

		elif atom == "2H":
			return "H2"

		return atom

	def findNonbonded( self, atom ):
		for i in range( len( self.ffnonbonded ) ):
			if self.ffnonbonded[i][0] == atom:
				return copy.deepcopy( self.ffnonbonded[i] )

		return None

	def conversor( self, atom, aminoAcid ):
		key = aminoAcid
		for i in range( len( self.aminoacids[key][0] ) ):
			if self.aminoacids[key][0][i][0] == atom:
				return self.aminoacids[key][0][i][1]

	def findDihedral( self, atom1, atom2, atom3, atom4 ):
		for i in range( len( self.ffbonded['dihedraltypes'] ) ):
			if self.ffbonded['dihedraltypes'][i][0] == atom1 and self.ffbonded['dihedraltypes'][i][1] == atom2\
				and self.ffbonded['dihedraltypes'][i][2] == atom3 and self.ffbonded['dihedraltypes'][i][3] == atom4:
				return copy.deepcopy( self.ffbonded['dihedraltypes'][i] )

	def findAngle( self, atom1, atom2, atom3 ):
		for i in range( len( self.ffbonded['angletypes'] ) ):
			if self.ffbonded['angletypes'][i][0] == atom1 and self.ffbonded['angletypes'][i][1] == atom2 and self.ffbonded['angletypes'][i][2] == atom3:
				return copy.deepcopy( self.ffbonded['angletypes'][i] )

	def findBond( self, atom1, atom2 ):
		for i in range( len( self.ffbonded['bondtypes'] ) ):
			if ( self.ffbonded['bondtypes'][i][0] == atom1 and self.ffbonded['bondtypes'][i][1] == atom2 ) or ( self.ffbonded['bondtypes'][i][0] == atom2 and self.ffbonded['bondtypes'][i][1] == atom1 ):
				return copy.deepcopy( self.ffbonded['bondtypes'][i] )

	def covalentBond( self, bond, dist ):
		r = dist * 0.1
		value = bond[3] * ( r - bond[2] )**2
		return value ## kJ / mol

	def bondAngles( self, atom1, atom2, atom3, b ):
		#lig a1 - a2 - a3
		#b: informacao da ligacao [atom1, atom2, atom3, theta, kt]

		A = self.distance( atom1, atom2 )
		B = self.distance( atom2, atom3 )
		C = self.distance( atom1, atom3 )

		aux = ( A**2 + B**2 - C**2 ) / ( 2*A*B )
		theta = math.acos( aux ) # Calcula angulo de ligacao 
		theta = math.degrees( theta )

		value = b[4] * ( theta - b[3] )**2
		return value ## kJ / mol

	def calcDihedralAngle( self, atom1, atom2, atom3, atom4 ):
		vector1 = [( atom2[0] - atom1[0] ), ( atom2[1] - atom1[1] ), ( atom2[2] - atom1[2] )]
		vector2 = [( atom3[0] - atom2[0] ), ( atom3[1] - atom2[1] ), ( atom3[2] - atom2[2] )]
		vector3 = [( atom4[0] - atom3[0] ), ( atom4[1] - atom3[1] ), ( atom4[2] - atom3[2] )]

		normal1 = [( vector1[1] * vector2[2] - vector1[2] * vector2[1] ), ( vector1[2] * vector2[0] - vector1[0] * vector2[2] ), ( vector1[0] * vector2[1] - vector1[1] * vector2[0] )]
		normal2 = [( vector2[1] * vector3[2] - vector2[2] * vector3[1] ), ( vector2[2] * vector3[0] - vector2[0] * vector3[2] ), ( vector2[0] * vector3[1] - vector2[1] * vector3[0] )]

		n1 = math.sqrt( ( normal1[0]**2) + ( normal1[1]**2 ) + ( normal1[2]**2 ) )
		n2 = math.sqrt( ( normal2[0]**2) + ( normal2[1]**2 ) + ( normal2[2]**2 ) )

		normal1 = [( normal1[0]/n1 ), ( normal1[1]/n1 ), ( normal1[2]/n1 )]
		normal2 = [( normal2[0]/n2 ), ( normal2[1]/n2 ), ( normal2[2]/n2 )]

		v2 = math.sqrt( ( vector2[0]**2 ) + ( vector2[1]**2 ) + ( vector2[2]**2 ) )
		vector2 = [( vector2[0]/v2 ), ( vector2[1]/v2 ), ( vector2[2]/v2 )]

		m1 = [( ( normal1[1] * normal2[2] ) - ( normal1[2] * normal2[1] ) ), ( ( normal1[2] * normal2[0] ) - ( normal1[0] * normal2[2] ) ), ( ( normal1[0] * normal2[1] ) - ( normal1[1] * normal2[0] ) )]

		x = ( normal1[0] * normal2[0] ) + ( normal1[1] * normal2[1] ) + ( normal1[2] * normal2[2] )
		y = ( m1[0] * vector2[0] ) + ( m1[1] * vector2[1] ) + ( m1[2] * vector2[2] )

		angle = round( math.degrees( math.atan2( y, x ) ), 2 )
		return angle

	def dihedralAngles( self, at1, at2, at3, at4, b ):
		a1 = [at1[3], at1[4], at1[5]]
		a2 = [at2[3], at2[4], at2[5]]
		a3 = [at3[3], at3[4], at3[5]]
		a4 = [at4[3], at4[4], at4[5]]
		angle = self.calcDihedralAngle( a1, a2, a3, a4 )
		phase = math.radians(b[4])
		
		Vda = ( b[5] / 2 ) * ( 1 + math.cos( ( b[6] * angle ) - phase ) )
		
		return Vda ## kJ / mol 

	def areDihedral( self, i, j, k, m ):
		total = 0.0
		
		atomI = self.protein[i][0]
		atomJ = self.protein[j][0]
		aminoI = self.protein[i][1]
		aminoJ = self.protein[j][1]
		seqAminoI = self.protein[i][2]
		seqAminoJ = self.protein[j][2]
		posAtomI = self.protein[i][3:6] 
		posAtomJ = self.protein[j][3:6]
		atomK = self.protein[k][0]
		aminoK = self.protein[k][1]
		seqAminoK = self.protein[k][2]
		posAtomK = self.protein[k][3:6] 
		atomM = self.protein[m][0]
		aminoM = self.protein[m][1]
		seqAminoM = self.protein[m][2]
		posAtomM = self.protein[m][3:6] 

		while 1:
			if self.areBonded( aminoK, atomK, atomM ) == True:
				a1 = self.conversor( atomI, aminoI )
				a2 = self.conversor( atomJ, aminoJ )
				a3 = self.conversor( atomK, aminoK )
				a4 = self.conversor( self.protein[m][0], self.protein[m][1] )
				die = self.findDihedral( a1, a2, a3, a4 )

				if die != None:
					Vda = self.dihedralAngles( self.protein[i], self.protein[j], self.protein[k], self.protein[m], die )
					total += Vda
				else:
					die = self.findDihedral( 'X', a2, a3, a4 )
					if die != None:
						Vda = self.dihedralAngles( self.protein[i], self.protein[j], self.protein[k], self.protein[m], die )
						total += Vda
					else: 
						die = self.findDihedral( 'X', 'X', a3, a4, ffbonded )
						if die != None:
							Vda = self.dihedralAngles( self.protein[i], self.protein[j], self.protein[k], self.protein[m], die )
							total += Vda
						else:
							die = self.findDihedral( 'X', a2, a3, 'X', ffbonded )
							if die != None:
								Vda = self.dihedralAngles( self.protein[i], self.protein[j], self.protein[k], self.protein[m], die )
								total += Vda
			m += 1
			if m == ( len( self.protein ) - 1 ) or self.protein[k][2] == self.protein[m+1][2]:
				break

		return total		

	def calcLennardJones( self, atomA, atomB, distance ):
	    sigma = ( atomA[1] + atomB[1] ) / 2 #media aritmetica
	    epsilon = math.sqrt( atomA[2]*atomB[2] ) #media geometrica
	    r = distance * 0.1 ##transformar em nm
	    C12 = ( sigma / r )**12
	    C6 = ( sigma / r )**6
	    LJ = 4*epsilon * ( C12 - C6 )
	    return LJ ##kJ / mol

	def calcCoulumb( self, chargeA, chargeB, distance ):
	    ke = 138.935485 #kJ mol-1 nm e-2
	    r = distance * 0.1 ##transformar em nm
	    VC = ( ke * chargeA * chargeB ) / r
	    return VC ##kJ/mol

	def getEnergy( self ):
		total = 0
		for i in range( len( self.protein ) - 1 ):
			for j in range( i + 1, len( self.protein ) ):
				#print i, j
				atomI = self.protein[i][0]
				atomJ = self.protein[j][0]
				aminoI = self.protein[i][1]
				aminoJ = self.protein[j][1]
				seqAminoI = self.protein[i][2]
				seqAminoJ = self.protein[j][2]
				posAtomI = self.protein[i][3:6] 
				posAtomJ = self.protein[j][3:6]
				dist = self.distance( posAtomI, posAtomJ )

				if not( atomI == "C" and atomJ == "N" and seqAminoI - seqAminoJ == -1 ):
					if ( ( ( seqAminoI == seqAminoJ and self.areBonded( aminoI, atomI, atomJ ) == False ) or ( seqAminoI != seqAminoJ ) ) and dist <= self.cutoff ):
						charge1, a1 = self.findCharge( atomI, aminoI ) 
						#print charge1, a1, aminoI
						charge2, a2 = self.findCharge( atomJ, aminoJ )						
						#print charge2, a2, aminoJ
						#print self.areBonded( aminoI, atomI, atomJ ), atomI, "(", aminoI, ")", atomJ, "(", aminoJ, ")"
					
						atomI = self.renameAtom( atomI )
						atomj = self.renameAtom( atomJ )

						atom1 = self.findNonbonded( a1 )
						atom2 = self.findNonbonded( a2 )

						ljValue = self.calcLennardJones( atom1, atom2, dist )
						cValue = self.calcCoulumb( charge1, charge2, dist )

						total += ljValue + cValue

					if ( seqAminoI == seqAminoJ and self.areBonded( aminoI, atomI, atomJ ) == True ):
						a1 = self.conversor( atomI, aminoI )
						a2 = self.conversor( atomJ, aminoJ )

						bond = self.findBond( a1, a2 )

						if bond != None:
							bondValue = self.covalentBond( bond, dist )
							#total += bondValue

						k = j
						atomK = self.protein[k][0]
						aminoK = self.protein[k][1]
						posAtomK = self.protein[k][3:6]

						while 1:
							if self.areBonded( aminoJ, atomJ, atomK ) == True:
								a3 = self.conversor( atomK, aminoK )
								ang = self.findAngle( a1, a2, a3 )
								if ang != None:
									va = self.bondAngles( posAtomI, posAtomJ, posAtomK, ang )
									#total += va

								vda = self.areDihedral( i, j, k, k+1 )
								#total += vda

							k += 1
							#print k, len( self.protein )
							if k >= ( len( self.protein ) - 1 ) or self.protein[j][2] == self.protein[k][2]:
								break
				else:
					bond = self.findBond( atomI, atomJ )
					if bond != None:
						Vb = self.covalentBond( bond, dist )
						#total += Vb
					for l in range( 0, 10 ):
						if self.protein[l][0] == 'CA' or self.protein[l][0] == 'H':
							if self.areBonded( aminoJ, atomJ, self.protein[l][0] ) == True:
								a3 = self.conversor( self.protein[l][0], self.protein[l][1] )
								ang = self.findAngle( self.protein[i][0], self.protein[j][0], a3 )
								if ang != None:
									Va = self.bondAngles( posAtomI, posAtomJ, self.protein[l][3:6], ang )
									#total += Va
							Vda = Vda = self.areDihedral( i, j, l, l+1 )
							#total += Vda

		return total

	def readProtein( self ):
		'''file = open( self.filename, 'r' )
		content = file.readlines()
		file.close()
	    
		lines = []
		for i in range( len( content ) ):
			lines.append( content[i].split() )

		for i in range( len( lines ) ):
			if lines[i][0] == 'ATOM':
				if lines[i][5] == '1':
					self.protein.append( [copy.deepcopy( lines[i][2] ), "N" + copy.deepcopy( lines[i][3] ), int( lines[i][4] ),float( lines[i][5] ), float( lines[i][6] ), float( lines[i][7] )] )
				else:
					if lines[i][5] == lines[len( lines )-2][5] :
						self.protein.append( [copy.deepcopy( lines[i][2]), "C" + copy.deepcopy( lines[i][3] ), int( lines[i][4] ),float( lines[i][5] ), float( lines[i][6] ), float( lines[i][7] )] )
					else:
						self.protein.append( [copy.deepcopy( lines[i][2]), copy.deepcopy( lines[i][3] ), int( lines[i][4] ),float( lines[i][5] ), float( lines[i][6] ), float( lines[i][7] )] )
		
		print self.protein'''
		#print len( self.pdb.atoms )
		#print self.pdb.posAtoms
		for i in range( len( self.pdb.posAtoms ) ):
			if self.pdb.content[i].getTag() == 'ATOM':
				print i, self.pdb.content[i].seqResidue
				if self.pdb.content[i].seqResidue == 1:					
					print copy.deepcopy( self.pdb.content[i].getAtom() )
					self.protein.append( [copy.deepcopy( self.pdb.content[i].getAtom() ), "N" + copy.deepcopy( self.pdb.content[i].getResidue() ), self.pdb.content[i].seqResidue, float( self.pdb.posAtoms[i][0] ), float( self.pdb.posAtoms[i][1] ), float( self.pdb.posAtoms[i][2] )] )
				else:
					if self.pdb.content[i].seqResidue == self.pdb.content[len( self.pdb.content )-2].seqResidue:
						self.protein.append( [copy.deepcopy( self.pdb.content[i].getAtom() ), "C" + copy.deepcopy( self.pdb.content[i].getResidue() ), self.pdb.content[i].seqResidue, float( self.pdb.posAtoms[i][0] ), float( self.pdb.posAtoms[i][1] ), float( self.pdb.posAtoms[i][2] )] )
					else:
						self.protein.append( [copy.deepcopy( self.pdb.content[i].getAtom() ), copy.deepcopy( self.pdb.content[i].getResidue() ), self.pdb.content[i].seqResidue, float( self.pdb.posAtoms[i][0] ), float( self.pdb.posAtoms[i][1] ), float( self.pdb.posAtoms[i][2] )] )

		print self.protein

	def readAmberFiles( self ):
		# aminoAcid:
		file = open( self.aminoFile, 'r' )
		content = file.readlines()
		file.close()
		
		index = 0
		while True:
			if content[index] == "; Next are non-terminal AA's\n":
				break
			else:
				index += 1
		index += 2
		
		atom = False
		bond = False
		
		key = ""
		
		for i in range( index, len( content ) ):
			if content[i].split() != []:
				if content[i][0] == '[':
					key = content[i].split()[1]
					self.aminoacids[key] = [[],[]] # [atoms, bonds]
				elif content[i].split()[1] == "atoms":
					atom = True
					bond = False
				elif content[i].split()[1] == "bonds":
					atom = False
					bond = True
				elif content[i].split()[1] == "impropers":
					atom = False
					bond = False
				elif atom == True:
					self.aminoacids[key][0].append( [content[i].split()[0], content[i].split()[1], float( content[i].split()[2] )] ) # atoms
				elif bond == True:
					self.aminoacids[key][1].append( [content[i].split()[0], content[i].split()[1]] ) # bonds

		#print self.aminoacids

		# bondend:
		file = open( self.bondedFile, 'r' )
		content = file.readlines()
		file.close()
		
		lines = []
		for i in range(len(content)):
			content[i] = content[i].split()
			if len(content[i]) != 0:
				lines.append(content[i])	
		
		for i in range( len( lines ) ):
			if lines[i][0] == '[':
				key = lines[i][1]
				if key == 'bondtypes':
					self.ffbonded[key] = []
					while 1:
						i += 1
						if lines[i][0] != ';': 
							if lines[i][0] == '[': 
								break
							self.ffbonded[key].append( [lines[i][0], lines[i][1], float( lines[i][3] ), float( lines[i][4] )] )
				
				if key == 'angletypes':
					self.ffbonded[key] = []
					while 1:
						i += 1
						if lines[i][0] != ';': 
							if lines[i][0] == '[': 
								break
							self.ffbonded[key].append( [lines[i][0], lines[i][1], lines[i][2], float( lines[i][4] ), float( lines[i][5] )] )
				
				if key == 'dihedraltypes':
					if self.ffbonded.get( key ) is None:
						self.ffbonded[key] = []
					while 1:
						i += 1
						if i != len( lines )-1:
							if lines[i][0] != ';i': 
								if lines[i][0] == '[': 
									break
								self.ffbonded[key].append( [lines[i][0], lines[i][1], lines[i][2], lines[i][3], float( lines[i][5] ), float( lines[i][6] ), float( lines[i][7] )] )
						else:
							self.ffbonded[key].append( [lines[i][0], lines[i][1], lines[i][2], lines[i][3], float( lines[i][5] ), float( lines[i][6] ), float( lines[i][7] )] )
							break

		#print self.ffbonded

		# nonbonded:
		file = open( self.nonbondedFile, 'r' )
		content = file.readlines()
		file.close()
		
		lines = []
		for i in range( len( content ) ):
			lines.append( content[i].split() )

		for i in range( len( lines ) ):
			if lines[i][0][0] != '[' and lines[i][0][0] != ';' :
				self.ffnonbonded.append( [copy.deepcopy( lines[i][0] ), float( lines[i][5] ), float( lines[i][6] )] )

		#print self.ffnonbonded