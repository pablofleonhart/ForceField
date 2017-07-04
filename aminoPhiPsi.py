from pdbReader import PDBReader
import copy
import math
import matplotlib.pyplot as plt
import numpy as np

class AminoPhiPsi:
	pdbPattern = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s} {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}"
	filename = None
	ATOM_TAG = "ATOM"
	END_TAG = "TER"
	ALPHA_TAG = " CA "
	CARBON_TAG = " C  "
	NITROGEN_TAG = " N  "
	OC_ATOMS = ("C", "O", "OC", "HOC", "HC", "HO", "OXT")
	NH_ATOMS = ("N", "H", "1H", "H1", "2H", "H2", "3H", "H3")
	NHC_ATOMS = (" N  ", " H  ", "1H  ", " H1 ", "2H  ", " H2 ", "3H  ", " H3 ", " CA ")

	phi = []
	psi = []
	angles = []
	omega = []
	pdb = None

	def __init__( self, filename ):
		self.filename = filename
		self.readFile()
		self.calcAngles()

	def degrees( self, x ):
		return x*180/math.pi

	def rad( self, x ):
		return math.pi*x/180

	def readFile( self ):
		self.pdb = PDBReader( self.filename )

	def calcAngles( self ):
		file = open( "aminoPhiPsi.txt", "w" )
		file.write( "{:5s}  {:>7s}  {:>7s}  {:>7s}".format( "Amino", "Phi", "Psi", "Omega" ) + "\n" )
		self.phi = []
		self.psi = []
		self.angles = []
		self.omega = []

		for i in range ( 0, len( self.pdb.dicContent ) ):
			phiValue = 360.00
			psiValue = 360.00
			omegaValue = 360.00

			if i > 0:
				phiValue = self.calcDihedralAngle( self.pdb.dicContent.get( str( i-1 ) ).getPosC(), self.pdb.dicContent.get( str( i ) ).getPosN(), self.pdb.dicContent.get( str( i ) ).getPosCA(), self.pdb.dicContent.get( str( i ) ).getPosC() )
				omegaValue = self.calcDihedralAngle( self.pdb.dicContent.get( str( i-1 ) ).getPosCA(), self.pdb.dicContent.get( str( i-1 ) ).getPosC(), self.pdb.dicContent.get( str( i ) ).getPosN(), self.pdb.dicContent.get( str( i ) ).getPosCA() )
			
			if i < len( self.pdb.dicContent )-1:
				psiValue = self.calcDihedralAngle( self.pdb.dicContent.get( str( i ) ).getPosN(), self.pdb.dicContent.get( str( i ) ).getPosCA(), self.pdb.dicContent.get( str( i ) ).getPosC(), self.pdb.dicContent.get( str( i+1 ) ).getPosN() )				

			self.phi.append( math.radians( phiValue ) )
			self.psi.append( math.radians( psiValue ) )
			self.angles.append( math.radians( phiValue ) )
			self.angles.append( math.radians( psiValue ) )
			self.omega.append( math.radians( omegaValue ) )

			file.write( "{:5s}  {:7.2f}  {:7.2f}  {:7.2f}".format( self.pdb.dicContent.get( str( i ) ).getAminoAcid(), phiValue, psiValue, omegaValue ) + "\n" )

		file.close()

	def writeAngles( self ):
		file = open( "aminoPhiPsi.txt", "w" )
		for i in range ( 0, len( self.pdb.dicContent ) ):
			file.write( "{:5s}  {:7.2f}  {:7.2f}  {:7.2f}".format( self.pdb.dicContent.get( str( i ) ).getAminoAcid(), \
						self.phi[i], self.psi[i], self.omega[i] ) + "\n" )

		file.close()

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

	def adjustOmegas( self, angles = [] ):
		sizeAminoAcids = len( self.pdb.dicContent.keys() )

		if len( angles ) == 0:
			minIndex = min( self.pdb.aminoAcids )
			angles = [math.pi]*sizeAminoAcids
			for i in xrange( sizeAminoAcids ):
				if i + minIndex < max( self.pdb.aminoAcids ):
					carbonIndex = zip( self.pdb.atoms, self.pdb.aminoAcids ).index( ( self.CARBON_TAG,  i + minIndex ) )		
					nitrogenIndex = zip( self.pdb.atoms, self.pdb.aminoAcids ).index( ( self.NITROGEN_TAG, i + 1 + minIndex ) )
					domega = math.atan2( math.sin( angles[i] - self.omega[i+1] ), math.cos( angles[i] - self.omega[i+1] ) )
					carbonPos = self.pdb.posAtoms[carbonIndex]
					nitrogenPos = self.pdb.posAtoms[nitrogenIndex]
					index = 0

					for atom in zip( self.pdb.atoms, self.pdb.aminoAcids ):
						if ( atom[1] > i + 1 + minIndex or ( atom[1] == i + 1 + minIndex and ( atom[0] != self.NITROGEN_TAG ) ) ):
							self.pdb.posAtoms[index] = self.rotateAtomsBond( domega, self.pdb.posAtoms[index], carbonPos, nitrogenPos )
						index += 1

	def adjustPhiPsi( self, angles = [] ):
		size = max( self.pdb.aminoAcids )
		if len( angles ) == 0:
			angles = [[math.radians( 180.0 ) for x in range( 2 )] for y in range( size )]

		#print "########################"
		#print self.pdb.posAtoms[15]
		'''print "angles", angles
		print self.phi
		print self.psi'''
		sizeAminoAcids = len( self.pdb.dicContent.keys() )
		minIndex = min( self.pdb.aminoAcids )
		maxIndex = max( self.pdb.aminoAcids )
		#print minIndex, maxIndex

		for i in xrange( sizeAminoAcids ):
			'''if i + minIndex <= maxIndex:
				nitrogenIndex = zip( self.pdb.atoms, self.pdb.aminoAcids ).index( ( self.NITROGEN_TAG, i + minIndex ) )   
				alphaIndex = zip( self.pdb.atoms, self.pdb.aminoAcids ).index( ( self.ALPHA_TAG, i + minIndex ) )
				dphi = math.atan2( math.sin( angles[i][0] - self.phi[i] ), math.cos( angles[i][0] - self.phi[i] ) )
				#print "dphi", dphi
				nitrogenPos = self.pdb.posAtoms[nitrogenIndex]
				posCA = self.pdb.posAtoms[alphaIndex]
				#print "N", nitrogenPos, nitrogenIndex
				#print "CA", posCA, alphaIndex
				idx = 0
				for atom in zip( self.pdb.atoms, self.pdb.aminoAcids ):
					if ( i > 0 ) and ( atom[1] > i + minIndex or ( atom[1] == i + minIndex and ( atom[0].strip() not in self.NHC_ATOMS ) ) ):
						#print atom
						self.pdb.posAtoms[idx] = self.rotateAtomsBond( dphi, self.pdb.posAtoms[idx], nitrogenPos, posCA )
					idx += 1

				carbonIndex = zip( self.pdb.atoms, self.pdb.aminoAcids ).index( ( self.CARBON_TAG,  i + minIndex ) )  
				alphaIndex = zip( self.pdb.atoms, self.pdb.aminoAcids ).index( ( self.ALPHA_TAG, i + minIndex ) )
				dpsi = math.atan2( math.sin( angles[i][1] - self.psi[i] ), math.cos( angles[i][1] - self.psi[i] ) )              
				carbonPos = self.pdb.posAtoms[carbonIndex] 
				posCA = self.pdb.posAtoms[alphaIndex]
				idx = 0
				for atom in zip( self.pdb.atoms, self.pdb.aminoAcids ):
					if ( i + minIndex < maxIndex ) and ( atom[1] > i + minIndex or ( atom[1] == i + minIndex and ( atom[0] == " O  " ) ) ): 
						self.pdb.posAtoms[idx] = self.rotateAtomsBond( dpsi, self.pdb.posAtoms[idx], posCA, carbonPos )
					idx += 1

				#print self.pdb.posAtoms
				self.writePDBFile( "test" + str(i) + ".pdb" )'''

			if i > 0:
				pc_i  = zip( self.pdb.atoms, self.pdb.aminoAcids).index((" C  ",  i-1 + min(self.pdb.aminoAcids)))
				n_i   = zip( self.pdb.atoms, self.pdb.aminoAcids).index((" N  ",  i + min(self.pdb.aminoAcids)))   
				ca_i  = zip( self.pdb.atoms, self.pdb.aminoAcids).index((" CA ", i + min(self.pdb.aminoAcids)))
				c_i   = zip( self.pdb.atoms, self.pdb.aminoAcids).index((" C  ",  i + min(self.pdb.aminoAcids)))
				pc_pos = self.pdb.posAtoms[pc_i]
				n_pos  = self.pdb.posAtoms[n_i]
				ca_pos = self.pdb.posAtoms[ca_i]
				c_pos  = self.pdb.posAtoms[c_i] 
				current_angle = math.radians( self.calcDihedralAngle( pc_pos, n_pos, ca_pos, c_pos) )
				#print self.phi[i]
				dphi = self.angle_diff(self.phi[i], angles[i][0])
				ia = 0
				for atom in zip( self.pdb.atoms, self.pdb.aminoAcids):
					if (atom[1] > i + min(self.pdb.aminoAcids) or (atom[1] == i + min(self.pdb.aminoAcids) and (atom[0] not in self.NHC_ATOMS))):
						print atom[0]
						self.pdb.posAtoms[ia] = self.rotateAtomsBond(dphi, self.pdb.posAtoms[ia], n_pos, ca_pos)   
					ia += 1        

			if i + min(self.pdb.aminoAcids) < max(self.pdb.aminoAcids):
				n_i  = zip( self.pdb.atoms, self.pdb.aminoAcids).index((" N  ",  i + min(self.pdb.aminoAcids)))
				ca_i = zip( self.pdb.atoms, self.pdb.aminoAcids).index((" CA ", i + min(self.pdb.aminoAcids)))
				c_i  = zip( self.pdb.atoms, self.pdb.aminoAcids).index((" C  ",  i + min(self.pdb.aminoAcids)))
				nn_i = zip( self.pdb.atoms, self.pdb.aminoAcids).index((" N  ",  i+1 + min(self.pdb.aminoAcids)))  
				n_pos  = self.pdb.posAtoms[n_i]
				ca_pos = self.pdb.posAtoms[ca_i]
				c_pos  = self.pdb.posAtoms[c_i]
				nn_pos = self.pdb.posAtoms[nn_i] 
				current_angle = math.radians( self.calcDihedralAngle(n_pos, ca_pos, c_pos, nn_pos) )
				#print self.psi[i]
				dpsi = self.angle_diff(self.psi[i], angles[i][1])
				ia = 0
				for atom in zip( self.pdb.atoms, self.pdb.aminoAcids):
					if (atom[1] > i+min(self.pdb.aminoAcids) or (atom[1] == i+min(self.pdb.aminoAcids) and (atom[0]==" O  "))): 
						self.pdb.posAtoms[ia] = self.rotateAtomsBond(dpsi, self.pdb.posAtoms[ia], ca_pos, c_pos)         
					ia += 1

			self.writePDBFile( "test" + str(i) + ".pdb" )

	def angle_diff(self, angle_source, angle_target):
		a = angle_target - angle_source
		if a > math.pi:
		    a -= math.pi * 2.0
		if a < -math.pi:
		    a += math.pi * 2.0
		return a

	def normalize( self, v ):
		norm = np.linalg.norm( v )
		if norm == 0: 
			return v
		return v/norm  

	def rotateAtomsBond( self, theta, posAtom, bond_start, bond_end ):
		v = np.array( posAtom ) - np.array( bond_start )
		k = np.array( bond_end ) - np.array( bond_start )
		k = self.normalize( k )
		rot_pos = v * np.cos( theta ) + ( np.cross( k, v ) ) * np.sin( theta ) + k * ( np.dot( k, v ) ) * ( 1.0 - np.cos( theta ) )
		return list( rot_pos + np.array( bond_start ) )

	def getOmegas( self ):
		angles = []
		ca = self.pdb.getCAInfo()
		n  = self.pdb.getNInfo()
		c  = self.pdb.getCInfo()

		for index in xrange(len(ca)):
			if index > 0:
				name = ca[index-1][1]
				posCA = ca[index-1][2]
				carbonPos = c[index-1][2]                
				nex_nitrogenPos = n[index][2]
				nex_ca_pos = ca[index][2]
				omega = self.calcDihedralAngle( posCA, carbonPos, nex_nitrogenPos, nex_ca_pos )
				angles.append( omega )
			else:
				angles.append( 360.00 )

		self.omega = angles
		return angles

	def getPhiPsi( self ):
		angles = []   
		ca = self.pdb.getCAInfo()
		n  = self.pdb.getNInfo()
		c  = self.pdb.getCInfo()
		self.phi = []
		self.psi = []

		for index in xrange( len( ca ) ):
			name = ca[index][1]
			if index > 0:
				pre_c_pos = c[index-1][2]
			nitrogenPos = n[index][2]
			posCA = ca[index][2]
			carbonPos = c[index][2]
			if index < len(ca) - 1:
				nex_nitrogenPos = n[index+1][2]       
			if index > 0: 
				phiValue = self.calcDihedralAngle( pre_c_pos, nitrogenPos, posCA, carbonPos )
			else:
				phiValue = 360.00  
			if index < len(ca) - 1: 
				psiValue = self.calcDihedralAngle( nitrogenPos, posCA, carbonPos, nex_nitrogenPos )
			else:
				psiValue = 360.00

			self.phi.append( phiValue )
			self.psi.append( psiValue )

			angles.append( phiValue )
			angles.append( psiValue )

		return angles

	def calc_angle_3(self, pos1, posC, pos2):
		pos1 = np.array(pos1)
		posC = np.array(posC)
		pos2 = np.array(pos2)
		bond1C = self.normalize(pos1 - posC)
		bond2C = self.normalize(pos2 - posC)
		dp = np.dot(bond1C, bond2C)
		angle = np.arccos(dp)
		return angle  

	def get_peptide_bond_angles(self):
		angles = []
		ca = self.pdb.getCAInfo()
		n  = self.pdb.getNInfo()
		c  = self.pdb.getCInfo()
		for index in xrange( len( ca ) ):
			if index < len(ca) - 1:
				name = ca[index][1]
				carbonPos = c[index][2]                
				nex_nitrogenPos = n[index+1][2]
				nex_ca_pos = ca[index+1][2]
				alpha = self.calc_angle_3( carbonPos, nex_nitrogenPos, nex_ca_pos )
				angles.append( alpha )
			else:
				angles.append( 2.0 * math.pi )

		return angles 

	def adjustPeptideBonds( self, angles = [] ):
		sizeAminoAcids = len( self.pdb.aminoAcids )
		if len( angles ) == 0:
			angles = [math.radians( 120.0 )] * sizeAminoAcids

		for i in xrange( sizeAminoAcids ):
			if i + min( self.pdb.aminoAcids ) < max( self.pdb.aminoAcids ):
				carbonIndex   = zip(self.pdb.atoms, self.pdb.aminoAcids).index(( self.CARBON_TAG,  i + min(self.pdb.aminoAcids)))
				nitrogenIndex  = zip(self.pdb.atoms, self.pdb.aminoAcids).index(( self.NITROGEN_TAG,  i + 1 + min(self.pdb.aminoAcids)))
				nh_i  = -1.0
				for a in self.NH_ATOMS:
					if a != "N":                                                 
						nh_i  = zip(self.pdb.atoms, self.pdb.aminoAcids).index((a,  i + 1 + min(self.pdb.aminoAcids))) if (a,  i + 1 + min(self.pdb.aminoAcids)) in zip(self.pdb.atoms, self.pdb.aminoAcids) else -1
						if nh_i >= 0:
							break
					    
				nca_i = zip(self.pdb.atoms, self.pdb.aminoAcids).index(( self.ALPHA_TAG, i + 1 + min(self.pdb.aminoAcids)))
				current_alphas = self.get_peptide_bond_angles()
				dalpha = math.atan2(math.sin(angles[i] - current_alphas[i]), math.cos(angles[i] - current_alphas[i]))
				carbonPos   = self.pdb.posAtoms[carbonIndex]
				nitrogenPos  = self.pdb.posAtoms[nitrogenIndex]
				if nh_i >= 0:
					nh_pos  = self.pdb.posAtoms[nh_i]
					current_alphaH = self.calc_angle_3(carbonPos, nitrogenPos, nh_pos)
					dalphaH = math.atan2(math.sin(angles[i] - current_alphaH), math.cos(angles[i] - current_alphaH)) 
				nca_pos = self.pdb.posAtoms[nca_i]
				index = 0
				for atom in zip(self.pdb.atoms, self.pdb.aminoAcids):
					if (atom[1] > i + 1 + min(self.pdb.aminoAcids) or (atom[1] == i + 1 + min(self.pdb.aminoAcids) and (atom[0] not in self.NH_ATOMS))): 
						self.pdb.posAtoms[index] = self.bend_bonds(dalpha, self.pdb.posAtoms[index], carbonPos, nitrogenPos, nca_pos)
					elif nh_i >= 0 and atom[1] == i + 1 + min(self.pdb.aminoAcids) and atom[0] in self.NH_ATOMS and atom[0] != " N  ":
						self.pdb.posAtoms[index] = self.bend_bonds(dalphaH, self.pdb.posAtoms[index], carbonPos, nitrogenPos, nh_pos)    
					index += 1

	def bend_bonds( self, theta, posAtom, pos1, posC, pos2 ):
		posC = np.array( posC )
		posAtom = np.array( posAtom ) - posC
		pos1 = np.array( pos1 ) - posC
		pos2 = np.array( pos2 ) - posC
		bond1C = self.normalize( pos1 )
		bond2C = self.normalize( pos2 )
		ortho = self.normalize( np.cross( bond1C, bond2C ) )

		c = np.cos( theta )
		s = np.sin( theta )
		t = 1.0 - c

		rot = np.matrix( [[c + ortho[0] * ortho[0] * t, ortho[0] * ortho[1] * t - ortho[2] * s, ortho[0] * ortho[2] * t + ortho[1] * s], 
		                 [ortho[0] * ortho[1] * t + ortho[2] * s, c + ortho[1] * ortho[1] * t, ortho[1] * ortho[2] * t - ortho[0] * s],
		                 [ortho[2] * ortho[0] * t - ortho[1] * s, ortho[2] * ortho[1] * t + ortho[0] * s, c + ortho[2] * ortho[2] * t]] )
		                 
		new_pos = np.matrix.tolist( np.matrix( posAtom ) * rot.transpose() )[0]               
		new_pos = list( np.array( new_pos ) + posC )
		return new_pos

	def writePDBFile( self, name ):
		pdbNew = open( name, "w" )
		countTotal = 1
		acid = 0
		index = None
		for z in range( 0, len( self.pdb.atoms ) ):
			key = z

			if self.pdb.content[key].seqResidue != index:
				index = self.pdb.content[key].seqResidue
				acid += 1
			pdbNew.write( self.pdbPattern.format( self.pdb.content[key].tag, countTotal, self.pdb.content[key].atom, self.pdb.content[key].locIndicator, self.pdb.content[key].residue, self.pdb.content[key].chainID, \
						  acid, self.pdb.content[key].insResidue, float( self.pdb.posAtoms[z][0] ), float( self.pdb.posAtoms[z][1] ), float( self.pdb.posAtoms[z][2] ), \
						  float( self.pdb.content[key].occupancy ), float( self.pdb.content[key].temperature ), self.pdb.content[key].symbol, self.pdb.content[key].chargeAtom ) + "\n" )

			countTotal += 1

		pdbNew.write( "TER\n" )
		pdbNew.close()

	def plotRamanchandran( self ):
		plt.plot( self.phi, self.psi, 'ro', color = "green", ms = 7.0 )
		plt.xlim( -180, 180 )
		plt.ylim( -180, 180 )

		plt.xticks( np.arange( -180.1, 180.1, 60 ) )
		plt.yticks( np.arange( -180.1, 180.1, 60 ) )
		plt.xlabel( "Phi(deg)" )
		plt.ylabel( "Psi(deg)" )
		plt.arrow( -180, 0, 360, 0 )
		plt.arrow( 0, -180, 0, 360 )

		fig = plt.gcf()
		fig.set_size_inches( 6.0, 6.0 )
		fig.savefig( 'ramachandran.png', dpi=300, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None )