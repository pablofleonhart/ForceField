from atom import Atom
import copy

class AminoAcids( object ):

	pdbPattern = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
	ATOM_TAG = "ATOM"
	END_TAG = "TER"
	ATOM_REF = "CA"
	OXIGEN_CARBOXYL = "OC"
	HYDROGEN_CARBOXYL = ( "HC", "HOC", "HO" )
	HYDROGEN_AMINO = ( "H", "1H" )
	NITROGEN = "N"

	dic = { "A": "files/alanine.pdb",
			"R": "files/arginine.pdb",
			"N": "files/asparagine.pdb",
			"D": "files/aspartic_acid.pdb",
			"C": "files/cysteine.pdb",
			"E": "files/glutamic_acid.pdb",
			"Q": "files/glutamine.pdb",
			"G": "files/glycine.pdb",
			"H": "files/histidine.pdb",
			"I": "files/isoleucine.pdb",
			"L": "files/leucine.pdb",
			"K": "files/lysine.pdb",
			"M": "files/methionine.pdb",
			"F": "files/phenalalanine.pdb",
			"P": "files/proline.pdb",
			"S": "files/serine.pdb",
			"T": "files/threonine.pdb",
			"W": "files/tryptophan.pdb",
			"Y": "files/tyrosine.pdb",
			"V": "files/valine.pdb" }

	asnNames = { "H" : " H2 ",
				 "1HD2" : "HD21",
				 "2HD2" : "HD22",
				 "2HB" : " HB2",
				 "1HB" : " HB3",
				 "2H" : " H2 ",
				 "1H" : " H1 "}

	leuNames = { "2H" : " H  ",
				 "2HB" : " HB2",
				 "1HB" : " HB3",
				 "1HD1" : "HD11",
				 "2HD1" : "HD12",
				 "3HD1" : "HD13",
				 "1HD2" : "HD21",
				 "2HD2" : "HD22",
				 "3HD2" : "HD23" }

	tyrNames = { "2H" : " H  ",
				 "2HB" : " HB2",
				 "1HB" : " HB3" }

	ileNames = { "2H" : " H  ",
				 "1HG1" : "HG13",
				 "2HG1" : "HG12",
				 "1HG2" : "HG21",
				 "2HG2" : "HG22",
				 "3HG2" : "HG23",
				 "1HD1" : "HD11",
				 "2HD1" : "HD12",
				 "3HD1" : "HD13" }

	glnNames = { "2H" : " H  ",
			     "2HB" : " HB2",
				 "1HB" : " HB3",
				 "2HG" : " HG2",
				 "1HG" : " HG3",
				 "1HE2" : "HE21",
				 "2HE2" : "HE22" }

	trpNames = { "2H" : " H  ",
			     "2HB" : " HB2",
				 "1HB" : " HB3" }

	lysNames = { "2HB" : " HB2",
				 "1HB" : " HB3",
				 "2HG" : " HG2",
				 "1HG" : " HG3",
				 "2HD" : " HD2",
				 "1HD" : " HD3",
				 "2HE" : " HE2",
				 "1HE" : " HE3",
				 "1HZ" : " HZ1",
				 "2HZ" : " HZ2",
				 "3HZ" : " HZ3" }

	aspNames = { "2H" : " H  ",
			     "2HB" : " HB2",
				 "1HB" : " HB3" }

	glyNames = { "2H" : " H  ",
			     "2HA" : " HA2",
				 "1HA" : " HA3" }

	proNames = { "2HB" : " HB2",
				 "1HB" : " HB3",
				 "2HG" : " HG2",
				 "1HG" : " HG3",
				 "2HD" : " HD2",
				 "1HD" : " HD3" }

	serNames = { "2H" : " H  ",
			     "2HB" : " HB2",
				 "1HB" : " HB3",
				 "OC" : " OXT" }

	argNames = { "2H" : " H  ",
				 "2HB" : " HB2",
				 "1HB" : " HB3",
				 "2HG" : " HG2",
				 "1HG" : " HG3",
				 "2HD" : " HD2",
				 "1HD" : " HD3",
				 "1HH1" : "HH11",
				 "2HH1" : "HH12",
				 "1HH2" : "HH21",
				 "2HH2" : "HH22" }

	dicNames = { "ASN" : asnNames,
				 "LEU" : leuNames,
				 "TYR" : tyrNames,
				 "ILE" : ileNames,
				 "GLN" : glnNames,
				 "TRP" : trpNames,
				 "LYS" : lysNames,
				 "ASP" : aspNames,
				 "GLY" : glyNames,
				 "PRO" : proNames,
				 "SER" : serNames,
				 "ARG" : argNames }

	def __init__( self, sequence, fileName ):
		self.sequence = sequence
		self.fileName = fileName
		self.dicContent = {}
		self.dicResults = {}

		self.readSequence()
		self.matchAminoAcids()
		self.writePDBFile()

	def updatePositions( self, original, newValues ):
		for j in range( 0, len( original ) ):
			original[j].xCor = newValues[j][0]
			original[j].yCor = newValues[j][1]
			original[j].zCor = newValues[j][2]

	def readSequence( self ):
		for i in range( 0, len( self.sequence ) ):
			fileName = None
			
			if self.sequence[i] is not None:
				fileName = self.dic.get( self.sequence[i] )

			if fileName is not None:
				if self.dicContent.get( self.sequence[i] ) is None:
					posAtoms = []
					translation = []
					atom = None
					pdb = open( fileName, "r" )

					stop = False
					content = []

					while not stop:
						line = pdb.readline()
						if not line:
							stop = True
						else:
							atom = Atom( line )
							if atom.getTag() == self.END_TAG:
								stop = True
							elif atom.getTag() == self.ATOM_TAG:
								content.append( atom )
								pos = map( float, atom.getPos() )
								if i == 0 and atom.getAtom() == self.ATOM_REF:
									translation = copy.deepcopy( pos )
									for j in xrange( 3 ):
										translation[j] *= -1

								posAtoms.append( pos )

					pdb.close()

					if i == 0:
						self.translateAtoms( posAtoms, translation )
						self.updatePositions( content, posAtoms )						

					self.dicContent[self.sequence[i]] = []
					self.dicContent[self.sequence[i]] = content

	def translateAtoms( self, posAtoms, translation ):
		for i in range( len( posAtoms ) ):
			for j in range( 3 ):
				posAtoms[i][j] += translation[j]

	def matchAminoAcids( self ):
		key = 0
		for i in range( 0, len( self.sequence ) - 1 ):
			j = i + 1
			keyI = self.sequence[i]
			keyJ = self.sequence[j]
			iPosOC = []
			jPosN = []
			iIndexOC = 0
			iIndexHC = 0
			jIndexHA = 0
			positionsJ = []

			if i == 0:
				keyContentI = copy.deepcopy( self.dicContent.get( keyI ) )

			keyContentJ = copy.deepcopy( self.dicContent.get( keyJ ) )

			if keyContentI is not None:
				for k in range( 0, len( keyContentI ) ):
					if keyContentI[k].getAtom() == self.OXIGEN_CARBOXYL:
						iIndexOC = k
						iPosOC = map( float, keyContentI[k].getPos() )

					elif keyContentI[k].getAtom() in self.HYDROGEN_CARBOXYL:
						iIndexHC = k

			if keyContentJ is not None:
				for k in range( 0, len( keyContentJ ) ):
					if keyContentJ[k].getAtom() in self.HYDROGEN_AMINO:
						jIndexHA = k

					if keyContentJ[k].getAtom() == self.NITROGEN:
						jPosN = map( float, keyContentJ[k].getPos() )

					positionsJ.append( map( float, keyContentJ[k].getPos() ) )

			if iIndexOC > iIndexHC:
				keyContentI.pop( iIndexOC )
				keyContentI.pop( iIndexHC )
			else:
				keyContentI.pop( iIndexHC )
				keyContentI.pop( iIndexOC )

			translation = [ 0, 0, 0 ]
			for k in xrange( 3 ):
				translation[k] = iPosOC[k] - jPosN[k]

			self.translateAtoms( positionsJ, translation )
			self.updatePositions( keyContentJ, positionsJ )
			self.dicResults[str( key )] = []
			self.dicResults[str( key )] = keyContentI

			keyContentJ.pop( jIndexHA )
			keyContentI = keyContentJ
			key += 1

		self.dicResults[str( key )] = []
		self.dicResults[str( key )] = keyContentI

		for z in range( 0, len( self.dicResults ) ):
			key = str( z )
			for y in range( 0, len( self.dicResults.get( key ) ) ):
				self.dicResults.get( key )[y].atom = self.renameAtom( self.dicResults.get( key )[y].atom, self.dicResults.get( key )[y].residue )

	def renameAtom( self, atom, amino ):
		dicAmino = self.dicNames.get( amino )
		if dicAmino is not None:
			if dicAmino.get( atom.strip() ) is not None:
				return dicAmino.get( atom.strip() )

		return atom

	def writePDBFile( self ):
		pdbNew = open( self.fileName, "w" )
		countTotal = 1
		print 'len', len( self.dicResults )
		for z in range( 0, len( self.dicResults ) ):
			key = str( z )
			for y in range( 0, len( self.dicResults.get( key ) ) ):
				self.dicResults.get( key )[y].serial = countTotal
				self.dicResults.get( key )[y].seqResidue = z+1				

				pdbNew.write( self.pdbPattern.format( self.dicResults.get( key )[y].tag, self.dicResults.get( key )[y].serial, self.dicResults.get( key )[y].atom, \
							  self.dicResults.get( key )[y].locIndicator, self.dicResults.get( key )[y].residue, self.dicResults.get( key )[y].chainID, self.dicResults.get( key )[y].seqResidue, \
							  self.dicResults.get( key )[y].insResidue, float( self.dicResults.get( key )[y].xCor ), float( self.dicResults.get( key )[y].yCor ), \
							  float( self.dicResults.get( key )[y].zCor ), float( self.dicResults.get( key )[y].occupancy ), float( self.dicResults.get( key )[y].temperature ), \
							  self.dicResults.get( key )[y].symbol, self.dicResults.get( key )[y].chargeAtom ) + "\n" )

				countTotal += 1

		pdbNew.write( "TER\n" )
		pdbNew.close()