class Atom( object ):
	tag = ''
	serial = ''
	atom = ''
	locIndicator = ''
	residue = ''
	chainID = ''
	seqResidue = ''
	insResidue = ''
	xCor = ''
	yCor = ''
	zCor = ''
	occupancy = ''
	temperature = ''
	symbol = ''
	chargeAtom = ''

	def __init__( self, line ):
		self.tag = line[0:6]
		self.serial = line[6:11]
		self.atom = line[12:16]
		self.locIndicator = line[16:17]
		self.residue = line[17:20]
		self.chainID = line[21:22]
		self.seqResidue = line[22:26]
		self.insResidue = line[26:27]
		self.xCor = line[30:38]
		self.yCor = line[38:46]
		self.zCor = line[46:54]
		self.occupancy = line[54:60]
		self.temperature = line[60:66]
		self.symbol = line[76:78]
		self.chargeAtom = line[78:80]

	def getTag( self ):
		return self.tag.strip()

	def getAtom( self ):
		return self.atom.strip()

	def getResidue( self ):
		return self.residue.strip()

	def getPos( self ):
		return map( float, ( self.xCor, self.yCor, self.zCor ) )