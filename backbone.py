import copy

class Backbone:
	NITROGEN_TAG = "N"
	ALPHA_TAG = "CA"
	CARBON_TAG = "C"
	posN = []
	posCA = []
	posC = []
	aminoAcid = None

	def setPosAtom( self, atom, aminoAcid, positions ):
		valid = False
		if atom == self.NITROGEN_TAG:
			self.posN = positions
			valid = True
		elif atom == self.ALPHA_TAG:
			self.posCA = positions
			valid = True
		elif atom == self.CARBON_TAG:
			self.posC = positions
			valid = True

		if valid:
			self.aminoAcid = aminoAcid

	def getPosN( self ):
		return copy.deepcopy( self.posN )

	def getPosCA( self ):
		return copy.deepcopy( self.posCA )

	def getPosC( self ):
		return copy.deepcopy( self.posC )

	def getAminoAcid( self ):
		return self.aminoAcid

	def _print( self ):
		print self.aminoAcid, self.posN, self.posCA, self.posC