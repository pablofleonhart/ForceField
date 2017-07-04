from aminoPhiPsi import AminoPhiPsi
from mapAminoAcids import AminoAcids
from pdbReader import PDBReader
from pdbAligner import PDBAligner
from acor import ACOR
import math
import numpy as np
import rmsd
import sys

class Builder( object ):
	experimental = None
	modified = None
	mod = None

	def __init__( self ):
		sequence = raw_input( "- Enter the desired aminoacid sequence or use the default (NLYIQWLKDGGPSSGRPPPS) by pressing 'Enter': " )

		if len( sequence ) == 0:
			sequence = "NLYIQWLKDGGPSSGRPPPS"

		if len( sequence ) > 1:
			fileName = "1L2Y-P.pdb"
			aminoAcids = AminoAcids( sequence, fileName )
			print "OK - The file '1L2Y-P.pdb' with the peptide bonds was generated."
			name = raw_input( "- Enter the PDB filename to calc the dihedral angles: Phi and Psi, or use the '" + "1L2Y-P.pdb" + "' by pressing 'Enter':" )

			if len( name ) == 0:
				name = fileName

			aminoPhiPsi = AminoPhiPsi( name )

			aminoPhiPsi.adjustOmegas()
			print "Omega", aminoPhiPsi.getOmegas()
			aminoPhiPsi.adjustPeptideBonds()
			aminoPhiPsi.adjustPhiPsi()
			print "Phi e Psi", aminoPhiPsi.getPhiPsi()
			aminoPhiPsi.writeAngles()
			aminoPhiPsi.plotRamanchandran()
			print "OK - The file 'aminoPhiPsi.txt' with the dihedral angles by amino acid was generated."
			print "OK - The file 'ramachandran.png' with the ramachandran map was generated."
			aminoPhiPsi.writePDBFile( "1L2Y-P.pdb" )

			self.readFiles( "files/1L2Y.pdb", "1L2Y-P.pdb" )
			self.calcRMSD()
			self.calcKabschRMSD()

			params = ['psi1', 'phi2', 'psi2', 'phi3', 'psi3', 'phi4', 'psi4', 'phi5', 'psi5', 'phi6', 'psi6', 'phi7', 'psi7', 'phi8',\
					  'psi8', 'phi9', 'psi9', 'phi10', 'psi10', 'phi11', 'psi11', 'phi12', 'psi12', 'phi13', 'psi13', 'phi14', 'psi14',\
					  'phi15', 'psi15', 'phi16', 'psi16', 'phi17', 'psi17', 'phi18', 'psi18', 'phi19', 'psi19', 'phi20']
			acor = ACOR( self.experimental, self.modified, params, False, 10 )
			acor.evolve()

		else:
			print "You must need specify at least two amino acids!"

	def readFiles( self, fileA, fileB ):
		self.experimental = PDBReader( fileA )
		self.modified = PDBReader( fileB )

		self.modified.adjustAtoms( self.experimental.atoms, self.experimental.aminoAcids )
		self.experimental.adjustAtoms( self.modified.atoms, self.modified.aminoAcids )

		'''print self.experimental.atoms
		print self.experimental.posAtoms
		print self.modified.atoms
		print self.modified.posAtoms'''

		self.experimental.calcBackbonePos()
		self.modified.calcBackbonePos()
		self.experimental.calcCaPos()
		self.modified.calcCaPos()

	def calcKabschRMSD( self ):
		P = np.array( self.experimental.posAtoms )
		Q = np.array( self.modified.posAtoms )
		#print rmsd.kabsch_rmsd( P, Q )
		P -= rmsd.centroid( P )
		Q -= rmsd.centroid( Q )
		print "{:15s} {:6.2f}".format( "Kabsch RMSD:", rmsd.kabsch_rmsd( P, Q ) )

	def calcRMSD( self ):
		print( len( self.experimental.atoms ), len( self.modified.atoms ) )
		#print( self.experimental.atoms )
		#print( self.modified.atoms )
		print( len( self.experimental.backbone ), len( self.modified.backbone ) )
		print( len( self.experimental.alpha ), len( self.modified.alpha ) )

		aligner = PDBAligner()
		print "{:15s} {:6.2f}".format( "CA RMSD:", aligner.calcRMSD( self.experimental.alpha, self.modified.alpha ) )
		print "{:15s} {:6.2f}".format( "Backbone RMSD:", aligner.calcRMSD( self.experimental.backbone, self.modified.backbone ) )
		print "{:15s} {:6.2f}".format( "All atoms RMSD:", aligner.calcRMSD( self.experimental.posAtoms, self.modified.posAtoms ) )

builder = Builder()