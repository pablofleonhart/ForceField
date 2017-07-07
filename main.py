from aminoPhiPsi import AminoPhiPsi
from mapAminoAcids import AminoAcids
from pdbReader import PDBReader
from pdbAligner import PDBAligner
from acor import ACOR
import math
import numpy as np
import rmsd
import sys
from funcEnergy import EnergyFunction

class Builder( object ):
	def __init__( self ):
		self.experimental = None
		self.modified = None
		self.mod = None
		#sequence = raw_input( "- Enter the desired aminoacid sequence or use the default (NLYIQWLKDGGPSSGRPPPS) by pressing 'Enter': " )

		#if len( sequence ) == 0:
		sequence = "NLYIQWLKDGGPSSGRPPPS"

		if len( sequence ) > 1:
			fileName = "1L2Y-P.pdb"
			aminoAcids = AminoAcids( sequence, fileName )
			#print "OK - The file '1L2Y-P.pdb' with the peptide bonds was generated."
			#name = raw_input( "- Enter the PDB filename to calc the dihedral angles: Phi and Psi, or use the '" + "1L2Y-P.pdb" + "' by pressing 'Enter':" )

			#if len( name ) == 0:
			name = fileName

			aminoPhiPsi = AminoPhiPsi( name )

			aminoPhiPsi.adjustOmegas()
			print "Omega", aminoPhiPsi.getOmegas()
			aminoPhiPsi.adjustPeptideBonds()
			aminoPhiPsi.adjustPhiPsi()
			print "Phi e Psi", aminoPhiPsi.getPhiPsi()
			aminoPhiPsi.writeAngles()
			#aminoPhiPsi.plotRamanchandran()
			#print "OK - The file 'aminoPhiPsi.txt' with the dihedral angles by amino acid was generated."
			#print "OK - The file 'ramachandran.png' with the ramachandran map was generated."
			aminoPhiPsi.writePDBFile( "1L2Y-P.pdb" )

			'''ang = [[math.pi*2, 1.2741482132786572], 
				[2.2283201680746201, 2.6267633751735717], 
				[1.0634146573400245, 0.52223129258622736], 
				[1.6050099640130078, 0.30649803705741041], 
				[-2.2040803357176495, -0.90176292448857343], 
				[0.079934578824690128, 0.27756899349685993], 
				[0.18108129015152885, 2.1142009127176769], 
				[-0.41623437621905968, -0.20717239853881342], 
				[2.9469258394311248, 1.733020396724231], 
				[0.47278485681107796, 1.3929740117443048], 
				[-3.045035654137688, 1.2461557148347264], 
				[1.142121268084795, 1.7031436014818488],
				 [-1.6539304901151612, -2.9477680057792339], 
				 [1.4419871869212022, -2.9491294143075724], 
				 [0.76509826318785779, 2.3151450233906496], 
				 [0.87020864291200439, -1.9501523674608274], 
				 [-1.4609214206122387, 2.4248772347632457], 
				 [1.8412165754693706, 0.53752326661797145], 
				 [0.62560371156870076, 2.9850669485863666], 
				 [2.8679314685005073, math.pi*2]]

			#ang = [[math.pi/2.9 for x in range( 2 )] for y in range( 20 )]

			aminoPhiPsi.adjustPhiPsi( ang )
			aminoPhiPsi.writePDBFile( "1L2Y-F.pdb" )'''

			self.readFiles( "files/1L2Y.pdb", "1L2Y-F.pdb" )
			self.calcRMSD()
			self.calcKabschRMSD()

			#print( len( self.experimental.atoms ), len( self.modified.atoms ) )
			params = ['psi1', 'phi2', 'psi2', 'phi3', 'psi3', 'phi4', 'psi4', 'phi5', 'psi5', 'phi6', 'psi6', 'phi7', 'psi7', 'phi8',\
					  'psi8', 'phi9', 'psi9', 'phi10', 'psi10', 'phi11', 'psi11', 'phi12', 'psi12', 'phi13', 'psi13', 'phi14', 'psi14',\
					  'phi15', 'psi15', 'phi16', 'psi16', 'phi17', 'psi17', 'phi18', 'psi18', 'phi19', 'psi19', 'phi20']
			acor = ACOR( self.experimental, self.modified, params, False, 200 )
			acor.evolve()

			'''fe = EnergyFunction( "1PLX-P.pdb" )
			print fe.getEnergy()'''
			'''app = AminoPhiPsi( "files/1L2Y.pdb" )
			fa = EnergyFunction( app.pdb )
			print fa.getEnergy()

			app = AminoPhiPsi( "1L2Y-P.pdb" )
			fa = EnergyFunction( app.pdb )
			print fa.getEnergy()

			app = AminoPhiPsi( "1L2Y-F.pdb" )
			fa = EnergyFunction( app.pdb )
			print fa.getEnergy()'''

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
		P = np.array( self.experimental.alpha )
		Q = np.array( self.modified.alpha )
		#print rmsd.kabsch_rmsd( P, Q )
		P -= rmsd.centroid( P )
		Q -= rmsd.centroid( Q )
		print "{:15s} {:6.2f}".format( "Kabsch RMSD:", rmsd.kabsch_rmsd( P, Q ) )

		P = np.array( self.experimental.backbone )
		Q = np.array( self.modified.backbone )
		#print rmsd.kabsch_rmsd( P, Q )
		P -= rmsd.centroid( P )
		Q -= rmsd.centroid( Q )
		print "{:15s} {:6.2f}".format( "Kabsch RMSD:", rmsd.kabsch_rmsd( P, Q ) )

		P = np.array( self.experimental.posAtoms )
		Q = np.array( self.modified.posAtoms )
		#print rmsd.kabsch_rmsd( P, Q )
		P -= rmsd.centroid( P )
		Q -= rmsd.centroid( Q )
		print "{:15s} {:6.2f}".format( "Kabsch RMSD:", rmsd.kabsch_rmsd( P, Q ) )

	def calcRMSD( self ):
		'''print( len( self.experimental.atoms ), len( self.modified.atoms ) )
		print( self.experimental.atoms )
		print( self.modified.atoms )
		print( len( self.experimental.backbone ), len( self.modified.backbone ) )
		print( len( self.experimental.alpha ), len( self.modified.alpha ) )'''

		aligner = PDBAligner()
		print "{:15s} {:6.2f}".format( "CA RMSD:", aligner.calcRMSD( self.experimental.alpha, self.modified.alpha ) )
		print "{:15s} {:6.2f}".format( "Backbone RMSD:", aligner.calcRMSD( self.experimental.backbone, self.modified.backbone ) )
		print "{:15s} {:6.2f}".format( "All atoms RMSD:", aligner.calcRMSD( self.experimental.posAtoms, self.modified.posAtoms ) )

builder = Builder()
