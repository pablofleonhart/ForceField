from aminoPhiPsi import AminoPhiPsi
from mapAminoAcids import AminoAcids
from pdbReader import PDBReader
import sys

sequence = raw_input( "- Enter the desired aminoacid sequence or use the default (NLYIQWLKDGGPSSGRPPPS) by pressing 'Enter': " )

if len( sequence ) == 0:
	sequence = "NLYIQWLKDGGPSSGRPPPS"

if len( sequence ) > 1:
	fileName = "1L2Y-P2.pdb"
	aminoAcids = AminoAcids( sequence, fileName )
	print "OK - The file '1L2Y-P.pdb' with the peptide bonds was generated."
	name = raw_input( "- Enter the PDB filename to calc the dihedral angles: Phi and Psi, or use the '" + "1L2Y-P.pdb" + "' by pressing 'Enter':" )

	if len( name ) == 0:
		name = fileName

	aminoPhiPsi = AminoPhiPsi( name )

	print "Phi e Psi", aminoPhiPsi.getPhiPsi()
	print "Omega", aminoPhiPsi.getOmegas()
	aminoPhiPsi.adjustOmegas()
	print "Omega", aminoPhiPsi.getOmegas()
	aminoPhiPsi.set_peptide_bond_angles()
	#print aminoPhiPsi.get_peptide_bond_angles()		
	aminoPhiPsi.rotateStructure()
	print "Phi e Psi", aminoPhiPsi.getPhiPsi()
	aminoPhiPsi.calcAngles()
	aminoPhiPsi.plotRamanchandran()
	print "OK - The file 'aminoPhiPsi.txt' with the dihedral angles by amino acid was generated."
	print "OK - The file 'ramachandran.png' with the ramachandran map was generated."
	aminoPhiPsi.writePDBFile( "1L2Y-P2.pdb" )
else:
	print "You must need specify at least two amino acids!"