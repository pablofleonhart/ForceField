from mapAminoAcids import AminoAcids
from pdbReader import PDBReader
import sys

sequence = raw_input( "- Enter the desired aminoacid sequence or use the default (VSCEDCPEHCSTQKAQAKCDNDKCVCEPI) by pressing 'Enter': " )

if len( sequence ) == 0:
	sequence = "VSCEDCPEHCSTQKAQAKCDNDKCVCEPI"

if len( sequence ) > 1:
	filename = "result.pdb"
	
	aminoAcids = AminoAcids( sequence, filename )
	pdb = PDBReader( filename )
	print pdb.atoms
	print pdb.posAtoms
	print pdb.aminoAcids
	print pdb.aAcids

	for k in pdb.dicContent.keys():
		pdb.dicContent.get( k )._print()
else:
	print "You must need specify at least two amino acids!"