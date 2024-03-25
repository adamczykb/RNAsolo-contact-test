import Bio
import sys
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import os
import copy







Select = Bio.PDB.Select
parser = PDBParser(PERMISSIVE=1, QUIET=True)
structure = parser.get_structure("str", str(sys.argv[1]))





io = PDBIO()
io.set_structure(structure)
if len(rna_chains) > 0:
    chains = copy.deepcopy(rna_chains)
    io.save(append_id(str(sys.argv[1]), "RNA"), ChainsSelect(models, chains))
if len(dna_chains) > 0:
    chains = copy.deepcopy(dna_chains)
    io.save(append_id(str(sys.argv[1]), "DNA"), ChainsSelect(models, chains))
if len(protein_chains) > 0:
    chains = copy.deepcopy(protein_chains)
    io.save(append_id(str(sys.argv[1]), "Protein"), ChainsSelect(models, chains))
