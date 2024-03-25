import Bio
import sys
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
import os

Select = Bio.PDB.Select
parser = PDBParser(PERMISSIVE=1, QUIET=True)
structure = parser.get_structure('str', str(sys.argv[1]))

if len(sys.argv) == 4: 
    models = str(sys.argv[3]).split(",")
else:
    models = None
    
residue_descs = []
models_to_remove = []
for model in structure:
    if str(model.serial_num) not in models:
        for idx, residue in reversed(list(enumerate(model.get_residues()))):
            chain_id = residue.get_parent().id
            model[chain_id].detach_child(residue.id)
        models_to_remove.append(model.id)

for model_id in reversed(models_to_remove):
    structure.detach_child(model_id)

io = PDBIO()
io.set_structure(structure)
io.save(str(sys.argv[2]))