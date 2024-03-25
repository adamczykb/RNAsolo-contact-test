import sys
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO

# Select = Bio.PDB.Select
parser = PDBParser(PERMISSIVE=1, QUIET=True)
structure = parser.get_structure('str', str(sys.argv[1]))

selected_residues = str(sys.argv[3]).split(',')
#print(selected_residues)
    
residue_descs = []
for model in structure:
    for idx, residue in reversed(list(enumerate(model.get_residues()))):
        chain_id = residue.get_parent().id
        a,serial,icode = residue.id
        res_id = chain_id + str(serial)
        if (icode != ' '):
            res_id = res_id + icode
        if res_id not in selected_residues: 
            model[chain_id].detach_child(residue.id)
        else:
            residue_descs.insert(0,'{}:{}'.format(res_id,residue.get_resname()))
          
io = PDBIO()
io.set_structure(structure)
io.save(str(sys.argv[2]))
