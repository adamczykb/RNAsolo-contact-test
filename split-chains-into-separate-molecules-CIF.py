import Bio
import sys
from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
import os
import copy

dna_dict = [ "DT", "DA", "DC", "DG" ]
rna_dict = [ "A", "C", "G", "U" ]
protein_dict = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]

def append_id(filename, molecule_type):
    name, ext = os.path.splitext(filename)
    return "{name}_{molecule_type}{ext}".format(name=name, molecule_type=molecule_type, ext=ext)

Select = Bio.PDB.Select
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("str", str(sys.argv[1]))

distribution = {}
rna_chains = {}
dna_chains = {}
protein_chains = {}
chains = {}
for model in structure:
    distribution[model.id] = {}
    for chain in model:
        distribution[model.id][chain.id] = {'RNA':0,'DNA':0,'Protein':0}
        for residue in chain:
            res = residue.get_resname().strip()
            if res in protein_dict:
                distribution[model.id][chain.id]['Protein'] = distribution[model.id][chain.id].pop('Protein') + 1
            elif res in dna_dict:
                distribution[model.id][chain.id]['DNA'] = distribution[model.id][chain.id].pop('DNA') + 1
            elif res in rna_dict:
                distribution[model.id][chain.id]['RNA'] = distribution[model.id][chain.id].pop('RNA') + 1
        molecule_type = max(distribution[model.id][chain.id], key=distribution[model.id][chain.id].get)
        if 'RNA' == molecule_type or 'DNA' == molecule_type and distribution[model.id][chain.id]['RNA'] > 0:
            if model.id not in rna_chains:
                rna_chains[model.id] = []
            rna_chains[model.id].append(chain.id)
        elif 'RNA' == molecule_type or 'DNA' == molecule_type and distribution[model.id][chain.id]['DNA'] > 0:
            if model.id not in dna_chains:
                dna_chains[model.id] = []
            dna_chains[model.id].append(chain.id)
        elif 'Protein' == molecule_type and distribution[model.id][chain.id]['Protein'] > 0:
            if model.id not in protein_chains:
                protein_chains[model.id] = []
            protein_chains[model.id].append(chain.id)

if len(sys.argv) == 3:
    models = str(sys.argv[2]).split(",")
else:
    models = None

class ChainsSelect(Select):
    models
    chains
    def __init__(self, m, c):
      self.models = m
      self.chains = c
    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        model_id = residue.get_parent().get_parent().id
        if (len(chains)>0) and ((models == None) or (models != None and str(model_id + 1) in models)) and (residue.get_parent().id in self.chains[model_id] and residue.id[0] == ' '):
            return True
        return False

io=MMCIFIO()
io.set_structure(structure)
if len(rna_chains)>0:
    chains = copy.deepcopy(rna_chains)
    io.save(append_id(str(sys.argv[1]),'RNA'), ChainsSelect(models, chains))
if len(dna_chains)>0:
    chains = copy.deepcopy(dna_chains)
    io.save(append_id(str(sys.argv[1]),'DNA'), ChainsSelect(models, chains))
if len(protein_chains)>0:
    chains = copy.deepcopy(protein_chains)
    io.save(append_id(str(sys.argv[1]),'Protein'), ChainsSelect(models, chains))