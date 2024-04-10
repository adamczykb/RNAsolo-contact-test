import Bio
import sys
from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO

dna_dict = [ "DT", "DA", "DC", "DG" ]
rna_dict = [ "A", "C", "G", "U" ]
protein_dict = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" ]

Select = Bio.PDB.Select
def calc(path):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("str", path)

    distribution = {}
    molecule_chains = {}
    for model in structure:
        distribution[model.id] = {}
        molecule_chains[model.id] = {'RNA':[],'DNA':[],'Protein':[],'Ligand':[]}
        for chain in model:
            distribution[model.id][chain.id] = {'RNA':0,'DNA':0,'Protein':0,'Ligand':0}
            for residue in chain:
                res = residue.get_resname().strip()
                if res in protein_dict:
                    distribution[model.id][chain.id]['Protein'] = distribution[model.id][chain.id].pop('Protein') + 1
                elif res in dna_dict:
                    distribution[model.id][chain.id]['DNA'] = distribution[model.id][chain.id].pop('DNA') + 1
                elif res in rna_dict:
                    distribution[model.id][chain.id]['RNA'] = distribution[model.id][chain.id].pop('RNA') + 1
                elif res[0]!=" " and res[0]!="W" and len(res)>1:
                    distribution[model.id][chain.id]['Ligand'] = distribution[model.id][chain.id].pop('Ligand') + 1

            molecule_type = max(distribution[model.id][chain.id], key=distribution[model.id][chain.id].get)
            molecule_chains[model.id][molecule_type].append(chain.id)

    return structure,molecule_chains

if __name__ == "__main__":
    structure,molecule_chains=calc(sys.argv[1])
    for model in structure:
        print('Model no:{}'.format(model.serial_num))
        print('DNA:{}'.format(','.join(molecule_chains[model.id]['DNA'])))
        print('RNA:{}'.format(','.join(molecule_chains[model.id]['RNA'])))
        print('Protein:{}'.format(','.join(molecule_chains[model.id]['Protein'])))