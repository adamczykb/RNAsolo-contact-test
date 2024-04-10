import copy
import Bio
import numpy as np
import requests


from tempfile import NamedTemporaryFile
from Bio.PDB import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser

from typing import Tuple
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBIO
from enum import Enum
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

class ProcessingException(Exception):
    pass


class FileType(Enum):
    PDB = 0
    CIF = 1


DNA_DICT = ["DT", "DA", "DC", "DG"]
RNA_DICT = ["A", "C", "G", "U"]
PROTEIN_DICT = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]


def molecule_chain_assigment(structure: Bio.PDB.Structure) -> Tuple[dict, dict]:
    """
    Recoginse and store cardinality of chain molecule affiliation

    Parameters:
            structure (Bio.PDB.Structure): Structure to process
    Returns:
            distribution (dict): Dict of molecule cardinality for each model and chain
            result_molecule_type (dict): Dict of molecule affilation for each model and chain
    """

    distribution = {}
    result_molecule_type = {}
    for model in structure:
        distribution[model.id] = {}
        result_molecule_type[model.id] = {}
        for chain in model:
            distribution[model.id][chain.id] = {"RNA": 0, "DNA": 0, "Protein": 0}
            for residue in chain:
                res = residue.get_resname().strip()
                if res in PROTEIN_DICT:
                    distribution[model.id][chain.id]["Protein"] = (
                        distribution[model.id][chain.id].pop("Protein") + 1
                    )
                elif res in DNA_DICT:
                    distribution[model.id][chain.id]["DNA"] = (
                        distribution[model.id][chain.id].pop("DNA") + 1
                    )
                elif res in RNA_DICT:
                    distribution[model.id][chain.id]["RNA"] = (
                        distribution[model.id][chain.id].pop("RNA") + 1
                    )
            result_molecule_type[model.id][chain.id] = max(
                distribution[model.id][chain.id],
                key=distribution[model.id][chain.id].get,
            )
    return result_molecule_type, distribution


def structure_chain_molecules(structure):
    """
    Structure chain molecule affiliaton

    Parameters:
            structure (Bio.PDB.Structure): Structure to process
    Returns:
            result_molecule_type (dict): Dict of molecule affilation for each model and chain
    """
    result_molecule_type, distribution = molecule_chain_assigment(structure)

    chain_assigment = {"RNA": {}, "DNA": {}, "PROT": {}}
    for model in structure:
        for chain in model:
            if (
                "RNA" == result_molecule_type[model.id][chain.id]
                or "DNA" == result_molecule_type[model.id][chain.id]
                and distribution[model.id][chain.id]["RNA"] > 0
            ):
                if model.id not in chain_assigment["RNA"]:
                    chain_assigment["RNA"][model.id] = []
                chain_assigment["RNA"][model.id].append(chain.id)
            elif (
                "RNA" == result_molecule_type[model.id][chain.id]
                or "DNA" == result_molecule_type[model.id][chain.id]
                and distribution[model.id][chain.id]["DNA"] > 0
            ):
                if model.id not in chain_assigment["DNA"]:
                    chain_assigment["DNA"][model.id] = []
                chain_assigment["DNA"][model.id].append(chain.id)
            elif (
                "Protein" == result_molecule_type[model.id][chain.id]
                and distribution[model.id][chain.id]["Protein"] > 0
            ):
                if model.id not in chain_assigment["PROT"]:
                    chain_assigment["PROT"][model.id] = []
                chain_assigment["PROT"][model.id].append(chain.id)
    return chain_assigment


class ModelsSelect(Select):
    def __init__(self, m, save_ions=False,ligands =True):
        self.models = int(m)
        self.save_ions = save_ions
        self.ligands = ligands

    def accept_residue(self, residue):     
        return (residue.get_parent().get_parent().id + 1 == self.models and residue.id[0] == " ") or  \
            (residue.get_parent().get_parent().id + 1 == self.models and residue.id[0] != "W" and self.ligands) or \
                (residue.get_parent().get_parent().id + 1 == self.models and self.save_ions and residue.id[0] != "W" and residue.id[0]!=" ")


class ChainsSelect(Select):
    def __init__(self, c, save_ions=False,ligands =True):
        self.chains = set(c)
        self.save_ions = save_ions
        self.ligands = ligands

    # def accept_residue(self, residue):
    #     return (residue.id[0] == " " and self.modified_atoms) or not self.modified_atoms

    # def accept_chain(
    #     self,
    #     chain,
    # ):
    #     # do I have to leave modified atoms
    #     return chain.id in self.chains 
    def accept_residue(self, residue):
        return  (residue.get_parent().id in self.chains and residue.id[0] != "W") or  \
            (residue.get_parent().id in self.chains and residue.id[0] != "W" and self.ligands and len(residue)>1) or \
                (self.save_ions and residue.id[0] != "W" and residue.id[0]!=" " and len(residue)==1)

class ResiduesSelect(Select):
    def __init__(self, r,ligands =True):
        self.residues = set(r)
        self.ligands=ligands
    def accept_residue(self, residue):
        res_id = residue.get_parent().id + str(residue.id[1])
        if residue.id[2] != " ":
            res_id = res_id + residue.id[2]
        return (res_id in self.residues and self.ligands) or (residue.id[0] == " " and not self.ligands and res_id in self.residues)


def structure_filter(pdb_file_path: str, selector: Select):
    structure = PDBParser(PERMISSIVE=0, QUIET=True).get_structure("str", pdb_file_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file_path, selector)

def structure_cif_filter(cif_file_path: str, selector: Select):
    structure = MMCIFParser(QUIET=True).get_structure("str", cif_file_path)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(cif_file_path, selector)

def get_file(pdb_id: str, model_id: int=0, file_type: FileType = None,save_ions=False,ligands =True) -> FileType:
    """
    Get structure file from rcsb.org

    Parameters:
            pdb_id (str): 4 literals structure ID
            model_id (int): Structure model to filter
            file_path (str): File path to store downloaded structure
    Returns:
            result_molecule_type (dict): Dict of molecule affilation for each model and chain
    """
    def get_cif():
        response = requests.get(
                f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=5000
            )
        if response.status_code != 200:
            raise ProcessingException("Cannot get structure")
        with NamedTemporaryFile(suffix=".cif", delete=False) as cif_file:
                with open(cif_file.name, "wb") as result_file:
                    result_file.write(response.content)

                # io = MMCIFIO()
                # parser = MMCIFParser(QUIET=True)
                # struct_dict = MMCIF2Dict(cif_file.name)
                # struct_dict["_atom_site.label_seq_id"]=struct_dict["_atom_site.auth_seq_id"]
                # io.set_dict(struct_dict)
                # io.save(cif_file.name)

                if model_id>0:
                    io = MMCIFIO()
                    parser = MMCIFParser(QUIET=True)
                    io.set_structure(parser.get_structure("str", cif_file.name))
                    io.save(cif_file.name, ModelsSelect(model_id,save_ions=save_ions,ligands=ligands))
                return cif_file.name, FileType.CIF
    def get_pdb():  
        response = requests.get(
            f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=5000
        )
        with NamedTemporaryFile(suffix=".pdb", delete=False) as pdb_file:
                with open(pdb_file.name, "wb") as result_file:
                    result_file.write(response.content)
                structure_filter(pdb_file.name, ModelsSelect(model_id,save_ions=save_ions,ligands=ligands))
                return pdb_file.name, FileType.PDB
    if file_type is None:
        response = requests.get(
            f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=5000
        )
        if response.status_code != 200:
            return get_cif()
        else:
            with NamedTemporaryFile(suffix=".pdb", delete=False) as pdb_file:
                with open(pdb_file.name, "wb") as result_file:
                    result_file.write(response.content)
                structure_filter(pdb_file.name, ModelsSelect(model_id,save_ions=save_ions,ligands=ligands))
                return pdb_file.name, FileType.PDB
    else:
        match file_type:
            case FileType.PDB:
                return get_pdb()
            case FileType.CIF:
                return get_cif()


def split_structure_to_hermetic_chains(cif_file_path, files_paths,model):
    structure = MMCIFParser(QUIET=True).get_structure("str", cif_file_path)
    distribution = structure_chain_molecules(structure)

    for idx, molecule in enumerate(["RNA", "DNA", "PROT"]):
        local_struct = copy.deepcopy(structure)
        io = MMCIFIO()
        io.set_structure(local_struct)
        try:
            io.save(files_paths[idx], ChainsSelect(distribution[molecule][model-1],ligands=True))
        except:
            pass
