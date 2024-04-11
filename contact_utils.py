import copy
from dataclasses import dataclass
import Bio
import requests


from tempfile import NamedTemporaryFile
from Bio.PDB import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser

from typing import Tuple
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBIO
from enum import Enum


class ProcessingException(Exception):
    pass


class MoleculeType(Enum):
    DNA = 0
    RNA = 1
    LIGAND = 2
    ION = 3
    PROTEIN = 4


@dataclass
class ContactResidue:
    resid: str
    value: str
    mtype: MoleculeType

    def __hash__(self) -> int:
        hash(self.resid)

    def __getitem__(self, indx):
        match indx:
            case 0:
                return self.resid
            case 1:
                return self.value
            case 2:
                return self.mtype


@dataclass
class InContactResult:
    hb_path_tsv: str = ""
    full_motif_cif: str = ""
    full_motif_fasta: str = ""
    in_contact = set()
    in_contact_desc = dict()
    dot_bracket_file_path: str = ""
    dssp_file_path: str = ""
    protein_file_cif: str = ""
    protein_file_fasta: str = ""
    dna_file_cif: str = ""
    dna_file_fasta: str = ""
    rna_file_cif: str = ""
    rna_file_fasta: str = ""
    ion_file_cif: str = ""
    ion_file_fasta: str = ""


CIF_HEADER = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "auth_seq_id",
    "auth_asym_id",
    "pdbx_PDB_model_num",
]
TSV_COLUMNS = [
    "donor",
    "donor_atom",
    "acceptor",
    "acceptor_atom",
    "distance",
    "atom_categories",
    "donor_acceptor_groups_gap",
    "CA_atoms_donor_acceptor_distance",
    "hydrogen_donor_acceptor_angle",
    "hydrogen_acceptor_distance",
    "acceptor_hydrogen_antecedent_angle",
    "donor_acceptor_antecedent_angle",
    "hydrogen_bonds_no",
]


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
            distribution[model.id][chain.id] = {
                "RNA": 0,
                "DNA": 0,
                "Protein": 0,
                "ION": 0,
            }
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
                elif (
                    residue.id[0] != " " and residue.id[0] != "W" and len(residue) == 1
                ):
                    distribution[model.id][chain.id]["ION"] = (
                        distribution[model.id][chain.id].pop("ION") + 1
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

    chain_assigment = {"RNA": {}, "DNA": {}, "PROT": {}, "ION": {}}
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
            elif (
                "ION" == result_molecule_type[model.id][chain.id]
                and distribution[model.id][chain.id]["ION"] > 0
            ):
                chain_assigment["ION"][model.id].append(chain.id)
    return chain_assigment


class ModelsSelect(Select):
    def __init__(self, m, save_ions=False, ligands=True):
        self.models = int(m)
        self.save_ions = save_ions
        self.ligands = ligands

    def accept_residue(self, residue):
        return (
            (
                residue.get_parent().get_parent().id + 1 == self.models
                and residue.id[0] == " "
            )
            or (
                residue.get_parent().get_parent().id + 1 == self.models
                and residue.id[0] != "W"
                and self.ligands
            )
            or (
                residue.get_parent().get_parent().id + 1 == self.models
                and self.save_ions
                and residue.id[0] != "W"
                and residue.id[0] != " "
            )
        )


class ChainsSelect(Select):
    def __init__(self, c, save_ions=False, ligands=True):
        self.chains = set(c)
        self.save_ions = save_ions
        self.ligands = ligands

    def accept_residue(self, residue):
        return (
            (residue.get_parent().id in self.chains and residue.id[0] != "W")
            or (
                residue.get_parent().id in self.chains
                and self.ligands
                and residue.id[0] != "W"
                and residue.id[0] != " "
                and len(residue) > 1
            )
            or (
                self.save_ions
                and residue.id[0] != "W"
                and residue.id[0] != " "
                and len(residue) == 1
            )
        )


class ResiduesSelect(Select):
    def __init__(self, r):
        self.residues = set(r)

    def accept_residue(self, residue):
        res_id = f"{residue.get_parent().id}.{residue.id[1]}"
        if residue.id[2] != " ":
            res_id = res_id + residue.id[2]
        return res_id in self.residues and residue.id[0] != "W"


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


def get_file(
    pdb_id: str,
    model_id: int = 0,
    file_type: FileType = None,
    save_ions=False,
    ligands=True,
) -> FileType:
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

            if model_id > 0:
                io = MMCIFIO()
                parser = MMCIFParser(QUIET=True)
                io.set_structure(parser.get_structure("str", cif_file.name))
                io.save(
                    cif_file.name,
                    ModelsSelect(model_id, save_ions=save_ions, ligands=ligands),
                )
            return cif_file.name, FileType.CIF

    def get_pdb():
        response = requests.get(
            f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=5000
        )
        with NamedTemporaryFile(suffix=".pdb", delete=False) as pdb_file:
            with open(pdb_file.name, "wb") as result_file:
                result_file.write(response.content)
            structure_filter(
                pdb_file.name,
                ModelsSelect(model_id, save_ions=save_ions, ligands=ligands),
            )
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
                structure_filter(
                    pdb_file.name,
                    ModelsSelect(model_id, save_ions=save_ions, ligands=ligands),
                )
                return pdb_file.name, FileType.PDB
    else:
        match file_type:
            case FileType.PDB:
                return get_pdb()
            case FileType.CIF:
                return get_cif()


def split_structure_to_hermetic_chains(
    cif_file_path, files_paths, model, in_contact_desc
):
    structure = MMCIFParser(QUIET=True).get_structure("str", cif_file_path)
    # distribution = structure_chain_molecules(structure)
    residue_assigment = {
        "RNA": set(),
        "DNA": set(),
        "PROT": set(),
        "ION": set(),
        "LIG": set(),
    }
    for key, value in in_contact_desc.items():
        for res in value:
            match res[2]:
                case MoleculeType.DNA:
                    residue_assigment["DNA"].add(res[0])
                case MoleculeType.LIGAND:
                    residue_assigment["LIG"].add(res[0])
                case MoleculeType.PROTEIN:
                    residue_assigment["PROT"].add(res[0])
                case MoleculeType.ION:
                    residue_assigment["ION"].add(res[0])
    residue_assigment["RNA"].update(list(in_contact_desc.keys()))
    for idx, molecule in enumerate(list(residue_assigment.keys())):
        local_struct = copy.deepcopy(structure)
        io = MMCIFIO()
        io.set_structure(local_struct)
        try:
            io.save(
                files_paths[idx],
                ResiduesSelect(residue_assigment[molecule]),
            )
        except:
            pass
