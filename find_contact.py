import copy
from tempfile import NamedTemporaryFile
from Bio.PDB import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser
from pymol import cmd

import subprocess
import requests
import tarfile
import os

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


def molecule_chain_assigment(structure):
    distribution = {}
    result_molecule_type = {}
    for model in structure:
        distribution[model.id] = {}
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


def molecule_distribution(structure):
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


class ProcessingException(Exception):
    pass


class ModelsSelect(Select):
    def __init__(self, m):
        self.models = int(m)

    def accept_residue(self, residue):
        # do I have to leave modified atoms
        return (
            residue.get_parent().get_parent().id + 1
            == self.models
            # and residue.id[0] == " "
        )


class ChainsSelect(Select):
    def __init__(self, c, modified_atoms=False):
        self.chains = set(c)
        self.modified_atoms = modified_atoms

    def accept_residue(self, residue):
        return (residue.id[0] == " " and self.modified_atoms) or not self.modified_atoms

    def accept_chain(
        self,
        chain,
    ):
        # do I have to leave modified atoms
        return chain in self.chains


class ResiduesSelect(Select):
    def __init__(self, r):
        self.residues = set(r)

    def accept_residue(self, residue):
        # do I have to leave modified atoms
        res_id = residue.get_parent().id + str(residue.id[1])
        if residue.id[2] != " ":
            res_id = res_id + residue.id[2]
        return res_id in self.residues


def structure_filter(pdb_file_path: str, selector: Select):
    structure = PDBParser(PERMISSIVE=0, QUIET=True).get_structure("str", pdb_file_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file_path, selector)


def get_pdb_file_maping(pdb_id) -> dict:  # dict of mapping
    response = requests.get(
        f"https://ftp.ebi.ac.uk/pub/databases/rcsb/pdb-remediated/compatible/pdb_bundle/v9/{pdb_id}/{pdb_id}-pdb-bundle.tar.gz",
        timeout=5000,
    )
    file_mapping = dict()
    if response.status_code != 200:
        raise ProcessingException("Cannot get file")
    with NamedTemporaryFile(suffix=".tar.gz", delete=False) as response_archive:
        response_archive.write(response.content)
        response_archive.close()
        tar = tarfile.open(response_archive.name, "r:gz")
        with NamedTemporaryFile(suffix=".txt", delete=False) as mapping_file:
            mapping_file.write(tar.extractfile(f"{pdb_id}-chain-id-mapping.txt"))
            mapping_file.close()
            with open(mapping_file, "r") as mapping_file_read:
                current_file = ""
                for line in mapping_file_read.readlines()[2:]:
                    row = line.split()
                    if len(row) == 1:
                        if "txt" in row[0] and current_file != row[0]:
                            current_file = str(row[0][:-1])
                    if len(row) == 2:
                        file_mapping[row[1]] = [current_file, row[0]]
            os.remove(mapping_file.name)
        os.remove(response_archive.name)
    return file_mapping


def get_pdb_file_from_bundle(pdb_id, file_from_bundle, save_to_path):  # dict of mapping
    response = requests.get(
        f"https://ftp.ebi.ac.uk/pub/databases/rcsb/pdb-remediated/compatible/pdb_bundle/v9/{pdb_id}/{pdb_id}-pdb-bundle.tar.gz",
        timeout=5000,
    )
    if response.status_code != 200:
        raise ProcessingException("Cannot get file")
    with NamedTemporaryFile(suffix=".tar.gz", delete=False) as response_archive:
        response_archive.write(response.content)
        response_archive.close()
        tar = tarfile.open(response_archive.name, "r:gz")
        with open(save_to_path, "wb") as mapping_file:
            mapping_file.write(tar.extractfile(file_from_bundle))
        os.remove(response_archive.name)


def get_pdb_file(pdb_id, chain_id, file_path):

    response = requests.get(
        f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=5000
    )
    if response.status_code != 200:
        try:
            file_mapping = get_pdb_file_maping(pdb_id)
            files = set()
            for idx, chain in enumerate(chain_id):
                chain_id[idx] = file_mapping[chain][1]
                files.add(file_mapping[chain][0])
            if len(files) != 1:
                raise ProcessingException("Too many files to process")

            get_pdb_file_from_bundle(pdb_id, list(files)[0], file_path)
        except ProcessingException:
            response = requests.get(
                f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=5000
            )
            with NamedTemporaryFile(suffix=".cif") as cif_file:
                with open(file_path, "wb") as cif_file:
                    cif_file.write(response.content)
                cmd.delete("all")
                cmd.load(cif_file.name)
                cmd.save(file_path)
    else:
        with open(file_path, "wb") as result_file:
            result_file.write(response.content)


def find_contact_ion():
    chains = ["B"]
    model = "1"
    pdb_id = "3d2v"
    with NamedTemporaryFile(suffix=".pdb") as pdb_file:
        get_pdb_file(pdb_id, chains, pdb_file.name)
        structure_filter(pdb_file.name, ModelsSelect(model))
        interaction_identificaton = subprocess.Popen(
            f"/home/solo/fingernat/code/fingeRNAt.py -r {pdb_file.name} -f SIMPLE",
            stdout=subprocess.PIPE,
            shell=True,
        )
        print(interaction_identificaton.communicate()[0].decode("utf-8"))


def split_structure_to_hermetic_chains(
    pdb_file_path, files_paths: tuple(str, str, str)
):
    structure = PDBParser(PERMISSIVE=0, QUIET=True).get_structure("str", pdb_file_path)
    distribution = molecule_distribution(structure)
    for idx, molecule in ["RNA", "DNA", "PROT"]:
        local_struct = copy.deepcopy(structure)
        io = PDBIO()
        io.set_structure(local_struct)
        io.save(files_paths[idx], ChainsSelect(distribution[molecule]))


def find_contact_hybrid_ligand_protein():
    chains = ["B"]
    model = "1"
    pdb_id = "6I0V".lower()
    with NamedTemporaryFile(suffix=".pdb") as pdb_file:
        get_pdb_file(pdb_id, chains, pdb_file.name)
        structure_filter(pdb_file.name, ModelsSelect(model))
        hydrogen_bound_extract = subprocess.Popen(
            f"hbplus {pdb_file.name}", shell=True, stdout=subprocess.DEVNULL
        )
        hydrogen_bound_extract.wait()

        if hydrogen_bound_extract.returncode != 0:
            raise ProcessingException("hbplus error")

        filter_hydrogen_bound = subprocess.Popen(
            "awk -vrna_chains="
            + (",".join(chains))
            + f" -f /opt/filter-hbonds.awk < {pdb_file.name.split('.')[0]}.hb2",
            stdout=subprocess.PIPE,
            shell=True,
        )
        filter_hydrogen_bound.communicate()[0].decode("utf-8")  # save?

        p1 = subprocess.Popen(
            f"awk -vrna_chains='{','.join(chains)}' -f /opt/filter-residues.awk < '{pdb_file.name.split('.')[0]}.hb2' ",
            stdout=subprocess.PIPE,
            shell=True,
        )
        p2 = subprocess.Popen(
            "sort -k1,1 -k2,2n ", stdin=p1.stdout, stdout=subprocess.PIPE, shell=True
        )
        p3 = subprocess.Popen(
            "awk '{a=$1 $2; if ($3!=\"-\") a = a $3; print a;}'" "",
            stdin=p2.stdout,
            stdout=subprocess.PIPE,
            shell=True,
        )
        p4 = subprocess.Popen(
            """awk '{printf("%s,",$0);}' """,
            stdin=p3.stdout,
            stdout=subprocess.PIPE,
            shell=True,
        )
        p5 = subprocess.Popen(
            """awk '{printf("%s",substr($0,1,length($0)-1));}'""",
            stdin=p4.stdout,
            stdout=subprocess.PIPE,
            shell=True,
        )
        # get pairing and then filter that which we are interested in

        residues_in_contact = p5.communicate()[0].decode("utf-8").split(",")
        # parse to fasta
        rna_file = NamedTemporaryFile(suffix="_RNA.pdb")
        dna_file = NamedTemporaryFile(suffix="_DNA.pdb")
        protein_file = NamedTemporaryFile(suffix="PROTEIN_.pdb")
        # if RNA protein
        split_structure_to_hermetic_chains(
            pdb_file.name, (rna_file.name, dna_file.name, protein_file.name)
        )
        structure_filter(pdb_file.name, ResiduesSelect(residues_in_contact))
        os.remove(f"{pdb_file.name.split('.')[0]}.hb2")


if __name__ == "__main__":
    molecule = "CR"
    match molecule:
        case "CR":
            find_contact_hybrid_ligand_protein()
        case "H":
            find_contact_hybrid_ligand_protein()
        case "A":
            # find_contact_hybrid_ligand_protein()
            find_contact_ion()
