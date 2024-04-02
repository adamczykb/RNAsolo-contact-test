# def into_molecule_files():
#     # parse to fasta
#     rna_file = NamedTemporaryFile(suffix="_RNA.cif")
#     dna_file = NamedTemporaryFile(suffix="_DNA.cif")
#     protein_file = NamedTemporaryFile(suffix="_PROTEIN.cif")
#     # if RNA protein
#     split_structure_to_hermetic_chains(
#         pdb_path, (rna_file.name, dna_file.name, protein_file.name)
#     )

#     structure_filter(pdb_path, ResiduesSelect(residues_in_contact))


import os
import subprocess
import sys
from tempfile import NamedTemporaryFile

from contact_utils import FileType, ProcessingException, get_file
from hbplus_process_cif import get_hbplus_result_for_large_structure
from hbplus_process_pdb import get_hbplus_result


def find_contact_ion():
    chains = ["B"]
    model = 1
    pdb_id = "3d2v"
    
    structure_file_path, file_type = get_file(pdb_id, model)
    match file_type:
        case FileType.CIF:
            os.remove(structure_file_path)
            raise ProcessingException("Cannot get PDB file for {pdb_id}")

    interaction_identificaton = subprocess.Popen(
        f"/home/solo/fingernat/code/fingeRNAt.py -r {structure_file_path} -f SIMPLE",
        stdout=subprocess.PIPE,
        shell=True,
    )
    print(interaction_identificaton.communicate()[0].decode("utf-8"))


def find_contact_hybrid_ligand_protein(pdb_id,chains):
    # chains = ["B"]
    model = 1
    # pdb_id = "6i0v".lower()
    hb_path_tsv = ""
    in_contact = set()
    structure_file_path, file_type = get_file(pdb_id, model)
    if file_type == FileType.CIF:
        hb_path_tsv, in_contact = get_hbplus_result_for_large_structure(
            structure_file_path, chains
        )
    elif file_type == FileType.PDB:
        hb_path_tsv, in_contact = get_hbplus_result(structure_file_path, chains)
    print(hb_path_tsv, in_contact)


if __name__ == "__main__":

    find_contact_hybrid_ligand_protein(sys.argv[1],sys.argv[2].split(','))
    # find_contact_ion()
