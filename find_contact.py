import os
import subprocess
import sys
from tempfile import NamedTemporaryFile
from typing import List

import numpy as np
import pandas as pd
from pymol import cmd
from rnapolis import parser, annotator, tertiary
from contact_utils import (
    FileType,
    ProcessingException,
    ResiduesSelect,
    get_file,
    split_structure_to_hermetic_chains,
    structure_cif_filter,
)
from hbplus_process_cif import get_hbplus_result_for_large_structure
from hbplus_process_pdb import get_hbplus_result

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


def get_dot_bracket_and_filter(rna_file_path, in_contact):
    output_dot_bracket = NamedTemporaryFile(suffix=".dot", delete=False)
    with NamedTemporaryFile(suffix="_changed.cif") as cif_rnapolis:
        p1 = subprocess.Popen(
            "awk ' {if (NF > 3) {$7 = $17; $9 = $16; print; } else { print; }}' <"
            + rna_file_path,
            stdout=subprocess.PIPE,
            shell=True,
        )
        cif_rnapolis.write(p1.communicate()[0])
        structure3d = parser.read_3d_structure(open(cif_rnapolis.name, "r"))
        structure2d = annotator.extract_secondary_structure(structure3d)
        mapping = tertiary.Mapping2D3D(
            structure3d, structure2d.basePairs, structure2d.stackings, False
        )
        structure_cif_filter(rna_file_path, ResiduesSelect(in_contact))
        struct3d = dict()  # chain,seq_id,pairing
        for i, nucleotide in enumerate(
            [
                i.chain + "." + str(i.number) + "." + i.name
                for i in structure3d.residues
                if i.is_nucleotide
            ]
        ):
            if "".join(nucleotide.split(".")[:2]) in in_contact:
                if nucleotide.split(".")[0] not in struct3d.keys():
                    struct3d[nucleotide.split(".")[0]] = []
                struct3d[nucleotide.split(".")[0]].append(
                    [
                        nucleotide.split(".")[0],
                        nucleotide.split(".")[1],
                        nucleotide.split(".")[2],
                        mapping.bpseq.dot_bracket.structure[i],
                    ]
                )
        for chain in struct3d.keys():
            output_dot_bracket.write((f">strand_{chain}\n").encode())
            temp = np.array(struct3d[chain])
            output_dot_bracket.write(("".join(temp[:, 2]) + "\n").encode())
            output_dot_bracket.write(("".join(temp[:, 3]) + "\n").encode())
    output_dot_bracket.close()
    return output_dot_bracket.name


def get_dssp_and_filter(cif_prot_file_path, in_contact):
    with NamedTemporaryFile(suffix=".dssp", delete=False) as dssp:
        p1 = subprocess.Popen(
            f"""mkdssp {cif_prot_file_path} {dssp.name}""",
            stdout=subprocess.PIPE,
            shell=True,
        )
        p1.wait()
        p2 = subprocess.Popen(
            f"""awk -v residues="{','.join(in_contact)}" -f /opt/filter-secondary-structure.awk < {dssp.name}""",
            stdout=subprocess.PIPE,
            shell=True,
        )
        v = p2.communicate()[0].decode("utf-8")
        dssp.seek(0)
        dssp.write(v)
    structure_cif_filter(cif_prot_file_path, ResiduesSelect(in_contact))
    return dssp.name


def find_contact_hybrid_ligand_protein(pdb_id: str, model: int, chains: List[str]):
    hb_path_tsv = ""
    in_contact = set()
    structure_file_path, file_type = get_file(pdb_id, model, FileType.CIF)
    if file_type == FileType.CIF:
        hb_path_tsv, in_contact = get_hbplus_result_for_large_structure(
            structure_file_path, chains
        )
    elif file_type == FileType.PDB:
        hb_path_tsv, in_contact = get_hbplus_result(structure_file_path, chains)
    if os.path.exists(structure_file_path):
        os.remove(structure_file_path)
    structure_file_path_cif, _ = get_file(pdb_id, model, FileType.CIF)
    rna_file = NamedTemporaryFile(suffix="_RNA.cif", delete=False)
    dna_file = NamedTemporaryFile(suffix="_DNA.cif", delete=False)
    protein_file = NamedTemporaryFile(suffix="_PROTEIN.cif", delete=False)
    split_structure_to_hermetic_chains(
        structure_file_path_cif,
        [rna_file.name, dna_file.name, protein_file.name],
        model,
    )

    dot_file_path = get_dot_bracket_and_filter(rna_file.name, in_contact)
    cmd.delete("all")
    cmd.load(rna_file.name, f"{pdb_id.upper()}_RNA", quiet=1)
    cmd.save(rna_file.name.split(".")[0] + ".fasta", quiet=1)
    try:
        structure_cif_filter(dna_file.name, ResiduesSelect(in_contact))
        cmd.delete("all")
        cmd.load(dna_file.name, f"{pdb_id.upper()}_DNA", quiet=1)
        cmd.save(dna_file.name.split(".")[0] + ".fasta", quiet=1)
    except ValueError:
        os.remove(dna_file.name)
    try:
        dssp_file_path = get_dssp_and_filter(protein_file.name, in_contact)
        cmd.delete("all")
        cmd.load(protein_file.name, f"{pdb_id.upper()}_PROTEIN", quiet=1)
        cmd.save(protein_file.name.split(".")[0] + ".fasta", quiet=1)
    except ValueError:
        os.remove(protein_file.name)
    print(hb_path_tsv, in_contact)


if __name__ == "__main__":

    find_contact_hybrid_ligand_protein(
        sys.argv[1], int(sys.argv[2]), sys.argv[3].split(",")
    )
    # find_contact_ion()
