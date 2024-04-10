from dataclasses import dataclass
import os
import subprocess
import sys
from tempfile import NamedTemporaryFile
from typing import List, Set
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.mmcifio import MMCIFIO

import numpy as np
import pandas as pd
import requests
from pymol import cmd
from rnapolis import parser, annotator, tertiary
from contact_utils import (
    ChainsSelect,
    FileType,
    ModelsSelect,
    ProcessingException,
    ResiduesSelect,
    get_file,
    split_structure_to_hermetic_chains,
    structure_cif_filter,
)
from hbplus_process_cif import MoleculeType, get_hbplus_result_for_large_structure
from hbplus_process_pdb import get_hbplus_result


@dataclass
class DNAProtLigandResult:
    hb_path_tsv: str = ""
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


def find_contact_ion(pdb_id: str, model: int, chains: List[str]):
    output_path = ""
    motif_residues_clean = set()
    in_contact_desc=dict()
    in_contact=set()
    structure_file_path, file_type = get_file(pdb_id, model, FileType.CIF,save_ions=True,ligands=False)
    for chain in chains:
        with NamedTemporaryFile(suffix=".cif") as structure_file_copy_cif:
            structure_file_copy_cif.write(open(structure_file_path,'rb').read())
            output_path=''
            with NamedTemporaryFile(suffix=".pdb") as structure_file_path_pdb:
                
                structure_cif_filter(
                    structure_file_copy_cif.name, ChainsSelect(chain, save_ions=True,ligands=False)
                )

                response = requests.post(
                    "http://tomek:8080",
                    headers={"Content-Type": "text/plain"},
                    data=open(structure_file_copy_cif.name, "rb"),
                    timeout=10000,
                )
                if response.status_code == 200:
                    with open(
                        structure_file_path_pdb.name, "wb"
                    ) as pdb_file_output:
                        pdb_file_output.write(response.content)
                
                interaction_identificaton = subprocess.Popen(
                    f"/opt/fingernat/code/fingeRNAt.py -r {structure_file_path_pdb.name} -f SIMPLE",
                    stdout=subprocess.PIPE,
                    stderr=subprocess.DEVNULL,
                    shell=True,
                )
                output_path = (
                    interaction_identificaton.communicate()[0].decode("utf-8").strip()
                )
            if output_path=='':
                continue
            filter_awk = subprocess.Popen(
                f"""
                    awk -f /opt/infer-residues.awk < {output_path}  
                """,
                stdout=subprocess.PIPE,
                shell=True,
            )
            filter_awk2 = subprocess.Popen(
                "sort -k1,1 -k2,2n",
                stdout=subprocess.PIPE,
                stdin=filter_awk.stdout,
                shell=True,
            )
            filter_awk3 = subprocess.Popen(
                'awk \'{a=$1"."$2; if ($3!="-") a = a $3; print a;}\'',
                stdout=subprocess.PIPE,
                stdin=filter_awk2.stdout,
                shell=True,
            )
            filter_awk4 = subprocess.Popen(
                "awk '{printf(\"%s,\",$0);}'",
                stdout=subprocess.PIPE,
                stdin=filter_awk3.stdout,
                shell=True,
            )
            filter_awk5 = subprocess.Popen(
                "awk '{printf(\"%s\",substr($0,1,length($0)-1));}'",
                stdout=subprocess.PIPE,
                stdin=filter_awk4.stdout,
                shell=True,
            )
            motif_residues = filter_awk5.communicate()[0].decode("utf-8")

            for residue in motif_residues.split(","):
                if residue.split(".")[0] in chains:
                    motif_residues_clean.add(residue)
            tb = pd.read_csv(output_path,sep='\t')
            tb= tb.set_axis([tb.columns[0]]+['.'.join(name.split('#')[1].split(':')[::-1]) for name in tb.columns[1:]],axis=1)
            tb['Ligand_name']=(tb['Ligand_name'].str.split(':').str[1:3].str[::-1]).apply('.'.join) + ":" +tb['Ligand_name'].str.split(':').str[0]
            tb=tb.rename({'Ligand_name':'ion'},axis=1)
            tb=tb.T
            tb.columns=tb.iloc[0]
            tb=tb[1:]
            for ion in  list(tb.columns):
                ion_name=ion.split(':')
                for residue in list(tb[tb[ion]==1].index):
                    if residue not in in_contact_desc.keys():
                        in_contact_desc[residue]=[]
                    in_contact_desc[residue].append([ion_name[0],ion_name[1],MoleculeType.ION])
            in_contact.update(list(in_contact_desc.keys()))
    structure_cif_filter(
        structure_file_path, ResiduesSelect(list(motif_residues_clean))
    )
    
    return structure_file_path, in_contact,in_contact_desc


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
        structure_cif_filter(rna_file_path, ResiduesSelect(in_contact,ligands=False))
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
    with NamedTemporaryFile(suffix="_dssp.cif") as dssp:
        p1 = subprocess.Popen(
            f"""mkdssp --output-format mmcif {cif_prot_file_path} {dssp.name}""",
            stdout=subprocess.PIPE,
            shell=True,
        )
        p1.wait()

        struct_dict = MMCIF2Dict(dssp.name)
        keys_from_cif = list(
            filter(
                lambda k: k.startswith("_struct_conf.")
                or k.startswith("_struct_conf_type."),
                struct_dict.keys(),
            )
        )
        keys_from_cif_to_save = list(
            filter(
                lambda k: k.startswith("_struct_conf.")
                or k.startswith("_struct_conf_type.")
                or k == "data_",
                struct_dict.keys(),
            )
        )
        for i in list(struct_dict.keys()):
            if i not in keys_from_cif_to_save:
                del struct_dict[i]

        for ind in range(len(struct_dict["_struct_conf.id"]) - 1, -1, -1):
            if (
                struct_dict["_struct_conf.end_auth_asym_id"][ind]
                + struct_dict["_struct_conf.end_auth_seq_id"][ind]
                in in_contact
                or struct_dict["_struct_conf.beg_auth_asym_id"][ind]
                + struct_dict["_struct_conf.beg_auth_seq_id"][ind]
                in in_contact
            ):
                pass
            else:
                for key in list(
                    filter(lambda k: k.startswith("_struct_conf."), keys_from_cif)
                ):
                    del struct_dict[key][ind]
        for ind in range(len(struct_dict["_struct_conf_type.id"]) - 1, -1, -1):
            if (
                struct_dict["_struct_conf_type.id"][ind]
                not in struct_dict["_struct_conf.conf_type_id"]
            ):
                del struct_dict["_struct_conf_type.id"][ind]
                del struct_dict["_struct_conf_type.criteria"][ind]
        io = MMCIFIO()
        io.set_dict(struct_dict)
        io.save(dssp.name)
        structure_cif_filter(cif_prot_file_path, ResiduesSelect(in_contact))
        with open(cif_prot_file_path, "ab") as cif_prot:
            cif_prot.writelines(dssp.readlines()[2:])


def find_contact_hybrid_ligand_protein(
    pdb_id: str, model: int, chains: List[str]
) -> DNAProtLigandResult:
    result = DNAProtLigandResult()
    hb_path_tsv = ""
    in_contact = set()
    structure_file_path, file_type = get_file(pdb_id, model, FileType.CIF,save_ions=False,ligands=True)
    hb_path_tsv, in_contact, in_contact_desc = get_hbplus_result_for_large_structure(
        structure_file_path, chains
    )
    result.hb_path_tsv = hb_path_tsv
    result.in_contact = in_contact
    result.in_contact_desc = in_contact_desc

    structure_file_path_cif, _ = get_file(pdb_id, model, FileType.CIF)
    rna_file = NamedTemporaryFile(suffix="_RNA.cif", delete=False)
    dna_file = NamedTemporaryFile(suffix="_DNA.cif", delete=False)
    protein_file = NamedTemporaryFile(suffix="_PROTEIN.cif", delete=False)
    split_structure_to_hermetic_chains(
        structure_file_path_cif,
        [rna_file.name, dna_file.name, protein_file.name],
        model,
    )

    result.dot_bracket_file_path = get_dot_bracket_and_filter(rna_file.name, in_contact)
    cmd.delete("all")
    cmd.load(rna_file.name, f"{pdb_id.upper()}_RNA", quiet=1)
    cmd.save(rna_file.name.split(".")[0] + ".fasta", quiet=1)

    result.rna_file_cif = rna_file.name
    result.rna_file_fasta = rna_file.name.split(".")[0] + ".fasta"

    try:
        structure_cif_filter(dna_file.name, ResiduesSelect(in_contact))
        cmd.delete("all")
        cmd.load(dna_file.name, f"{pdb_id.upper()}_DNA", quiet=1)
        cmd.save(dna_file.name.split(".")[0] + ".fasta", quiet=1)
        result.dna_file_cif = dna_file.name
        result.dna_file_fasta = dna_file.name.split(".")[0] + ".fasta"
    except ValueError:
        os.remove(dna_file.name)
    try:
        get_dssp_and_filter(protein_file.name, in_contact)
        cmd.delete("all")
        cmd.load(protein_file.name, f"{pdb_id.upper()}_PROTEIN", quiet=1)
        cmd.save(protein_file.name.split(".")[0] + ".fasta", quiet=1)
        result.protein_file_cif = protein_file.name
        result.protein_file_fasta = protein_file.name.split(".")[0] + ".fasta"
    except ValueError:
        os.remove(protein_file.name)
    return result


if __name__ == "__main__":

    print(
        find_contact_hybrid_ligand_protein(
            sys.argv[1], int(sys.argv[2]), sys.argv[3].split(",")
        )
    )
    print(find_contact_ion(sys.argv[1], int(sys.argv[2]), sys.argv[3].split(",")))
