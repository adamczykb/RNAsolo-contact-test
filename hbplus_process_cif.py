from enum import Enum
import os
import shutil
import subprocess
import time
import numpy as np
import os
import pandas as pd
import io as inputoutput
import uuid
from typing import List, Tuple
from tempfile import NamedTemporaryFile

import requests
from contact_utils import (
    DNA_DICT,
    PROTEIN_DICT,
    RNA_DICT,
    ChainsSelect,
    ProcessingException,
)
from classify_chain_molecule import calc
from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from multiprocessing import Lock, Pool


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


class MoleculeType(Enum):
    DNA = 0
    RNA = 1
    LIGAND = 2
    ION = 3
    PROTEIN = 4


def check_molecule(residue):
    hetflag, resseq, icode = residue.get_id()
    if hetflag != " " and hetflag != "W" and len(residue) > 1:
        return MoleculeType.LIGAND

    if hetflag != " " and hetflag != "W" and len(residue) == 1:
        return MoleculeType.ION

    if residue.resname in DNA_DICT:
        return MoleculeType.DNA

    if residue.resname in RNA_DICT:
        return MoleculeType.RNA

    if residue.resname in PROTEIN_DICT:
        return MoleculeType.PROTEIN


def init_processing_locks(lfile):
    global lockfile
    lockfile = lfile


def process_append_result(
    rna_residue, dna_prot_lig_residue, output_file_path, stucture_file_path
):
    """
    Single instance process function for async multiprocessing
    [Necessary locks]
    - lockpymol - lock access to pymol cif -> pdb converter
    - lockfile - lock access to file descriptor for writing

    Parameters:
            rna_residue (str): Residue of RNA to find hydrogen bond
            dna_prot_lig_residue (str): Residue of DNA/Protein/Ligand to find hydrogen bond
            output_file_path (str): Global path to store tsv result
            stucture_file_path (str): Global path to cif file

    Returns:
            local_in_contact (set): Set of residues in contact both RNA and DNA/PROT
    """
    in_contact_description = dict()
    io = MMCIFIO()
    parser = MMCIFParser(QUIET=True)
    io.set_structure(parser.get_structure("str", stucture_file_path)[0])
    op_dir = f"/tmp/{str(uuid.uuid4()).split('-',maxsplit=1)[0]}_{rna_residue}_{dna_prot_lig_residue}_processing"
    shutil.rmtree(op_dir, ignore_errors=True)
    os.mkdir(op_dir)
    os.chdir(op_dir)
    with NamedTemporaryFile(suffix=f"_{dna_prot_lig_residue}.pdb") as temp_file_pdb:
        with NamedTemporaryFile(
            suffix=f"_{dna_prot_lig_residue}.cif", delete=False
        ) as temp_file_cif:
            io.save(
                temp_file_cif.name, ChainsSelect([rna_residue, dna_prot_lig_residue])
            )

            struct_dict = MMCIF2Dict(temp_file_cif.name)
            temp_list = np.array(struct_dict["_atom_site.auth_asym_id"])
            struct_dict["_atom_site.label_seq_id"] = struct_dict[
                "_atom_site.auth_seq_id"
            ]

            if dna_prot_lig_residue != "R" and rna_residue != "D":
                temp_list = np.where(
                    temp_list == rna_residue,
                    "R",
                    np.where(temp_list == dna_prot_lig_residue, "D", temp_list),
                )
            else:
                temp_list = np.where(
                    temp_list == rna_residue,
                    "RNA",
                    np.where(temp_list == dna_prot_lig_residue, "DNA_PROT", temp_list),
                )
                temp_list = np.where(
                    temp_list == "RNA",
                    "R",
                    np.where(temp_list == "DNA_PROT", "D", temp_list),
                )

            struct_dict["_atom_site.auth_asym_id"] = temp_list
            io.set_dict(struct_dict)
            io.save(temp_file_cif.name)
            response = requests.post(
                "http://tomek:8080",
                headers={"Content-Type": "text/plain"},
                data=open(temp_file_cif.name, "rb"),
                timeout=10000,
            )
            if response.status_code == 200:
                with open(temp_file_pdb.name, "wb") as pdb_file_output:
                    pdb_file_output.write(response.content)

            else:
                raise ProcessingException("Parsing cif error")
            try:
                os.remove("hbdebug.dat")
            except Exception:
                pass

            hydrogen_bound_extract = subprocess.Popen(
                f"hbplus {temp_file_pdb.name}",
                shell=True,
                stdout=subprocess.DEVNULL,
            )
            hydrogen_bound_extract.wait()

            if hydrogen_bound_extract.returncode != 0:
                raise ProcessingException("hbplus error")

            filter_hydrogen_bound = subprocess.Popen(
                "awk -vrna_chains="
                + ("R")
                + f" -f /opt/filter-hbonds.awk < {op_dir}/{temp_file_pdb.name.split('/')[-1].split('.')[0]}.hb2",
                stdout=subprocess.PIPE,
                shell=True,
            )

            stream_output = (
                filter_hydrogen_bound.communicate()[0].decode("utf-8").split("\n")[3:]
            )
            if len(stream_output) > 3:
                result = "\n".join(
                    [
                        ",".join(row.split())
                        for row in np.array(
                            [
                                [
                                    (result[0:9]).replace(" ", "")
                                    + " "
                                    + result[10:13]
                                    + " "
                                    + result[14:23].replace(" ", "")
                                    + " "
                                    + result[24:]
                                ]
                                for result in stream_output
                            ]
                        )[:, 0]
                    ]
                )
                df = pd.read_csv(
                    inputoutput.StringIO(result),
                    names=TSV_COLUMNS,
                ).replace(-1, np.nan)
                structure_to_analyze = parser.get_structure("str", stucture_file_path)[
                    0
                ]
                for index, row in df.iterrows():

                    if row["donor"][0] == "R" and row["acceptor"][0] != "R":
                        key = f'{rna_residue}.{int("".join(row["donor"][1:]).split("-", maxsplit=1)[0])}'
                        if rna_residue not in in_contact_description.keys():
                            in_contact_description[key] = []
                        res_id = int(
                            "".join(row["acceptor"][1:]).split("-", maxsplit=1)[0]
                        )
                        try:
                            residue = structure_to_analyze[dna_prot_lig_residue][res_id]
                        except KeyError:
                            residue = structure_to_analyze[dna_prot_lig_residue][
                                (
                                    f'H_{"".join(row["acceptor"][1:]).split("-", maxsplit=1)[1]}',
                                    res_id,
                                    " ",
                                )
                            ]
                        in_contact_description[key].append(
                            (
                                f"{dna_prot_lig_residue}.{res_id}",
                                residue.resname,
                                check_molecule(residue),
                            )
                        )
                    elif row["donor"][0] != "R" and row["acceptor"][0] == "R":
                        key = f'{rna_residue}.{int("".join(row["acceptor"][1:]).split("-", maxsplit=1)[0])}'
                        if rna_residue not in in_contact_description.keys():
                            in_contact_description[key] = []
                        res_id = int(
                            "".join(row["donor"][1:]).split("-", maxsplit=1)[0]
                        )
                        try:
                            residue = structure_to_analyze[dna_prot_lig_residue][res_id]
                        except KeyError:
                            residue = structure_to_analyze[dna_prot_lig_residue][
                                (
                                    f'H_{"".join(row["donor"][1:]).split("-", maxsplit=1)[1]}',
                                    res_id,
                                    " ",
                                )
                            ]
                        in_contact_description[key].append(
                            (
                                f"{dna_prot_lig_residue}.{res_id}",
                                residue.resname,
                                check_molecule(residue),
                            )
                        )

                df["donor"] = (
                    df["donor"]
                    .str.replace(r"^R", rna_residue, regex=True)
                    .str.replace(r"^D", dna_prot_lig_residue, regex=True)
                )
                df["acceptor"] = (
                    df["acceptor"]
                    .str.replace(r"^R", rna_residue, regex=True)
                    .str.replace(r"^D", dna_prot_lig_residue, regex=True)
                )

                with lockfile:
                    with NamedTemporaryFile() as tmp_file:
                        with open(output_file_path, "ab") as otsv:
                            df.to_csv(tmp_file.name, sep="\t")
                            otsv.write(
                                "".join(
                                    [i.decode() for i in tmp_file.readlines()[1:]]
                                ).encode()
                            )
            p1 = subprocess.Popen(
                f"awk -vrna_chains='R' -f /opt/filter-residues.awk < "
                + f" {op_dir}/{temp_file_pdb.name.split('/')[-1].split('.')[0]}.hb2 ",
                stdout=subprocess.PIPE,
                shell=True,
            )
            p2 = subprocess.Popen(
                "sort -k1,1 -k2,2n ",
                stdin=p1.stdout,
                stdout=subprocess.PIPE,
                shell=True,
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
            v = p5.communicate()[0].decode("utf-8")
            try:

                local_in_contact = np.array(v.split(","))
                local_in_contact = np.delete(
                    local_in_contact, np.where(local_in_contact == "")
                )
                local_in_contact = np.char.replace(
                    local_in_contact, "D", f'{dna_prot_lig_residue}.'
                )
                local_in_contact = np.char.replace(local_in_contact, "R",f'{rna_residue}.' )
            except:
                pass
    shutil.rmtree(op_dir, ignore_errors=True)
    return local_in_contact, in_contact_description


def get_hbplus_result_for_large_structure(
    cif_path: str, rna_chains: List[str]
) -> Tuple[str, set]:
    """
    Process structure given cif format by parsing each pair [RNA,DNA/PROT]

    Parameters:
            cif_path (str): Global path to cif file
            rna_chains List[str]: List of RNA chains to analyze

    Returns:
            tsv_file_path (str): Hydrogen bonds tsv file path
            residues_in_contact (set): Set of residues in contact both RNA and DNA/PROT
    """
    structure, molecule_chains = calc(cif_path)

    lock_hb = Lock()
    jobs = []

    in_contact = set()
    in_contact_desc = dict()
    output_file = NamedTemporaryFile(
        suffix=f"{cif_path.split('/')[-1]}.tsv", delete=False
    )
    df = pd.DataFrame(
        columns=TSV_COLUMNS,
    )
    df.to_csv(output_file.name, sep="\t")
    with Pool(
        processes=14,
        initializer=init_processing_locks,
        initargs=(lock_hb,),
    ) as pool:
        for rna in rna_chains:
            for dna_prot in (
                molecule_chains[list(structure.get_models())[0].id]["DNA"]
                + molecule_chains[list(structure.get_models())[0].id]["Protein"]
                + molecule_chains[list(structure.get_models())[0].id]["Ligand"]
            ):
                p = pool.apply_async(
                    process_append_result,
                    args=(rna, dna_prot, output_file.name, cif_path),
                )
                jobs.append(p)
        while True:
            time.sleep(1)
            # catch exception if results are not ready yet
            try:
                ready = [result.ready() for result in jobs]
                successful = [result.successful() for result in jobs]
            except Exception:
                continue
            # exit loop if all tasks returned success
            if all(successful):
                break
            # raise exception reporting exceptions received from workers
            if all(ready) and not all(successful):
                raise Exception(
                    f"Workers raised following exceptions {[result._value for result in jobs if not result.successful()]}"
                )
        for p in jobs:
            res = p.get()
            in_contact.update(res[0])
            for key, value in res[1].items():
                if key in in_contact_desc:
                    in_contact_desc[key].extend(value)
                else:
                    in_contact_desc[key] = value

        df = pd.read_csv(output_file.name, sep="\t")
        df.drop(["hydrogen_bonds_no", "Unnamed: 0"], axis=1, inplace=True)
        df.reset_index(inplace=True)
        df.rename(columns={"index": "hydrogen_bonds_no"}, inplace=True)
        df["hydrogen_bonds_no"] += 1
        df.to_csv(output_file.name, sep="\t", index=False)
        return output_file.name, in_contact, in_contact_desc
