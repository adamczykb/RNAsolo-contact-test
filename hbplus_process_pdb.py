import os
import shutil
import subprocess
import io as inputoutput

from typing import List

import numpy as np
import pandas as pd

from tempfile import NamedTemporaryFile

from contact_utils import ProcessingException

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


def get_filtered_hydrogen_bonds(
    pdb_path: str, rna_chains: List[str], output_tsv_path: str, op_dir: str
):
    hydrogen_bound_extract = subprocess.Popen(
        f"hbplus {pdb_path}", shell=True, stdout=subprocess.DEVNULL
    )
    hydrogen_bound_extract.wait()

    if hydrogen_bound_extract.returncode != 0:
        raise ProcessingException("hbplus error")

    filter_hydrogen_bound = subprocess.Popen(
        "awk -vrna_chains="
        + (",".join(rna_chains))
        + f" -f /opt/filter-hbonds.awk < {op_dir}/{pdb_path.split('/')[-1].split('.')[0]}.hb2",
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

        with NamedTemporaryFile() as tmp_file:
            with open(output_tsv_path, "ab") as otsv:
                df.to_csv(tmp_file.name, sep="\t")
                otsv.write(
                    "".join([i.decode() for i in tmp_file.readlines()[1:]]).encode()
                )
        df.drop(
            ["hydrogen_bonds_no", "Unnamed: 0"], axis=1, inplace=True, errors="ignore"
        )
        df.reset_index(inplace=True)
        df.rename(columns={"index": "hydrogen_bonds_no"}, inplace=True)
        df["hydrogen_bonds_no"] += 1
        df.to_csv(output_tsv_path, sep="\t", index=False)


def get_residue_in_contact(pdb_path: str, rna_chains: List[str], op_dir: str):
    p1 = subprocess.Popen(
        f"awk -vrna_chains='{','.join(rna_chains)}' -f /opt/filter-residues.awk < {op_dir}//{pdb_path.split('/')[-1].split('.')[0]}.hb2 ",
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

    return set(p5.communicate()[0].decode("utf-8").split(","))


def get_hbplus_result(pdb_path: str, rna_chains: List[str]):
    """
    Process structure given cif format by parsing each pair [RNA,DNA/PROT]

    Parameters:
            pdb_path (str): Global path to cif file
            rna_chains List[str]: List of RNA chains to analyze

    Returns:
            tsv_file_path (str): Hydrogen bonds tsv file path
            residues_in_contact (set): Set of residues in contact both RNA and DNA/PROT
    """
    op_dir = f"/tmp/{pdb_path.split('/')[-1].split('.')[0]}_{'_'.join(rna_chains)}_processing"
    shutil.rmtree(op_dir, ignore_errors=True)
    os.mkdir(op_dir)
    os.chdir(op_dir)

    output_file = NamedTemporaryFile(suffix=".tsv", delete=False)
    get_filtered_hydrogen_bonds(pdb_path, rna_chains, output_file.name, op_dir)
    residues_in_contact = get_residue_in_contact(pdb_path, rna_chains, op_dir)
    shutil.rmtree(op_dir, ignore_errors=True)
    return output_file.name, residues_in_contact
