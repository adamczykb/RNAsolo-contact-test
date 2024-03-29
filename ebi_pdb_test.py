import os
import shutil
import subprocess
import sys
import time
import traceback
from tempfile import NamedTemporaryFile
import Bio
import numpy as np
from find_contact import ChainsSelect, ModelsSelect, ProcessingException
from show_chain_molecule_type_CIF import calc
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBIO
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pymol import cmd
import multiprocessing
from multiprocessing import Lock, Pool, Process
import os
import pandas as pd
import io as inputoutput

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


def process_append_result(rna, dna_prot, output_file_path, stucture_file_path):
    local_in_contact = []
    io = MMCIFIO()
    parser = MMCIFParser(QUIET=False)
    io.set_structure(parser.get_structure("str", stucture_file_path)[0])

    shutil.rmtree(f"/tmp/{rna}_{dna_prot}_processing", ignore_errors=True)
    os.mkdir(f"/tmp/{rna}_{dna_prot}_processing")
    os.chdir(f"/tmp/{rna}_{dna_prot}_processing")
    with NamedTemporaryFile(suffix=f"_{dna_prot}.pdb") as temp_file_pdb:
        with NamedTemporaryFile(suffix=f"_{dna_prot}.cif") as temp_file_cif:
            io.save(temp_file_cif.name, ChainsSelect([rna, dna_prot]))
            struct_dict = MMCIF2Dict(temp_file_cif.name)
            temp_list = np.array(struct_dict["_atom_site.auth_asym_id"])

            if dna_prot != "R" and rna != "D":
                temp_list = np.where(
                    temp_list == rna,
                    "R",
                    np.where(temp_list == dna_prot, "D", temp_list),
                )
            else:
                temp_list = np.where(
                    temp_list == rna,
                    "RNA",
                    np.where(temp_list == dna_prot, "DNA_PROT", temp_list),
                )
                temp_list = np.where(
                    temp_list == "RNA",
                    "R",
                    np.where(temp_list == "DNA_PROT", "D", temp_list),
                )

            struct_dict["_atom_site.auth_asym_id"] = temp_list
            io.set_dict(struct_dict)
            io.save(temp_file_cif.name, ChainsSelect([rna, dna_prot]))

            with lockpymol:
                cmd.delete("all")
                cmd.load(temp_file_cif.name, quiet=1)
                cmd.save(temp_file_pdb.name, quiet=1)
            try:
                os.remove("hbdebug.dat")
            except:
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
                + (",".join(rna))
                + f" -f /opt/filter-hbonds.awk < /tmp/{rna}_{dna_prot}_processing/{temp_file_pdb.name.split('/')[-1].split('.')[0]}.hb2",
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

                df["donor"] = (
                    df["donor"]
                    .str.replace(r"^R", rna, regex=True)
                    .str.replace(r"^D", dna_prot, regex=True)
                )
                df["acceptor"] = (
                    df["acceptor"]
                    .str.replace(r"^R", rna, regex=True)
                    .str.replace(r"^D", dna_prot, regex=True)
                )

                with lockfile:
                    with NamedTemporaryFile() as tmpFile:
                        with open(output_file_path, "ab") as otsv:
                            df.to_csv(tmpFile.name, sep="\t")
                            otsv.write(
                                "".join(
                                    [i.decode() for i in tmpFile.readlines()[1:]]
                                ).encode()
                            )
            p1 = subprocess.Popen(
                f"awk -vrna_chains='{','.join(rna)}' -f /opt/filter-residues.awk < '/tmp/{rna}_{dna_prot}_processing/{temp_file_pdb.name.split('/')[-1].split('.')[0]}.hb2' ",
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
            try:
                local_in_contact = np.array(p5.communicate()[0].decode("utf-8").split(","))
                local_in_contact = np.delete(
                    local_in_contact, np.where(local_in_contact == "")
                )
                local_in_contact = np.char.replace(local_in_contact, "D", dna_prot)
                local_in_contact = np.char.replace(local_in_contact, "R", rna)
            except:
                pass
    shutil.rmtree(f"/tmp/{rna}_{dna_prot}_processing", ignore_errors=True)
    return local_in_contact


def init(lpymol, lfile):
    global lockpymol
    global lockfile
    lockpymol = lpymol
    lockfile = lfile


if __name__ == "__main__":
    structure, molecule_chains = calc(sys.argv[1])

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    lock = Lock()
    lock_hb = Lock()
    jobs = []

    residue_dictionary = dict()
    in_contact = set()
    with NamedTemporaryFile(suffix=f"{sys.argv[1].split('/')[-1]}.tsv") as output_file:

        df = pd.DataFrame(
            columns=TSV_COLUMNS,
        )
        df.to_csv(output_file.name, sep="\t")
        
        with Pool(processes=1, initializer=init, initargs=(lock, lock_hb,)) as pool:

            for model in structure:
                # io.set_structure(model)
                for rna in molecule_chains[model.id]["RNA"]:
                    for dna_prot in (
                        molecule_chains[model.id]["DNA"]
                        + molecule_chains[model.id]["Protein"]
                    ):
                        p = pool.apply_async(
                            process_append_result,
                            args=(
                                rna,
                                dna_prot,
                                output_file.name,
                                sys.argv[1]

                            ),
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
                in_contact.update(p.get())

            print(in_contact)
            df = pd.read_csv(output_file.name, sep="\t")
            df.reset_index(inplace=True)
            df = df.rename(columns={"index": "hydrogen_bonds_no"})
            # df.to_csv(
            #     "/".join(sys.argv[1].split("/")[:-1])
            #     + sys.argv[1].split(".")[0]
            #     + ".tsv",
            #     sep="\t",
            # )
            df.to_csv('/tmp/elo.tsv')
