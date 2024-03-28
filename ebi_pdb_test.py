import os
import shutil
import subprocess
import sys
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
from multiprocessing import Lock
import os
import pandas as pd
import io as inputoutput

if __name__ == "__main__":
    structure, molecule_chains = calc(sys.argv[1])
    parser = MMCIFParser(QUIET=False)
    io = MMCIFIO()
    io.set_structure(parser.get_structure("str", sys.argv[1]))

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    lock = Lock()
    lock_hb = Lock()
    jobs = []

    residue_dictionary = dict()
    in_contact = set()
    output_file = NamedTemporaryFile(suffix=f"{sys.argv[1]}.tsv", delete=False)
    output_file.write(
        """<---DONOR---> <-ACCEPTOR-->    atom                        ^               
h    n   atom  resd res      DA  || num        DHA   H-A  angle D-A-AA Bond
n    s   type  num  typ     dist DA aas  dist angle  dist       angle   num
""".encode()
    )
    for model in structure:
        io.set_structure(model)
        for rna in molecule_chains[model.id]["RNA"]:
            for dna_prot in (
                molecule_chains[model.id]["DNA"] + molecule_chains[model.id]["Protein"]
            ):

                shutil.rmtree(f"/tmp/{rna}_{dna_prot}_processing", ignore_errors=True)
                os.mkdir(f"/tmp/{rna}_{dna_prot}_processing")
                os.chdir(f"/tmp/{rna}_{dna_prot}_processing")

                # try:
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

                        with lock:
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
                            filter_hydrogen_bound.communicate()[0]
                            .decode("utf-8")
                            .split("\n")[3:]
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
                                names=[
                                    "donor",
                                    "donor_atom",
                                    "acceptor",
                                    "acceptor_atom",
                                    "distance",
                                    "atom_categories",
                                    "donor_acceptor_groups_gap",
                                    "CA_atoms_donor_acceptor_dostance",
                                    "hydrogen_donor_acceptor_angle",
                                    "hydrogen_acceptor_distance",
                                    "acceptor_hydrogen_antecedent_angle",
                                    "donor_acceptor_antecedent_angle",
                                    "#hydrogen_bonds",
                                ],
                            ).replace(-1,np.nan)
                            print(df)
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
                            df.to_csv("/tmp/zabawa.tsv", sep="\t")
                            with lock_hb:
                                output_file.write(result.encode())
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
                        local_in_contact = np.array(
                            p5.communicate()[0].decode("utf-8").split(",")
                        )
                        np.char.replace(local_in_contact, "D", dna_prot)
                        np.char.replace(local_in_contact, "R", rna)
                        in_contact.add(local_in_contact)
                shutil.rmtree(f"/tmp/{rna}_{dna_prot}_processing", ignore_errors=True)

                # except Bio.PDB.PDBExceptions.PDBIOException:
                #     print(traceback.print_exc())

    output_file.close()
    print(in_contact)
