import os
import subprocess
import sys
import traceback
from tempfile import NamedTemporaryFile
import Bio
import numpy as np
from tqdm import tqdm
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

def generate_possible_chain_names(chain_to_skip):
    output_list = []
    for asci_code_low in range(97, 123):
        if ord(asci_code_low) not in chain_to_skip:
            output_list.append(ord(asci_code_low))
    for asci_code_capital in range(65, 91):
        if ord(asci_code_capital) not in chain_to_skip:
            output_list.append(ord(asci_code_capital))
    output_list.extend([str(i) for i in range(0, 10)])


if __name__ == "__main__":
    structure, molecule_chains = calc(sys.argv[1])
    parser = MMCIFParser(QUIET=False)
    io = MMCIFIO()
    io.set_structure(parser.get_structure("str", sys.argv[1]))
    
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    lock = Lock()
    jobs = []

    residue_dictionary = dict()
    in_contact = set()

    for model in structure:
            io.set_structure(model)
            for rna in molecule_chains[model.id]["RNA"]:
                for dna_prot in (
                    molecule_chains[model.id]["DNA"] + molecule_chains[model.id]["Protein"]
                ):
                    os.mkdir(f"/tmp/{rna}_{dna_prot}_processing")
                    os.chdir(f"/tmp/{rna}_{dna_prot}_processing")

                    try:
                        with NamedTemporaryFile(suffix=f"_{dna_prot}.pdb") as temp_file_pdb:
                            with NamedTemporaryFile(
                                suffix=f"_{dna_prot}.cif"
                            ) as temp_file_cif:
                                io.save(temp_file_cif.name, ChainsSelect([rna, dna_prot]))
                                struct_dict = MMCIF2Dict(parser.get_structure("str", temp_file_cif.name))
                                temp_list = np.array(struct_dict['_atom_site.auth_asym_id'])
                                struct_dict['_atom_site.auth_asym_id']
                                if dna_prot!='R' and rna!='D':
                                    temp_list=np.where(temp_list==rna,'R',np.where(temp_list==dna_prot,'D',temp_list))
                                else:
                                    temp_list=np.where(temp_list==rna,'RNA',np.where(temp_list==dna_prot,'DNA_PROT',temp_list))
                                    temp_list=np.where(temp_list=='RNA','R',np.where(temp_list=='DNA_PROT','D',temp_list))

                                
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
                                    + f" -f /opt/filter-hbonds.awk < {temp_file_pdb.name.split('.')[0]}.hb2",
                                    stdout=subprocess.PIPE,
                                    shell=True,
                                )

                                result = filter_hydrogen_bound.communicate()[0].decode(
                                    "utf-8"
                                )

                                p1 = subprocess.Popen(
                                    f"awk -vrna_chains='{','.join(rna)}' -f /opt/filter-residues.awk < '{temp_file_pdb.name.split('.')[0]}.hb2' ",
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
                                print(result)
                                print(p5.communicate()[0].decode("utf-8"))

                    except Bio.PDB.PDBExceptions.PDBIOException:
                        print(traceback.print_exc())
    print(res)
    print("")
    print(res_n)
