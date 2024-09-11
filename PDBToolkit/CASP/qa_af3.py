"""
Calculate QA scores and rank for CASP models.
"""
import os
import subprocess
import shutil
import zipfile
from pathlib import Path
import json
import pandas as pd
import argparse
import numpy as np
from multiprocessing import Pool
import logging
from Bio import PDB

logging.basicConfig(level=logging.INFO)
cif2pdb_path = Path(__file__).parent.parent / "PDBOps" / "cif2pdb.py"


def extract_files(zip_file, target_dir):
    with zipfile.ZipFile(zip_file, "r") as zip_ref:
        for file in zip_ref.namelist():
            if file.endswith(".cif") or ("summary" in file and file.endswith(".json")):
                zip_ref.extract(file, target_dir)

def calc_qa(directory, only_ptm = False):
    ptm_list = []
    iptm_list = []
    has_clash_list = []
    pdb_list = []
    qa_list = []
    for file in os.listdir(directory):
        name, ext = os.path.splitext(file)
        if ext == ".json":
            pdb_file = name.replace("summary_confidences", "model") + ".pdb"
            pdb_list.append(pdb_file)
            with open(os.path.join(directory, file), "r") as f:
                metrics = json.load(f)
                iptm = metrics["iptm"]
                ptm = metrics["ptm"]
                has_clash = metrics["has_clash"]
                if only_ptm:
                    qa = ptm
                else:
                    qa = iptm * 0.8 + ptm * 0.2
                has_clash_list.append(has_clash)
                iptm_list.append(iptm)
                ptm_list.append(ptm)
                qa_list.append(qa)

    data = pd.DataFrame({
        "file": pdb_list, 
        "iptm": iptm_list, 
        "ptm": ptm_list, 
        "has_clash": has_clash_list,
        "qa": qa_list, 
    })
    data["qa"] = data["qa"].map('{:.3f}'.format)
    data["ptm"] = data["ptm"].map('{:.2f}'.format)
    if only_ptm:
        data.drop(columns=["iptm"], inplace=True)
    else:
        data["iptm"] = data["iptm"].map('{:.2f}'.format)
    data.sort_values(by="qa", ascending=False, inplace=True)
    data["rank"] = [f"rank_{i}.pdb" for i in range(1, len(data) + 1)]

    return data
    
def calc_plddt(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    model = parser.get_structure('structure', pdb_file)[0]

    b_factors = [
        atom.bfactor 
        for chain in model
        for residue in chain
        for atom in residue
    ]

    return np.mean(b_factors) / 100

def calc_plddt_wrapper(file):
    plddt = calc_plddt(file)
    return file, plddt

def qa_pipeline(input_dir, output_dir, renumber = True, no_clash = False, only_ptm = False, n_cpu = 1):
    # unzip
    for file in os.listdir(input_dir):
        if file.endswith(".zip"):
            extract_files(os.path.join(input_dir, file), output_dir)
    # cif to pdb
    logging.info("Running cif2pdb.py")
    cif2pdb_args = [
        "python", str(cif2pdb_path), 
        output_dir, output_dir, 
        "--n_cpu", str(n_cpu)
    ]
    if not renumber:
        cif2pdb_args.append("--no_renumber")
    subprocess.run(cif2pdb_args, check=True)
    # qa
    data = calc_qa(output_dir, only_ptm=only_ptm)
    # plddt
    with Pool(n_cpu) as pool:
        plddts = pool.map(calc_plddt_wrapper, [os.path.join(output_dir, f) for f in data["file"].to_list()])
    plddt_dict = {os.path.basename(file): plddt for file, plddt in plddts}
    data["plddt"] = data["file"].map(f"{plddt_dict:.4f}")
    # rank
    for pdb_file, rank in zip(data["file"], data["rank"]):
        shutil.copy(os.path.join(output_dir, pdb_file), os.path.join(output_dir, rank))
    # clash
    if no_clash:
        clash_files = data[data["has_clash"] == 1.0]["rank"]
        for file in clash_files:
            os.remove(os.path.join(output_dir, file))
        data = data[data["has_clash"] == 0.0].copy()
    
    data.to_csv(os.path.join(output_dir, "qa.csv"), index=False, sep="\t")

def main(args):
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    qa_pipeline(input_dir, output_dir, not args.no_renumber, args.no_clash, args.only_ptm, args.n_cpu)
    logging.info("QA calculation completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="The directory containing the zip files to be processed")
    parser.add_argument("output_dir", help="The directory where the processed files will be saved")
    parser.add_argument('--no_renumber', action='store_true', help='Do not renumber atoms in the structure.')
    parser.add_argument('--no_clash', action='store_true', help=
                        'Do not include structures with clashes in the final ranking.')
    parser.add_argument('--only_ptm', action='store_true', help='Only calculate the ptm score.')
    parser.add_argument('--n_cpu', type=int, default=1, help='Number of CPUs to use for processing.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)
    
    main(args)
