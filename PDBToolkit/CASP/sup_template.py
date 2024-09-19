"""
Superpose models to template and calculate tmscore.
"""
import os
import re
import subprocess
import pandas as pd
import logging
from multiprocessing import Pool
import argparse

from PDBToolkit.config import USALIGN_PATH

logging.basicConfig(level=logging.INFO)


def run_usalign(model, reference, output_prefix = None, extra_args = None):
    command = [
        USALIGN_PATH,
        model, 
        reference
    ]
    if output_prefix:
        command.extend(["-o", output_prefix])
    if extra_args:
        command.extend(extra_args)
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        logging.error(f"Error processing {model}: {result.stderr.strip()}")
        return None
    
    tmscore_pattern = re.compile(r'TM-score\s*=\s*([0-9.]+)')
    tmscore = float(tmscore_pattern.findall(result.stdout)[1])
    logging.info(f"TM-score for {model}: {tmscore}")

    return tmscore


def wrapper(model, reference, output_prefix, extra_args):
    tmscore = run_usalign(model, reference, output_prefix, extra_args)
    return model, tmscore


def process_in_parallel(model_dir, reference_file, sup_dir = None, extra_args = None, n_cpu = 1):
    model_list = [os.path.join(model_dir, model) for model in os.listdir(model_dir) if model.endswith('.pdb')]
    reference_list = [reference_file for model in model_list]
    if sup_dir:
        output_prefix_list = [os.path.join(sup_dir, os.path.basename(model).replace(".pdb", "_sup")) for model in model_list]
    else:
        output_prefix_list = [None for model in model_list]
    extra_args_list = [extra_args for model in model_list]
    total_args = zip(model_list, reference_list, output_prefix_list, extra_args_list)
    
    with Pool(n_cpu) as pool:
        results = pool.starmap(wrapper, total_args)
    
    logging.info(f"Processed {len(results)} models.")
    
    return dict(results)


def main(args):
    model_dir = os.path.abspath(args.model_dir)
    reference_file = os.path.abspath(args.reference)
    if args.sup_dir:
        sup_dir = os.path.abspath(args.sup_dir)
        os.makedirs(sup_dir, exist_ok=True)
    else:
        sup_dir = None
    tmscore_dict = process_in_parallel(model_dir, reference_file, sup_dir, args.extra_args, args.n_cpu)
    tmscore_dict = {os.path.basename(k): v for k, v in tmscore_dict.items()}
    tmscore_df = pd.DataFrame({"model": tmscore_dict.keys(), "tmscore": tmscore_dict.values()})
    tmscore_df.to_csv(os.path.join(sup_dir, "tmscore.csv"), index=False, sep="\t")
    
    for filename in os.listdir(sup_dir):
        if filename.endswith(".pml"):
            file_path = os.path.join(sup_dir, filename)
            os.remove(file_path)

    logging.info(f"Results saved to {os.path.join(sup_dir, 'tmscore.csv')}")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Superpose models to template and calculate TM-score.")
    parser.add_argument("model_dir", help="Directory containing model files to process.")
    parser.add_argument("reference", help="Reference PDB file.")
    parser.add_argument("--sup_dir", help="Directory to save the output PDB files.")
    parser.add_argument('--extra_args', nargs='*', default=None, 
                        help='Additional arguments for USalign.')
    parser.add_argument("--n_cpu", type=int, default=1, help="Number of CPUs to use.")
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)

    main(args)
