"""
Calculate phenix clashscore.
"""
import os
import subprocess
from multiprocessing import Pool
import json
import logging
import re
import argparse

from PDBToolkit.config import PHENIX_CLASHSCORE_PATH

logging.basicConfig(level=logging.INFO)


def calc_clashscore(file):
    command = [
        PHENIX_CLASHSCORE_PATH, file,
        'nuclear=True',
        'keep_hydrogens=True'
    ]
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        logging.error(f"Error processing {file}: {result.stderr.strip()}")
        return None
    
    match = re.search(r'clashscore\s*=\s*([\d.]+)', result.stdout)
    clashscore = float(match.group(1))
    logging.info(f"Clashscore for {file}: {clashscore}")
    
    return clashscore

def wrapper(file):
    clashscore = calc_clashscore(file)
    return file, clashscore


def process_in_parallel(file_list, output_path, n_cpu):
    with Pool(n_cpu) as pool:
        results = pool.map(wrapper, file_list)
    results = dict(results)

    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)


def main(args):
    if args.file:
        files = [args.file]
    if args.directory:
        files = [os.path.join(args.directory, f) for f in os.listdir(args.directory) if f.endswith('.pdb')]
    if args.list:
        with open(args.list, 'r') as file_list:
            files = [line.strip() for line in file_list]

    output_path = os.path.abspath(args.output_path)
    dirname = os.path.dirname(output_path)
    os.makedirs(dirname, exist_ok=True)

    process_in_parallel(files, output_path, args.n_cpu)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate phenix clashscore.')
    parser.add_argument('-f', '--file', type=str, help='Single PDB file to process.')
    parser.add_argument('-d', '--directory', type=str, help='Directory containing PDB files.')
    parser.add_argument('-l', '--list', type=str, help='File containing list of PDB files.')
    parser.add_argument('output_path', type=str, help='Path to the output file.')
    parser.add_argument('--n_cpu', type=int, default=1, help='Number of CPUs to use for parallel processing.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)

    options = [args.file, args.directory, args.list]
    if options.count(None) != 2:
        raise ValueError("You must specify exactly one of --file, --directory, or --list.")

    main(args)
