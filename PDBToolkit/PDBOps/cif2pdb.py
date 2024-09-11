"""
Convert CIF to PDB.
"""
import os
from Bio import PDB
import logging
import argparse
from multiprocessing import Pool
from renumber_atom import renumber_atom

logging.basicConfig(level=logging.INFO)


def cif_to_pdb(input_path, output_path, renumber = False):
    directory = os.path.dirname(output_path)
    os.makedirs(directory, exist_ok=True)

    parser = PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure('structure', input_path)
    
    if renumber:
        renumber_atom(structure, output_path)
    else:
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_path)


def cif_to_pdb_in_parallel(input_dir, output_dir, renumber = False, n_cpu = 1):
    os.makedirs(output_dir, exist_ok=True)
    total_args = []
    for filename in os.listdir(input_dir):
        if filename.lower().endswith('.cif'):
            input_path = os.path.join(input_dir, filename)
            output_filename = os.path.splitext(filename)[0] + '.pdb'
            output_path = os.path.join(output_dir, output_filename)
            total_args.append((input_path, output_path, renumber))

    with Pool(n_cpu) as pool:
        pool.starmap(cif_to_pdb, total_args)


def main(args):
    input_path = os.path.abspath(args.input_path)
    output_path = os.path.abspath(args.output_path)
    if os.path.isfile(input_path):
        cif_to_pdb(input_path, output_path, args.renumber)
    else:
        cif_to_pdb_in_parallel(input_path, output_path, args.renumber, args.n_cpu)
    logging.info("Done.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert CIF files to PDB files.')
    parser.add_argument('input_path', type=str, help='Path to the input CIF file or directory.')
    parser.add_argument('output_path', type=str, help='Path to the output PDB file or directory.')
    parser.add_argument('--renumber', action='store_true', help='Renumber atoms in the structure.')
    parser.add_argument('--n_cpu', type=int, default=1, help='Number of CPUs to use for parallel processing.')
    args = parser.parse_args()

    print("-----------------------------------------------------------------------------", flush=True)
    print("User settings:", flush=True)
    for key, value in vars(args).items():
        print(f"{key}: {value}", flush=True)
    print("-----------------------------------------------------------------------------", flush=True)

    main(args)
