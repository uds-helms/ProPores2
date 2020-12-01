import os
import platform
import argparse as ap
import subprocess as sp
from typing import List
import multiprocessing as mp


# PROPORES call function
def run_propores(args: List[str]):
    # cast numeric input to string
    args = [str(arg) for arg in args]
    # run PROPORES for the given PDB file
    try:
        sp.run(args, shell=platform.system() == 'Windows', check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
        print('+', end='')
    # if an error occurred, write the input file name into the log file
    except sp.CalledProcessError as e:
        print(e)
        return


def batch(in_dir: str, out_dir: str, propores_path: str, cores: int, args: List[str]):
    # generate input file list and output (sub-directory)
    file_list = []
    # iterate over all files and potential sub-directories in the PDB input directory
    for root, dirs, files in os.walk(in_dir):
        for file in files:
            # skip files that are not marked as PDB files
            if not file.lower().endswith('.pdb'):
                continue
            # PDB file path
            file_path = os.path.abspath(os.path.join(root, file))
            # sub-directory of the PDB file
            sub_directory = os.path.abspath(os.path.join(out_dir, root.split(in_dir)[1].strip(os.sep)))

            file_list.append((file_path, sub_directory))

    # create a list of arguments for the propores(log, pdb_file_path, args) function above
    # ADD ADDITIONAL PROPORES ARGUMENTS TO THE END OF THE TUPLE
    propores_args = tuple(args)
    tasks = [(propores_path, '-i', pdb_path, '-o', output_sub_directory_path) + propores_args
             for pdb_path, output_sub_directory_path in file_list]

    print('PROPORES 2.0 Batch Processing')
    print('PDB files to process: {0:,}'.format(len(tasks)))

    # run the task list in parallel with the given number of cores
    with mp.Pool(processes=cores, maxtasksperchild=1) as pool:
        pool.map(run_propores, tasks)

    print('\nFinished!')


# the script needs to be wrapped in this guard for multiprocessing to work on Windows
if __name__ == '__main__':
    # command line parsing
    parser = ap.ArgumentParser(description='Run PROPORES 2.0 on a directory of one or more PDB files. PDB files '
                                           'contained in sub-directories are included as well. The directory '
                                           'hierarchy in the input directory is maintained in the result directory.')
    parser.add_argument('input', type=str, help='path to a directory with one or more PDBs')
    parser.add_argument('output', type=str, help='path to to the output directory')
    parser.add_argument('propores', type=str, help='path to the PROPORES executable')
    parser.add_argument('--cores', type=int, default=1, help='number of cores to use in parallel')
    parser.add_argument('--args', nargs=ap.REMAINDER, default=[], help='command line arguments that should get passed '
                                                                       'to PROPORES')
    cmd_args = parser.parse_args()

    cmd_args.propores = os.path.abspath(cmd_args.propores)

    # executing batch processing
    batch(cmd_args.input, cmd_args.output, cmd_args.propores, cmd_args.cores, cmd_args.args)

