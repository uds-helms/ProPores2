import os
import platform
import argparse as ap
import subprocess as sp
import multiprocessing as mp


# PROPORES call function
def propores(args):
    # cast numeric input to string
    args = [str(arg) for arg in args]
    # run PROPORES for the given PDB file
    try:
        sp.run(args, shell=platform.system() == 'Windows', check=True, stdout=sp.DEVNULL)
    # if an error occurred, write the input file name into the log file
    except sp.CalledProcessError as e:
        return


# the script needs to be wrapped in this guard for multiprocessing to work on Windows
if __name__ == '__main__':
    # COMMAND LINE PARSING
    parser = ap.ArgumentParser(description='Run PROPORES 2.0 on a directory of one or more PDB files. PDB files '
                                           'contained in sub-directories are included as well. The directory '
                                           'hierarchy in the input directory is maintained in the result directory.')
    parser.add_argument('input-directory', type=str, help='path to a directory with one or more PDBs')
    parser.add_argument('output-directory', type=str, help='path to to the output directory')
    parser.add_argument('propores', type=str, help='path to the PROPORES executable')
    parser.add_argument('--cores', type=int, default=1, help='number of cores to use in parallel')
    parser.add_argument('--args', nargs=ap.REMAINDER, default=[], help='command line arguments that should get passed '
                                                                       'to PROPORES')
    args = parser.parse_args()

    # SCRIPT
    # generate input file list and output (sub-directory)
    file_list = []
    # iterate over all files and potential sub-directories in the PDB input directory
    for root, dirs, files in os.walk(args.input_directory):
        for file in files:
            # skip files that are not marked as PDB files
            if not file.lower().endswith('.pdb'):
                continue
            # PDB file path
            file_path = os.path.join(root, file)
            # sub-directory of the PDB file
            sub_directory = os.path.join(output_directory_path, root.split(PDB_file_directory_path)[1].strip(os.sep))

            file_list.append((file_path, sub_directory))

    # create a list of arguments for the propores(log, pdb_file_path, args) function above
    # ADD ADDITIONAL PROPORES ARGUMENTS TO THE END OF THE TUPLE
    propores_args = tuple(args.args)
    tasks = [(args.propores, '-i', pdb_path, '-o', output_sub_directory_path) + propores_args
             for pdb_path, output_sub_directory_path in file_list]

    # run the task list in parallel with the given number of cores
    with mp.Pool(processes=args.cores, maxtasksperchild=1) as pool:
        pool.starmap(propores, tasks)
