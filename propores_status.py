import os
import yaml
import platform
import argparse as ap
import subprocess as sp
import multiprocessing as mp

from collections import defaultdict


def get_value(d, k):
    if k in d:
        return str(d[k])
    return '?'

# the script needs to be wrapped in this guard for multiprocessing to work on Windows
if __name__ == '__main__':
    # COMMAND LINE PARSING
    parser = ap.ArgumentParser(description='Compute the run status (not yet started, in progress, finished) from '
                                           'an directory that contains the output of one or more PROPORES 2.0 runs. '
                                           'This includes start time, used parameters and runtime.')
    parser.add_argument('directory', type=str, help='path to a directory containing the output of one or more '
                                                    'PROPORES runs')
    parser.add_argument('overview', type=str, help='path to the file that is supposed to contain the status overview')
    parser.add_argument('-i', '--input-directory', type=str, default='',
                        help='(optional) path to the input directory that was used to determine which PROPORES runs '
                             'have not been started yet')
    args = parser.parse_args()

    # SCRIPT
    runs = defaultdict(lambda: defaultdict(lambda: '?'))
    # generate input file list if the input directory was given (and exists)
    if args.input_directory and os.path.isdir(args.input_directory):
        # iterate over all files and potential sub-directories in the PDB input directory
        for root, dirs, files in os.walk(args.input_directory):
            for file in files:
                # skip files that are not marked as PDB files
                if not file.lower().endswith('.pdb'):
                    continue
                # add PDB file path
                runs['.'.join(file.split('.')[:-1])]['PDB.path'] = os.path.join(root, file)

        print('PDBs in the input directory: {0:,}'.format(len(runs)))

    for root, dirs, files in os.walk(args.directory):
        # skip sub-directory if it does not belong to a PROPORES run
        if 'axes_trace_log.yaml' not in files and 'gate_open_log.yaml' not in files and 'pore_id_log.yaml' not in files:
            continue

        # extract pore ID information
        if 'pore_id_log.yaml' in files:
            # load the parameter and log file into a dictionary
            with open(os.path.join(root, 'pore_id_log.yaml'), 'r') as file:
                log_dict = yaml.load(file, Loader=yaml.SafeLoader)
            params = log_dict['parameters'] if 'parameters' in log_dict else '?'    # type: Dict
            log = log_dict['log'] if 'log' in log_dict else '?'                     # type: Dict
            name = get_value(params, 'PDB name')

            if name == '?' or params == '?':
                continue

            runs[name]['PDB.path'] = os.path.abspath(params['PDB path'])
            runs[name]['output.directory'] = os.path.abspath(root)
            runs[name]['skip.H'] = get_value(params, 'skip hydrogen atoms')
            runs[name]['skip.HETATM'] = get_value(params, 'skip hetero atoms')
            runs[name]['skip.non.std.amino.acids'] = get_value(params, 'skip non-standard amino acids in ATOM records')
            runs[name]['keep.alternative.locations'] = get_value(params, 'keep alternative atom locations')
            runs[name]['run.axis.prep'] = get_value(params, 'run axes trace preparation')
            runs[name]['run.gate.prep'] = get_value(params, 'run gate open preparation')
            runs[name]['resolution'] = get_value(params, 'resolution')
            runs[name]['solvent.radius'] = get_value(params, 'solvent radius')
            runs[name]['probe.radius'] = get_value(params, 'probe radius')
            runs[name]['volume.threshold'] = get_value(params, 'volume thresholds')
            runs[name]['computation.mode'] = get_value(params, 'computation mode')
            runs[name]['ID.filter'] = get_value(params, 'filter')
            runs[name]['ID.start.date'] = get_value(log, 'start time')
            if 'PDB parsing stats' in log:
                runs[name]['ID.atoms'] = get_value(log['PDB parsing stats'], 'atoms')
                runs[name]['ID.removed.atoms'] = get_value(log['PDB parsing stats'], 'total skipped atoms')
            runs[name]['ID.n.grid.boxes'] = get_value(log, 'number of grid boxes')
            runs[name]['ID.atom.pairs'] = get_value(log, 'atom pairs')
            runs[name]['ID.used.computation.mode'] = get_value(log, 'used computation mode')
            runs[name]['ID.identified.pores'] = get_value(log, 'identified pores')
            runs[name]['ID.end.date'] = get_value(log, 'end time')
            runs[name]['ID.total.runtime'] = get_value(log, 'total runtime')

        if 'axis_trace_log.yaml' in files:
            with open(os.path.join(root, 'axis_trace_log.yaml'), 'r') as file:
                log_dict = yaml.load(file, Loader=yaml.SafeLoader)
            params = log_dict['parameters']
            log = log_dict['log']
            name = params['PDB name']

            runs[name]['PDB.path'] = os.path.abspath(params['PDB path'])
            runs[name]['output.directory'] = os.path.abspath(root)
            runs[name]['surface.patch.threshold'] = get_value(params, 'surface patch threshold')
            runs[name]['axis.start.date'] = get_value(log, 'start time')
            runs[name]['ID.identified.pores'] = get_value(log, 'pores')
            runs[name]['axis.end.date'] = get_value(log, 'end time')
            runs[name]['axis.total.runtime'] = get_value(log, 'total runtime')

        if 'gate_open_log.yaml' in files:
            with open(os.path.join(root, 'gate_open_log.yaml'), 'r') as file:
                log_dict = yaml.load(file, Loader=yaml.SafeLoader)
            params = log_dict['parameters']
            log = log_dict['log']
            name = params['PDB name']

            runs[name]['PDB.path'] = os.path.abspath(params['PDB path'])
            runs[name]['output.directory'] = os.path.abspath(root)
            runs[name]['perturb.value'] = get_value(params, 'perturb.value')
            runs[name]['clash.tolerance'] = get_value(params, 'clash tolerance')
            runs[name]['difficulty.threshold'] = get_value(params, 'gate difficulty threshold')
            runs[name]['re.estimate.difficulty'] = get_value(params, 're-estimate gate difficulty')
            runs[name]['gate.start.date'] = get_value(log, 'start time')
            runs[name]['gates'] = get_value(log, 'gates')
            runs[name]['gate.end.date'] = get_value(log, 'end time')
            runs[name]['gate.total.runtime'] = get_value(log, 'total runtime')

    # create a list of arguments for the propores(log, pdb_file_path, args) function above
    # ADD ADDITIONAL PROPORES ARGUMENTS TO THE END OF THE TUPLE
