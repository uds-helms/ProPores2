import os
import yaml
import argparse as ap
from collections import defaultdict
from typing import Dict, Union, Tuple, List


# extract a value from a dictionary if the key exists, return a default value if it does not
def get_value(dictionary: Dict[str, str], key: str) -> str:
    if key in dictionary:
        return str(dictionary[key])
    return '?'


# get the parameters and log dictionaries from a YAML configuration file
def get_log(file_path: str) -> Tuple[Union[Dict, str], Union[Dict, str], str]:
    with open(file_path, 'r') as log_file:
        d = yaml.load(log_file, Loader=yaml.SafeLoader)                                 # type: Dict

    # extract the PDB name if the parameter log is given
    log_dict = d['log'] if 'log' in d else dict()                                       # type: Dict
    param_dict = d['parameters'] if 'parameters' in d else dict()                       # type: Dict
    pdb_name = get_value(param_dict, 'PDB name')                                        # type: str

    return param_dict, log_dict, pdb_name


# compute the run status
def run_state(start: str, end: str, total_runtime: str) -> str:
    if start == '?':
        return 'NOT STARTED'
    if start != '?' and end == '?':
        return 'IN PROGRESS'
    return total_runtime


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
                name = '.'.join(file.split('.')[:-1])
                runs[name]['PDB.path'] = os.path.join(root, file)
                runs[name]['PDB.name'] = name

        print('PDBs in the input directory: {0:,}'.format(len(runs)))

    for root, dirs, files in os.walk(args.directory):
        # skip sub-directory if it does not belong to a PROPORES run
        if 'axes_trace_log.yaml' not in files and 'gate_open_log.yaml' not in files and 'pore_ID_log.yaml' not in files:
            continue

        # extract pore ID information
        if 'pore_ID_log.yaml' in files:
            # load the parameter and log file into a dictionary
            params, log, name = get_log(os.path.join(root, 'pore_ID_log.yaml'))

            # there is some error with the run if the parameters are not given
            if name == '?' or params == '?':
                continue

            runs[name]['PDB.name'] = name
            runs[name]['PDB.path'] = os.path.abspath(params['PDB path']) if 'PDB path' in params else '?'
            runs[name]['output.directory'] = os.path.abspath(root)
            runs[name]['ID.skip.H'] = get_value(params, 'skip hydrogen atoms')
            runs[name]['ID.skip.HETATM'] = get_value(params, 'skip hetero atoms')
            runs[name]['ID.skip.non.std.amino.acids'] = get_value(params, 'skip non-standard amino acids in ATOM records')
            runs[name]['ID.keep.alternative.locations'] = get_value(params, 'keep alternative atom locations')
            runs[name]['ID.run.axis.prep'] = get_value(params, 'run axes trace preparation')
            runs[name]['ID.run.gate.prep'] = get_value(params, 'run gate open preparation')
            runs[name]['ID.resolution'] = get_value(params, 'resolution')
            runs[name]['ID.solvent.radius'] = get_value(params, 'solvent radius')
            runs[name]['ID.probe.radius'] = get_value(params, 'probe radius')
            runs[name]['ID.volume.threshold'] = get_value(params, 'volume threshold')
            runs[name]['ID.selected.computation.mode'] = get_value(params, 'computation mode')
            runs[name]['ID.filter'] = get_value(params, 'filter')
            runs[name]['ID.start.date'] = get_value(log, 'start time')
            if 'PDB parsing stats' in log:
                runs[name]['atoms'] = get_value(log['PDB parsing stats'], 'atoms')
                runs[name]['removed.atoms'] = get_value(log['PDB parsing stats'], 'total skipped atoms')
            runs[name]['ID.n.grid.boxes'] = get_value(log, 'number of grid boxes')
            runs[name]['ID.atom.pairs'] = get_value(log, 'atom pairs')
            runs[name]['ID.used.computation.mode'] = get_value(log, 'used computation mode')
            runs[name]['ID.identified.pores'] = get_value(log, 'identified pores')
            runs[name]['ID.end.date'] = get_value(log, 'end time')
            runs[name]['ID.total.runtime'] = get_value(log, 'total runtime')

        # extract axis trace information
        if 'axes_trace_log.yaml' in files:
            # load the parameter and log file into a dictionary
            params, log, name = get_log(os.path.join(root, 'axes_trace_log.yaml'))

            # there is some error with the run if the parameters are not given
            if name == '?' or params == '?':
                continue

            runs[name]['PDB.name'] = name
            runs[name]['PDB.path'] = os.path.abspath(params['PDB path']) if 'PDB path' in params else '?'
            runs[name]['output.directory'] = os.path.abspath(root)
            runs[name]['axis.enabled'] = str(True)
            runs[name]['axis.surface.patch.threshold'] = get_value(params, 'surface patch threshold')
            runs[name]['axis.start.date'] = get_value(log, 'start time')
            runs[name]['ID.identified.pores'] = get_value(log, 'pores')
            runs[name]['axis.end.date'] = get_value(log, 'end time')
            runs[name]['axis.total.runtime'] = get_value(log, 'total runtime')

        # extract gate open information
        if 'gate_open_log.yaml' in files:
            # load the parameter and log file into a dictionary
            params, log, name = get_log(os.path.join(root, 'gate_open_log.yaml'))

            # there is some error with the run if the parameters are not given
            if name == '?' or params == '?':
                continue

            runs[name]['PDB.name'] = name
            runs[name]['PDB.path'] = os.path.abspath(params['PDB path']) if 'PDB path' in params else '?'
            runs[name]['output.directory'] = os.path.abspath(root)
            runs[name]['gate.enabled'] = str(True)
            runs[name]['gate.perturb.value'] = get_value(params, 'perturb value')
            runs[name]['gate.clash.tolerance'] = get_value(params, 'clash tolerance')
            runs[name]['gate.difficulty.threshold'] = get_value(params, 'gate difficulty threshold')
            runs[name]['gate.re.estimate.difficulty'] = get_value(params, 're-estimate gate difficulty')
            runs[name]['gate.start.date'] = get_value(log, 'start time')
            runs[name]['gates'] = get_value(log, 'gates')
            runs[name]['gate.end.date'] = get_value(log, 'end time')
            runs[name]['gate.total.runtime'] = get_value(log, 'total runtime')

    # inform the users how many runs had at least some output
    print('Runs in the output directory: {0:,}'.format(sum(1 for k, d in runs.items()
                                                           if d['ID.start.date'] != '?' or d['axis.start.date'] != '?'
                                                           or d['gate.start.date'] != '?')))

    # compute the status
    for key, run in runs.items():
        run['ID.status'] = run_state(run['ID.start.date'], run['ID.end.date'], run['ID.total.runtime'])
        run['axis.status'] = run_state(run['axis.start.date'], run['axis.end.date'], run['axis.total.runtime'])
        run['gate.status'] = run_state(run['gate.start.date'], run['gate.end.date'], run['gate.total.runtime'])
            
    # generate and output the overview
    with open(args.overview, 'w') as file:
        header = ['PDB.name', 'ID.status', 'axis.status', 'gate.status', 'atoms', 'removed.atoms', 'ID.n.grid.boxes',
                  'ID.atom.pairs', 'ID.used.computation.mode', 'ID.identified.pores', 'gates',
                  'ID.resolution', 'ID.solvent.radius', 'ID.probe.radius', 'ID.volume.threshold',
                  'ID.skip.H', 'ID.skip.HETATM', 'ID.skip.non.std.amino.acids', 'ID.keep.alternative.locations',
                  'ID.selected.computation.mode', 'ID.filter',
                  'axis.surface.patch.threshold', 'gate.perturb.value', 'gate.clash.tolerance',
                  'gate.re.estimate.difficulty', 'PDB.path', 'output.directory']
        fmt = '\t'.join('{' + str(i) + '}' for i in range(len(header))) + '\n'

        file.write(fmt.format(*header))

        for name, run in sorted(runs.items()):
            file.write(fmt.format(*[run[key] for key in header]))


