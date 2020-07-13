# import standard or third party modules
import os
import yaml
import traceback
from shutil import which
from sys import platform
from copy import deepcopy
from tkinter import messagebox
from collections import defaultdict

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2020 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


""" HELPER FUNCTIONS """


def check_dir(name, path, files=False):
    """
    Checks if the directory exists and is a directory, and shows an error message if it is not.
    :param name: name of the directory
    :param path: directory path
    :param files: True if files in directory are to be checked as well, False otherwise
    :return: True if the directory is valid, False otherwise
    """
    path = adjust_dir_path(path)

    if not os.path.isdir(path):
        if not files:
            os.makedirs(path)
            return True
        else:
            messagebox.showerror('Error', 'The {0} directory does not exist.\n\n'
                                          'Directory path: {1}'.format(name, path))
            return False

    # if specified, check if the files in the directory are accessible for reading
    if files:
        for fn in os.listdir(path):
            if os.path.isfile(fn) and not os.access(path + fn, os.R_OK):
                messagebox.showerror('Error', 'A file in the {1} directory is not readable.\n\n'
                                              'File name: {0}\n\n'
                                              'Directory path: {2}'.format(fn, name, path))
                return False

    return True


def setup_directory(name, directory):
    """
    Sets up the specified directory if it does not exist yet and shows and error message if that fails.
    :param name: name of the directory
    :param directory: directory path
    :return: True if the directory already existed/could be created and is accessible, False otherwise
    """
    dir = adjust_dir_path(directory)

    # try to create the directory if it does not exist
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
            return True
        except OSError:
            messagebox.showerror('Error', 'The {0} directory did not exist and could not be created.\n\n'
                                          'Directory path: {1}'.format(name, dir))
            return False

    # check if it is accessible for reading
    if not os.access(dir, os.R_OK):
        messagebox.showerror('Error', 'The {0} directory is not readable.\n\n'
                                      'Directory path: {1}'.format(name, dir))
        return False

    # check if it is accessible for writing
    if not os.access(dir, os.W_OK):
        messagebox.showerror('Error', 'The {0} directory is not writable.\n\n'
                                      'Directory path: {1}'.format(name, dir))
        return False

    return True


def clean_directory(name, directory, subdirs=False):
    """
    Deletes all files and subdirectories.
    :param name: name of the directory
    :param directory: directory path
    :return: True if all files/subdirectories could be delete, False otherwise
    """
    directory = adjust_dir_path(directory)

    # do nothing if it is not a valid directory
    if not check_dir(name, directory):
        return False

    failed = set()

    # try to delete all files and subdirectories unless they are exempt
    for fn in os.listdir(directory):
        # delete files
        if os.path.isfile(directory + fn):
            try:
                os.remove(directory + fn)
            except OSError:
                failed.add(fn)
        # delete subdirectories if specified
        elif os.path.isdir(directory + fn) and subdirs:
            if not clean_directory(name, directory + fn, subdirs):
                failed.add(fn)
            else:
                os.rmdir(directory + fn)

    if failed:
        msg = 'The following files in the {0} directory could not be deleted:\n- {2}\n\n' \
              'Directory path: {1}'
        messagebox.showerror('Error', msg.format(name, directory, '\n- '.join(sorted(failed))))
        return False

    return True


def find_exe(names, path):
    """
    Tries to find the specified executable by first checking if the given path points to a valid executable. If that
    fails, try to see if the executable name on the path is in the system environment. If that also fails,
    try to see if the name is in the system environment. If all that fails, return the given input path.
    :param names: list of names name of the executable, i.e. propores, Propores, propores64
    :param path: given path of the executable
    :return: actual path to the executable
    """
    # test if the given path points to a valid executable and if so, return it
    given_exe_path = which(path, mode=os.X_OK)

    # found it
    if given_exe_path:
        return adjust_file_path(given_exe_path)

    # otherwise test if executable at the end of the given path is in the system environment and if so, return it
    ending_exe_path = which(adjust_file_path(path).split('/')[-1], mode=os.X_OK)

    # found it
    if ending_exe_path:
        return adjust_file_path(ending_exe_path)

    # otherwise test if the given name of the executable is in the system environment and if so, return it
    for name in names:
        name_exe_path = which(name, mode=os.X_OK)
        # found the executable
        if name_exe_path:
            return adjust_file_path(name_exe_path)

    return adjust_file_path(given_exe_path)


def adjust_dir_path(path, abs=False, rel=False):
    """
    Replaces backlashes with forward slashes and adds a forward slash at the end if missing.
    :param path: directory path
    :return: adjusted directory path
    """
    if not path:
        return path

    path = path.replace('\\', '/')

    if path[-1] != '/':
        path += '/'

    if abs:
        return os.path.abspath(path).replace('\\', '/')

    if rel:
        return os.path.relpath(path).replace('\\', '/')

    return path.replace('\\', '/')


def adjust_file_path(path, abs=False, rel=False):
    """
    Replaces backlashes with forward slashes.
    :param path: file path
    :return:  adjusted file path
    """
    if not path:
        return path

    if abs:
        return os.path.abspath(path).replace('\\', '/')

    if rel:
        return os.path.relpath(path).replace('\\', '/')

    return path.replace('\\', '/')


""" GUI CONFIGURATION """

class Options:
    """
    Constants for option menus.
    """
    # pore ID computation mode
    mode_autodetect = 'autodetect'                      # type: str
    mode_standalone = 'standalone'                      # type: str
    mode_ray_trace = 'ray-trace'                        # type: str
    # pore ID cavity/pore type selection
    pore_all = 'all'                                        # type: str
    pore_only = 'only pores'                                # type: str
    cavity_only = 'only cavities'                           # type: str
    # gate open difficulty modes
    difficulty_all = 'run all'                              # type: str
    difficulty_medium = 'skip medium and hard gates'        # type: str
    difficulty_hard = 'skip hard gates'                     # type: str
    # H atom processing
    h_atom_keep_all = 'keep all'                                    # type: str
    h_atom_remove_all = 'remove all'                                # type: str
    h_atom_remove_only_protein = 'remove only in ATOM records'      # type: str
    h_atom_remove_only_hetero = 'remove only in HETATM records'     # type: str
    # hetero atoms
    hetero_keep_all = 'keep all'                                    # type: str
    hetero_remove_all = 'remove all'                                # type: str
    hetero_remove_except_dummy = 'remove all except dummy atoms'    # type: str
    hetero_remove_only_dummy = 'remove only dummy atoms'            # type: str
    # preparation
    prep_both = 'prepare both'                                      # type: str
    prep_only_axis = 'only axis trace preparation'                  # type: str
    prep_only_gate = 'only gate open preparation'                   # type: str
    prep_none = 'no preparation'                                    # type: str
    prep_auto = 'autodetect'                                        # type: str


class GuiEntity:
    """
    Represents an entity in the GUI: name + label + description
    """
    def __init__(self, label='', desc=''):
        self.label = label      # type: str
        self.desc = desc        # type: str


class GuiKeys:
    """
    Constant keys for the GUI.
    """
    """ GENERAL / REQUIRED """
    gr_title = 'General'
    out_dir = 'Output directory'
    pdb_single = 'Single Protein PDB file'
    out_name = 'Out name'
    pdb_batch = 'Protein PDB batch'
    cores = 'Number of cores'
    input_selection = 'Input selection'
    h_atoms = 'H-atoms'
    keep_alt = 'Keep alternative atom locations'
    hetero = 'HETATM entries (hetero atoms)'
    skip_non_std = 'Skip non-standard amino acids in ATOM records'
    exe = 'PROPORES executable'

    """ PORE ID """
    pore_title = 'Pore ID'
    res = 'Resolution'
    solvent = 'Solvent radius'
    probe = 'Probe radius'
    volume = 'Volume threshold'
    preparation = 'Preparation'
    computation_mode = 'Computation mode'
    pore_filter = 'Pore filter'

    """ AXIS """
    axis_title = 'Axis determination'
    surface = 'Surface patch'
    axis_selection = 'Axis selection'
    axis_single = 'Axis single'
    axis_dir = 'Axis directory'

    """ GATE """
    gate_title = 'Gate Opening'
    clash = 'Clash tolerance'
    re_estimate = 'Re-estimate difficulty'
    difficulty = 'Gate difficulty'
    gate_selection = 'Gate selection'
    gate_single = 'Gate single'
    gate_dir = 'Gate directory'
    rotamer = 'Rotamer library'

    """ MENU BAR """
    # settings
    ms_title = 'Settings'
    ms_save = 'Save current settings'                   # type: str
    ms_def = 'Restore default settings'                 # type: str
    ms_save_def = 'Restore and save default settings'   # type: str
    ms_ana = 'Advanced analysis settings'               # type: str
    ms_ngs = 'Advanced NGS settings'                    # type: str
    # help
    mh_title = 'Help'                       # type: str
    mh_man = 'Open user manual'             # type: str
    mh_inst = 'Open installation guide'     # type: str
    mh_install = 'Install NGS'              # type: str

    """ OTHERS """
    g_desc = 'Description'                  # type: str
    run_button = 'Run button'               # type: str


""" CONFIGURATION COLLECTIONS """


class UserConfig:
    """
    Analysis, NGS and converter input paths, settings and output directories.
    """
    def __init__(self, cfg_dict):
        self.cfg_dict = deepcopy(cfg_dict)                                                     # type: dict

        # general
        self.pdb_path = adjust_file_path(cfg_dict['General']['PDB path'])                   # type: str
        self.out_dir = adjust_dir_path(cfg_dict['General']['output directory'])             # type: str
        self.out_name = cfg_dict['General']['output name']                                  # type: str
        self.h_atom = str(cfg_dict['General']['H atoms'])                                   # type: str
        self.keep_alternative = bool(cfg_dict['General']['keep alternative locations'])     # type: bool
        self.hetero = str(cfg_dict['General']['hetero atoms'])                              # type: str
        self.skip_non_std = bool(cfg_dict['General']['skip non-std amino acids'])           # type: bool
        self.exe = adjust_file_path(cfg_dict['General']['PROPORES executable'])             # type: str
        self.pdb_batch_dir = adjust_dir_path(cfg_dict['General']['batch input'])            # type: str
        self.pdb_batch_cores = int(cfg_dict['General']['cores'])                            # type: int
        self.input_selection = int(cfg_dict['General']['input selection'])                  # type: int

        # pore ID
        self.run_id = bool(cfg_dict['Pore ID']['run'])                                      # type: bool
        self.resolution = float(cfg_dict['Pore ID']['resolution'])                          # type: float
        self.solvent = float(cfg_dict['Pore ID']['solvent radius'])                         # type: float
        self.probe = float(cfg_dict['Pore ID']['probe radius'])                             # type: float
        self.volume = float(cfg_dict['Pore ID']['volume threshold'])                        # type: float
        self.computation_mode = cfg_dict['Pore ID']['computation mode']                     # type: str
        self.prep = str(cfg_dict['Pore ID']['preparation'])                                 # type: str
        self.pore_filter = cfg_dict['Pore ID']['pore filter']                               # type: str

        # axis trace
        self.run_axis = bool(cfg_dict['Axis trace']['run'])                                 # type: bool
        self.axis_single = adjust_file_path(cfg_dict['Axis trace']['single input'])         # type: str
        self.axis_dir = adjust_dir_path(cfg_dict['Axis trace']['directory input'])          # type: str
        self.surface_patch = float(cfg_dict['Axis trace']['surface patch threshold'])       # type: float
        self.axis_from_input = bool(cfg_dict['Axis trace']['run from input'])               # type: bool
        self.axis_selection = int(cfg_dict['Axis trace']['input selection'])                # type: int

        # gate-open
        self.run_gate = bool(cfg_dict['Gate open']['run'])                                  # type: bool
        self.gate_single = adjust_file_path(cfg_dict['Gate open']['single input'])          # type: str
        self.gate_dir = adjust_dir_path(cfg_dict['Gate open']['directory input'])           # type: str
        self.clash = float(cfg_dict['Gate open']['clash tolerance'])                        # type: float
        self.reestimate = bool(cfg_dict['Gate open']['re-estimate difficulty'])             # type: bool
        self.skip_difficulty = cfg_dict['Gate open']['difficulty']                          # type: str
        self.gate_selection = int(cfg_dict['Gate open']['input selection'])                 # type: int
        self.rotamer = adjust_dir_path(cfg_dict['Gate open']['rotamer library'])            # type: str

    def update_cfg_dict(self):
        """
        Updates and returns the configuration dictionary.
        :return: updated configuration dictionary
        """
        # general
        self.cfg_dict['General']['PDB path'] = adjust_file_path(self.pdb_path)
        self.cfg_dict['General']['output directory'] = adjust_dir_path(self.out_dir)
        self.cfg_dict['General']['output name'] = self.out_name
        self.cfg_dict['General']['H atoms'] = self.h_atom
        self.cfg_dict['General']['keep alternative locations'] = bool(self.keep_alternative)
        self.cfg_dict['General']['PROPORES executable'] = adjust_file_path(self.exe)
        self.cfg_dict['General']['hetero atoms'] = self.hetero
        self.cfg_dict['General']['skip non-std amino acids'] = bool(self.skip_non_std)
        self.cfg_dict['General']['batch input'] = adjust_dir_path(self.pdb_batch_dir)
        self.cfg_dict['General']['cores'] = int(self.pdb_batch_cores)
        self.cfg_dict['General']['input selection'] = int(self.input_selection)
        
        # pore ID
        self.cfg_dict['Pore ID']['run'] = bool(self.run_id)
        self.cfg_dict['Pore ID']['resolution'] = float(self.resolution)
        self.cfg_dict['Pore ID']['solvent radius'] = float(self.solvent)
        self.cfg_dict['Pore ID']['probe radius'] = float(self.probe)
        self.cfg_dict['Pore ID']['volume threshold'] = float(self.volume)
        self.cfg_dict['Pore ID']['computation mode'] = self.computation_mode
        self.cfg_dict['Pore ID']['preparation'] = self.prep
        self.cfg_dict['Pore ID']['pore types'] = self.pore_filter

        # axis trace
        self.cfg_dict['Axis trace']['run'] = bool(self.run_axis)
        self.cfg_dict['Axis trace']['single input'] = adjust_file_path(self.axis_single)
        self.cfg_dict['Axis trace']['directory input'] = adjust_dir_path(self.axis_dir)
        self.cfg_dict['Axis trace']['surface patch threshold'] = float(self.surface_patch)
        self.cfg_dict['Axis trace']['run from input'] = bool(self.axis_from_input)
        self.cfg_dict['Axis trace']['input selection'] = int(self.axis_selection)

        # gate-open
        self.cfg_dict['Gate open']['run'] = bool(self.run_gate)
        self.cfg_dict['Gate open']['single input'] = adjust_file_path(self.gate_single)
        self.cfg_dict['Gate open']['directory input'] = adjust_dir_path(self.gate_dir)
        self.cfg_dict['Gate open']['clash tolerance'] = float(self.clash)
        self.cfg_dict['Gate open']['re-estimate difficulty'] = self.reestimate
        self.cfg_dict['Gate open']['difficulty'] = self.skip_difficulty
        self.cfg_dict['Gate open']['input selection'] = int(self.gate_selection)
        self.cfg_dict['Gate open']['rotamer library'] = adjust_dir_path(self.rotamer)

        return deepcopy(self.cfg_dict)


class GuiConfig:
    def __init__(self, cfg_dict):
        self.cfg_dict = cfg_dict                                                # type: dict

        """ COLOURS """
        # foreground
        self.st_colour = cfg_dict['colours']['section title']                   # type: str
        self.sst_colour = cfg_dict['colours']['sub-section title']              # type: str
        self.dt_colour = cfg_dict['colours']['description title']               # type: str
        self.dd_colour = cfg_dict['colours']['description text']                # type: str

        # background

        """ FONTS """
        self.st_font = tuple(cfg_dict['fonts']['section title'])
        self.sst_font = tuple(cfg_dict['fonts']['sub-section title'])
        self.sl_font = tuple(cfg_dict['fonts']['status line'])
        self.el_font = tuple(cfg_dict['fonts']['entry field label'])
        self.ef_font = tuple(cfg_dict['fonts']['entry field text'])
        self.dt_font = tuple(cfg_dict['fonts']['description title'])
        self.dd_font = tuple(cfg_dict['fonts']['description text'])
        self.cb_font = tuple(cfg_dict['fonts']['check button'])
        self.tcb_font = tuple(cfg_dict['fonts']['title check button'])
        self.fb_font = tuple(cfg_dict['fonts']['file chooser button'])
        self.rb_font = tuple(cfg_dict['fonts']['run button'])
        self.om_font = tuple(cfg_dict['fonts']['option menu'])
        self.rdb_font = tuple(cfg_dict['fonts']['radiobutton'])

        """ STYLES """
        self.theme = cfg_dict['styles']['theme']                                # type: str
        self.st_style = 'SectionTitle.TLabel'                                   # type: str
        self.sst_style = 'SubSectionTitle.TLabel'                               # type: str
        self.sl_style = 'StatusLine.TLabel'                                     # type: str
        self.el_style = 'EntryLabel.TLabel'                                     # type: str
        self.dt_style = 'DescriptionTitle.TLabel'                               # type: str
        self.cb_style = 'Custom.TCheckbutton'                                   # type: str
        self.tcb_style = 'SectionTitle.TCheckbutton'                            # type: str
        self.fb_style = 'FileChooser.TButton'                                   # type: str
        self.rb_style = 'RunButton.TButton'                                     # type: str
        self.om_style = 'Options.TMenubutton'                                   # type: str
        self.rdb_style = 'Selector.TRadiobutton'                                # type: str

        """ WIDTHS """
        # entry fields
        self.fe_width = cfg_dict['widths']['file/directory entry field']        # type: int
        self.ne_width = cfg_dict['widths']['number entry field']                # type: int
        # buttons
        self.fb_width = cfg_dict['widths']['file/directory chooser button']     # type: int
        self.rb_width = cfg_dict['widths']['run button']                        # type: int
        self.rsb_width = cfg_dict['widths']['right settings button']            # type: int
        self.lsb_width = cfg_dict['widths']['left settings button']             # type: int
        # borders
        self.b_width = cfg_dict['widths']['status box']                         # type: int
        self.sb_width = cfg_dict['widths']['status box border']                 # type: int
        self.sp_width = cfg_dict['widths']['separator line width']              # type: int
        # labels/messages
        self.dd_width = cfg_dict['widths']['description text']                  # type: int
        # frames
        self.df_width = cfg_dict['widths']['description frame']                 # type: int

        """ HEIGHT """
        # frames
        self.df_height = cfg_dict['heights']['description frame']               # type: int

        """ X-MARGINS """
        # labels
        self.st_padx = tuple(cfg_dict['x-margins']['section title'])            # type: (int, int)
        self.sst_padx = tuple(cfg_dict['x-margins']['subsection title'])        # type: (int, int)
        # window borders
        self.wl_padx = cfg_dict['x-margins']['window left']                     # type: (int, int)
        self.wr_padx = cfg_dict['x-margins']['window right']                    # type: (int, int)
        self.col_padx = cfg_dict['x-margins']['column']                         # type: (int, int)
        # separators
        self.sp_padx = tuple(cfg_dict['x-margins']['separator'])                # type: (int, int)
        self.tss_padx = tuple(cfg_dict['x-margins']['text separator'])          # type: (int, int)
        # buttons
        self.rb_padx = tuple(cfg_dict['x-margins']['run button'])               # type: (int, int)
        self.lsb_padx = tuple(cfg_dict['x-margins']['left settings button'])    # type: (int, int)
        self.rsb_padx = tuple(cfg_dict['x-margins']['right settings button'])   # type: (int, int)
        self.fb_padx = tuple(cfg_dict['x-margins']['file chooser button'])      # type: (int, int)
        # check buttons
        self.tcb_padx = tuple(cfg_dict['x-margins']['text check button'])       # type: (int, int)
        self.ttcb_padx = tuple(cfg_dict['x-margins']['title check button'])     # type: (int, int)
        self.ecb_padx = tuple(cfg_dict['x-margins']['empty check button'])      # type: (int, int)
        # mouse-over description
        self.dt_padx = tuple(cfg_dict['x-margins']['description title'])        # type: (int, int)
        self.dd_padx = tuple(cfg_dict['x-margins']['description text'])         # type: (int, int)
        # entry fields
        self.el_padx = tuple(cfg_dict['x-margins']['entry field label'])        # type: (int, int)
        self.ef_padx = tuple(cfg_dict['x-margins']['entry field text'])         # type: (int, int)
        # status box
        self.sf_padx = tuple(cfg_dict['x-margins']['status box frame'])         # type: (int, int)

        """ Y-MARGINS """
        # labels
        self.st_pady = tuple(cfg_dict['y-margins']['section title'])            # type: (int, int)
        self.sst_pady = tuple(cfg_dict['y-margins']['subsection title'])        # type: (int, int)
        # window borders
        self.wb_pady = cfg_dict['y-margins']['window bottom']                   # type: (int, int)
        # separators
        self.sp_pady = tuple(cfg_dict['y-margins']['separator'])                # type: (int, int)
        self.tss_pady = tuple(cfg_dict['y-margins']['text separator'])          # type: (int, int)
        # buttons
        self.rb_pady = tuple(cfg_dict['y-margins']['run button'])               # type: (int, int)
        self.lsb_pady = tuple(cfg_dict['y-margins']['left settings button'])    # type: (int, int)
        self.rsb_pady = tuple(cfg_dict['y-margins']['right settings button'])   # type: (int, int)
        self.fb_pady = tuple(cfg_dict['y-margins']['file chooser button'])      # type: (int, int)
        # check buttons
        self.tcb_pady = tuple(cfg_dict['y-margins']['text check button'])       # type: (int, int)
        self.ttcb_pady = tuple(cfg_dict['y-margins']['title check button'])     # type: (int, int)
        self.ecb_pady = tuple(cfg_dict['y-margins']['empty check button'])      # type: (int, int)
        # mouse-over description
        self.dt_pady = tuple(cfg_dict['y-margins']['description title'])        # type: (int, int)
        self.dd_pady = tuple(cfg_dict['y-margins']['description text'])         # type: (int, int)
        # entry fields
        self.el_pady = tuple(cfg_dict['y-margins']['entry field label'])        # type: (int, int)
        self.ef_pady = tuple(cfg_dict['y-margins']['entry field text'])         # type: (int, int)
        # status box
        self.sf_pady = tuple(cfg_dict['y-margins']['status box frame'])         # type: (int, int)

        """ CONSTANTS """
        self.key = GuiKeys()                # type: GuiKeys
        self.desc = self.descriptions()     # type: dict
        # application title
        self.title = 'PROPORES'             # type: str

    def descriptions(self):
        """
        Sets up the labels and descriptions of different GUI elements.
        :return: dictionary with labels and descriptions
        """
        d = defaultdict(GuiEntity)

        """ GENERAL """
        d[self.key.gr_title] = GuiEntity(label='General',
                                         desc='This section specifies the required input and output behaviour shared '
                                              'by all PROPORES components.')
        d[self.key.out_dir] = GuiEntity(label='Result directory',
                                        desc='Directory that is going to contain a sub-directory of the current '
                                             'PROPORES run. Output files can include files for identified '
                                             'pores (pseudo-PDB), their axes (pseudo-PDB) and lining residues (tab-'
                                             'separated), as well as versions of the input protein (PDB) with '
                                             'open gates between neighbouring pores.'
                                             '\n\nPDB and pseudo-PDB files can be visualised with PDB-viewers such as '
                                             'Chimera or PyMOL.')
        d[self.key.pdb_single] = GuiEntity(label='Single PDB input',
                                           desc='Required file with protein information such as IDs, names and '
                                                'coordinates of atoms and residues.'
                                                '\n\nFormat: Protein Data Bank (.pdb)')
        d[self.key.out_name] = GuiEntity(label='Result name',
                                         desc='Name of the result sub-directory as well as prefix of all output files. '
                                              'If the field is left blank, the name of the input protein PDB file is '
                                              'automatically used. '
                                              '\n\nExample: If the input PDB file is "example_data/1EA5.pdb", '
                                              'then this defaults to "1EA5". If the result directory is "results/", '
                                              'all results will be written to "results/1EA5/".')
        d[self.key.pdb_batch] = GuiEntity(label='Batch input',
                                          desc='Directory that contains one or more protein PDB files. If this option '
                                               'is selected, PROPORES will be run on all PDBs in this directory and '
                                               'all sub-directories. The directory hierarchy of the input directory '
                                               'will be maintained in the output directory.')
        d[self.key.cores] = GuiEntity(label='Number of cores',
                                      desc='Number of cores to use when running PROPORES on a batch of PDB files in '
                                           'parallel. If this value is set to 1, all PDB files in the input directory '
                                           'will be run sequentially.')
        d[self.key.input_selection] = GuiEntity(label='Input selection',
                                                desc='Run PROPORES on a single PDB file or on a batch of several PDB '
                                                     'files in parallel.')
        d[self.key.exe] = GuiEntity(label='PROPORES executable',
                                    desc='Path to the PROPORES executable. This defaults to the standard location in '
                                         'the same directory as this graphical user interface. If it cannot be found '
                                         'there, RPOPORES will try to find it in the operating system environment. '
                                         '\n\nThis settings can be used to set the path manually, in case the '
                                         'executable was saved in a different location or was renamed, or PROPORES is '
                                         'unable to find the executable for other reasons.')
        d[self.key.h_atoms] = GuiEntity(label='Hydrogen processing',
                                        desc='This setting allows to fine-tune how hydrogen atoms should be dealt with '
                                             'when reading in the input protein file, and can distinguish between '
                                             'hydrogen atoms in PDB protein atom records (ATOM) and hydrogen atoms in '
                                             'non-protein atom records (HETATM).')
        d[self.key.hetero] = GuiEntity(label='Hetero atom processing',
                                       desc='Hetero atoms are non-protein atoms listed under "HETATM" in a PDB file. '
                                            'Some PDB files also contain dummy atoms as HETATM records, for example '
                                            'membrane dummy atoms in PDBs from the OPM database.'
                                            '\n\nThis settings can be used to specify how hetero atoms should be dealt '
                                            'with when reading in the input protein file.')
        d[self.key.keep_alt] = GuiEntity(label='Keep alternative locations',
                                         desc='Some atoms have multiple possible locations in the input protein file. '
                                              'By default, only the primary location is considered and the other '
                                              'locations are ignored. If enabled, all locations will be kept, which '
                                              'means that all alternative locations will be loaded as additional '
                                              'atoms.')
        d[self.key.skip_non_std] = GuiEntity(label='Skip non-standard amino acids',
                                             desc='Protein atoms are listed under "ATOM" in a PDB file. If enabled, '
                                                  'all "ATOM" records with non-standard resiue names (ALA, ARG, ASN, '
                                                  'ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, '
                                                  'THR, TRP, TYR, VAl) will be ignored when loading the PDB file. '
                                                  'Atoms in "HETATM" records are not affected by this.')

        # PORE ID
        d[self.key.pore_title] = GuiEntity(label='Pore ID',
                                           desc='Identifies all pores in a given protein PDB file, and for each pore '
                                                'generates a tab-separated file with lining residues and a pseudo PDB '
                                                'file with the empty space within the pore. The pseudo PDB file can be '
                                                'used for visualising the pore together with the protein.')
        d[self.key.res] = GuiEntity(label='Resolution',
                                    desc='Pore ID moves the protein on a 3D grid. The resolution specifies the '
                                         'length of a single grid box. The smaller the resolution value, the more '
                                         'fine-grained the computed pores, but also the potentially longer the '
                                         'required runtime.'
                                         '\n\nUnit: Ångström (Å)')
        d[self.key.solvent] = GuiEntity(label='Solvent radius',
                                        desc='To identify pores on the surface, Pore ID checks if enclosed areas in '
                                             'the protein are accessible for the surrounding solvent, i.e. water '
                                             'molecules with a radius of 1.4 Å. If the value is increased, hollow '
                                             'areas in the protein must also be accessible for these larger molecules '
                                             'to be classified as pores and not cavities.'
                                             '\n\nUnit: Ångström (Å)')
        d[self.key.probe] = GuiEntity(label='Probe radius',
                                      desc='After the solvent accessible surface (SAS) of the protein is computed (see '
                                           'above), Pore ID "rolls" a probe across this surface and removes shallow '
                                           'surface regions that are not deeper than the probe. The higher this value, '
                                           'the deeper surface areas need to be before they are considered to be part '
                                           'of a pore.'
                                           '\n\nUnit: Ångström (Å)')
        d[self.key.volume] = GuiEntity(label='Volume threshold',
                                       desc='Minimum volume (size) of pores and cavities. Smaller potential pores and '
                                            'cavities are not included in the results or further analysis steps.'
                                            '\n\nUnit: cubic Ångström (Å³)')
        d[self.key.preparation] = GuiEntity(label='Preparation',
                                            desc='By default, Pore ID prepares files for gate open or axis trace when '
                                                 'they are not currently enabled so that they can be run at a later '
                                                 'point without having to re-run Pore ID. This behaviour can be '
                                                 'controlled with this option.')
        d[self.key.computation_mode] = GuiEntity(label='Computation mode',
                                                 desc='Enclosed spaces within the protein or on its surface can be '
                                                      'computed in several ways. By default, Pore ID tries to guess '
                                                      'the most efficient version for each run (autodetect).'
                                                      '\n\nRay-trace: Faster for smaller proteins or higher resolution '
                                                      'values; increased RAM usage.'
                                                      '\n\nStandalone: Potentially faster for larger proteins or lower '
                                                      'resolution values; reduced RAM usage.')
        d[self.key.pore_filter] = GuiEntity(label='Pore type',
                                            desc='By default, Pore ID generates results for pores and cavities. This '
                                                 'can be restricted to only pores (accessible from the protein '
                                                 'surface) or to only cavities (completely enclosed by the protein).')

        # AXIS TRACE
        d[self.key.axis_title] = GuiEntity(label='Axis Trace',
                                           desc='Determines the axis of a pore or cavity and writes it into a pseudo '
                                                'PDB file for visualisation. Some pores have multiple axes, in which '
                                                'cases several output files are generated.')
        d[self.key.surface] = GuiEntity(label='Surface patch threshold',
                                        desc='Minimum area of a pore surface patch for it to count as '
                                             'a potential pore axis origin.'
                                             '\n\nUnit: square Ångström (Å²)')
        d[self.key.axis_selection] = GuiEntity(label='Run from prepared input',
                                               desc='Run axis trace for one or more previously computed pores or '
                                                    'cavities, which requires that Pore ID was run with "axis trace '
                                                    'preparation" before. If Pore ID is currently enabled, this '
                                                    'setting is ignored and internally generated data is used instead '
                                                    'to compute the axes of all pores or cavities.')
        d[self.key.axis_single] = GuiEntity(label='Single input file',
                                            desc='Path to a single axis trace input file from a previous Pore ID run '
                                                 'with enabled "axis trace preparation".')
        d[self.key.axis_dir] = GuiEntity(label='Input directory',
                                         desc='Path to a directory with one or more axis trace input files from a '
                                              'previous Pore ID run with enabled "axis trace preparation".')

        # GATE OPEN
        d[self.key.gate_title] = GuiEntity(label='Gate Opening',
                                           desc='Rotates the shared lining residues of two neighbouring pores in an '
                                                'effort to open the gate between them as much as possible. The result '
                                                'is a PDB file of the entire protein with the rotated residues.')
        d[self.key.clash] = GuiEntity(label='Clash tolerance',
                                      desc='Van der Waals radius overlap tolerance for rotamer generation. Lowering '
                                           'the clash tolerance can substantially speed up the gate opening runtime, '
                                           'however, too low values can also eliminate all potential rotamers.'
                                           '\n\nUnit: Ångström (Å)')
        d[self.key.re_estimate] = GuiEntity(label='Re-estimate difficulty',
                                            desc='If enabled, the estimated difficulty of potential gates is '
                                                 're-assessed after all rotamers have been generated and clashes with '
                                                 'the rest of the protein have been eliminated, thus leaving a '
                                                 'potentially much reduced number of possible rotamer combinations.')
        d[self.key.difficulty] = GuiEntity(label='Gate difficulty',
                                           desc='Some amino acids only have a small number of possible rotamers, '
                                                'while others have over 10,000. Depending on the number of residues '
                                                'in the gate and the number of possible rotamers of each gate residue, '
                                                'opening the gate can be very fast or take a long time. This setting '
                                                'allows to skip gates that are likely to take a long time to compute.')
        d[self.key.gate_single] = GuiEntity(label='Single input file',
                                            desc='Path to a single gate open input file from a previous Pore ID run '
                                                 'with enabled "gate opening preparation".')
        d[self.key.gate_dir] = GuiEntity(label='Input directory',
                                         desc='Path to a directory with one or more gate open input files from a '
                                              'previous Pore ID run with enabled "gate opening preparation".')
        d[self.key.gate_selection] = GuiEntity(label='Run from prepared input',
                                               desc='Run gate opening for one or more pair of previously computed '
                                                    'pores or cavities, which requires that Pore ID was run with '
                                                    '"gate opening preparation" before. If Pore ID is currently '
                                                    'enabled, this setting is ignored and internally generated data is '
                                                    'used instead to open all potential gates of sufficiently low '
                                                    'difficulty.')
        d[self.key.rotamer] = GuiEntity(label='Rotamer library',
                                        desc='Running gate opening requires a library of viable rotamer conformations '
                                             'for each standard amino acid. PROPORES comes with a rotamer library that '
                                             'it will try to automatically load.'
                                             '\n\nThis settings can be used to set the path manually, in case the '
                                             'library was saved in a different location or was renamed, or PROPORES is '
                                             'unable to find it for other reasons.')

        """ MENU BAR """
        # settings
        d[self.key.ms_save] = GuiEntity(label='Save')
        d[self.key.ms_def] = GuiEntity(label='Default')
        d[self.key.ms_save_def] = GuiEntity(label='Default+Save')

        """ OTHERS """
        d[self.key.g_desc] = GuiEntity(label='Description',
                                       desc='Place the mouse on a text, button or input field to obtain a brief '
                                            'description. More information can be found on github.com/Markus-'
                                            'Hollander/PROPORES, where you can also ask questions or report issues.'
                                            '\n\n\n\n\n\n')
        d[self.key.run_button] = GuiEntity(label='Run',
                                           desc='Run all selected PROPORES components and write the results into the '
                                                'specified result directory. This window might stop responding while '
                                                'PROPORES is running in the background; runtime information will '
                                                'appear in the command line window.'
                                                '\n\nWarning: If the "results/result_name/" directory already exists, '
                                                'this will delete previously generated content in that directory!')

        return d


class Config:
    def __init__(self, cfg_path):
        """
        Initialises the application configuration depending on the operating system.
        """
        """ OPERATING SYSTEM """
        # OS names
        self.windows = 'Windows'    # type: str
        self.linux = 'Linux'        # type: str
        self.mac = 'Mac OS X'       # type: str
        self.cygwin = 'Cygwin'      # type: str
        # True if the OS is supported
        self.os_supported = True    # type: bool
        # determine OS
        self.os = self._get_os()    # type: str

        """ LOAD CONFIGURATION FROM FILE """
        # path to the configuration file
        self.cfg_path = cfg_path                                                            # type: str
        # load the configuration file
        self.cfg = yaml.load(open(self.cfg_path, 'r'), Loader=yaml.Loader)                  # type: dict
        # get the user, default and GUI configuration depending on the operating system
        self.user = UserConfig(self.cfg[self.os]['User'])                                   # type: UserConfig
        self.gui = GuiConfig(self.cfg[self.os]['GUI'])                                      # type: GuiConfig
        self.options = Options()                                                            # type: Options

    def _get_os(self):
        """
        Determines the operating system running the program.
        :return: operating system
        """
        if platform.startswith('linux'):
            return self.linux
        elif platform == 'darwin':
            return self.mac
        elif platform == 'cygwin':
            return self.cygwin
        elif platform == 'win32':
            return self.windows
        else:
            self.os_supported = False
            return ''

    def restore_default(self):
        """
        Restores the default settings.
        """
        # restore the default configuration
        self.cfg[self.os]['User'] = self.cfg[self.os]['Default']
        self.user = UserConfig(self.cfg[self.os]['User'])

    def write_config(self):
        """
        Writes the configuration dictionary into the configuration.yaml file.
        """
        # write the configuration to disk
        with open(self.cfg_path, 'w') as cfg_file:
            # prevent aliases in the YAML file
            dumper = yaml.SafeDumper
            dumper.ignore_aliases = lambda dump, data: True
            cfg_file.write(yaml.dump(self.cfg, default_flow_style=False, Dumper=dumper))

    def save_config(self):
        """
        Saves the current settings to the configuration file.
        """
        # update the user configuration dict
        self.cfg[self.os]['User'] = self.user.update_cfg_dict()
        # write the configuration to disk
        self.write_config()


""" INITIALISATION """
try:
    # load the configuration
    cfg = Config('config.yaml')
    gui = cfg.gui
except Exception:
    print(os.getcwd())
    messagebox.showerror('Configuration Loading Error',
                         'An unhandled exception occurred while loading the user configuration.\n\n'
                         '{0}'.format(traceback.format_exc()))