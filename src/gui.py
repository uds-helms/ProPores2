# import standard or third party modules
import os
import re
import traceback
import tkinter as tk
from shlex import quote
from shutil import which
from subprocess import Popen, run
from tkinter import ttk, filedialog, messagebox, font

# import own modules
from src.configuration import cfg, gui, adjust_dir_path, adjust_file_path, find_exe, check_dir, clean_directory, \
    setup_directory

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2017 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'


def out_warn(out_dir):
    """
    Warn the user that the old output is going to be deleted.
    """
    return messagebox.askokcancel('Warning',
                                  'This is going to delete all files in the specified output directory.\n\n'
                                  'Directory path: {0}\n\n'
                                  'Press OK to continue, or Cancel to abort.'.format(out_dir))


class Key:
    """ Can be used to give other classes support for mouse-over. """
    def __init__(self, key):
        self.key = key
        self.label = gui.desc[key].label
        self.desc = gui.desc[key].desc
        
    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour for the widget.
        :param desc: description box
        """
        self.bind('<Enter>', lambda _: desc.describe(self.label, self.desc))
        self.bind('<Leave>', desc.clear)


class TLabel(ttk.Label, Key):
    """ A title label for sections that supports mouse-over. """
    def __init__(self, master, key):
        Key.__init__(self, key)
        ttk.Label.__init__(self, master, text=self.label, style=gui.st_style)

    def place(self, row, column=0):
        self.grid(row=row, column=column, columnspan=5 - column, pady=gui.st_pady, padx=gui.st_padx, sticky=tk.W)


class STLabel(ttk.Label, Key):
    """ A title label for subsections that supports mouse-over. """
    def __init__(self, master, key):
        Key.__init__(self, key)
        ttk.Label.__init__(self, master, text=self.label, style=gui.sst_style)

    def place(self, row):
        self.grid(row=row, column=1, pady=gui.sst_pady, padx=gui.sst_padx, sticky=tk.W)


class SLabel(ttk.Label, Key):
    """ A status label that supports mouse-over. """
    def __init__(self, master, key, var, pos='l'):
        Key.__init__(self, key)
        ttk.Label.__init__(self, master, textvariable=var, style=gui.sl_style)
        if pos == 'c':
            self.configure(anchor=tk.CENTER, width=gui.b_width)
        self.pos = pos

    def place(self, row):
        if self.pos == 'l':
            self.grid(row=row, column=0, sticky=tk.W)
        elif self.pos == 'r':
            self.grid(row=row, column=1, sticky=tk.E)
        else:
            self.grid(row=row, column=0, columnspan=2, sticky=tk.EW)


class ELabel(ttk.Label, Key):
    """ A label naming an entry field that supports mouse-over. """
    def __init__(self, master, key):
        Key.__init__(self, key)
        ttk.Label.__init__(self, master, text=self.label, style=gui.el_style)

    def place(self, row):
        self.grid(row=row, column=1, pady=gui.el_pady, padx=gui.el_padx, sticky=tk.W)


class RButton(ttk.Button, Key):
    """ Main run button that supports mouse-over. """
    def __init__(self, master, key, cmd):
        Key.__init__(self, key)
        ttk.Button.__init__(self, master, text=self.label, command=cmd, width=gui.rb_width, style=gui.rb_style)

    def place(self, row):
        self.grid(row=row, column=3, columnspan=2, sticky=tk.E, pady=gui.rb_pady, padx=gui.rb_padx)


class FButton(ttk.Button, Key):
    """ File/directory chooser button that supports mouse-over. """
    def __init__(self, master, key, cmd):
        Key.__init__(self, key)
        ttk.Button.__init__(self, master, text='...', command=cmd, style=gui.fb_style, width=gui.fb_width)


class CEntry(ttk.Entry, Key):
    """ Custom entry field that supports mouse-over. """
    def __init__(self, master, key, var, width):
        Key.__init__(self, key)
        ttk.Entry.__init__(self, master, textvariable=var, width=width, font=gui.ef_font)


class LCheckbutton(ttk.Checkbutton, Key):
    """ Checkbutton with label that supports mouse-over. """
    def __init__(self, master, key, default):
        Key.__init__(self, key)
        self.var = tk.BooleanVar(value=default)
        ttk.Checkbutton.__init__(self, master, text=self.label, variable=self.var, style=gui.cb_style)

    def place(self, row):
        self.grid(row=row, column=1, columnspan=3, stick=tk.W, pady=gui.tcb_pady, padx=gui.tcb_padx)


class ERadiobutton(ttk.Radiobutton, Key):
    def __init__(self, master, key, var, value):
        Key.__init__(self, key)
        ttk.Radiobutton.__init__(self, master, text=self.label, var=var, value=value, style=gui.rdb_style)

    def show(self, row):
        self.grid(row=row, column=1, sticky=tk.W)


class TCheckbutton(ttk.Checkbutton, Key):
    """ Checkbutton with label that supports mouse-over. """
    def __init__(self, master, key, var):
        Key.__init__(self, key)
        ttk.Checkbutton.__init__(self, master, text='', variable=var, style=gui.tcb_style)

    def place(self, row):
        self.grid(row=row, column=0, stick=tk.W, pady=gui.ttcb_pady, padx=gui.ttcb_padx)


class EOptionMenu(ttk.OptionMenu, Key):
    def __init__(self, master, key, var, options):
        Key.__init__(self, key)
        ttk.OptionMenu.__init__(self, master, var, var.get(), *options, style=gui.om_style)


class Sep(ttk.Separator, Key):
    def __init__(self, master, key):
        Key.__init__(self, key)
        ttk.Separator.__init__(self, master, orient=tk.HORIZONTAL)


""" AGGREGATES """


class TSep(ttk.Frame):
    """ Frame that contains an (optional) title and a separator line that supports mouse-over. """
    def __init__(self, master, key):
        ttk.Frame.__init__(self, master)
        self.columnconfigure(1, weight=1)

        if gui.desc[key]:
            self.title = STLabel(self, key)
        else:
            self.title = None

        self.sep = Sep(self, key)

    def place(self, row, col=1):
        sep_col = 0
        if self.title:
            self.title.grid(row=0, column=0, sticky=tk.W, pady=gui.sst_pady, padx=gui.sst_padx)
            sep_col += 1
        self.sep.grid(row=0, column=sep_col, sticky=tk.EW, pady=cfg.gui.tss_pady, padx=gui.tss_padx)
        self.grid(row=row, column=col, columnspan=5 - col, sticky=tk.EW)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        if self.title:
            self.title.mouse_over(desc)
        self.sep.mouse_over(desc)


class Options:
    """
    Option menu that supports mouse-over: title - option menu
    """
    def __init__(self, master, key, default, options):
        self.var = tk.StringVar(value=default)
        self.title = ELabel(master, key)
        self.menu = EOptionMenu(master, key, self.var, options)

    def show(self, row):
        self.title.grid(row=row, column=1, sticky=tk.W, pady=gui.el_pady, padx=gui.el_padx)
        self.menu.grid(row=row, column=2, columnspan=2, sticky=tk.W, pady=gui.ef_pady, padx=gui.ef_padx)

    def mouse_over(self, desc):
        self.title.mouse_over(desc)
        self.menu.mouse_over(desc)


class NInput:
    """
    Input line for numbers that supports mouse-over: title - entry field.
    """
    def __init__(self, master, key, default, pos=False, neg=False):
        self.key = key
        self.title = ELabel(master, key)
        self.var = tk.StringVar(value=default)
        self.field = CEntry(master, key, self.var, gui.ne_width)
        self.pos = pos
        self.neg = neg

    def show(self, row):
        """
        Shows the input fields in the application window.
        """
        self.title.grid(row=row, column=1, sticky=tk.W, pady=gui.el_pady, padx=gui.el_padx)
        self.field.grid(row=row, column=2, sticky=tk.W, pady=gui.ef_pady, padx=gui.ef_padx, columnspan=2)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.title.mouse_over(desc)
        self.field.mouse_over(desc)

    def validate(self):
        """
        Checks if the value is a number.
        :return: True if the value is a number, False otherwise
        """
        try:
            val = float(self.var.get())

            if self.pos and val < 0:
                messagebox.showerror('Input Error', '{0} needs to be >= 0.'.format(gui.desc[self.key].label))
                return False

            if self.neg and val > 0:
                messagebox.showerror('Input Error', '{0} needs to be <= 0.'.format(gui.desc[self.key].label))
                return False

            return True

        except ValueError:
            messagebox.showerror('Input Error', '{0} needs to be a number.'.format(gui.desc[self.key].label))
            return False


class TextInput:
    """
    Input line for text.
    """
    def __init__(self, master, key, default):
        self.key = key
        self.title = ELabel(master, key)
        self.var = tk.StringVar(value=default)
        self.field = CEntry(master, key, self.var, gui.fe_width)

    def show(self, row, label=True):
        """
        Shows the input fields in the application window.
        """
        if label:
            self.title.grid(row=row, column=1, sticky=tk.W, pady=gui.el_pady, padx=gui.el_padx)
        self.field.grid(row=row, column=2, columnspan=2, sticky=tk.W, pady=gui.ef_pady, padx=gui.ef_padx)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.title.mouse_over(desc)
        self.field.mouse_over(desc)

    def validate(self):
        """
        Checks if the input only contains allowed characters.
        :return: True if input is valid, False otherwise
        """
        self.var.set(' '.join(self.var.get().split()).strip())

        if not re.match(r'[\w -]*$', self.var.get()):
            messagebox.showerror('Input Error', '{0} may only contain letters, numbers, underscores, dashes and '
                                                'spaces.'.format(gui.desc[self.key].label))
            return False
        return True


class FDInput:
    """
    Input line for file/directory paths that supports mouse-over: title - entry field - chooser button.
    """
    def __init__(self, master, key, default):
        self.title = ELabel(master, key)
        self.var = tk.StringVar(value=default)
        self.field = CEntry(master, key, self.var, gui.fe_width)
        self.chooser = FButton(master, key, self.ask)
        self.header = gui.desc[key].label
        self.titlefy()

    def show(self, row, label=True):
        """
        Shows the input fields in the application window.
        """
        if label:
            self.title.grid(row=row, column=1, sticky=tk.W, pady=gui.el_pady, padx=gui.el_padx)
        self.field.grid(row=row, column=2, columnspan=2, sticky=tk.W, pady=gui.ef_pady, padx=gui.ef_padx)
        self.chooser.grid(row=row, column=4, sticky=tk.E, pady=gui.fb_pady, padx=gui.fb_padx)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.title.mouse_over(desc)
        self.field.mouse_over(desc)
        self.chooser.mouse_over(desc)

    def titlefy(self):
        """
        Transforms the name to title case.
        """
        words = self.header.split()

        for i in range(len(words)):
            if not words[i].isupper():
                words[i] = words[i].title()

        self.header = ' '.join(words)

    def ask(self):
        pass

    def validate(self):
        pass


class DirInput(FDInput):
    """
    Input field for a directory.
    """
    def __init__(self, master, key, default):
        FDInput.__init__(self, master, key, default)
        self.name = gui.desc[key].label.lower().replace('directory', '').strip()

    def ask(self):
        """
        Open the directory chooser dialog.
        """
        dir_path = filedialog.askdirectory(initialdir=(os.path.abspath(self.var.get())),
                                           title='Select {0}'.format(self.header))
        if dir_path:
            self.var.set(adjust_dir_path(dir_path))

    def open(self):
        """
        Opens the directory in the OS file explorer.
        """
        path = os.path.abspath(self.var.get())

        # open command depending on the operating system
        if cfg.os == cfg.linux:
            Popen(["xdg-open", path])
        elif cfg.os == cfg.mac:
            Popen(["open", path])
        elif cfg.os == cfg.windows:
            os.startfile(path)
        elif cfg.os == cfg.cygwin:
            Popen(["xdg-open", path])

    def validate(self):
        """
        Checks if the directory exists and is a directory, and shows an error message if it is not.
        :return: True if the directory is valid, False otherwise
        """
        return check_dir(self.name, self.var.get())

    def setup(self):
        """
        Sets up the specified directory if it does not exist yet and shows and error message if that fails.
        :return: True if the directory already existed/could be created and is accessible, False otherwise
        """
        return setup_directory(self.name, self.var.get())

    def clean(self):
        """
        Deletes all files and subdirectories unless they are exempt.
        :return: True if all files/subdirectories could be delete, False otherwise
        """
        return clean_directory(self.name, self.var.get(), True)


class ReadFileInput(FDInput):
    """
    Input field for files.
    """
    def __init__(self, master, key, default):
        FDInput.__init__(self, master, key, default)
        self.name = gui.desc[key].label.lower().replace('file', '').strip()

    def ask(self):
        """
        Open the file chooser dialog.
        """
        directory, file_name = os.path.split(os.path.abspath(self.var.get()))
        file_path = filedialog.askopenfilename(initialfile=file_name,
                                               initialdir=directory,
                                               title='Select {0}'.format(self.header))
        if file_path:
            self.var.set(adjust_file_path(file_path))

    def validate(self):
        """
        Checks if the file exists, is a file and if it is readable and shows an error message if it is not.
        :return: True if the file is valid, False otherwise
        """
        path = adjust_file_path(self.var.get())

        if not os.path.isfile(path):
            messagebox.showerror('File Error', 'The {0} file does not exist.\n\n'
                                               'File: {1}'.format(self.name, path))
            return False
        if not os.access(path, os.R_OK):
            messagebox.showerror('File Error', 'The {0} file is not readable.\n\n'
                                               'File {1}'.format(self.name, path))
            return False
        return True

    def open(self):
        """
        Open the directory containing the file.
        """
        path = os.path.dirname(os.path.abspath(self.var.get()))

        # open command depending on the operating system
        if cfg.os == cfg.linux:
            Popen(["xdg-open", path])
        elif cfg.os == cfg.mac:
            Popen(["open", path])
        elif cfg.os == cfg.windows:
            os.startfile(path)
        elif cfg.os == cfg.cygwin:
            Popen(["xdg-open", path])


class ExeInput(FDInput):
    """
    Input field for executables.
    """
    def __init__(self, master, key, default, names):
        FDInput.__init__(self, master, key, default)
        self.names = names

    def ask(self):
        """
        Open the file chooser dialog.
        """
        self.var.set(adjust_file_path(filedialog.askopenfilename()))
        # check if the file path points to an executable, if not try to find it
        if not os.access(self.var.get(), os.X_OK):
            self.var.set(find_exe(self.names, self.var.get()))

    def validate(self):
        """
        Checks if a path points to a valid executable and shows an error message if not.
        :return: True if it is valid, False otherwise
        """
        path = find_exe(self.names, adjust_file_path(self.var.get()))

        if not which(path):
            msg = 'The {0} path does not point to an executable. Please select a valid executable.\n\n' \
                  'Given executable path: {1}. '
            messagebox.showerror('Executable Path Error',
                                 msg.format(self.title, path))
            return False

        self.var.set(adjust_file_path(path))

        return True


class FileSubSection:
    """
    Check button with title and separator line, as well as associated input fields that can be toggled
    via the check button.
    """
    def __init__(self, master, key, inputs, label=True):
        self.sep = TSep(master, key)
        self.inputs = inputs
        self.label = label

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.sep.mouse_over(desc)

        for button, input_field in self.inputs:
            input_field.mouse_over(desc)
            button.mouse_over(desc)

    def show(self, row):
        """
        Places the selector and the associated input fields in the application window.
        """
        self.sep.place(row, 1)

        for i in range(len(self.inputs)):
            button, input_field = self.inputs[i]
            button.show(row + i + 1)
            input_field.show(row + i + 1, self.label)


class TSelector:
    """
    Check button with title and separator line, as well as associated input fields that can be toggled
    via the check button.
    """
    def __init__(self, master, key, default, label=True):
        self.var = tk.BooleanVar(value=default)
        self.check = TCheckbutton(master, key, self.var)
        self.title = TLabel(master, key)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour.
        """
        self.check.mouse_over(desc)
        self.title.mouse_over(desc)

    def show(self, row):
        """
        Places the selector and the associated input fields in the application window.
        """
        self.check.place(row)
        self.title.place(row, 1)


class Description(ttk.Frame):
    """
    A section showing descriptions of widgets.
    """
    def __init__(self, master, row):
        ttk.Frame.__init__(self, master, width=gui.dd_width)
        # default values
        self.title_text = gui.desc[gui.key.g_desc].label
        self.desc_text = gui.desc[gui.key.g_desc].desc
        # text variables
        self.title_var = tk.StringVar(value=self.title_text)
        self.desc_var = tk.StringVar(value=self.desc_text)
        # title and description
        self.title = ttk.Label(self, textvariable=self.title_var, style=gui.dt_style)
        self.desc = ttk.Label(self, textvariable=self.desc_var)
        self.desc.configure(font=gui.dd_font, wraplength=gui.dd_width, foreground=gui.dd_colour, anchor=tk.W)

        # place in the application window
        self.show(row)

    def show(self, row):
        """
        Places the section and its widgets in the window.
        """
        self.title.grid(row=0, column=0, sticky=tk.W, pady=gui.dt_pady, padx=gui.dt_padx)
        self.desc.grid(row=1, column=0, sticky=tk.EW, pady=gui.dd_pady, padx=gui.dd_padx)
        self.grid(row=row, column=0, columnspan=5, sticky=(tk.N, tk.W, tk.E))

    def describe(self, title, desc):
        """
        Shows the specified title and description.
        """
        self.desc_var.set(desc)
        self.title_var.set(title)

    def clear(self, _):
        """
        Shows the default title and description.
        """
        self.desc_var.set(self.desc_text)
        self.title_var.set(self.title_text)


class AxisSection:
    """
    Creates and sets up a section for axis trace with input fields.
    """
    def __init__(self, master, desc, row):
        self.master = master
        self.title = TSelector(master, gui.key.axis_title, cfg.user.run_axis)
        # surface patch size input (float >= 0)
        self.surface = NInput(master, gui.key.surface, cfg.user.surface_patch, pos=True)

        # previously compute axis trace input
        # one input field for a file path and one for a directory path
        self.single_input = ReadFileInput(master, gui.key.axis_single, cfg.user.axis_single)
        self.directory_input = DirInput(master, gui.key.axis_dir, cfg.user.axis_dir)
        # mutually exclusive selection buttons
        self.selection = tk.IntVar(value=cfg.user.axis_selection)
        self.single_button = ERadiobutton(master, gui.key.axis_single, self.selection, 0)
        self.dir_button = ERadiobutton(master, gui.key.axis_dir, self.selection, 1)
        # combine the buttons and the input fields
        self.input_sel = FileSubSection(master, gui.key.axis_selection,
                                        [(self.single_button, self.single_input), (self.dir_button, self.directory_input)],
                                        label=False)

        self.place(row)
        self.mouse_over(desc)

    def place(self, row):
        """
        Places the axis trace section widgets in the main application window.
        :param row: start row of the section
        """
        self.title.show(row)
        self.surface.show(row + 1)
        self.input_sel.show(row + 2)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour to each widget of the axis trace section.
        :param desc: description box
        """
        self.title.mouse_over(desc)
        self.surface.mouse_over(desc)
        self.input_sel.mouse_over(desc)

    def save(self):
        """
        Updates the user configuration with values from the entry fields.
        """
        cfg.user.run_axis = bool(self.title.var.get())
        if self.surface.validate():
            cfg.user.surface_patch = float(self.surface.var.get())
        cfg.user.axis_dir = adjust_dir_path(self.directory_input.var.get())
        cfg.user.axis_single = adjust_file_path(self.single_input.var.get())

    def update(self):
        """
        Updates the entry fields with the current values from the user configuration.
        """
        self.title.var.set(bool(cfg.user.run_axis))
        self.surface.var.set(cfg.user.surface_patch)
        self.directory_input.var.set(adjust_dir_path(cfg.user.axis_dir))
        self.single_input.var.set(adjust_file_path(cfg.user.axis_single))
        self.master.update()


class GateSection:
    """
    Creates and sets up a section for gate opening with input fields.
    """
    def __init__(self, master, desc, row):
        self.master = master
        self.title = TSelector(master, gui.key.gate_title, cfg.user.run_gate)

        # rotamer clash tolerance (float >= 0)
        self.clash = NInput(master, gui.key.clash, cfg.user.clash, pos=True)

        # previously compute axis trace input
        # one input field for a file path and one for a directory path
        self.single_input = ReadFileInput(master, gui.key.gate_single, cfg.user.gate_single)
        self.directory_input = DirInput(master, gui.key.gate_dir, cfg.user.gate_dir)
        # mutually exclusive selection buttons
        self.selection = tk.IntVar(value=cfg.user.gate_selection)
        self.single_button = ERadiobutton(master, gui.key.gate_single, self.selection, 0)
        self.dir_button = ERadiobutton(master, gui.key.gate_dir, self.selection, 1)
        # combine the buttons and the input fields
        self.input_sel = FileSubSection(master, gui.key.gate_selection,
                                        [(self.single_button, self.single_input),
                                         (self.dir_button, self.directory_input)],
                                        label=False)

        # option menu for the difficulty settings
        self.difficulty = Options(master, gui.key.difficulty, cfg.user.skip_difficulty,
                                  [cfg.options.difficulty_all, cfg.options.difficulty_medium,
                                   cfg.options.difficulty_hard])
        # checkbutton for re-estimating the difficulty
        self.re_estimate = LCheckbutton(master, gui.key.re_estimate, cfg.user.reestimate)

        self.place(row)
        self.mouse_over(desc)

    def place(self, row):
        """
        Places the gate opening section widgets in the main application window.
        :param row: start row of the section
        """
        self.title.show(row)
        self.clash.show(row + 1)
        self.input_sel.show(row + 2)

        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 5, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)

        self.difficulty.show(row + 6)
        self.re_estimate.place(row + 7)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour to each widget of the gate opening section.
        :param desc: description box
        """
        self.title.mouse_over(desc)
        self.clash.mouse_over(desc)
        self.input_sel.mouse_over(desc)
        self.re_estimate.mouse_over(desc)
        self.difficulty.mouse_over(desc)

    def save(self):
        """
        Updates the user configuration with values from the entry fields.
        """
        cfg.user.run_gate = bool(self.title.var.get())
        if self.clash.validate():
            cfg.user.clash = float(self.clash.var.get())
        cfg.user.gate_dir = adjust_dir_path(self.directory_input.var.get())
        cfg.user.gate_single = adjust_file_path(self.single_input.var.get())
        cfg.user.skip_difficulty = self.difficulty.var.get()
        cfg.user.re_estimate = bool(self.re_estimate.var.get())

    def update(self):
        """
        Updates the entry fields with the current values from the user configuration.
        """
        self.title.var.set(bool(cfg.user.run_gate))
        self.clash.var.set(float(cfg.user.clash))
        self.single_input.var.set(adjust_file_path(cfg.user.gate_single))
        self.directory_input.var.set(adjust_dir_path(cfg.user.gate_dir))
        self.difficulty.var.set(cfg.user.skip_difficulty)
        self.re_estimate.var.set(bool(cfg.user.reestimate))
        self.master.update()


class PoreIDSection:
    """
    Creates and sets up a section for pore identification with input fields.
    """
    def __init__(self, master, desc, row):
        self.master = master
        self.title = TSelector(master, gui.key.pore_title, cfg.user.run_id)

        # grid and cluster settings
        self.resolution = NInput(master, gui.key.res, cfg.user.resolution, pos=True)
        self.solvent = NInput(master, gui.key.solvent, cfg.user.solvent, pos=True)
        self.probe = NInput(master, gui.key.probe, cfg.user.probe, pos=True)
        self.volume = NInput(master, gui.key.volume, cfg.user.volume, pos=True)

        # run flags
        self.gate_prep = LCheckbutton(master, gui.key.gate_prep, cfg.user.gate_prep)
        self.axis_prep = LCheckbutton(master, gui.key.axis_prep, cfg.user.axis_prep)

        # computation mode and pore type selection
        self.computation_mode = Options(master, gui.key.cylinder_mode, cfg.user.cylinder_mode,
                                        [cfg.options.cylinder_autodetect, cfg.options.cylinder_ray_trace,
                                          cfg.options.cylinder_standalone])
        self.pore_type = Options(master, gui.key.pore_type, cfg.user.pore_type,
                                 [cfg.options.pore_all, cfg.options.pore_only, cfg.options.cavity_only])

        self.place(row)
        self.mouse_over(desc)

    def place(self, row):
        """
        Places the analysis section widgets in the main application window.
        :param row: start row of the section
        """
        self.title.show(row)
        self.resolution.show(row + 1)
        self.solvent.show(row + 2)
        self.probe.show(row + 3)
        self.volume.show(row + 4)

        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 5, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)

        self.computation_mode.show(row + 6)
        self.pore_type.show(row + 7)

        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 8, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)

        self.gate_prep.place(row + 9)
        self.axis_prep.place(row + 10)

        # separator between the input fields and the run button
        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 11, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour to each widget of the analysis section.
        :param desc: description box
        """
        self.title.mouse_over(desc)
        self.resolution.mouse_over(desc)
        self.solvent.mouse_over(desc)
        self.probe.mouse_over(desc)
        self.volume.mouse_over(desc)
        self.gate_prep.mouse_over(desc)
        self.axis_prep.mouse_over(desc)
        self.computation_mode.mouse_over(desc)
        self.pore_type.mouse_over(desc)

    def save(self):
        """
        Updates the user configuration with values from the entry fields.
        """
        if self.resolution.validate():
            cfg.user.resolution = float(self.resolution.var.get())
        if self.solvent.validate():
            cfg.user.solvent = float(self.solvent.var.get())
        if self.probe.validate():
            cfg.user.probe = float(self.probe.var.get())
        if self.volume.validate():
            cfg.user.volume = float(self.volume.var.get())
        cfg.user.gate_prep = bool(self.gate_prep.var.get())
        cfg.user.axis_prep = bool(self.axis_prep.var.get())
        cfg.user.cylinder_mode = self.computation_mode.var.get()
        cfg.user.pore_type = self.pore_type.var.get()

    def update(self):
        """
        Updates the entry fields with the current values from the user configuration.
        """
        self.resolution.var.set(cfg.user.resolution)
        self.solvent.var.set(cfg.user.solvent)
        self.probe.var.set(cfg.user.probe)
        self.volume.var.set(cfg.user.volume)

        self.axis_prep.var.set(cfg.user.axis_prep)
        self.gate_prep.var.set(cfg.user.gate_prep)
        self.computation_mode.var.set(cfg.user.cylinder_mode)
        self.pore_type.var.set(cfg.user.pore_type)

        self.master.update()


class GeneralSection:
    """
    Creates and sets up a section for general settings with input fields and run button.
    """
    def __init__(self, master, desc, row):
        self.master = master
        self.title = TLabel(master, gui.key.gr_title)

        # PDB input file and output directory
        self.pdb = ReadFileInput(master, gui.key.pdb, cfg.user.pdb_path)
        self.out_dir = DirInput(master, gui.key.out_dir, cfg.user.out_dir)
        self.out_name = TextInput(master, gui.key.out_name, cfg.user.out_name)

        # PROPORES executable
        self.exe = ExeInput(master, gui.key.exe, cfg.user.exe, 'propores')
        self.skip_h = LCheckbutton(master, gui.key.h_atoms, cfg.user.skip_h)
        self.keep_alternative = LCheckbutton(master, gui.key.keep_alt, cfg.user.keep_alternative)

        self.mouse_over(desc)
        self.place(row)

    def place(self, row):
        """
        Places the general section widgets in the main application window.
        :param row: start row of the section
        """
        self.title.place(row)

        self.exe.show(row + 1)
        self.pdb.show(row + 2)
        self.out_dir.show(row + 3)
        self.out_name.show(row + 4)

        sep = ttk.Separator(self.master, orient=tk.HORIZONTAL)
        sep.grid(row=row + 5, column=1, columnspan=4, pady=gui.sp_pady, padx=gui.sp_padx, sticky=tk.EW)
        self.skip_h.place(row + 6)
        self.keep_alternative.place(row + 7)

    def mouse_over(self, desc):
        """
        Activates mouse-over behaviour to each widget of the general section.
        :param desc: description box
        """
        self.title.mouse_over(desc)
        self.pdb.mouse_over(desc)
        self.out_dir.mouse_over(desc)
        self.exe.mouse_over(desc)
        self.keep_alternative.mouse_over(desc)
        self.skip_h.mouse_over(desc)
        self.out_name.mouse_over(desc)

    def save(self):
        """
        Updates the user configuration with values from the entry fields.
        """
        if self.pdb.validate():
            cfg.user.pdb_path = adjust_file_path(self.pdb.var.get())
        cfg.user.out_dir = adjust_file_path(self.out_dir.var.get())
        if self.out_name.validate():
            cfg.user.out_name = self.out_name.var.get()
        if self.exe.validate():
            cfg.user.exe = adjust_file_path(self.exe.var.get())
        cfg.user.keep_alternative = bool(self.keep_alternative.var.get())
        cfg.user.skip_h = bool(self.skip_h.var.get())
        cfg.save_config()

    def update(self):
        """
        Updates the entry fields with the current values from the user configuration.
        """
        self.pdb.var.set(adjust_file_path(cfg.user.pdb_path))
        self.out_dir.var.set(adjust_dir_path(cfg.user.out_dir))
        self.out_name.var.set(cfg.user.out_name)
        self.exe.var.set(adjust_file_path(cfg.user.exe))
        self.keep_alternative.var.set(cfg.user.keep_alternative)
        self.skip_h.var.set(cfg.user.skip_h)
        self.master.update()


class GUI:
    def __init__(self, app):
        self.app = app
        # set up the widget styles
        self.style = ttk.Style()
        self.styles()
        # set up the frames
        self.left, self.right = self.setup()
        # set up the sections_verbose
        self.desc = Description(self.right, 0)
        self.desc.configure(height=gui.df_height, width=gui.df_width)
        self.pore_id = PoreIDSection(self.left, self.desc, 10)                  # type: PoreIDSection
        self.general = GeneralSection(self.left, self.desc, 0)                  # type: GeneralSection
        self.axis = AxisSection(self.right, self.desc, 2)                       # type: AxisSection
        self.gate = GateSection(self.right, self.desc, 10)                      # type: GateSection
        # button for running PROPORES
        self.run_button = RButton(self.left, gui.key.run_button, lambda: self.run())
        self.run_button.mouse_over(self.desc)
        self.run_button.place(22)
        # set up the menu bar
        self.menubar()
        # centre the application window
        self.centre()

    def menubar(self):
        """
        Sets up the menu bar with sub-menus for settings, the file parser and help.
        """
        menubar = tk.Menu(self.app)
        # set the menu as the menu bar of the application window
        self.app.configure(menu=menubar)

        # add settings options to the menu bar
        menubar.add_command(label=gui.desc[gui.key.ms_save].label, command=self.save)
        menubar.add_command(label=gui.desc[gui.key.ms_def].label, command=self.default)
        menubar.add_command(label=gui.desc[gui.key.ms_save_def].label, command=self.save_default)

        return menubar

    def setup(self):
        """
        Sets up window borders and section frames.
        :return: left frame, right frame
        """
        # set the title
        self.app.title(gui.title)
        # set up the window border frames and the separating frame
        ttk.Frame(self.app).grid(row=0, column=0, sticky=tk.NSEW, ipadx=gui.wl_padx)
        ttk.Frame(self.app).grid(row=0, column=4, sticky=tk.NSEW, ipadx=gui.wr_padx)
        ttk.Frame(self.app).grid(row=0, column=2, sticky=tk.NSEW, ipadx=gui.col_padx)
        # don't show bottom padding on Mac OS
        if cfg.os != cfg.mac:
            ttk.Frame(self.app).grid(row=1, column=0, ipady=gui.wb_pady)
        self.app.rowconfigure(0, weight=1)
        # set up the left and right section frame
        left = ttk.Frame(self.app)
        right = ttk.Frame(self.app)
        left.grid(row=0, column=1, sticky=tk.NSEW)
        right.grid(row=0, column=3, sticky=tk.NSEW)
        # change the expanding behaviour of the right frame
        right.rowconfigure(1, weight=1)
        left.rowconfigure(8, weight=1)

        return left, right

    def centre(self):
        """
        Places the application window in the middle of the screen.
        """
        self.app.update()
        # compute the current height and width
        h = self.app.winfo_reqheight()
        w = self.app.winfo_reqwidth()
        # compute the x and y coordinates of the window using the screen width and height
        x = (self.app.winfo_screenwidth() - w) // 2
        y = (self.app.winfo_screenheight() - h) // 2

        # set the window position (x coordinate, y coordinate) on the screen
        self.app.geometry('{0}x{1}+{2}+{3}'.format(w, h, x, y))
        self.app.update()

    def fonts(self, style, fnt, fg=None):
        """
        Build the font by trying to convert it into a Tkinter font.
        :param style: base_style name
        :param fnt: tuple of font name, and optionally size and weight
        :param fg: text colour
        """
        # no font specified
        if not fnt:
            return

        if not fnt[0]:
            return
        # see if the font is a TKinter font
        try:
            temp = (font.nametofont(fnt[0]), *fnt[1:])
        except Exception:
            temp = fnt

        if fg:
            self.style.configure(style, foreground=fg, font=temp)
        else:
            self.style.configure(style, font=temp)

    def styles(self):
        """
        Sets up the widget styles.
        """
        if cfg.gui.theme in self.style.theme_names():
            self.style.theme_use(cfg.gui.theme)

        self.fonts(gui.dt_style, gui.dt_font, gui.dt_colour)
        self.fonts(gui.st_style, gui.st_font, gui.st_colour)
        self.fonts(gui.sst_style, gui.sst_font, gui.sst_colour)
        self.fonts(gui.fb_style, gui.fb_font)
        self.fonts(gui.el_style, gui.el_font)
        self.fonts(gui.sl_style, gui.sl_font)
        self.fonts(gui.cb_style, gui.cb_font)
        self.fonts(gui.rb_style, gui.rb_font)
        self.fonts(gui.om_style, gui.om_font)
        self.fonts(gui.tcb_style, gui.tcb_font, gui.st_colour)
        self.fonts(gui.rdb_style, gui.rdb_font)

    def default(self):
        """
        Restore the default settings.
        """
        cfg.restore_default()
        self.general.update()
        self.pore_id.update()
        self.axis.update()
        self.gate.update()

    def save(self):
        """
        Save the current settings.
        """
        self.general.save()
        self.pore_id.save()
        self.axis.save()
        self.gate.save()
        cfg.save_config()

    def save_default(self):
        """
        First restore default settings and then save them.
        """
        self.default()
        self.save()

    @staticmethod
    def open(file):
        """
        Opens the user manual PDF file.
        """
        try:
            # open command depending on the operating system
            if cfg.os == cfg.linux:
                Popen(["xdg-open", file])
            elif cfg.os == cfg.mac:
                Popen(["open", file])
            elif cfg.os == cfg.windows:
                os.startfile(file)
            elif cfg.os == cfg.cygwin:
                os.startfile(file)
        except Exception:
            messagebox.showerror('Open Error',
                                 'An unhandled exception occurred when opening the following file:\n{0}'
                                 '\n\n{1}'.format(file, traceback.format_exc()))

    def validate(self):
        """
        Check if user input of enabled components is valid.
        :return: True if all enabled inputs are valid, False otherwise
        """
        # General: Propores exectuable, input PDB file and output directory
        if not self.general.exe.validate():
            return False

        if not self.general.pdb.validate():
            return False

        if not self.general.out_name.validate():
            return False

        # Pore ID: resolution, solvent and probe radius, volume threshold
        if self.pore_id.title.var.get():
            if not self.pore_id.resolution.validate():
                return False

            if not self.pore_id.solvent.validate():
                return False

            if not self.pore_id.probe.validate():
                return False

            if not self.pore_id.volume.validate():
                return False

        # Axis Trace: surface area threshold and input
        if self.axis.title.var.get():
            if not self.axis.surface.validate():
                return False

            # only check axis input if pore ID is not enabled
            if not self.pore_id.title.var.get():
                # single file input
                if self.axis.selection.get() == 0:
                    if not self.axis.single_input.validate():
                        return False
                # directory input
                elif self.axis.selection.get() == 1:
                    if not self.axis.directory_input.validate():
                        return False

        # Gate Opening: clash tolerance and input
        if self.gate.title.var.get():
            if not self.gate.clash.validate():
                return False

            # only check gate input if pore ID is not enabled
            if not self.pore_id.title.var.get():
                # single file input
                if self.gate.selection.get() == 0:
                    if not self.gate.single_input.validate():
                        return False
                # directory input
                elif self.axis.selection.get() == 1:
                    if not self.gate.directory_input.validate():
                        return False

        return True

    def run(self):
        """
        Run all enabled PROPORES components.
        """
        # at least one component needs to be enabled, show an error message if that is not the case
        if not self.pore_id.title.var.get() and not self.axis.title.var.get() and not self.gate.title.var.get():
            messagebox.showerror('Input Error',
                                 'Please enable at least one of Pore ID, Axis Trace or Gate Opening.')
            return

        # check if the input of all enabled components is valid, show an error message if that is not the case
        if not self.validate():
            messagebox.showerror('Input Error', 'PROPORES could not be run due to erroneous input.')
            return

        # setup the command line arguments for PROPORES, starting with the executable, PDB input and output directory
        args = [adjust_file_path(self.general.exe.var.get(), abs=True),
                '-i', quote(adjust_file_path(self.general.pdb.var.get(), abs=True)),
                '-o', quote(adjust_dir_path(self.general.out_dir.var.get(), abs=True))]

        # set the output name, if specified
        if self.general.out_name.var.get():
            args += ['--name', self.general.out_name.var.get()]

        # set the flag for skipping H-atoms during PDB parsing, if specified
        if self.general.skip_h.var.get():
            args += ['--skip-H-atoms']

        # set the flag for keeping during PDB parsing, if specified
        if self.general.keep_alternative.var.get():
            args += ['--keep-alternative']

        # check if Pore ID is enabled for running
        if self.pore_id.title.var.get():
            # add Pore ID flag, resolution, solvent and probe radius, volume threshold
            args += ['pore-id',
                     '-b', quote(self.pore_id.resolution.var.get()),
                     '-s', quote(self.pore_id.solvent.var.get()),
                     '-p', quote(self.pore_id.probe.var.get()),
                     '-v', quote(self.pore_id.volume.var.get())]

            # check preparation flags for axis trace and gate opening
            if self.pore_id.gate_prep.var.get():
                args += ['--axis-preparation']
            if self.pore_id.axis_prep.var.get():
                args += ['--gate-preparation']

            # set mutually exclusive computation mode flag, if auto-detect is disabled
            if self.pore_id.computation_mode.var.get() == cfg.options.cylinder_ray_trace:
                args += ['--cylinder-ray-trace']
            elif self.pore_id.computation_mode.var.get() == cfg.options.cylinder_standalone:
                args += ['--cylinder-standalone']

            # set pore type flag, if only a specific type is selected
            if self.pore_id.pore_type.var.get() == cfg.options.pore_only:
                args += ['--only-pores']
            elif self.pore_id.pore_type.var.get() == cfg.options.cavity_only:
                args += ['--only-cavities']

        # check if axis trace is enabled
        if self.axis.title.var.get():
            # add axis trace flag and surface area threshold
            args += ['axis-trace',
                     '-spt', quote(self.axis.surface.var.get())]

            # if pore ID is not enabled, check for single file or directory input
            if not self.pore_id.title.var.get():
                # single file input
                if self.axis.selection.get() == 0:
                    args += ['-ts', quote(adjust_file_path(self.axis.single_input.var.get(), abs=True))]
                # directory input
                elif self.axis.selection.get() == 1:
                    args += ['-td', quote(adjust_dir_path(self.axis.directory_input.var.get(), abs=True))]

        # check if gate open is enabled
        if self.gate.title.var.get():
            # add gate open flag and clash tolerance
            args += ['gate-open',
                     '-ct', quote(self.gate.clash.var.get())]

            # if pore ID is not enabled, check for single file or directory input
            if not self.pore_id.title.var.get():
                # single file input
                if self.gate.selection.get() == 0:
                    args += ['-gs', quote(adjust_file_path(self.gate.single_input.var.get(), abs=True))]
                # directory input
                elif self.gate.selection.get() == 1:
                    args += ['-gd', quote(adjust_dir_path(self.gate.directory_input.var.get(), abs=True))]

            # check if difficulty restrictions are applied and if yes, if difficulty re-estimation is desired
            if self.gate.difficulty.var.get() == cfg.options.difficulty_medium:
                args += ['--skip-medium-gates']
                if self.gate.re_estimate.var.get():
                    args += ['--re-estimate']
            elif self.gate.difficulty.var.get() == cfg.options.difficulty_hard:
                args += ['--skip-hard-gates']
                if self.gate.re_estimate.var.get():
                    args += ['--re-estimate']

        # run PROPORES with the command line options and report errors
        try:
            run(args=' '.join(args), shell=True)
            self.open(adjust_dir_path(self.general.out_dir.var.get(), abs=True))
        except Exception:
            messagebox.showerror('PROPORES Error',
                                 'An unhandled exception occurred while trying to run PROPORES.\n\n'
                                 '{0}'.format(traceback.format_exc()))
