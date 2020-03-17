import os
import warnings
import tkinter as tk
import src.gui as gui

__author__ = 'Markus Hollander'
__copyright__ = 'Copyright (C) 2020 Markus Hollander'
__license__ = 'GPL V3'
__version__ = '2.0'

"""
Main function that creates the graphical user interface (GUI), which controls the program flow.
"""
if __name__ == '__main__':
    # don't show warnings in std-out
    warnings.filterwarnings('ignore')

    # extract and set the proper working directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # build and run the GUI
    root = tk.Tk()
    app_gui = gui.GUI(root)
    root.mainloop()
    # extract and set the proper working directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # build and run the GUI
    root = tk.Tk()
    app_gui = gui.GUI(root)
    root.mainloop()
