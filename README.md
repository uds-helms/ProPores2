# ProPores2

## Introduction
ProPores2 is a C++ command line tool for analysing pores (and cavities) in proteins. It provides three components that can be run together or individually:

**Pore Identification:** Identifies all pores in a given protein PDB file, and for each pore generates a tab-separated file with lining residues and a pseudo PDB file with the empty space within the pore. The pseudo PDB file can be used for visualising the pore together with the protein. 

**Axis Trace:** Determines the axis of a pore and writes it into a pseudo PDB file for visualisation. Some pores have multiple axes, in which cases several output files are generated.

**Gate Opening:** Rotates the shared lining residues of two neighbouring pores in an effort to open the gate between them as much as possible. The result is a PDB file of the entire protein with the rotated residues.

## Step 1: Installation
The installation requires a C++17 compiler and CMake 3.13 or newer.

### Windows 7 and 10
1. Install a C++17 compiler, or update an existing C++ compiler if necessary. [Microsoft provides a tutorial for installing the Visual C++ compiler.](https://docs.microsoft.com/en-us/cpp/build/vscpp-step-0-installation?view=vs-2019) Make sure to tick "Desktop development in C++" during the isntallation of Visual Studio.
2. Install CMake 3.13 or newer, or update an existing CMake if necessary. [CMake.org provides a guide for installing CMake on your computer.](https://cmake.org/install) This should add `cmake` to the environment variables but might require a restart of the computer to take effect.
3. Download this repository as a ZIP file and extract it to a location of your choice, or use Git to clone the repository.
4. Open the ProPores2 folder where you have extracted or cloned it, right-click on `install_windows.bat` and select `Run as administrator`. This creates the executable `propores.exe` in the ProPores2 folder.

   If `Run as administrator` does not work, try simply double-clicking `install_windows.bat`. 

   If that does not work either, open the Windows command line, e.g. via the search field in the task bar, and navigate to the PROPORES folder with the `cd` command ([tutorial](https://docs.microsoft.com/en-us/windows-server/administration/windows-commands/cd)). Then execute the following commands:
   ```
   if exist "build\" rmdir /Q /S build
   cmake -S . -B build
   cmake --build build --config Release --target install
   if exist "build\" rmdir /Q /S build
   ```

### MacOS and Linux
1. Install a C++17 compiler, or update an existing C++ compiler if necessary. If you are using the GCC compiler, you need GCC-9 (and G++9) or newer. 

   Some older Ubuntu long-term support versions, e.g. Ubuntu 18.04 LTS, do not provide GCC-9 in their package manager. In that case, you can install GCC-9 as follows:
   ```
   sudo apt update
   sudo apt install build-essential
   sudo apt-get install manpages-dev
   sudo apt install software-properties-common
   sudo add-apt-repository ppa:ubuntu-toolchain-r/test
   sudo apt install gcc-9 g++-9
   sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 90 --slave /usr/bin/g++ g++ /usr/bin/g++-9 --slave /usr/bin/gcov gcov /usr/bin/gcov-9
   sudo update-alternatives --config gcc
   ```
2. Install CMake 3.13 or newer, or update an existing CMake if necessary. [CMake.org provides a guide for installing CMake on your computer.](https://cmake.org/install)
3. Download this repository as a ZIP file and extract it to a location of your choice, or use Git to clone the repository.
4. Open the ProPores2 folder in the terminal and run 
   ```
   sh install_linux_macOS.sh
   ```
   This creates the executable `propores` in the ProPores2 folder.
5. (Optional) If you want to use the user interface, install Python3 with the packages `tkinter` and `pyyaml`. Open the ProPores2 folder in the terminal and run
   ```
   sh install_python_ubuntu_debian.sh
   ```
   or
   ```
   install_python_macOS.sh
   ```


## Step 2: Running ProPores2 on the Command Line
Open the terminal (or command prompt) and navigate to the ProPores2 folder with the `propores` executable. The general command setup on Windows is
```
propores.exe <run command(s)> -i <path> -o <path> [options...]
```
and on MacOS or Linux:
```
./propores <run command(s)> -i <path> -o <path> [options...]
```
The three run commands are `pore-id`, `axis-trace` and `gate-open`, and at least one of them has to be given. The input protein PDB file `-i` and the output directory `-o` are always required. To see a full list of available options, run `propores.exe -h` or `propores.exe --help`. Some use case examples are listed below.

### Example 1: Only Pore Identification
Running only pore identification with 0.5 Angstrom resolution, with the option to run axis trace and gate opening at a later point.
```
propores.exe pore-id -i input/1EA5_R.pdb -o output --name example -b 0.5 --axis-preparation --gate-preparation
```
This writes the pore identification output, as well as files for independent axis trace and gate opening, to the folder `output/example/` in the ProPores2 folder.

### Example 2: Only Axis Trace
Running only axis trace based on previously generated files. For a single pore:
```
propores.exe axis-trace -i input/1EA5_R.pdb -o output --name example -ts output/example/axis_trace_input/0_pore.tsv
```
and for all pores in a directory:
```
propores.exe axis-trace -i input/1EA5_R.pdb -o output --name example -td output/example/axis_trace_input
```
This writes the axes to the folder `output/example/axes`.

### Example 3: Only Gate Opening
Running only opening based on previously generated files. For a single pore:
```
propores.exe gate-open -i input/1EA5_R.pdb -o output --name example -gs output/example/gate_open_input/0_pore.tsv
```
and for all pores in a directory:
```
propores.exe gate-open -i input/1EA5_R.pdb -o output --name example -td output/example/gate_open_input
```
While most gates can be opened within seconds, some gates have a very large number of possible residue configurations and can therefore take a lot longer. ProPores2 estimates the difficulty of each gate based on the residue composition, and gates can be skipped based on that difficulty:
```
propores.exe gate-open -i input/1EA5_R.pdb -o output --name example -td output/example/gate_open_input --skip-hard-gates
```
This writes the protein versions with open gates to the folder `output/example/gate_open`.

### Example 4: Everything at Once
```
propores.exe pore-id axis-trace gate-open -i input/1EA5_R.pdb -o output --name example -b 0.5 --skip-hard-gates
```

## Step 2: Running ProPores2 via the Graphical User Interface
ProPores2 comes with a graphical user interface that can be used instead of the command line. After installing ProPores2 as described in Step 1, the user interface can be opened as follows:

### Windows 7 and 10
Double-click on `ProporesGUI_32.exe` or `ProporesGUI_64.exe` depending on whether you have a 32-bit or 64-bit Windows version. If in doubt, select the 32-bit version.

### MacOS and Linux
Open the terminal and navigate to the ProPores2 folder with the `propores` executable. Run
```
python3 propores_gui.py
```

### Avoiding Errors
Make sure that `ProporesGUI_32.exe`, `ProporesGUI_64.exe` or `propores_gui.py` are in the same folder as the `config.yaml` file and the `src` folder, otherwise the user interface will not open. Check if the automatically detected path to the ProPores2 executable is correct, otherwise set it manually under `General → ProPores executable`.

## Attribution
The approach was designed by Po-Hsien Lee and published as:

>Lee, PH, Helms, V (2012). Identifying continuous pores in protein structures with PROPORES by computational repositioning of gating residues. Proteins, 80, 2:421-32. https://www.ncbi.nlm.nih.gov/pubmed/22095919

The original version was implemented by Po-Hsien Lee (2011) in Perl. ProPores2 is a faster and more memory-efficient C++ implementation by Markus Hollander and Moomal Aziz.
