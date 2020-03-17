# PPROPORS 2.0

## Introduction
PROPORES 2.0 is a C++ command line tool for analysing pores (and cavities) in proteins. It provides three components that can be run together or individually:

**Pore Identification:** Identifies all pores in a given protein PDB file, and for each pore generates a tab-separated file with lining residues and a pseudo PDB file with the empty space within the pore. The pseudo PDB file can be used for visualising the pore together with the protein. 

**Axis Trace:** Determines the axis of a pore and writes it into a pseudo PDB file for visualisation. Some pores have multiple axes, in which cases several output files are generated.

**Gate Opening:** Rotates the shared lining residues of two neighbouring pores in an effort to open the gate between them as much as possible. The result is a PDB file of the entire protein with the rotated residues.

## Step 1: Installation
The installation requires a C++17 compiler and CMake 3.15 or newer.

### Windows 7 and 10
1. Install a C++17 compiler, or update an existing C++ compiler if necessary. [Microsoft provides a tutorial for installing the Visual C++ compiler.](https://docs.microsoft.com/en-us/cpp/build/vscpp-step-0-installation?view=vs-2019)
2. Install CMake 3.15 or newer, or update an existing CMake if necessary. [CMake.org provides a guide for installing CMake on your computer.](https://cmake.org/install) This should add `cmake` to the environment variables but might require a restart of the computer to take effect.
3. Download this repository as a ZIP file and extract it to a location of your choice, or use Git to clone the repository.
4. Open the PROPORES folder where you have extracted or cloned it, right-click on `install_windows.bat` and select `Run as administrator`. This creates the executable `propores.exe` in the PROPORES folder.

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
2. Install CMake 3.15 or newer, or update an existing CMake if necessary. [CMake.org provides a guide for installing CMake on your computer.](https://cmake.org/install)
3. Download this repository as a ZIP file and extract it to a location of your choice, or use Git to clone the repository.
4. Open the PROPORES folder in the terminal and run 
   ```
   sh install_linux_macOS.sh
   ```
   This creates the executable `propores` in the PROPORES folder.


## Step 2: Running PROPORES
Open the terminal (or command prompt) and navigate to the PROPORES folder with the `propores` executable. The general command setup on Windows is
```
propores.exe <run command(s)> -i <path> -o <path> [options...]
```
and on MacOS and Linux:
```
./propores <run command(s)> -i <path> -o <path> [options...]
```
The three run commands are `pore-id`, `axis-trace` and `gate-open`, and at least one of them has to be given. The input protein PDB file `-i` and the output directory `-o` are always required. To see a full list of available options, run `propores.exe -h` or `propores.exe --help`. Some use case examples are listed below.

### Example 1: Only Pore Identification
Running only pore identification with 0.5 Angstrom resolution, with the option to run axis trace and gate opening at a later point.
```
propores.exe pore-id -i input/1EA5_R.pdb -o output --name example -b 0.5 --axis-preparation --gate-preparation
```
This writes the pore identification output, as well as files for independent axis trace and gate opening, to the folder `output/example/` in the PROPORES folder.

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
While most gates can be opened within seconds, some gates have a very large number of possible residue configurations and can therefore take a lot longer. PROPORES estimates the difficulty of each gate based on the residue composition, and gates can be skipped based on that difficulty:
```
propores.exe gate-open -i input/1EA5_R.pdb -o output --name example -td output/example/gate_open_input --skip-hard-gates
```
This writes the protein versions with open gates to the folder `output/example/gate_open`.

### Example 4: Everything at Once
```
propores.exe pore-id axis-trace gate-open -i input/1EA5_R.pdb -o output --name example -b 0.5 --skip-hard-gates
```

## Attribution
The approach was designed by Po-Hsien Lee and published as:

>Lee, PH, Helms, V (2012). Identifying continuous pores in protein structures with PROPORES by computational repositioning of gating residues. Proteins, 80, 2:421-32. https://www.ncbi.nlm.nih.gov/pubmed/22095919

The original version was implemented by Po-Hsien Lee (2011) in Perl. PROPORES 2.0 is a faster and more memory-efficient C++ implementation by Markus Hollander and Moomal Aziz.
