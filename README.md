# PPROPORS 2.0

## Introduction


## Step 1: Installation
The installation requires a C++17 compiler and CMake 3.15 or newer.

### Windows 7 and 10
1. Install a C++17 compiler, or update an existing C++ compiler if necessary. [Microsoft provides a tutorial for installing the Visual C++ compiler.](https://docs.microsoft.com/en-us/cpp/build/vscpp-step-0-installation?view=vs-2019)
2. Install CMake 3.15 or newer, or update an existing CMake if necessary. [CMake.org provides a guide for installing CMake on your computer.](https://cmake.org/install)
3. Download this repository as a ZIP file and extract it to a location of your choice, or use Git to clone the repository.
4. Open the PROPORES folder where you have extracted or cloned it, right-click on `install_windows.bat` and select `Run as administrator`. This creates the executable `propores.exe` in the PROPORES folder.

### MacOS and Linux
1. Install a C++17 compiler, or update an existing C++ compiler if necessary. If you are using the GCC compiler, you need GCC-9 (and G++9) or newer. Some older Ubuntu long-term support versions, e.g. Ubuntu 18.04 LTS, do not provide GCC-9 in their package manager. In that case, you can install GCC-9 as follows:
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
   sh install_macos_linux.sh
   ```
   This creates the executable `propores` in the PROPORES folder.


## Step 2: Running PROPORES

## Attribution

The approach was designed by Po-Hsien Lee and published as:

>Lee, PH, Helms, V (2012). Identifying continuous pores in protein structures with PROPORES by computational repositioning of gating residues. Proteins, 80, 2:421-32. https://www.ncbi.nlm.nih.gov/pubmed/22095919

The original version was implemented by Po-Hsien Lee (2011) in Perl. PROPORES 2.0 is a faster and more memory-efficient C++ implementation by Markus Hollander and Moomal Aziz.

## FAQ
