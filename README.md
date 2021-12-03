# ProPores2

## Introduction
ProPores2 is a C++ command line tool for analysing pores (and cavities) in proteins. It provides three components that can be run together or individually:

**Pore Identification:** Identifies all pores in a given protein PDB file, and for each pore generates a tab-separated file with lining residues and a pseudo PDB file with the empty space within the pore. The pseudo PDB file can be used for visualising the pore together with the protein. 

**Axis Trace:** Determines the axis of a pore and writes it into a pseudo PDB file for visualisation. Some pores have multiple axes, in which cases several output files are generated.

**Gate Opening:** Rotates the shared lining residues of two neighbouring pores in an effort to open the gate between them as much as possible. The result is a PDB file of the entire protein with the rotated residues.
