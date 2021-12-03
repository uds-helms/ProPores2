# ProPores2

## Introduction
ProPores2 is a C++ command line tool for analysing pores (and cavities) in proteins. It provides three components that can be run together or individually:

**Pore Identification:** Identifies all pores in a given protein PDB file, and for each pore generates a tab-separated file with lining residues and a pseudo PDB file with the empty space within the pore. The pseudo PDB file can be used for visualising the pore together with the protein.

**Axis Trace:** Determines the axis of a pore and writes it into a pseudo PDB file for visualisation. Some pores have multiple axes, in which cases several output files are generated.

**Gate Opening:** Rotates the shared lining residues of two neighbouring pores in an effort to open the gate between them as much as possible. The result is a PDB file of the entire protein with the rotated residues.

## Webservice
ProPores2 is also available as a webservice: https://service.bioinformatik.uni-saarland.de/propores/index.html

## Citation

>Hollander, M., Rasp, D., Aziz, M., Helms, V. (2021). ProPores2: Web Service and Stand-Alone Tool for Identifying, Manipulating, and Visualizing Pores in Protein Structures. Journal of Chemical Information and Modeling 2021 61 (4), 1555-1559. DOI: 10.1021/acs.jcim.1c00154. https://pubs.acs.org/doi/10.1021/acs.jcim.1c00154

## Previous Publication
The approach was designed by Po-Hsien Lee and published as:

>Lee, P.H., Helms, V. (2012). Identifying continuous pores in protein structures with PROPORES by computational repositioning of gating residues. Proteins, 80, 2:421-32. https://www.ncbi.nlm.nih.gov/pubmed/22095919

The original version was implemented by Po-Hsien Lee (2011) in Perl. ProPores2 is a faster and more memory-efficient C++ implementation by Markus Hollander and Moomal Aziz.
