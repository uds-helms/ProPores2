/*
 * Copyright (C) 2020 Markus Hollander (markus.hollander@bioinformatik.uni-saarland.de). All rights reserved.
 *
 * This file is part of PROPORES 2.0. PROPORES was originally developed by Po-Hsien Lee (2011) in Perl under the terms
 * of the GNU General Public License version 3 or higher. The approach was published as:
 * Lee, PH, Helms, V (2012). Identifying continuous pores in protein structures with PROPORES by computational
 * repositioning of gating residues. Proteins, 80, 2:421-32. https://www.ncbi.nlm.nih.gov/pubmed/22095919
 *
 * PROPORES 2.0 is a C++ implementation of the approach outlined in Lee (2012) by Markus Hollander and Moomal Aziz.
 * It is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the license, or (at your option) any later version.
 *
 * PROPORES 2.0 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with this program. If not, see
 * <https://www.gnu.org./licenses/>.
*/

#ifndef PROPORES_READER_H
#define PROPORES_READER_H

#include <vector>
#include "atom.h"

// extract all valid atom entries from a given PDB input file, allow to skip H-atoms and to load alternative locations
// of atoms in addition to the primary location
void parse_PDB(const std::string &filename, std::vector<std::shared_ptr<Atom>> &atoms, bool skip_H,
               bool keep_alternative);

// format an index and coordinate in pseudo PDB format, can be used to output pore/cavity boxes for visualisation
std::string pseudo_pdb(size_t id, const Vec<double> &coord);

// format information of a single atom in PDB format
std::string PDB_atom_entry(const std::shared_ptr<Atom> &atom);

// output the information of all atoms in a PDB file
void write_PDB(const std::string &file_path, const std::vector<std::shared_ptr<Atom>> &atoms);

#endif //PROPORES_READER_H