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

#ifndef PROPORES_GATE_OPEN_H
#define PROPORES_GATE_OPEN_H

#include <stack>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include "grid.h"
#include "gate.h"
#include "settings.h"

// estimate the difficulty of a gate by considering the potential number of rotamer combinations
size_t estimate_difficulty(const Gate &gate);

// extract potential gates between neighbouring pores/cavities from the protein grid
void gates_from_grid(ProteinGrid &grid, std::vector<Gate> &gates);

// extract a single gate from a single input file
void gate_from_file(const std::string &file_path, std::vector<Gate> &gates);

// extract gates from a directory of input files
void gates_from_directory(const std::string &dir_path, std::vector<Gate> &gates);

// create an output file for each gate with the pore/cavity IDs, estimated difficulty and gating residues
void output_gates(const Settings &settings, const std::vector<Gate> &gates);

// generate the library with valid rotamer angles for each residue
void load_rotamer_library(const Settings &cfg, std::unordered_map<ResidueType, std::vector<std::vector<int>>> &lib);

// compute the dihedral angles of a given gate residue
void compute_dihedral_angles(const Gate &gate, Residue &res);

// rotate the side chain of a given residue to generate the specified rotamer
void side_chain_rotate(const Residue &res, Rotamer &rotamer);

// check if a residue side chain rotamer clashes substantially with the rest of the protein
bool serious_clash(const Gate &gate, const Residue &res, Rotamer &rotamer, double tolerance, double gamma_tolerance);

// generate valid rotamers for a given gate residue
void generate_rotamers(const Settings &settings, const Gate &gate, Residue &res,
                       const std::vector<std::vector<int>> &rotamer_angles);

// check how many residues do not have valid rotamers, and add to them the original residue configuration
size_t fix_residues(Gate &gate);

// sort residues by the number of rotamers in ascending order and adjust their indices accordingly
void sort_residues(Gate &gate);

// compute the interaction energy of two rotamers based on the distance between them
double interaction_energy(const Rotamer &rot_1, const Rotamer &rot_2);

// eliminates rotamers that are dead-ends in the effort to minimise the gate energy
void dead_end_elimination(Gate &gate);

// remove rotamers that have been eliminated
void remove_eliminated(Gate &gate);

// check all combinations of the remaining rotamers to determine the combination that yields the lowest gate energy
std::vector<Rotamer> exhaustive_search(const Gate &gate);

// main gate-open routine
void gate_open(Settings &settings, std::vector<Gate> &gates);

#endif //PROPORES_GATE_OPEN_H
