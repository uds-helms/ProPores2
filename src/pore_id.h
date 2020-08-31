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

#ifndef PROPORES_PORE_ID_H
#define PROPORES_PORE_ID_H

#include <vector>
#include "atom.h"
#include "grid.h"
#include "gate.h"
#include "vector.h"

// scan for empty grid boxes involved in PSP events on the given axis
std::vector<uint8_t> psp_scan(ProteinGrid &grid);

// check if the connecting vector between atom pairs is interrupted by occupied grid boxes
void collision_detection(ProteinGrid &grid);

// determine which cylinder function should be run, if true => run cylinder standalone
bool switch_to_standalone(const ProteinGrid &grid, const Settings &settings);

// form a cylinder around the connecting vector between atom pairs
void cylinder_trace(ProteinGrid &grid, const std::shared_ptr<AtomPair> &atom_pair, double r_sq, int r_grid,
                    double threshold, double search);

// form a cylinder around the connecting vector between atom pairs
void cylinder_completion(ProteinGrid &grid);

// compute which atom pair cylinders are perpendicular to each other
void perpendicularity_cylinder(ProteinGrid &grid);

// compute which atom pairs are perpendicular to each other without computing PSP events or cylinders
void perpendicularity_standalone(ProteinGrid &grid, size_t cores);

// mark grid boxes as part of the background if they are accessible from the solvent without passing through occupied
// boxes or boxes that are (potentially) part of a pore
void seal_gaps(ProteinGrid &grid);

// remove shallow pore regions on the protein surface
void remove_shallow_regions(ProteinGrid &grid, double shallow_radius);

// determine the largest ball radius that can be centred in each empty grid box
void identify_pore_nuclei(ProteinGrid &grid);

// build initial clusters of connected potential pore boxes
void initial_pore_assembly(ProteinGrid &grid, double tolerance);

// after the clusters are established, add the surrounding potential pore/cavity boxes to the nearest cluster
void complete_pores(ProteinGrid &grid);

// for each empty grid box involved in a pore/cavity, compute if it's at the protein surface, in touch with
// residues, in touch with another cluster or in the cluster centre
void label_pore_boxes(ProteinGrid &grid);

// classify clusters as pores or as cavities
void classify_clusters(ProteinGrid &grid);

// for each pore/cavity cluster add the surrounding residues
void lining_residues(ProteinGrid &grid);

// create files for pore/cavity boxes and lining residues
void output_pores(ProteinGrid &grid, const Settings &settings);

// main pore/cavity identification routine
void pore_ID(ProteinGrid &grid, const Settings &settings);

#endif //PROPORES_PORE_ID_H
