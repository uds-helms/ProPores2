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

#ifndef PROPORES_TRACE_H
#define PROPORES_TRACE_H

#include "grid.h"

// extract a single pore or cavity grid box cluster from a single input file
void cluster_from_file(const std::string &file_path, std::vector<PoreCluster> &clusters);

// extract pore or cavity grid box clusters from a directory of input files
void clusters_from_directory(const std::string &dir_path, std::vector<PoreCluster> &clusters);

// create an output file for each pore or cavity grid box cluster
void output_trace_cluster(const Settings &settings, const ProteinGrid &grid);

// computes starting points of pore/cavity axes
void compute_starting_points(PoreGrid &grid, size_t min_size);

// computes the not yet visited pore/cavity box with the lowest score
int min_index(const std::vector<double> &scores, const std::vector<uint8_t> &visited);

// use Dijkstra to compute the pore/cavity axes from the given starting point
void trace(PoreGrid &grid, size_t start);

// generates an output file for each axis of a pore/cavity
void write_axes(const PoreGrid &grid, const std::string &path_base);

void axis_trace(const Settings &settings, const std::vector<PoreCluster> &clusters);

#endif //PROPORES_TRACE_H
