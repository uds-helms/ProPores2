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
void trace(PoreGrid &grid, const std::string &output_path_base, size_t start);

void axis_trace(const Settings &settings, const std::vector<PoreCluster> &clusters);

#endif //PROPORES_TRACE_H
