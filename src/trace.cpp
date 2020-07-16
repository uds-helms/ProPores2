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

#include <stack>
#include <cmath>
#include <algorithm>
#include "trace.h"
#include "grid.h"
#include "vector.h"
#include "reader.h"

namespace fs = std::filesystem;

// extract a single pore or cavity grid box cluster from a single input file
void cluster_from_file(const std::string &file_path, std::vector<PoreCluster> &clusters) {
    std::ifstream file(file_path);
    std::string line;
    // get the first line, which contains the ID of the cluster, and whether it is a pore or a cavity
    std::getline(file, line);
    std::vector<std::string> cols = split(r_strip(line));
    // the ID line is incomplete
    if (cols.size() < 3) {
        file.close();
        return;
    }
    // ID_pore or ID_cavity
    cols = split(cols[2], '_');
    size_t id = stoi(cols[0]);
    bool pore = cols[1] == "pore";
    // extract the grid box length that was used to generate the pore/cavity
    std::getline(file, line);
    cols = split(r_strip(line));
    // the grid box length line is incomplete
    if (cols.size() < 7) {
        file.close();
        return;
    }
    double box_length = stod(cols[4]);
    // initialise the grid box cluster and set its type
    PoreCluster cluster = PoreCluster(id, box_length);
    cluster.pore = pore;
    // extract the grid box information
    std::vector<PoreBox> boxes;
    while (std::getline(file, line)) {
        cols = split(r_strip(line), '\t');
        if (cols.size() != 6) continue;
        PoreBox box = PoreBox(Vec<int>(stoi(cols[0]), stoi(cols[1]), stoi(cols[2])));
        box.state = to_box_state(cols[3]);
        box.distance_to_closest = stod(cols[4]);
        box.distance_to_centre = stod(cols[5]);
        boxes.push_back(box);
    }
    file.close();
    // make sure that there is at least one grid box in the input file
    cluster.boxes = std::move(boxes);
    if (!cluster.empty()) clusters.push_back(cluster);
}

// extract pore or cavity grid box clusters from a directory of input files
void clusters_from_directory(const std::string &dir_path, std::vector<PoreCluster> &clusters) {
    for (auto &p: fs::directory_iterator(dir_path)) {
        if (fs::is_directory(p)) continue;
        cluster_from_file(p.path().string(), clusters);
    }
}

// create an output file for each pore or cavity grid box cluster
void output_trace_cluster(const Settings &settings, const ProteinGrid &grid) {
    for (const PoreCluster &cluster: grid.clusters) {
        // generate the output path and create the output file
        fs::path output_path = settings.axis_preparation_dir / fs::path(cluster.name() + ".tsv");
        std::ofstream file(output_path.string());
        // write a header with the ID and type (pore or cavity) of the cluster, as well as the used grid box length
        file << "# Name: " << cluster.name() << std::endl;
        file << "# Grid box length: " << cluster.box_length << " in Angstrom" << std::endl;
        // write the grix box information
        for (const PoreBox &box: cluster.boxes) {
            file << box.coord.x << "\t" << box.coord.y << "\t" << box.coord.z << "\t" << to_str(box.state) << "\t";
            file << box.distance_to_closest << "\t" << box.distance_to_centre << std::endl;
        }
    }
}

// computes starting points of pore/cavity axes
void compute_starting_points(PoreGrid &grid, const size_t min_size) {
    std::vector<uint8_t> in_patch(grid.size());
    Vec<int> max_dist_vec = grid.min;
    double max_dist = 0.0;
    size_t count = 1;

    for (int x = grid.min.x; x < grid.max.x; x++) {
        for (int y = grid.min.y; y < grid.max.y; y++) {
            for (int z = grid.min.z; z < grid.max.z; z++) {
                // 1D index of the box
                size_t index = grid.index(x, y, z);
                // if the box belongs to the pore/cavity, check if it's further away from the protein centre of mass
                // than the current maximum. if so, set it as the new maximum
                if (grid.at(index) != UNCLASSIFIED && grid.distance_to_centre[index] > max_dist) {
                    max_dist = grid.distance_to_centre[index];
                    max_dist_vec = Vec<int>(x, y, z);
                }
                // skip boxes that are not exposed to the solvent or are already in a surface patch
                if (grid.at(index) != TOUCHES_BACKGROUND || in_patch[index]) continue;

                // use a stack to create a surface patch connected to the current box
                std::stack<Vec<int>> to_check;
                to_check.push(Vec<int>(x, y, z));
                std::vector<Vec<int>> patch;
                // iteratively add all connected boxes that are exposed to the solvent to the surface patch
                while (!to_check.empty()) {
                    Vec<int> box = to_check.top();
                    to_check.pop();
                    // add the current box to the current surface patch
                    patch.push_back(box);
                    in_patch[grid.index(box)] = 1;
                    // search all surrounding boxes
                    for (const RelativeSearchBox &rel_box: grid.search_spaces[1]) {
                        Vec<int> search_box = box + rel_box.grid;
                        // skip search boxes that are not in the grid
                        if (!grid.valid(search_box)) continue;
                        size_t idx = grid.index(search_box);
                        // if the search box touches the solvent and is not yet in a surface patch, add it to the stack
                        if (grid.at(idx) == TOUCHES_BACKGROUND && !in_patch[idx]) to_check.push(search_box);
                    }
                }
                // skip surface patches that do not pass the minimum size threshold
                if (patch.size() * grid.box_length_sq < min_size) continue;

                // compute the box in the surface patch with the largest possible radius
                Vec<int> max_radius_vec = patch[0];
                double max_radius = grid.radii[grid.index(max_radius_vec)];
                for (const Vec<int> &vec: patch) {
                    double radius = grid.radii[grid.index(vec)];
                    if (radius > max_radius) {
                        max_radius = radius;
                        max_radius_vec = vec;
                    }
                }
                // add the box as an axes starting point
                AxisPoint point = AxisPoint(grid.index(max_radius_vec), max_radius_vec,
                        "surface patch " + std::to_string(count));
                grid.starting_points.push_back(point);
                count++;
            }
        }
    }
    // if there were no sufficiently large surface patches, use the box with the largest distance from the protein
    // centre of mass as the only axes starting point
    if (grid.starting_points.empty()) {
        AxisPoint point = AxisPoint(grid.index(max_dist_vec), max_dist_vec, "furthest away from centre of mass");
        grid.starting_points.push_back(point);
    }
}

// computes the not yet visited pore/cavity box with the lowest score
int min_index(const std::vector<double> &scores, const std::vector<uint8_t> &visited) {
    // initialise the minimum score and index
    int min_idx = -1;
    double min_score = INFINITY;

    for (size_t i = 0; i < scores.size(); i++) {
        // skip boxes that have already been visited or have a score of infinity (= not (yet) reachable)
        if (visited[i] || std::isinf(scores[i])) continue;
        // if the current box has a lower score, set it as the new minimum
        if (scores[i] < min_score) {
            min_idx = i;
            min_score = scores[i];
        }
    }

    return min_idx;
}

// use Dijkstra to compute the pore/cavity axes from the given starting point
void trace(PoreGrid &grid, const size_t start_idx) {
    // DIJKSTRA DATA
    // initialise the scores to infinity
    std::vector<double> scores(grid.size(), INFINITY);
    // mark all boxes as unvisited
    std::vector<uint8_t> visited(grid.size());
    // initialise the trace back to the 1D box index outside the grid
    std::vector<size_t> trace_back(grid.size(), grid.size());
    // initialise the path lengths to infinity
    std::vector<double> lengths(grid.size(), INFINITY);

    // STARTING POINT INITIALISATION
    AxisPoint start = grid.starting_points[start_idx];
    // the starting point has the smallest score of 0.0
    scores[start.index] = 0.0;
    // the starting point has no previous box in the trace back
    trace_back[start.index] = start.index;
    lengths[start.index] = 0;

    // mark grid boxes that are not part of the pore/cavity als already visited, so that their scores stay at infinity
    // and they will not get considered as part of the axis
    for (size_t i = 0; i < visited.size(); i++) {
        if (grid.at(i) == UNCLASSIFIED) visited[i] = 1;
    }

    // start with the starting point (int and not size_t because min_index returns -1 if there is no next box)
    int current_index = start.index;
    // perform Dijkstra's algorithm as long as there are still not yet visited boxes left
    while (current_index >= 0) {
        // get the 3D coordinates, score and path length of the current box
        Vec<int> current_box = grid.coordinates[current_index];
        double current_score = scores[current_index];
        size_t current_path_length = lengths[current_index];
        // process all direct neighbours of the current box (including diagonal)
        for (const RelativeSearchBox &rel_box: grid.search_spaces[1]) {
            // skip neighbouring boxes that are not in the grid
            if (!grid.valid(current_box + rel_box.grid)) continue;
            // get the 1D index of the neighbour box
            size_t neighbour = grid.index(current_box + rel_box.grid);
            // skip already visited neighbours
            if (visited[neighbour]) continue;
            // check if the neighbour box has a lower score if reached from the current box
            double new_neighbour_score = current_score + grid.individual_score(neighbour);
            // if yes, adjust the score and set the current box as the new predecessor in the trace back
            if (new_neighbour_score < scores[neighbour]) {
                scores[neighbour] = new_neighbour_score;
                trace_back[neighbour] = current_index;
                lengths[neighbour] = current_path_length + 1;
            }
        }
        // mark the current box as visited
        visited[current_index] = 1;
        // get the index of the next not yet visited box (= -1 if no such box exists)
        current_index = min_index(scores, visited);
    }

    // find the box on the outside of the pore/cavity with the lowest score as the end point of the axis
    // initialise the minimum score as infinity
    double min_score = INFINITY;
    int min_index = -1;
    // iterate over all boxes that are part of the pore/cavity and lay on its outer layer. make sure they are not
    // the start point and do not have a score of infinity
    for (size_t i = 0; i < grid.size(); i++) {
        if (grid.at(i) == UNCLASSIFIED || grid.at(i) == ON_CLUSTER_INSIDE || i == start.index
            || std::isinf(scores[i])) {
            continue;
        }
        // adjusted the box score by dividing it by the length of the path to the power of 1.5 to favour longer
        // paths, which avoids selecting boxes directly adjacent to the start point
        double adjusted_score = scores[i] / pow(lengths[i], 1.5);
        // if the adjusted score is lower than the current minimum, set the box as the candidate for the end point
        if (adjusted_score < min_score) {
            min_score = adjusted_score;
            min_index = i;
        }
    }

    // add all not yet processed start points from surface patches to the list of possible axis end points
    std::vector<AxisPoint> end_points;
    for (size_t i = start_idx + 1; i < grid.starting_points.size(); i++) {end_points.push_back(grid.starting_points[i]);}
    // add the end point with the lowest score to the list of end points if such a point could be found
    if (min_index > -1) {
        AxisPoint end = AxisPoint(min_index, grid.coordinates[min_index], "lowest score point");
        end_points.push_back(end);
    }

    // compute the axes between the current start point and all identified possible end points
    for (const AxisPoint &end: end_points) {
        Axis axis = Axis(start.type, end.type);
        axis.distance = distance(grid.centre(start.grid), grid.centre(end.grid));
        axis.length = lengths[end.index] * grid.box_length;
        axis.score = scores[end.index];
        size_t index = end.index;
        // start the trace back with end point and stop when the start point is reached
        while (index != start.index) {
            axis.axis.push_back(grid.centre(grid.coordinates[index]));
            index = trace_back[index];
        }
        // finally, add the start point
        axis.axis.push_back(grid.centre(grid.coordinates[start.index]));
        // add the axis to the pore/cavity
        grid.axes.push_back(axis);
    }
}

// generates an output file for each axis of a pore/cavity
void write_axes(PoreGrid &grid, const std::string &path_base) {
    // sort axis by distance between start and end point in descending order
    sort(grid.axes.begin(), grid.axes.end(),
            [](const Axis &lhs, const Axis &rhs) {return lhs.distance >= rhs.distance;});
    // generate an output file for each axis
    for (size_t i = 0; i < grid.axes.size(); i++) {
        // get the current axis
        const Axis &axis = grid.axes[i];
        // open the output file for the current axis
        std::ofstream file(path_base + std::to_string(i) + ".pdb");
        // write the header
        file << "# start (last point in the list): " << axis.start << std::endl;
        file << "# end (first point in the list): " << axis.end << std::endl;
        file << "# axis length: " << axis.length << " Angstrom" << std::endl;
        file << "# axis score: " << axis.score << std::endl;
        file << "# distance between axis start and end: " << axis.distance << " Angstrom" << std::endl;
        // write the boxes constituting the current axis
        size_t id = 0;
        for (const Vec<double> &coordinate: grid.axes[i].axis) {
            file << pseudo_pdb(id, coordinate) << std::endl;
            id++;
        }
        file.close();
    }
}

// main axis-trace routine
void axis_trace(const Settings &settings, const std::vector<PoreCluster> &clusters) {
    // compute the axis or axes of the given pores/cavities
    for (const PoreCluster &cluster: clusters) {
        auto start = std::chrono::high_resolution_clock::now();
        // construct the pore/cavity grid
        PoreGrid pore = PoreGrid(cluster, settings.perturb_value);
        // determine the starting point(s) of the axis (or axes)
        compute_starting_points(pore, settings.surface_patch_threshold);
        // trace the axis or axes
        for (size_t i = 0; i < pore.starting_points.size(); i++) {
            trace(pore, i);
        }
        // construct the base of the output file(s)
        std::string base = (settings.axis_dir / fs::path(settings.output_name + cluster.name() + "_axis_")).string();
        // generate output for each axis
        write_axes(pore, base);
        print(1, start, "> " + cluster.name() + ": computed axes in");
    }
}