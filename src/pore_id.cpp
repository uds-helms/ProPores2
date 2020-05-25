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

#include <cmath>
#include <stack>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <filesystem>
#include "grid.h"
#include "atom.h"
#include "basics.h"
#include "vector.h"
#include "reader.h"

namespace fs = std::filesystem;

// scan for empty grid boxes involved in PSP events on the given axis
std::vector<uint8_t> psp_scan(ProteinGrid &grid) {
    std::vector<uint8_t> psp = std::vector<uint8_t>(grid.size());

    for (const char axis: {'x', 'y', 'z'}) {
        // axes = {x index, y index, z index}
        // 0 = outer most loop, 1 = middle loop, 2 = internal loop
        // the given axis of interest is assigned index 2 and the other 0 and 1
        // default: x-axis
        Vec<int> axes = Vec<int>(2, 1, 0);
        if (axis == 'y') {
            axes = Vec<int>(0, 2, 1);
        } else if (axis == 'z') {
            axes = Vec<int>(1, 0, 2);
        }
        // maps the current index and last index
        std::vector<int> index(3);
        std::vector<int> start(3);
        std::vector<int> end(3);
        start[axes[0]] = grid.min.x;
        start[axes[1]] = grid.min.y;
        start[axes[2]] = grid.min.z;
        end[axes[0]] = grid.max.x;
        end[axes[1]] = grid.max.y;
        end[axes[2]] = grid.max.z;

        // iterate over the other two axes, hold them constant while checking for PSP events on the given axis
        for (index[0] = start[0]; index[0] < end[0]; index[0]++) {
            for (index[1] = start[1]; index[1] < end[1]; index[1]++) {
                // find the first non-empty box on the given axis
                int first_encounter = -1;
                for (index[2] = start[2]; index[2] < end[2]; index[2]++) {
                    if (grid.at(index[axes[0]], index[axes[1]], index[axes[2]]) == OCCUPIED) {
                        first_encounter = index[2];
                        break;
                    }
                }
                // no PSP event if there are only empty boxes or only a single non-empty box right at the end
                if (first_encounter == -1 || first_encounter >= end[2] - 1) continue;

                // start from the other side to see if there is at least one more non-empty box
                int second_encounter = -1;
                for (index[2] = end[2] - 1; index[2] > first_encounter; index[2]--) {
                    if (grid.at(index[axes[0]], index[axes[1]], index[axes[2]]) == OCCUPIED) {
                        second_encounter = index[2];
                        break;
                    }
                }
                // no PSP-event if there was no second encounter or the outer most non-empty boxes are next to each other
                if (second_encounter == -1 || second_encounter == first_encounter + 1) continue;

                // increase the PSP event counter of all empty boxes between first and second encounter
                for (index[2] = first_encounter + 1; index[2] < second_encounter; index[2]++) {
                    // get the 1D index of the current box and corresponding PSP counter and box state
                    size_t i = grid.index(index[axes[0]], index[axes[1]], index[axes[2]]);
                    uint8_t &psp_counter = psp[i];
                    BoxState &state = grid.at(i);
                    // increase the PSP counter if the box is empty
                    if (state != OCCUPIED) psp_counter++;
                    // mark the box as a potential pore as soon as the PSP counter reaches 2
                    if (psp_counter >= 2) state = POTENTIAL_PORE;
                }
            }
        }
    }
    return psp;
}

// check if the connecting vector between atom pairs is interrupted by occupied grid boxes
void collision_detection(ProteinGrid &grid) {
    double step = grid.box_radius * 0.25;
    // the smallest object that could fit between two atoms: H-atom, solvent radius or a grid box
    double const min = std::min(vdw_radius(H), grid.solvent_radius);

    // iterate over all atom pairs
    for (std::shared_ptr<Atom> &atom_1: grid.atoms) {
        for (std::shared_ptr<Atom> &atom_2: grid.atoms) {
            // make sure to only look at each atom pair once
            if (atom_1->id >= atom_2->id) continue;
            // vector connecting atom 1 and 2:   atom 1 -----------. atom 2
            Vec<double> v = atom_2->coord - atom_1->coord;
            // skip atoms pairs that are too close to each other
            // = the two vdW-radii + the minimum of the box length, prove radius and vdW-radius of H
            if (v.length() <= (atom_1->vdw + atom_2->vdw + min)) continue;

            // check if the vector from atom 1 to 2 passes through occupied grid boxes (= collisions)

            //                         full distance between atom 1 and 2
            //                  a1|----------------------------------------->|a2
            //                    |--->|------------------------------->|--->|
            //                     ^                  ^                   ^
            //   vdW-radius of atom 1 | length to search for collisions | vdW-radius of atom 2

            // The vdW radius of atom 1 or 2 can mark a grid box as occupied as long as it touches the centre point of
            // the box. As a result, starting directly outside of the vdW-radius could cause a collision detection
            // with atom 1 or 2, which is not intended. Adding/subtracting the box-radius makes sure that boxes
            // partly occupied by atom 1 or 2 are skipped.
            double search = atom_1->vdw + grid.box_radius;
            double search_end = v.length() - atom_2->vdw - grid.box_radius;

            // obtain the direction between atom 1 and 2 (= the connecting vector normed to unit length)
            // to move along the connecting vector in search of collisions
            Vec<double> v_unit = v.unit();

            bool collision = false;
            while (search < search_end) {
                // compute vector pointing from atom 1 to the current search point along the connecting vector
                // a1----.......a2
                // if the box is occupied, a collision was detected
                if (grid.at(atom_1->coord + v_unit * search) == OCCUPIED) {
                    collision = true;
                    break;
                }
                // advance the collision search by half a box if no collision was detected so far
                search += step;
            }
            // if a collision was detected, the atom pair is not involved in a pore
            if (collision) continue;
            // if there was no collision in between, store the atom pair for further information
            grid.add_atom_pair(atom_1, atom_2, v, v_unit);
        }
    }
}


// determine which cylinder function should be run, if true => run cylinder standalone
bool switch_to_standalone(const ProteinGrid &grid, const Settings &settings) {
    bool autodetect = grid.occupied.size() > settings.box_threshold;
    autodetect = autodetect || grid.pairs.size() < settings.atom_pair_threshold;
    autodetect = autodetect && settings.cylinder_tag == AUTODETECT;
    // run standalone cylinder perpendicularity if enabled
    return autodetect || settings.cylinder_tag == STANDALONE;
}

// form a cylinder around the connecting vector between atom pairs
void cylinder_trace(ProteinGrid &grid, const std::shared_ptr<AtomPair> &atom_pair, const double r_sq, const int r_grid,
                    const double threshold, const double search) {
    // location of the current search point
    Vec<double> search_point = atom_pair->atom_1->coord + atom_pair->unit * search;
    // grid box of the search point
    Vec<int> search_grid = grid.box(search_point);
    // search around the current search point
    for (const RelativeSearchBox &rel_box: grid.search_spaces[r_grid]) {
        Vec<int> box = search_grid + rel_box.grid;
        // make sure the search box is within the grid
        if (!grid.valid(box)) continue;
        // skip boxes when the box centre is not within range of the cylinder search radius AND they are not safe to add
        // do not mark this box as scanned yet, as moving along the search vector could bring the box centre
        // into range of the search radius
        if ((search_point - grid.centre(box)).length_squared() > r_sq && rel_box.length_sq > threshold) continue;
        // when a grid box is already part of a perpendicularly crossing we don't need to check it again
        // if the box is occupied, it cannot be part of the pore and can be skipped
        size_t index = grid.index(box);
        BoxState state = grid.at(index);
        if (state == POTENTIAL_PORE || state == OCCUPIED) continue;
        // trace the atom pair
        grid.trace[index].insert(atom_pair->id);
    }
}

// form a cylinder around the connecting vector between atom pairs
void cylinder_completion(ProteinGrid &grid) {
    grid.trace = std::vector<std::unordered_set<size_t>>(grid.size());
    for (const std::shared_ptr<AtomPair> &atom_pair: grid.pairs) {
        // get the search radius on the grid (atom_pair->r = smaller vdW-radius of atom_1 and atom_2)
        int r_grid = int(ceil(atom_pair->r / grid.box_length));
        // square the radius for easier comparisons in cylinder_trace
        double r_sq = pow(atom_pair->r, 2);
        //
        double threshold = pow(r_grid - grid.box_radius_sqrt3, 2) / grid.box_length_sq;
        // build the cylinder around the vector from atom 1 to 2

        //                         full distance between atom 1 and 2
        //                  a1|----------------------------------------->|a2
        //                    |-->|------------------------------------->|
        //                     ^                  ^
        //             box length |   length to trace the cylinder       |

        // start one box away from atom 1, search in box length increments along the vector until atom 2 is reached
        double search = grid.box_length;
        double search_length = atom_pair->vec.length();
        // compute the cylinder from near atom 1 to near atom 2
        while (search < search_length) {
            // check all boxes within the search radius of the search point if they are suitable for tracing
            cylinder_trace(grid, atom_pair, r_sq, r_grid, threshold, search);
            // advance the search
            search += grid.box_length;
        }
    }
}


// compute which atom pair cylinders are perpendicular to each other
void perpendicularity_cylinder(ProteinGrid &grid) {
    for (size_t i = 0; i < grid.size(); i++) {
        // skip boxes that are occupied or already marked as a potential pore
        BoxState &state = grid.at(i);
        if (state == POTENTIAL_PORE || state == OCCUPIED) continue;
        // skip boxes where less than two atom pairs cross
        const std::unordered_set<size_t> &trace = grid.trace[i];
        if (trace.size() < 2) continue;

        // check all combinations of atom pairs whose cylinders cross in this box
        for (const size_t id_1: trace) {
            for (const size_t id_2: trace) {
                // only look at each combination of pairs once
                if (id_1 >= id_2) continue;
                // if the atom pairs are almost perpendicular, mark the box as a potential pore and continue with
                // the next box
                if (almost_perpendicular(grid.pairs[id_1]->vec, grid.pairs[id_2]->vec)) {
                    state = POTENTIAL_PORE;
                    goto finished;
                }
            }
        }
        finished:
        continue;
    }
}

// compute which atom pair cylinders are perpendicular to each other
void perpendicularity_standalone(ProteinGrid &grid) {
    for (const std::shared_ptr<AtomPair> &pair_1: grid.pairs) {
        for (const std::shared_ptr<AtomPair> &pair_2: grid.pairs) {
            // only consider combinations once
            if (pair_1->id >= pair_2->id) continue;
            // only consider almost perpendicular combinations
            if (!almost_perpendicular(pair_1->vec, pair_2->vec)) continue;
            // compute the closest distance between the two vectors, as well as the corresponding points on each vector
            auto[point_1, point_2, dist] = closest(pair_1->atom_1->coord, pair_1->atom_2->coord,
                                                   pair_2->atom_1->coord, pair_2->atom_2->coord);

            // check if the two vectors are close enough that their cylinders overlap
            if (dist >= pair_1->r + pair_2->r) continue;

            // make sure the points are on the vectors between the atoms, not just on the line going through them
            if (!vec_in_interval(pair_1->atom_1->coord, point_1, pair_1->atom_2->coord)) continue;
            if (!vec_in_interval(pair_2->atom_1->coord, point_2, pair_2->atom_2->coord)) continue;

            // point_1 is the point on the vector of atom pair 1 that is the closest to the vector between atom pair 2,
            // where the closest point is point_2

            //                           pair_2->r
            //                   |<---------------------------|
            // point_1 --------------------------------------> point_2
            //                   |<-------p-------->|
            //                        r       r
            //        |---------------------------->|
            //           pair_1->r

            // the goal is to identify point p on the vector between point_1 and point_2 that is at the centre of the
            // overlap between the two cylinders
            double r;
            Vec<int> p;
            if (dist == 0) {
                r = std::min(pair_1->r, pair_2->r);
                // when the distance is 0, the two vectors cross exactly, which means point_1 == point_2
                p = grid.box(point_1);
            } else {
                Vec<double> unit = (point_2 - point_1).unit();
                r = (dist - (dist - pair_1->r) - (dist - pair_2->r)) / 2;
                p = grid.box(point_1 + unit * (pair_1->r - r));
            }

            for (const RelativeSearchBox &rel_box: grid.search_spaces[int(floor(r / grid.box_length))]) {
                Vec<int> search_box = p + rel_box.grid;
                // make sure the search box is within the grid bounds
                if (!grid.valid(search_box)) continue;
                // mark empty boxes in the overlap region as potential pore boxes
                BoxState &state = grid.at(search_box);
                if (state != OCCUPIED) state = POTENTIAL_PORE;
            }
        }
    }
}

// mark grid boxes as part of the background if they are accessible from the solvent without passing through occupied
// boxes or boxes that are (potentially) part of a pore
void seal_gaps(ProteinGrid &grid) {
    // step 1: mark grid boxes as part of the background if they are accessible from the solvent without passing
    //        through occupied boxes or boxes that are (potentially) part of a pore
    Vec<int> start_vec = grid.min;
    // use a stack to avoid recursion depth problems
    std::stack<Vec<int>> stack = std::stack<Vec<int>>();
    stack.push(start_vec);

    while (!stack.empty()) {
        Vec<int> coord = stack.top();
        stack.pop();
        // only consider boxes that are directly next to the current box, not diagonal
        for (const Vec<int> &search_box: grid.neighbours) {
            Vec<int> search = coord + search_box;
            // skip boxes that are not in the grid
            if (!grid.valid(search)) continue;

            BoxState &state = grid.at(search);
            if (state == UNCLASSIFIED) {
                state = BACKGROUND;
                stack.push(search);
            }
        }
    }
    // step 2: label all not yet classified boxes as potential pore boxes
    for (BoxState &state: grid.states) {
        if (state == UNCLASSIFIED) state = POTENTIAL_PORE;
    }
}

// remove shallow pore regions on the protein surface
void remove_shallow_regions(ProteinGrid &grid, const double shallow_radius) {
    // squared shallow radius length in grid coordinates
    double shallow_sq = pow(shallow_radius / grid.box_length, 2);
    // store shallow regions that should be made part of the background, if they are marked as background boxes
    // right away, they are included in the next iterations and thus a lot more than just shallow regions get
    // removed
    std::vector<uint8_t> make_background = std::vector<uint8_t>(grid.size());

    for (int x = grid.min.x; x < grid.max.x; x++) {
        for (int y = grid.min.y; y < grid.max.y; y++) {
            for (int z = grid.min.z; z < grid.max.z; z++) {
                Vec<int> box(x, y, z);
                // skip boxes that are not part of the background on the protein surface
                if (grid.at(x, y, z) != BACKGROUND) continue;
                // check if the box is a direct neighbour of an empty, non-background box
                bool remove = false;
                // the relative search box is diagonal, not directly next to the current box
                for (const Vec<int> &rel_box: grid.neighbours) {
                    Vec<int> search = box + rel_box;
                    // skip search boxes that are not in the grid
                    if (!grid.valid(search)) continue;
                    uint8_t state = grid.at(search);
                    if (state != BACKGROUND && state != OCCUPIED) {
                        remove = true;
                        break;
                    }
                }
                // do nothing if the box is only in contact with other background boxes or occupied boxes
                if (!remove) continue;

                std::stack<Vec<int>> stack;
                stack.push(box);

                while (!stack.empty()) {
                    Vec<int> current_box = stack.top();
                    stack.pop();
                    // otherwise search all surrounding boxes up within the shallow criterion to remove them
                    for (const Vec<int> &rel_box: grid.neighbours) {
                        Vec<int> search_box = current_box + rel_box;
                        if (!grid.valid(search_box)) continue;
                        size_t index = grid.index(search_box);
                        if (grid.at(index) == OCCUPIED) continue;
                        if (make_background[index]) continue;
                        if (distance_squared(box, search_box) <= shallow_sq) {
                            make_background[index] = 1;
                            stack.push(search_box);
                        }
                    }
                }
            }
        }
    }
    // make the removed shallow regions part of the background
    for (size_t i = 0; i < make_background.size(); i++) {
        if (make_background[i]) grid.at(i) = BACKGROUND;
    }
}

// determine the largest ball radius that can be centred in each empty grid box
void identify_pore_nuclei(ProteinGrid &grid) {
    grid.ball_radii = std::vector<double>(grid.size());
    for (int x = grid.min.x; x < grid.max.x; x++) {
        for (int y = grid.min.y; y < grid.max.y; y++) {
            for (int z = grid.min.z; z < grid.max.z; z++) {
                size_t index = grid.index(x, y, z);
                BoxState &state = grid.at(index);
                if (state != POTENTIAL_PORE) continue;
                // determine the distance of the box centre to the vdW-radius of the closest atom
                // ball radius is initialised to the largest possible distance in the grid
                auto[dist, _] = grid.closest(x, y, z);
                grid.ball_radii[index] = std::max(dist, 0.0);
                if (grid.ball_radii[index] >= grid.solvent_radius) state = NUCLEUS;
            }
        }
    }
}

// build initial clusters of connected potential pore boxes
void initial_pore_assembly(ProteinGrid &grid, const double tolerance) {
    auto box_tolerance = size_t(ceil(tolerance / grid.box_length));
    for (int x = grid.min.x; x < grid.max.x; x++) {
        for (int y = grid.min.y; y < grid.max.y; y++) {
            for (int z = grid.min.z; z < grid.max.z; z++) {
                // start a new pore/cavity cluster if the current box is not part of another cluster and is a
                // potential pore nucleus (= empty, marked as a potential pore and with a sufficiently large distance
                // to the van der Waals radius of the closest atom)
                if (grid.at(x, y, z) != NUCLEUS) continue;
                std::stack<std::pair<Vec<int>, size_t>> to_check;
                to_check.push(std::make_pair(Vec<int>(x, y, z), 0));
                PoreCluster &cluster = grid.new_cluster();

                while (!to_check.empty()) {
                    std::pair<Vec<int>, size_t> pair = to_check.top();
                    to_check.pop();
                    Vec<int> box = pair.first;
                    size_t away_from_nucleus = pair.second;
                    // add the current box to the current cluster
                    size_t index = grid.index(box);
                    cluster.add_box(index, box);
                    grid.at(index) = IN_CLUSTER;
                    // search all direct neighbours, not diagonal
                    for (const Vec<int> &rel_box: grid.neighbours) {
                        Vec<int> neighbour = box + rel_box;
                        // skip search boxes that are not in the grid
                        if (!grid.valid(neighbour)) continue;
                        // do not add boxes that are part of the background, already in a cluster or not a potential
                        // pore nucleus (= empty, marked as a potential pore and with a sufficiently large distance to
                        // the closest atom)
                        // otherwise add the box to the cluster and recursively continue the pore completion
                        BoxState state = grid.at(neighbour);
                        if (state == NUCLEUS) {
                            to_check.push(std::make_pair(neighbour, 0));
                        } else if (state == POTENTIAL_PORE && away_from_nucleus < box_tolerance) {
                            to_check.push(std::make_pair(neighbour, away_from_nucleus + 1));
                        }
                    }
                }
            }
        }
    }
}

// after the clusters are established, add the surrounding potential pore/cavity boxes to the nearest cluster
void complete_pores(ProteinGrid &grid) {
    double radius_sq = ceil(grid.solvent_radius_sq / grid.box_length_sq) + 1;

    for (PoreCluster &cluster: grid.clusters) {
        // create a copy of the boxes contained in the current pore/cavity cluster for iteration, since the vector is
        // going to be adjusted during the while loop, which would otherwise lead to unpredictable behaviour
        std::vector<PoreBox> boxes = cluster.boxes;
        // keep track of which boxes were already added to the current pore/cavity cluster
        std::vector<uint8_t> add_to_cluster(grid.size());
        // iteratively add
        for (const PoreBox &box: boxes) {
            std::stack<Vec<int>> stack;
            stack.push(box.coord);

            while (!stack.empty()) {
                Vec<int> current_box = stack.top();
                stack.pop();
                // otherwise search all surrounding boxes within the shallow radius to remove them
                for (const Vec<int> &rel_box: grid.neighbours) {
                    Vec<int> neighbour = current_box + rel_box;
                    // skip neighbours that are not in the grid
                    if (!grid.valid(neighbour)) continue;
                    size_t index = grid.index(neighbour);
                    // skip neighbours that are not potentially part of a pore or cavity
                    if (grid.at(index) != POTENTIAL_PORE) continue;
                    // skip neighbours that were already added to the cluster
                    if (add_to_cluster[index]) continue;
                    //
                    if (distance_squared(box.coord, neighbour) <= radius_sq) {
                        add_to_cluster[index] = 1;
                        grid.at(index) = IN_CLUSTER;
                        cluster.add_box(index, neighbour);
                        stack.push(neighbour);
                    }
                }
            }
        }
        for (PoreBox &box: cluster.boxes) { box.distance_to_closest = grid.ball_radii[grid.index(box.coord)]; }
    }
    // the complete ball radius vector is no longer needed
    grid.ball_radii.clear();
}

// classify clusters as pores or as cavities
void classify_clusters(ProteinGrid &grid) {
    for (PoreCluster &cluster: grid.clusters) {
        // iterate over each box in the current pore/cavity
        for (const PoreBox &box: cluster.boxes) {
            // check the surrounding boxes; only consider search boxes that are directly adjacent, not diagonal
            for (const Vec<int> &rel_box: grid.neighbours) {
                Vec<int> neighbour_box = box.coord + rel_box;
                // skip search boxes that are not in the grid
                if (!grid.valid(neighbour_box)) continue;
                // if the neighbouring box touches the background (= solvent), this takes precedent since it is
                // potentially part of a pore surface patch and thus a potential start/end point of a pore axis
                if (grid.at(neighbour_box) == BACKGROUND) {
                    cluster.pore = true;
                    break;
                }
            }
        }
    }
}

// for each empty grid box involved in a pore/cavity, compute if it is at the protein surface, in touch with
// the protein, another pore/cavity cluster, unclustered potential pore/cavity boxes or inside the pore/cavity cluster
void label_pore_boxes(ProteinGrid &grid) {
    for (PoreCluster &cluster: grid.clusters) {
        // iterate over each box in the current pore/cavity
        for (PoreBox &box: cluster.boxes) {
            // compute the distance of its centre to the protein's centre of mass
            box.distance_to_centre = (grid.centre_of_mass - grid.centre(box.coord)).length();
            // initially label the box as being located on the inside of the pore/cavity cluster
            BoxState label = ON_CLUSTER_INSIDE;
            // check the surrounding boxes; only consider search boxes that are directly adjacent, not diagonal
            for (const Vec<int> &rel_box: grid.neighbours) {
                Vec<int> neighbour_box = box.coord + rel_box;
                // skip search boxes that are not in the grid
                if (!grid.valid(neighbour_box)) continue;
                // obtain the state of the neighbouring box
                BoxState state = grid.at(neighbour_box);
                // if the neighbouring box touches the background (= solvent), this takes precedent since it is
                // potentially part of a pore surface patch and thus a potential start/end point of a pore axis
                if (state == BACKGROUND) {
                    cluster.pore = true;
                    label = TOUCHES_BACKGROUND;
                    break;
                    // check for other potential types of neighbour boxes in order of their priority
                } else if (state == OCCUPIED) {
                    label = TOUCHES_PROTEIN;
                } else if (state == IN_CLUSTER && label != TOUCHES_PROTEIN) {
                    if (!cluster.has_index(grid.index(neighbour_box))) label = TOUCHES_OTHER_CLUSTER;
                } else if ((state == POTENTIAL_PORE || state == NUCLEUS) && label != TOUCHES_PROTEIN
                           && label != TOUCHES_OTHER_CLUSTER) {
                    label = TOUCHES_POTENTIAL_PORE;
                }
            }
            // finally set the label of the current box
            box.state = label;
        }
    }
}

// for each pore/cavity cluster add the surrounding residues
void lining_residues(ProteinGrid &grid) {
    for (PoreCluster &cluster: grid.clusters) {
        for (const PoreBox &box: cluster.boxes) {
            for (const Vec<int> &rel_box: grid.neighbours) {
                Vec<int> search_box = box.coord + rel_box;
                // skip search boxes that are not in the grid
                if (!grid.valid(search_box)) continue;
                for (const size_t atom_id: grid.occupied[grid.index(search_box)]) {
                    cluster.add_lining_residue(grid.atoms[atom_id]);
                }
            }
        }
    }
}

// create files for pore/cavity boxes and lining residues
void output_pores(ProteinGrid &grid, const Settings &settings) {
    // for each pore/cavity, create a PDB file with pore/cavity grid boxes and a file for lining residues
    for (const PoreCluster &cluster: grid.clusters) {
        // file name basis
        std::string name = settings.output_name + cluster.name();
        // PDB
        // pore/cavity volume in Angstrom
        double volume = cluster.size() * grid.box_volume;
        // create the PDB file path, open the empty file and write the pore/cavity volume
        std::ofstream pdb_file((settings.pore_dir / fs::path(name + ".pdb")).string());
        pdb_file << "REMARK\tPore Volume = " << volume << " (cubic Angstrom)" << std::endl;
        // write each pore/cavity grid box as an "atom" in pseudo PDB format so that it can be visualised
        for (size_t i = 0; i < cluster.size(); i++) {
            pdb_file << pseudo_pdb(i + 1, grid.centre(cluster.boxes[i].coord)) << std::endl;
        }
        pdb_file.close();
        // LINING RESIDUES
        // create the lining residues file path and open the empty file
        std::ofstream lining((settings.lining_dir / fs::path(name + ".tsv")).string());
        // for each lining residue, write the chain, residue ID and residue type in tab-separated format
        for (const auto &[chain, res_id, res_type]: cluster.lining_residues) {
            lining << chain << "\t" << res_id << "\t" << to_str(res_type) << std::endl;
        }
        lining.close();
    }
}

// main pore/cavity identification routine
void pore_ID(ProteinGrid &grid, const Settings &settings) {
    auto start = std::chrono::high_resolution_clock::now();
    fs::path log = settings.pore_log;

    // PROTEIN-SOLVENT-PROTEIN IDENTIFICATION
    add_comment(log, 1, "identifying enclosed, empty spaces");
    // scan the x, y and z axis for protein-solvent-protein (PSP) events
    psp_scan(grid);
    add_entry(log, 1, "PSP scan", start);
    // first run collision detection
    auto sub_start = std::chrono::high_resolution_clock::now();
    collision_detection(grid);
    add_entry(log, 1, "atom pairs", grid.pairs.size());
    add_entry(log, 1, "collision detection", sub_start);
    // compute
    sub_start = std::chrono::high_resolution_clock::now();
    if (switch_to_standalone(grid, settings)) {
        add_entry(log, 1, "used computation mode", to_str(STANDALONE));
        perpendicularity_standalone(grid);
    // otherwise run cylinder trace
    } else {
        add_entry(log, 1, "used computation mode", to_str(RAY_TRACE));
        // first generate and trace cylinders in the grid
        cylinder_completion(grid);
        // then compute which boxes are crossed by perpendicular cylinders
        perpendicularity_cylinder(grid);
        // release the trace memory
        grid.trace = std::vector<std::unordered_set<size_t>>();
    }
    add_entry(log, 1, "cylinder computation", sub_start);
    // release the atom pair memory
    grid.pairs = std::vector<std::shared_ptr<AtomPair>>();

    // PROTEIN SURFACE
    add_comment(log, 1, "processing protein surface");
    // 1. seal gaps on the protein surface to compute the background surrounding the protein
    // 2. remove shallow regions on the protein surface
    sub_start = std::chrono::high_resolution_clock::now();
    seal_gaps(grid);
    add_entry(log, 1, "gap sealing", sub_start);
    sub_start = std::chrono::high_resolution_clock::now();
    remove_shallow_regions(grid, settings.probe_radius);
    add_entry(log, 1, "shallow region removal", sub_start);
    print(1, start, "> marked potential pore or cavity areas in");

    // PORE IDENTIFICATION AND CLASSIFICATION
    add_comment(log, 1, "identifying pores");
    // 1. identify nuclei that could form the core of a pore or cavity
    // 2. merge all connected pore/cavity nucleus boxes into clusters (this avoids merging pores/cavities that are
    //    close to each other and connected by narrow stretches of empty boxes)
    // 3. after the clusters are established, add the surrounding potential pore/cavity boxes to the nearest cluster
    start = std::chrono::high_resolution_clock::now();
    sub_start = std::chrono::high_resolution_clock::now();
    identify_pore_nuclei(grid);
    add_entry(log, 1, "pore nuclei identification", sub_start);
    sub_start = std::chrono::high_resolution_clock::now();
    initial_pore_assembly(grid, 1.0);
    add_entry(log, 1, "initial pore assembly", sub_start);
    sub_start = std::chrono::high_resolution_clock::now();
    complete_pores(grid);
    add_entry(log, 1, "pore completion", sub_start);

    // POST-PROCESSING
    add_comment(log, 1, "pore post-processing");
    // 1. for each empty grid box involved in a pore/cavity, compute if it's at the protein surface, in touch with
    //    residues, in touch with another cluster or in the cluster centre
    // 2. if the cluster touches the background, it is a pore, otherwise it is a cavity
    sub_start = std::chrono::high_resolution_clock::now();
    if (settings.run_axis_trace || settings.run_axis_preparation) label_pore_boxes(grid);
    classify_clusters(grid);
    add_entry(log, 1, "pore classification", sub_start);
    // remove pores/cavities that are too small, or only keep pores or cavities
    sub_start = std::chrono::high_resolution_clock::now();
    grid.filter_pores(settings.volume_threshold, settings.remove_tag);
    add_entry(log, 1, "pore filter", sub_start);
    // compute the residues surrounding each pore/cavity
    sub_start = std::chrono::high_resolution_clock::now();
    lining_residues(grid);
    grid.occupied = std::vector<std::unordered_set<size_t>>();
    add_entry(log, 1, "lining residue", sub_start);
    // create files for the pores/cavities and corresponding lining residues
    sub_start = std::chrono::high_resolution_clock::now();
    output_pores(grid, settings);
    add_entry(log, 1, "pore output", sub_start);
    print(1, start, "> identified " + std::to_string(grid.clusters.size()) + " pores/cavities in");
    add_comment(log, 1, "results");
    add_entry(log, 1, "identified pores", grid.clusters.size());
}
