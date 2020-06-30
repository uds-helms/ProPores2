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

#include "grid.h"
#include <set>
#include <fstream>
#include "gate.h"
#include "settings.h"
#include "basics.h"
#include "atom.h"
#include <stack>
#include "reader.h"
#include <filesystem>
#include <unordered_map>
#include <iostream>
#include "angles.h"

namespace fs = std::filesystem;

// estimate the difficulty of a gate by considering the potential number of rotamer combinations
GateDifficulty estimate_difficulty(const Gate &gate, const bool check_for_empty) {
    size_t residues_with_rotamers = gate.residues.size();
    size_t hard = 0;
    size_t medium = 0;
    unsigned long long int combinations = 1;
    // compute the potential number of rotamer combinations
    for (const Residue &res: gate.residues) {
        bool valid = true;
        switch (res.type) {
            case ARG:
                combinations *= 14086;
                hard++;
                break;
            case ASN:
                combinations *= 130;
                break;
            case ASP:
                combinations *= 64;
                break;
            case CYS:
                combinations *= 12;
                break;
            case GLN:
                combinations *= 1316;
                medium++;
                break;
            case GLU:
                combinations *= 641;
                medium++;
                break;
            case HIS:
                combinations *= 142;
                break;
            case ILE:
                combinations *= 124;
                break;
            case LEU:
                combinations *= 104;
                break;
            case LYS:
                combinations *= 14741;
                hard++;
                break;
            case MET:
                combinations *= 1297;
                medium++;
                break;
            case PHE:
                combinations *= 64;
                break;
            case SER:
                combinations *= 12;
                break;
            case THR:
                combinations *= 12;
                break;
            case TRP:
                combinations *= 127;
                break;
            case TYR:
                combinations *= 58;
                break;
            case VAL:
                combinations *= 12;
                break;
            default:
                residues_with_rotamers--;
                valid = false;
        }
        // if the residue is a valid rotating residue but does not have any rotamers left, it is also invalid
        if (valid && check_for_empty && res.rotamers.empty()) residues_with_rotamers--;
    }
    // if there is no residue that can be rotated, there can be no gate between the two pores/cavities
    if (residues_with_rotamers == 0) return DIFFICULTY_ERROR;
    // trivial if there is only one residue with (potential) rotamers
    if (residues_with_rotamers == 1) return EASY;
    // likely very long runtime if there is more than one hard residue
    if (hard > 1) return HARD;
    // if there is only one hard residue, the difficulty class can vary a fair amount
    if (hard == 1) {
        // if there is only one other residue that is neither hard nor medium, then the runtime is probably low
        if (medium == 0 && residues_with_rotamers < 6) return EASY;
        // if there is only one medium residue and at most two easy residues, then the runtime is probably still fine
        if (medium == 1 && residues_with_rotamers < 4) return MEDIUM;
        // if the number of easy or medium residues is higher, the runtime is probably very long
        return HARD;
    }
    // check different cases with two or more medium residues and no hard residues
    if (medium > 4) return HARD;
    if (medium >= 2) {
        if (residues_with_rotamers < 6) return MEDIUM;
        return HARD;
    }
    // if there are only easy residues and at most one medium residue, the difficulty depends on the number of residues
    if (residues_with_rotamers < 6) return EASY;
    if (residues_with_rotamers < 9) return MEDIUM;
    return HARD;
}

// extract potential gates between neighbouring pores/cavities from the protein grid
void gates_from_grid(ProteinGrid &grid, std::vector<Gate> &gates) {
    // check all unique combinations of pores/cavities
    for (const PoreCluster &cluster_1: grid.clusters) {
        for (const PoreCluster &cluster_2: grid.clusters) {
            // skip combinations that were already examined
            if (cluster_1.id >= cluster_2.id) continue;
            // identify all lining residues that are shared by the two pores/cavities
            std::vector<std::tuple<std::string, size_t, std::string>> shared;
            for (const std::tuple<std::string, size_t, std::string> &residue: cluster_1.lining_residues) {
                ResidueType res = to_residue_type(std::get<2>(residue));
                // skip lining residues that are not standard amino acids or cannot be rotated
                if (res == INVALID_RESIDUE || res == PRO || res == GLY || res == ALA) continue;
                // check if the residue is shared
                if (cluster_2.lining_residues.find(residue) != cluster_2.lining_residues.end()) {
                    shared.push_back(residue);
                }
            }
            // the two pores/cavities do not share lining residues and therefore do not have a potential gate
            if (shared.empty()) continue;
            // create the gate and estimate its difficulty
            Gate gate = Gate(cluster_1.id, cluster_2.id, shared);
            gate.difficulty = estimate_difficulty(gate, false);
            // DIFFICULTY_ERROR means that there is no way of opening the gate, e.g. because none of the gating
            // residues have rotamers like in the case of glycine
            if (gate.difficulty < DIFFICULTY_ERROR) gates.push_back(std::move(gate));
        }
    }
}

// extract a single gate from a single input file
void gate_from_file(const std::string &file_path, std::vector<Gate> &gates) {
    std::ifstream file(file_path);
    std::string line;
    // get the first line, which contains the IDs of the two neighbouring pores/cavities
    std::getline(file, line);
    std::vector<std::string> cols = split(r_strip(line));
    // the ID line is incomplete
    if (cols.size() < 6) {
        file.close();
        return;
    }
    size_t id_1 = stoi(cols[3]);
    size_t id_2 = stoi(cols[5]);
    // extract the estimated gate difficulty and check if it is below the difficulty threshold
    std::getline(file, line);
    cols = split(r_strip(line));
    // the difficulty line is incomplete
    if (cols.size() < 3) {
        file.close();
        return;
    }
    GateDifficulty difficulty = to_gate_difficulty(cols[3]);
    // extract the gating residues
    std::vector<std::tuple<std::string, size_t, std::string>> residues;
    while (std::getline(file, line)) {
        cols = split(r_strip(line), '\t');
        if (cols.size() != 3) continue;
        residues.emplace_back(std::make_tuple(cols[0], stoi(cols[1]), cols[2]));
    }
    file.close();
    // make sure that there is at least one residue in the input file
    if (residues.empty()) return;
    // create the gate
    Gate gate = Gate(id_1, id_2, residues);
    // set the gate difficulty and add it to the list of gates
    gate.difficulty = difficulty;
    gates.push_back(gate);
}

// extract gates from a directory of input files
void gates_from_directory(const std::string &dir_path, std::vector<Gate> &gates) {
    for (auto &p: fs::directory_iterator(dir_path)) {
        if (fs::is_directory(p)) continue;
        gate_from_file(p.path().string(), gates);
    }
}

// create an output file for each gate with the pore/cavity IDs, estimated difficulty and gating residues
void output_gates(const Settings &settings, const std::vector<Gate> &gates) {
    for (const Gate &gate: gates) {
        // do not output gates that do not have valid gating residues or where the difficulty could not be estimated
        if (gate.difficulty > HARD || gate.residues.empty()) continue;
        // generate the output path and create the output file
        fs::path output_path = settings.gate_preparation_dir / fs::path("between_" + gate.name() + ".txt");
        std::ofstream file(output_path.string());
        // write a header with the IDs of the two pores/cavities and the estimated difficulty
        file << "# Pore/cavity IDs: " << gate.pore_1 << " and " << gate.pore_2 << std::endl;
        file << "# Estimated difficulty: " << to_str(gate.difficulty) << std::endl;
        // write the residue chain, ID and type for each gating residue
        for (const Residue &res: gate.residues) {
            file << res.chain << "\t" << res.pdb_id << "\t" << to_str(res.type) << std::endl;
        }
        file.close();
    }
}

// generate the library with valid rotamer angles for each residue
void load_rotamer_library(const Settings &cfg, std::unordered_map<ResidueType, std::vector<std::vector<int>>> &lib) {
    // read in the rotamer file of each residue type
    for (const std::string &res_name: {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
                                       "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL"}) {
        // construct the file path and check if it's pointing towards an existing file
        std::string file_path = cfg.rotamer_file_path(res_name);
        if (file_path.empty()) continue;
        // open the file and read the lines, each representing a rotamer angle configuration for this residue
        std::ifstream file(file_path);
        std::string line;
        std::vector<std::vector<int>> rotamer_definitions;
        while (std::getline(file, line)) {
            std::vector<int> angles;
            // extract the rotamer angles
            for (std::string angle: split(line, '\t')) {
                angle = remove_spaces(angle);
                if (!angle.empty()) angles.push_back(stoi(angle));
            }
            rotamer_definitions.push_back(angles);
        }
        // insert the residue rotamer information into the library
        lib.insert(std::make_pair(to_residue_type(res_name), rotamer_definitions));
        file.close();
    }
}

// compute the dihedral angles of a given gate residue
void compute_dihedral_angles(const Gate &gate, Residue &res) {
    // res.defs = indices of the main backbone and side chain components (e.g. N-CA-CB-CG-CD)
    // the loop computes the dihedral angles along these main components
    for (size_t i = 0; i < res.defs.size() - 3; i++) {
        double dihedral = dihedral_angle(res.atoms[res.defs[i].index].coord,
                                         res.atoms[res.defs[i + 1].index].coord,
                                         res.atoms[res.defs[i + 2].index].coord,
                                         res.atoms[res.defs[i + 3].index].coord);
        // store the dihedral angle, the angle index and the indices of the relevant components for later steps
        res.angles.emplace_back(i, res.defs[i + 1].index, res.defs[i + 2].index, radian_to_degree(dihedral));
        // continue with the next section B-C-D-E and the angle between B ---> C and D ---> E
    }
}

// rotate the side chain of a given residue to generate the specified rotamer
void side_chain_rotate(const Residue &res, Rotamer &rotamer) {
    // residues have up to 4 dihedral angles along the main backbone and side chain components (e.g. N-CA-CB-CG-CD)
    // the loop iterates over these angles starting with the angle of the entire side chain and then the angles of every
    // smaller parts of the side chain
    // the first rotation point is C-beta (CB), the next is C-gamma (CG),...
    // atoms before the current rotation point in the residue chain are not rotated
    for (const Dihedral &chi: res.angles) {
        // chi.dihedral is the current dihedral angle and rotamer.angle is the desired dihedral angle
        // at this point of the side chain
        // => the difference is by how much the side chain starting at this point has to be rotated
        double rotation_angle = rotamer.angles[chi.index] - chi.angle;
        // do not bother to rotate if the difference is negligible
        if (abs(rotation_angle) < 1) continue;
        // the current rotation point, e.g. C-beta in the first iteration (N-CA-CB-CG)
        Vec<double> rotation_point = rotamer.atoms[chi.rotation_point].coord;
        // compute the rotation axis
        Vec<double> rotation_axis = rotation_point - rotamer.atoms[chi.rotation_axis_start].coord;
        // construct the quaternion rotation matrix
        auto[row_1, row_2, row_3] = quaternion_rotation_matrix(rotation_angle, rotation_axis);
        // use the rotation matrix and rotation point to rotate the position of all atoms subsequent to the current
        // rotation centre atom in the side chain
        for (size_t i = chi.rotation_point + 1; i < rotamer.atoms.size(); i++) {
            BasicAtom &atom = rotamer.atoms[i];
            // rotation point -------> atom position
            Vec<double> temp = atom.coord - rotation_point;
            // new position of the atom relative to the rotation point
            Vec<double> rotated(dot_product(row_1, temp), dot_product(row_2, temp), dot_product(row_3, temp));
            // => add the rotation point to the new position for the new position of the atom
            atom.coord = rotated + rotation_point;
        }
    }
}

// check if a residue side chain rotamer clashes substantially with the rest of the protein
bool serious_clash(const Gate &gate, const Residue &res, Rotamer &rotamer, const double tolerance,
                   const double gamma_tolerance) {
    double energy = 0;
    // check each rotamer side chain atom for a clash with the rest of the protein
    for (BasicAtom &atom: rotamer.atoms) {
        // skip all atoms that do not get rotated, i.e. components of the backbone or the C-beta atom
        if (!atom.is_side_chain || atom.type == "CB") continue;
        // find the closest atoms and compute the distance
        auto[dist, closest] = KD_nearest_neighbour(gate.root, atom.coord);

        // vdW radius |          | vdw radius
        // atom------>|--------->|<------closest
        //    |----------------->|<-------
        //                ^ dist

        // case: no overlap between the two van der Waals radii
        //       => dist > 0 and dist - atom.vdw > 0
        // case: atom within the van der Waals radius of the closest atom
        //       => dist < 0 and dist - atom.vdw < 0
        // case: overlapping van der Waals radii
        //       => dist > 0 and dist - atom.vdw < 0

        // tolerate a certain amount of overlap that can be set by the user
        double clash_tolerance = tolerance;
        // if the closest atom is the C-alpha atom of the same residue AND the current atom is in the gamma position
        // of the side chain, use a fixed tolerance instead that better matches the increased overlap of this case
        if (res.pdb_id == closest->residue_id && closest->full_type == "CA" && atom.type.size() >= 2 &&
            atom.type[1] == 'G')
            clash_tolerance = gamma_tolerance;
        // the more negative dist - atom.vdw, the bigger the overlap/clash
        // => if the clash is more negative than the negative clash tolerance it is too severe to accept the rotamer
        dist -= atom.vdw;
        if (dist < -clash_tolerance) {
            return true;
        }
        // the self-energy of a rotamer is comprised of the distance to the pore axis and the distance to the
        // centre of the gate: E_self = E_sphere - E_CA
        // E_sphere = sum of the distance to the closest van der Waals radius for each rotamer side chain atom
        //      > the larger E_sphere, the closer the rotamer is to the pore axis and the more it is blocking the gate
        //      > the smaller E_sphere, the further away the rotamer is from the pore axis and the less it is blocking
        //        (the same applies the more negative the distance, i.e. the bigger the overlap between the two
        //         van der Waals radii)
        // E_CA = sum of the distances to the centre of the gate (= C_alpha centre) for each rotamer side chain atom
        //      > the smaller E_CA, the closer the rotamer is to the gate centre and the more it is blocking the gate
        //      > the larger E_CA, the further away the rotamer is from the gate centre and the less it is blocking
        // => the smaller E_self, the less it is blocking the gate and thus the better the rotamer
        energy += dist;
        if (gate.residues.size() > 1) energy -= distance(gate.centre, atom.coord);
    }
    rotamer.energy = energy;
    return false;
}

// generate valid rotamers for a given gate residue
void generate_rotamers(const Settings &settings, const Gate &gate, Residue &res,
                       const std::vector<std::vector<int>> &rotamer_angles) {
    size_t id = 0;
    for (const std::vector<int> &angles: rotamer_angles) {
        Rotamer rotamer(id, res.atoms, angles);
        // the number of desired dihedral angles must match the number of actual dihedral angles
        if (rotamer.angles.size() != res.angles.size()) continue;
        // rotate the side chain to the desired dihedral angles
        side_chain_rotate(res, rotamer);
        // discard rotamers that lead to a serious clash with the rest of the protein
        if (serious_clash(gate, res, rotamer, settings.clash_tolerance, settings.clash_tolerance_gamma)) continue;
        // add valid rotamers to the residue
        res.rotamers.push_back(rotamer);
        id++;
    }
}

// check how many residues do not have valid rotamers, and add to them the original residue configuration
size_t fix_residues(Gate &gate) {
    size_t empty = 0;
    for (Residue &res: gate.residues) {
        if (res.rotamers.empty()) {
            res.rotamers.emplace_back(0, res.atoms, std::vector<int>(res.angles.size(), 0));
            empty++;
        }
    }
    return empty;
}

// sort residues by the number of rotamers in ascending order and adjust their indices accordingly
void sort_residues(Gate &gate) {
    sort(gate.residues.begin(), gate.residues.end(),
         [](const Residue &lhs, const Residue &rhs) { return lhs.rotamers.size() < rhs.rotamers.size(); });
}

// compute the interaction energy of two rotamers based on the distance between them
double interaction_energy(const Rotamer &rot_1, const Rotamer &rot_2) {
    double energy = 0;
    // go over all combinations of side chain atoms between the two rotamers
    for (const BasicAtom &atom_1: rot_1.atoms) {
        if (!atom_1.is_side_chain) continue;
        for (const BasicAtom &atom_2: rot_2.atoms) {
            if (!atom_2.is_side_chain) continue;

            //   vdW radius |          | vdw radius
            // atom 1------>|--------->|<------atom 2
            //      |------------------------->
            //                ^ dist

            // case: no overlap between the two van der Waals radii
            //       => dist > 0 and dist - atom_1.vdw - atom_2.vdw > 0
            // case: overlapping van der Waals radii
            //       => dist >= 0 and dist - atom_1.vdw - atom_2.vdw < 0

            // E_interaction(a, b) = sum of distances between all pairs of atoms in rotamer 1 and 2
            //      > overlaps lead to an increase of the energy, and the bigger the overlap, the bigger the increase
            //      > no overlap leads to a decrease of the energy, and the bigger the distance, the bigger the decrease
            // => the further away the two rotamers, the more open the gate and the smaller the energy
            energy -= distance(atom_1.coord, atom_2.coord) - atom_1.vdw - atom_2.vdw;
        }
    }
    return energy;
}

// eliminates rotamers that are dead-ends in the effort to minimise the gate energy
void dead_end_elimination(Gate &gate) {
    // iterate over all residues in the gate
    for (Residue &current_residue: gate.residues) {
        // iterate over all unique pairs of rotamers in the current residue
        for (Rotamer &rot_1: current_residue.rotamers) {
            // skip already eliminated rotamers
            if (rot_1.eliminated) continue;
            for (Rotamer &rot_2: current_residue.rotamers) {
                // skip pairs that were already examined
                if (rot_2.id <= rot_1.id) continue;
                // skip already eliminated rotamers
                if (rot_2.eliminated) continue;
                // determine which of the two rotamers cannot have a lower energy than the other and eliminate it
                double min_difference_sum_1 = 0;
                double min_difference_sum_2 = 0;
                // go over all other residues and their rotamers to check which of the rotamers has a lower energy
                for (const Residue &other_residue: gate.residues) {
                    // skip the current residue
                    if (other_residue.id == current_residue.id) continue;
                    double min_difference_1 = INFINITY;
                    double min_difference_2 = INFINITY;
                    for (const Rotamer &other_rot: other_residue.rotamers) {
                        // skip eliminated rotamers
                        if (other_rot.eliminated) continue;
                        // compute the interaction energy of the two rotamers with the rotamer of the other residue
                        // then compute the difference to see which of the two rotamers has the lower interaction energy
                        double difference = interaction_energy(rot_1, other_rot) - interaction_energy(rot_2, other_rot);
                        // determine the minimum difference for the other residue
                        min_difference_1 = std::min(difference, min_difference_1);
                        min_difference_2 = std::min(-difference, min_difference_2);
                    }
                    // add the minimum difference for the other residue to the sum
                    min_difference_sum_1 += min_difference_1;
                    min_difference_sum_2 += min_difference_1;
                }
                // compute the elimination criterion based on the difference between the rotamer self-energies and
                // the minimum interaction energy difference
                double elimination_criterion_1 = rot_1.energy - rot_2.energy + min_difference_sum_1;
                double elimination_criterion_2 = rot_2.energy - rot_1.energy + min_difference_sum_2;
                // check which of the two rotamers can be eliminated
                if (elimination_criterion_1 > 0) {
                    rot_1.eliminated = true;
                } else if (elimination_criterion_2 > 0) {
                    rot_2.eliminated = true;
                }
            }
        }
    }
}

// remove rotamers that have been eliminated
void remove_eliminated(Gate &gate) {
    for (Residue &res: gate.residues) {
        std::vector<Rotamer> not_eliminated;
        for (Rotamer &rot: res.rotamers) { if (!rot.eliminated) not_eliminated.push_back(std::move(rot)); }
        std::swap(res.rotamers, not_eliminated);
    }
}

// check all combinations of the remaining rotamers to determine the combination that yields the lowest gate energy
std::vector<Rotamer> exhaustive_search(const Gate &gate) {
    // initialise the gate energy and stack, since an iterative approach circumvents potential recursion limit issues
    double gate_energy = INFINITY;
    std::vector<Rotamer> best_combination;
    std::stack<std::vector<Rotamer>> stack;
    stack.push(std::vector<Rotamer>());
    // keep adding rotamer combinations to the stack until all possible combinations have been evaluated
    while (!stack.empty()) {
        // retrieve the current combination from the top of the stack
        std::vector<Rotamer> current_combination = stack.top();
        stack.pop();
        // iterate over all rotamers of the next residue in the gate to extend the current combination
        for (const Rotamer &rotamer: gate.residues[current_combination.size()].rotamers) {
            // extend the current combination with the current rotamer
            std::vector<Rotamer> copy = current_combination;
            copy.push_back(rotamer);
            // if the combination is not fully extended yet, add it to the stack
            if (copy.size() != gate.residues.size()) {
                stack.push(copy);
                continue;
            }
            // if the combination is fully extended, compute its energy
            double current_energy = 0;
            // compute the interaction energy of the rotamers in the current combination
            for (size_t i = 0; i < copy.size(); i++) {
                for (size_t j = i + 1; j < copy.size(); j++) {
                    current_energy += interaction_energy(copy[i], copy[j]);
                }
            }
            // compute the self energies of the rotamers in the current combination
            for (Rotamer &rot: copy) { current_energy += rot.energy; }
            // the rotamer combination is more optimal if its energy is less than the current best combination
            if (current_energy < gate_energy) {
                gate_energy = current_energy;
                best_combination = copy;
            }
        }
    }
    return best_combination;
}

// main gate-open routine
void gate_open(Settings &settings, std::vector<Gate> &gates) {
    // load the rotamer library
    auto library_start = std::chrono::high_resolution_clock::now();
    std::unordered_map<ResidueType, std::vector<std::vector<int>>> rotamer_lib;
    load_rotamer_library(settings, rotamer_lib);
    print(1, library_start, "> loaded the rotamer library in");
    // load the PDB to log all changed atom coordinates and changed residues
    std::vector<std::shared_ptr<Atom>> global_atoms;
    PDBReaderStats stats = PDBReaderStats();
    parse_PDB(settings, global_atoms, stats);
    std::set<std::string> global_residues;
    // process each gate
    for (Gate &gate: gates) {
        // skip the gate if the difficulty is too high and a re-assessment is not desired
        if (!settings.re_estimate && gate.difficulty >= settings.difficulty_threshold) {
            print(1, "> gate between pores " + std::to_string(gate.pore_1) + " and "
                     + std::to_string(gate.pore_2) + ": difficulty higher than threshold");
            continue;
        }
        // load the protein atoms into the gate
        print(1, "> gate between pores " + std::to_string(gate.pore_1) + " and "
                 + std::to_string(gate.pore_2) + ": " + std::to_string(gate.residues.size()) + " residue(s)");
        auto start_gate = std::chrono::high_resolution_clock::now();
        gate.load_atoms(settings);
        // make sure that there is at least one gate residue with valid rotamers
        if (gate.residues.empty()) {
            print(2, "> no residues with valid rotamers");
            continue;
        }
        // compute dihedral angles in the side chains of the gate residues and generate the rotamers
        for (Residue &res: gate.residues) {
            compute_dihedral_angles(gate, res);
            generate_rotamers(settings, gate, res, rotamer_lib.at(res.type));
        }
        // re-estimate the difficulty now that rotamers have been eliminated and skip the ones that are more difficult
        // than the difficulty threshold
        if (settings.re_estimate) {
            gate.difficulty = estimate_difficulty(gate, true);
            if (gate.difficulty == DIFFICULTY_ERROR) {
                print(2, "> no residues with valid rotamers");
                continue;
            }
            if (gate.difficulty >= settings.difficulty_threshold) {
                print(2, "> re-assessed difficulty still higher than threshold");
                continue;
            }
        }
        // make sure that there is at least one gate residue with valid rotamers
        if (gate.residues.size() == fix_residues(gate)) {
            print(2, "> no residues with valid rotamers");
            continue;
        }
        // eliminate rotamers that cannot be part of an optimal solution with dead-end elimination
        // sort the gating residues according to the number of rotamers first to improve the runtime
        sort_residues(gate);
        dead_end_elimination(gate);
        // remove eliminated rotamers from the gate and sort the residues again by the number of rotamers
        remove_eliminated(gate);
        // compute the total energy of all possible gate configurations based on the remaining rotamers
        std::vector<Rotamer> best_combination = exhaustive_search(gate);
        // apply the positions of the best rotamer combination and write the newly positioned atoms in a PDB file
        // also extract and store the chain and ID of the changed residues for visualisation on the web server
        std::set<std::string> gate_residues;
        for (const Rotamer &rot: best_combination) {
            for (const BasicAtom &atom: rot.atoms) {
                // identify the original atom
                std::shared_ptr<Atom> gate_atom = gate.atoms[atom.index];
                std::shared_ptr<Atom> global_atom = global_atoms[atom.index];
                // update atom coordinates
                gate_atom->coord = atom.coord;
                global_atom->coord = atom.coord;
                // store the chain and ID of changed residues
                gate_residues.insert(gate_atom->chain + ":" + std::to_string(gate_atom->residue_id));
                global_residues.insert(global_atom->chain + ":" + std::to_string(global_atom->residue_id));
            }
        }
        // write the altered PDB and add a comment line with the altered residues
        std::string file_name = settings.output_name + gate.name() + ".pdb";
        write_PDB((settings.gate_dir / fs::path(file_name)).string(), gate.atoms);
        add_comment((settings.gate_dir / fs::path(file_name)).string(), gate_residues);
        // output the total run time needed for this gate
        print(2, start_gate, "=> total gate runtime:");
    }
    // write the altered PDB with all opened gates and add a comment line with the altered residues
    write_PDB(settings.gate_merged_pdb.string(), global_atoms);
    add_comment(settings.gate_merged_pdb.string(), global_residues);
}