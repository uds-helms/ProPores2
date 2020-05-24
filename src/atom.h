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

#ifndef PROPORES_ATOM_H
#define PROPORES_ATOM_H

#include <string>
#include <vector>
#include <memory>
#include "enums.h"
#include "basics.h"
#include "vector.h"

// ATOM PROCESSING
// parses the atom type column in a PDB line
std::pair<AtomType, std::string> parse_atom_type(const std::string &str, RecordType record);

// van der Waals radius of different atom types
double vdw_radius(AtomType type);

// atomic mass of different atom types
double atomic_mass(AtomType type);

struct Atom {
    // atom ID in the program and the PDB file
    size_t id;
    size_t PDB_id;
    // record type, i.e. ATOM or HETATM
    RecordType record;
    // atom type, i.e. H, C, S, O,...
    AtomType type;
    // full type, i.e. CA, CB
    std::string full_type;
    // side chain mobility can lead to several possible locations of an atom, this is the corresponding label
    std::string alternate_location;
    // corresponding residue information
    ResidueType residue_type;
    std::string residue_name;
    std::string chain;
    size_t residue_id;
    std::string insertion_code;
    // x, y, z coordinates in Angstrom
    Vec<double> original;
    Vec<double> coord;
    // fraction of proteins where the atom is in this position
    std::string atom_occupancy;
    // movement of the atom
    std::string temperature;
    // element symbol
    AtomType element;
    // charge
    std::string charge;
    // vdW-radius
    double vdw;
    // atomic mass
    double mass;

    Atom(const size_t atom_id, const std::string &pdb_line) :
            id(atom_id),
            record(to_record_type(pdb_substr(pdb_line, 0, 6))),
            PDB_id(size_t(pdb_substr(pdb_line, 6, 5, 100000))),
            alternate_location(pdb_substr(pdb_line, 16, 1)),
            residue_type(to_residue_type(pdb_substr(pdb_line, 17, 3))),
            residue_name(pdb_substr(pdb_line, 17, 3)),
            chain(pdb_substr(pdb_line, 21, 1)),
            residue_id(size_t(pdb_substr(pdb_line, 22, 4, 10000))),
            insertion_code(pdb_substr(pdb_line, 26, 1)),
            coord(Vec<double>(stod(pdb_substr(pdb_line, 30, 8)),
                              stod(pdb_substr(pdb_line, 38, 8)),
                              stod(pdb_substr(pdb_line, 46, 8)))),
            atom_occupancy(pdb_substr(pdb_line, 54, 6)),
            temperature(pdb_substr(pdb_line, 60, 6)),
            element(to_atom_type(pdb_substr(pdb_line, 76, 2))),
            charge(pdb_substr(pdb_line, 78, 2)) {
        // extract the long and short form of the atom type
        if (element != INVALID_ATOM) {
            type = element;
            full_type = pdb_substr(pdb_line, 12, 4);
        } else {
            std::pair<AtomType, std::string> pair = parse_atom_type(pdb_line.substr(12, 4), record);
            type = pair.first;
            full_type = pair.second;
            element = type;
        }
        // compute the van der Waals radius and the atomic mass
        vdw = vdw_radius(type);
        mass = atomic_mass(type);
        // save the input coordinates in case the coordinates get altered and the original one's are still needed
        original = coord;
    }

    // KD-TREE
    // pointers to the left and right child nodes
    std::shared_ptr<Atom> KD_left, KD_right;
    // index of the splitting axis, i.e. 0 for the x-axis, 1 for the y-axis and 2 for the z-axis
    size_t KD_axis = 0;
    // true if the atom is part of a KD-tree, false otherwise
    bool in_tree = false;

    // checks
    [[nodiscard]] bool has_left() const { return KD_left != nullptr; }

    [[nodiscard]] bool has_right() const { return KD_right != nullptr; }

    [[nodiscard]] bool is_leaf() const { return in_tree && !has_left() && !has_right(); }
};

struct AtomPair {
    // pair ID
    size_t id;
    // the two atoms
    std::shared_ptr<Atom> atom_1;
    std::shared_ptr<Atom> atom_2;
    // vector connecting atom 1 and 2 and its direction
    Vec<double> vec;
    Vec<double> unit;
    // smaller van der Waals radius
    double r;

    AtomPair(const size_t i, std::shared_ptr<Atom> &a_1, std::shared_ptr<Atom> &a_2,
             const Vec<double> &vector, const Vec<double> &vector_unit) :
    // pair, atom 1 and atom 2 ID
            id(i),
            atom_1(a_1),
            atom_2(a_2),
            // vector connecting atom 1 and 2
            vec(vector),
            unit(vector_unit),
            r(std::min(a_1->vdw, a_2->vdw)) {}
};


// ALGORITHMS
// sort atoms by their ID in ascending order
void sort_atoms(std::vector<std::shared_ptr<Atom>> &atoms);

// recursively build a KD-tree of atom coordinates
std::shared_ptr<Atom> &KD_tree(int depth, std::vector<std::shared_ptr<Atom>>::iterator left,
                               std::vector<std::shared_ptr<Atom>>::iterator right);

// iteratively compute the distance to the closest van der Waals radius for an input coordinate (q)
std::tuple<double, std::shared_ptr<Atom>> KD_nearest_neighbour(const std::shared_ptr<Atom> &root, const Vec<double> &q);

#endif //PROPORES_ATOM_H
