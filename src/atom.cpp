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
#include <cctype>
#include <string>
#include <iostream>
#include "atom.h"
#include "enums.h"
#include "vector.h"

// ATOM PROCESSING
// parses the atom type column in a PDB line
std::pair<AtomType, std::string> parse_atom_type(const std::string &input_str, const RecordType record) {
    // get a version of the atom type column without surrounding spaces
    std::string str = remove_spaces(input_str);

    // check if the atom type column contains 1-4 characters
    if (str.empty() || str.size() > 4) return std::make_pair(INVALID_ATOM, "");

    // if there is only one character given, it must belong to the atom type
    if (str.size() == 1) return std::make_pair(to_atom_type(str), str);

    // for hydrogen atoms, the numbering may be in front of the atom abbreviation 'H'
    if (std::isdigit(str.at(0))) {
        for (const char &c: str) {
            if (c == 'H') return std::make_pair(H, str);
            if (!std::isdigit(c)) break;
        }
        return std::make_pair(INVALID_ATOM, "");
    }

    // the atom type column consists of two columns of two chars each, where the first sub-column contains the type
    // and the second sub-column contains numbering
    // there can be cases where the atom type is only 1 character long and the numbering is 3 characters long,
    // e.g. HG11, which is the first H-atom attached to the gamma C-atom CG1
    std::string left = remove_spaces(input_str.substr(0, 2));
    std::string right = remove_spaces(input_str.substr(2, 2));

    // if the left sub-column is empty, then there is an alignment problem
    if (left.empty()) {
        left = right;
        right = "";
    }

    // if there is only one char in the left sub-column, it must belong to the a one-letter atom type
    if (left.size() == 1) return std::make_pair(to_atom_type(left), str);

    // if there are two characters in the left sub-column, they can either belong both to the atom type or the right one
    // could belong to the numbering (see example above)
    AtomType left_type = to_atom_type(left);
    AtomType leftmost_type = to_atom_type(left.at(0));

    // if the entire left column cannot be converted to a valid atom type, then only the first char can belong to
    // the atom type
    if (left_type == INVALID_ATOM) return std::make_pair(leftmost_type, str);
    // if the entire left column can be converted to a valid atom type but only the left most character cannot,
    // then the entire left sub-column must be the atom type
    if (leftmost_type == INVALID_ATOM) return std::make_pair(left_type, str);

    // this is a case where both the entire left sub-column with two characters or only the leftmost character can
    // form valid atom type
    // the second character might belong to a 3-character long numbering or to a misaligned 2-character numbering

    // ATOM records contain atoms of standard amino acids (C, H, N, O, S), if that is not the case, there is an error
    // (the hydrogen case with numbering at the beginning is taken care of above)
    if (record == ATOM) {
        if (leftmost_type == C || leftmost_type == H || leftmost_type == N || leftmost_type == O || leftmost_type == S)
            return std::make_pair(leftmost_type, str);

        return std::make_pair(INVALID_ATOM, "");
    }

    // the alternative are HETATM records where the situation is not as clear cut
    // OXT is the carboxy oxygen atom in some residues
    if (str.size() >= 3 && str.substr(0, 3) == "OXT") return std::make_pair(O, str);

    return std::make_pair(left_type, str);
}

// van der Waals radius of different atom types in Ångström
double vdw_radius(const AtomType type) {
    // sources:
    // 1.  https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
    //     which uses data from:
    //     - A. Bondi (1964). "van der Waals Volumes and Radii". J. Phys. Chem. 68: 441. doi:10.1021/j100785a001
    //     - M. Mantina; A.C. Chamberlin; R. Valero; C.J. Cramer; D.G. Truhlar (2009). "Consistent van der
    //       Waals Radii for the Whole Main Group". J. Phys. Chem. A. 113 (19): 5806–12
    // 2.  Batsanov, S.S., Van der Waals Radii Evaluated from
    //     Structural Parameters of Metals, Zh. Fiz. Khim., 2000,
    //     vol. 74, no. 7, pp. 1273–1276
    // 3.  https://www.internetchemie.info/chemische-elemente/
    // 4.  Chemical Elements Pocket Guide: Detailed Summary of the Periodic Table
    //     By Coventry House Publishing

    // protein atoms
    if (type == C) return 1.7;
    if (type == H) return 1.2;
    if (type == N) return 1.55;
    if (type == O) return 1.52;
    if (type == S) return 1.8;
    // hetero atoms
    if (type == Ac) return 2.0;
    if (type == Ag) return 1.72;
    if (type == Al) return 1.84;
    if (type == Am) return 2.28;
    if (type == Ar) return 1.88;
    if (type == As) return 1.85;
    if (type == At) return 2.02;
    if (type == Au) return 1.66;
    if (type == B) return 1.92;
    if (type == Ba) return 2.68;
    if (type == Be) return 1.53;
    if (type == Bh) return 2.0;
    if (type == Bi) return 2.07;
    if (type == Bk) return 2.0;
    if (type == Br) return 1.85;
    if (type == Ca) return 2.31;
    if (type == Cd) return 1.58;
    if (type == Ce) return 2.0;
    if (type == Cf) return 2.0;
    if (type == Cl) return 1.75;
    if (type == Cm) return 2.0;
    if (type == Co) return 1.95;
    if (type == Cr) return 1.97;
    if (type == Cs) return 3.43;
    if (type == Cu) return 1.4;
    if (type == Db) return 2.0;
    if (type == Ds) return 1.86;
    if (type == Dy) return 1.55;
    if (type == Er) return 2.02;
    if (type == Es) return 2.0;
    if (type == Eu) return 2.0;
    if (type == F) return 1.47;
    if (type == Fe) return 1.96;
    if (type == Fm) return 2.0;
    if (type == Fr) return 3.48;
    if (type == Ga) return 1.87;
    if (type == Gd) return 1.72;
    if (type == Ge) return 2.11;
    if (type == He) return 1.4;
    if (type == Hf) return 2.24;
    if (type == Hg) return 1.55;
    if (type == Ho) return 1.96;
    if (type == Hs) return 2.0;
    if (type == I) return 1.98;
    if (type == In) return 1.93;
    if (type == Ir) return 2.03;
    if (type == K) return 2.75;
    if (type == Kr) return 2.02;
    if (type == La) return 2.51;
    if (type == Li) return 1.82;
    if (type == Lr) return 2.0;
    if (type == Lu) return 2.0;
    if (type == Md) return 2.0;
    if (type == Mg) return 1.73;
    if (type == Mn) return 1.96;
    if (type == Mo) return 2.06;
    if (type == Mt) return 2.0;
    if (type == Na) return 2.27;
    if (type == Nb) return 2.13;
    if (type == Nd) return 2.0;
    if (type == Ne) return 1.54;
    if (type == Ni) return 1.63;
    if (type == No) return 2.0;
    if (type == Np) return 2.0;
    if (type == Os) return 2.03;
    if (type == P) return 1.8;
    if (type == Pa) return 2.0;
    if (type == Pb) return 2.02;
    if (type == Pd) return 1.63;
    if (type == Pm) return 2.0;
    if (type == Po) return 1.97;
    if (type == Pr) return 2.0;
    if (type == Pt) return 1.75;
    if (type == Pu) return 2.0;
    if (type == Ra) return 2.83;
    if (type == Rb) return 3.03;
    if (type == Re) return 2.05;
    if (type == Rf) return 2.0;
    if (type == Rh) return 2.02;
    if (type == Rn) return 2.2;
    if (type == Ru) return 2.02;
    if (type == Sb) return 2.06;
    if (type == Sc) return 2.11;
    if (type == Se) return 1.9;
    if (type == Sg) return 2.0;
    if (type == Si) return 2.1;
    if (type == Sm) return 2.0;
    if (type == Sn) return 2.17;
    if (type == Sr) return 2.49;
    if (type == Ta) return 2.13;
    if (type == Tb) return 1.66;
    if (type == Tc) return 2.04;
    if (type == Te) return 2.06;
    if (type == Th) return 2.43;
    if (type == Ti) return 2.14;
    if (type == Tl) return 1.96;
    if (type == Tm) return 2.0;
    if (type == Ts) return 1.38;
    if (type == U) return 1.86;
    if (type == V) return 2.03;
    if (type == W) return 2.07;
    if (type == Xe) return 2.16;
    if (type == Y) return 2.45;
    if (type == Yb) return 2.0;
    if (type == Zn) return 1.39;
    if (type == Zr) return 2.25;
    // default
    return 0.0;
}

// atomic mass of different atom types
double atomic_mass(const AtomType type) {
    // protein atoms
    if (type == C) return 12.01;
    if (type == H) return 1.008;
    if (type == O) return 16.00;
    if (type == N) return 14.01;
    if (type == S) return 32.07;
    // hetero atoms
    if (type == Ac) return 227;
    if (type == Ag) return 107.87;
    if (type == Al) return 26.98;
    if (type == Am) return 243;
    if (type == Ar) return 39.95;
    if (type == As) return 74.92;
    if (type == At) return 210;
    if (type == Au) return 196.97;
    if (type == B) return 10.81;
    if (type == Ba) return 137.33;
    if (type == Be) return 9.01;
    if (type == Bh) return 264;
    if (type == Bi) return 208.98;
    if (type == Bk) return 247;
    if (type == Br) return 79.9;
    if (type == Ca) return 40.08;
    if (type == Cd) return 112.41;
    if (type == Ce) return 140.12;
    if (type == Cf) return 251;
    if (type == Cl) return 35.45;
    if (type == Cm) return 247;
    if (type == Cn) return 285.18;
    if (type == Co) return 58.93;
    if (type == Cr) return 52.0;
    if (type == Cs) return 132.91;
    if (type == Cu) return 63.55;
    if (type == Db) return 262;
    if (type == Ds) return 271;
    if (type == Dy) return 162.5;
    if (type == Er) return 167.26;
    if (type == Es) return 252;
    if (type == Eu) return 151.97;
    if (type == F) return 19;
    if (type == Fe) return 55.85;
    if (type == Fm) return 257;
    if (type == Fr) return 223;
    if (type == Ga) return 69.72;
    if (type == Gd) return 157.25;
    if (type == Ge) return 72.64;
    if (type == He) return 4.003;
    if (type == Hf) return 178.49;
    if (type == Hg) return 200.59;
    if (type == Ho) return 164.93;
    if (type == Hs) return 277;
    if (type == I) return 126.9;
    if (type == In) return 114.82;
    if (type == Ir) return 192.22;
    if (type == K) return 39.1;
    if (type == Kr) return 83.8;
    if (type == La) return 138.91;
    if (type == Li) return 6.94;
    if (type == Lr) return 262;
    if (type == Lu) return 174.97;
    if (type == Lv) return 293.21;
    if (type == Mc) return 290.2;
    if (type == Md) return 258;
    if (type == Mg) return 24.31;
    if (type == Mn) return 54.94;
    if (type == Mo) return 95.94;
    if (type == Mt) return 268;
    if (type == Na) return 22.99;
    if (type == Nb) return 92.91;
    if (type == Nd) return 144.24;
    if (type == Ne) return 20.18;
    if (type == Nh) return 286.18;
    if (type == Ni) return 58.69;
    if (type == No) return 259;
    if (type == Np) return 237;
    if (type == Og) return 294.21;
    if (type == Os) return 190.23;
    if (type == P) return 30.97;
    if (type == Pa) return 231.04;
    if (type == Pb) return 207.2;
    if (type == Pd) return 106.42;
    if (type == Pm) return 145;
    if (type == Po) return 210;
    if (type == Pr) return 140.91;
    if (type == Pt) return 195.08;
    if (type == Pu) return 244;
    if (type == Ra) return 226;
    if (type == Rb) return 85.47;
    if (type == Re) return 186.21;
    if (type == Rf) return 261;
    if (type == Rg) return 272;
    if (type == Rh) return 102.91;
    if (type == Rn) return 220;
    if (type == Ru) return 101.07;
    if (type == Sb) return 121.76;
    if (type == Sc) return 44.96;
    if (type == Se) return 78.96;
    if (type == Sg) return 266;
    if (type == Si) return 28.09;
    if (type == Sm) return 150.36;
    if (type == Sn) return 118.71;
    if (type == Sr) return 87.62;
    if (type == Ta) return 180.95;
    if (type == Tb) return 158.93;
    if (type == Tc) return 98;
    if (type == Te) return 127.6;
    if (type == Th) return 232.04;
    if (type == Ti) return 47.87;
    if (type == Tl) return 204.38;
    if (type == Tm) return 168.93;
    if (type == Ts) return 294.211;
    if (type == U) return 238.03;
    if (type == V) return 50.94;
    if (type == W) return 183.84;
    if (type == Xe) return 131.29;
    if (type == Y) return 88.91;
    if (type == Yb) return 173.04;
    if (type == Zn) return 65.41;
    if (type == Zr) return 91.22;
    // default
    return 0.0;
}

// ALGORITHMS
// sort atoms by their ID in ascending order
void sort_atoms(std::vector<std::shared_ptr<Atom>> &atoms) {
    sort(atoms.begin(), atoms.end(),
         [](const std::shared_ptr<Atom> &a_1, const std::shared_ptr<Atom> &a_2) { return a_1->id < a_2->id; });
}

// recursively build a KD-tree of atom coordinates
std::shared_ptr<Atom> &KD_tree(const int depth, std::vector<std::shared_ptr<Atom>>::iterator left,
                               std::vector<std::shared_ptr<Atom>>::iterator right) {
    // determine the axis at this depth
    int const axis = depth % 3;
    // sort the atoms in the current range by the coordinates on the current axis
    // it has to be right + 1 since in sort the end is not included, but right is included
    sort(left, right + 1,
         [&axis](const std::shared_ptr<Atom> &a_1, const std::shared_ptr<Atom> &a_2) {
             return a_1->coord[axis] < a_2->coord[axis];
         });
    // determine the median atom
    auto median = left + int(floor((right - left) / 2));
    std::shared_ptr<Atom> &atom = (*median);
    // store the axis and depth of the atom
    atom->KD_axis = axis;
    atom->in_tree = true;

    // recursively add the left child
    if (left < median) atom->KD_left = KD_tree(depth + 1, left, median - 1);
    // recursively add the right child
    if (right > median) atom->KD_right = KD_tree(depth + 1, median + 1, right);

    return atom;
}

// iteratively compute the distance to the closest van der Waals radius for an input coordinate (q)
std::tuple<double, std::shared_ptr<Atom>>
KD_nearest_neighbour(const std::shared_ptr<Atom> &root, const Vec<double> &q) {
    double closest_dist = INFINITY;
    std::shared_ptr<Atom> closest = root;
    std::stack<std::shared_ptr<Atom>> stack = std::stack<std::shared_ptr<Atom>>();
    stack.push(root);

    while (!stack.empty()) {
        // DON'T use &node, otherwise it messes up all the pointers
        std::shared_ptr<Atom> node = stack.top();
        stack.pop();
        // update the currently closest node if the distance is shorter than the current minimum
        // need to incorporate the vdW-radius since grid boxes within the radii are considered to be part of the
        // atom
        double dist = (node->coord - q).length() - node->vdw;
        if (dist < closest_dist) {
            closest_dist = dist;
            closest = node;
        }
        // reached a leaf node and cannot recurse further
        if (node->is_leaf()) continue;

        // distance between the coordinate of the current node and the query point on the current axis
        double axis_distance = node->coord[node->KD_axis] - q[node->KD_axis];
        // closest_dist = dist(query, closest node) - vdW-radius of the closest node
        // => when the axis distance of the current node to the query is smaller than the closest distance + the biggest
        //    vdW-radius, then search in the left AND right sub-tree
        bool check_both = abs(axis_distance) < closest_dist + vdw_radius(S);
        // also check in the respective sub-tree when the query point on the current axis is smaller (left) or
        // bigger (right) than the current splitting value
        bool check_left = check_both || axis_distance >= 0;
        bool check_right = check_both || axis_distance <= 0;

        // when the query value on this axis is bigger than the current node's splitting value, then the query is
        // more likely to be located in the right sub-tree => search there first
        if (axis_distance >= 0) {
            // due to the iterative approach the node to be checked first has to be pushed last
            if (node->has_right() && check_right) stack.push(node->KD_right);
            if (node->has_left() && check_left) stack.push(node->KD_left);
            // otherwise search first in the right-subtree
        } else {
            if (node->has_left() && check_left) stack.push(node->KD_left);
            if (node->has_right() && check_right) stack.push(node->KD_right);
        }
    }
    return std::make_tuple(closest_dist, closest);
}


