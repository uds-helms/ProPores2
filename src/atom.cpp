#include <stack>
#include <string>
#include <iostream>
#include "atom.h"
#include "enums.h"
#include "vector.h"

// ATOM PROCESSING
// parses the atom type column in a PDB line
std::pair<AtomType, std::string> parse_atom_type(const std::string &str) {
    if (str.empty() || str.size() > 4) return std::make_pair(to_atom_type(""), "");
    for (const char &c: str) {
        if (c == 'H') return std::make_pair(to_atom_type("H"), str);
    }
    return std::make_pair(to_atom_type(str.substr(0, 1)), str);
}

// van der Waals radius of different atom types
double vdw_radius(const AtomType type) {
    if (type == N) return 1.55;
    if (type == C) return 1.7;
    if (type == O) return 1.52;
    if (type == H) return 1.2;
    if (type == S) return 1.8;
    // default
    return 0.0;
}

// atomic mass of different atom types
double atomic_mass(const AtomType type) {
    if (type == N) return 14.01;
    if (type == C) return 12.01;
    if (type == O) return 16.00;
    if (type == H) return 1.008;
    if (type == S) return 32.07;
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


