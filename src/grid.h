#ifndef PROPORES_GRID_H
#define PROPORES_GRID_H

#include <set>
#include <map>
#include <stack>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <iostream>
#include <exception>
#include <unordered_set>
#include <unordered_map>
#include "atom.h"
#include "gate.h"
#include "enums.h"
#include "vector.h"
#include "basics.h"
#include "reader.h"
#include "settings.h"

// represents a grid box in a (potential) pore or cavity
struct PoreBox {
    Vec<int> coord;
    BoxState state = IN_CLUSTER;
    // corresponding states and distance to the protein centre of mass
    double distance_to_centre = 0.0;
    // maximum sphere radius fitting at each box without touching the closest van der Waals radius
    double distance_to_closest = 0.0;

    // constructor
    explicit PoreBox(const Vec<int> &box_coordinates) : coord(box_coordinates) {}
    // equality (for testing)
    bool operator==(const Vec<int> &vec) const { return coord == vec;}
};

// a cluster of boxes that forms a pore or cavity
struct PoreCluster {
    size_t id;
    // 3D grid coordinates of boxes in the cluster, as well as a set of their 1D indices for duplicate elimination
    std::vector<PoreBox> boxes;
    std::unordered_set<size_t> indices;
    // residues directly surrounding the pore or cavity
    std::set<std::tuple<std::string, size_t, ResidueType>> lining_residues;
    // true if the cluster is a pore, false if it's a cavity
    bool pore = false;
    // grid box length
    double box_length;

    // constructor
    PoreCluster(const size_t identifier, const double length) :
            id(identifier),
            box_length(length) {}

    // add a box to the cluster
    void add_box(const size_t index, const Vec<int> &box) {
        size_t previous_size = indices.size();
        indices.insert(index);
        // make sure to only add box coordinates once
        if (indices.size() > previous_size) boxes.emplace_back(box);
    }

    // extract residue information from the given atom and add it to the lining residue set
    void add_lining_residue(const std::shared_ptr<Atom> &atom) {
        lining_residues.emplace(std::make_tuple(atom->chain, atom->residue_id, atom->residue_type));
    }

    // check if a 1D index is part of the cluster
    [[nodiscard]] bool has_index(const size_t index) const { return indices.find(index) != indices.end(); }

    // number of boxes in the pore or cavity cluster
    [[nodiscard]] size_t size() const { return boxes.size(); }

    [[nodiscard]] bool empty() const { return boxes.empty(); }

    // return the output name of the cluster
    [[nodiscard]] std::string name() const { return std::to_string(id) + (pore ? "_pore" : "_cavity"); }
};

struct Grid {
    // CONSTANTS
    // grid box dimensions
    const double box_length;
    const double box_length_sq;
    const double box_radius;
    const double box_volume;

    // DIMENSIONS
    // minimum and maximum box coordinates in each dimension (can be negative, minimum is included, maximum not)
    Vec<int> min;
    Vec<int> max;
    // number of boxes in each dimension of the grid
    size_t width = 0;
    size_t height = 0;
    size_t depth = 0;
    // number of boxes in the grid
    size_t boxes = 0;


    // search bubbles of relative coordinates for different radii
    RelativeSearchSpaces search_spaces;
    // direct neighbour boxes
    std::vector<Vec<int>> neighbours = {Vec<int>(-1, 0, 0), Vec<int>(1, 0, 0), Vec<int>(0, -1, 0),
                                        Vec<int>(0, 1, 0), Vec<int>(0, 0, -1), Vec<int>(0, 0, 1)};

    // state of each grid box (e.g. occupied, part of the background, potential pore,...)
    std::vector<BoxState> states;

    explicit Grid(const double length) :
            box_length(length),
            box_length_sq(pow(box_length, 2)),
            box_radius(box_length / 2),
            box_volume(pow(box_length, 3)) {}

    // sets and initialises grid dimensions, number of grid boxes and box states
    void set_dimensions() {
        // set the grid dimensions
        width = abs(max.x - min.x);
        height = abs(max.y - min.y);
        depth = abs(max.z - min.z);
        // compute the number of grid boxes
        boxes = height * width * depth;
        // initialise the grid box states
        states = std::vector<BoxState>(boxes);
    }

    size_t size() { return boxes; }

    // COORDINATE CONVERSIONS and GRID ACCESS
    // get box coordinates confined to the grid
    [[nodiscard]] Vec<int> box(const Vec<double> &vec) const {
        return Vec<int>(int(floor(vec.x / box_length)), int(floor(vec.y / box_length)), int(floor(vec.z / box_length)));
    }

    // checks if a box is within the bounds of the grid
    bool valid(const Vec<int> &vec) {
        if (vec.x < min.x || vec.x >= max.x) return false;
        if (vec.y < min.y || vec.y >= max.y) return false;
        return vec.z >= min.z && vec.z < max.z;
    }

    bool valid(const size_t i) { return i >= 0 && i < boxes; }

    // get box centre coordinates from grid coordinates
    [[nodiscard]] Vec<double> centre(const int x, const int y, const int z) const {
        return Vec<double>(x * box_length + box_radius, y * box_length + box_radius, z * box_length + box_radius);
    }

    [[nodiscard]] Vec<double> centre(const Vec<int> &vec) const { return centre(vec.x, vec.y, vec.z); }

    // convert 3D coordinates to 1D vector indices
    [[nodiscard]] size_t index(const int x, const int y, const int z) const {
        Vec<int> v = Vec<int>(x, y, z) - min;
        return v.x + width * (v.y + height * v.z);
    }

    [[nodiscard]] size_t index(const Vec<double> &vec) const {
        return index(int(floor(vec.x / box_length)), int(floor(vec.y / box_length)), int(floor(vec.z / box_length)));
    }

    [[nodiscard]] size_t index(const Vec<int> &vec) const { return index(vec.x, vec.y, vec.z); }

    // get the state/label of a box
    BoxState &at(const size_t i) { return states[i]; }

    BoxState &at(const int x, const int y, const int z) { return states[index(x, y, z)]; }

    BoxState &at(const Vec<double> &vec) { return states[index(vec)]; }

    BoxState &at(const Vec<int> &vec) { return states[index(vec)]; }
};

struct ProteinGrid : Grid {
    // CONSTANTS
    // grid box dimensions
    const double box_radius_sqrt3;
    // maximum van der Waals radius
    const double max_vdw;
    // radius of the solvent surrounding the protein
    const double solvent_radius;
    const double solvent_radius_sq;
    // margin on each side of the grid
    double margin;

    // centre of mass
    Vec<double> centre_of_mass;
    double mass = 0;

    // DATA
    // atoms in the input molecule
    std::vector<std::shared_ptr<Atom>> atoms;
    // atom pairs with empty space between them
    std::vector<std::shared_ptr<AtomPair>> pairs;
    // maximum sphere radius fitting at each box without touching the closest van der Waals radius
    std::vector<double> ball_radii;
    // which boxes are occupied by which atoms and van der Waals radii
    std::vector<std::unordered_set<size_t>> occupied;
    // which boxes are part of cylinders between atom pairs
    std::vector<std::unordered_set<size_t>> trace;
    // clusters of empty grid boxes that potentially constitute pores or cavities
    std::vector<PoreCluster> clusters;

    // HELPERS
    size_t next_pair_id = 0;
    size_t next_cluster_id = 0;
    // KD-tree for efficiently identifying the closest distance to van der Waals radii
    bool KD_constructed = false;
    std::shared_ptr<Atom> KD_root;

    // CONSTRUCTORS
    // constructor interface for the actual program with user command line input
    explicit ProteinGrid(const Settings &settings) :
            ProteinGrid(settings.pdb_path.string(), settings.box_length, settings.solvent_radius, settings.skip_H,
                        settings.keep_alternative) {}

    // detailed constructor for testing
    ProteinGrid(const std::string &pdb_path, const double length, const double solvent,
                const bool skip_H, const bool keep_alternative) :
            Grid(length),
            box_radius_sqrt3(box_radius * sqrt(3)),
            solvent_radius(solvent),
            solvent_radius_sq(pow(solvent, 2)),
            max_vdw(vdw_radius(S)) {
        // INPUT PARSING
        // parse atom information from the input file
        parse_PDB(pdb_path, atoms, skip_H, keep_alternative);
        // make sure that there are valid atoms in the input file
        if (atoms.empty()) throw std::runtime_error("There are no valid atoms in the input PDB file.");

        // COORDINATE SHIFTING
        // compute the grid margin from the solvent radius and the maximum vdW-radius
        margin = max_vdw + solvent_radius * 2;

        // initialise the minimum and maximum with the coordinates of the first atom
        Vec<double> atom_min(atoms[0]->coord.x, atoms[0]->coord.y, atoms[0]->coord.z);
        Vec<double> atom_max(atoms[0]->coord.x, atoms[0]->coord.y, atoms[0]->coord.z);

        // update the minimum and maximum coordinates
        for (const std::shared_ptr<Atom> &atom: atoms) {
            if (atom->coord.x < atom_min.x) atom_min.x = atom->coord.x;
            if (atom->coord.y < atom_min.y) atom_min.y = atom->coord.y;
            if (atom->coord.z < atom_min.z) atom_min.z = atom->coord.z;

            if (atom->coord.x > atom_max.x) atom_max.x = atom->coord.x;
            if (atom->coord.y > atom_max.y) atom_max.y = atom->coord.y;
            if (atom->coord.z > atom_max.z) atom_max.z = atom->coord.z;
        }

        // GRID CONSTRUCTION
        // compute the maximum coordinates including the margin on each side
        Vec<double> max_coord = atom_max + margin;
        Vec<double> min_coord = atom_min - margin;
        // compute the number of grid boxes in each dimension
        max = Vec<int>(int(ceil(max_coord.x / box_length)),
                       int(ceil(max_coord.y / box_length)),
                       int(ceil(max_coord.z / box_length)));

        min = Vec<int>(int(floor(min_coord.x / box_length)),
                       int(floor(min_coord.y / box_length)),
                       int(floor(min_coord.z / box_length)));

        // initialise the grid dimensions based on the minimum and maximum
        set_dimensions();

        // initialise the managers
        occupied = std::vector<std::unordered_set<size_t>>(boxes);

        // FILL THE GRID (and compute CENTRE OF MASS)
        for (std::shared_ptr<Atom> &atom: atoms) {
            // add the atom to the corresponding grid box
            Vec<int> grid_box = box(atom->coord);
            size_t i = index(grid_box);
            occupied[i].insert(atom->id);
            states[i] = OCCUPIED;

            // also include boxes with centre points inside the vdW-radius of the atom
            double vdw_radius_squared = pow(atom->vdw, 2);
            // number of boxes to check around the atom box
            int check_radius = int(ceil(atom->vdw / box_length));

            for (int x = grid_box.x - check_radius; x <= grid_box.x + check_radius; x++) {
                for (int y = grid_box.y - check_radius; y <= grid_box.y + check_radius; y++) {
                    for (int z = grid_box.z - check_radius; z <= grid_box.z + check_radius; z++) {
                        Vec<int> search_box = Vec<int>(x, y, z);
                        // skip search boxes that are outside the grid bounds
                        if (!valid(search_box)) continue;
                        // check if the center of the box is within the vdW sphere
                        // (the vdW-radius is squared because we use the squared distance from the atom in order not to
                        // compute the square root)
                        if (distance_squared(centre(search_box), atom->coord) <= vdw_radius_squared) {
                            i = index(x, y, z);
                            occupied[i].insert(atom->id);
                            states[i] = OCCUPIED;
                        }
                    }
                }
            }
            // accumulate the atom coordinates scaled by the respective atomic mass, as well as the total atomic mass
            // of non-H atoms
            if (atom->type != H) {
                centre_of_mass += atom->coord * atom->mass;
                mass += atom->mass;
            }
        }
        // compute the centre of mass
        if (mass > 0) centre_of_mass /= mass;
    }

    // KD-TREE NEAREST NEIGHBOUR SEARCH
    // construct a KD-tree of the atoms in the grid
    void build_KD_tree() {
        // recursively build the tree
        KD_root = KD_tree(0, atoms.begin(), atoms.end() - 1);
        // bring the atoms back into the original order, meaning index 0 is atom with ID 0
        sort_atoms(atoms);
        KD_constructed = true;
    }

    // use the KD-tree to identify the closest distance to atom van der Waals radii
    std::tuple<double, std::shared_ptr<Atom>> closest(const Vec<double> &query) {
        if (!KD_constructed) build_KD_tree();
        return KD_nearest_neighbour(KD_root, query);
    }

    std::tuple<double, std::shared_ptr<Atom>> closest(const int x, const int y, const int z) {
        return closest(centre(x, y, z));
    }

    // PORE-ID
    // create and add an atom pair to the list
    void add_atom_pair(std::shared_ptr<Atom> &atom_1, std::shared_ptr<Atom> &atom_2, const Vec<double> &v,
                       const Vec<double> &v_unit) {
        pairs.emplace_back(new AtomPair{next_pair_id, atom_1, atom_2, v, v_unit});
        next_pair_id++;
    }

    // create a new cluster and add it to the cluster vector
    PoreCluster &new_cluster() {
        clusters.emplace_back(next_cluster_id, box_length);
        next_cluster_id++;
        return clusters[clusters.size() - 1];
    }

    // filters and sorts pores/cavities
    void filter_pores(const double volume_threshold, const RemoveTag remove_tag) {
        // extract all valid pores/cavities
        std::vector<PoreCluster> valid_pores;
        for (PoreCluster &cluster: clusters) {
            // if specified, remove cavities
            if (remove_tag == REMOVE_CAVITIES && !cluster.pore) continue;
            // if specified, remove pores
            if (remove_tag == REMOVE_PORES && cluster.pore) continue;
            // remove pores or cavities that are too small
            if (cluster.size() * box_volume < volume_threshold) continue;
            valid_pores.push_back(std::move(cluster));
        }
        // sort the remaining pores/cavities by size
        std::sort(valid_pores.begin(), valid_pores.end(),
                  [](PoreCluster &lhs, PoreCluster &rhs) { return lhs.size() > rhs.size(); });
        // adjust the IDs of the remaining pores/cavities
        for (size_t i = 0; i < valid_pores.size(); i++) { valid_pores[i].id = i; }
        next_cluster_id = valid_pores.size();
        // only keep the valid pores/cavities
        clusters = valid_pores;
    }
};

struct PoreGrid : Grid {
    // DATA
    // maximum sphere radius fitting at each box without touching the closest van der Waals radius
    std::vector<double> radii;
    // distance of the box to the protein centre of mass
    std::vector<double> distance_to_centre;
    // starting points for axes determination with Dijkstra
    std::vector<Vec<int>> starting_points;
    // 3D coordinates
    std::vector<Vec<int>> coordinates;

    // HELP
    // small number to avoid divisions by 0
    const double perturb;
    // true if the pore has at least one large enough surface patch exposed to the solvent, false otherwise
    bool has_surface_patch = false;
    // number of boxes belonging to the pore/cavity
    size_t n_pore_boxes;

    // constructor
    PoreGrid(const PoreCluster &cluster, const double perturb) :
            Grid(cluster.box_length),
            perturb(perturb) {
        // INPUT VALIDATION
        // make sure there is at least one input box
        if (cluster.boxes.empty()) throw std::invalid_argument("PoreGrid: no pore boxes given.");

        // GRID CONSTRUCTION
        // number of pore boxes in the grid
        n_pore_boxes = cluster.size();
        // determine the minimum and maximum box coordinates; since the input coordinates are already grid coordinates,
        // the minimum and maximum can be compute directly
        min = cluster.boxes[0].coord;
        max = cluster.boxes[0].coord;

        for (const PoreBox &box: cluster.boxes) {
            if (box.coord.x < min.x) min.x = box.coord.x;
            if (box.coord.y < min.y) min.y = box.coord.y;
            if (box.coord.z < min.z) min.z = box.coord.z;

            if (box.coord.x > max.x) max.x = box.coord.x;
            if (box.coord.y > max.y) max.y = box.coord.y;
            if (box.coord.z > max.z) max.z = box.coord.z;
        }

        // max is exclusive, thus it has to be incremented by 1
        max += 1;
        // compute the grid dimensions based on the minimum and maximum
        set_dimensions();
        // initialise
        radii = std::vector<double>(boxes);
        distance_to_centre = std::vector<double>(boxes);
        coordinates = std::vector<Vec<int>>(boxes);

        // add the input information to the 1D representation
        for (const PoreBox &box: cluster.boxes) {
            size_t idx = index(box.coord);
            radii[idx] = box.distance_to_closest;
            distance_to_centre[idx] = box.distance_to_centre;
            states[idx] = box.state;
        }
        // store the 3D coordinates in the 1D representation
        for (int x = min.x; x < max.x; x++) {
            for (int y = min.y; y < max.y; y++) {
                for (int z = min.z; z < max.z; z++) {
                    coordinates[index(x, y, z)] = Vec<int>(x, y, z);
                }
            }
        }
    }

    // Dijkstra score
    double individual_score(const size_t i) {
        // if the pore/cavity has at least two starting points at surface patches or no surface patch at all,
        // then the score is 1 divided by the maximum radius fitting into the box
        if (starting_points.size() > 1 || !has_surface_patch) return 1 / (radii[i] + perturb);
        // otherwise if the pore/cavity currently only has one starting point from a surface patch, also
        // consider the distance of the box to the protein centre of mass
        return (distance_to_centre[i] + perturb) / (radii[i] + perturb);
    }
};

#endif //PROPORES_GRID_H
