#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "atom.h"

// extract all valid atom entries from a given PDB input file, allow to skip H-atoms and to load alternative locations
// of atoms in addition to the primary location
void parse_PDB(const std::string &filename, std::vector<std::shared_ptr<Atom>> &atoms, const bool skip_H,
               const bool keep_alternative) {
    std::ifstream file;
    file.open(filename);

    std::string line;
    size_t id = 0;

    while (getline(file, line)) {
        // record type: ATOM (protein atoms) or HETATM (non-protein atoms, including solvent molecules)
        // only consider protein atoms
        if (remove_spaces(line.substr(0, 6)) != "ATOM") continue;
        // it can happen that the z-coordinate column is not given in full
        // => make sure that the z-column is at least 1 long and determine its length dynamically
        if (line.size() < 47) continue;
        // create the atom and a shared pointer to it
        std::shared_ptr<Atom> atom = std::make_shared<Atom>(id, line);
        // make sure the atom is valid before adding it to the atom vector
        if (atom->type == INVALID_ATOM || atom->residue_type == INVALID_RESIDUE) continue;
        // skip H atoms if specified
        if (skip_H && atom->type == H) continue;
        // skip alternative locations of atoms if keep-alternatives is not enabled
        if (!keep_alternative && !atom->alternate_location.empty() && atom->alternate_location != "A") continue;
        // add the atom to the list
        atoms.push_back(atom);
        id++;
    }
    file.close();
}

// format an index and coordinate in pseudo PDB format, can be used to output pore/cavity boxes for visualisation
std::string pseudo_pdb(const size_t id, const Vec<double> &coord) {
    std::stringstream ss;
    ss << "ATOM  ";
    // due to width limitations of this column, higher IDs can only be represented with asterisks
    ss << std::setw(5) << std::setfill(' ') << std::right;
    (id <= 99999) ? ss << id : ss << "*****";
    ss << "  " << std::setw(3) << std::setfill(' ') << std::left << "C";
    ss << std::setw(1) << std::setfill(' ') << std::right << "";
    ss << std::setw(3) << std::setfill(' ') << std::right << "AAA";
    ss << std::setw(2) << std::setfill(' ') << std::right << "A";
    ss << std::setw(4) << std::setfill(' ') << std::right << "1";
    ss << std::setw(1) << std::setfill(' ') << std::right << "";
    ss << "   ";
    ss << std::setw(8) << std::setprecision(3) << std::setfill(' ') << std::right << std::fixed << coord.x;
    ss << std::setw(8) << std::setprecision(3) << std::setfill(' ') << std::right << std::fixed << coord.y;
    ss << std::setw(8) << std::setprecision(3) << std::setfill(' ') << std::right << std::fixed << coord.z;
    return ss.str();
}

// format information of a single atom in PDB format
std::string PDB_atom_entry(const std::shared_ptr<Atom> &atom) {
    if (std::isinf(atom->coord.x) || std::isinf(atom->coord.y) || std::isinf(atom->coord.z)) return "";
    std::stringstream ss;
    ss << "ATOM  ";
    // due to width limitations of this column, higher IDs can only be represented with asterisks
    if (atom->PDB_id <= 99999) {
        ss << std::setw(5) << std::setfill(' ') << std::right << atom->PDB_id;
    } else {
        ss << "*****";
    }
    ss << "  " << std::setw(3) << std::setfill(' ') << std::left << atom->full_type;
    ss << std::setw(1) << std::setfill(' ') << std::right << atom->alternate_location;
    ss << std::setw(3) << std::setfill(' ') << std::right << to_str(atom->residue_type);
    ss << std::setw(2) << std::setfill(' ') << std::right << atom->chain;
    // due to width limitations of this column, higher IDs can only be represented with asterisks
    if (atom->residue_id <= 9999) {
        ss << std::setw(4) << std::setfill(' ') << std::right << atom->residue_id;
    } else {
        ss << "****";
    }
    ss << std::setw(1) << std::setfill(' ') << std::right << atom->insertion_code;
    ss << "   ";
    ss << std::setw(8) << std::setprecision(3) << std::setfill(' ') << std::right << std::fixed << atom->coord.x;
    ss << std::setw(8) << std::setprecision(3) << std::setfill(' ') << std::right << std::fixed << atom->coord.y;
    ss << std::setw(8) << std::setprecision(3) << std::setfill(' ') << std::right << std::fixed << atom->coord.z;
    ss << std::setw(6) << std::setprecision(2) << std::setfill(' ') << atom->atom_occupancy;
    ss << std::setw(6) << std::setprecision(2) << std::setfill(' ') << atom->temperature;
    ss << std::setw(12) << std::setfill(' ') << atom->element;
    ss << std::setw(2) << std::setfill(' ') << atom->charge;
    return ss.str();
}

// output the information of all atoms in a PDB file
void write_PDB(const std::string &file_path, const std::vector<std::shared_ptr<Atom>> &atoms) {
    std::ofstream file;
    file.open(file_path);
    for (const std::shared_ptr<Atom> &atom: atoms) { file << PDB_atom_entry(atom) << "\n"; }
    file.close();
}