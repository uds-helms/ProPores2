#ifndef PROPORES_READER_H
#define PROPORES_READER_H

#include <vector>
#include "atom.h"

// extract all valid atom entries from a given PDB input file, allow to skip H-atoms and to load alternative locations
// of atoms in addition to the primary location
void parse_PDB(const std::string &filename, std::vector<std::shared_ptr<Atom>> &atoms, bool skip_H,
               bool keep_alternative);

// format an index and coordinate in pseudo PDB format, can be used to output pore/cavity boxes for visualisation
std::string pseudo_pdb(size_t id, const Vec<double> &coord);

// format information of a single atom in PDB format
std::string PDB_atom_entry(const std::shared_ptr<Atom> &atom);

// output the information of all atoms in a PDB file
void write_PDB(const std::string &file_path, const std::vector<std::shared_ptr<Atom>> &atoms);

#endif //PROPORES_READER_H