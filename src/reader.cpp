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

#include <vector>
#include <cctype>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "atom.h"
#include "reader.h"
#include "settings.h"

// extract all valid atom entries from a given PDB input file, allow to skip H-atoms and to load alternative locations
// of atoms in addition to the primary location
void parse_PDB(const std::string &pdb_path, const std::string &kept_path, const std::string &skipped_path,
               std::vector<std::shared_ptr<Atom>> &atoms, PDBReaderStats &stats,
               const bool skip_H,
               const bool keep_alternative,
               const bool skip_hetero_atoms,
               const bool skip_non_standard_amino_acids) {
    // open the input PDB file, the file for logging the processed and kept atom lines and the file for skipped atoms
    std::ifstream pdb_file(pdb_path);
    std::ofstream kept_file(kept_path);
    std::ofstream skipped_file(skipped_path);
    // for extracting the content of the input PDB file line by line
    std::string line;
    // next atom ID
    size_t id = 0;

    while (getline(pdb_file, line)) {
        // record type: ATOM (protein atoms) or HETATM (non-protein atoms, including solvent molecules)
        RecordType record_type = to_record_type(remove_spaces(pdb_substr(line, 0, 6)));
        // skip all non-atom entries
        if (record_type != ATOM && record_type != HETATOM) continue;
        // skip hetero atom entries if the corresponding flag is set
        if (skip_hetero_atoms && record_type == HETATOM) {
            skipped_file << line << std::endl;
            stats.skipped_hetero++;
            continue;
        }
        // it can happen that the z-coordinate column is not given in full
        // => make sure that the z-column is at least 1 long and determine its length dynamically
        if (line.size() < 47) {
            skipped_file << line << std::endl;
            stats.line_too_short++;
            continue;
        }
        // create the atom and a shared pointer to it
        std::shared_ptr<Atom> atom = std::make_shared<Atom>(id, line);
        // make sure the atom is valid before adding it to the atom vector
        if (atom->type == INVALID_ATOM) {
            skipped_file << line << std::endl;
            stats.atom_type_parse_error++;
            continue;
        }
        // make sure that ATOM records have a valid amino acid residue
        if (atom->record == ATOM && skip_non_standard_amino_acids && atom->residue_type == INVALID_RESIDUE) {
            skipped_file << line << std::endl;
            stats.skipped_non_standard_amino_acids++;
            continue;
        }
        // make sure that van der Waals radius and atomic mass information is available
        if (atom->vdw == 0.0) {
            skipped_file << line << std::endl;
            stats.no_vdw_radius++;
            continue;
        }
        if (atom->mass == 0.0) {
            skipped_file << line << std::endl;
            stats.no_atomic_mass++;
            continue;
        }
        // skip H atoms if specified
        if (skip_H && atom->type == H) {
            skipped_file << line << std::endl;
            stats.skipped_H++;
            continue;
        }
        // skip alternative locations of atoms if keep-alternatives is not enabled
        if (!keep_alternative && !atom->alternate_location.empty() && atom->alternate_location != "A") {
            skipped_file << line << std::endl;
            stats.skipped_alternative++;
            continue;
        }
        // add the atom to the list and log it in the kept file
        kept_file << line << std::endl;
        stats.valid++;
        atoms.push_back(atom);
        id++;
    }
    pdb_file.close();
    kept_file.close();
    skipped_file.close();
}

// wrapper that allows calling parse_PDB with the program settings
void parse_PDB(const Settings &settings, std::vector<std::shared_ptr<Atom>> &atoms, PDBReaderStats &stats) {
    parse_PDB(settings.pdb_path.string(), settings.kept_atoms.string(), settings.skipped_atoms.string(), atoms, stats,
              settings.skip_H, settings.keep_alternative, settings.skip_hetero_atoms,
              settings.skip_non_standard_amino_acids);
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
    ss << std::setw(6) << std::setprecision(2) << std::setfill(' ') << std::right << std::fixed << 1.0;
    ss << std::setw(6) << std::setprecision(2) << std::setfill(' ') << std::right << std::fixed << 0.0;
    ss << std::setw(12) << std::setfill(' ') << std::right << std::fixed << "C";
    ss << "  ";
    return ss.str();
}

// format information of a single atom in PDB format
std::string PDB_atom_entry(const std::shared_ptr<Atom> &atom) {
    if (std::isinf(atom->coord.x) || std::isinf(atom->coord.y) || std::isinf(atom->coord.z)) return "";
    std::stringstream ss;
    ss << std::setw(6) << std::setfill(' ') << std::left << std::fixed << to_str(atom->record);
    // due to width limitations of this column, higher IDs can only be represented with asterisks
    if (atom->PDB_id <= 99999) {
        ss << std::setw(5) << std::setfill(' ') << std::right << atom->PDB_id;
    } else {
        ss << "*****";
    }
    ss << " ";

    if (atom->full_type.size() == 4 || to_str(atom->type).size() == 2 || std::isdigit(atom->full_type.at(0))) {
        ss << std::setw(4) << std::setfill(' ') << std::left << atom->full_type;
    } else {
        ss << " " << std::setw(3) << std::setfill(' ') << std::left << atom->full_type;
    }
    ss << std::setw(1) << std::setfill(' ') << std::right << atom->alternate_location;
    ss << std::setw(3) << std::setfill(' ') << std::right << atom->residue_name;
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
    ss << std::setw(12) << std::setfill(' ') << std::right << std::fixed << to_str(atom->element);
    ss << std::setw(2) << std::setfill(' ') << std::left << std::fixed << atom->charge;
    return ss.str();
}

// output the information of all atoms in a PDB file
void write_PDB(const std::string &file_path, const std::vector<std::shared_ptr<Atom>> &atoms) {
    std::ofstream file;
    file.open(file_path);
    for (const std::shared_ptr<Atom> &atom: atoms) { file << PDB_atom_entry(atom) << "\n"; }
    file.close();
}