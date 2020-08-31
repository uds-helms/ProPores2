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

#ifndef PROPORES_READER_H
#define PROPORES_READER_H

#include <set>
#include <map>
#include <vector>
#include "atom.h"
#include "enums.h"
#include "basics.h"
#include "settings.h"

struct PDBReaderStats {
    size_t valid = 0;
    size_t skipped_H = 0;
    size_t skipped_alternative = 0;
    size_t skipped_hetero = 0;
    size_t skipped_non_standard_amino_acids = 0;
    size_t line_too_short = 0;
    size_t atom_type_parse_error = 0;
    size_t no_vdw_radius = 0;
    size_t no_atomic_mass = 0;

    [[nodiscard]] size_t total_skipped() const {
        return skipped_H + skipped_alternative + skipped_hetero + skipped_non_standard_amino_acids + line_too_short + atom_type_parse_error + no_vdw_radius + no_atomic_mass;
    }
    
    void write(const std::string &file_path) const {
        add_entry(file_path, 1, "PDB parsing stats");
        add_comment(file_path, 2, "number of atoms that were successfully extracted from the PDB file");
        add_entry(file_path, 2, "atoms", valid);
        add_comment(file_path, 2, "number of atoms that were skipped for various reasons");
        add_entry(file_path, 2, "total skipped atoms", total_skipped());
        add_comment(file_path, 2, "reasons why atoms were skipped");
        add_entry(file_path, 2, "skipped HETATM entries", skipped_hetero);
        add_entry(file_path, 2, "record line too short", line_too_short);
        add_entry(file_path, 2, "atom type could not be determined", atom_type_parse_error);
        add_entry(file_path, 2, "skipped non-standard amino acids in ATOM records", skipped_non_standard_amino_acids);
        add_entry(file_path, 2, "no van der Waals radius information", no_vdw_radius);
        add_entry(file_path, 2, "no atomic mass information", no_atomic_mass);
        add_entry(file_path, 2, "skipped hydrogen atoms", skipped_H);
        add_entry(file_path, 2, "skipped alternative atom locations", skipped_alternative);
    }
};

struct PoreInfo {
    size_t id = 0;
    bool is_pore = true;
    double volume = 0.0;
    std::map<ResidueType, size_t> lining;
    size_t axis_count = 0;
    double axis_length = 0.0;
    double axis_distance = 0.0;
    double axis_score = 0.0;

};

// extract all valid atom entries from a given PDB input file, allow to skip H-atoms and to load alternative locations
// of atoms in addition to the primary location
void parse_PDB(const std::string &pdb_path, const std::string &kept_path, const std::string &skipped_path,
               std::vector<std::shared_ptr<Atom>> &atoms, PDBReaderStats &stats,
               HAtomTag h_atom,
               HeteroAtomTag hetero,
               bool keep_alternative,
               bool skip_non_standard_amino_acids);

// wrapper that allows calling parse_PDB with the program settings
void parse_PDB(const Settings &settings, std::vector<std::shared_ptr<Atom>> &atoms, PDBReaderStats &stats);

// format an index and coordinate in pseudo PDB format, can be used to output pore/cavity boxes for visualisation
std::string pseudo_pdb(size_t id, const Vec<double> &coord);

// format information of a single atom in PDB format
std::string PDB_atom_entry(const std::shared_ptr<Atom> &atom);

// output the information of all atoms in a PDB file
void write_PDB(const std::string &file_path, const std::vector<std::shared_ptr<Atom>> &atoms);

// add a comment with comma separated elements to a file
void add_comment(const std::string &file_path, const std::set<std::string> &elements);

//void generate_overview(const Settings &settings);

#endif //PROPORES_READER_H