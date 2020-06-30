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

#ifndef PROPORES_GATE_H
#define PROPORES_GATE_H

#include <set>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <exception>
#include <unordered_map>
#include "atom.h"
#include "enums.h"
#include "vector.h"
#include "reader.h"
#include "settings.h"

// residues have a backbone and side chain definition, e.g. ARG = N-CA-CB-CG-CD-NE-CZ
// this represents a component of this definition, e.g. N or CA (C-alpha)
struct DefinitionComponent {
    // full atom type and atom index in the residue atom list
    // e.g. in the above example, type = N and index = 0, or type = CA and index = 3
    size_t index = 0;
    std::string type;

    // constructor
    explicit DefinitionComponent(std::string atom_type) : type(std::move(atom_type)) {}
};

// represents a dihedral angle in a residue
// for example if the residue backbone and side chain definition is ARG = N-CA-CB-CG-CD-NE-CZ, the first dihedral side
// chain angle would be between the plane spanned by N-CA and CA-CB and the plane spanned by CA-CB and CB-CG
// in this example, the side chain would be rotated starting with CB to change the angle, making CB the rotation centre
// the rotation axis would be formed by the vector between CA and CB, making CA the start of the rotation axis
struct Dihedral {
    // index in the dihedral angle list of the residue
    size_t index;
    // indices of the rotation axis start and rotation point atoms in the residue atom list
    size_t rotation_point;
    size_t rotation_axis_start;
    // the dihedral angle in degrees
    double angle;

    // constructor
    Dihedral(const size_t angle_index, const size_t previous, const size_t current, const double angle) :
            index(angle_index),
            rotation_axis_start(previous),
            rotation_point(current),
            angle(angle) {}
};

struct BasicAtom {
    // index of the full-atom in the gate's atom list
    size_t index;
    // full atom type, e.g. CA for C-alpha
    std::string type;
    // atom coordinates
    Vec<double> coord;
    // atom van der Waals radius
    double vdw;
    // true if the atom is part of the residue chain, false if it is part of the residue backbone
    bool is_side_chain = false;

    // constructor
    explicit BasicAtom(const std::shared_ptr<Atom> &atom) :
            index(atom->id),
            type(atom->full_type),
            coord(atom->coord),
            vdw(atom->vdw) {}
};

// represents a single rotamer of a residue
struct Rotamer {
    // unique rotamer ID within the residue
    size_t id;
    // basic information on atoms in this rotamer
    std::vector<BasicAtom> atoms;
    // desired dihedral angles in degrees
    std::vector<int> angles;
    // true if the rotamer was eliminated during dead-end elimination and should not be considered further
    bool eliminated = false;
    // self-energy of the rotamer
    double energy = 0;

    // constructor
    Rotamer(const size_t identifier, std::vector<BasicAtom> atom_list, std::vector<int> desired_angles) :
            id(identifier),
            angles(std::move(desired_angles)),
            atoms(std::move(atom_list)) {}
};


struct Residue {
    // ID in the gate residue list
    size_t id;
    // residue ID in the protein PDB
    size_t pdb_id;
    // residue type, e.g. MET
    ResidueType type;
    // protein chain to which the residue belongs
    std::string chain;
    // basic information on atoms in this residue
    std::vector<BasicAtom> atoms;
    // backbone and side chain definition, e.g. ARG = N-CA-CB-CG-CD-NE-CZ
    std::vector<DefinitionComponent> defs;
    // index of the first side chain atom in 'atoms'
    size_t side_chain_start = 0;
    // list of dihedral angles
    std::vector<Dihedral> angles;
    // collision free rotamers of this residue
    std::vector<Rotamer> rotamers;

    // constructor
    Residue(const size_t id_in_gate, std::string res_chain, const size_t res_id, const ResidueType res_type) :
            id(id_in_gate),
            chain(std::move(res_chain)),
            pdb_id(res_id),
            type(res_type) {}

    // add an atom to the atom list
    void add_atom(const std::shared_ptr<Atom> &atom) { atoms.emplace_back(atom); }

    // set up the backbone and side chain atom indices
    bool setup_indices() {
        // load the main backbone and side chain component definitions
        for (const std::string &def: residue_definition(type)) { defs.emplace_back(def); }
        // in ILE, rotating the side chain starting with CG1 should not rotate the branch with CG2
        // => need to move CG2 in front of CG1, if that is not already the case
        if (type == ILE) {
            size_t cg1_index = 0;
            size_t cg2_index = 0;
            for (size_t i = 0; i < atoms.size(); i++) {
                if (atoms[i].type == "CG1") cg1_index = i;
                if (atoms[i].type == "CG2") cg2_index = i;
            }
            if (cg1_index < cg2_index) std::iter_swap(atoms.begin() + cg1_index, atoms.begin() + cg2_index);
        }
        // build the indices of the backbone and side chain definition
        for (size_t i = 0; i < atoms.size(); i++) {
            BasicAtom &atom = atoms[i];
            // side chains start with the C beta (= CB) atom
            if (atom.type == "CB") side_chain_start = i;
            // since in PDBs the atoms are ordered by their position from the amino group, the side chain cannot
            // start at index 0 and thus an index of 0 means the side chain start has not been found yet
            if (i >= side_chain_start && side_chain_start > 0) atom.is_side_chain = true;
            // check if the current atom is the start of a backbone or side chain chunk and set the index of the
            // corresponding definition accordingly
            for (DefinitionComponent &def: defs) { if (atom.type == def.type) { def.index = i; }}
        }
        // validate that the indices are strictly increasing, if not, processing this residue becomes problematic
        for (size_t i = 1; i < defs.size(); i++) {
            if (defs[i].index <= defs[i - 1].index) return false;
            DefinitionComponent &def = defs[i];
        }
        return true;
    }

    [[nodiscard]] size_t size() const { return rotamers.size(); }

    [[nodiscard]] std::string str() const {
        return chain + "-" + std::to_string(pdb_id) + "-" + to_str(type) + "-" + std::to_string(size()) + "rot";
    }
};

// represents a (potential) gate between two neighbouring pores/cavities
struct Gate {
    // IDs of the two pores/cavities joined by this gate
    size_t pore_1;
    size_t pore_2;
    // all atoms in the protein, extracted from the input PDB file
    std::vector<std::shared_ptr<Atom>> atoms;
    // gating residues between the two pores/cavities
    std::vector<Residue> residues;
    // estimation of how difficult (= how slow) computing the best rotamer combination is likely going to be
    GateDifficulty difficulty = DIFFICULTY_ERROR;
    // coordinates of the gate centre, in terms of the geometric centre of C-alpha atoms in the gating residues
    Vec<double> centre;
    // KD-tree root, starting point of KD-nearest neighbour search
    std::shared_ptr<Atom> root;

    // constructor, no protein atoms are loaded initially
    Gate(const size_t id_1, const size_t id_2, std::vector<std::tuple<std::string, size_t, std::string>> &gating_res) :
            pore_1(id_1),
            pore_2(id_2) {
        // add the gating residues to the gate
        size_t next_id = 0;
        for (const auto &[chain, id, res]: gating_res) {
            // do not add invalid residues
			ResidueType type = to_residue_type(res);
            if (type == INVALID_RESIDUE) continue;
            residues.emplace_back(next_id, chain, id, type);
            next_id++;
        }
    }

    // load the protein atoms from the input PDB file, and match them to the gating residues
    void load_atoms(Settings &settings) {
        // parse the input PDB file and extract the protein atoms
        PDBReaderStats stats = PDBReaderStats();
        parse_PDB(settings, atoms, stats);
        // log PDB parsing stats if that was not already done for a previous gate
        if (!settings.gate_already_logged_pdb_parsing_stats) {
            stats.write(settings.gate_log.string());
            settings.gate_already_logged_pdb_parsing_stats = true;
        }
        // associate atoms with gating residues and determine the position of the C-alpha atoms in the gating residues
        size_t number_c_alpha = 0;
        for (const std::shared_ptr<Atom> &atom: atoms) {
            for (Residue &res: residues) {
                // the atom is not part of a gating residue
                if (atom->chain != res.chain || atom->residue_id != res.pdb_id) continue;
                // accumulate the position of C-alpha atoms
                if (atom->full_type == "CA") {
                    centre += atom->coord;
                    number_c_alpha++;
                }
                // do not need to store the atoms of alanine, glycine or proline since they are going to be removed
                // from the gate later anyway due to missing side chains or lack of rotamers
                if (res.type == ALA || res.type == GLY || res.type == PRO || res.type == INVALID_RESIDUE) continue;
                res.add_atom(atom);
            }
        }
        // compute the centre of the gate in terms of C-alpha gate atoms
        if (number_c_alpha > 0) centre /= double(number_c_alpha);

        // only keep residues that have proper side chains/rotamers and whose indices could be setup
        std::vector<Residue> valid_residues;
        for (Residue &residue: residues) {
            // alanine, glycine and proline do not have side chains or none that can be converted to rotamers
            if (residue.type == ALA || residue.type == GLY || residue.type == PRO || residue.type == INVALID_RESIDUE) continue;
            // some error occurred while the side chains were set up
            if (!residue.setup_indices()) continue;
            // set the coordinates of side chain atoms to infinity to (temporarily) remove them from the protein
            // and nearest neighbour search to avoid collisions with rotamers
            for (const BasicAtom &atom: residue.atoms) {
                if (atom.is_side_chain) atoms[atom.index]->coord = Vec<double>(INFINITY, INFINITY, INFINITY);
            }
            valid_residues.push_back(std::move(residue));
        }
        residues = valid_residues;

        // recursively build the KD-tree for nearest neighbour search
        root = KD_tree(0, atoms.begin(), atoms.end() - 1);
        // bring the atoms back into the original order, meaning index 0 is atom with ID 0
        sort_atoms(atoms);
    }

    // generate the name of the gate based on the two pore/cavity IDs
    [[nodiscard]] std::string name() const {
        return std::to_string(std::min(pore_1, pore_2)) + "_and_" + std::to_string(std::max(pore_1, pore_2));
    }

    // string representation of all residues in the gate
    [[nodiscard]] std::string residue_str() const {
        if (residues.empty()) return "no residues";
        std::string str;
        for (size_t i = 0; i < residues.size(); i++) {
            const Residue &res = residues[i];
            if (i > 0) str += ", ";
            str += res.str();
        }
        return str;
    }
};

#endif //PROPORES_GATE_H
