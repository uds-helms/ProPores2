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

#include <chrono>
#include <iostream>
#include "atom.h"
#include "grid.h"
#include "trace.h"
#include "pore_id.h"
#include "settings.h"
#include "gate_open.h"


int main(int argc, char *argv[]) {
    // parse and check input parameters
    Settings settings = Settings(argc, argv);
    // initialisation of gates and pore/cavity box clusters
    std::vector<Gate> gates;
    std::vector<PoreCluster> clusters;
    // header
    std::string header = "\nPROPORES 2.0 Copyright (C) 2020 Markus Hollander";
    // start the run time clock
    auto initial = std::chrono::high_resolution_clock::now();
    // PORE-ID
    // run pore/cavity identification, if specified
    if (settings.run_pore_id) {
        print(header);
        auto pore_id_start = std::chrono::high_resolution_clock::now();
        print("\nPore-ID");
        add_entry(settings.pore_log, 1, "started", true);
        // PDB parsing and grid construction
        ProteinGrid grid = ProteinGrid(settings);
        print(1, pore_id_start, "> loaded " + settings.pdb_name + " with " + std::to_string(grid.atoms.size()) + " atoms in");

        add_entry(settings.pore_log, 1, "number of grid boxes", grid.boxes);
        add_entry(settings.pore_log, 1, "PDB parsing runtime", pore_id_start);
        // run pore/cavity identification
        pore_ID(grid, settings);

        if (settings.run_axis_preparation || settings.run_gate_preparation || settings.run_gate_open) {
            add_comment(settings.pore_log, 1, "preparation");
        }

        // compute a list of potential gate between two neighbouring pores/cavities for subsequent gate-opening and/or
        // as preparation for later manual or automated analysis
        if (settings.run_gate_preparation || settings.run_gate_open) {
            // extract the gates from the grid
            auto gate_preparation_start = std::chrono::high_resolution_clock::now();
            std::string gate_steps = "> prepared";
            gates_from_grid(grid, gates);
            // if gate preparation is enabled, write a file with shared lining residues for each potential gate
            if (settings.run_gate_preparation) {
                output_gates(settings, gates);
                gate_steps += " and wrote";
            }
            print(1, gate_preparation_start,
                  gate_steps + " information of " + std::to_string(gates.size()) + " gate(s) in");
            add_entry(settings.pore_log, 1, "potential gates", gates.size());
            add_entry(settings.pore_log, 1, "gate preparation", gate_preparation_start);
        }

        // write a file with grid boxes of each pore/cavity cluster for subsequent gate-opening and/or as preparation
        // for later manual or automated analysis
        if (settings.run_axis_preparation) {
            auto axis_preparation = std::chrono::high_resolution_clock::now();
            output_trace_cluster(settings, grid);
            print(1, axis_preparation,
                  "> wrote information of " + std::to_string(grid.clusters.size()) + " pore(s) in");
            add_entry(settings.pore_log, 1, "axis trace preparation", axis_preparation);
        }
        // save the pore/cavity grid box clusters for axis-determination
        if (settings.run_axis_trace) clusters = std::move(grid.clusters);
        print(1, pore_id_start, "=> total pore-ID runtime:");
        add_comment(settings.pore_log, 1, "finished");
        add_entry(settings.pore_log, 1, "finished", true);
        add_entry(settings.pore_log, 1, "total runtime", pore_id_start);
    }

    // determine the axes of pores and cavities
    if (settings.run_axis_trace) {
        if (!settings.run_pore_id) print(header);
        print("\nAxis-Trace");
        auto axis_start = std::chrono::high_resolution_clock::now();

        add_entry(settings.axis_log, 1, "started", true);
        // if pore-ID was not performed at the start of this run, obtain the cluster(s) from user provided file(s)
        if (!settings.run_pore_id) {
            if (settings.load_cluster_from_single_file)
                cluster_from_file(settings.axis_single_input.string(), clusters);
            if (settings.load_cluster_from_directory)
                clusters_from_directory(settings.axis_directory_input.string(), clusters);
            print(1, axis_start, "> loaded " + std::to_string(clusters.size()) + " pore(s) in");
        }

        // if at least one cluster could be loaded, try to open all clusters in the list
        if (!clusters.empty()) {
            add_entry(settings.axis_log, 1, "pores", clusters.size());
            axis_trace(settings, clusters);
            print(1, axis_start, "=> total axis_trace runtime:");
        }

        add_comment(settings.axis_log, 1, "finished");
        add_entry(settings.axis_log, 1, "finished", true;
        add_entry(settings.axis_log, 1, "total runtime", axis_start);
    }

    // open gates between neighbouring pores/cavities
    if (settings.run_gate_open) {
        if (!settings.run_pore_id && !settings.run_axis_trace) print(header);
        print("\nGate-Open");
        auto gate_open_start = std::chrono::high_resolution_clock::now();

        add_entry(settings.gate_log, 1, "started", true);
        // if pore-ID was not performed at the start of this run, obtain the potential gates from user provided file(s)
        if (!settings.run_pore_id) {
            if (settings.load_gate_from_single_file)
                gate_from_file(settings.gate_single_input.string(), gates);
            if (settings.load_gates_from_directory)
                gates_from_directory(settings.gate_directory_input.string(), gates);
            print(1, gate_open_start, "> loaded " + std::to_string(gates.size()) + " gate(s) in");
        }
        // if at least one gate could be loaded, try to open all gates in the list
        if (!gates.empty()) {
            add_entry(settings.gate_log, 1, "gates", gates.size());
            gate_open(settings, gates);
            print(1, gate_open_start, "=> total gate-open runtime:");
        }

        add_comment(settings.gate_log, 1, "finished");
        add_entry(settings.gate_log, 1, "finished", true);
        add_entry(settings.gate_log, 1, "total runtime", gate_open_start);
    }
    print(0, initial, "\n=> total runtime:");
}