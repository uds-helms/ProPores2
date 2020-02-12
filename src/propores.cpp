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
        // PDB parsing and grid construction
        ProteinGrid grid = ProteinGrid(settings);
        print(1, pore_id_start, "> loaded protein with " + std::to_string(grid.atoms.size()) + " atoms in");
        // run pore/cavity identification
        pore_ID(grid, settings);

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
        }

        // write a file with grid boxes of each pore/cavity cluster for subsequent gate-opening and/or as preparation
        // for later manual or automated analysis
        if (settings.run_axis_trace_preparation) {
            auto axis_preparation = std::chrono::high_resolution_clock::now();
            output_trace_cluster(settings, grid);
            print(1, axis_preparation,
                  "> wrote information of " + std::to_string(grid.clusters.size()) + " pore(s) in");
        }
        // save the pore/cavity grid box clusters for axis-determination
        if (settings.run_axis_trace) clusters = std::move(grid.clusters);
        print(1, pore_id_start, "=> total pore-ID runtime:");
    }

    // determine the axes of pores and cavities
    if (settings.run_axis_trace) {
        if (!settings.run_pore_id) print(header);
        print("\nAxis-Trace");
        // if pore-ID was not performed at the start of this run, obtain the cluster(s) from user provided file(s)
        if (!settings.run_pore_id) {
            auto axis_loading_start = std::chrono::high_resolution_clock::now();
            if (settings.load_cluster_from_single_file)
                cluster_from_file(settings.axis_single_input_file.string(), clusters);
            if (settings.load_cluster_from_directory)
                clusters_from_directory(settings.axis_directory_input_path.string(), clusters);
            print(1, axis_loading_start, "> loaded " + std::to_string(clusters.size()) + " pore(s) in");
        }
        // if at least one cluster could be loaded, try to open all clusteers in the list
        if (!clusters.empty()) {
            auto axis_start = std::chrono::high_resolution_clock::now();
            axis_trace(settings, clusters);
            print(1, axis_start, "=> total axis_trace runtime:");
        }
    }

    // open gates between neighbouring pores/cavities
    if (settings.run_gate_open) {
        if (!settings.run_pore_id && !settings.run_axis_trace) print(header);
        print("\nGate-Open");
        // if pore-ID was not performed at the start of this run, obtain the potential gates from user provided file(s)
        if (!settings.run_pore_id) {
            auto gate_loading_start = std::chrono::high_resolution_clock::now();
            if (settings.load_gate_from_single_file)
                gate_from_file(settings.gate_single_input_file_path.string(), gates);
            if (settings.load_gates_from_directory)
                gates_from_directory(settings.gate_input_directory_path.string(), gates);
            print(1, gate_loading_start, "> loaded " + std::to_string(gates.size()) + " gate(s) in");
        }
        // if at least one gate could be loaded, try to open all gates in the list
        if (!gates.empty()) {
            auto gate_open_start = std::chrono::high_resolution_clock::now();
            gate_open(settings, gates);
            print(1, gate_open_start, "=> total gate-open runtime:");
        }
    }
    print(0, initial, "\n=> total runtime:");
}