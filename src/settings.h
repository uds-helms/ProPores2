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

#ifndef PROPORES_SETTINGS_H
#define PROPORES_SETTINGS_H

#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include "enums.h"
#include "basics.h"

namespace fs = std::filesystem;

struct Settings {
    // INPUT and OUTPUT
    // PDB parsing
    fs::path pdb_path;
    std::string pdb_name;
    // main output specifications
    fs::path out_dir;
    std::string output_name;
    // pore-ID paths
    fs::path pore_dir;
    fs::path lining_dir;
    // axis trace paths
    fs::path axis_dir;
    fs::path axis_preparation_dir;
    fs::path axis_single_input;
    fs::path axis_directory_input;
    // gate open paths
    fs::path gate_dir;
    fs::path gate_preparation_dir;
    fs::path gate_single_input;
    fs::path gate_directory_input;
    fs::path gate_merged_pdb;
    // logging
    fs::path kept_atoms;
    fs::path skipped_atoms;
    fs::path pore_log;
    fs::path axis_log;
    fs::path gate_log;
    // static input
    fs::path rotamer_lib = fs::path("rotamer_lib");

    // PARAMETERS
    // PORE-ID and AXIS-TRACE
    // grid dimensions and pore-ID parameters
    double box_length = 1.0;
    double solvent_radius = 1.2;
    double probe_radius = 1.4;
    // PORE ID
    // threshold for the minimum volume of a pore/cavity in Angstrom
    double volume_threshold = 50;
    // highest number of grid boxes before switching from cylinder trace to cylinder standalone
    size_t box_threshold = 1500000;
    // highest number of atoms before switching from cylinder trace to cylinder standalone
    size_t atom_pair_threshold = 30000;
    // number of cores for standalone pore-ID
    size_t cores = 1;
    // AXIS-TRACE
    // threshold for the minimum number of grid boxes in a pore surface patch
    size_t surface_patch_threshold = 30;
    double perturb_value = 0.001;
    // GATE-OPEN
    // acceptable van der Waals radius overlap in Angstrom
    double clash_tolerance = 0.75;
    // acceptable wan der Waals radius overlap in Angstrom for the gamma-atom in a gate residue side chain
    double clash_tolerance_gamma = 1.0;

    // FLAGS
    // main program components
    bool run_pore_id = false;
    bool run_axis_trace = false;
    bool run_gate_open = false;
    // PDB-PARSING
    // false if only the main location of an atom should be considered, true if the alternative location should be kept
    bool keep_alternative = false;
    // true if H-atoms should be discarded from the input PDB file
    HAtomTag h_atom_tag = KEEP_ALL_H_ATOMS;
    // true if HETATM records should be ignored
    HeteroAtomTag hetero_atom_tag = KEEP_ALL_HETERO_ATOMS;
    // true if non-standard amino acids should be ignored
    bool skip_non_standard_amino_acids = false;
    // PORE-ID
    // whether to generate output files that can be used for axis-trace or gate-open
    PreparationTag preparation = AXIS_AND_GATE;
    // whether to remove pores, cavities or keep both
    RemoveTag pore_filter = REMOVE_NOTHING;
    // whether to run cylinder trace, cylinder standalone or try to auto-detect the best approach
    ComputationMode computation_mode = AUTODETECT;
    // whether to load gates (gate-open) or pore/cavity clusters (axis-trace) from a single file or directory
    bool load_gate_from_single_file = false;
    bool load_gates_from_directory = false;
    bool load_cluster_from_single_file = false;
    bool load_cluster_from_directory = false;
    // GATE-OPEN
    // true if the gate-open difficulty should be re-assessed after rotamers have been generated and check for collision
    bool re_estimate = false;
    // gates with this difficulty or higher are skipped
    GateDifficulty difficulty_threshold = DIFFICULTY_ERROR;
    // for each gate the PDB is parsed again, true if the PDB parsing stats were already logged to avoid duplicates
    bool gate_already_logged_pdb_parsing_stats = false;

    // CONSTRUCTORS
    Settings() = default;

    // from command line input
    Settings(int argc, char *argv[]) {
        // COMMAND LINE PARSING
        // get the command line arguments
        args = std::vector<std::string>(argv, argv + argc);
        // extract the program run flags
        run_pore_id = option_exists("pore-id");
        run_axis_trace = option_exists("axis-trace");
        run_gate_open = option_exists("gate-open");
        // extract file paths and output specifications
        output_name = get_option("--name");
        std::string input_pdb = get_option("-i");
        std::string output_dir = get_option("-o");
        std::string axis_input_single = get_option("-ts");
        std::string axis_input_dir = get_option("-td");
        std::string gate_input_single = get_option("-gs");
        std::string gate_input_dir = get_option("-gd");
        std::string custom_rotamer_path = get_option("--rotamers");
        // extract simple flags
        keep_alternative = option_exists("--keep-alternative");
        skip_non_standard_amino_acids = option_exists("--skip-non-std-amino-acids");
        re_estimate = option_exists("--re-estimate");
        bool help_short = option_exists("-h");
        bool help_long = option_exists("--help");
        // extract complex flags
        std::string h_atom_flag = get_option("--h-atom");
        std::string hetero_flag = get_option("--hetero");
        std::string mode_flag = get_option("--mode");
        std::string pore_filter_flag = get_option("--pore-filter");
        std::string preparation_flag = get_option("--preparation");
        std::string difficulty_flag = get_option("--difficulty");

        if (help_short || help_long) print_help("");

        if (!run_pore_id && !run_axis_trace && !run_gate_open) print_help("No program part(s) specified.");

        if (input_pdb.empty()) print_help("No protein PDB file (-i) was given.");

        try {
            pdb_path = fs::absolute(std::filesystem::path(input_pdb));
            if (!fs::is_regular_file(pdb_path)) {
                print_help("The given protein PDB file path (-i) does not point to an existing file. "
                           "Make sure to wrap file paths with spaces in quotation marks.");
            }
        } catch (...) { print_help("The given protein PDB file path (-i) is invalid."); }

        if (output_dir.empty()) print_help("No output directory (-o) was given.");

        try {
            out_dir = fs::absolute(std::filesystem::path(output_dir));
        } catch (...) { print_help("The given output directory (-o) is invalid."); }

        // extract numerical parameters
        set_option("-b", box_length, "grid box length (-b)");
        set_option("-s", solvent_radius, "solvent radius (-s)");
        set_option("-p", probe_radius, "probe radius (-p)");
        set_option("-v", volume_threshold, "volume threshold (-v)");
        set_option("-spt", surface_patch_threshold, "surface patch threshold (-spt)");
        set_option("-ct", clash_tolerance, "clash tolerance (-ct)");
        set_option("--cores", cores, "number of multi-threading cores (--cores)");
        // check if the numerical parameters are within the expected value range
        if (box_length < 0) print_help("The option '-b <decimal>' was negative.");
        if (solvent_radius < 0) print_help("The option '-s <decimal>' was negative.");
        if (probe_radius < 0) print_help("The option '-p <decimal>' was negative.");
        if (volume_threshold < 0) print_help("The option '-v <decimal>' was negative.");
        if (surface_patch_threshold < 0) print_help("The option '-spt <decimal>' was negative.");
        if (clash_tolerance < 0) print_help("The option '-ct <decimal>' was negative.");

        // check and set complex flags
        if (!h_atom_flag.empty()) {
            h_atom_tag = to_h_atom_tag(h_atom_flag);
            if (h_atom_tag == INVALID_H_ATOM_TAG) {print_help("The \"--h-atom\" option was invalid.");}
        }

        if (!hetero_flag.empty()) {
            hetero_atom_tag = to_hetero_atom_tag(hetero_flag);
            if (hetero_atom_tag == INVALID_HETERO_ATOM_TAG) {print_help("The \"--hetero\" option was invalid.");}
        }

        if (!mode_flag.empty()) {
            computation_mode = to_computation_mode(mode_flag);
            if (computation_mode == INVALID_COMPUTATION_MODE) {print_help("The \"--mode\" option was invalid.");}
        }

        if (!pore_filter_flag.empty()) {
            pore_filter = to_remove_tag(pore_filter_flag);
            if (pore_filter == INVALID_REMOVE_TAG) {print_help("The \"--pore-filter\" option was invalid.");}
        }

        if (!preparation_flag.empty()) {
            preparation = to_preparation_tag(preparation_flag);
            if (preparation == INVALID_PREPARATION_TAG) {print_help("The \"--preparation\" option was invalid.");}
        }

        if (!difficulty_flag.empty()) {
            difficulty_threshold = to_gate_difficulty(difficulty_flag);
            if (difficulty_threshold == DIFFICULTY_ERROR) {print_help("The \"--difficulty\" option was invalid.");}
        }

        // extract the PDB name
        pdb_name = pdb_path.stem().string();
        // set the output name to the PDB name if none was given on the command line
        if (output_name.empty()) output_name = pdb_name;
        // construct the output directory path
        out_dir /= fs::path(output_name);
        if (!output_name.empty() && output_name.back() != '_') output_name += '_';

        // generate output sub-directory paths
        pore_dir = out_dir / fs::path("pores");
        lining_dir = out_dir / fs::path("lining");
        axis_dir = out_dir / fs::path("axes");
        gate_dir = out_dir / fs::path("gate_open");
        axis_preparation_dir = out_dir / fs::path("axis_trace_input");
        gate_preparation_dir = out_dir / fs::path("gate_open_input");
        gate_merged_pdb = gate_dir / fs::path(output_name + "all_in_one.pdb");

        // process axis-trace file paths
        if (run_axis_trace) {
            load_cluster_from_single_file = !axis_input_single.empty();
            load_cluster_from_directory = !axis_input_dir.empty();
            // both a single axis-trace file and an axis-trace directory were provided
            if (load_cluster_from_directory && load_cluster_from_single_file) {
                print_help("The options '-ts <path>' and '-td <path>' are mutually exclusive.");
            }
            // pore-ID is not enabled and there was no custom axis-trace input
            if (!run_pore_id && axis_input_single.empty() && axis_input_dir.empty()) {
                // check if maybe the axis trace preparation directory exists and contains input
                if (fs::exists(axis_preparation_dir) && !fs::is_empty(axis_preparation_dir)) {
                    load_cluster_from_directory = true;
                    axis_directory_input = axis_preparation_dir;
                } else {
                    print_help("When pore-ID is not run, axis-trace input must be provided "
                               "(either '-ts <path>' or '-td <path>').");
                }
            }
            // pore-ID is enabled and there was at least one custom axis-trace input
            else if (run_pore_id && (load_cluster_from_directory || load_cluster_from_single_file)) {
                print_help("When pore-ID is run, do not provide axis determination input files ('-ts' or '-td'). "
                           "For running axis-trace with user defined input, do not run pore-ID.");
            }
            // at this point we know that if the single axis-trace file is provided, then pore-ID is not enabled
            // and there is no axis-trace directory given either, so we just have to check the file path
            else if (load_cluster_from_single_file) {
                axis_single_input = fs::path(axis_input_single);
                if (!fs::is_regular_file(axis_single_input)) {
                    print_help("The axis trace input file path (-ts) does not point to an existing file. "
                               "Make sure to wrap file paths with spaces in quotation marks.");
                } else if (fs::is_empty(axis_single_input)) {
                    print_help("The axis trace input file path (-ts) points to an empty file.");
                }
            }
            // at this point we know that if the axis-trace directory is provided, then pore-ID is not enabled
            // and there is no single axis-trace file given either, so we just have to check the directory path
            else if (load_cluster_from_directory) {
                axis_directory_input = fs::path(axis_input_dir);
                if (!fs::is_directory(axis_directory_input)) {
                    print_help("The axis trace input directory path (-td) does not point to an existing directory. "
                               "Make sure to wrap directory paths with spaces in quotation marks.");
                } else if (fs::is_empty(axis_directory_input)) {
                    print_help("The axis trace input directory path (-td) points to an empty directory.");
                }
            }
        }

        // process gate-open file paths
        if (run_gate_open) {
            load_gate_from_single_file = !gate_input_single.empty();
            load_gates_from_directory = !gate_input_dir.empty();
            // both a single gate-open file and an gate-open directory were provided
            if (load_gate_from_single_file && load_gates_from_directory) {
                print_help("The options '-gs <path>' and '-gd <path>' are mutually exclusive.");
            }
            // pore-ID is not enabled and there was no custom gate-open input
            if (!run_pore_id && gate_input_single.empty() && gate_input_dir.empty()) {
                if (fs::exists(gate_preparation_dir) && !fs::is_empty(gate_preparation_dir)) {
                    load_gates_from_directory = true;
                    gate_directory_input = gate_preparation_dir;
                } else {
                    print_help("When pore-ID is not run, gate-open input must be provided "
                               "(either '-gs <path>' or '-gd <path>').");
                }
            }
            // pore-ID is enabled and there was at least one custom gate-open input
            else if (run_pore_id && (load_gates_from_directory || load_gate_from_single_file)) {
                print_help("When pore-ID is run, do not provide gate-open input files ('-gs' or '-gd'). "
                           "For running gate-open with user defined input, do not run pore-ID.");
            }
            // at this point we know that if the single gate-open file is provided, then pore-ID is not enabled
            // and there is no gate-open directory given either, so we just have to check the file path
            else if (load_gate_from_single_file) {
                gate_single_input = fs::path(gate_input_single);
                if (!fs::is_regular_file(gate_single_input)) {
                    print_help("The gate-open input file path (-gs) does not point to an existing file. "
                               "Make sure to wrap file paths with spaces in quotation marks.");
                } else if (fs::is_empty(gate_single_input)) {
                    print_help("The gate-open input file path (-gs) points to an empty file.");
                }
            }
            // at this point we know that if the gate-open directory is provided, then pore-ID is not enabled
            // and there is no single gate-open file given either, so we just have to check the directory path
            else if (load_gates_from_directory) {
                gate_directory_input = fs::path(gate_input_dir);
                if (!fs::is_directory(gate_directory_input)) {
                    print_help("The gate-open input directory path (-gd) does not point to an existing directory. "
                               "Make sure to wrap directory paths with spaces in quotation marks.");
                } else if (fs::is_empty(gate_directory_input)) {
                    print_help("The gate-open input directory path (-gd) points to an empty directory.");
                }
            }
        }

        // try to delete the previous output depending on which program parts are enabled
        if (fs::exists(out_dir)) {
            // for pore-ID all previous output is deleted
            if (run_pore_id) {
                try {
                    fs::remove_all(out_dir);
                } catch (...) {
                    print_help("The content of the already existing output directory at (" + out_dir.string()
                               + ") could not be deleted. Check if it or the content is open in another program.");
                }
            } else {
                // if only axis-trace and/or gate-open are enabled, only delete previous axis-trace/gate-open output
                if (run_axis_trace && fs::exists(axis_dir)) {
                    try {
                        fs::remove_all(axis_dir);
                    } catch (...) {
                        print_help("The content of the already existing axis-trace directory at '"
                                   + axis_dir.string() + "' could not be deleted. Check if it or the content is "
                                                          "open in another program.");
                    }
                }
                if (run_gate_open && fs::exists(gate_dir)) {
                    try {
                        fs::remove_all(gate_dir);
                    } catch (...) {
                        print_help("The content of the already existing gate-open directory at '"
                                   + gate_dir.string() + "' could not be deleted. Check if it or the content is "
                                                          "open in another program.");
                    }
                }
            }
        }
        // create the main output directory
        fs::create_directories(out_dir);
        // create the sub-directories for enabled program parts
        if (run_pore_id) {
            fs::create_directories(pore_dir);
            fs::create_directories(lining_dir);
        }
        if (run_axis_trace) fs::create_directories(axis_dir);
        if (run_gate_open) fs::create_directories(gate_dir);
        if (preparation == AXIS_AND_GATE || preparation == ONLY_AXIS) fs::create_directories(axis_preparation_dir);
        if (preparation == AXIS_AND_GATE || preparation == ONLY_GATE) fs::create_directories(gate_preparation_dir);

        // if gate-open is enabled, check if the rotamer library is where it should be
        if (run_gate_open) {
            std::string flavour;
            // check if a custom rotamer library path was specified
            if (!custom_rotamer_path.empty()) {
                rotamer_lib = fs::path(custom_rotamer_path);
                flavour = "custom ";
            }
            if (!fs::is_directory(rotamer_lib)) {
                print_help("The " + flavour + "rotamer library could not be found at: " + rotamer_lib.string());
            } else if (fs::is_empty(rotamer_lib)) {
                print_help("The " + flavour + " rotamer library at '" + rotamer_lib.string() + "' is empty.");
            }
        }

        // generate logging paths
        kept_atoms = out_dir / fs::path("used_atoms.pdb");
        skipped_atoms = out_dir / fs::path("skipped_atoms.pdb");
        pore_log = out_dir / fs::path("pore_ID_log.yaml");
        axis_log = out_dir / fs::path("axes_trace_log.yaml");
        gate_log = out_dir / fs::path("gate_open_log.yaml");

        // log the program component parameters
        if (run_pore_id) write_pore_id_parameters();
        if (run_axis_trace) write_axes_trace_parameters();
        if (run_gate_open) write_gate_open_parameters();

        // copy input PDB to output folder for convenience, if it's not already in there
        fs::path pdb_copy = out_dir / fs::path(pdb_name + ".pdb");
        if (!fs::exists(pdb_copy)) fs::copy(pdb_path, out_dir / fs::path(pdb_name + ".pdb"));

    }

    // generates the file path with all rotamer definitions for the given residue type
    [[nodiscard]] std::string rotamer_file_path(const std::string &residue_type) const {
        fs::path file_path = rotamer_lib / fs::path(residue_type + "_trimed.txt");
        if (fs::exists(file_path)) return file_path.string();
        return "";
    }

private:
    std::vector<std::string> args;
    const size_t col = 26;
    const std::string sep = "  ";
    const size_t line_width = 100;

    [[nodiscard]] bool option_exists(const std::string &opt) const {
        return std::find(args.begin(), args.end(), opt) != args.end();
    }

    [[nodiscard]] std::string get_option(const std::string &opt) const {
        auto it = std::find(args.begin(), args.end(), opt);
        if (it == args.end() || it == args.end() - 1) return "";
        return *(it + 1);
    }

    void set_option(const std::string &opt, double &member, const std::string &msg) {
        if (!option_exists(opt)) return;
        try {
            member = stod(get_option(opt));
        } catch (...) { print_help("Given " + msg + " is not a number"); }
        if (member < 0) { print_help("Given " + msg + " is < 0."); }
    }

    void set_option(const std::string &opt, size_t &member, const std::string &msg) {
        if (!option_exists(opt)) return;
        try {
            member = stoi(get_option(opt));
        } catch (...) { print_help("Given " + msg + " is not a number"); }
        if (member < 0) { print_help("Given " + msg + " is < 0."); }
    }

    void print_help(const std::string &msg) const {
        print_text("PROPORES 2.0 Copyright (C) 2020 Markus Hollander [GNU General Public License v3]");
        print_text("Identification of protein pores, cavities and channels with options for axis determination and "
                   "opening connections between neighbouring pores.");
        std::cout << std::endl;

        if (!msg.empty()) {
            print_text("ERROR: " + msg + " See the following usage description.");
            std::cout << std::endl;
        }

        std::cout << "Usage:" << std::endl;
        std::cout << sep << "propores <run command(s)> -i <path> -o <path> [options...]" << std::endl << std::endl;
        std::cout << "Example: " << std::endl;
        std::cout << sep << "propores pore-id axis-trace -i input/1EA5_R.pdb -o output" << std::endl << std::endl;

        std::cout << "Run commands (at least one, space separated):" << std::endl;
        print_option("pore-id", "Run pore/cavity identification.");
        print_option("axis-trace", "Run pore/cavity axis determination.");
        print_option("gate-open", "Run gate-open for neighbouring pores/cavities.");
        std::cout << std::endl;

        std::cout << "Required:" << std::endl;
        print_option("-i <path>", "Path to the input PDB file.");
        print_option("-o <path>", "Path to the output directory. A sub-directory will be created for the "
                                  "current run. If the sub-directory already exists and pore-ID is run, its content "
                                  "will get deleted!");
        std::cout << std::endl;

        std::cout << "General options:" << std::endl;
        print_option("-h, --help", "Show this help message.");
        print_option("--name <string>", "Name for the output sub-directory and  all output files. "
                                        "[default: PDB name]");
        print_option("--h-atom <number>", "Filter H atoms in the PDB:+ "
                                          "0: keep all H atoms+ "
                                          "1: remove all H atoms+ "
                                          "2: remove only H atoms in ATOM records+ "
                                          "3: remove only H atoms in HETATM records+ "
                                          "[default: keep all H atoms]");
        print_option("--hetero <number>", "Filter HETATM records in the PDB:+ "
                                          "0: keep all hetero atoms+ "
                                          "1: remove all hetero atoms+ "
                                          "2: remove all hetero atoms except dummy atoms+ "
                                          "3: remove only dummy hetero atoms+ "
                                          "[default: keep all hetero atoms]");
        print_option("--keep-alternative", "Load alternative locations of atoms in addition to the primary location. "
                                           "[default: false]");
        print_option("--skip-non-std-amino-acids", "Do not load atoms of non-standard residues (= not one "
                                                   "of the 20 'classical' amino acids) from ATOM records. HETATM "
                                                   "records are unaffected by this setting. [default: false]");

        std::cout << "Pore-ID options:" << std::endl;
        print_option("-b <decimal>", "Grid box length in Angstrom (= resolution). [default: 1.0]");
        print_option("-s <decimal>", "Solvent radius in Angstrom. [default: 1.2]");
        print_option("-p <decimal>", "Probe radius in Angstrom. The higher the value, the deeper pore regions "
                                     "have to be on the protein surface. [default: 1.4]");
        print_option("-v <decimal>", "Volume threshold for pores/cavities in cubic Angstrom. [default: 50]");
        print_option("--mode <number>", "Computation mode for pore identification:+ "
                                        "0: autodetect+ "
                                        "1: ray trace, which is faster for smaller proteins or larger grid box lengths, "
                                        "but comes with increased RAM usage.+ "
                                        "2: standalone, which has a low RAM usage but can be potentially a lot slower+ "
                                        "[default: autodetect]");
        print_option("--pore-filter <number>", "Generate output and analysis for+ "
                                               "0: pores and cavities+ "
                                               "1: only pores+ "
                                               "2: only cavities+ "
                                               "[default: all]");
        print_option("--preparation <number>", "Output files that can later be used for+ "
                                               "0: axis-trace and gate-open+ "
                                               "1: only axis-trace+ "
                                               "2: only gate-open+ "
                                               "3: no preparation+ "
                                               "[default: axis-trace and gate-open]");

        std::cout << std::endl << "Axis-trace options:" << std::endl;
        print_option("-spt <number>", "Minimum area of a pore surface patch  in Angstrom for it to count as a potential "
                                      "pore axis origin. [default: 30]");
        print_option("-ts <path>", "Path to a file with information for determining the axis of a single "
                                   "pore/cavity. Only needed if axis trace is not run together with pore ID. "
                                   "[default: none]");
        print_option("-td <path>", "Path to a directory with file(s) with information for determining the axis of a "
                                   "single pore/cavity each. Only needed if axis trace is not run together with pore "
                                   "ID. [default: none]");

        std::cout << std::endl << "Gate-open options:" << std::endl;
        print_option("-ct <decimal>", "Van der Waals radius overlap tolerance (clash tolerance) for rotamer generation "
                                      "in Angstrom. Lowering the clash tolerance can substantially speed up the "
                                      "gate-open runtime, but too low values can also eliminate all potential "
                                      "rotamers. [default: 0.75]");
        print_option("-gs <path>", "Path to a file with information for opening the gate between a single pair of "
                                   "two pores/cavities. Only needed if gate-open is not run together with pore ID. "
                                   "[default: none]");
        print_option("-gd <path>", "Path to a directory with file(s) with information for opening the gate between a "
                                   "single pair of two pores/cavities. Only needed if gate-open is not run together "
                                   "with pore ID. [default: none]");
        print_option("--rotamers <path>", "Custom path to the rotamer library. [default: none]");
        print_option("--difficulty <number>", "Filter gates by their (potential) difficulty+ "
                                              "0: try to open all gates+ "
                                              "1: skip very time intensive gates+ "
                                              "2: skip potentially time intensive gates+ "
                                              "[default: try top open all gates]");
        print_option("--re-estimate", "Re-estimate the difficulty after rotamers have been checked for "
                                      "collisions before applying the difficulty threshold. "
                                      "[default: false]");

        if (msg.empty()) exit(0);
        exit(1);
    }

    void print_option(const std::string &opt, const std::string &exp) const {
        std::vector<std::string> rows = wrap(exp, line_width - col - sep.length() - indent(1).length());
        for (size_t i = 0; i < rows.size(); i++) {
            if (i == 0) {
                std::cout << indent(1) << std::setw(col) << std::left << opt << sep << rows[i] << std::endl;
            } else {
                std::cout << indent(1) << std::setw(col) << " " << sep << rows[i] << std::endl;
            }
        }
    }

    void print_text(const std::string &text) const {
        std::vector<std::string> rows = wrap(text, line_width);
        for (const std::string &row: rows) { std::cout << row << std::endl; }
    }

    void write_pore_id_parameters() const {
        // create the log file or overwrite the old one
        std::ofstream file(pore_log);
        file.close();
        // parameters
        add_comment(pore_log, 0, "given and inferred command line parameters");
        add_entry(pore_log, 0, "parameters");
        add_comment(pore_log, 1, "input and output");
        add_entry(pore_log, 1, "PDB name", pdb_name);
        add_entry(pore_log, 1, "PDB path", pdb_path.string());
        add_entry(pore_log, 1, "output directory", out_dir.string());
        add_entry(pore_log, 1, "output name", output_name);
        add_comment(pore_log, 1, "PDB parsing");
        add_entry(pore_log, 1, "H-atoms", to_str(h_atom_tag));
        add_entry(pore_log, 1, "HETATM entries", to_str(hetero_atom_tag));
        add_entry(pore_log, 1, "skip non-standard amino acids in ATOM records", skip_non_standard_amino_acids);
        add_entry(pore_log, 1, "keep alternative atom locations", keep_alternative);
        add_comment(pore_log, 1, "run flags");
        add_entry(pore_log, 1, "output preparation", to_str(preparation));
        add_comment(pore_log, 1, "identification parameters");
        add_entry(pore_log, 1, "resolution", box_length);
        add_entry(pore_log, 1, "solvent radius", solvent_radius);
        add_entry(pore_log, 1, "probe radius", probe_radius);
        add_entry(pore_log, 1, "volume threshold", volume_threshold);
        add_entry(pore_log, 1, "computation mode", to_str(computation_mode));
        add_entry(pore_log, 1, "filter", to_str(pore_filter));
        // initialise log
        add_comment(pore_log, 0, "information about the pore ID run");
        add_entry(pore_log, 0, "log");
    }

    void write_axes_trace_parameters() const {
        // create the log file or overwrite the old one
        std::ofstream file(axis_log);
        file.close();
        // parameters
        add_comment(axis_log, 0, "given and inferred command line parameters");
        add_entry(axis_log, 0, "parameters");
        add_comment(axis_log, 1, "input and output");
        add_entry(axis_log, 1, "PDB name", pdb_name);
        add_entry(axis_log, 1, "output directory", axis_dir.string());
        add_entry(axis_log, 1, "load from pore ID", run_pore_id);
        add_entry(axis_log, 1, "load from single file", load_cluster_from_single_file);
        add_entry(axis_log, 1, "load from directory", load_cluster_from_directory);

        if (run_pore_id) {
            add_entry(axis_log, 1, "input path", "None");
        } else if (load_cluster_from_single_file) {
            add_entry(axis_log, 1, "input path", axis_single_input.string());
        } else {
            add_entry(axis_log, 1, "input path", axis_directory_input.string());
        }
        add_comment(axis_log, 1, "axes trace parameters");
        add_entry(axis_log, 1, "surface patch threshold", surface_patch_threshold);
        // initialise log
        add_comment(axis_log, 0, "information about the axes trace run");
        add_entry(axis_log, 0, "log");
    }

    void write_gate_open_parameters() const {
        // create the log file or overwrite the old one
        std::ofstream file(gate_log);
        file.close();
        // parameters
        add_comment(gate_log, 0, "given and inferred command line parameters");
        add_entry(gate_log, 0, "parameters");
        add_comment(gate_log, 1, "input and output");
        add_entry(gate_log, 1, "PDB name", pdb_name);
        add_entry(gate_log, 1, "PDB path", pdb_path.string());
        add_entry(gate_log, 1, "output directory", gate_dir.string());
        add_comment(gate_log, 1, "PDB parsing");
        add_entry(pore_log, 1, "H-atoms", to_str(h_atom_tag));
        add_entry(pore_log, 1, "HETATM entries", to_str(hetero_atom_tag));
        add_entry(pore_log, 1, "skip non-standard amino acids in ATOM records", skip_non_standard_amino_acids);
        add_entry(gate_log, 1, "keep alternative atom locations", keep_alternative);
        add_comment(gate_log, 1, "gate input");
        add_entry(gate_log, 1, "load from pore ID", run_pore_id);
        add_entry(gate_log, 1, "load from single file", load_gate_from_single_file);
        add_entry(gate_log, 1, "load from directory", load_gates_from_directory);
        if (run_pore_id) {
            add_entry(gate_log, 1, "input path", "None");
        } else if (load_cluster_from_single_file) {
            add_entry(gate_log, 1, "input path", gate_single_input.string());
        } else {
            add_entry(gate_log, 1, "input path", gate_directory_input.string());
        }
        add_comment(gate_log, 1, "gate open parameters");
        add_entry(gate_log, 1, "perturb value", perturb_value);
        add_entry(gate_log, 1, "clash tolerance", clash_tolerance);
        add_entry(gate_log, 1, "gate difficulty threshold", to_str(difficulty_threshold));
        add_entry(gate_log, 1, "re-estimate gate difficulty", re_estimate);
        // initialise log
        add_comment(gate_log, 0, "information about the gate open run");
        add_entry(gate_log, 0, "log");
    }
};

#endif //PROPORES_SETTINGS_H
