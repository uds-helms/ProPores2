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
    fs::path pdb_path;
    fs::path out_dir_path;
    fs::path pore_path;
    fs::path lining_path;
    fs::path axes_path;
    fs::path axes_preparation_path;
    fs::path gate_path;
    fs::path gate_preparation_path;
    fs::path axis_single_input_file;
    fs::path axis_directory_input_path;
    fs::path gate_single_input_file_path;
    fs::path gate_input_directory_path;
    fs::path rotamer_lib_path = fs::path("rotamer_lib");

    std::string output_name;

    // PARAMETERS
    // PORE-ID and AXIS-TRACE
    // grid dimensions and pore-ID parameters
    double box_length = 1.0;
    double solvent_radius = 1.2;
    double probe_radius = 1.4;
    // PORE ID
    // threshold for the minimum volume of a pore/cavity in Angstrom
    double volume_threshold = 50;
    // TODO
    // highest number of grid boxes before switching from cylinder trace to cylinder standalone
    size_t box_threshold = 1500000;
    // highest number og atoms before switching from cylinder trace to cylinder standalone
    size_t atom_pair_threshold = 30000;
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
    bool skip_H = false;
    // PORE-ID
    // whether to generate output files that can be used for axis-trace or gate-open
    bool run_axis_trace_preparation = false;
    bool run_gate_preparation = false;
    // whether to remove pores, cavities or keep both
    RemoveTag remove_tag = REMOVE_NOTHING;
    // whether to run cylinder trace, cylinder standalone or try to auto-detect the best approach
    CylinderTag cylinder_tag = AUTODETECT;
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
        // extract flags
        run_axis_trace_preparation = option_exists("--axis-preparation");
        run_gate_preparation = option_exists("--gate-preparation");
        keep_alternative = option_exists("--keep-alternative");
        re_estimate = option_exists("--re-estimate");
        bool help_short = option_exists("-h");
        bool help_long = option_exists("--help");
        bool only_cavities = option_exists("--only-cavities");
        bool only_pores = option_exists("--only-pores");
        bool standalone = option_exists("--cylinder-standalone");
        bool ray_trace = option_exists("--cylinder-ray-trace");
        bool skip_hard_gates = option_exists("--skip-hard-gates");
        bool skip_medium_gates = option_exists("--skip-medium-gates");
        // extract file paths and output specifications
        output_name = get_option("--name");
        std::string input_pdb = get_option("-i");
        std::string output_dir = get_option("-o");
        std::string axis_input_single = get_option("-ts");
        std::string axis_input_dir = get_option("-td");
        std::string gate_input_single = get_option("-gs");
        std::string gate_input_dir = get_option("-gd");

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
            out_dir_path = fs::absolute(std::filesystem::path(output_dir));
        } catch (...) { print_help("The given output directory (-o) is invalid."); }

        // extract numerical parameters
        set_option("-b", box_length, "grid box length (-b)");
        set_option("-s", solvent_radius, "solvent radius (-s)");
        set_option("-p", probe_radius, "probe radius (-p)");
        set_option("-v", volume_threshold, "volume threshold (-v)");
        set_option("-spt", surface_patch_threshold, "surface patch threshold (-spt)");
        set_option("-ct", clash_tolerance, "clash tolerance (-ct)");
        // check if the numerical parameters are within the expected value range
        if (box_length < 0) print_help("The option '-b <decimal>' was negative.");
        if (solvent_radius < 0) print_help("The option '-s <decimal>' was negative.");
        if (probe_radius < 0) print_help("The option '-p <decimal>' was negative.");
        if (volume_threshold < 0) print_help("The option '-v <decimal>' was negative.");
        if (surface_patch_threshold < 0) print_help("The option '-spt <decimal>' was negative.");
        if (clash_tolerance < 0) print_help("The option '-ct <decimal>' was negative.");

        // check if mutually exclusive options are used
        if (only_pores && only_cavities) {
            print_help("The options '--only-pores' and '--only-cavities' are mutually exclusive.");
        } else if (only_pores) {
            remove_tag = REMOVE_CAVITIES;
        } else if (only_cavities) { remove_tag = REMOVE_PORES; }

        if (ray_trace && standalone) {
            print_help("The options '--cylinder-ray-trace' and '--cylinder-standalone' are mutually exclusive.");
        } else if (ray_trace) {
            cylinder_tag = RAY_TRACE;
        } else if (standalone) { cylinder_tag = STANDALONE; }

        if (skip_hard_gates && skip_medium_gates) {
            print_help("The options '--skip-medium-gates' and '--skip-hard-gates' are mutually exclusive.");
        } else if (skip_hard_gates) {
            difficulty_threshold = HARD;
        } else if (skip_medium_gates) { difficulty_threshold = MEDIUM; }

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
                print_help("When pore-ID is not run, axis-trace input must be provided "
                           "(either '-ts <path>' or '-td <path>').");
            }
            // pore-ID is enabled and there was at least one custom axis-trace input
            if (run_pore_id && (load_cluster_from_directory || load_cluster_from_single_file)) {
                print_help("When pore-ID is run, do not provide axis determination input files ('-ts' or '-td'). "
                           "For running axis-trace with user defined input, do not run pore-ID.");
            }
            // at this point we know that if the single axis-trace file is provided, then pore-ID is not enabled
            // and there is no axis-trace directory given either, so we just have to check the file path
            if (load_cluster_from_single_file) {
                axis_single_input_file = fs::path(axis_input_single);
                if (!fs::is_regular_file(axis_single_input_file)) {
                    print_help("The axis trace input file path (-ts) does not point to an existing file. "
                               "Make sure to wrap file paths with spaces in quotation marks.");
                } else if (fs::is_empty(axis_single_input_file)) {
                    print_help("The axis trace input file path (-ts) points to an empty file.");
                }
            }
            // at this point we know that if the axis-trace directory is provided, then pore-ID is not enabled
            // and there is no single axis-trace file given either, so we just have to check the directory path
            if (load_cluster_from_directory) {
                axis_directory_input_path = fs::path(axis_input_dir);
                if (!fs::is_directory(axis_directory_input_path)) {
                    print_help("The axis trace input directory path (-td) does not point to an existing directory. "
                               "Make sure to wrap directory paths with spaces in quotation marks.");
                } else if (fs::is_empty(axis_directory_input_path)) {
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
                print_help("When pore-ID is not run, gate-open input must be provided "
                           "(either '-gs <path>' or '-gd <path>').");
            }
            // pore-ID is enabled and there was at least one custom gate-open input
            if (run_pore_id && (load_gates_from_directory || load_gate_from_single_file)) {
                print_help("When pore-ID is run, do not provide gate-open input files ('-gs' or '-gd'). "
                           "For running gate-open with user defined input, do not run pore-ID.");
            }
            // at this point we know that if the single gate-open file is provided, then pore-ID is not enabled
            // and there is no gate-open directory given either, so we just have to check the file path
            if (load_gate_from_single_file) {
                gate_single_input_file_path = fs::path(gate_input_single);
                if (!fs::is_regular_file(gate_single_input_file_path)) {
                    print_help("The gate-open input file path (-gs) does not point to an existing file. "
                               "Make sure to wrap file paths with spaces in quotation marks.");
                } else if (fs::is_empty(gate_single_input_file_path)) {
                    print_help("The gate-open input file path (-gs) points to an empty file.");
                }
            }
            // at this point we know that if the gate-open directory is provided, then pore-ID is not enabled
            // and there is no single gate-open file given either, so we just have to check the directory path
            if (load_gates_from_directory) {
                gate_input_directory_path = fs::path(gate_input_dir);
                if (!fs::is_directory(gate_input_directory_path)) {
                    print_help("The gate-open input directory path (-gd) does not point to an existing directory. "
                               "Make sure to wrap directory paths with spaces in quotation marks.");
                } else if (fs::is_empty(gate_input_directory_path)) {
                    print_help("The gate-open input directory path (-gd) points to an empty directory.");
                }
            }
        }

        // set the output name to the PDB name if none was given on the command line
        if (output_name.empty()) output_name = pdb_path.stem().string();
        // construct the output directory path
        out_dir_path /= fs::path(output_name);
        if (!output_name.empty() && output_name.back() != '_') output_name += '_';

        // generate output sub-directory paths
        pore_path = out_dir_path / fs::path("pores");
        lining_path = out_dir_path / fs::path("lining");
        axes_path = out_dir_path / fs::path("axes");
        gate_path = out_dir_path / fs::path("gate_open");
        axes_preparation_path = out_dir_path / fs::path("axis_trace_input");
        gate_preparation_path = out_dir_path / fs::path("gate_open_input");

        // try to delete the previous output depending on which program parts are enabled
        if (fs::exists(out_dir_path)) {
            // for pore-ID all previous output is deleted
            if (run_pore_id) {
                try {
                    fs::remove_all(out_dir_path);
                } catch (...) {
                    print_help("The content of the already existing output directory at (" + out_dir_path.string()
                               + ") could not be deleted. Check if it or the content is open in another program.");
                }
            } else {
                // if only axis-trace and/or gate-open are enabled, only delete previous axis-trace/gate-open output
                if (run_axis_trace && fs::exists(axes_path)) {
                    try {
                        fs::remove_all(axes_path);
                    } catch (...) {
                        print_help("The content of the already existing axis-trace directory at '"
                                   + axes_path.string() + "' could not be deleted. . Check if it or the content is "
                                                          "open in another program.");
                    }
                }
                if (run_gate_open && fs::exists(gate_path)) {
                    try {
                        fs::remove_all(gate_path);
                    } catch (...) {
                        print_help("The content of the already existing gate-open directory at '"
                                   + gate_path.string() + "' could not be deleted. . Check if it or the content is "
                                                          "open in another program.");
                    }
                }
            }
        }
        // create the main output directory
        fs::create_directories(out_dir_path);
        // create the sub-directories for enabled program parts
        if (run_pore_id) {
            fs::create_directories(pore_path);
            fs::create_directories(lining_path);
        }
        if (run_axis_trace) fs::create_directories(axes_path);
        if (run_gate_open) fs::create_directories(gate_path);
        if (run_axis_trace_preparation) fs::create_directories(axes_preparation_path);
        if (run_gate_preparation) fs::create_directories(gate_preparation_path);
        // if gate-open is enabled, check if the rotamer library is where it should be
        if (run_gate_open) {
            if (!fs::is_directory(rotamer_lib_path)) {
                print_help("The rotamer library could not be found at: " + rotamer_lib_path.string());
            } else if (fs::is_empty(rotamer_lib_path)) {
                print_help("The rotamer library at '" + rotamer_lib_path.string() + "' is empty.");
            }
        }

    }

    // generates the file path with all rotamer definitions for the given residue type
    [[nodiscard]] std::string rotamer_file_path(const std::string &residue_type) const {
        fs::path file_path = rotamer_lib_path / fs::path(residue_type + "_trimed.txt");
        if (fs::exists(file_path)) return file_path.string();
        return "";
    }

private:
    std::vector<std::string> args;
    const size_t col = 21;
    const std::string indent = "  ";
    const std::string sep = "  ";
    const size_t line_width = 90;

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
        std::cout << "PROPORES 2.0" << std::endl;
        print_text("Identification of protein pores, cavities and channels with options for axis determination and "
                   "opening connections between neighbouring pores.");
        std::cout << std::endl;

        if (!msg.empty()) {
            print_text("ERROR: " + msg + " See the following usage description.");
            std::cout << std::endl;
        }

        std::cout << "Usage:" << std::endl;
        std::cout << sep << "propores <run command(s)> -i <path> -o <path> [options...]" << std::endl << std::endl;

        std::cout << "Run commands (at least one, space separated):" << std::endl;
        print_option("pore-id", "Run pore/cavity identification.");
        print_option("axis-trace", "Run pore/cavity axis determination.");
        print_option("gate-open", "Run gate-open for neighbouring pores/cavities.");
        std::cout << std::endl;

        std::cout << "Required:" << std::endl;
        print_option("-i <path>", "Path to the input PDB file.");
        print_option("-o <path>", "Path to the output directory. A sub-directory will be created for the current run. "
                                  "If the sub-directory already exists and pore-ID is run, its content will get "
                                  "deleted!");
        std::cout << std::endl;

        std::cout << "General options:" << std::endl;
        print_option("-h, --help", "Show this help message.");
        print_option("--name <string>", "Name for the output sub-directory and  all output files. [default: PDB name]");
        print_option("--skip-H-atoms", "Ignore H-atoms in the input PDB file. [default: false]");
        print_option("--keep-alternative", "Load alternative locations of atoms in addition to the primary location. "
                                           "[default: false");

        std::cout << "Pore-ID options:" << std::endl;
        print_option("-b <decimal>", "Grid box length in Angstrom (= resolution). [default: 1.0]");
        print_option("-s <decimal>", "Solvent radius in Angstrom. [default: 1.0]");
        print_option("-p <decimal>", "Probe radius in Angstrom. [default: 1.0]");
        print_option("-v <decimal>", "Volume threshold for pores/cavities in cubic Angstrom. [default: 50]");
        print_option("--cylinder-ray-trace", "Faster for smaller proteins or larger grid box lengths. "
                                             "Increased RAM usage. [default: autodetect]");
        print_option("--cylinder-standalone", "Potentially faster for larger proteins or smaller grid box lengths. "
                                              "Reduced RAM usage. [default: autodetect]");
        print_option("--only-pores", "Only generate output and analyses for pores, not cavities. [default: both]");
        print_option("--only-cavities", "Only generate output and analyses for cavities. [default: both]");
        print_option("--axis-preparation", "Output files that can be used to trace the axis of individual "
                                           "pores/cavities. [default: false]");
        print_option("--gate-preparation", "Output files that can be used to open individual gates between "
                                           "neighbouring pores/cavities. [default: false]");

        std::cout << std::endl << "Axis-trace options:" << std::endl;
        print_option("-spt <number>", "Threshold for the minimum number of grid boxes in pore surface patches. "
                                      "[default: 30]");
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
                                      "rotamers. [default: 0.5]");
        print_option("-gs <path>", "Path to a file with information for opening the gate between a single pair of "
                                   "two pores/cavities. Only needed if gate-open is not run together with pore ID. "
                                   "[default: none]");
        print_option("-gd <path>", "Path to a directory with file(s) with information for opening the gate between a "
                                   "single pair of two pores/cavities. Only needed if gate-open is not run together "
                                   "with pore ID. [default: none]");
        print_option("--skip-hard-gates", "Skip very time intensive gates. [default: false]");
        print_option("--skip-medium-gates", "Skip potentially time intensive gates. This will automatically enable "
                                            "'--skip-hard-gates' as well. [default: false]");
        print_option("--re-estimate", "Re-estimate the difficulty after rotamers have been checked for "
                                      "collisions before applying the difficulty threshold. "
                                      "[default: false]");

        if (msg.empty()) exit(0);
        exit(1);
    }

    void print_option(const std::string &opt, const std::string &exp) const {
        std::vector<std::string> rows = wrap(exp, line_width - col - sep.length() - indent.length());
        for (size_t i = 0; i < rows.size(); i++) {
            if (i == 0) {
                std::cout << indent << std::setw(col) << std::left << opt << sep << rows[i] << std::endl;
            } else {
                std::cout << indent << std::setw(col) << " " << sep << rows[i] << std::endl;
            }
        }
    }

    void print_text(const std::string &text) const {
        std::vector<std::string> rows = wrap(text, line_width);
        for (const std::string &row: rows) { std::cout << row << std::endl; }
    }
};

#endif //PROPORES_SETTINGS_H
