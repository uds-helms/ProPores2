#ifndef PROPORES_PORE_ID_H
#define PROPORES_PORE_ID_H

#include <vector>
#include "atom.h"
#include "grid.h"
#include "gate.h"
#include "vector.h"

// scan for empty grid boxes involved in PSP events on the given axis
std::vector<uint8_t> psp_scan(ProteinGrid &grid);

// check if the connecting vector between atom pairs is interrupted by occupied grid boxes
void collision_detection(ProteinGrid &grid);

// form a cylinder around the connecting vector between atom pairs
void cylinder_trace(ProteinGrid &grid, const std::shared_ptr<AtomPair> &atom_pair, double r_sq, int r_grid,
                    double threshold, double search);

// form a cylinder around the connecting vector between atom pairs
void cylinder_completion(ProteinGrid &grid);

// compute which atom pair cylinders are perpendicular to each other
void perpendicularity_cylinder(ProteinGrid &grid);

// compute which atom pairs are perpendicular to each other without computing PSP events or cylinders
void perpendicularity_standalone(ProteinGrid &grid);

// mark grid boxes as part of the background if they are accessible from the solvent without passing through occupied
// boxes or boxes that are (potentially) part of a pore
void seal_gaps(ProteinGrid &grid);

// remove shallow pore regions on the protein surface
void remove_shallow_regions(ProteinGrid &grid, double shallow_radius);

// determine the largest ball radius that can be centred in each empty grid box
void identify_pore_nuclei(ProteinGrid &grid);

// build initial clusters of connected potential pore boxes
void initial_pore_assembly(ProteinGrid &grid, double tolerance);

// after the clusters are established, add the surrounding potential pore/cavity boxes to the nearest cluster
void complete_pores(ProteinGrid &grid);

// for each empty grid box involved in a pore/cavity, compute if it's at the protein surface, in touch with
// residues, in touch with another cluster or in the cluster centre
void label_pore_boxes(ProteinGrid &grid);

// classify clusters as pores or as cavities
void classify_clusters(ProteinGrid &grid);

// for each pore/cavity cluster add the surrounding residues
void lining_residues(ProteinGrid &grid);

// create files for pore/cavity boxes and lining residues
void output_pores(ProteinGrid &grid, const Settings &settings);

// main pore/cavity identification routine
void pore_ID(ProteinGrid &grid, const Settings &settings);

#endif //PROPORES_PORE_ID_H
