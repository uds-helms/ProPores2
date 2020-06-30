/*
 * Copyright (C) 2020 Markus Hollander (markus.hollander@bioinformatik.uni-saarland.de). All rights reserved.
 *
 * This file is part of PROPORES 2.0. PROPORES was originally developed by Po-Hsien Lee (2011) in Perl under the terms
 * of the GNU General Public License version 3 or higher. The approach was published as: * Lee, PH, Helms, V (2012). Identifying continuous pores in protein structures with PROPORES by computational
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

#include <string>
#include <vector>
#include "enums.h"

std::string to_str(const RecordType record) {
    switch (record) {
        case ATOM_RECORD: return "ATOM";
        case HETERO_RECORD: return "HETATM";
        default: return "INVALID_RECORD_TYPE";
    }
}

RecordType to_record_type(const std::string &str) {
    if (str == "ATOM") return ATOM_RECORD;
    if (str == "HETATM") return HETERO_RECORD;
    return INVALID_RECORD;
}

std::string to_str(const AtomType atom) {
    switch (atom) {
        // protein atoms
        case C: return "C";
        case H: return "H";
        case N: return "N";
        case O: return "O";
        case S: return "S";
        // hetero atoms
        case Ac: return "AC";
        case Ag: return "AG";
        case Al: return "AL";
        case Am: return "AM";
        case Ar: return "AR";
        case As: return "AS";
        case At: return "AT";
        case Au: return "AU";
        case B: return "B";
        case Ba: return "BA";
        case Be: return "BE";
        case Bh: return "BH";
        case Bi: return "BI";
        case Bk: return "BK";
        case Br: return "BR";
        case Ca: return "CA";
        case Cd: return "CD";
        case Ce: return "CE";
        case Cf: return "CF";
        case Cl: return "CL";
        case Cm: return "CM";
        case Cn: return "CN";
        case Co: return "CO";
        case Cr: return "CR";
        case Cs: return "CS";
        case Cu: return "CU";
        case Db: return "DB";
        case Ds: return "DS";
        case Dy: return "DY";
        case Er: return "ER";
        case Es: return "ES";
        case Eu: return "EU";
        case F: return "F";
        case Fe: return "FE";
        case Fl: return "FL";
        case Fm: return "FM";
        case Fr: return "FR";
        case Ga: return "GA";
        case Gd: return "GD";
        case Ge: return "GE";
        case He: return "HE";
        case Hf: return "HF";
        case Hg: return "HG";
        case Ho: return "HO";
        case Hs: return "HS";
        case I: return "I";
        case In: return "IN";
        case Ir: return "IR";
        case K: return "K";
        case Kr: return "KR";
        case La: return "LA";
        case Li: return "LI";
        case Lr: return "LR";
        case Lu: return "LU";
        case Lv: return "LV";
        case Mc: return "MC";
        case Md: return "MD";
        case Mg: return "MG";
        case Mn: return "MN";
        case Mo: return "MO";
        case Mt: return "MT";
        case Na: return "NA";
        case Nb: return "NB";
        case Nd: return "ND";
        case Ne: return "NE";
        case Nh: return "NH";
        case Ni: return "NI";
        case No: return "NO";
        case Np: return "NP";
        case Og: return "OG";
        case Os: return "OS";
        case P: return "P";
        case Pa: return "PA";
        case Pb: return "PB";
        case Pd: return "PD";
        case Pm: return "PM";
        case Po: return "PO";
        case Pr: return "PR";
        case Pt: return "PT";
        case Pu: return "PU";
        case Ra: return "RA";
        case Rb: return "RB";
        case Re: return "RE";
        case Rf: return "RF";
        case Rg: return "RG";
        case Rh: return "RH";
        case Rn: return "RN";
        case Ru: return "RU";
        case Sb: return "SB";
        case Sc: return "SC";
        case Se: return "SE";
        case Sg: return "SG";
        case Si: return "SI";
        case Sm: return "SM";
        case Sn: return "SN";
        case Sr: return "SR";
        case Ta: return "TA";
        case Tb: return "TB";
        case Tc: return "TC";
        case Te: return "TE";
        case Th: return "TH";
        case Ti: return "TI";
        case Tl: return "TL";
        case Tm: return "TM";
        case Ts: return "TS";
        case U: return "U";
        case V: return "V";
        case W: return "W";
        case Xe: return "XE";
        case Y: return "Y";
        case Yb: return "YB";
        case Zn: return "ZN";
        case Zr: return "ZR";
        default: return "INVALID_ATOM";
    }
}

AtomType to_atom_type(const std::string &str) {
    // protein atoms
    if (str == "C") return C;
    if (str == "H") return H;
    if (str == "N") return N;
    if (str == "O") return O;
    if (str == "S") return S;
    // hetero atoms
    if (str == "AC" || str == "Ac") return Ac;
    if (str == "AG" || str == "Ag") return Ag;
    if (str == "AL" || str == "Al") return Al;
    if (str == "AM" || str == "Am") return Am;
    if (str == "AR" || str == "Ar") return Ar;
    if (str == "AS" || str == "As") return As;
    if (str == "AT" || str == "At") return At;
    if (str == "AU" || str == "Au") return Au;
    if (str == "B") return B;
    if (str == "BA" || str == "Ba") return Ba;
    if (str == "BE" || str == "Be") return Be;
    if (str == "BH" || str == "Bh") return Bh;
    if (str == "BI" || str == "Bi") return Bi;
    if (str == "BK" || str == "Bk") return Bk;
    if (str == "BR" || str == "Br") return Br;
    if (str == "CA" || str == "Ca") return Ca;
    if (str == "CD" || str == "Cd") return Cd;
    if (str == "CE" || str == "Ce") return Ce;
    if (str == "CF" || str == "Cf") return Cf;
    if (str == "CL" || str == "Cl") return Cl;
    if (str == "CM" || str == "Cm") return Cm;
    if (str == "CN" || str == "Cn") return Cn;
    if (str == "CO" || str == "Co") return Co;
    if (str == "CR" || str == "Cr") return Cr;
    if (str == "CS" || str == "Cs") return Cs;
    if (str == "CU" || str == "Cu") return Cu;
    if (str == "DB" || str == "Db") return Db;
    if (str == "DS" || str == "Ds") return Ds;
    if (str == "DY" || str == "Dy") return Dy;
    if (str == "ER" || str == "Er") return Er;
    if (str == "ES" || str == "Es") return Es;
    if (str == "EU" || str == "Eu") return Eu;
    if (str == "F") return F;
    if (str == "FE" || str == "Fe") return Fe;
    if (str == "FL" || str == "Fl") return Fl;
    if (str == "FM" || str == "Fm") return Fm;
    if (str == "FR" || str == "Fr") return Fr;
    if (str == "GA" || str == "Ga") return Ga;
    if (str == "GD" || str == "Gd") return Gd;
    if (str == "GE" || str == "Ge") return Ge;
    if (str == "HE" || str == "He") return He;
    if (str == "HF" || str == "Hf") return Hf;
    if (str == "HG" || str == "Hg") return Hg;
    if (str == "HO" || str == "Ho") return Ho;
    if (str == "HS" || str == "Hs") return Hs;
    if (str == "I") return I;
    if (str == "IN" || str == "In") return In;
    if (str == "IR" || str == "Ir") return Ir;
    if (str == "K") return K;
    if (str == "KR" || str == "Kr") return Kr;
    if (str == "LA" || str == "La") return La;
    if (str == "LI" || str == "Li") return Li;
    if (str == "LR" || str == "Lr") return Lr;
    if (str == "LU" || str == "Lu") return Lu;
    if (str == "LV" || str == "Lv") return Lv;
    if (str == "MC" || str == "Mc") return Mc;
    if (str == "MD" || str == "Md") return Md;
    if (str == "MG" || str == "Mg") return Mg;
    if (str == "MN" || str == "Mn") return Mn;
    if (str == "MO" || str == "Mo") return Mo;
    if (str == "MT" || str == "Mt") return Mt;
    if (str == "NA" || str == "Na") return Na;
    if (str == "NB" || str == "Nb") return Nb;
    if (str == "ND" || str == "Nd") return Nd;
    if (str == "NE" || str == "Ne") return Ne;
    if (str == "NH" || str == "Nh") return Nh;
    if (str == "NI" || str == "Ni") return Ni;
    if (str == "NO" || str == "No") return No;
    if (str == "NP" || str == "Np") return Np;
    if (str == "OG" || str == "Og") return Og;
    if (str == "OS" || str == "Os") return Os;
    if (str == "P") return P;
    if (str == "PA" || str == "Pa") return Pa;
    if (str == "PB" || str == "Pb") return Pb;
    if (str == "PD" || str == "Pd") return Pd;
    if (str == "PM" || str == "Pm") return Pm;
    if (str == "PO" || str == "Po") return Po;
    if (str == "PR" || str == "Pr") return Pr;
    if (str == "PT" || str == "Pt") return Pt;
    if (str == "PU" || str == "Pu") return Pu;
    if (str == "RA" || str == "Ra") return Ra;
    if (str == "RB" || str == "Rb") return Rb;
    if (str == "RE" || str == "Re") return Re;
    if (str == "RF" || str == "Rf") return Rf;
    if (str == "RG" || str == "Rg") return Rg;
    if (str == "RH" || str == "Rh") return Rh;
    if (str == "RN" || str == "Rn") return Rn;
    if (str == "RU" || str == "Ru") return Ru;
    if (str == "SB" || str == "Sb") return Sb;
    if (str == "SC" || str == "Sc") return Sc;
    if (str == "SE" || str == "Se") return Se;
    if (str == "SG" || str == "Sg") return Sg;
    if (str == "SI" || str == "Si") return Si;
    if (str == "SM" || str == "Sm") return Sm;
    if (str == "SN" || str == "Sn") return Sn;
    if (str == "SR" || str == "Sr") return Sr;
    if (str == "TA" || str == "Ta") return Ta;
    if (str == "TB" || str == "Tb") return Tb;
    if (str == "TC" || str == "Tc") return Tc;
    if (str == "TE" || str == "Te") return Te;
    if (str == "TH" || str == "Th") return Th;
    if (str == "TI" || str == "Ti") return Ti;
    if (str == "TL" || str == "Tl") return Tl;
    if (str == "TM" || str == "Tm") return Tm;
    if (str == "TS" || str == "Ts") return Ts;
    if (str == "U") return U;
    if (str == "V") return V;
    if (str == "W") return W;
    if (str == "XE" || str == "Xe") return Xe;
    if (str == "Y") return Y;
    if (str == "YB" || str == "Yb") return Yb;
    if (str == "ZN" || str == "Zn") return Zn;
    if (str == "ZR" || str == "Zr") return Zr;
    return INVALID_ATOM;
}

AtomType to_atom_type(const char &str) {
    return to_atom_type(std::string(1, str));
}

std::string to_str(const ResidueType residue) {
    switch (residue) {
        case ALA: return "ALA";
        case ARG: return "ARG";
        case ASN: return "ASN";
        case ASP: return "ASP";
        case CYS: return "CYS";
        case GLN: return "GLN";
        case GLU: return "GLU";
        case GLY: return "GLY";
        case HIS: return "HIS";
        case ILE: return "ILE";
        case LEU: return "LEU";
        case LYS: return "LYS";
        case MET: return "MET";
        case PHE: return "PHE";
        case PRO: return "PRO";
        case SER: return "SER";
        case THR: return "THR";
        case TRP: return "TRP";
        case TYR: return "TYR";
        case VAL: return "VAL";
        default: return "INVALID_RESIDUE";
    }
}

ResidueType to_residue_type(const std::string &str) {
    if (str == "ALA") return ALA;
    if (str == "ARG") return ARG;
    if (str == "ASN") return ASN;
    if (str == "ASP") return ASP;
    if (str == "CYS") return CYS;
    if (str == "GLN") return GLN;
    if (str == "GLU") return GLU;
    if (str == "GLY") return GLY;
    if (str == "HIS") return HIS;
    if (str == "ILE") return ILE;
    if (str == "LEU") return LEU;
    if (str == "LYS") return LYS;
    if (str == "MET") return MET;
    if (str == "PHE") return PHE;
    if (str == "PRO") return PRO;
    if (str == "SER") return SER;
    if (str == "THR") return THR;
    if (str == "TRP") return TRP;
    if (str == "TYR") return TYR;
    if (str == "VAL") return VAL;
    return INVALID_RESIDUE;
}

std::string to_str(const BoxState state) {
    switch (state) {
        case OCCUPIED: return "OCCUPIED";
        case BACKGROUND: return "BACKGROUND";
        case POTENTIAL_PORE: return "POTENTIAL_PORE";
        case NUCLEUS: return "NUCLEUS";
        case IN_CLUSTER: return "IN_CLUSTER";
        case ON_CLUSTER_INSIDE: return "ON_CLUSTER_INSIDE";
        case TOUCHES_BACKGROUND: return "TOUCHES_BACKGROUND";
        case TOUCHES_PROTEIN: return "TOUCHES_PROTEIN";
        case TOUCHES_OTHER_CLUSTER: return "TOUCHES_OTHER_CLUSTER";
        case TOUCHES_POTENTIAL_PORE: return "TOUCHES_POTENTIAL_PORE";
        default: return "UNCLASSIFIED";
    }
}

BoxState to_box_state(const std::string &str) {
    if (str == "OCCUPIED") return OCCUPIED;
    if (str == "BACKGROUND") return BACKGROUND;
    if (str == "POTENTIAL_PORE") return POTENTIAL_PORE;
    if (str == "NUCLEUS") return NUCLEUS;
    if (str == "IN_CLUSTER") return IN_CLUSTER;
    if (str == "ON_CLUSTER_INSIDE") return ON_CLUSTER_INSIDE;
    if (str == "TOUCHES_BACKGROUND") return TOUCHES_BACKGROUND;
    if (str == "TOUCHES_PROTEIN") return TOUCHES_PROTEIN;
    if (str == "TOUCHES_OTHER_CLUSTER") return TOUCHES_OTHER_CLUSTER;
    if (str == "TOUCHES_POTENTIAL_PORE") return TOUCHES_POTENTIAL_PORE;
    return UNCLASSIFIED;
}

std::string to_str(HAtomTag tag) {
    switch (tag) {
        case KEEP_ALL_H_ATOMS: return "KEEP_ALL_H_ATOMS";
        case REMOVE_ALL_H_ATOMS: return "REMOVE_ALL_H_ATOMS";
        case REMOVE_ONLY_PROTEIN_H_ATOMS: return "REMOVE_ONLY_PROTEIN_H_ATOMS";
        case REMOVE_ONLY_HETERO_H_ATOMS: return "REMOVE_ONLY_HETERO_H_ATOMS";
        default: return "INVALID_H_ATOM_TAG";
    }
}

HAtomTag to_h_atom_tag(const std::string &str) {
    if (str == "0" || str == "KEEP_ALL_H_ATOMS") return KEEP_ALL_H_ATOMS;
    if (str == "1" || str == "REMOVE_ALL_H_ATOMS") return REMOVE_ALL_H_ATOMS;
    if (str == "2" || str == "REMOVE_ONLY_PROTEIN_H_ATOMS") return REMOVE_ONLY_PROTEIN_H_ATOMS;
    if (str == "3" || str == "REMOVE_ONLY_HETERO_H_ATOMS") return REMOVE_ONLY_HETERO_H_ATOMS;
    return INVALID_H_ATOM_TAG;
}

std::string to_str(HeteroAtomTag tag) {
    switch (tag) {
        case KEEP_ALL_HETERO_ATOMS: return "KEEP_ALL_HETERO_ATOMS";
        case REMOVE_ALL_HETERO_ATOMS: return "REMOVE_ALL_HETERO_ATOMS";
        case REMOVE_HETERO_ATOMS_EXCEPT_DUMMY: return "REMOVE_HETERO_ATOMS_EXCEPT_DUMMY";
        case REMOVE_ONLY_DUMMY_HETERO_ATOMS: return "REMOVE_ONLY_DUMMY_HETERO_ATOMS";
        default: return "INVALID_HETERO_ATOM_TAG";
    }
}

HeteroAtomTag to_hetero_atom_tag(const std::string &str) {
    if (str == "0" || str == "KEEP_ALL_HETERO_ATOMS") return KEEP_ALL_HETERO_ATOMS;
    if (str == "1" || str == "REMOVE_ALL_HETERO_ATOMS") return REMOVE_ALL_HETERO_ATOMS;
    if (str == "2" || str == "REMOVE_HETERO_ATOMS_EXCEPT_DUMMY") return REMOVE_HETERO_ATOMS_EXCEPT_DUMMY;
    if (str == "3" || str == "REMOVE_ONLY_DUMMY_HETERO_ATOMS") return REMOVE_ONLY_DUMMY_HETERO_ATOMS;
    return INVALID_HETERO_ATOM_TAG;
}

std::string to_str(const GateDifficulty difficulty) {
    switch (difficulty) {
        case EASY: return "EASY";
        case MEDIUM: return "MEDIUM";
        case HARD: return "HARD";
        default: return "DIFFICULTY_ERROR";
    }
}

GateDifficulty to_gate_difficulty(const std::string &str) {
    if (str == "EASY" || str == "0") return EASY;
    if (str == "MEDIUM" || str == "1") return MEDIUM;
    if (str == "HARD" || str == "2") return HARD;
    return DIFFICULTY_ERROR;
}

std::string to_str(ComputationMode tag) {
    switch (tag) {
        case AUTODETECT: return "AUTODETECT";
        case RAY_TRACE: return "RAY_TRACE";
        case STANDALONE: return "STANDALONE";
        default: return "INVALID_COMPUTATION_MODE";
    }
}

ComputationMode  to_computation_mode(const std::string &str) {
    if (str == "AUTODETECT" || str == "0") return AUTODETECT;
    if (str == "RAY_TRACE" || str == "1") return RAY_TRACE;
    if (str == "STANDALONE" || str == "2") return STANDALONE;
    return INVALID_COMPUTATION_MODE;
}

std::string to_str(RemoveTag tag) {
    switch (tag) {
        case REMOVE_NOTHING: return "REMOVE_NOTHING";
        case REMOVE_PORES: return "REMOVE_PORES";
        case REMOVE_CAVITIES: return "REMOVE_CAVITIES";
        default: return "INVALID_REMOVE_TAG";
    }
}

RemoveTag to_remove_tag(const std::string &str) {
    if (str == "REMOVE_NOTHING" || str == "0") return REMOVE_NOTHING;
    if (str == "REMOVE_PORES" || str == "1") return REMOVE_PORES;
    if (str == "REMOVE_CAVITIES" || str == "2") return REMOVE_CAVITIES;
    return INVALID_REMOVE_TAG;
}

std::string to_str(PreparationTag tag) {
    switch (tag) {
        case AXIS_AND_GATE: return "AXIS_AND_GATE";
        case ONLY_AXIS: return "ONLY_AXIS";
        case ONLY_GATE: return "ONLY_GATE";
        case NO_PREPARATION: return "NO_PREPARATION";
        default: return "INVALID_PREPARATION_TAG";
    }
}

PreparationTag to_preparation_tag(const std::string &str) {
    if (str == "AXIS_AND_GATE" || str == "0") return AXIS_AND_GATE;
    if (str == "ONLY_AXIS" || str == "1") return ONLY_AXIS;
    if (str == "ONLY_GATE" || str == "2") return ONLY_GATE;
    if (str == "NO_PREPARATION" || str == "3") return NO_PREPARATION;
    return INVALID_PREPARATION_TAG;
}

// map residue types to backbone definitions
std::vector<std::string> residue_definition(const ResidueType type) {
    if (type == ARG) return {"N", "CA", "CB", "CG", "CD", "NE", "CZ"};
    if (type == ASN) return {"N", "CA", "CB", "CG", "OD1"};
    if (type == ASP) return {"N", "CA", "CB", "CG", "OD1"};
    if (type == CYS) return {"N", "CA", "CB", "SG"};
    if (type == GLN) return {"N", "CA", "CB", "CG", "CD", "OE1"};
    if (type == GLU) return {"N", "CA", "CB", "CG", "CD", "OE1"};
    if (type == HIS) return {"N", "CA", "CB", "CG", "ND1"};
    if (type == ILE) return {"N", "CA", "CB", "CG1", "CD1"};
    if (type == LEU) return {"N", "CA", "CB", "CG", "CD1"};
    if (type == LYS) return {"N", "CA", "CB", "CG", "CD", "CE", "NZ"};
    if (type == MET) return {"N", "CA", "CB", "CG", "SD", "CE"};
    if (type == PHE) return {"N", "CA", "CB", "CG", "CD1"};
    if (type == PRO) return {"N", "CA", "CB", "CG", "CD"};
    if (type == SER) return {"N", "CA", "CB", "OG"};
    if (type == THR) return {"N", "CA", "CB", "OG1"};
    if (type == TRP) return {"N", "CA", "CB", "CG", "CD1"};
    if (type == TYR) return {"N", "CA", "CB", "CG", "CD1"};
    if (type == VAL) return {"N", "CA", "CB", "CG1"};
    return {};
}