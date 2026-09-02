import streamlit as st
import numpy as np
from Bio.PDB.vectors import refmat, calc_angle

def nr_atoms(structure):
    num_atoms = len(list(structure.get_atoms()))
    return num_atoms


def nr_residues(structure):
    num_residues = len(list(structure.get_residues()))
    return num_residues


def nr_chains(structure):
    num_chains = len(list(structure.get_chains()))
    return num_chains

def is_dna_residue(residue):
    return residue.resname in ["DA", "DT", "DG", "DC", "A", "T", "G","C"]

def is_hbond_donor(atom):
    residue = atom.get_parent()
    resname = residue.resname
    atom_name = atom.name

    # =========================
    # PROTEINE
    # =========================

    # Azotul din backbone
    if resname != "PRO" and atom_name == "N":
        return True

    # Donori din sidechain
    if resname == "SER" and atom_name == "OG":
        return True

    if resname == "THR" and atom_name == "OG1":
        return True

    if resname == "TYR" and atom_name == "OH":
        return True

    if resname == "CYS" and atom_name == "SG":
        return True

    # Donori N din sidechain
    if atom.element == "N":
        return True

    # =========================
    # ADN
    # =========================

    if resname in ["DA", "A"]:
        if atom_name == "N6":
            return True

    if resname in ["DT", "T"]:
        if atom_name == "N3":
            return True

    if resname in ["DG", "G"]:
        if atom_name in ["N1", "N2"]:
            return True

    if resname in ["DC", "C"]:
        if atom_name == "N4":
            return True

    return False


def is_hbond_acceptor(atom):
    residue = atom.get_parent()
    resname = residue.resname
    atom_name = atom.name

    # =========================
    # ADN
    # =========================

    if resname in ["DA", "A"]:
        if atom_name in ["N1", "N7"]:
            return True

    if resname in ["DT", "T"]:
        if atom_name in ["O2", "O4"]:
            return True

    if resname in ["DG", "G"]:
        if atom_name in ["O6", "N3", "N7"]:
            return True

    if resname in ["DC", "C"]:
        if atom_name in ["O2", "N3"]:
            return True

    # =========================
    # PROTEINE
    # =========================

    if atom.element == "O":
        return True

    if atom.element == "N":
        return True

    return False


def is_backbone_atom(atom):
    return atom.name in ["N", "CA", "C", "O"]


def constr_afis_h_bonds(distance, cutoff, hydrogen_bonds, donor, acceptor):
    if distance <= cutoff:
        hydrogen_bonds.append({
            "donor_atom": donor.serial_number,
            "donor_residue": donor.get_parent().resname,
            "donor_residue_number": donor.get_parent().id[1],
            "donor_chain": donor.get_parent().get_parent().id,

            "acceptor_atom": acceptor.serial_number,
            "acceptor_residue": acceptor.get_parent().resname,
            "acceptor_residue_number": acceptor.get_parent().id[1],
            "acceptor_chain": acceptor.get_parent().get_parent().id,

            "distance": float(distance)
        })



def constr_h_poz(atom_n, atom_ca, atom_c):
    c = atom_c.get_vector()
    n = atom_n.get_vector()
    ca = atom_ca.get_vector()
    h = n - ca
    c = c - n
    mirror = refmat(h, c)
    h = h.left_multiply(mirror)
    h_position = atom_n.get_vector() + h
    return h_position



def find_hydrogen_bonds(structure,distance_cutoff=3.5,angle_cutoff=120):
    donors = []
    acceptors = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:

                    if is_hbond_donor(atom):
                        donors.append(atom)

                    if is_hbond_acceptor(atom):
                        acceptors.append(atom)

    hydrogen_bonds_sec = []
    hydrogen_bonds_tert = []

    for donor in donors:
        for acceptor in acceptors:

            if donor == acceptor:
                continue

            if donor.get_parent() == acceptor.get_parent():
                continue

            distance = np.linalg.norm(donor.coord - acceptor.coord )
            if distance > distance_cutoff:
                continue

            angle = None
            if donor.name == "N":

                residue = donor.get_parent()
                chain = residue.get_parent()

                if "CA" in residue:
                    atom_ca = residue["CA"]
                    residues = list(chain)
                    try:
                        residue_index = residues.index(residue)
                    except ValueError:
                        continue
                    if residue_index > 0:
                        previous_residue = (residues[residue_index - 1])
                        if "C" in previous_residue:
                            atom_c = (previous_residue["C"])
                            h_position = constr_h_poz(donor,atom_ca,atom_c)
                            angle = np.degrees(calc_angle(donor.get_vector(),h_position, acceptor.get_vector() ) )

            if angle is None:
                continue

            if angle < angle_cutoff:
                continue

            donor_residue = donor.get_parent()
            acceptor_residue = acceptor.get_parent()

            donor_residue = donor.get_parent()
            acceptor_residue = acceptor.get_parent()

            bond = {
                "donor_atom": donor.serial_number,
                "acceptor_atom": acceptor.serial_number,

                "donor_residue": donor_residue.resname,
                "acceptor_residue": acceptor_residue.resname,

                "donor_residue_number": donor_residue.id[1],
                "acceptor_residue_number": acceptor_residue.id[1],

                "donor_chain": donor_residue.get_parent().id,
                "acceptor_chain": acceptor_residue.get_parent().id,

                "distance": float(distance),
                "angle": float(angle)
            }

            if (
                    is_dna_residue(donor_residue)
                    and is_dna_residue(acceptor_residue)
            ):

                # ADN → structură secundară
                hydrogen_bonds_sec.append(bond)

            elif (
                    is_backbone_atom(donor)
                    and is_backbone_atom(acceptor)
            ):

                # Proteină → backbone → structură secundară
                hydrogen_bonds_sec.append(bond)

            else:

                # Proteină → sidechain → structură terțiară
                hydrogen_bonds_tert.append(bond)

    return {
        "secondary": hydrogen_bonds_sec,
        "tertiary": hydrogen_bonds_tert
    }


def find_ss_bridges(structure, cutoff=3.0):
    ss_bridges = []
    s_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == "CYS":
                    s_atoms.append(residue["SG"])

    for i in range(len(s_atoms)):
        for j in range(i + 1, len(s_atoms)):
            atom1 = s_atoms[i]
            atom2 = s_atoms[j]

            distance = np.linalg.norm(
                atom1.coord - atom2.coord
            )

            if distance <= cutoff:
                residue1 = atom1.get_parent()
                residue2 = atom2.get_parent()

                chain1 = residue1.get_parent()
                chain2 = residue2.get_parent()

                ss_bridges.append({
                    "atom1": atom1.serial_number,
                    "atom2": atom2.serial_number,

                    "residue1": residue1.resname,
                    "residue2": residue2.resname,

                    "residue_number1": residue1.id[1],
                    "residue_number2": residue2.id[1],

                    "chain1": chain1.id,
                    "chain2": chain2.id,

                    "distance": float(distance)
                })

    return ss_bridges


def get_statistics(structure,distance_cutoff=3.5,angle_cutoff=120):
    return {
        "atoms": nr_atoms(structure),
        "residues": nr_residues(structure),
        "chains": nr_chains(structure),
        "hydrogen_bonds": find_hydrogen_bonds(structure,distance_cutoff=distance_cutoff,angle_cutoff=angle_cutoff),
        "ss_bridges": find_ss_bridges(structure)
    }
