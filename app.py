import streamlit as st
from stmol import showmol

from pdb_utils import download_pdb, parse_pdb, get_pdb_text
from statistics import get_statistics
from view import create_protein_view

st.set_page_config(layout="wide")

st.title("Protein 3D Viewer")

st.write(
    "Aplicație pentru vizualizarea structurilor proteinelor"
)

datadir = st.text_input(
    "Introdu locația pentru download"
)

pdb_id = st.text_input(
    "Introdu PDB ID"
)

if st.button("Download"):

    pdb_file = download_pdb(
        pdb_id,
        datadir
    )

    if pdb_file is None:
        st.error("Fișierul PDB nu a fost descărcat.")
        st.stop()

    pdb_text = get_pdb_text(
        pdb_file
    )

    st.session_state["pdb_file"] = pdb_file
    st.session_state["pdb_text"] = pdb_text
    st.session_state["pdb_id"] = pdb_id


if "pdb_text" in st.session_state:

    pdb_file = st.session_state["pdb_file"]
    pdb_text = st.session_state["pdb_text"]


    st.sidebar.subheader("Hydrogen bond criteria")

    distance_cutoff = st.sidebar.slider(
        "Maximum distance between heavy atoms (Å)",
        min_value=1.0,
        max_value=3.5,
        value=3.5,
        step=0.1
    )

    angle_cutoff = st.sidebar.slider(
        "Minimum D-H-A angle (°)",
        min_value=120,
        max_value=180,
        value=120,
        step=1
    )

    st.sidebar.write(
        f"Current criteria: "
        f"distance ≤ {distance_cutoff:.1f} Å, "
        f"angle ≥ {angle_cutoff}°"
    )

    structure = parse_pdb(pdb_file )

    stats = get_statistics(
        structure,
        distance_cutoff=distance_cutoff,
        angle_cutoff=angle_cutoff
    )


    st.sidebar.title("Protein statistics")

    st.write("Number of atoms:", stats["atoms"])

    st.write("Number of residues:", stats["residues"])

    st.write("Number of chains:", stats["chains"])

    st.sidebar.subheader("Hydrogen bonds")

    hydrogen_bonds = stats["hydrogen_bonds"]

    secondary_bonds = hydrogen_bonds["secondary"]
    tertiary_bonds = hydrogen_bonds["tertiary"]


    st.sidebar.write("Hydrogen bonds contributing to secondary structure:",len(secondary_bonds))

    for i, bond in enumerate(secondary_bonds):
        st.sidebar.write(
            f"{i + 1}. "
            f"Atom {bond['donor_atom']} "
            f"({bond['donor_residue']} "
            f"{bond['donor_residue_number']}, "
            f"chain {bond['donor_chain']}) "
            f"→ "
            f"Atom {bond['acceptor_atom']} "
            f"({bond['acceptor_residue']} "
            f"{bond['acceptor_residue_number']}, "
            f"chain {bond['acceptor_chain']}) "
            f"— distance: {bond['distance']:.2f} Å "
            f"— angle: {bond['angle']:.1f}°"
        )


    st.sidebar.write("Hydrogen bonds contributing to tertiary structure:", len(tertiary_bonds))

    for i, bond in enumerate(tertiary_bonds):
        st.sidebar.write(
            f"{i + 1}. "
            f"Atom {bond['donor_atom']} "
            f"({bond['donor_residue']} "
            f"{bond['donor_residue_number']}, "
            f"chain {bond['donor_chain']}) "
            f"→ "
            f"Atom {bond['acceptor_atom']} "
            f"({bond['acceptor_residue']} "
            f"{bond['acceptor_residue_number']}, "
            f"chain {bond['acceptor_chain']}) "
            f"— distance: {bond['distance']:.2f} Å "
            f"— angle: {bond['angle']:.1f}°"
        )

    st.sidebar.subheader("Disulphide bridges")
    ss_bonds = stats["ss_bridges"]
    st.sidebar.write("Number of disulphide bridges:",len(ss_bonds))

    for i, bond in enumerate(ss_bonds):
        st.sidebar.write(
            f"{i + 1}. "
            f"Atom {bond['atom1']} "
            f"({bond['residue1']} "
            f"{bond['residue_number1']}, "
            f"chain {bond['chain1']}) "
            f"↔ "
            f"Atom {bond['atom2']} "
            f"({bond['residue2']} "
            f"{bond['residue_number2']}, "
            f"chain {bond['chain2']}) "
            f"— {bond['distance']:.2f} Å"
        )


    st.subheader("3D Protein Viewer")

    representation = st.selectbox(
        "Representation",
        ["Cartoon",
         "Ball and stick",
         "Sticks",
         "Spheres",
         "Surface"]
    )

    color_scheme = st.selectbox(
        "Color scheme",
        ["Spectrum",
         "Chain",
         "Amino acid"]
    )

    representation_map = {
        "Cartoon": "cartoon",
        "Ball and stick": "ball and stick",
        "Sticks": "sticks",
        "Spheres": "spheres",
        "Surface": "surface"
    }


    color_map = {
        "Spectrum": "spectrum",
        "Chain": "chain",
        "Amino acid": "amino"
    }


    representation_value = (representation_map[representation])

    color_value = (color_map[color_scheme])

    view = create_protein_view(
        pdb_text=pdb_text,
        representation=representation_value,
        color_scheme=color_value,
        width=1500,
        height=750
    )

    showmol(
        view,
        width=1500,
        height=750
    )