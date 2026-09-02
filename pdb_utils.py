import os
import sys
import urllib.request

from Bio.PDB import PDBParser


def download_pdb(pdbcode, datadir):
    pdbcode = pdbcode.upper()
    filename = pdbcode + ".pdb"

    url = f"https://files.rcsb.org/download/{filename}"
    filepath = os.path.join(datadir, filename)

    try:
        urllib.request.urlretrieve(url, filepath)
        print(f"Fișier descărcat: {filepath}")
        return filepath

    except Exception as err:
        print(f"Eroare la descărcare: {err}", file=sys.stderr)
        return None



def parse_pdb(filename):
    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("protein", filename)

    return structure



def get_pdb_text(filename):
    with open(filename, "r") as file:
        return file.read()


