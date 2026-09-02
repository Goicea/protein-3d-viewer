import py3Dmol
from py3Dmol import view
from stmol import showmol


def create_protein_view(
    pdb_text,
    representation="cartoon",
    color_scheme="spectrum",
    width=1100,
    height=750
):

    view = py3Dmol.view(
        width=width,
        height=height
    )

    # Adaugă structura
    view.addModel(
        pdb_text,
        "pdb"
    )

    # Alegerea reprezentării
    if representation == "cartoon":

        view.setStyle({
            "cartoon": {
                "color": color_scheme
            }
        })

    elif representation == "ball and stick":

        view.setStyle({
            "stick": {
                "colorscheme": color_scheme,
                "radius": 0.15
            },
            "sphere": {
                "colorscheme": color_scheme,
                "scale": 0.3
            }
        })

    elif representation == "sticks":

        view.setStyle({
            "stick": {
                "colorscheme": color_scheme
            }
        })

    elif representation == "spheres":

        view.setStyle({
            "sphere": {
                "colorscheme": color_scheme
            }
        })

    elif representation == "surface":

        view.setStyle({
            "cartoon": {
                "color": color_scheme
            }
        })

        view.addSurface(
            py3Dmol.VDW,
            {
                "opacity": 0.7
            }
        )

    # Centrează și afișează toată proteina
    view.zoomTo()

    return view