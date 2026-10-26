"""Building of ``pw.x`` input decks from a template plus parameter overrides.

The generators used to copy a single bare-atom/cluster template file and patch
``nat`` / ``ntyp`` / ``starting_magnetization`` by regex.  The same regex work
is duplicated across :mod:`lco.qe_input` (shared helpers) and the thin CLI
wrappers :mod:`Adsorbed_QE_Pt` and :mod:`Adsorbed_QE_lico`.

A critical correctness fix lives here: the templates hard-code
``prefix = 'pt4-tetrahedron'`` and ``outdir = './tmp/'``.  Every generated
input inherited both, so several variants in the same directory wrote to the
same scratch -- clobbering one another's saved charge density.  That made
``startingpot='file'`` restarts dangerous (another variant's density could be
read silently) and made the convergence experiments non-reproducible.  We now
give **each variant a unique prefix and a private outdir**, derived from the
system, site and variant tag.
"""
import re

#: Elements understood by the campaign templates, with mass and pseudopotential.
PSEUDO_DB = {
    "Li": {"mass": 6.941, "upf": "li_pbe_v1.4.uspp.F.UPF"},
    "Co": {"mass": 58.933, "upf": "co_pbe_v1.2.uspp.F.UPF"},
    "O": {"mass": 15.999, "upf": "o_pbe_v1.2.uspp.F.UPF"},
    "Pt": {"mass": 195.084, "upf": "pt_pbe_v1.4.uspp.F.UPF"},
    "Ni": {"mass": 58.693, "upf": "ni_pbe_v1.4.uspp.F.UPF"},
    "Cu": {"mass": 63.546, "upf": "cu_pbe_v1.2.uspp.F.UPF"},
}

#: Cards that are regenerated wholesale and must be stripped from a template.
CARDS_TO_STRIP = ("ATOMIC_SPECIES", "K_POINTS", "HUBBARD", "CELL_PARAMETERS",
                  "ATOMIC_POSITIONS", "PREFIX", "OUTDIR")

#: Hubbard U block appended to every generated input.
HUBBARD_BLOCK = "HUBBARD (ortho-atomic)\nU Co-3d 3.2\n"

#: k-points card used by all slab surface calculations.
SURFACE_KPOINTS = "K_POINTS (automatic)\n4 3 1 0 0 0\n"

def build_atomic_species_block(elements):
    """Generate the ``ATOMIC_SPECIES`` block for the given element list."""
    lines = ["ATOMIC_SPECIES"]
    for elem in elements:
        if elem not in PSEUDO_DB:
            raise KeyError(f"Element '{elem}' missing in PSEUDO_DB configuration.")
        info = PSEUDO_DB[elem]
        lines.append(f"  {elem:<2s}  {info['mass']:7.3f}  {info['upf']}")
    return "\n".join(lines) + "\n"


def strip_old_cards(template_content):
    """Remove template cards that are regenerated wholesale.

    ``prefix``/``outdir`` are deliberately stripped so each generated input
    gets a unique, collision-free pair (``prefix``/``outdir`` are then injected
    by :func:`make_prefix_directive`).
    """
    cleaned = []
    skip = False
    for line in template_content.splitlines():
        stripped = line.strip()
        upper = stripped.upper()
        if any(upper.startswith(card) for card in CARDS_TO_STRIP):
            skip = True
            continue
        if skip and (stripped.startswith("&") or stripped.upper().startswith("OCCUPATIONS")):
            skip = False
        if not skip and stripped:
            cleaned.append(line)
    return "\n".join(cleaned).strip()


def make_prefix_directive(prefix, outdir):
    """Return the ``prefix = ...`` / ``outdir = ...`` lines for a variant."""
    return f"   prefix = '{prefix}'\n   outdir = '{outdir}',"


def variant_prefix(system, site, tag):
    """A unique, collision-free QE prefix for one experiment variant.

    Components: ``lco`` (project marker), the system tag (e.g. ``Pt1`` or
    ``Pt4_tetra``), the adsorption site, and an experiment tag (e.g. the
    template stem ``44.1_template``).
    """
    return f"lco_{system}_{site}_{tag}"
