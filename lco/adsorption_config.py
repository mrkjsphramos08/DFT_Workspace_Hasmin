"""Slab profiles for the LCO adsorption campaign.

Each profile names the relaxed slab CIF to build on and the output namespace its
adsorbed structures live under. ``get_slab_profile`` and ``slab_output_root`` are
the single source of truth used by the structure and job generators.
"""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

SLAB_PROFILES = {
    "3_layers": {
        "slab_path": PROJECT_ROOT / "04_slab" / "3_layers" / "3x2_3_layers.cif",
        "output_name": "3_layers",
    },
    "4_layers": {
        "slab_path": PROJECT_ROOT / "04_slab" / "4_layers" / "3x2_4_layers.cif",
        "output_name": "4_layers",
    },
    "3_layers_90": {
        "slab_path": PROJECT_ROOT / "04_slab" / "3_layers_90" / "3x2_3_layers_90.cif",
        "output_name": "3_layers_90",
    },
    "3_layers_new": {
        "slab_path": PROJECT_ROOT / "04_slab" / "3_layers_new" / "3x2_3_layers_new.cif",
        "output_name": "3_layers_new",
    },
}

DEFAULT_SLAB = "3_layers"
SLAB_CHOICES = tuple(SLAB_PROFILES)


def get_slab_profile(slab_name):
    try:
        return SLAB_PROFILES[slab_name]
    except KeyError as exc:
        choices = ", ".join(SLAB_PROFILES)
        raise ValueError(f"Unknown slab '{slab_name}'. Choose from: {choices}") from exc


def slab_output_root(metal_dir, slab_name):
    """Return the slab-specific adsorption directory for one cluster family."""
    profile = get_slab_profile(slab_name)
    return metal_dir / "Adsorbed" / profile["output_name"]
