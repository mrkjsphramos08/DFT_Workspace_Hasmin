"""Canonical adsorption-site definitions for the LCO cluster campaign.

``SITE_SPECS`` is the single source of truth. Each site has a ``kind``:
  - "atop":   above the surface atom of ``elements[0]`` closest to the slab center
  - "bridge": midpoint of the pair of surface atoms (``elements``) closest to the
              slab center; pairs farther apart than ``max_dist`` (Å) are ignored
              (``None`` = no cutoff)
  - "center": center of the xy cell

Dict order is the generation order, so new sites go at the end to keep
existing results reproducible.
"""
SITE_SPECS = {
    "li":       {"kind": "atop",   "elements": ("Li",)},
    "co":       {"kind": "atop",   "elements": ("Co",)},
    "lico":     {"kind": "bridge", "elements": ("Li", "Co"), "max_dist": None},
    "centered": {"kind": "center", "elements": ()},
    "oxygen":   {"kind": "atop",   "elements": ("O",)},
    "lili":     {"kind": "bridge", "elements": ("Li", "Li"), "max_dist": 4.5},
    "coco":     {"kind": "bridge", "elements": ("Co", "Co"), "max_dist": 4.5},
    # new sites
    "lio":      {"kind": "bridge", "elements": ("Li", "O"),  "max_dist": 3.0},
    "coo":      {"kind": "bridge", "elements": ("Co", "O"),  "max_dist": 3.0},
}

SITE_ORDER = tuple(SITE_SPECS)
VALID_SITES = set(SITE_ORDER)
