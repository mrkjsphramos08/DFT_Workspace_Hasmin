from pathlib import Path
from ase import Atoms
from ase.io import write

out_dir = Path("05_metals/Pt1")
out_dir.mkdir(parents=True, exist_ok=True)

pt = Atoms("Pt", positions=[[7.5, 7.5, 7.5]], cell=[15.0, 15.0, 15.0], pbc=True)
write(out_dir / "Pt1_latest.cif", pt)