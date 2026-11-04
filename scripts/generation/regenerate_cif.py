"""Rebuild ``*_latest.cif`` structure files from Quantum ESPRESSO outputs.

Walks ``root_dir`` for ``OUTPUT_*.out`` and rewrites the CIF found in the same
directory with the last ionic step. Run from the repository root, because the
default ``root_dir`` is relative.
"""
from pathlib import Path
from ase.io import read, write

def regenerate_cifs(root_dir="05_metals"):
    root = Path(root_dir)
    out_files = list(root.rglob("OUTPUT_*.out"))
    
    if not out_files:
        print("No OUTPUT_*.out files found.")
        return

    for out_path in out_files:
        folder = out_path.parent
        # Find existing .cif file in the immediate directory to preserve exact naming
        cif_files = [f for f in folder.glob("*.cif") if f.is_file()]
        
        if cif_files:
            cif_path = cif_files[0]
        else:
            cif_path = folder / f"{out_path.stem.replace('OUTPUT_', '')}.cif"

        try:
            # Read the final relaxed ionic step from the QE output file
            atoms = read(out_path, index="-1", format="espresso-out")
            write(cif_path, atoms, format="cif")
            print(f"[SUCCESS] Updated: {cif_path}")

        except Exception as e:
            # Fallback to schema XML if the .out text stream is truncated or broken
            xml_files = list(folder.rglob("*.xml"))
            if xml_files:
                try:
                    atoms = read(xml_files[0], format="espresso-xml")
                    write(cif_path, atoms, format="cif")
                    print(f"[RECOVERED VIA XML] Updated: {cif_path}")
                    continue
                except Exception:
                    pass

            print(f"[FAILED] Could not process {out_path}: {e}")

if __name__ == "__main__":
    regenerate_cifs()