#!/usr/bin/env python3
"""Analyze VASP calculation outputs (vasprun.xml), print summary, and export to Markdown or text."""

import argparse
from pathlib import Path
import numpy as np
from pymatgen.io.vasp.outputs import Vasprun

# Energy conversion constant
RY_TO_EV = 13.605693122994  # 1 Ry in eV


def ev_to_ry(energy_ev: float) -> float:
    """Convert energy from electronvolts to Rydberg."""
    return energy_ev / RY_TO_EV


def parse_vasp(vasprun_path: Path) -> dict:
    """Parse vasprun.xml and extract key relaxation and convergence data."""
    v = Vasprun(
        str(vasprun_path),
        parse_dos=False,
        parse_eigen=False,
        parse_potcar_file=False,
    )

    structure = v.final_structure

    # Forces and stress
    if v.ionic_steps:
        last_step = v.ionic_steps[-1]
        forces = np.asarray(last_step.get("forces", []))
        stress = last_step.get("stress")
    else:
        forces = np.array([])
        stress = None

    max_force = float(np.linalg.norm(forces, axis=1).max()) if len(forces) > 0 else 0.0
    max_stress = float(np.max(np.abs(stress))) if stress is not None else None

    # Final energy
    energy_ev = float(v.final_energy)
    energy_ry = ev_to_ry(energy_ev)

    # Formula units in the simulation cell
    reduced_formula, n_formula_units = (
        structure.composition.get_reduced_composition_and_factor()
    )

    energy_fu_ev = energy_ev / n_formula_units
    energy_fu_ry = ev_to_ry(energy_fu_ev)

    # Lattice parameters
    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles

    return {
        "vasprun_path": str(vasprun_path.resolve()),
        "converged_electronic": v.converged_electronic,
        "converged_ionic": v.converged_ionic,
        "converged": v.converged,
        "n_ionic_steps": len(v.ionic_steps),
        "energy_ev": energy_ev,
        "energy_ry": energy_ry,
        "reduced_formula": reduced_formula,
        "n_formula_units": n_formula_units,
        "energy_fu_ev": energy_fu_ev,
        "energy_fu_ry": energy_fu_ry,
        "a": a,
        "b": b,
        "c": c,
        "alpha": alpha,
        "beta": beta,
        "gamma": gamma,
        "volume": structure.volume,
        "max_force": max_force,
        "max_stress_kb": max_stress,
        "structure": structure,
    }


def format_text(data: dict) -> str:
    """Format extracted data as plain text."""
    lines = [
        "=== VASP BULK RELAXATION ANALYSIS ===",
        "",
        "--- Convergence ---",
        f"Electronic converged: {data['converged_electronic']}",
        f"Ionic converged:      {data['converged_ionic']}",
        f"Overall converged:    {data['converged']}",
        f"Ionic steps:          {data['n_ionic_steps']}",
        "",
        "--- Total Energy ---",
        f"Final energy (eV):    {data['energy_ev']:.8f}",
        f"Final energy (Ry):    {data['energy_ry']:.8f}",
        "",
        f"--- Energy per {data['reduced_formula']} ---",
        f"Formula units:        {data['n_formula_units']:g}",
        f"Energy (eV/f.u.):     {data['energy_fu_ev']:.8f}",
        f"Energy (Ry/f.u.):     {data['energy_fu_ry']:.8f}",
        "",
        "--- Optimized Structure ---",
        f"a (Å):                {data['a']:.8f}",
        f"b (Å):                {data['b']:.8f}",
        f"c (Å):                {data['c']:.8f}",
        f"Alpha (°):            {data['alpha']:.6f}",
        f"Beta (°):             {data['beta']:.6f}",
        f"Gamma (°):            {data['gamma']:.6f}",
        f"Volume (Å³):          {data['volume']:.8f}",
        "",
        "--- Forces & Stress ---",
        f"Maximum force (eV/Å): {data['max_force']:.8e}",
    ]
    if data["max_stress_kb"] is not None:
        lines.append(f"Max stress (kB):      {data['max_stress_kb']:.8f}")

    lines.extend([
        "",
        "=== ANALYSIS COMPLETE ===",
        "",
    ])
    return "\n".join(lines)


def format_markdown(data: dict) -> str:
    """Format extracted data as GitHub-Flavored Markdown."""
    lines = [
        "# VASP Bulk Relaxation Analysis",
        "",
        f"- **Calculation File:** `{data['vasprun_path']}`",
        f"- **Reduced Formula:** {data['reduced_formula']}",
        f"- **Formula Units in Cell:** {data['n_formula_units']:g}",
        "",
        "## Convergence",
        "",
        "| Metric | Status / Value |",
        "| :--- | :--- |",
        f"| Electronic converged | `{data['converged_electronic']}` |",
        f"| Ionic converged | `{data['converged_ionic']}` |",
        f"| Overall converged | `{data['converged']}` |",
        f"| Ionic steps | {data['n_ionic_steps']} |",
        "",
        "## Energy Summary",
        "",
        "| Quantity | Energy (eV) | Energy (Ry) |",
        "| :--- | :--- | :--- |",
        f"| Total Final Energy | {data['energy_ev']:.8f} | {data['energy_ry']:.8f} |",
        f"| Energy per {data['reduced_formula']} | {data['energy_fu_ev']:.8f} | {data['energy_fu_ry']:.8f} |",
        "",
        "## Optimized Lattice Parameters",
        "",
        "| Parameter | Value | Unit |",
        "| :--- | :--- | :--- |",
        f"| *a* | {data['a']:.8f} | Å |",
        f"| *b* | {data['b']:.8f} | Å |",
        f"| *c* | {data['c']:.8f} | Å |",
        f"| $\\alpha$ | {data['alpha']:.6f} | ° |",
        f"| $\\beta$ | {data['beta']:.6f} | ° |",
        f"| $\\gamma$ | {data['gamma']:.6f} | ° |",
        f"| Cell Volume | {data['volume']:.8f} | Å³ |",
        "",
        "## Residual Forces & Stress",
        "",
        "| Metric | Value | Unit |",
        "| :--- | :--- | :--- |",
        f"| Maximum Atomic Force | {data['max_force']:.8e} | eV/Å |",
    ]
    if data["max_stress_kb"] is not None:
        lines.append(f"| Maximum Residual Stress | {data['max_stress_kb']:.8f} | kB |")

    structure = data.get("structure")
    if structure is not None:
        lines.extend([
            "",
            "## Final Atomic Positions (Fractional)",
            "",
            "| # | Species | *x* | *y* | *z* |",
            "| :--- | :--- | :--- | :--- | :--- |",
        ])
        for idx, site in enumerate(structure):
            lines.append(
                f"| {idx + 1} | {site.specie.symbol} | "
                f"{site.frac_coords[0]:.6f} | {site.frac_coords[1]:.6f} | {site.frac_coords[2]:.6f} |"
            )

    lines.append("")
    return "\n".join(lines)


def analyze_vasp(
    vasprun_path: str = "vasprun.xml",
    output_path: str | None = None,
    output_format: str = "auto",
    quiet: bool = False,
) -> dict:
    """Analyze VASP calculation, print summary to console, and write MD or TXT file."""
    path = Path(vasprun_path)
    if not path.exists():
        # Fallback relative to script location
        script_dir_path = Path(__file__).resolve().parent / vasprun_path
        if script_dir_path.exists():
            path = script_dir_path
        else:
            raise FileNotFoundError(f"VASP output not found: {vasprun_path}")

    data = parse_vasp(path)

    # Format text for stdout and text file
    txt_content = format_text(data)
    md_content = format_markdown(data)

    if not quiet:
        print(txt_content)

    # Determine output format and paths to save
    calc_dir = path.parent
    files_to_write = []

    if output_path is not None:
        out_file = Path(output_path)
        ext = out_file.suffix.lower()
        if output_format == "auto":
            fmt = "txt" if ext == ".txt" else "md"
        else:
            fmt = output_format

        content = md_content if fmt == "md" else txt_content
        files_to_write.append((out_file, content))
    else:
        # Default behavior:
        # If output_format is 'auto' or 'md' -> vasp_analysis.md
        # If 'txt' -> vasp_analysis.txt
        # If 'both' -> both vasp_analysis.md and vasp_analysis.txt
        if output_format in ("auto", "md"):
            files_to_write.append((calc_dir / "vasp_analysis.md", md_content))
        elif output_format == "txt":
            files_to_write.append((calc_dir / "vasp_analysis.txt", txt_content))
        elif output_format == "both":
            files_to_write.append((calc_dir / "vasp_analysis.md", md_content))
            files_to_write.append((calc_dir / "vasp_analysis.txt", txt_content))

    for target_file, content in files_to_write:
        target_file.parent.mkdir(parents=True, exist_ok=True)
        target_file.write_text(content, encoding="utf-8")
        if not quiet:
            print(f"\n[Saved report: {target_file}]")

    return data


def main():
    parser = argparse.ArgumentParser(
        description="Analyze VASP calculation (vasprun.xml), print summary, and export MD/TXT report."
    )
    parser.add_argument(
        "vasprun",
        nargs="?",
        default="vasprun.xml",
        help="Path to vasprun.xml (default: vasprun.xml)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Custom output file path (e.g. analysis.md or analysis.txt). Defaults to vasp_analysis.md in calculation dir.",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["auto", "md", "txt", "both"],
        default="auto",
        help="Output report format: 'md' (Markdown), 'txt' (plain text), 'both', or 'auto' (default: auto).",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress terminal printing (only write output file).",
    )

    args = parser.parse_args()
    analyze_vasp(
        vasprun_path=args.vasprun,
        output_path=args.output,
        output_format=args.format,
        quiet=args.quiet,
    )


if __name__ == "__main__":
    main()
