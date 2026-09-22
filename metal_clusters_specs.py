#!/usr/bin/env python3
"""Export fresh summary CSV from the current bare-cluster QE outputs."""
import csv, math, re
from pathlib import Path

ROOT = Path('/home/mark/LCO/01_bulk')
METALS = ROOT / '05_metals'
OUT = ROOT / 'metal_cluster_summary.csv'

def rgx(text, pat, default=None):
    m = re.findall(pat, text, re.I | re.M)
    return m[-1] if m else default

def parse_out(out_path):
    text = out_path.read_text(errors='replace')
    final = rgx(text, r'Final energy\s*=\s*([-+0-9.Ee]+)\s+Ry')
    bang = rgx(text, r'^\s*!\s*total energy\s*=\s*([-+0-9.Ee]+)\s+Ry')
    fermi = rgx(text, r'the Fermi energy is\s*([-+0-9.Ee]+)\s+ev', default=None)
    mag = rgx(text, r'total magnetization\s*=\s*([-+0-9.Ee]+)\s+Bohr mag/cell', default=None)
    abs_mag = rgx(text, r'absolute magnetization\s*=\s*([-+0-9.Ee]+)\s+Bohr mag/cell', default=None)
    total_force = rgx(text, r'Total force\s*=\s*([-+0-9.Ee]+)')
    scf_corr = rgx(text, r'Total SCF correction\s*=\s*([-+0-9.Ee]+)')
    conv_iters = rgx(text, r'convergence has been achieved in\s+(\d+)\s+iterations', default=None)
    # Try the explicit "number of scf cycles = X" format first
    scf_cycles = rgx(text, r'number of scf cycles\s*=\s*(\d+)', default=None)
    # Fallback: parse from "bfgs converged in X scf cycles and Y bfgs steps"
    if not scf_cycles:
        scf_cycles = rgx(text, r'bfgs converged in\s+(\d+)\s+scf cycles', default=None)
    bfgs_steps = rgx(text, r'number of bfgs steps\s*=\s*(\d+)', default=None)
    if not bfgs_steps:
        bfgs_steps = rgx(text, r'bfgs converged in\s+\d+\s+scf cycles and\s+(\d+)\s+bfgs steps', default=None)
    # Match CPU/WALL time: supports both "5m 9.85s" and "1h 5m 9.85s" formats
    cpu = rgx(text, r'PWSCF\s*:.*?(\d+h\s*\d+m\s*[\d.]*s)\s*CPU')
    if not cpu:
        cpu = rgx(text, r'PWSCF\s*:.*?(\d+m\s*[\d.]*s)\s*CPU')
    wall = rgx(text, r'PWSCF\s*:.*?(\d+h\s*\d+m\s*[\d.]*s)\s*WALL')
    if not wall:
        wall = rgx(text, r'PWSCF\s*:.*?(\d+m\s*[\d.]*s)\s*WALL')
    return {
        'final_energy': final, 'scf_energy': bang, 'fermi': fermi,
        'mag': mag, 'abs_mag': abs_mag, 'total_force': total_force,
        'scf_correction': scf_corr, 'scf_iters': conv_iters,
        'scf_cycles': scf_cycles, 'bfgs_steps': bfgs_steps,
        'cpu': cpu, 'wall': wall,
    }

def parse_positions(text):
    begin = text.find('Begin final coordinates')
    end = text.find('End final coordinates', begin)
    block = text[begin:end] if begin != -1 and end != -1 else ''
    atoms = []
    for line in block.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[0][0].isalpha():
            atoms.append({'symbol': parts[0], 'x': float(parts[1]), 'y': float(parts[2]), 'z': float(parts[3])})
    return atoms

def bonds(atoms):
    dists = []
    for i in range(len(atoms)):
        for j in range(i+1, len(atoms)):
            dx = atoms[i]['x'] - atoms[j]['x']
            dy = atoms[i]['y'] - atoms[j]['y']
            dz = atoms[i]['z'] - atoms[j]['z']
            dists.append(math.sqrt(dx*dx + dy*dy + dz*dz))
    dists.sort()
    out = {'min': dists[0], 'max': dists[-1], 'avg': sum(dists)/len(dists)}
    out['pairs'] = dists
    return out

def classify(atoms):
    b = bonds(atoms)
    edge_avg = sum(b['pairs'][:3]) / 3
    off = b['pairs'][3:]
    if off and max(off) > 1.5 * edge_avg:
        return 'Rhombus'
    return 'Tetrahedron'

def parse_input(in_path):
    text = in_path.read_text(errors='replace')
    vals = {}
    vals['calc'] = rgx(text, r"calculation\s*=\s*'([^']+)'")
    vals['ecutwfc'] = rgx(text, r'ecutwfc\s*=\s*([\d.]+)')
    vals['ecutrho'] = rgx(text, r'ecutrho\s*=\s*([\d.]+)')
    vals['smearing'] = rgx(text, r"smearing\s*=\s*'([^']+)'")
    vals['degauss'] = rgx(text, r'degauss\s*=\s*([\d.Ee+-]+)')
    vals['nspin'] = rgx(text, r'nspin\s*=\s*(\d+)')
    vals['start_mag'] = rgx(text, r'starting_magnetization\(\d+\)\s*=\s*([\d.+-]+)')
    vals['vdw'] = rgx(text, r"vdw_corr\s*=\s*'([^']+)'")
    vals['mix_mode'] = rgx(text, r"mixing_mode\s*=\s*'([^']+)'")
    vals['beta'] = rgx(text, r'mixing_beta\s*=\s*([\d.]+)')
    vals['forc_conv'] = rgx(text, r'forc_conv_thr\s*=\s*([\d.De+-]+)')
    vals['nat'] = rgx(text, r'\bnat\s*=\s*(\d+)')
    vals['ntyp'] = rgx(text, r'\bntyp\s*=\s*(\d+)')
    vals['ibrav'] = rgx(text, r'\bibrav\s*=\s*(\d+)')
    return vals

def main():
    rows = []
    for metal in sorted(p.name for p in METALS.iterdir() if p.is_dir()):
        for shape in ['Rhom', 'Tetra']:
            base = METALS / metal / shape
            inp = next(base.glob('INPUT_*_ecut40.4_500.in'), None)
            out = next(base.glob('OUTPUT_*_ecut40.4_500.out'), None)
            if not inp or not out:
                continue
            pos = parse_positions(out.read_text(errors='replace'))
            b = bonds(pos)
            geo = classify(pos)
            spec = parse_out(out)
            params = parse_input(inp)
            pairs = b['pairs']
            rows.append({
                'System': f'{metal}_{shape}', 'Metal': metal, 'Geometry': geo,
                'Status': 'COMPLETE' if spec['final_energy'] else 'PARSE_ERROR',
                'Calculation': params['calc'], 'Ecutwfc_Ry': params['ecutwfc'],
                'Ecutrho_Ry': params['ecutrho'], 'Smearing': params['smearing'],
                'Degauss_Ry': params['degauss'], 'Nspin': params['nspin'],
                'Start_Magnetization': params['start_mag'], 'VDW': params['vdw'],
                'Force_Conv': params['forc_conv'],
                'Max_Force': spec['total_force'], 'SCF_Correction': spec['scf_correction'],
                'SCF_Iterations': spec['scf_iters'], 'SCF_Cycles': spec['scf_cycles'],
                'BFGS_Steps': spec['bfgs_steps'],
                'Final_Energy_Ry': spec['final_energy'], 'SCF_Energy_Ry': spec['scf_energy'],
                'Fermi_eV': spec['fermi'], 'Magnetization': spec['mag'],
                'Abs_Magnetization': spec['abs_mag'],
                'CPU_Time': spec['cpu'], 'Wall_Time': spec['wall'],
                'Num_Atoms': len(pos), 'Cell_A': pos[0]['x']*2 if pos else None,
                'Min_Bond_A': b['min'], 'Max_Bond_A': b['max'], 'Avg_Edge_A': sum(pairs[:3])/3,
                'Overall_Avg_A': b['avg'],
                'Bond_M1M2_A': pairs[0], 'Bond_M1M3_A': pairs[1], 'Bond_M1M4_A': pairs[2],
                'Bond_M2M3_A': pairs[3], 'Bond_M2M4_A': pairs[4], 'Bond_M3M4_A': pairs[5],
                'Short_Diag_A': min(pairs[3:]) if len(pairs) > 3 else None,
                'Long_Diag_A': max(pairs[3:]) if len(pairs) > 3 else None,
                'Bond_Range_A': b['max'] - b['min'],
                'Homogeneity': 'High' if (b['max'] - b['min']) < 0.1 else 'Medium' if (b['max'] - b['min']) < 0.5 else 'Low',
            })
    fieldnames = list(rows[0].keys()) if rows else []
    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f'Wrote {OUT} with {len(rows)} rows')
    for r in rows:
        print(f"{r['System']:<12} {r['Geometry']:<10} E={r['Final_Energy_Ry']:<14} M={r['Magnetization']:<6} Bond={r['Min_Bond_A']}-{r['Max_Bond_A']}")

if __name__ == '__main__':
    main()
