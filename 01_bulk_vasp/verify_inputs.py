#!/usr/bin/env python3
"""Read-only validation for this fixed LiCoO2 conversion; never reads POTCAR."""
from pathlib import Path
from decimal import Decimal
import hashlib
import json
import re

root = Path(__file__).resolve().parent
qe_bytes = (root / 'reference/LCO_bulk.in').read_bytes()
qe = qe_bytes.decode('utf-8-sig').splitlines()

def require(condition, label):
    if not condition:
        raise SystemExit('FAIL: ' + label)
    print('PASS: ' + label)

def block(prefix, count):
    i = next(i for i, line in enumerate(qe) if line.strip().startswith(prefix))
    return [line.split() for line in qe[i+1:i+1+count]]

provenance = json.loads((root / 'provenance.json').read_text())
require(hashlib.sha256(qe_bytes).hexdigest() == provenance['source_sha256'],
        'reference QE bytes match recorded SHA256')
pos = (root / 'POSCAR').read_text().splitlines()
require(len(pos) == 20 and len(pos[0]) <= 40, 'POSCAR has 20 records and short header')
require(pos[1] == '1.0' and pos[5].split() == ['Li', 'Co', 'O']
        and pos[6].split() == ['3', '3', '6'] and pos[7] == 'Direct',
        'POSCAR scale, species, counts and Direct convention')
lattice = block('CELL_PARAMETERS', 3)
atoms = block('ATOMIC_POSITIONS', 12)
require([p.split() for p in pos[2:5]] == lattice,
        'all nine lattice numeric tokens exactly equal QE')
require([p.split() for p in pos[8:]] == [row[1:] for row in atoms],
        'all 36 position numeric tokens exactly equal QE, including signs/boundaries')
require([row[0] for row in atoms] == ['Li']*3 + ['Co']*3 + ['O']*6,
        'original atom sequence matches Li3 Co3 O6')
require(all(len(p.split()) == 3 for p in pos[2:5]+pos[8:]),
        'no extra fields or inline comments in POSCAR numeric records')
require(all('#' not in p and '!' not in p for p in pos[1:]),
        'POSCAR comments confined to header')
matrix = [[Decimal(x) for x in row] for row in lattice]
a,b,c = matrix
det = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
       + a[2]*(b[0]*c[1]-b[1]*c[0]))
require(det > 0, 'cell is nonsingular and right-handed')
print('Cell volume (angstrom^3):', det)

kp = (root / 'KPOINTS').read_text().splitlines()
require(len(kp)==5 and kp[1:] == ['0','Gamma','12 12 4','0 0 0'],
        'KPOINTS is exactly the requested automatic Gamma 12x12x4 zero-shift format')
require(block('K_POINTS',1)[0] == ['12','12','4','0','0','0'],
        'mesh dimensions and zero shift match QE source')

tags = {}
incar_lines = (root / 'INCAR').read_text().splitlines()
for i,line in enumerate(incar_lines):
    data = re.split('[#!]',line,maxsplit=1)[0].strip()
    if not data:
        continue
    require('=' in data and ';' not in data, 'INCAR record syntax line '+str(i+1))
    key,value = [x.strip() for x in data.split('=',1)]
    require(key not in tags, 'INCAR unique tag '+key)
    j = i-1
    while j>=0 and not incar_lines[j].strip():
        j-=1
    require(j>=0 and incar_lines[j].lstrip().startswith('#'),
            'explanatory comment for '+key)
    tags[key] = value

expected = {
 'SYSTEM':'LiCoO2_bulk_PBE_U3.2_D3','GGA':'PE','ENCUT':'550','PREC':'Accurate',
 'EDIFF':'1E-6','NELM':'500','ALGO':'Normal','LREAL':'.FALSE.','LASPH':'.TRUE.',
 'ISMEAR':'0','SIGMA':'0.13606','ISPIN':'2','MAGMOM':'3*0.0 3*1.0 6*0.0',
 'LDAU':'.TRUE.','LDAUTYPE':'2','LDAUL':'-1 2 -1','LDAUU':'0 3.2 0',
 'LDAUJ':'0 0 0','LMAXMIX':'4','IVDW':'11','IBRION':'2','ISIF':'3',
 'PSTRESS':'0','NSW':'100','EDIFFG':'-0.0103','KPAR':'4','NCORE':'1',
 'ISTART':'0','ICHARG':'2','POMASS':'6.941 58.933 15.999',
}
require(tags == expected, 'all 31 explicit INCAR values match documented conversion')
require(len(tags['LDAUL'].split())==3 and len(tags['LDAUU'].split())==3
        and len(tags['LDAUJ'].split())==3, 'Hubbard arrays align with three species')
require(Decimal(tags['LDAUU'].split()[1])-Decimal(tags['LDAUJ'].split()[1])
        == Decimal('3.2'), 'Co Ueff is 3.2 eV')
require(sum(int(x.split('*')[0]) for x in tags['MAGMOM'].split()) == 12,
        'MAGMOM expands to twelve atom moments')

script = (root / 'run_vasp.sh').read_text().splitlines()
directives = [line for line in script if line.startswith('#SBATCH')]
require(all('#' not in line[len('#SBATCH'):] for line in directives),
        'SBATCH directives contain no trailing inline comments')
for directive in ['--ntasks=40','--nodes=1','--cpus-per-task=1','--partition=high']:
    require('#SBATCH '+directive in directives, 'Slurm '+directive)
first_command = next(i for i,line in enumerate(script)
                     if line.strip() and not line.lstrip().startswith('#'))
require(all(i<first_command for i,line in enumerate(script) if line.startswith('#SBATCH')),
        'all scheduler directives precede executable shell statements')
require(40 % int(tags['KPAR']) == 0 and 40//int(tags['KPAR'])==10,
        '40 MPI ranks split into four groups of ten')

for path in root.rglob('*'):
    if not path.is_file():
        continue
    require(path.name != 'POTCAR' and path.suffix.upper() != '.UPF',
            'no proprietary dataset file: '+path.relative_to(root).as_posix())
    if path.relative_to(root).as_posix() == 'reference/LCO_bulk.in':
        continue
    data = path.read_bytes()
    require(not data.startswith(b'\xef\xbb\xbf') and b'\r' not in data,
            'LF/no BOM: '+path.relative_to(root).as_posix())

manifest = root / 'MANIFEST.sha256'
if manifest.is_file():
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ',1)
        require(hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,
                'manifest checksum: '+name)
print('All package checks passed. No POTCAR contents read; no VASP execution.')
