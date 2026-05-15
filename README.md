# Toolkit

Theoretical calculation toolkit for structure manipulation and analysis. Supports VASP POSCAR, XYZ, and CIF structure files.

## Installation

```bash
pip install -e .
```

Requires Python >= 3.9. Dependencies: numpy, pandas, scipy, click, questionary. Optional: matplotlib, seaborn (visualization), ase (adsorption analysis).

## Quick Start

Toolkit provides three ways to use it:

### 1. Interactive Menu

```bash
tk
```

Arrow-key navigation, no arguments needed:

```
  Toolkit
  ========
  > Structure Conversion
    Structure Manipulation
    Analysis
    Quit
```

### 2. Command-Line Arguments

```bash
# Convert structure to POSCAR
tk convert tovasp input.xyz -o output.vasp --direct

# Move atoms
tk manipulate move POSCAR --atoms "1 23-24" --displacement "1 1 1"

# Join two structures
tk manipulate joint POSCAR1 POSCAR2 --axis z --gap 0.5

# Create supercell
tk manipulate extend POSCAR --nx 3 --ny 3 --nz 1

# Density profile from trajectory
tk analyze density pos_*.xyz --axis z --element O --interval 0.2

# Bond angle analysis
tk analyze angle pos_*.xyz --period xy --main O --sub H --cutoff 1.1

# CP2K to DeepMD raw format
tk convert cp2raw -i cp2k.out -p PROJECT_NAME
```

### 3. Python Library

```python
from toolkit import Structure, Writefile, parse_cp2k_inp
from toolkit.analysis import density_profile, bond_angle_analysis

# Read structure
s = Structure("POSCAR")
print(s.lattice())
print(s.coord())

# Write structure
coords = s.coord()
wf = Writefile(coords, s.lattice())
wf.tovasp("output.vasp", direct=True)

# Read lattice from CP2K input
params = parse_cp2k_inp("cp2k.inp")

# Density analysis
counts, density = density_profile(s.coord(), s.lattice(), axis='z', element='O')

# Bond angle analysis
normal, abnormal = bond_angle_analysis(
    s.coord(), s.lattice(), period='xy', main_ele='O', sub_ele='H', cutoff=1.1
)
```

## Commands

| Group | Command | Description |
|-------|---------|-------------|
| `convert` | `tovasp` | Convert structure to VASP POSCAR |
| `convert` | `c2d` | Cartesian to direct coordinates |
| `convert` | `cp2raw` | CP2K output to DeepMD raw files |
| `analyze` | `density` | Density profile along axis |
| `analyze` | `angle` | Bond angle and surface analysis |
| `manipulate` | `move` | Move selected atoms |
| `manipulate` | `joint` | Join two structures |
| `manipulate` | `extend` | Create supercell |

## Supported Formats

- **VASP POSCAR** — `.vasp`, `POSCAR`
- **XYZ** — `.xyz` (requires lattice parameters)
- **CIF** — `.cif`
- **CP2K** — `.inp` (lattice extraction), energy/force/coordinate files

## License

MIT
