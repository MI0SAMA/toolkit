import click
import numpy as np
import pandas as pd

from toolkit.structure import Structure
from toolkit.io import Writefile


@click.group()
def manipulate():
    """Manipulate structures (move atoms, join, extend)."""
    pass


def _parse_atom_range(spec):
    nums = []
    for part in spec.split():
        if '-' in part:
            lo, hi = part.split('-')
            nums.extend(range(int(lo), int(hi) + 1))
        else:
            nums.append(int(part))
    return nums


@manipulate.command()
@click.argument('infile')
@click.option('--atoms', '-a', default=None,
              help='Atom numbers to move, e.g. "1 23-24"')
@click.option('--displacement', '-d', default=None,
              help='Displacement vector, e.g. "1 1 1"')
@click.option('--direct/--cartesian', default=False,
              help='Input displacement in direct (fractional) coordinates')
@click.option('-o', '--output', default='POSCAR_moved', help='Output file')
def move(infile, atoms, displacement, direct, output):
    """Move selected atoms in a structure."""
    s = Structure(infile)

    if not atoms:
        atoms = click.prompt('Atom numbers to move (e.g. 1 23-24)')
    if not displacement:
        displacement = click.prompt('Displacement vector (e.g. 1 1 1)')

    atom_nums = _parse_atom_range(atoms)
    disp_vec = [float(x) for x in displacement.split()]

    moved_coord = s.move(atom_nums, disp_vec, direct=direct)
    wf = Writefile(moved_coord, s.lattice())
    wf.tovasp(output, direct=False)
    click.echo(f'Written to {output}')


@manipulate.command()
@click.argument('infile1')
@click.argument('infile2')
@click.option('--axis', '-a', default='z',
              type=click.Choice(['x', 'y', 'z']),
              help='Axis along which to join structures')
@click.option('--gap', '-g', default=0.5, type=float,
              help='Gap between structures (Angstrom)')
@click.option('-o', '--output', default='POSCAR_joint', help='Output file')
def joint(infile1, infile2, axis, gap, output):
    """Join two structures along a given axis."""
    s1 = Structure(infile1)
    s2 = Structure(infile2)

    lattice1 = s1.lattice()
    lattice2 = s2.lattice()
    coord1 = s1.coord()
    coord2 = s2.coord()

    axis_idx = {'x': 0, 'y': 1, 'z': 2}
    idx = axis_idx[axis]
    other_axes = [0, 1, 2]
    other_axes.pop(idx)

    if not (np.allclose(lattice1[other_axes[0]], lattice2[other_axes[0]]) and
            np.allclose(lattice1[other_axes[1]], lattice2[other_axes[1]])):
        click.echo('Warning: lattices do not match in the joining plane')

    shift = lattice1[idx][idx] + gap
    coord2[axis] = coord2[axis] + shift

    new_coord = pd.concat([coord1, coord2])
    new_lattice = lattice1.copy()
    new_lattice[idx] = lattice1[idx] + lattice2[idx]
    new_lattice[idx][idx] += gap

    wf = Writefile(new_coord, new_lattice)
    wf.tovasp(output, direct=True)
    click.echo(f'Written to {output}')


@manipulate.command()
@click.argument('infile')
@click.option('--nx', type=int, default=1, help='Repeat in x')
@click.option('--ny', type=int, default=1, help='Repeat in y')
@click.option('--nz', type=int, default=1, help='Repeat in z')
@click.option('-o', '--output', default='POSCAR_ext', help='Output file')
def extend(infile, nx, ny, nz, output):
    """Create a supercell by repeating the structure."""
    s = Structure(infile)
    ext_lattice, ext_coord = s.extend(nx, ny, nz)
    wf = Writefile(ext_coord, ext_lattice)
    wf.tovasp(output, direct=True)
    click.echo(f'Written to {output}')
