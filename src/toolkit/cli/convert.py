import click
import numpy as np
import os
import glob as glob_mod

from toolkit.structure import Structure
from toolkit.io import Writefile, parse_cp2k_inp


@click.group()
def convert():
    """Convert between structure file formats."""
    pass


@convert.command()
@click.argument('infile')
@click.option('-o', '--outfile', default=None,
              help='Output file (default: infile base + .vasp)')
@click.option('--lattice', '-l', default=None,
              help='Lattice params: comma-separated 6 or 9 values for xyz files')
@click.option('--direct/--cartesian', default=True,
              help='Write direct (fractional) or Cartesian coordinates')
def tovasp(infile, outfile, lattice, direct):
    """Convert a structure file to VASP POSCAR format."""
    lattice_params = None
    if lattice:
        lattice_params = [float(x) for x in lattice.split(',')]
    else:
        lattice_params = parse_cp2k_inp()
        if lattice_params:
            click.echo(f'Read lattice from cp2k.inp: {lattice_params}')

    s = Structure(infile, lattice_params)
    wf = Writefile(s.coord(), s.lattice(), infile)

    if outfile is None:
        outfile = os.path.splitext(infile)[0] + '.vasp'
    wf.tovasp(outfile, direct=direct)
    click.echo(f'Written to {outfile}')


@convert.command()
@click.argument('infile')
@click.option('-o', '--outfile', default=None, help='Output file')
@click.option('--lattice', '-l', default=None,
              help='Lattice params: comma-separated 6 or 9 values for xyz files')
def c2d(infile, outfile, lattice):
    """Convert Cartesian to Direct coordinates (writes POSCAR)."""
    lattice_params = None
    if lattice:
        lattice_params = [float(x) for x in lattice.split(',')]
    else:
        lattice_params = parse_cp2k_inp()
    if lattice_params is None:
        raise click.UsageError(
            'Lattice required. Use --lattice or place cp2k.inp in current directory.'
        )
    s = Structure(infile, lattice_params)
    wf = Writefile(s.coord(), s.lattice(), infile)
    if outfile is None:
        outfile = infile + '_d'
    wf.tovasp(outfile, direct=True)
    click.echo(f'Written to {outfile}')


@convert.command()
@click.option('--infile', '-i', default='cp2k.out',
              help='CP2K output file (default: cp2k.out)')
@click.option('--project', '-p', default=None,
              help='Project name (auto-detected if omitted)')
@click.option('--step-start', 'step_start', type=int, default=0,
              help='Start step index (default: 0)')
@click.option('--step-end', 'step_end', type=int, default=None,
              help='End step index (default: last)')
@click.option('--output-dir', '-o', default='raw',
              help='Output directory (default: raw/)')
def cp2raw(infile, project, step_start, step_end, output_dir):
    """Convert CP2K output to DeepMD raw files."""
    if project is None:
        with open(infile, 'r') as f:
            for line in f:
                parts = line.rstrip().split()
                if 'Project' in parts and 'name' in parts:
                    project = parts[-1]
                    break
        if project is None:
            raise click.UsageError('Cannot detect project name. Use --project.')

    element_type, element_num, elements = _read_cp2k_types(project)

    ener_file = f'{project}-1.ener'
    with open(ener_file) as f:
        total_steps = sum(1 for _ in f)
    if step_end is None:
        step_end = total_steps - 1

    energy = _read_cp2k_energy(project)
    force = _read_cp2k_force(project, element_num)
    coord = _read_cp2k_coord(project, element_num)
    box = _read_cp2k_box(infile)

    os.makedirs(output_dir, exist_ok=True)
    n_steps = step_end - step_start

    with open(os.path.join(output_dir, 'box.raw'), 'w') as fb, \
         open(os.path.join(output_dir, 'energy.raw'), 'w') as fe, \
         open(os.path.join(output_dir, 'force.raw'), 'w') as ff, \
         open(os.path.join(output_dir, 'coord.raw'), 'w') as fc:
        for i in range(n_steps):
            idx = step_start + i
            fb.write(' '.join(str(v) for v in box) + '\n')
            fe.write(f'{energy[idx]}\n')
            ff.write(' '.join(str(v) for v in force[idx]) + '\n')
            fc.write(' '.join(str(v) for v in coord[idx]) + '\n')

    with open(os.path.join(output_dir, 'type.raw'), 'w') as ft:
        for elem in elements:
            ft.write(f'{element_type[elem]}\n')
    with open(os.path.join(output_dir, 'type_map.raw'), 'w') as fm:
        for etype in element_type:
            fm.write(f'{etype}\n')

    click.echo(f'Written {n_steps} frames to {output_dir}/')


# CP2K reading helpers

def _read_cp2k_types(project):
    elements = []
    with open(f'{project}-pos-1.xyz') as f:
        lines = f.readlines()
        element_num = int(lines[0].split()[0])
        for line in lines[2:element_num + 2]:
            parts = line.rstrip().split()
            if len(parts) == 4:
                elements.append(parts[0])
    unique = list(dict.fromkeys(elements))
    element_type = dict(zip(unique, range(len(unique))))
    return element_type, element_num, elements


def _read_cp2k_energy(project):
    energy = []
    with open(f'{project}-1.ener') as f:
        for line in f:
            parts = line.rstrip().split()
            try:
                energy.append(float(parts[-2]) * 27.211386)
            except (ValueError, IndexError):
                pass
    return energy


def _read_cp2k_force(project, element_num):
    force_file = f'{project}-force-1.xyz'
    force = []
    factor = 27.211386 / 0.529177
    try:
        with open(force_file) as f:
            frame = []
            count = 0
            for line in f:
                parts = line.rstrip().split()
                if len(parts) == 6:
                    frame.extend(float(parts[i]) * factor for i in (3, 4, 5))
                    count += 1
                    if count == element_num:
                        force.append(frame)
                        frame = []
                        count = 0
    except FileNotFoundError:
        click.echo(f'Warning: {force_file} not found, forces will be empty', err=True)
    return force


def _read_cp2k_coord(project, element_num):
    coord = []
    with open(f'{project}-pos-1.xyz') as f:
        frame = []
        count = 0
        for line in f:
            parts = line.rstrip().split()
            if len(parts) == 4:
                frame.extend(float(parts[i]) for i in (1, 2, 3))
                count += 1
                if count == element_num:
                    coord.append(frame)
                    frame = []
                    count = 0
    return coord


def _read_cp2k_box(infile):
    box = []
    with open(infile, 'r') as f:
        for line in f:
            parts = line.rstrip().split()
            if 'CELL_TOP|' in line and 'Vector' in line:
                box.extend(float(parts[i]) for i in (4, 5, 6))
            if 'CELL_TOP|' in line and 'Angle' in line:
                break
    return box
