import click
import pandas as pd
import numpy as np
import glob as glob_mod
import os

from toolkit.structure import Structure
from toolkit.io import parse_cp2k_inp
from toolkit.analysis.density import density_profile
from toolkit.analysis.angle import bond_angle_analysis


@click.group()
def analyze():
    """Analyze structures and trajectories."""
    pass


@analyze.command()
@click.argument('files', nargs=-1)
@click.option('--axis', default='z', type=click.Choice(['x', 'y', 'z']))
@click.option('--element', '-e', default='O', help='Element to analyze')
@click.option('--interval', '-i', default=0.2, type=float, help='Bin width')
@click.option('--min', 'vmin', type=float, default=None)
@click.option('--max', 'vmax', type=float, default=None)
@click.option('--lattice', '-l', default=None,
              help='Lattice params for xyz files (comma-separated)')
@click.option('-o', '--output', default=None, help='Output CSV file')
@click.option('--plot/--no-plot', default=True, help='Show density plot')
def density(files, axis, element, interval, vmin, vmax, lattice, output, plot):
    """Calculate density profile from xyz trajectory files."""
    files = _expand_files(files)
    lattice_params = _get_lattice(lattice)

    all_counts, all_density = None, None
    valid = 0
    for fp in files:
        s = Structure(fp, lattice_params)
        counts, dens = density_profile(
            s.coord(), s.lattice(), axis=axis, element=element,
            interval=interval, vmin=vmin, vmax=vmax
        )
        if all_counts is None:
            all_counts = counts
            all_density = dens
        else:
            all_counts = all_counts + counts
            all_density = all_density + dens
        valid += 1
        click.echo(f'  {fp}')

    all_density = all_density / valid * 18 / 6.02 * 10
    click.echo(f'Processed {valid} files')

    if output:
        all_density.to_csv(output, header=False)
        click.echo(f'Saved to {output}')

    if plot:
        import matplotlib.pyplot as plt
        all_density.plot(kind='line', marker='o', color='r')
        out_png = output.replace('.csv', '.png') if output else 'density.png'
        plt.savefig(out_png)
        click.echo(f'Plot saved to {out_png}')


@analyze.command()
@click.argument('files', nargs=-1)
@click.option('--period', '-p', type=click.Choice(['xy', 'yz', 'xz', 'xyz']),
              required=True, help='Periodicity direction')
@click.option('--main', '-m', default='O', help='Central element')
@click.option('--sub', '-s', default='H', help='Bonded element')
@click.option('--cutoff', '-c', default=1.1, type=float, help='Bond cutoff (A)')
@click.option('--num-coord', '-n', default=2, type=int,
              help='Expected coordination number')
@click.option('--lattice', '-l', default=None,
              help='Lattice params (comma-separated)')
@click.option('--output', '-o', default=None, help='Output prefix for CSV files')
def angle(files, period, main, sub, cutoff, num_coord, lattice, output):
    """Analyze bond angles from xyz trajectory files."""
    files = _expand_files(files)
    lattice_params = _get_lattice(lattice)

    all_normal = pd.DataFrame()
    all_abnormal = pd.DataFrame()

    for fp in files:
        s = Structure(fp, lattice_params)
        norm, abnorm = bond_angle_analysis(
            s.coord(), s.lattice(), period=period,
            main_ele=main, sub_ele=sub, cutoff=cutoff, num_coord=num_coord
        )
        all_normal = pd.concat([all_normal, norm], axis=0)
        all_abnormal = pd.concat([all_abnormal, abnorm], axis=0)
        click.echo(f'  {fp}')

    prefix = output or 'angle'
    all_normal.index.name = 'atom_index'
    all_abnormal.index.name = 'atom_index'
    n_file = f'{prefix}_normal.csv'
    a_file = f'{prefix}_abnormal.csv'
    all_normal.to_csv(n_file)
    all_abnormal.to_csv(a_file)
    click.echo(f'Saved to {n_file}, {a_file}')


def _expand_files(patterns):
    if not patterns:
        raise click.UsageError('At least one file or pattern required')
    result = []
    for pat in patterns:
        if '*' in pat or '?' in pat:
            result.extend(sorted(glob_mod.glob(pat)))
        else:
            result.append(pat)
    if not result:
        raise click.UsageError(f'No files matched: {patterns}')
    return result


def _get_lattice(lattice_arg):
    if lattice_arg:
        return [float(x) for x in lattice_arg.split(',')]
    return parse_cp2k_inp()
