"""Interactive menu shell using questionary."""

import questionary
import click
import pandas as pd
import os
import glob as glob_mod

from toolkit.structure import Structure
from toolkit.io import Writefile, parse_cp2k_inp
from toolkit.analysis.density import density_profile
from toolkit.analysis.angle import bond_angle_analysis


def run_menu():
    """Launch the main interactive menu loop."""
    while True:
        choice = questionary.select(
            'Toolkit',
            choices=[
                'Structure Conversion',
                'Structure Manipulation',
                'Analysis',
                'Quit',
            ]
        ).ask()

        if choice is None or choice == 'Quit':
            break
        elif choice == 'Structure Conversion':
            _convert_menu()
        elif choice == 'Structure Manipulation':
            _manipulate_menu()
        elif choice == 'Analysis':
            _analysis_menu()


def _convert_menu():
    choice = questionary.select(
        'Structure Conversion',
        choices=[
            'Convert to VASP (POSCAR)',
            'Cartesian to Direct coordinates',
            'CP2K to DeepMD raw',
            'Back',
        ]
    ).ask()

    if choice is None or choice == 'Back':
        return
    elif 'VASP' in choice:
        infile = questionary.text('Input file:').ask()
        if not infile:
            return
        fmt = questionary.select(
            'Coordinate format:',
            choices=['Direct (fractional)', 'Cartesian'],
            default='Direct (fractional)',
        ).ask()
        direct = 'Direct' in fmt
        outfile = questionary.text('Output file:', default=os.path.splitext(infile)[0] + '.vasp').ask()
        lattice_params = parse_cp2k_inp()
        s = Structure(infile, lattice_params)
        wf = Writefile(s.coord(), s.lattice(), infile)
        wf.tovasp(outfile, direct=direct)
        click.echo(f'Written to {outfile}')
    elif 'Direct' in choice:
        infile = questionary.text('Input file:').ask()
        lattice_params = _prompt_lattice()
        if not lattice_params:
            return
        s = Structure(infile, lattice_params)
        wf = Writefile(s.coord(), s.lattice(), infile)
        outfile = questionary.text('Output file:', default=infile + '_d').ask()
        wf.tovasp(outfile, direct=True)
        click.echo(f'Written to {outfile}')
    elif 'CP2K' in choice:
        infile = questionary.text('CP2K output file:', default='cp2k.out').ask()
        project = questionary.text('Project name:').ask()
        click.echo(f'Run: tk convert cp2raw -i {infile} -p {project}')


def _manipulate_menu():
    choice = questionary.select(
        'Structure Manipulation',
        choices=[
            'Move atoms',
            'Join two structures',
            'Extend supercell',
            'Back',
        ]
    ).ask()

    if choice is None or choice == 'Back':
        return

    if 'Move' in choice:
        infile = questionary.text('Input file:').ask()
        atoms_spec = questionary.text('Atom numbers (e.g. 1 23-24):').ask()
        disp = questionary.text('Displacement vector (e.g. 1 1 1):').ask()
        outfile = questionary.text('Output file:', default='POSCAR_moved').ask()

        s = Structure(infile)
        atom_nums = []
        for part in atoms_spec.split():
            if '-' in part:
                lo, hi = part.split('-')
                atom_nums.extend(range(int(lo), int(hi) + 1))
            else:
                atom_nums.append(int(part))
        disp_vec = [float(x) for x in disp.split()]
        moved = s.move(atom_nums, disp_vec)
        wf = Writefile(moved, s.lattice())
        wf.tovasp(outfile, direct=False)
        click.echo(f'Written to {outfile}')

    elif 'Join' in choice:
        f1 = questionary.text('First structure file:').ask()
        f2 = questionary.text('Second structure file:').ask()
        axis = questionary.select('Join axis:', choices=['z', 'x', 'y'], default='z').ask()
        gap = questionary.text('Gap (Angstrom):', default='0.5').ask()
        outfile = questionary.text('Output file:', default='POSCAR_joint').ask()

        s1 = Structure(f1)
        s2 = Structure(f2)
        l1, l2 = s1.lattice(), s2.lattice()
        c1, c2 = s1.coord(), s2.coord()

        aidx = {'x': 0, 'y': 1, 'z': 2}
        idx = aidx[axis]
        shift = l1[idx][idx] + float(gap)
        c2[axis] = c2[axis] + shift

        new_lat = l1.copy()
        new_lat[idx] = l1[idx] + l2[idx]
        new_lat[idx][idx] += float(gap)

        wf = Writefile(pd.concat([c1, c2]), new_lat)
        wf.tovasp(outfile, direct=True)
        click.echo(f'Written to {outfile}')

    elif 'Extend' in choice:
        infile = questionary.text('Input file:').ask()
        nx = questionary.text('Repeat X:', default='1').ask()
        ny = questionary.text('Repeat Y:', default='1').ask()
        nz = questionary.text('Repeat Z:', default='1').ask()
        outfile = questionary.text('Output file:', default='POSCAR_ext').ask()

        s = Structure(infile)
        ext_lat, ext_coord = s.extend(int(nx), int(ny), int(nz))
        wf = Writefile(ext_coord, ext_lat)
        wf.tovasp(outfile, direct=True)
        click.echo(f'Written to {outfile}')


def _analysis_menu():
    choice = questionary.select(
        'Analysis',
        choices=[
            'Density profile',
            'Bond angle analysis',
            'Back',
        ]
    ).ask()

    if choice is None or choice == 'Back':
        return

    if 'Density' in choice:
        pattern = questionary.text('File pattern (e.g. pos_*.xyz):').ask()
        element = questionary.text('Element:', default='O').ask()
        axis = questionary.select('Axis:', choices=['z', 'x', 'y'], default='z').ask()
        interval = questionary.text('Interval:', default='0.2').ask()
        lattice_params = _prompt_lattice()

        files = sorted(glob_mod.glob(pattern))
        all_counts, all_dens = None, None
        valid = 0
        for fp in files:
            s = Structure(fp, lattice_params)
            counts, dens = density_profile(
                s.coord(), s.lattice(), axis=axis, element=element,
                interval=float(interval)
            )
            if all_counts is None:
                all_counts, all_dens = counts, dens
            else:
                all_counts += counts
                all_dens += dens
            valid += 1

        all_dens = all_dens / valid * 18 / 6.02 * 10
        outfile = questionary.text('Output CSV:', default='density.csv').ask()
        all_dens.to_csv(outfile, header=False)
        click.echo(f'Processed {valid} files -> {outfile}')

    elif 'angle' in choice:
        pattern = questionary.text('File pattern (e.g. pos_*.xyz):').ask()
        period = questionary.select('Periodicity:', choices=['xy', 'yz', 'xz', 'xyz'], default='xy').ask()
        main_ele = questionary.text('Main element:', default='O').ask()
        sub_ele = questionary.text('Sub element:', default='H').ask()
        cutoff = questionary.text('Cutoff (A):', default='1.1').ask()
        lattice_params = _prompt_lattice()

        files = sorted(glob_mod.glob(pattern))
        all_norm, all_abnorm = pd.DataFrame(), pd.DataFrame()
        for fp in files:
            s = Structure(fp, lattice_params)
            n, a = bond_angle_analysis(
                s.coord(), s.lattice(), period=period,
                main_ele=main_ele, sub_ele=sub_ele, cutoff=float(cutoff)
            )
            all_norm = pd.concat([all_norm, n], axis=0)
            all_abnorm = pd.concat([all_abnorm, a], axis=0)
            click.echo(f'  {fp}')

        prefix = questionary.text('Output prefix:', default='angle').ask()
        all_norm.to_csv(f'{prefix}_normal.csv')
        all_abnorm.to_csv(f'{prefix}_abnormal.csv')
        click.echo(f'Saved to {prefix}_normal.csv, {prefix}_abnormal.csv')


def _prompt_lattice():
    """Prompt for lattice parameters, trying cp2k.inp first."""
    from_cp2k = parse_cp2k_inp()
    if from_cp2k:
        click.echo(f'Found cp2k.inp: {from_cp2k}')
        use = questionary.confirm('Use these?', default=True).ask()
        if use:
            return from_cp2k
    inp = questionary.text(
        'Lattice params (6 or 9 comma-separated values):'
    ).ask()
    if inp:
        return [float(x) for x in inp.split(',')]
    return None
