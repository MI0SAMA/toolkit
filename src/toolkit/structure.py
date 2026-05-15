import numpy as np
import pandas as pd
import itertools


def _is_num(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False


_MASS_TO_ELEM = {
    1.008: 'H', 4.0026: 'He', 6.94: 'Li', 9.012: 'Be', 10.81: 'B',
    12.011: 'C', 14.007: 'N', 15.999: 'O', 18.998: 'F', 20.180: 'Ne',
    22.990: 'Na', 24.305: 'Mg', 26.982: 'Al', 28.085: 'Si', 30.974: 'P',
    32.06: 'S', 35.45: 'Cl', 39.098: 'K', 40.078: 'Ca', 44.956: 'Sc',
    47.867: 'Ti', 50.942: 'V', 51.996: 'Cr', 54.938: 'Mn', 55.845: 'Fe',
    58.933: 'Co', 58.693: 'Ni', 63.546: 'Cu', 65.38: 'Zn', 69.723: 'Ga',
    72.630: 'Ge', 74.922: 'As', 78.971: 'Se', 79.904: 'Br', 83.798: 'Kr',
    85.468: 'Rb', 87.62: 'Sr', 88.906: 'Y', 91.224: 'Zr', 92.906: 'Nb',
    95.95: 'Mo', 97.91: 'Tc', 101.07: 'Ru', 102.91: 'Rh', 106.42: 'Pd',
    107.87: 'Ag', 112.41: 'Cd', 114.82: 'In', 118.71: 'Sn', 121.76: 'Sb',
    127.60: 'Te', 126.90: 'I', 131.29: 'Xe', 132.91: 'Cs', 137.33: 'Ba',
    138.91: 'La', 140.12: 'Ce', 140.91: 'Pr', 144.24: 'Nd', 144.91: 'Pm',
    150.36: 'Sm', 151.96: 'Eu', 157.25: 'Gd', 158.93: 'Tb', 162.50: 'Dy',
    164.93: 'Ho', 167.26: 'Er', 168.93: 'Tm', 173.05: 'Yb', 174.97: 'Lu',
    178.49: 'Hf', 180.95: 'Ta', 183.84: 'W', 186.21: 'Re', 190.23: 'Os',
    192.22: 'Ir', 195.08: 'Pt', 196.97: 'Au', 200.59: 'Hg', 204.38: 'Tl',
    207.2: 'Pb', 208.98: 'Bi', 208.98: 'Po', 210.0: 'At', 222.0: 'Rn',
}


def _guess_element(mass):
    """Guess element symbol from atomic mass."""
    best, best_elem = float('inf'), 'X'
    for m, elem in _MASS_TO_ELEM.items():
        diff = abs(m - mass)
        if diff < best:
            best, best_elem = diff, elem
    return best_elem


def _cell_p2v(A, B, C, alpha, beta, gamma):
    """Convert cell parameters (lengths + angles) to 3x3 lattice vectors."""
    alpha, beta, gamma = np.deg2rad([alpha, beta, gamma])
    lx = A
    xy = B * np.cos(gamma)
    xz = C * np.cos(beta)
    ly = B * np.sin(gamma)
    if not ly > 0.1:
        raise RuntimeError(f"ly:=B* sin(gamma)=={ly}, must be greater than 0.1")
    yz = (B * C * np.cos(alpha) - xy * xz) / ly
    lz2 = C ** 2 - xz ** 2 - yz ** 2
    if not lz2 > 0.01:
        raise RuntimeError(f"lz^2:=C**2-xz**2-yz**2=={lz2}, must be greater than 0.01")
    lz = np.sqrt(lz2)
    return np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]], dtype=float)


class Structure:
    """Read a structure file (VASP POSCAR, CIF, XYZ) and provide lattice + coordinates."""

    def __init__(self, filename, lattice_params=None):
        self._filename = str(filename)
        if '.' in self._filename:
            self._filetype = self._filename.split('.')[-1].lower()
        elif 'poscar' in self._filename.lower():
            self._filetype = 'vasp'
        with open(self._filename, 'r') as f:
            self._filecontent = [x.rstrip().split() for x in f]
        self._lattice_params = lattice_params

    @property
    def filecontent(self):
        return self._filecontent

    @property
    def filename(self):
        return self._filename

    @property
    def filetype(self):
        if self._filetype in ('vasp', 'cif', 'xyz', 'data'):
            return self._filetype
        raise ValueError('Structure input type not recognized (vasp/cif/xyz/data)')

    def lattice(self):
        """Return 3x3 lattice matrix as np.ndarray."""
        if self.filetype == 'vasp':
            lattice = np.array(self.filecontent[2:5], dtype=float)
        elif self.filetype == 'cif':
            dic = {}
            for i in self.filecontent:
                if len(i) == 2:
                    dic[i[0]] = i[1]
            lg = np.array([dic['_cell_length_a'], dic['_cell_length_b'],
                           dic['_cell_length_c']], dtype=float)
            ag = np.array([dic['_cell_angle_alpha'], dic['_cell_angle_beta'],
                           dic['_cell_angle_gamma']], dtype=float)
            lattice = _cell_p2v(lg[0], lg[1], lg[2], ag[0], ag[1], ag[2])
        elif self.filetype == 'xyz':
            if self._lattice_params is None:
                raise ValueError('Lattice parameters required for .xyz files')
            if len(self._lattice_params) == 6:
                A, B, C, alpha, beta, gamma = self._lattice_params
                lattice = _cell_p2v(A, B, C, alpha, beta, gamma)
            elif len(self._lattice_params) == 9:
                lattice = np.array(self._lattice_params).reshape(3, 3)
            else:
                raise ValueError(
                    'Lattice params must be 6 values (A,B,C,alpha,beta,gamma) '
                    'or 9 values (xx,xy,xz,yx,yy,yz,zx,zy,zz)'
                )
        elif self.filetype == 'data':
            box = self._parse_lammps_sections()['box']
            xlo, xhi = box['xlo'], box['xhi']
            ylo, yhi = box['ylo'], box['yhi']
            zlo, zhi = box['zlo'], box['zhi']
            xy, xz, yz = box.get('xy', 0.0), box.get('xz', 0.0), box.get('yz', 0.0)
            lattice = np.array([
                [xhi - xlo, 0, 0],
                [xy, yhi - ylo, 0],
                [xz, yz, zhi - zlo],
            ], dtype=float)
        else:
            raise ValueError('Unknown file type')
        lattice[np.logical_and(lattice < 0.0001, lattice > 0)] = 0
        return lattice

    def coord(self, direct=False):
        """Return atomic coordinates as pd.DataFrame with element-indexed rows."""
        if self.filetype == 'vasp':
            at_num = list(map(int, self.filecontent[6]))
            at_type = dict(zip(self.filecontent[5], at_num))
            for i in (8, 9):
                if _is_num(self.filecontent[i][0]):
                    coord_arr = np.array(
                        self.filecontent[i:i + sum(at_num)]
                    )[:, :3].astype(float)
                    if self.filecontent[i - 1][0][0] == 'D':
                        coord_arr = np.dot(coord_arr, self.lattice())
                    break
            coord = pd.DataFrame(
                coord_arr, columns=['x', 'y', 'z'],
                index=[val for val in at_type for _ in range(at_type[val])]
            )
        elif self.filetype == 'xyz':
            at_num = int(self.filecontent[0][0])
            coord_arr = np.array(self.filecontent[2:2 + at_num])[:, 0:4]
            coord = pd.DataFrame(
                coord_arr[:, 1:4].astype(float), columns=['x', 'y', 'z'],
                index=coord_arr[:, 0]
            )
        elif self.filetype == 'cif':
            count = 0
            for i in self.filecontent:
                if i != [] and i[0] == '_atom_site_type_symbol':
                    break
                count += 1
            coord_arr = np.array(self.filecontent[count + 1:])
            coord = pd.DataFrame(
                np.dot(coord_arr[:, 2:5].astype(float), self.lattice()),
                columns=['x', 'y', 'z'],
                index=coord_arr[:, -1]
            )
        elif self.filetype == 'data':
            sections = self._parse_lammps_sections()
            atoms_data = sections['atoms']
            mass_map = sections['mass_map']
            ncols = atoms_data.shape[1]
            if ncols == 5:
                cols = ['id', 'type', 'x', 'y', 'z']
            elif ncols == 6:
                cols = ['id', 'type', 'q', 'x', 'y', 'z']
            elif ncols == 7:
                cols = ['id', 'mol', 'type', 'q', 'x', 'y', 'z']
            elif ncols >= 8:
                cols = ['id', 'mol', 'type', 'q', 'x', 'y', 'z'] + [f'extra_{k}' for k in range(ncols - 8)]
            else:
                cols = ['id', 'type', 'x', 'y', 'z']
            df = pd.DataFrame(atoms_data, columns=cols)
            type_col = 'type' if 'type' in cols else cols[1]
            elem_labels = [mass_map.get(int(t), 'X') for t in df[type_col]]
            coord = pd.DataFrame(
                df[['x', 'y', 'z']].astype(float).to_numpy(),
                columns=['x', 'y', 'z'],
                index=elem_labels
            )
        if direct:
            coord_d = np.dot(coord, np.linalg.inv(self.lattice()))
            return pd.DataFrame(coord_d, columns=['x', 'y', 'z'], index=coord.index)
        return coord

    def _parse_lammps_sections(self):
        """Parse LAMMPS data file sections: box, masses, atoms."""
        content = self.filecontent
        lines = content

        # Find section boundaries
        masses_start = None
        masses_end = None
        atoms_start = None
        n_atoms = None
        box_vals = {}

        for idx, tokens in enumerate(lines):
            line_str = ' '.join(tokens).lower()
            # Detect counts line
            if 'atoms' in line_str and len(tokens) > 0 and tokens[0].isdigit():
                n_atoms = int(tokens[0])
            # Box bounds
            if len(tokens) >= 2 and _is_num(tokens[0]) and _is_num(tokens[1]):
                if 'xlo' in line_str and 'xhi' in line_str:
                    box_vals['xlo'] = float(tokens[0])
                    box_vals['xhi'] = float(tokens[1])
                    if len(tokens) >= 3 and _is_num(tokens[2]):
                        box_vals['xy'] = float(tokens[2])
                elif 'ylo' in line_str and 'yhi' in line_str:
                    box_vals['ylo'] = float(tokens[0])
                    box_vals['yhi'] = float(tokens[1])
                    if len(tokens) >= 3 and _is_num(tokens[2]):
                        box_vals['xz'] = float(tokens[2])
                elif 'zlo' in line_str and 'zhi' in line_str:
                    box_vals['zlo'] = float(tokens[0])
                    box_vals['zhi'] = float(tokens[1])
                    if len(tokens) >= 3 and _is_num(tokens[2]):
                        box_vals['yz'] = float(tokens[2])
                # Standalone tilt line: "xy xz yz" at end of box section
                elif 'xy' in line_str and 'xz' in line_str and 'yz' in line_str:
                    box_vals['xy'] = float(tokens[0])
                    box_vals['xz'] = float(tokens[1])
                    box_vals['yz'] = float(tokens[2])
            # Masses section start
            if not tokens:
                continue
            if tokens[0].lower() == 'masses':
                masses_start = idx + 1
                continue
            # Atoms section start
            if tokens[0].lower() == 'atoms':
                if masses_start is not None:
                    masses_end = idx
                atoms_start = idx + 1
                continue
            # Detect end of masses section (first line after masses that isn't a "type mass" pair)
            if masses_start is not None and masses_end is None and atoms_start is None:
                if idx > masses_start:
                    if not _is_num(tokens[0]):
                        masses_end = idx

        # Auto-detect masses_end if it was never set
        if masses_start is not None and masses_end is None:
            for idx in range(masses_start, len(lines)):
                tokens = lines[idx]
                if not tokens or not _is_num(tokens[0]):
                    masses_end = idx
                    break
            if masses_end is None:
                masses_end = atoms_start if atoms_start else len(lines)

        # Parse mass_map: type_num -> element symbol
        mass_map = {}
        if masses_start is not None:
            for idx in range(masses_start, masses_end or len(lines)):
                tokens = lines[idx]
                if not tokens or len(tokens) < 2 or not _is_num(tokens[0]):
                    if len(tokens) == 0 or not _is_num(tokens[0]):
                        continue
                    break
                type_num = int(tokens[0])
                mass = float(tokens[1])
                mass_map[type_num] = _guess_element(mass)

        # Parse atoms
        atoms_rows = []
        if atoms_start is not None:
            for idx in range(atoms_start, len(lines)):
                tokens = lines[idx]
                if not tokens or len(tokens) < 4:
                    continue
                if not _is_num(tokens[0]):
                    break
                atoms_rows.append([float(t) for t in tokens])

        return {
            'atoms': np.array(atoms_rows) if atoms_rows else np.empty((0, 5)),
            'mass_map': mass_map,
            'box': box_vals,
        }

    def extend(self, x, y, z):
        """Create supercell. Returns (ext_lattice, ext_coord)."""
        org_coord = self.coord()
        org_lattice = self.lattice()
        ext_lattice = org_lattice * [x, y, z]

        ranges = [list(range(1, n + 1)) for n in (x, y, z)]
        extend_list = np.array(list(itertools.product(*ranges)))

        ext_coord = pd.DataFrame()
        for idx in extend_list:
            shifted = org_coord + np.dot(idx - [1, 1, 1], org_lattice)
            ext_coord = pd.concat([ext_coord, shifted])
        return ext_lattice, ext_coord

    def move(self, atom_nums, displacement, direct=False):
        """Move specified atoms by displacement vector.

        Args:
            atom_nums: list of 1-indexed atom numbers to move
            displacement: [dx, dy, dz]
            direct: whether input coordinates are in direct (fractional) space
        """
        coord = self.coord(direct=direct)
        for idx in atom_nums:
            coord.iloc[idx - 1] = coord.iloc[idx - 1] + displacement
        if direct:
            cart = np.dot(coord, self.lattice())
            return pd.DataFrame(cart, columns=['x', 'y', 'z'], index=coord.index)
        return coord
