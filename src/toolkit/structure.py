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
        if self._filetype in ('vasp', 'cif', 'xyz'):
            return self._filetype
        raise ValueError('Structure input type not recognized (vasp/cif/xyz)')

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
        if direct:
            coord_d = np.dot(coord, np.linalg.inv(self.lattice()))
            return pd.DataFrame(coord_d, columns=['x', 'y', 'z'], index=coord.index)
        return coord

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
