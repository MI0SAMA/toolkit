import pytest
import numpy as np
from toolkit.structure import Structure, _guess_element
from conftest import data_path

LATTICE_CUBIC_10 = np.array([
    [10.0, 0.0, 0.0],
    [0.0, 10.0, 0.0],
    [0.0, 0.0, 10.0],
])


class TestGuessElement:
    def test_hydrogen(self):
        assert _guess_element(1.008) == 'H'

    def test_oxygen(self):
        assert _guess_element(15.999) == 'O'

    def test_carbon(self):
        assert _guess_element(12.011) == 'C'

    def test_iron(self):
        assert _guess_element(55.845) == 'Fe'

    def test_unknown_returns_closest(self):
        result = _guess_element(14.0)
        assert result == 'N'  # nitrogen, mass 14.007


class TestStructureVASP:
    def test_read_lattice_vasp(self):
        s = Structure(data_path('water.vasp'))
        np.testing.assert_array_almost_equal(s.lattice(), LATTICE_CUBIC_10)

    def test_read_coord_vasp(self):
        s = Structure(data_path('water.vasp'))
        coord = s.coord()
        assert list(coord.columns) == ['x', 'y', 'z']
        assert len(coord.loc[['O']]) == 1
        assert len(coord.loc[['H']]) == 2
        # O at (5,5,5)
        np.testing.assert_array_almost_equal(coord.loc[['O']].iloc[0], [5.0, 5.0, 5.0])

    def test_read_direct_coord_vasp(self):
        s = Structure(data_path('water_direct.vasp'))
        coord = s.coord(direct=False)
        # Direct coords * lattice = Cartesian
        np.testing.assert_array_almost_equal(coord.loc[['O']].iloc[0], [5.0, 5.0, 5.0])

    def test_coord_to_direct(self):
        s = Structure(data_path('water.vasp'))
        d = s.coord(direct=True)
        np.testing.assert_array_almost_equal(d.loc['O'].iloc[0], [0.5, 0.5, 0.5])

    def test_filetype_detection(self):
        s = Structure(data_path('water.vasp'))
        assert s.filetype == 'vasp'

    def test_invalid_filetype_raises(self):
        s = Structure(data_path('water.vasp'))
        s._filetype = 'unknown'
        with pytest.raises(ValueError, match='not recognized'):
            _ = s.filetype


class TestStructureXYZ:
    LATTICE_9 = [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0]

    def test_read_with_9_lattice_params(self):
        s = Structure(data_path('water.xyz'), lattice_params=self.LATTICE_9)
        np.testing.assert_array_almost_equal(s.lattice(), LATTICE_CUBIC_10)

    def test_read_with_6_lattice_params(self):
        s = Structure(data_path('water.xyz'), lattice_params=[10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
        np.testing.assert_array_almost_equal(s.lattice(), LATTICE_CUBIC_10)

    def test_read_xyz_without_lattice_raises(self):
        with pytest.raises(ValueError, match='Lattice parameters required'):
            Structure(data_path('water.xyz')).lattice()

    def test_read_coord_xyz(self):
        s = Structure(data_path('water.xyz'), lattice_params=self.LATTICE_9)
        coord = s.coord()
        assert coord.loc[['O']].shape[0] == 1
        assert coord.loc[['H']].shape[0] == 2


class TestStructureLAMMPS:
    def test_read_lattice_data(self):
        s = Structure(data_path('water.data'))
        np.testing.assert_array_almost_equal(s.lattice(), LATTICE_CUBIC_10)

    def test_read_coord_data(self):
        s = Structure(data_path('water.data'))
        coord = s.coord()
        assert 'O' in coord.index
        assert 'H' in coord.index
        assert coord.loc[['O']].shape[0] == 1
        assert coord.loc[['H']].shape[0] == 2

    def test_triclinic_lattice(self):
        """LAMMPS data with xy tilt factor."""
        s = Structure(data_path('triclinic.data'))
        lattice = s.lattice()
        assert lattice[1, 0] != 0  # xy tilt is non-zero


class TestExtend:
    LATTICE_9 = [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0]

    def test_extend_2x2x2(self):
        s = Structure(data_path('water.xyz'), lattice_params=self.LATTICE_9)
        ext_lattice, ext_coord = s.extend(2, 2, 2)
        assert len(ext_coord) == 3 * 8  # 3 atoms * 2*2*2
        np.testing.assert_array_equal(ext_lattice, LATTICE_CUBIC_10 * [2, 2, 2])

    def test_extend_1x1x1(self):
        s = Structure(data_path('water.xyz'), lattice_params=self.LATTICE_9)
        ext_lattice, ext_coord = s.extend(1, 1, 1)
        assert len(ext_coord) == 3  # same as original


class TestMove:
    LATTICE_9 = [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0]

    def test_move_single_atom(self):
        s = Structure(data_path('water.xyz'), lattice_params=self.LATTICE_9)
        moved = s.move([1], [1.0, 0.0, 0.0])
        orig = s.coord()
        assert moved.iloc[0]['x'] == orig.iloc[0]['x'] + 1.0

    def test_move_direct(self):
        s = Structure(data_path('water.xyz'), lattice_params=self.LATTICE_9)
        moved = s.move([1], [0.1, 0.0, 0.0], direct=True)
        orig = s.coord()
        # 0.1 fractional * 10.0 lattice = 1.0 Angstrom
        np.testing.assert_almost_equal(moved.iloc[0]['x'], orig.iloc[0]['x'] + 1.0)
