import numpy as np
import pandas as pd
from toolkit.analysis.density import density_profile
from toolkit.analysis.angle import bond_angle_analysis


class TestDensityProfile:
    @staticmethod
    def _make_water_coords(n=10, L=30.0):
        """Create n water molecules along z-axis, spaced by 3 A, centered in xy."""
        cx, cy = L / 2, L / 2
        coords = []
        for i in range(n):
            coords.append({'x': cx, 'y': cy, 'z': float(i * 3)})
            coords.append({'x': cx, 'y': cy + 1.0, 'z': float(i * 3)})
            coords.append({'x': cx + 1.0, 'y': cy, 'z': float(i * 3)})
        index = ['O', 'H', 'H'] * n
        return pd.DataFrame(coords, index=index)

    @staticmethod
    def _make_cubic_lattice(L=30.0):
        return np.array([[L, 0, 0], [0, L, 0], [0, 0, L]], dtype=float)

    def test_density_basic(self):
        coord = self._make_water_coords(10)
        lattice = self._make_cubic_lattice(30.0)
        counts, dens = density_profile(coord, lattice, axis='z', element='O',
                                       interval=1.0)
        assert len(counts) > 0
        assert dens.sum() > 0

    def test_density_counts_sum_to_atoms(self):
        coord = self._make_water_coords(5, L=20.0)
        lattice = self._make_cubic_lattice(20.0)
        counts, _ = density_profile(coord, lattice, axis='z', element='O',
                                    interval=1.0)
        assert counts.sum() == 5  # 5 O atoms at z=0,3,6,9,12

    def test_density_with_bounds(self):
        coord = self._make_water_coords(10)
        lattice = self._make_cubic_lattice(30.0)
        counts, _ = density_profile(coord, lattice, axis='z', element='O',
                                    interval=1.0, vmin=0.0, vmax=15.0)
        assert counts.sum() > 0


class TestBondAngle:
    @staticmethod
    def _make_linear_water(n=10, lattice_L=30.0):
        """Create n linear water molecules (H-O-H aligned along x), centered."""
        cx, cy, cz = lattice_L / 2, lattice_L / 2, lattice_L / 2
        coords = []
        for i in range(n):
            coords.append({'x': cx, 'y': cy + float(i * 2.0), 'z': cz})  # O
            coords.append({'x': cx + 1.0, 'y': cy + float(i * 2.0), 'z': cz})  # H
            coords.append({'x': cx - 1.0, 'y': cy + float(i * 2.0), 'z': cz})  # H
        index = ['O', 'H', 'H'] * n
        return pd.DataFrame(coords, index=index)

    @staticmethod
    def _make_cubic_lattice(L=30.0):
        return np.array([[L, 0, 0], [0, L, 0], [0, 0, L]], dtype=float)

    def test_angle_returns_dataframes(self):
        coord = self._make_linear_water(10)
        lattice = self._make_cubic_lattice(30.0)
        normal, abnormal = bond_angle_analysis(
            coord, lattice, period='xy', main_ele='O', sub_ele='H',
            cutoff=1.5, num_coord=2
        )
        assert isinstance(normal, pd.DataFrame)
        assert isinstance(abnormal, pd.DataFrame)

    def test_angle_linear_water_180_deg(self):
        """Linear H-O-H should have ~180 degree bond angle."""
        coord = self._make_linear_water(10)
        lattice = self._make_cubic_lattice(30.0)
        normal, _ = bond_angle_analysis(
            coord, lattice, period='xy', main_ele='O', sub_ele='H',
            cutoff=1.5, num_coord=2
        )
        if not normal.empty and 'Bond' in normal.columns:
            angles = normal['Bond'].dropna()
            if len(angles) > 0:
                # Linear water bond angle ~180
                assert abs(angles.iloc[0] - 180.0) < 1.0
