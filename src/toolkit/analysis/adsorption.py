import numpy as np
from scipy.spatial import Delaunay


def get_surface_by_range(atoms, z_range):
    """Get surface atom indices within a z-range from ASE Atoms object."""
    pos = atoms.get_positions()
    lo, hi = min(z_range), max(z_range)
    return [i for i in range(len(atoms)) if lo < pos[i][2] < hi]


def get_extended_atoms(atoms):
    """Extend surface 3x3x1 in xy for periodic analysis."""
    tmp = atoms.copy()
    tmp = atoms * [3, 3, 1]
    cell = atoms.get_cell()
    tmp.set_positions(tmp.get_positions() - cell[0] - cell[1])
    return tmp


def filter_outside_cell(points, cell):
    """Remove points whose fractional coords are outside [0,1]."""
    fcoords = np.dot(points, np.linalg.inv(cell))
    mask = (fcoords >= 0).all(axis=1) & (fcoords <= 1).all(axis=1)
    return points[mask]


def find_adsorption_sites(surface_atoms):
    """Find top, bridge, and hollow sites from surface atom positions.

    Args:
        surface_atoms: ASE Atoms object of surface layer

    Returns:
        top_sites, hollow_sites, bridge_sites: np.ndarrays of positions
    """
    surf_ext = get_extended_atoms(surface_atoms)
    pos_ext = surf_ext.get_positions()
    top_sites = surface_atoms.get_positions()
    cell = surface_atoms.get_cell()

    tri = Delaunay(pos_ext[:, :2])
    nodes = pos_ext[tri.simplices]

    hollow_sites = np.array([t.sum(axis=0) / 3 for t in nodes])
    bridge_sites = []
    for t in nodes:
        bridge_sites.append((t[0] + t[1]) / 2)
        bridge_sites.append((t[0] + t[2]) / 2)
        bridge_sites.append((t[1] + t[2]) / 2)
    bridge_sites = np.array(bridge_sites)

    hollow_sites = filter_outside_cell(hollow_sites, cell)
    bridge_sites = filter_outside_cell(bridge_sites, cell)

    return top_sites, hollow_sites, bridge_sites
