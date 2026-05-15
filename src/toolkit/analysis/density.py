import numpy as np
import pandas as pd


def density_profile(coord, lattice, axis='z', element=None,
                    interval=0.2, vmin=None, vmax=None):
    """Calculate density profile along a lattice axis.

    Args:
        coord: pd.DataFrame with columns ['x','y','z'], index=element labels
        lattice: 3x3 np.ndarray
        axis: 'x', 'y', or 'z'
        element: element label to filter (None = all atoms)
        interval: bin width
        vmin, vmax: range bounds (None = auto from data)

    Returns:
        counts: count in each bin (pd.Series)
        density: number density per area per interval (pd.Series)
    """
    if element is not None:
        coord_subset = coord.loc[[element], axis].sort_values()
    else:
        coord_subset = coord[axis].sort_values()

    if vmin is None:
        vmin = np.floor(coord_subset.min())
    if vmax is None:
        vmax = np.ceil(coord_subset.max()) + interval

    bins = np.arange(vmin, vmax, interval)
    bin_labels = [
        f'{round(b, 3)}~{round(b + interval, 3)}' for b in bins
    ]
    counts = pd.Series(0.0, index=bin_labels)
    for j, bin_start in enumerate(bins):
        bin_end = bin_start + interval
        mask = (coord_subset >= bin_start) & (coord_subset < bin_end)
        counts.iloc[j] = mask.sum()

    axis_idx = {'x': [1, 2], 'y': [0, 2], 'z': [0, 1]}
    a, b = axis_idx[axis]
    area = np.linalg.norm(np.cross(lattice[a], lattice[b]))
    density = counts / area / interval
    return counts, density
