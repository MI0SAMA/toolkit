import numpy as np
import pandas as pd
from scipy import spatial
import itertools


def _vec_angle(v1, v2):
    cos_a = np.inner(v1, v2)
    sin_a = np.linalg.norm(v1) * np.linalg.norm(v2)
    return np.degrees(np.arccos(cos_a / sin_a))


def _surface_angle(angle):
    return 90.0 - angle


def _half_coord(pt, period, lattice):
    mapping = {'xy': (0.5, 0.5, 1), 'yz': (1, 0.5, 0.5),
               'xz': (0.5, 1, 0.5), 'xyz': (0.5, 0.5, 0.5)}
    fcoord = np.dot(pt, np.linalg.inv(lattice))
    fcoord = fcoord * mapping[period]
    return np.dot(fcoord, lattice)


def bond_angle_analysis(coord, lattice, period, main_ele='O', sub_ele='H',
                        cutoff=1.1, num_coord=2):
    """Analyze bond angles and surface orientations.

    Args:
        coord: pd.DataFrame with columns ['x','y','z'], index=element labels
        lattice: 3x3 np.ndarray
        period: 'xy', 'yz', 'xz', or 'xyz'
        main_ele: central element label
        sub_ele: bonded element label
        cutoff: distance cutoff for bonds (Angstrom)
        num_coord: expected coordination number

    Returns:
        info_normal: pd.DataFrame of normal (matching coordination) sites
        info_abnormal: pd.DataFrame of abnormal (mismatched) sites
    """
    ext_map = {'xy': (3, 3, 1), 'yz': (1, 3, 3), 'xz': (3, 1, 3), 'xyz': (3, 3, 3)}
    nx, ny, nz = ext_map[period]

    org_coord = coord.loc[[main_ele, sub_ele]]
    ext_lattice = lattice * [nx, ny, nz]

    ranges = [list(range(1, n + 1)) for n in (nx, ny, nz)]
    extend_list = np.array(list(itertools.product(*ranges)))

    ext_coord = pd.DataFrame()
    for idx in extend_list:
        shifted = org_coord + np.dot(idx - [1, 1, 1], lattice)
        ext_coord = pd.concat([ext_coord, shifted])

    ext_coord_f = ext_coord.dot(np.linalg.inv(ext_lattice))
    ext_coord_f.columns = ['x', 'y', 'z']

    main_f = ext_coord_f.loc[main_ele]
    sub_f = ext_coord_f.loc[sub_ele]

    # Filter to central region
    main_mask = ((main_f['x'] > 1/3) & (main_f['x'] < 2/3) &
                 (main_f['y'] > 1/3) & (main_f['y'] < 2/3))
    sub_mask = ((sub_f['x'] > 0.3) & (sub_f['x'] < 0.7) &
                (sub_f['y'] > 0.3) & (sub_f['y'] < 0.7))

    main_cart = np.dot(main_f[main_mask].to_numpy(), ext_lattice)
    sub_cart = np.dot(sub_f[sub_mask].to_numpy(), ext_lattice)

    normal_records = {'Distance1': [], 'Distance2': [], 'Bond': [],
                      'Bond_Surface': [], 'Mid_Surface': [], 'Main_position': []}
    abnormal_records = {'Main_position': [], 'Sub_position': []}

    surface_vec = {'xy': [0, 0, 1], 'yz': [0, 1, 0], 'xz': [1, 0, 0]}

    for i in main_cart:
        dist = spatial.distance.cdist(sub_cart, [i], 'euclidean')
        dist_df = pd.DataFrame(dist, columns=['distance'])
        close = dist_df[dist_df['distance'] < cutoff]

        if len(close) == num_coord:
            idx0, idx1 = close.index[0], close.index[1]
            normal_records['Distance1'].append(float(close.iloc[0]))
            normal_records['Distance2'].append(float(close.iloc[1]))

            a1 = _vec_angle(sub_cart[idx0] - i, sub_cart[idx1] - i)
            normal_records['Bond'].append(a1)

            if period != 'xyz':
                sv = surface_vec[period]
                a2 = _surface_angle(_vec_angle(sub_cart[idx0] - i, sv))
                normal_records['Bond_Surface'].append(a2)
                a3 = _surface_angle(
                    _vec_angle((sub_cart[idx0] + sub_cart[idx1]) / 2 - i, sv)
                )
                normal_records['Mid_Surface'].append(a3)

            half = _half_coord(i, period, ext_lattice)
            normal_records['Main_position'].append(half)
        else:
            half = _half_coord(i, period, ext_lattice)
            abnormal_records['Main_position'].append(half)
            saved = []
            for j in sub_cart[close.index]:
                saved.append(_half_coord(j, period, ext_lattice))
            abnormal_records['Sub_position'].append(saved)

        sub_cart = np.delete(sub_cart, close.index, axis=0)

    info_normal = pd.DataFrame(normal_records)
    info_abnormal = pd.DataFrame(abnormal_records)
    if not info_normal.empty:
        info_normal.iloc[:, -1] = info_normal.iloc[:, -1].map(lambda x: x[0])
    if not info_abnormal.empty:
        info_abnormal.iloc[:, 0] = info_abnormal.iloc[:, 0].map(lambda x: x[0])
        info_abnormal.iloc[:, 1] = info_abnormal.iloc[:, 1].map(lambda x: x[0])
    return info_normal, info_abnormal
