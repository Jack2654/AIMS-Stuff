import numpy as np
import math
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from BasicBandOut import band_info


# made python file 01/12/2026, most updated angle information


# returns set of atomic coordinates with auxiliary information removed
def atoms_trimmed(file):
    atoms = []
    with open(file, "r") as f:
        for line in f:
            if "atom" in line and "#" not in line:
                atoms.append(line.strip())
    return [[float(x) for x in atom.split()[1:4]] for atom in atoms]


# returns set of lattice vectors
def lattice_vectors(file):
    lat_vecs = []
    with open(file, "r") as f:
        for line in f:
            if "lattice" in line:
                lat_vecs.append(line)
    return [[float(x) for x in lat.split()[1:4]] for lat in lat_vecs]


# computes angle in degrees, assumes np.array() inputs
def angle(pt1, center, pt2):
    ba = pt1 - center
    bc = pt2 - center
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    result = np.arccos(cosine_angle)
    return np.degrees(result)


# returns shifted set of separate bonds as defined by input parameters
def pull_data(readfile, atom_idxs, shiftmap):
    if len(atom_idxs) != 8 or len(shiftmap) != 8:
        raise ValueError("wrong dimension atom_idxs or shiftmap passed to pull_data")

    # pulls lattice vectors and atoms from the geometry.in file
    lv = lattice_vectors(readfile)
    at = atoms_trimmed(readfile)
    atom_positions = [at[int(atom_idx) - 1] for atom_idx in atom_idxs]

    # shifts the atoms as necessary
    shifted_pos = [[coord for coord in atom_pos] for atom_pos in atom_positions]
    for k, shift in enumerate(shiftmap):
        if shift is not None:
            for i in range(3):
                for j in range(3):
                    shifted_pos[k][i] += shift[j] * lv[j][i]
    return shifted_pos


# sets z=0 for metal #1, rotates the plane so z=0 for metal #3
# finally rotates once to minimize the average height of metal #2/#4
def align_atoms(atom_pos):
    recenter = np.array(atom_pos[0])
    shifted_pos = [np.array(pos) - recenter for pos in atom_pos]

    axis = shifted_pos[2] - shifted_pos[6]
    axis = axis / np.linalg.norm(axis)
    cur_flat_dist = math.sqrt(shifted_pos[4][0] ** 2 + shifted_pos[4][1] ** 2)
    cur_vert_dist = shifted_pos[4][2]
    cur_angle = np.arctan(cur_vert_dist / cur_flat_dist)
    rotation = R.from_rotvec(-1 * cur_angle * axis)  # -1 to "undo" the current rotation
    for i in range(len(shifted_pos)):
        shifted_pos[i] = rotation.apply(shifted_pos[i])

    axis = shifted_pos[4] - shifted_pos[0]
    axis = axis / np.linalg.norm(axis)
    cur_flat_dist_1 = math.sqrt(shifted_pos[2][0] ** 2 + shifted_pos[2][1] ** 2)
    cur_vert_dist_1 = shifted_pos[2][2]
    cur_flat_dist_2 = math.sqrt(shifted_pos[6][0] ** 2 + shifted_pos[6][1] ** 2)
    cur_vert_dist_2 = shifted_pos[6][2]
    a1 = (cur_vert_dist_1 / cur_flat_dist_1)
    a2 = (cur_vert_dist_2 / cur_flat_dist_2)
    cur_angle = np.arctan((a1 - a2) / 2)
    rotation = R.from_rotvec(cur_angle * axis)
    for i in range(len(shifted_pos)):
        shifted_pos[i] = rotation.apply(shifted_pos[i])

    return shifted_pos


# writes list of atom positions as FHI-aims geometry.in file format for easy testing
def write_as_aims(atom_pos):
    print("lattice_vector 20 0 0")
    print("lattice_vector 0 20 0")
    print("lattice_vector 0 0 20")
    for atom in atom_pos:
        print(f'atom {" ".join([str(x) for x in atom])} H')


# takes set of 8 positions and formats as easy-to-iterate bonds
def partition_to_bonds(atom_pos):
    bonds = []
    for i in range(4):
        cur_coords = atom_pos[2 * i: 2 * i + 3]
        if i == 3:
            cur_coords.append(atom_pos[0])
        bonds.append(cur_coords)
    return bonds


# from bonds compute 4 beta_ins by projecting and applying reference frame
# even bonds (0 and 2) should be inside and
# odd bonds (1 and 3) should be outside
def compute_beta_ins(bonds):
    flat_atoms = [[np.array([atom[0], atom[1], 0]) for atom in bond] for bond in bonds]
    centroid = (flat_atoms[0][0] + flat_atoms[1][0] + flat_atoms[2][0] + flat_atoms[3][0]) / 4
    beta_in = []
    for count, bond in enumerate(flat_atoms):
        m = (bond[0][1] - bond[2][1]) / (bond[0][0] - bond[2][0])
        b = bond[0][1] - m * bond[0][0]
        x_proj = (bond[1][0] + m * bond[1][1] - m * b) / (1 + m * m)
        y_proj = m * x_proj + b
        dist_actual = np.linalg.norm(bond[1] - centroid)
        dist_proj = np.linalg.norm(np.array([x_proj, y_proj, 0]) - centroid)

        cur_angle = angle(*bond)
        if dist_actual <= dist_proj:
            beta_in.append(cur_angle if (count % 2) == 0 else 360 - cur_angle)
        else:
            beta_in.append(cur_angle if (count % 2) == 1 else 360 - cur_angle)
    return beta_in


# from bonds compute 4 beta_ins by projecting and applying reference frame
# even bonds (0 and 1) should be down and
# odd bonds (2 and 3) should be up
def compute_beta_outs(bonds):
    # construct a plane connecting endpoints and parallel to (0, 0, 1)
    def construct_plane_normal(metal_1, metal_2):
        vec_thru_metals = np.array(metal_2) - np.array(metal_1)
        n = np.cross(vec_thru_metals, np.array([0, 0, 1]))
        return n

    # project a point onto a given plane
    def project_point_to_plane(p, p0, n):
        w = np.array(p) - np.array(p0)
        scalar = np.dot(n, w) / np.dot(n, n)
        proj = np.array(p) - scalar * n
        return proj

    beta_out = []
    for count, bond in enumerate(bonds):
        p1 = np.array(bond[0])
        p2 = np.array(bond[2])
        p_to_project = np.array(bond[1])
        normal = construct_plane_normal(p1, p2)
        p_proj = project_point_to_plane(p_to_project, p1, normal)
        beta_out_cur = angle(p1, p_proj, p2)

        new_normal = np.cross(np.array([1, 1, 0]), p2 - p1)
        new_normal = new_normal if new_normal[2] >= 0 else -1 * new_normal
        side = np.dot(new_normal, p_proj - p1)

        if side <= 0:
            beta_out_cur = beta_out_cur if count < 2 else 360 - beta_out_cur
        else:
            beta_out_cur = beta_out_cur if count >= 2 else 360 - beta_out_cur
        beta_out.append(beta_out_cur)

    return beta_out


# given any set of four ordered betas, compute and return two possible delta beta values
def compute_delta_betas(betas):
    db_1 = abs(betas[0] + betas[1] - betas[2] - betas[3]) / 2
    db_2 = abs(betas[0] + betas[3] - betas[2] - betas[1]) / 2
    return [db_1, db_2]


# computes:
#  - 4x beta -- no reference frame
#  - 4x beta_in -- make all atoms roughly coplanar (rotations), then project onto 2D plane and apply reference frame
#  - 4x beta_out -- same as beta_in but with different planes being used for projection
# Inputs:
#  - readfile: path to geometry.in to read for coordinates
#  - shiftmap: how to shift each atom (lv units) at same index to form one perfect copy of 8-atom "square"
#  - atom_idxs: list of atom indexes (1-indexed) as the following:
# 1-Ag  --  2-I  --  3-Bi
#  |                  |
# 8-I                4-I
#  |                  |
# 7-Bi  --  6-I  --  5-Ag
def new_delta_beta(readfile, atom_idxs, shiftmap):
    atom_pos = pull_data(readfile, atom_idxs, shiftmap)
    rotated_pos = align_atoms(atom_pos)
    if False and "IJAYUQ" in readfile:
        write_as_aims(rotated_pos)
    bonds = partition_to_bonds(rotated_pos)

    betas = [angle(*bond) for bond in bonds]
    beta_ins = compute_beta_ins(bonds)
    beta_outs = compute_beta_outs(bonds)

    delta_betas = compute_delta_betas(betas)
    delta_beta_ins = compute_delta_betas(beta_ins)
    delta_beta_outs = compute_delta_betas(beta_outs)
    return [delta_betas, delta_beta_ins, delta_beta_outs]


def plot_data(cur_x, cur_y, x_label):
    cur_x = np.array(cur_x)
    cur_y = np.array(cur_y)
    plt.clf()
    A = np.vstack([cur_x, np.ones_like(cur_x)]).T
    m, c = np.linalg.lstsq(A, cur_y, rcond=None)[0]
    y_pred = m * cur_x + c
    ss_res = np.sum((cur_y - y_pred) ** 2)
    ss_tot = np.sum((cur_y - np.mean(cur_y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    plt.scatter(cur_x, cur_y, color='r')
    plt.plot(cur_x, y_pred, 'k', label=f"$R^2$ = {r_squared:.3f}")
    plt.xlabel(x_label)
    plt.ylabel(r'$\Delta E^\pm$')
    plt.legend()
    plt.tight_layout()
    plt.show()


# plots for AgBi real structures
def AgBi_plots():
    # setup data
    all_idxs, all_shifts = {}, {}
    structures = ["VOXKIF/", "SELWOZ/", "SELWUF/", "IJAYUQ/"]
    all_idxs[structures[0]] = [2, 7, 1, 10, 2, 9, 1, 8]
    all_shifts[structures[0]] = [None, None, None, None, [0, 1, 0], None, [1, 0, 0], None]
    all_idxs[structures[1]] = [1, 2, 27, 4, 1, 7, 27, 3]
    all_shifts[structures[1]] = [None, None, None, [0, 1, 0], [0, 1, 0], [0, 1, 0], [-1, 0, 0], None]
    all_idxs[structures[2]] = [2, 9, 1, 3, 2, 6, 1, 8]
    all_shifts[structures[2]] = [None, None, None, None, [1, 0, 0], [0, -1, 0], [0, -1, 0], [0, -1, 0]]
    all_idxs[structures[3]] = [2, 10, 19, 3, 1, 11, 19, 6]
    all_shifts[structures[3]] = [None, None, None, None, [0, 0, -1], [0, 0, -1], [1, 0, 0], [0, 0, -1]]
    all_idxs["IJAYUQ_2"] = [20, 4, 2, 12, 20, 5, 1, 9]
    all_shifts["IJAYUQ_2"] = [[0, 0, -1], [0, 0, -1], None, None, [-1, 0, -1], None, [0, -1, -1], [0, 0, -1]]

    base_dir = "../../FHI-aims/Double_Perovskites/Real_Structures/input_geos/"
    db_res = []
    E_spin_splitting = []
    for structure in structures:
        folder = base_dir + "../calcs_and_outputs/bands/" + structure
        splits = [band_info(folder, f'band100{i}.out', band_gap=False, spin_splitting=True, verbose=False) for i in
                  range(1, 5)]
        E_spin_splitting.append(max([float(x) for x in splits]))
        cur_structure = base_dir + structure + "geometry.in"
        db_res.append(new_delta_beta(cur_structure, all_idxs[structure], all_shifts[structure]))
        if "IJAYUQ" in structure:
            addtl_db = new_delta_beta(cur_structure, all_idxs["IJAYUQ_2"], all_shifts["IJAYUQ_2"])
            for i in range(3):
                db_res[-1][i] += addtl_db[i]

    db_old_max = [max(structure[0]) for structure in db_res]
    db_old_avg = [sum(structure[0]) / len(structure[0]) for structure in db_res]
    db_in_max = [max(structure[1]) for structure in db_res]
    db_in_avg = [sum(structure[1]) / len(structure[1]) for structure in db_res]
    db_out_max = [max(structure[2]) for structure in db_res]
    db_out_avg = [sum(structure[2]) / len(structure[2]) for structure in db_res]

    db_new_max = [(db_in_max[i] + db_out_max[i]) / 2 for i in range(len(structures))]
    db_new_avg = [(db_in_avg[i] + db_out_avg[i]) / 2 for i in range(len(structures))]

    db_new_full_max = [max([db_in_max[i], db_out_max[i]]) for i in range(len(structures))]

    plot_data(db_old_max, E_spin_splitting, "DB old max")  # 0.590
    plot_data(db_old_avg, E_spin_splitting, "DB old avg")  # 0.862

    plot_data(db_in_max, E_spin_splitting, "DB in max")  # 0.973 !!!
    plot_data(db_in_avg, E_spin_splitting, "DB in avg")  # 0.921

    plot_data(db_out_max, E_spin_splitting, "DB out max")  # 0.271
    plot_data(db_out_avg, E_spin_splitting, "DB out avg")  # 0.557

    plot_data(db_new_max, E_spin_splitting, "DB new max")  # 0.914
    plot_data(db_new_avg, E_spin_splitting, "DB new avg")  # 0.875
    plot_data(db_new_full_max, E_spin_splitting, "DB full max")  # 0.973 !!!


AgBi_plots()
