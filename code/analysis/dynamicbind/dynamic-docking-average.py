import MDAnalysis as mda
import numpy as np
import os
import matplotlib.pyplot as plt
from MDAnalysis.lib.distances import distance_array

# Parameters
DISTANCE_CUTOFF = 9.5  # Cutoff threshold
FOLDER = "A6"  # Folder path
OUTPUT_AVG = os.path.join(FOLDER, "avg_matrix_dynamic.txt")
OUTPUT_DIFF = os.path.join(FOLDER, "diff_matrix.txt")
OUTPUT_DIFF_HEATMAP = os.path.join(FOLDER, "diff_matrix_heatmap.tiff")

# =============== Functions ==================

def compute_residue_distance_matrix(universe, n_residues=112):
    """Compute minimum distance matrix for protein residues (n_residues x n_residues)"""
    all_residues = universe.select_atoms(f'protein and resid 1:{n_residues}')
    if len(all_residues.residues) != n_residues:
        raise ValueError(f"Number of residues in PDB {universe.filename} ({len(all_residues.residues)}) != {n_residues}")
    matrix = np.full((n_residues, n_residues), np.inf)

    for i, res1 in enumerate(all_residues.residues):
        for j, res2 in enumerate(all_residues.residues):
            if i <= j:
                pos1 = res1.atoms.positions
                pos2 = res2.atoms.positions
                dist = distance_array(pos1, pos2)
                min_dist = np.min(dist)
                matrix[i, j] = min_dist
                matrix[j, i] = min_dist
    return matrix

def save_matrix_txt_only_values(matrix, filename):
    """Save numeric matrix only (without row/column labels)"""
    np.savetxt(filename, matrix, fmt="%.3f", delimiter="\t")
    print(f"Matrix saved to {filename} (values only, shape {matrix.shape[0]}x{matrix.shape[1]})")

def plot_diff_heatmap(matrix, filename, n_residues=112):
    """Plot heatmap of difference matrix"""
    plt.figure(figsize=(10, 8))
    vmax = np.max(np.abs(matrix))
    im = plt.imshow(
        np.flipud(matrix),
        cmap="bwr",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
        extent=[0, n_residues, 0, n_residues],
        aspect="equal"
    )

    tick_interval = 20
    tick_positions = np.arange(0, n_residues+1, tick_interval)
    tick_labels = [str(i) for i in tick_positions]
    plt.xticks(tick_positions, tick_labels, fontsize=16)
    plt.yticks(tick_positions, tick_labels, fontsize=16)
    plt.xlabel("Residue Index", fontsize=18, fontweight="bold")
    plt.ylabel("Residue Index", fontsize=18, fontweight="bold")

    # Divider lines (adjust based on your regions)
    plt.axvline(x=46, color='black', linestyle='--', linewidth=1.5)
    plt.axhline(y=46, color='black', linestyle='--', linewidth=1.5)
    plt.axvline(x=79, color='black', linestyle='--', linewidth=1.5)
    plt.axhline(y=79, color='black', linestyle='--', linewidth=1.5)

    # === Colorbar settings ===
    cbar = plt.colorbar(im, orientation='vertical')
    cbar.set_label("Distance Change (Å)", fontsize=18, labelpad=10)
    cbar.set_ticks([-5, 0, 5])
    cbar.ax.tick_params(labelsize=18)

    plt.tight_layout()
    plt.savefig(filename, format="tiff", dpi=300)
    plt.close()
    print(f"Diff heatmap saved to {filename}")

# =============== Main program ==================

def main():
    # Read rank*.pdb files
    pdb_files = sorted([os.path.join(FOLDER, f) for f in os.listdir(FOLDER) if f.startswith("rank") and f.endswith(".pdb")])
    if not pdb_files:
        raise FileNotFoundError(f"No rank*.pdb files found in {FOLDER}")

    print(f"Found {len(pdb_files)} PDB files")

    # Compute all rank matrices
    matrices = []
    for pdb in pdb_files:
        u = mda.Universe(pdb)
        matrix = compute_residue_distance_matrix(u)
        # Truncate distances greater than DISTANCE_CUTOFF
        matrix = np.where(matrix > DISTANCE_CUTOFF, DISTANCE_CUTOFF, matrix)
        matrices.append(matrix)
        print(f"Processed {pdb}")

    # Average matrix
    avg_matrix = np.mean(matrices, axis=0)
    save_matrix_txt_only_values(avg_matrix, OUTPUT_AVG)

    # Initial docking matrix
    docking_matrix = np.loadtxt(os.path.join(FOLDER, "avg_matrix_docking.txt"))
    if docking_matrix.shape != avg_matrix.shape:
        raise ValueError(f"Initial matrix {docking_matrix.shape} does not match average matrix {avg_matrix.shape}")
    docking_matrix = np.where(docking_matrix > DISTANCE_CUTOFF, DISTANCE_CUTOFF, docking_matrix)

    # Difference matrix
    diff_matrix = avg_matrix - docking_matrix
    save_matrix_txt_only_values(diff_matrix, OUTPUT_DIFF)

    # Visualize difference matrix
    plot_diff_heatmap(diff_matrix, OUTPUT_DIFF_HEATMAP, n_residues=diff_matrix.shape[0])

if __name__ == "__main__":
    main()

