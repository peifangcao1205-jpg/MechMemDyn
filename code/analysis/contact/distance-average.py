# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from matplotlib.colors import ListedColormap

# === Parameters ===
DISTANCE_CUTOFF = 9.5
LOWER_CUTOFF = 4.0  # Added lower cutoff value (Å)

# === Residue information ===
TREM2_LEN = 46
DAP12_LEN = 66
N_RESIDUES = TREM2_LEN + DAP12_LEN

# === Visualization function ===
def plot_submatrix_heatmap(matrix, tag):
    # Generate binary matrix: 3 Å < distance ≤ 9.5 Å set to 1, others set to NaN
    binary_matrix = np.where((matrix > LOWER_CUTOFF) & (matrix <= DISTANCE_CUTOFF), 1, np.nan)
    flipped_matrix = np.flipud(binary_matrix)

    plt.figure(figsize=(10, 8))

    # Use single-color colormap (gray), NaN displayed as white
    cmap = ListedColormap(['gray'])
    masked = np.ma.masked_invalid(flipped_matrix)

    # Plot (no gradient color, no colorbar)
    extent = [0, N_RESIDUES, 0, N_RESIDUES]
    plt.imshow(masked, cmap=cmap, interpolation='nearest', extent=extent, aspect='equal')
    
    # Axis tick settings
    tick_interval = 20
    tick_positions = np.arange(0, N_RESIDUES, tick_interval)
    tick_labels = [str(i) for i in tick_positions]

    plt.xticks(tick_positions, tick_labels, fontsize=18)
    plt.yticks(tick_positions, tick_labels, fontsize=18)

    plt.xlabel("Residue Index", fontsize=18, fontweight='bold')
    plt.ylabel("Residue Index", fontsize=18, fontweight='bold')

    # Add dashed dividing lines
    plt.axvline(x=46, color='black', linestyle='--', linewidth=1.5)
    plt.axhline(y=46, color='black', linestyle='--', linewidth=1.5)
    plt.axvline(x=79, color='black', linestyle='--', linewidth=1.5)
    plt.axhline(y=79, color='black', linestyle='--', linewidth=1.5)

    plt.tight_layout()
    out_file = os.path.join('average', f"distance_heatmap-{tag}.tif")
    plt.savefig(out_file, dpi=300, format='tiff')
    plt.close()
    print(f"Saved heatmap: {out_file}")

# === Main program ===
matrix_files = glob.glob("avg_matrix_last25000-*.txt")
target_files = matrix_files  # Removed "protein" not in f filter, include all files

for file in target_files:
    tag = os.path.splitext(os.path.basename(file))[0].replace("avg_matrix_last25000-", "")
    try:
        matrix = np.loadtxt(file, delimiter='\t')
        if matrix.shape != (N_RESIDUES, N_RESIDUES):
            raise ValueError(f"Invalid matrix shape {matrix.shape}, expected {N_RESIDUES}x{N_RESIDUES}")
        print(f"Loaded matrix for {tag}, shape: {matrix.shape}")

        # Filter distances: keep values ≤ DISTANCE_CUTOFF, set others to 0
        filtered_matrix = np.where(matrix <= DISTANCE_CUTOFF, matrix, 0)

        # Save as Excel file consistent with A1_contact_matrix.xlsx
        output_dir = 'average'
        os.makedirs(output_dir, exist_ok=True)
        df_matrix = pd.DataFrame(filtered_matrix, index=range(1, N_RESIDUES + 1), columns=range(1, N_RESIDUES + 1))
        df_matrix.to_excel(os.path.join(output_dir, f"{tag}_contact_matrix.xlsx"))

        # Continue generating heatmap
        plot_submatrix_heatmap(matrix, tag)
    except Exception as e:
        print(f"Error processing {file}: {e}")

