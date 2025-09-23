# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# === Parameters ===
FILES = ['A1', 'AT1', 'AT2', 'AT3', 'AT4', 'AT5', 'B1', 'BT1', 'BT2', 'BT3', 
         'C1', 'CT1', 'CT2', 'CT3', 'D1', 'DT1', 'DT2', 'DT3']
MATRIX_PROT = 'avg_matrix_last25000-protein.txt'
DISTANCE_CUTOFF = 9.5                # Distance cutoff (Å)
DPI = 300                            # Image resolution
OUTPUT_DIR = 'difference'            # Output directory

# === Create output directory if it doesn't exist ===
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# === Residue information ===
TREM2_LEN = 46
DAP12_LEN = 66
N_RESIDUES = TREM2_LEN + DAP12_LEN

# === Process each file ===
for tag in FILES:
    matrix_file = f'avg_matrix_last25000-{tag}.txt'
    interface_file = os.path.join(OUTPUT_DIR, f'interface_pairs-{tag}.xlsx')
    enhanced_file = os.path.join(OUTPUT_DIR, f'enhanced_pairs-{tag}.xlsx')
    heatmap_file = os.path.join(OUTPUT_DIR, f'diff_distance_matrix-{tag}.tif')

    # Read target matrix
    try:
        matrix = np.loadtxt(matrix_file)
    except Exception as e:
        print(f"❌ Failed to read {matrix_file}: {e}")
        continue

    # Read protein matrix
    try:
        matrix_prot = np.loadtxt(MATRIX_PROT)
    except Exception as e:
        print(f"❌ Failed to read {MATRIX_PROT}: {e}")
        continue

    # Limit maximum value to 9.5 Å
    matrix = np.where(matrix > DISTANCE_CUTOFF, DISTANCE_CUTOFF, matrix)
    matrix_prot = np.where(matrix_prot > DISTANCE_CUTOFF, DISTANCE_CUTOFF, matrix_prot)

    # Difference matrix
    diff_matrix = matrix - matrix_prot

    # Mask values beyond ±9.5 Å
    masked_matrix = np.ma.masked_where(np.abs(diff_matrix) > DISTANCE_CUTOFF, diff_matrix)

    # Flip matrix for Y-axis
    flipped_matrix = np.flipud(masked_matrix)

    # === Plot individual heatmap ===
    plt.figure(figsize=(10, 8))  # Square figure
    cmap = plt.cm.seismic.copy()
    cmap.set_bad(color='white')

    # Ensure square plot
    extent = [0, N_RESIDUES, 0, N_RESIDUES]
    im = plt.imshow(flipped_matrix, cmap=cmap, vmin=-DISTANCE_CUTOFF, vmax=DISTANCE_CUTOFF, 
                    interpolation='nearest', extent=extent, aspect='equal')

    # Add colorbar
    cbar = plt.colorbar(im, orientation='vertical')
    cbar.set_label("Distance Change (Å)", fontsize=18, labelpad=10)
    cbar.ax.tick_params(labelsize=18)

    # Add black dashed lines at x=46, y=46, x=79, y=79
    plt.axvline(x=46, color='black', linestyle='--', linewidth=1.5)  # Vertical at res 46
    plt.axhline(y=46, color='black', linestyle='--', linewidth=1.5)  # Horizontal at res 46
    plt.axvline(x=79, color='black', linestyle='--', linewidth=1.5)  # Vertical at res 79
    plt.axhline(y=79, color='black', linestyle='--', linewidth=1.5)  # Horizontal at res 79
    # Set ticks
    tick_interval = 20
    tick_positions = np.arange(0, N_RESIDUES, tick_interval)
    tick_labels = [str(i) for i in tick_positions]

    plt.xticks(tick_positions, tick_labels, fontsize=18)
    plt.yticks(tick_positions, tick_labels, fontsize=18)

    # Set labels
    plt.xlabel("Residue Index", fontsize=18, fontweight='bold')
    plt.ylabel("Residue Index", fontsize=18, fontweight='bold')

    plt.tight_layout()
    plt.savefig(heatmap_file, dpi=DPI, format='tiff')
    plt.close()

    # === Statistics for residue pairs ===
    all_pairs = []

    total_residues = diff_matrix.shape[0]
    for i in range(total_residues):
        for j in range(i + 1, total_residues):
            delta = diff_matrix[i, j]
            region = "TREM2–TREM2" if i < TREM2_LEN and j < TREM2_LEN else \
                     "DAP12–DAP12" if i >= TREM2_LEN and j >= TREM2_LEN else \
                     "TREM2–DAP12"
            all_pairs.append([i, j, delta, region])

    # Use all pairs for interface_pairs (include TREM2–TREM2, TREM2–DAP12, DAP12–DAP12)
    df_all_pairs = pd.DataFrame(all_pairs, columns=["Res_i", "Res_j", "ΔDistance", "Region"])
    df_interface = df_all_pairs[["Res_i", "Res_j", "ΔDistance", "Region"]]  # Include all regions
    df_enhanced = df_all_pairs[["Res_i", "Res_j", "ΔDistance", "Region"]]  # Unchanged

    df_interface.to_excel(interface_file, index=False)
    df_enhanced.to_excel(enhanced_file, index=False)

    print(f"Distance difference matrix and pair data for {tag} processed. Heatmap saved to {heatmap_file}")

# === Output summary ===
print("\n✅ All processing completed. Files generated in 'difference' folder:")
for tag in FILES:
    print(f" - {os.path.join(OUTPUT_DIR, f'interface_pairs-{tag}.xlsx')}")
    print(f" - {os.path.join(OUTPUT_DIR, f'enhanced_pairs-{tag}.xlsx')}")
    print(f" - {os.path.join(OUTPUT_DIR, f'diff_distance_matrix-{tag}.tif')}")
