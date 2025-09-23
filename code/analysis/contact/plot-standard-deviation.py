# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from matplotlib.colors import LinearSegmentedColormap

# === Parameters ===
DISTANCE_CUTOFF = 9.5
LOWER_CUTOFF = 4.0  # Lower cutoff value (Å)

# === Residue information ===
TREM2_LEN = 46
DAP12_LEN = 66
N_RESIDUES = TREM2_LEN + DAP12_LEN

# === Custom gradient colormap: soft pink → soft blue → deep blue ===
custom_cmap = LinearSegmentedColormap.from_list(
    "pink_blue_deepblue", [
        "#fbb4c4",  # Soft pink
        "#b3cde3",  # Soft blue
        "#1c4e80"   # Deep blue
    ]
)

# === Heatmap plotting function ===
def plot_std_heatmap_gradient(matrix, tag, out_dir):
    flipped_matrix = np.flipud(matrix)  # Flip vertically

    plt.figure(figsize=(10, 8))
    masked = np.ma.masked_invalid(flipped_matrix)

    extent = [0, N_RESIDUES, 0, N_RESIDUES]
    im = plt.imshow(masked, cmap=custom_cmap, interpolation='nearest', extent=extent, aspect='equal', vmin=0, vmax=5)

    tick_interval = 20
    tick_positions = np.arange(0, N_RESIDUES + 1, tick_interval)
    tick_labels = [str(i) for i in tick_positions]

    plt.xticks(tick_positions, tick_labels, fontsize=18)
    plt.yticks(tick_positions, tick_labels, fontsize=18)

    plt.xlabel("Residue Index", fontsize=18, fontweight='bold')
    plt.ylabel("Residue Index", fontsize=18, fontweight='bold')

    # Dashed dividing lines
    plt.axvline(x=46, color='black', linestyle='--', linewidth=1.5)
    plt.axhline(y=46, color='black', linestyle='--', linewidth=1.5)
    plt.axvline(x=79, color='black', linestyle='--', linewidth=1.5)
    plt.axhline(y=79, color='black', linestyle='--', linewidth=1.5)

    # Add colorbar and set ticks
    cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
    cbar.set_label('STD value', fontsize=16)
    cbar.set_ticks(np.arange(0, 5.1, 1))
    cbar.ax.tick_params(labelsize=14)

    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"distance_heatmap-{tag}.tif")
    plt.savefig(out_file, dpi=300, format='tiff')
    plt.close()
    print(f"✅ Saved heatmap: {out_file}")

# === Main program ===
std_files = glob.glob("std_matrix_last25000-*.txt")
output_dir = 'std_singlecolor'
os.makedirs(output_dir, exist_ok=True)

for std_file in std_files:
    tag = os.path.splitext(os.path.basename(std_file))[0].replace("std_matrix_last25000-", "")
    avg_file = std_file.replace("std_matrix", "avg_matrix")

    if not os.path.exists(avg_file):
        print(f"⚠️ Warning: Avg matrix file not found for {tag}, skipping.")
        continue

    try:
        std_matrix = np.loadtxt(std_file, delimiter='\t')
        avg_matrix = np.loadtxt(avg_file, delimiter='\t')

        if std_matrix.shape != (N_RESIDUES, N_RESIDUES):
            raise ValueError(f"Invalid shape in STD matrix for {tag}")
        if avg_matrix.shape != (N_RESIDUES, N_RESIDUES):
            raise ValueError(f"Invalid shape in AVG matrix for {tag}")

        print(f"📊 Processing {tag}...")

        mask = (avg_matrix > LOWER_CUTOFF) & (avg_matrix <= DISTANCE_CUTOFF)
        filtered_std_matrix = np.where(mask, std_matrix, np.nan)

        # Save Excel
        excel_path = os.path.join(output_dir, f"{tag}_std_matrix_filtered.xlsx")
        df_matrix = pd.DataFrame(filtered_std_matrix,
                                 index=range(1, N_RESIDUES + 1),
                                 columns=range(1, N_RESIDUES + 1)).fillna(0)
        df_matrix.to_excel(excel_path)
        print(f"💾 Saved Excel: {excel_path}")

        # Plot heatmap
        plot_std_heatmap_gradient(filtered_std_matrix, tag, output_dir)

    except Exception as e:
        print(f"❌ Error processing {tag}: {e}")

# === Merge Excel files ===
merged_excel_path = os.path.join(output_dir, "merged_std_matrix_filtered.xlsx")

groups = {
    'A': ['A1', 'AT1', 'AT2', 'AT3', 'AT4', 'AT5'],
    'B': ['B1', 'BT1', 'BT2', 'BT3'],
    'C': ['C1', 'CT1', 'CT2', 'CT3'],
    'D': ['D1', 'DT1', 'DT2', 'DT3']
}

all_tags = sum(groups.values(), [])
first_tag = all_tags[0]
first_file = os.path.join(output_dir, f"{first_tag}_std_matrix_filtered.xlsx")
df_first = pd.read_excel(first_file, index_col=0)

merged_data = []
for tag in all_tags:
    file_path = os.path.join(output_dir, f"{tag}_std_matrix_filtered.xlsx")
    if not os.path.exists(file_path):
        print(f"⚠️ Missing file for {tag}, filling NaN.")
        merged_data.append([np.nan] * (N_RESIDUES * N_RESIDUES))
    else:
        df = pd.read_excel(file_path, index_col=0)
        flattened = df.values.flatten()
        merged_data.append(flattened.tolist())

merged_df = pd.DataFrame(merged_data, columns=[f"{i}-{j}" for i in df_first.index for j in df_first.columns])
merged_df.insert(0, 'Sample', all_tags)

with pd.ExcelWriter(merged_excel_path, engine='openpyxl') as writer:
    header = ['Sample'] + [f"{i}-{j}" for i in df_first.index for j in df_first.columns]
    header_df = pd.DataFrame([header])
    header_df.to_excel(writer, index=False, header=False, startrow=0)
    merged_df.to_excel(writer, index=False, header=False, startrow=1)

print(f"📦 Merged matrix saved to: {merged_excel_path}")

