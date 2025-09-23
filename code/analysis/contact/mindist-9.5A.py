# -*- coding: utf-8 -*-
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.lib.distances import distance_array

# Parameters
DISTANCE_CUTOFF = 9.5
CONTACT_FILE = 'min_residue_distances.txt'          # Contact list
AVG_MATRIX_FILE = 'avg_matrix_last25000.txt'         # Average matrix
STD_MATRIX_FILE = 'std_matrix_last25000.txt'         # Standard deviation matrix

# Load structure and trajectory
u = mda.Universe('com_md.gro', 'com_md1000ns_noPBC.xtc')
n_frames = len(u.trajectory)

# Residue selection
protein1 = u.select_atoms('protein and resid 1:46')     # TREM2
protein2 = u.select_atoms('protein and resid 47:112')   # DAP12
all_residues = protein1.residues + protein2.residues
n_residues = len(all_residues)

print(f"TREM2 residues: {len(protein1.residues)}")
print(f"DAP12 residues: {len(protein2.residues)}")
print(f"Total frames: {n_frames}")

# Record minimum distance pairs for the last 2000 frames
min_distances = []
distance_matrices = []

for ts in u.trajectory[-25000:]:
    print(f"Analyzing frame {ts.frame}")
    frame_matrix = np.full((n_residues, n_residues), np.inf)

    for i, res1 in enumerate(all_residues):
        for j, res2 in enumerate(all_residues[i:], start=i):
            pos1 = res1.atoms.positions
            pos2 = res2.atoms.positions
            dist = distance_array(pos1, pos2)
            if i == j:
                np.fill_diagonal(dist, np.inf)
            min_dist = np.min(dist)
            frame_matrix[i, j] = min_dist
            frame_matrix[j, i] = min_dist
            if min_dist < DISTANCE_CUTOFF:
                min_distances.append((ts.frame, res1.resid, res2.resid, min_dist))

    distance_matrices.append(frame_matrix)

# Save contact list file
with open(CONTACT_FILE, 'w', encoding='utf-8') as f:
    f.write('Frame\tRes1\tRes2\tMin_Distance\n')
    for frame, r1, r2, d in min_distances:
        f.write(f"{frame}\t{r1}\t{r2}\t{d:.3f}\n")

print(f"Contact list saved to {CONTACT_FILE}")

# Calculate average distance matrix & standard deviation matrix
print("Calculating average and standard deviation matrices...")
distance_stack = np.stack(distance_matrices)  # shape: (2000, n_residues, n_residues)

avg_matrix = np.mean(distance_stack, axis=0)
np.savetxt(AVG_MATRIX_FILE, avg_matrix, fmt="%.3f", delimiter="\t")
print(f"Average matrix saved to {AVG_MATRIX_FILE}")

std_matrix = np.std(distance_stack, axis=0)
np.savetxt(STD_MATRIX_FILE, std_matrix, fmt="%.3f", delimiter="\t")
print(f"Standard deviation matrix saved to {STD_MATRIX_FILE}")

