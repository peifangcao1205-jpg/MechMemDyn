# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Read contact map file
def read_contact_map(file_path):
    try:
        # Load matrix from tab-separated file
        contact_matrix = np.loadtxt(file_path, delimiter='\t')
        print(f"Successfully loaded matrix from {file_path}, shape {contact_matrix.shape}")
        return contact_matrix
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None

# Plot contact map heatmap
def plot_contact_map(contact_matrix, file_name):
    plt.figure(figsize=(10, 8))
    # Create custom color map
    cmap = plt.cm.colors.ListedColormap(['green', 'white', 'red'])
    # Define color mapping boundaries and normalization
    bounds = [-1.5, -0.5, 0.5, 1.5]
    norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)
    # Set background to gray
    plt.gca().set_facecolor('gray')
    # Flip matrix vertically so y-axis starts from bottom
    flipped_matrix = np.flipud(contact_matrix)
    # Plot heatmap using custom color map
    ax = sns.heatmap(flipped_matrix, cmap=cmap, norm=norm, cbar=False, square=True)
    # Find and plot enlarged non-zero points
    indices = np.where(flipped_matrix != 0)
    x_coords, y_coords = indices[1], indices[0]
    x_coords = np.clip(x_coords, 0, contact_matrix.shape[1] - 1)
    y_coords = np.clip(y_coords, 0, contact_matrix.shape[0] - 1) + 1
    colors = ['red' if val > 0 else 'green' for val in flipped_matrix[indices]]
    ax.scatter(x_coords, y_coords, s=80, c=colors, edgecolors='none', marker='s')  # Enlarged square points

    #plt.title(f'Contact Map bound ({os.path.basename(file_name).replace(".map", "")})')
    plt.xlabel('Residue Index', fontsize=18, fontweight='bold')
    plt.ylabel('Residue Index', fontsize=18, fontweight='bold')
    
    # Add vertical and horizontal divider lines (adjust y-axis to match flip)
    ax.axvline(x=46, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(x=79, color='black', linestyle='--', linewidth=1.5)
    ax.axhline(y=contact_matrix.shape[0] - 46, color='black', linestyle='--', linewidth=1.5)  # Adjusted to 46
    ax.axhline(y=contact_matrix.shape[0] - 79, color='black', linestyle='--', linewidth=1.5)  # Adjusted to 79
    
    # Add border around heatmap, adjust right border to contact_matrix.shape[1]
    width = contact_matrix.shape[1]
    height = contact_matrix.shape[0]
    ax.plot([0, width], [0, 0], color='black', linewidth=2)  # Bottom border
    ax.plot([0, width], [height, height], color='black', linewidth=2)  # Top border
    ax.plot([0, 0], [0, height], color='black', linewidth=2)  # Left border
    ax.plot([width, width], [0, height], color='black', linewidth=2)  # Right border
    
    # Set x and y axis ticks (0, 10, ..., 110)
    tick_interval = 20
    x_tick_positions = np.arange(0, contact_matrix.shape[1], tick_interval)
    y_tick_positions = np.arange(0, contact_matrix.shape[0], tick_interval)
    x_tick_labels = [str(i) for i in x_tick_positions]
    y_tick_labels = [str(i) for i in y_tick_positions]
    
    # Adjust y-axis ticks to start from bottom (0)
    y_display_positions = contact_matrix.shape[0] - 1 - y_tick_positions
    plt.xticks(x_tick_positions, x_tick_labels, fontsize=18)
    plt.yticks(y_display_positions, y_tick_labels, fontsize=18)
    
    # Save image, modify file name
    output_file = f"{os.path.basename(file_name).replace('.map', '_bound.tif')}"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', format='tiff')
    # Display figure only in interactive mode
    if os.environ.get('DISPLAY', None):
        plt.show()
    plt.close()  # Close figure to free memory

# Main function
def main():
    # Get directory path
    dir_path = "bound"
    if not os.path.isdir(dir_path):
        print(f"Error: Directory {dir_path} does not exist.")
        exit(1)

    # Iterate over .map files in directory, exclude difference files
    map_files = [f for f in os.listdir(dir_path) if f.endswith('.map') and not f.startswith('contact_difference_')]
    if not map_files:
        print(f"Error: No original .map files found in directory {dir_path}.")
        exit(1)

    for file_name in map_files:
        file_path = os.path.join(dir_path, file_name)
        contact_matrix = read_contact_map(file_path)
        if contact_matrix is not None:
            plot_contact_map(contact_matrix, file_path)

if __name__ == '__main__':
    main()
