# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Check dependencies
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError as e:
    print(f"Error: Missing dependency - {e}")
    exit(1)

# Read the contact map from the file
def read_contact_map(file_path):
    try:
        # Load the matrix from the tab-separated file
        contact_matrix = np.loadtxt(file_path, delimiter='\t')
        return contact_matrix
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

# Perform subtraction and save the result
def subtract_contact_maps(bound_file, unbound_file, output_file):
    # Read the matrices
    bound_matrix = read_contact_map(bound_file)
    unbound_matrix = read_contact_map(unbound_file)

    if bound_matrix is None or unbound_matrix is None:
        print("Failed to load one or both matrices.")
        return None

    # Ensure bound_matrix is at least as large as unbound_matrix
    if bound_matrix.shape[0] < unbound_matrix.shape[0] or bound_matrix.shape[1] < unbound_matrix.shape[1]:
        print("bound_matrix must be larger than or equal to unbound_matrix in size.")
        return None

    # Take the common size (unbound_matrix size) for subtraction
    common_size = min(bound_matrix.shape[0], unbound_matrix.shape[0]), min(bound_matrix.shape[1], unbound_matrix.shape[1])
    bound_subset = bound_matrix[:common_size[0], :common_size[1]]
    difference_matrix = bound_subset - unbound_matrix

    # Save the result to a new file with the same format
    np.savetxt(output_file, difference_matrix, delimiter='\t', fmt='%.18e')
    print(f"Difference matrix saved to {output_file}.")
    return difference_matrix

# Plot the contact map as a heatmap
def plot_contact_map(contact_matrix, file_name):
    plt.figure(figsize=(10, 8))
    # Create a custom colormap
    cmap = plt.cm.colors.ListedColormap(['green', 'white', 'red'])
    # Define bounds and norm for the colormap to cover -2.5, -0.5, 0.5, 2.5
    bounds = [-2.5, -0.5, 0.5, 2.5]
    norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)
    # Set background to white
    plt.gca().set_facecolor('white')
    # Flip matrix vertically to make y-axis start from bottom
    flipped_matrix = np.flipud(contact_matrix)
    # Plot heatmap with custom colormap
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
    
    # Add vertical and horizontal lines adjusted for flipped y-axis
    ax.axvline(x=46, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(x=79, color='black', linestyle='--', linewidth=1.5)
    ax.axhline(y=contact_matrix.shape[0] - 46, color='black', linestyle='--', linewidth=1.5)  # Adjusted for 46
    ax.axhline(y=contact_matrix.shape[0] - 79, color='black', linestyle='--', linewidth=1.5)  # Adjusted for 79
    
    # Add border around heatmap, adjust right border to contact_matrix.shape[1]
    width = contact_matrix.shape[1]
    height = contact_matrix.shape[0]
    ax.plot([0, width], [0, 0], color='black', linewidth=2)  # Bottom border
    ax.plot([0, width], [height, height], color='black', linewidth=2)  # Top border
    ax.plot([0, 0], [0, height], color='black', linewidth=2)  # Left border
    ax.plot([width, width], [0, height], color='black', linewidth=2)  # Right border
    
    # Set x and y axis ticks to match (0, 10, ..., 110)
    tick_interval = 20
    x_tick_positions = np.arange(0, contact_matrix.shape[1], tick_interval)
    y_tick_positions = np.arange(0, contact_matrix.shape[0], tick_interval)
    x_tick_labels = [str(i) for i in x_tick_positions]
    y_tick_labels = [str(i) for i in y_tick_positions]
    
    # Adjust y-ticks to start from bottom (0)
    y_display_positions = contact_matrix.shape[0] - 1 - y_tick_positions
    plt.xticks(x_tick_positions, x_tick_labels, fontsize=18)
    plt.yticks(y_display_positions, y_tick_labels, fontsize=18)
    
    # Save the figure with a modified name
    output_file = f"{os.path.basename(file_name).replace('.map', '_difference.tif')}"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', format='tiff')
    # Show plot only in interactive mode
    if os.environ.get('DISPLAY', None):
        plt.show()
    plt.close()  # Close the figure to free memory

# Main function
def main():
    # Get directory paths
    bound_dir = "bound"
    unbound_dir = "unbound"
    if not os.path.isdir(bound_dir) or not os.path.isdir(unbound_dir):
        print(f"Error: Directory {bound_dir} or {unbound_dir} does not exist.")
        exit(1)

    # Traverse all .map files in the directories
    bound_files = [f for f in os.listdir(bound_dir) if f.endswith('.map')]
    unbound_files = [f for f in os.listdir(unbound_dir) if f.endswith('.map')]
    if not bound_files or not unbound_files:
        print(f"Error: No .map files found in {bound_dir} or {unbound_dir}.")
        exit(1)

    # Match files by name
    common_files = set(bound_files) & set(unbound_files)
    if not common_files:
        print("Error: No matching .map files in bound and unbound directories.")
        exit(1)

    for file_name in common_files:
        bound_file = os.path.join(bound_dir, file_name)
        unbound_file = os.path.join(unbound_dir, file_name)
        output_file = os.path.join(bound_dir, f"contact_difference_{file_name.replace('.map', '')}.map")

        # Subtract matrices and plot
        difference_matrix = subtract_contact_maps(bound_file, unbound_file, output_file)
        if difference_matrix is not None:
            plot_contact_map(difference_matrix, output_file)

if __name__ == '__main__':
    main()
