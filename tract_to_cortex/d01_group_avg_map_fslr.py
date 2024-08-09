import json
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from neuromaps.datasets import fetch_fslr
from neuromaps.images import construct_shape_gii, load_gifti, load_data
from neuromaps import images
from neuromaps.parcellate import Parcellater
import nibabel as nib
import numpy as np
import os
from os.path import join as ospj
import re
from surfplot import Plot
import sys
from tqdm import tqdm

"""
This script takes the binary, subject-level tract-to-cortex maps (in fsLR 32k) for each tract, and averages them across subjects. 
This produces a group-level tract-to-cortex map of the proportion of subjects that have a tract connecting to each cortical region.
Saves the maps out as giftis (no threshold applied)
You can also save out pngs (these visualize maps at different probability thresholds) - but code for figures would need to be run locally
on a mounted project -- haven't figured out rendering issues on cubic. 
"""

###########################
## Set variables and dirs #
###########################
dataset = sys.argv[1]
threshold = float(sys.argv[2])

config_file = f"/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_{dataset}.json"

with open(config_file, "rb") as f:
    config = json.load(f)

data_root = config['data_root']
derivs_dir = ospj(data_root, f"derivatives/{dataset}_vol_to_surf")
out_dir = ospj(derivs_dir, "group")

if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print(f"Directory {out_dir} created.")
else:
        print(f"Directory {out_dir} already exists.")

png_dir = ospj(out_dir, "figures")

if not os.path.exists(png_dir):
        os.makedirs(png_dir)
        print(f"Directory {png_dir} created.")
else:
        print(f"Directory {png_dir} already exists.")


###################
# Define functions 
###################
def load_tract_maps(subject_dir, filename):
    """Load the tract-to-cortex map for a given tract from a subject directory."""
    file_path = ospj(subject_dir, filename)
    if os.path.exists(file_path):
        tract_gifti = nib.load(file_path)
        tract_data = tract_gifti.darrays[0].data
        return tract_data
    else:
        raise FileNotFoundError(f'Tract file {file_path} not found.')

def average_tract_maps(derivs_dir):
    """Average tract-to-cortex maps across subjects to get population probability maps for each tract.
    Remember that each subject-level map is already binarized. 
    """
    subject_dirs = [os.path.join(derivs_dir, d, "fslr_32k") for d in os.listdir(derivs_dir) if d.startswith('sub-') and os.path.isdir(os.path.join(derivs_dir, d))]
    subject_ids = [os.path.basename(d) for d in os.listdir(derivs_dir) if d.startswith('sub-') and os.path.isdir(os.path.join(derivs_dir, d))]
    
    # Assume all subjects have the same tract names
    filenames = os.listdir(subject_dirs[0])
    tract_names = set()  # Use a set to avoid duplicates
    
    for filename in filenames:
        if filename.endswith('.shape.gii'):  # Ensure you're looking at the correct files
            parts = filename.split('.')
            if len(parts) > 2:
                tract_name_part = parts[1]
                tract_name = tract_name_part.split('_')[0]
                tract_names.add(tract_name)

    # Initialize dictionary to hold the sum of tract maps for each tract
    tract_sums = {tract: None for tract in tract_names}
    tract_counts = {tract: 0 for tract in tract_names}  # To keep track of how many subjects have each tract

    for subject_dir in tqdm(subject_dirs, desc="Processing subjects"):
        for filename in os.listdir(subject_dir):
            if filename.endswith('.shape.gii'):
                parts = filename.split('.')
                if len(parts) > 2:
                    tract_name_part = parts[1]
                    tract = tract_name_part.split('_')[0]
                    if tract in tract_names:
                        try:
                            tract_map = load_tract_maps(subject_dir, filename)
                            if tract_sums[tract] is None:
                                tract_sums[tract] = np.zeros_like(tract_map)
                            tract_sums[tract] += tract_map
                            tract_counts[tract] += 1
                        except FileNotFoundError:
                            continue

    # Average the tract maps
    print("Averaging across subjects")
    tract_averages = {tract: tract_sums[tract] / tract_counts[tract] if tract_counts[tract] > 0 else None for tract in tract_names}
    
    return tract_averages

def save_group_maps(tract_averages, output_dir):
    """Save the group-level tract-to-cortex maps."""
    os.makedirs(output_dir, exist_ok=True)
    
    for tract, avg_map in tract_averages.items():
        if avg_map is not None:
            filename = f"group_{tract}_1.25.shape.gii"
            output_file = ospj(output_dir, filename)
            data = avg_map.astype(np.float32)
            gii_data = nib.gifti.gifti.GiftiDataArray(data)
            gii_array = nib.gifti.gifti.GiftiImage(darrays=[gii_data])
            nib.save(gii_array, output_file)
            print(f"Saved GIFTI file for {tract}")


def save_tract_figs(tract_data, threshold, png_dir, hemi):
    """
   Threshold maps and save out
    Args:
    - tract_data (dict): Dictionary with tract names as keys and data arrays as values.
    - threshold (float): Threshold value to apply to the data arrays.
    - png_dir (str): Directory to save the output PNG files.
    - hemi (str): ('L' or 'R')  
    """
    for tract, data in tract_data.items():
        if tract.endswith("L"):
            thresholded_tract = np.where(data > threshold, data, 0)
            
            p = Plot(lh, views=['lateral', 'medial'], brightness = 0.7, zoom = 1.2)
            p.add_layer(thresholded_tract, cmap = aquamarine, cbar = True, color_range=(threshold,1)) # colors
            
            kws = dict(location='bottom', draw_border=False, aspect=10,
                    decimals=1, pad=0)
            fig = p.build(cbar_kws=kws)
            fig.axes[0].set_title(tract, pad= 20)
            fig.show()
            
            # Save the figure as PNG
            output_file = os.path.join(png_dir, f'{tract}_threshold{threshold}.png')
            fig.savefig(output_file, dpi=300)
            print(f'Saved {output_file}')
            
            # Close the figure to free memory
            plt.close(fig)
        else: 
            thresholded_tract = np.where(data > threshold, data, 0)
            
            p = Plot(rh, views=['lateral', 'medial'], brightness = 0.7, zoom = 1.2)
            p.add_layer(thresholded_tract, cmap = aquamarine, cbar = True, color_range=(threshold,1)) # colors
            
            kws = dict(location='bottom', draw_border=False, aspect=10,
                    decimals=1, pad=0)
            fig = p.build(cbar_kws=kws)
            fig.axes[0].set_title(tract, pad= 20)
            fig.show()
            
            # Save the figure as PNG
            output_file = os.path.join(png_dir, f'{tract}_threshold{threshold}.png')
            fig.savefig(output_file, dpi=300)
            print(f'Saved {output_file}')
            
            # Close the figure to free memory
            plt.close(fig)


#########################################################
# Compute population probability maps for each tract
#########################################################
tract_averages = average_tract_maps(derivs_dir)
save_group_maps(tract_averages, out_dir)
print("Group-level tract-to-cortex maps have been saved.")


#########################################################
# Save out images of these maps in fslr 32k 
#########################################################
# load my custom colormap :)
with open('/cbica/projects/luo_wm_dev/code/tract_profiles/colormaps/aquamarine.json', 'r') as f:
    hex_colors = json.load(f)
rgb_colors = [mcolors.hex2color(hex_color) for hex_color in hex_colors]
aquamarine = LinearSegmentedColormap.from_list("loaded_cmap", rgb_colors, N=256)


surfaces = fetch_fslr()
lh, rh = surfaces['veryinflated']

threshold = 0.3
# if want to save out figures, you'd need to run this code locally merp.
#save_tract_figs(tract_averages, threshold, png_dir, 'L')
#save_tract_figs(tract_averages, threshold, png_dir, 'R')