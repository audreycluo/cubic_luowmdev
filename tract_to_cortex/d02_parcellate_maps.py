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
import pandas as pd
import re
from surfplot import Plot
import sys

"""
This script takes the population probability tract-to-cortex maps and parcellates them to Glasser. 
Saves out csv's of the maps (no threshold applied).
Can also save out pngs (these visualize maps at different probability thresholds), 
but code for figures would need to be run locally
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
dataset = config['dataset']
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
 
def parcellate_tract_maps(parcellator_object,  tract):
    """
    Apply a cortical parcellation to tract-to-region probability maps
    """
    left_gii = tract_averages[f"{tract}L"]
    right_gii = tract_averages[f"{tract}R"]
    combined_gii = [left_gii, right_gii]
    parcellated_maps = parcellator_object.fit_transform(combined_gii, 'fsLR') 
    return parcellated_maps

def save_parcellated_maps(parcellator_object,  tract, atlas):
    parcellated_maps = parcellate_tract_maps(parcellator_object,  tract)
    print(f"{tract} parcellated")

    lh_parcellated = parcellated_maps[:180]  # first 180 parcels for lh
    rh_parcellated = parcellated_maps[180:]  # next 180 parcels for rh

    # save parcellated lh
    lh_labels = np.unique(lh_glasser.darrays[0].data)[1:]  # exclude 0 (medial wall)
    lh_data = np.column_stack((lh_parcellated, lh_labels))
    lh_df = pd.DataFrame(lh_data, columns=['probability', 'regionID'])
    lh_df.to_csv(ospj(out_dir, f'{tract}L_{atlas}.csv'), index=False)
    
    # save parcellated rh
    rh_labels = np.unique(rh_glasser.darrays[0].data)[1:]  
    rh_data = np.column_stack((rh_parcellated, rh_labels))
    rh_df = pd.DataFrame(rh_data, columns=['probability', 'regionID'])
    rh_df.to_csv(ospj(out_dir, f'{tract}R_{atlas}.csv'), index=False)

    # create vertex-level data from the parcellated data
    lh_vertex_data = np.zeros(32492)
    rh_vertex_data = np.zeros(32492)

    # assign the parcellated values to the corresponding vertices for plotting
    for idx, label in enumerate(np.unique(lh_glasser.darrays[0].data)):
        if label != 0:  # don't map the medial wall  
            lh_vertex_data[lh_glasser.darrays[0].data == label] = lh_parcellated[idx - 1] # identify which vertices belong to the current parcel

    for idx, label in enumerate(np.unique(rh_glasser.darrays[0].data)):
        if label != 0:   
            rh_vertex_data[rh_glasser.darrays[0].data == label] = rh_parcellated[idx - 1]

    vertex_data = np.concatenate((lh_vertex_data, rh_vertex_data))
    print(f"{tract} parcellated data saved")
    return vertex_data

def save_parc_map_figs(parcellator_object, tract, threshold, atlas):
    vertex_data = save_parcellated_maps(parcellator_object, tract, atlas)
    # threshold map
    vertex_data = np.where(vertex_data > threshold, vertex_data, 0)
   
    p = Plot(lh, rh, views=['lateral', 'medial'], brightness=0.7, zoom=1.2)
    p.add_layer(vertex_data, cmap=aquamarine, cbar=True, color_range=(threshold, 1))
    p.add_layer(vertex_data, cmap = aquamarine, as_outline = True, cbar = False, alpha=0.6) # outline
                    
    kws = dict(location='bottom', draw_border=False, aspect=10,
            decimals=1, pad=0)
    fig = p.build(cbar_kws=kws)
    fig.axes[0].set_title(tract, pad= 20)
    output_file = os.path.join(png_dir, f'{tract}_threshold{threshold}_{atlas}.png')
    fig.savefig(output_file, dpi=300)
    print(f'Saved {output_file}')
    

    plt.close(fig)


#########################################################
# Load vertex-level population probability tract maps
#########################################################
tract_averages = {}
for filename in os.listdir(out_dir):
    if filename.endswith('.shape.gii'):
        # extract tract name
        parts = filename.split('_')
        if len(parts) >= 2:
            tract_name_part = parts[1]
            tract_name = tract_name_part.split('_')[0]
            file_path = os.path.join(out_dir, filename)
            tract_averages[tract_name] = file_path # save the gifti filepath into a dictionary {tractname:tractmap}


print(tract_averages.keys())

tracts = {tract[:-1] for tract in tract_averages.keys()}


#########################################################
# Load glasser
#########################################################
lh_glasser = nib.load("/cbica/projects/luo_wm_dev/atlases/glasser/FSLRVertex/glasser_360_L.label.gii")
rh_glasser = nib.load("/cbica/projects/luo_wm_dev/atlases/glasser/FSLRVertex/glasser_360_R.label.gii")

# generate neuromaps fsLR based Glasser 360 parcellation (needs relabeling to have consecutive region IDs)
glasser = images.relabel_gifti((lh_glasser, rh_glasser), background=['Medial_wall'])

# let's assign these labels so that we can call them later for plotting
lh_glasser = glasser[0]
rh_glasser = glasser[1]

# create parcellater object
glasser_parc = Parcellater(glasser, 'fsLR').fit()
print(glasser_parc)


#########################################################
# Parcellate, plot, and save!
#########################################################
# load my custom colormap
with open('/cbica/projects/luo_wm_dev/code/tract_profiles/colormaps/aquamarine.json', 'r') as f:
    hex_colors = json.load(f)
rgb_colors = [mcolors.hex2color(hex_color) for hex_color in hex_colors]
aquamarine = LinearSegmentedColormap.from_list("loaded_cmap", rgb_colors, N=256)

surfaces = fetch_fslr()
lh, rh = surfaces['veryinflated']

#for tract in tracts:
   #save_parc_map_figs(glasser_parc, tract, threshold, "glasser") # this would need to be run locally. merp.

for tract in tracts:
    save_parcellated_maps(glasser_parc, tract, "glasser") 