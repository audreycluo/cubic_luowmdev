#!/usr/bin/env python
from dipy.stats.analysis import afq_profile, gaussian_weights
from dipy.io.streamline import load_trk
from dipy.tracking.streamline import transform_streamlines, set_number_of_points, values_from_volume
from fury import actor, window
from fury.colormap import create_colormap
import glob 
import json
import numpy as np
import nibabel as nib
import os
from os.path import join as ospj
import pandas as pd
from scipy.ndimage import map_coordinates
import SimpleITK as sitk
import sys

'''
This script computes the distance between each node and the closest cortical endpoint with the following steps:

1. downsample streamlines to 100 points on a core bundle (mean of streamline coordinates)
2. For each of the 100 points, find the coordinate from each streamline and make a 3d volume where the voxel containing that point has a 1. All other voxels have 0
3. Sample a gray matter mask into this same voxel grid
4. Calculate the distance between the closest GM voxel and node 0 and node 99, respectively. (These endpoint nodes are right next to GM in most tracts, with the expection of ATR and CST)
5. For node i along my tract, get the distance between node 0 and i, and between node 99 and i by summing the distances between consecutive nodes 
6. Take the minimum between these two values as the distance between node i and cortical endpoint.
'''

########################################
# Set directories
########################################
# Parse command-line arguments
if len(sys.argv) != 3:
   print("Usage: python compute_dist_to_cortex.py <config_file> <subject_id>")
   sys.exit(1)

config_file = sys.argv[1]
#config_file="/cbica/projects/luo_wm_dev/code/tract_profiles/config/config_HCPD.json"
subject_id = sys.argv[2]
#subject_id = "sub-0001305" # hcpd

# Read config from the specified file
with open(config_file, "rb") as f:
    config = json.load(f)

dataset = config['dataset']
data_root = config['data_root']
qsiprep_dir = ospj(data_root, "datalad_qsiprep")
study_path = ospj(data_root, "babs_qsirecon_pyafq/merge_ds/qsirecon/")
deriv_path = ospj(study_path, subject_id)
 
if(dataset=="HCPD"):
    afq_path = ospj(deriv_path, "ses-V1/dwi", f"{subject_id}_ses-V1_space-T1w_desc-preproc/")
if(dataset=="HBN"):
    afq_path = glob.glob(ospj(deriv_path, "ses*/dwi", f"{subject_id}_*/"))[0] # need globs due to ses-{site_name}
 
bundle_path = ospj(afq_path,'clean_bundles')
 
outputs_root = config['outputs_root']
outputs_dir = ospj(outputs_root, subject_id)

if not os.path.exists(outputs_dir):
    os.makedirs(outputs_dir)
    print(f"{outputs_dir} created")

################ 
# Load images  #
################ 
t1w_img = nib.load(ospj(qsiprep_dir, f"qsiprep/{subject_id}/anat/{subject_id}_desc-preproc_T1w.nii.gz"))
t1w = t1w_img.get_fdata()

md_img = nib.load(glob.glob(ospj(afq_path, f'{subject_id}_*DTI_desc-MD_dwi.nii.gz'))[0])
md = md_img.get_fdata()

print("T1 and MD images loaded")

########################################
# Define functions
########################################
def get_streamlines(bundle_name):
    print(bundle_name)
    fname = glob.glob(ospj(bundle_path, f"{subject_id}_*{bundle_name}_tractography.trk"))[0]
    sft = load_trk(fname, md_img)
    print("trk's loaded")

    # We transform the bundle coordinates, 
    # first into the RASMM common coordinate frame 
    # and then subsequently into the coordinate frame of the T1-weighted data
    sft.to_rasmm()
    tract_t1w = transform_streamlines(sft.streamlines, np.linalg.inv(t1w_img.affine))
    print("streamlines transformed")
    return tract_t1w

def calculate_tract_profile_coords_and_distances(streamlines, num_points=100):
    core_bundle = np.median(np.asarray(set_number_of_points(streamlines, 100)), axis=0) # not sure if i should use median or mean for core bundle
    # Calculate Euclidean distances between consecutive nodes
    distances = np.linalg.norm(np.diff(core_bundle, axis=0), axis=1)
    return core_bundle, distances

def create_3d_volume(points, volume_shape, voxel_size):
    # 'points' is a np array of shape (100, 3) containing 100 points
    # 'volume_shape' is a tuple representing the shape of 3D volume, e.g., (x_dim, y_dim, z_dim)
    # 'voxel_size' is the size of each voxel in my 3D volume
    volume = np.zeros(volume_shape, dtype=np.uint8)

    # convert points coordinates to voxel indices
    voxel_indices = np.round(points / voxel_size).astype(int)

    # if voxel contains a point on the core bundle, set it to 1
    volume[voxel_indices[:, 0], voxel_indices[:, 1], voxel_indices[:, 2]] = 1

    return volume
 
# Resample GM mask to the same voxel grid
def resample_to_voxel_grid(data, target_shape, source_affine, target_affine):
    # Create a grid of coordinates in the target image space
    coords = np.indices(target_shape).reshape(3, -1)
    coords = np.vstack((coords, np.ones((1, coords.shape[1]))))

    # Transform coordinates to the source image space
    voxel_coords = np.linalg.inv(source_affine).dot(target_affine.dot(coords))[:3]

    # Use map_coordinates to resample the source data at the transformed coordinates
    resampled_data = map_coordinates(data, voxel_coords, order=1, mode='nearest')
    return resampled_data.reshape(target_shape)

def find_closest_gm_voxel(node_coords, binarized_gm_mask, voxel_size):
    # get shape of binarized GM mask
    mask_shape = binarized_gm_mask.shape

    # initialize variables to store closest voxel coordinates and minimum distance
    closest_voxel_coords = None
    min_distance = float('inf')

    # iterate over all voxels in the GM mask
    for x in range(mask_shape[0]):
        for y in range(mask_shape[1]):
            for z in range(mask_shape[2]):
                # calculate Euclidean distance from current voxel to node 0 coordinates
                distance = np.sqrt((x - node_coords[0])**2 + (y - node_coords[1])**2 + (z - node_coords[2])**2)

                # update closest voxel coordinates and minimum distance if necessary
                if distance < min_distance and binarized_gm_mask[x, y, z] == 1:
                    min_distance = distance
                    closest_voxel_coords = (x, y, z)

    return closest_voxel_coords

def dist_closest_gm_voxel(node_coords, closest_voxel_coords):
    distance = np.linalg.norm(closest_voxel_coords - node_coords)
    return distance

# compute distances to cortical endpoint for each node along the tract 
def compute_min_distance_to_cortical_endpoint(bundle_name, distances_to_gm_data):
    node_coords = core_bundles[bundle_name][0]  # coordinates of all nodes in the bundle
    node_distances = core_bundles[bundle_name][1]  # distances between consecutive nodes
    distances_to_gm_data = distances_to_gm[bundle_name]
    num_nodes = len(node_coords) # computing distances for nodes 0-99

    distances_to_cortical_endpoints = np.zeros(num_nodes)
    cumulative_distance_from_0 = np.zeros(num_nodes)
    cumulative_distance_from_99 = np.zeros(num_nodes)

    if "ATR" in bundle_name:
        # compute cumulative distances from node 0
        for i in range(1, num_nodes):
            cumulative_distance_from_0[i] = cumulative_distance_from_0[i-1] + node_distances[i-1] # contains cumulative distances for nodes 0-99
        cumulative_distance_from_0 += distances_to_gm_data['node_0'] # add the distance between node_0 and GM 
        distances_to_cortical_endpoints = cumulative_distance_from_0
    
    elif "CST" in bundle_name:
        # compute cumulative distances from node 99
        for i in range(num_nodes-2, -1, -1):
            cumulative_distance_from_99[i] = cumulative_distance_from_99[i+1] + node_distances[i]
        cumulative_distance_from_99 += distances_to_gm_data['node_99'] # add the distance between node_99 and GM 
        distances_to_cortical_endpoints = cumulative_distance_from_99
    
    else:
        # compute cumulative distances from node 0
        for i in range(1, num_nodes):
            cumulative_distance_from_0[i] = cumulative_distance_from_0[i-1] + node_distances[i-1]  
        cumulative_distance_from_0 += distances_to_gm_data['node_0']  
         
        # compute cumulative distances from node 99
        for i in range(num_nodes-2, -1, -1):
            cumulative_distance_from_99[i] = cumulative_distance_from_99[i+1] + node_distances[i]
        cumulative_distance_from_99 += distances_to_gm_data['node_99']  
         
        # Compute the minimum distance to cortical endpoints for each node
        for i in range(num_nodes):
            distances_to_cortical_endpoints[i] = min(cumulative_distance_from_0[i], cumulative_distance_from_99[i])

    return distances_to_cortical_endpoints

 
#############################################################
# 1. Resample streamlines to 100 points and get core bundle #
#############################################################
# for each bundle:
bundles = ["ARCL", "ARCR", "ATRL", "ATRR", "CGCL", "CGCR",  "CSTL", "CSTR", "FA", "FP", "IFOL", "IFOR", "ILFL", "ILFR","pARCL", "pARCR","SLFL", "SLFR","UNCL", "UNCR","VOFL", "VOFR"]
# bundles = ["ARCL", "ATRL", "ATRR", "CSTL", "UNCR"] 
core_bundles = {}
streamlines = {}
for bundle_name in bundles:
    streamlines[bundle_name] = get_streamlines(bundle_name)
    core_bundles[bundle_name] = calculate_tract_profile_coords_and_distances(streamlines[bundle_name])
print("Streamlines resampled to 100 points")

# save out consecutive distances 
bundle_names = list(core_bundles.keys())
data = []

for bundle_name in bundle_names:
    distances = core_bundles[bundle_name][1]
    # Create column names for consecutive distances
    columns = [f'node{i}_node{i+1}' for i in range(len(distances))]
    data.append({'bundle_name': bundle_name, **dict(zip(columns, distances))})
df_consec_dist = pd.DataFrame(data)
consec_csv_file_path = ospj(outputs_dir, 'consecutive_distances.csv')
df_consec_dist.to_csv(consec_csv_file_path, index=False)


################################################## 
# 2. Make a binary 3d volume for each tract profile 
##################################################
t1w_sitk_img = sitk.ReadImage(ospj(qsiprep_dir, f"qsiprep/{subject_id}/anat/{subject_id}_desc-preproc_T1w.nii.gz"))
volume_shape = t1w_img.shape  # Volume shape
voxel_size = t1w_sitk_img.GetSpacing()[0] # voxel size
 
# i actually don't think i need this code anymore
##tract_prof_volumes = {}
# create the 3D volume
#for bundle_name, points_list in core_bundles.items():
 #   points = points_list[0]
 #   volume = create_3d_volume(points, t1w.shape, voxel_size)
 #   tract_prof_volumes[bundle_name] = volume
 
 
##############################################################
### 3. Sample a gray matter mask into this same voxel grid ###
##############################################################
gm_mask_img = nib.load(ospj(data_root, f"freesurfer_qsiprep_xfm/{subject_id}/ribbon_transformed.nii.gz"))
gm_mask_data = gm_mask_img.get_fdata()
gm_mask_affine = gm_mask_img.affine
gm_mask_data = np.where(gm_mask_data != 0, 1, 0) # make sure that mask is binarized

# Check if the affines are the same
if np.array_equal(t1w_img.affine, gm_mask_img.affine):
    # No transformation needed, directly resample GM mask
    resampled_gm_mask = map_coordinates(gm_mask_data, np.indices(volume_shape), order=1, mode='nearest')
else:
    # Calculate the transformation matrix to resample GM mask to the T1w space
    gm_to_t1w_transform = np.linalg.inv(t1w_img.affine).dot(gm_mask_img.affine)
    resampled_gm_mask = resample_to_voxel_grid(gm_mask_data, volume_shape, gm_mask_img.affine, t1w_img.affine)
 
print("GM mask loaded")


################################################################################################################
### 4. Find closest GM voxel to nodes 0 and 99 respectively, then compute the respective Euclidean distances ###
################################################################################################################
'''
Initially, I wanted to use SimpleITK's Maurer distance, but I encountered a limitation: 
for tracts lacking cortical endpoints (such as the CST or ATR), the closest GM voxel 
simpleITK identified was NOT a cortical endpoint, but rather a neighboring GM voxel. 
For instance, for the ATR, the thalamic node was very close to the cingulum, resulting 
in the computation of the distance between the cingulum and the thalamic tract endpoint 
instead of the frontal cortex.

Instead, I have adopted this approach: 
For cortico-cortical tracts, I first identify the closest GM voxel to nodes 0 and 
99 respectively. Then, I compute the Euclidean distances between node 0 and its 
closest GM voxel, and node 99 and its closest GM voxel.

In tracts where there was only 1 cortical endpoint, I did the following:
I treated the non-cortical terminal node similarly to the other nodes along the tract.
i.e. If the non-cortical temrinal node is node 99:
I added the consecutive Euclidean distances between node 0 to node 99 (computing the along-tract length). To this
summed distance, I added the Euclidean distance between node 0 and its closest GM voxel. 

This is approach less computationally efficient, but the code is more transparent and 
most importantly, it actually compute the distance that we want!!
''' 


distances_to_gm = {}

# Compute distance between GM mask and node 0 and 99 respectively
for bundle_name, bundle_data in core_bundles.items():
    print(f"Computing distances for bundle: {bundle_name}")
    if "ATR" in bundle_name: # note that we're using the older version of pyAFQ here where tract orientations are NOT consistent!!
        # may need to update code if using updated pyAFQ!!!
        node_0_coords = bundle_data[0][0]  # Coordinates of node 0 - frontal endpoint
        node_99_coords = bundle_data[0][99]  # Coordinates of node 99 - thalamic endpoint

        # Find the closest GM voxel to node 0 - frontal endpoint
        closest_gm_voxel_0 = find_closest_gm_voxel(node_0_coords, gm_mask_data, voxel_size)
        distance_to_gm_voxel_0 = dist_closest_gm_voxel(node_0_coords, closest_gm_voxel_0)

        # compute distance along tract between 0 and 99
        distance_along_tract_to_99 = sum(bundle_data[1])
        
        # total distance between cortical endpoint and node 99 (thalamus)
        distance_to_gm_voxel_99 = distance_along_tract_to_99 + distance_to_gm_voxel_0
    elif "CST" in bundle_name:
        node_0_coords = bundle_data[0][0]  # Coordinates of node 0 - brainstem endpoint
        node_99_coords = bundle_data[0][99]  # Coordinates of node 99 - motor cortex endpoint

        # find the closest GM voxel to node 99 - motor cortex endpoint
        closest_gm_voxel_99 = find_closest_gm_voxel(node_99_coords, gm_mask_data, voxel_size)
        distance_to_gm_voxel_99 = dist_closest_gm_voxel(node_99_coords, closest_gm_voxel_99)
        
        # compute distance along tract between 0 and 99
        distance_along_tract_to_0 = sum(bundle_data[1])

        # compute distance between cortical endpoint and node 0 (brainstem)
        distance_to_gm_voxel_0 = distance_along_tract_to_0 + distance_to_gm_voxel_99
    else: # for all other tracts    
        node_0_coords = bundle_data[0][0]  # Coordinates of node 0
        node_99_coords = bundle_data[0][99]  # Coordinates of node 99

        # find the closest GM voxel to node 0 
        closest_gm_voxel_0 = find_closest_gm_voxel(node_0_coords, gm_mask_data, voxel_size)
        distance_to_gm_voxel_0 = dist_closest_gm_voxel(node_0_coords, closest_gm_voxel_0)

        # find the closest GM voxel to node 99
        closest_gm_voxel_99 = find_closest_gm_voxel(node_99_coords, gm_mask_data, voxel_size)
        distance_to_gm_voxel_99 = dist_closest_gm_voxel(node_99_coords, closest_gm_voxel_99)

    # save out to dictionary
    distances_to_gm[bundle_name] = {
        'node_0': distance_to_gm_voxel_0,
        'node_99': distance_to_gm_voxel_99
    }
     

# compute minimum distance to cortical endpoints 
distances_to_cortical_endpoints = {} 
 
for bundle_name in bundles: 
    distances_to_cortical_endpoints[bundle_name] = compute_min_distance_to_cortical_endpoint(bundle_name, distances_to_gm[bundle_name])
 

################################ 
### Save out distances as csv!!
################################ 
df = pd.DataFrame(distances_to_cortical_endpoints).transpose().reset_index()
df = df.rename(columns={'index': 'bundle_name'})
df.columns = ['bundle_name'] + [f'node_{i}_dist' for i in range(100)]

# save out
csv_file_path = ospj(outputs_dir, 'distances_to_gm.csv')
df.to_csv(csv_file_path, index=False)