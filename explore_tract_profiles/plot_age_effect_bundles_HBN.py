from os.path import join as ospj
import nibabel as nib
import numpy as np
from dipy.io.streamline import load_trk
from dipy.tracking.streamline import transform_streamlines, set_number_of_points
from fury import actor, window
from fury.colormap import create_colormap
import pandas as pd
from matplotlib.cm import tab20
dataset = "HBN"
subject ='sub-NDARAC331VEH'
 

################ 
# Define paths #
################ 
study_path = ospj(f"/cbica/projects/luo_wm_dev/input/{dataset}/babs_qsirecon_pyafq/merge_ds/qsirecon/")
deriv_path = ospj(
    study_path, subject)
afq_path = ospj(
    deriv_path,
    'ses-HBNsiteCBIC', 'dwi',  f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc')
bundle_path = ospj(afq_path,
                      'clean_bundles')
out_folder = f"/Users/audluo/PennLINC/luowm_local/output/tract_profiles_testing/core_bundle_figs/{dataset}"


################ 
# Load images  #
################ 

# load FA and MD images
fa_img = nib.load(ospj(afq_path,
                          f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc_dwi_model-DTI_desc-FA_dwi.nii.gz'))
fa = fa_img.get_fdata()

md_img = nib.load(ospj(afq_path,
                          f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc_dwi_model-DTI_desc-MD_dwi.nii.gz'))
md = md_img.get_fdata()


# load T1w image
#t1w_img = nib.load(ospj(f'/cbica/projects/luo_wm_dev/input/HBN/HBN_qsiprep/{subject}/{subject}_desc-preproc_T1w.nii.gz'))
t1w_img = nib.load(ospj(f'/cbica/projects/luo_wm_dev/input/HBN/datalad_qsiprep/qsiprep/{subject}/anat/{subject}_desc-preproc_T1w.nii.gz'))

t1w = t1w_img.get_fdata()

# load age effects
gam_age_md = pd.read_csv("/cbica/projects/luo_wm_dev/output/HBN/tract_profiles/GAM/dti_md/GAMresults.dti_md.age_plotbundles.csv")

####################
# Define functions #
####################
def lines_as_tubes(sl, line_width, **kwargs): 
    line_actor = actor.line(sl, **kwargs)
    line_actor.GetProperty().SetRenderLinesAsTubes(1)
    line_actor.GetProperty().SetLineWidth(line_width)
    return line_actor


def slice_volume(data, x=None, y=None, z=None):
    # from https://yeatmanlab.github.io/pyAFQ/tutorials/tutorial_examples/plot_003_viz.html#sphx-glr-tutorials-tutorial-examples-plot-003-viz-py
    # The anatomical image is rendered using slicer actors. 
    # These are actors that visualize one slice of a three dimensional volume. 
    # Again, we create a helper function that will slice a volume along the x, y, and z dimensions. 
    # This function returns a list of the slicers we want to include in our visualization. 
    # This can be one, two, or three slicers, depending on how many of {x,y,z} are set. 
    slicer_actors = []
    slicer_actor_z = actor.slicer(data)
    if z is not None:
        slicer_actor_z.display_extent(
            0, data.shape[0] - 1,
            0, data.shape[1] - 1,
            z, z)
        slicer_actors.append(slicer_actor_z)
    if y is not None:
        slicer_actor_y = slicer_actor_z.copy()
        slicer_actor_y.display_extent(
            0, data.shape[0] - 1,
            y, y,
            0, data.shape[2] - 1)
        slicer_actors.append(slicer_actor_y)
    if x is not None:
        slicer_actor_x = slicer_actor_z.copy()
        slicer_actor_x.display_extent(
            x, x,
            0, data.shape[1] - 1,
            0, data.shape[2] - 1)
        slicer_actors.append(slicer_actor_x)

    return slicer_actors


def load_and_transform_tract(bundle_name, bundle_path, subject, fa_img, t1w_img):
    sft = load_trk(ospj(bundle_path, f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-{bundle_name}_tractography.trk'), fa_img)
    sft.to_rasmm()
    tract_t1w = transform_streamlines(sft.streamlines, np.linalg.inv(t1w_img.affine))
    return tract_t1w

def add_tract_to_scene(scene, tract_t1w, age_effects):
    """
    Add tracts to the scene with coloring based on age effects.

    Parameters:
    scene: fury.Scene
        The scene to add the tracts to.
    tract_t1w: list of arrays
        List of tractography streamlines.
    age_effects: list or array
        A vector of age effects for coloring the streamlines.
    """
    # Ensure age_effects is a numpy array
    age_effects = np.array(age_effects)
    
    # Normalize age effects to the range [0, 1] for colormap
    #normalized_age_effects = (age_effects - age_effects.min()) / (age_effects.max() - age_effects.min())
    
    # Create a colormap based on age effects
    #colors = create_colormap(age_effects, 'Spectral_r')
    
    for sl in tract_t1w:
        # Repeat the colormap for all points in the streamline
        #streamline_colors = np.tile(colors, (len(sl), 1))
        
        # Create the line actor with the streamline and its corresponding colors
        line_actor = actor.line([sl], opacity=0.3, colors = tab20.colors[1])
        
        # Add the line actor to the scene
        scene.add(line_actor)

    # add core bundles
    core = np.median(np.asarray(set_number_of_points(tract_t1w, 100)), axis=0)
     
    
    core_actor = lines_as_tubes(
        [core],
        40,
        colors=create_colormap(age_effects, 'inferno')
    )

    scene.add(core_actor)
   
         


#bundle_names = ['FA', 'FP', 'ARCL', 'ATRL', 'CGCL', 'CSTL', 'IFOL', 'ILFL', 'pARCL', 'SLFL', 'UNCL', 'VOFL', 
          #      'ARCR', 'ATRR', 'CGCR', 'CSTR', 'IFOR', 'ILFR', 'pARCR', 'SLFR', 'UNCR', 'VOFR']
bundle_names = ['CGCL', 'CSTL', 'SLFL','VOFL', 'ILFL']

# Load and transform tracts
tracts_t1w = {name: load_and_transform_tract(name, bundle_path, subject, fa_img, t1w_img) for name in bundle_names}

# Initialize scene
scene = window.Scene()
scene.set_camera(position=(779.66, 88.82, 62.28),
                 focal_point=(96.00, 114.00, 96.00),
                 view_up=(0, 0, 1))

# Add tracts to scene
for bundle_name, tract_t1w in tracts_t1w.items():
    age_effects = gam_age_md[gam_age_md['bundle_name'].str.contains(bundle_name)]['s_age.delta.adj.rsq_signed']
    add_tract_to_scene(scene, tract_t1w, age_effects)

# Add slicers to scene
slicers = slice_volume(t1w, x=t1w.shape[0] // 2, z=t1w.shape[-1] // 3)
for slicer in slicers:
    scene.add(slicer)

# Save the scene
output_path = ospj("/Users/audluo/PennLINC/luowm_local/output/tract_profiles_testing/both_datasets/", f'{subject}_all_bundles_age_effect_HBN.png')
window.record(scene, out_path=output_path, size=(1500, 1500))


 