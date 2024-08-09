from os.path import join as ospj
import nibabel as nib
import numpy as np
from dipy.io.streamline import load_trk
from dipy.tracking.streamline import transform_streamlines, set_number_of_points
from fury import actor, window
from fury.colormap import create_colormap
import pandas as pd

################################
# Set variables and directories#
################################ 
dataset = "HBN"
#subject = "sub-NDARZZ810LVF"  # age 8.5
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
t1w_img = nib.load(ospj(f'/cbica/projects/luo_wm_dev/input/{dataset}/datalad_qsiprep/qsiprep/{subject}/anat/{subject}_desc-preproc_T1w.nii.gz'))
t1w = t1w_img.get_fdata()

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


def view_bundles(L_t1w, R_t1w, sft_L, sft_R, L_profile, R_profile):

    slicers = slice_volume(t1w, x=t1w.shape[0] // 2, z=t1w.shape[-1] // 3)
    # from https://yeatmanlab.github.io/pyAFQ/tutorials/tutorial_examples/plot_003_viz.html#sphx-glr-tutorials-tutorial-examples-plot-003-viz-py
    # The next kind of fury object we will be working with is a window.Scene object. 
    # This is the (3D!) canvas on which we are drawing the actors. 
    # We initialize this object and call the scene.add method to add the actors.
    
    sft_L.to_vox()

 
    sl_L = list(L_t1w)
    for sl_L in L_t1w:
        L_actor = actor.line(
            [sl_L],
            colors=create_colormap(np.linspace(0, 100, len(sl_L)), 'Spectral_r'), opacity=0.3)
        scene.add(L_actor)

    sft_R.to_vox()
    sl_R = list(R_t1w)
    for sl_R in R_t1w:
        R_actor = actor.line(
            [sl_R],
            colors=create_colormap(np.linspace(0, 100, len(sl_R)), 'Spectral_r'), opacity=0.3)
        scene.add(R_actor)


    for slicer in slicers:
        scene.add(slicer)

    # add core bundles
    core_L = np.median(np.asarray(set_number_of_points(L_t1w, 100)), axis=0)
    core_R = np.median(np.asarray(set_number_of_points(R_t1w, 100)), axis=0)

    
    core_L_actor = lines_as_tubes(
        [core_L],
        40,
        colors=create_colormap(L_profile, 'inferno')
    )

    core_R_actor = lines_as_tubes(
        [core_R],
        40,
        colors=create_colormap(R_profile, 'inferno')
    )

    scene.add(core_L_actor)
    scene.add(core_R_actor)

    #window.show(scene)
    #exit()

def view_bundles_CC(t1w_img, sft, profile):


    slicers = slice_volume(t1w, x=t1w.shape[0] // 2, z=t1w.shape[-1] // 3)
  
    sl = list(t1w_img)
    for sl in t1w_img:
        CC_actor = actor.line(
            [sl],
            colors=create_colormap(np.linspace(0, 100, len(sl)), 'Spectral_r'), opacity=0.2)
        scene.add(CC_actor)

    
    for slicer in slicers:
        scene.add(slicer)

    # add core bundles
    core = np.median(np.asarray(set_number_of_points(t1w_img, 100)), axis=0)
 
    
    core_actor = lines_as_tubes(
        [core],
        30,
        colors=create_colormap(profile, 'inferno')
    )

    scene.add(core_actor)
    
 

def run_view_bundles(tract_name, scalar):
    if scalar == "dti_fa":
        scalar_img = fa_img
    else:
        scalar_img = md_img
    if tract_name == "FA" or tract_name == "FP":
        
        sft = load_trk(ospj(bundle_path,
                            f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-{tract_name}_tractography.trk'), scalar_img)
        print("trk's loaded")
        
        # We transform the bundle coordinates, 
        # first into the RASMM common coordinate frame 
        # and then subsequently into the coordinate frame of the T1-weighted data
        sft.to_rasmm()
        tract_t1w = transform_streamlines(sft.streamlines,
                                    np.linalg.inv(t1w_img.affine))
        
        print("streamlines transformed")

        # import already-computed tract profiles from pyAFQ babs project
        tract_profiles  = pd.read_csv(ospj(afq_path, f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-profiles_dwi.csv'))
        profile = tract_profiles[tract_profiles['tractID'] == tract_name][scalar]


        print("viewing bundles with window")
        view_bundles_CC(tract_t1w, sft, profile)
    
    else:
        sft_L = load_trk(ospj(bundle_path,
                            f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-{tract_name}L_tractography.trk'), scalar_img)
        sft_R = load_trk(ospj(bundle_path,
                            f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-{tract_name}R_tractography.trk'), scalar_img)
        print("trk's loaded")

        # We transform the bundle coordinates, 
        # first into the RASMM common coordinate frame 
        # and then subsequently into the coordinate frame of the T1-weighted data
        sft_L.to_rasmm()
        sft_R.to_rasmm()
        L_t1w = transform_streamlines(sft_L.streamlines,
                                    np.linalg.inv(t1w_img.affine))
        R_t1w = transform_streamlines(sft_R.streamlines,
                                    np.linalg.inv(t1w_img.affine))
        
        print("streamlines transformed")

        # import already-computed tract profiles from pyAFQ babs project
        tract_profiles  = pd.read_csv(ospj(afq_path, f'{subject}_ses-HBNsiteCBIC_acq-64dir_space-T1w_desc-preproc_dwi_space-RASMM_model-probCSD_algo-AFQ_desc-profiles_dwi.csv'))
        tract_name_L = tract_name + "_L"
        tract_name_R = tract_name + "_R"
        L_profile = tract_profiles[tract_profiles['tractID'] == tract_name_L][scalar]
        R_profile = tract_profiles[tract_profiles['tractID'] == tract_name_R][scalar]

        print("viewing bundles with window")
        view_bundles(L_t1w, R_t1w, sft_L, sft_R, L_profile, R_profile)




###########################
## View bundles and save ##
###########################

bundles = ["ARC", "ATR", "CGC",  "CST", "FA", "FP", "IFO", "ILF", "pARC", "SLF", "UNC", "VOF"]

 
# fractional anisotropy
for bundle in bundles:
    if bundle == "FA" or bundle == "FP":
        scene = window.Scene()
        scene.set_camera(position=(95.30, 75.69, 779.88),
                        focal_point=(96.00, 114.00, 96.00),
                        view_up=(0, -1, 0))
    
        run_view_bundles(bundle, "dti_fa")
        window.record(scene, out_path=ospj(out_folder, f'{subject}_{bundle}_dti_fa.png'), size=(1500, 1500))
        scene.clear()
    else: 
        # save left hemi
        scene = window.Scene()
    
        scene.set_camera(position=(779.66, 88.82, 62.28),
                        focal_point=(96.00, 114.00, 96.00),
                        view_up=(0, 0, 1))
        
        run_view_bundles(bundle, "dti_fa")
    
        window.record(scene, out_path=ospj(out_folder, f'{subject}_{bundle}_dti_fa_L.png'), size=(1500, 1500))
        scene.clear()

        # save right hemi
        scene = window.Scene()
        scene.set_camera(position=(-586.14, 64.17, 59.13),
                        focal_point=(96.00, 114.00, 96.00),
                        view_up=(0, 0, 1))
        
        run_view_bundles(bundle, "dti_fa")
        window.record(scene, out_path=ospj(out_folder, f'{subject}_{bundle}_dti_fa_R.png'), size=(1500, 1500))
        scene.clear()



# mean diffusivity
for bundle in bundles:
    if bundle == "FA" or bundle == "FP":
        scene = window.Scene()
        scene.set_camera(position=(95.30, 75.69, 779.88),
                        focal_point=(96.00, 114.00, 96.00),
                        view_up=(0, -1, 0))
    
        run_view_bundles(bundle, "dti_md")
        window.record(scene, out_path=ospj(out_folder, f'{subject}_{bundle}_dti_md.png'), size=(1500, 1500))
        scene.clear()
    else: 
        # save left hemi
        scene = window.Scene()
    
        scene.set_camera(position=(779.66, 88.82, 62.28),
                        focal_point=(96.00, 114.00, 96.00),
                        view_up=(0, 0, 1))
        
        run_view_bundles(bundle, "dti_md")
    
        window.record(scene, out_path=ospj(out_folder, f'{subject}_{bundle}_dti_md_L.png'), size=(1500, 1500))
        scene.clear()

        # save right hemi
        scene = window.Scene()
        scene.set_camera(position=(-586.14, 64.17, 59.13),
                        focal_point=(96.00, 114.00, 96.00),
                        view_up=(0, 0, 1))
        
        run_view_bundles(bundle, "dti_md")
        window.record(scene, out_path=ospj(out_folder, f'{subject}_{bundle}_dti_md_R.png'), size=(1500, 1500))
        scene.clear()

  