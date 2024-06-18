#! /bin/sh
# FIXEL-BASED ANALYSIS
# https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html
module load mrtrix
cd Desktop/3022060.01 # Project directory

# Create subjects list
subjects=$(find analysis/data -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort)
#subjects="sub-0001_HBU002657EC4DFC5EEB9"
echo "SUBJECTS LIST:"
echo "${subjects}"

# Create temporary directory to store individual processing scripts
rm -r scripts/dwi_prep/temp_scripts
mkdir -p scripts/dwi_prep/temp_scripts 

OUT=Desktop/3022060.01/analysis/fixel-based_analysis # Output directory

# Create individual processing scripts 
for subject in $subjects; do
echo "Processing $subject"

IN=Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep # Individual subject directory

#4. Computing (average) tissue response functions
step_4a="dwi2response dhollander $IN/dwi_unbiased.mif $IN/response_wm.txt $IN/response_gm.txt $IN/response_csf.txt"
step_4b="responsemean analysis/data/*/dwi/fba_dwi_prep/response_wm.txt $OUT/group_average_response_wm.txt"
step_4c="responsemean analysis/data/*/dwi/fba_dwi_prep/response_gm.txt $OUT/group_average_response_gm.txt"
step_4d="responsemean analysis/data/*/dwi/fba_dwi_prep/response_csf.txt $OUT/group_average_response_csf.txt"

#5. Upsampling DW images (~8min)
step_5="mrgrid $IN/dwi_unbiased.mif regrid -voxel 1.25 $IN/dwi_unbiased_upsampled.mif"

#6. Compute upsampled brain mask images (~3min)
step_6="dwi2mask $IN/dwi_unbiased_upsampled.mif $IN/dwi_unbiased_upsampled_mask.mif"
## Visually check masks; can contain non-brain

#7. Fibre Orientation Distribution estimation (multi-shell multi-tissue spherical deconvolution)
step_7="dwi2fod msmt_csd $IN/dwi_unbiased_upsampled.mif $OUT/group_average_response_wm.txt $IN/wmfod.mif $OUT/group_average_response_gm.txt $IN/gm.mif $OUT/group_average_response_csf.txt $IN/csf.mif -mask $IN/dwi_unbiased_upsampled_mask.mif -force"
## Visually check masks; brain masks must contain only brain for next step

#8. Joint bias field correction and intensity normalisation
step_8="mtnormalise $IN/wmfod.mif $IN/wmfod_norm.mif $IN/gm.mif $IN/gm_norm.mif $IN/csf.mif $IN/csf_norm.mif -mask $IN/dwi_unbiased_upsampled_mask.mif -force"

script_name="scripts/dwi_prep/temp_scripts/fba_4-8_${subject}"
echo "$step_4a" > $script_name
echo "$step_4b" >> $script_name
echo "$step_4c" >> $script_name
echo "$step_4d" >> $script_name
echo "$step_5" >> $script_name
echo "$step_6" >> $script_name
echo "$step_7" >> $script_name
echo "$step_8" >> $script_name
chmod a+rwx $script_name
# Run each script in computing cluster 
echo "Submitting to batch..."
qsub -N "fba_${subject}" -V -l 'nodes=1:ppn=4,walltime=24:00:00,mem=12gb' $script_name # Check resources request for batch processing
done
