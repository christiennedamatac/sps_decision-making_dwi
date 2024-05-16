#! /bin/sh
# Pre-process DWI data from the Healthy Brain Study with MRtrix protocol
# https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html
module load mrtrix
module load ANTs

cd project/3022060.01
data_dir=analysis/data
subjects=$(find "$data_dir" -maxdepth 1 -type d -exec basename {} \;)
#subjects="sub-0001_HBU002657EC4DFC5EEB9"

rm -r scripts/dwi_prep/temp_scripts
mkdir -p scripts/dwi_prep/temp_scripts 

for subject in $subjects; do
echo "Processing $subject"
source_dir=$data_dir/${subject}/dwi
out_dir=$source_dir/fba_dwi_prep
mkdir -p $out_dir

#1. Convert diffusion data
command_1="mrconvert $source_dir/eddy_output.nii.gz $out_dir/dwi.mif -fslgrad $source_dir/eddy_output.eddy_rotated_bvecs $source_dir/bval.bval -force"

#6. Remove bias field
command_6="dwibiascorrect ants $out_dir/dwi.mif $out_dir/dwi_unbiased.mif -bias $out_dir/dwi_bias_field.mif -mask $source_dir/b0_brain_mask.nii.gz -fslgrad $source_dir/eddy_output.eddy_rotated_bvecs $source_dir/bval.bval -scratch $source_dir/fba_dwi_prep -force"

#Submit to batch
#echo "Submitting to batch..."
#scriptname="/project/3022060.01/scripts/dwi_prep/temp_scripts/fba_dwi_prep_${subject}"
#echo "$command_1" > $scriptname
#echo "$command_6" >> $scriptname
#chmod a+rwx $scriptname
#sbatch --time=18:00:00 --mem=20gb $scriptname
done

for script in scripts/dwi_prep/temp_scripts/*; do
    bash "$script"
done

# NOTES
#fieldmap_folder=/project/3022043.01/DELTA/3T/converted/study/$subject/ses-mri01/fmap
# Not necessary because there will be no B0 pair creation 
#2. Perform dwidenoise
# If denoising and/or Gibbs ringing removal are performed as part of the preprocessing, they must be performed prior to any other processing steps: most other processing steps, in particular those that involve interpolation of the data, will invalidate the original properties of the image data that are exploited by dwidenoise and mrdegibbs at this stage, and would render the result prone to errors. Because the Healthy Brain Study already pre-processed the data, dwidenoise and mrdegibbs will not be included in this pre-processing protocol. 

#3. Perform mrdegibbs
# See comment under #2.

#4. Create B0 pair for distortion correction
# not necessary; Healthy Brain Study preprocessing already included FSL's TOPUP and EDDY

#5. Run dwipreproc
# not necessary; Healthy Brain Study preprocessing already included FSL's TOPUP and EDDY
