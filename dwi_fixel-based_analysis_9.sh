#! /bin/sh
# FIXEL-BASED ANALYSIS
# https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html
module load mrtrix
cd Desktop/3022060.01 # Project directory

#9. Generate a study-specific unbiased FOD template from 40 subjects
# https://mrtrix.readthedocs.io/en/latest/reference/commands/population_template.html	 

# Randomly pick a subset of 40 subjects (call external script and capture output)
template_subjects=$(python scripts/make_template_subjects_list.py)
echo "TEMPLATE SUBJECTS:"
echo "$template_subjects"

# Create temporary template-building directories
OUT=Desktop/3022060.01/analysis/fixel-based_analysis # Output directory
mkdir -p $OUT/template/fod_input $OUT/template/mask_input

# Symbolic link images and masks 
for subject in $template_subjects; do
ln -sr analysis/data/${subject}/dwi/fba_dwi_prep/wmfod_norm.mif $OUT/template/fod_input/${subject}.mif
ln -sr analysis/data/${subject}/dwi/fba_dwi_prep/dwi_unbiased_upsampled_mask.mif $OUT/template/mask_input/${subject}.mif
done

# Build template
# If process times out of computing cluster, use the -continue option when executing the command the next time, using the last successfully created file (last shared file between two warp dirs), for example: "-continue $OUT/template/population_template-tmp-XGC1KW/ warps_15/389_HBUF6A4AD5EE00B79606.mif"
# Needed only one run to complete when using 4 cores per node: cput=135:14:50,walltime=37:18:50,mem=68719493120b       
echo "population_template $OUT/template/fod_input -mask_dir $OUT/template/mask_input $OUT/template/wmfod_template.mif -voxel_size 1.25 -scratch $OUT/template/ -nocleanup -force > $OUT/template/population_template.log 2>&1" | qsub -N 'wmfod_template' -V -l 'nodes=1:ppn=4,walltime=48:00:00,mem=64gb'

# Optional: print job status every 5 minutes
while true; do qstat; sleep 300; 

# Remove temporary directory after checking population_template output
#rm -r analysis/fixel-based_analysis/template/population_template-tmp-*/
