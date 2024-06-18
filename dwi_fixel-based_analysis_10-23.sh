#! /bin/sh
# FIXEL-BASED ANALYSIS
# https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html
module load mrtrix
cd Desktop/3022060.01 # Project directory

# Create subjects list
subjects=$(find analysis/data -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort)
#subjects="sub-0001_HBU002657EC4DFC5EEB9"
echo "SUBJECTS LIST: ${subjects}"

# Create temporary directory to store individual processing scripts
rm -r scripts/dwi_prep/temp_scripts
mkdir -p scripts/dwi_prep/temp_scripts 

OUT=Desktop/3022060.01/analysis/fixel-based_analysis # Output directory

# Create individual processing scripts 
for subject in $subjects; do
echo "Processing $subject"

IN=Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep # Individual subject directory

#10. Register all subject FOD images to the FOD template
step_10="mrregister $IN/wmfod_norm.mif -mask1 $IN/dwi_unbiased_upsampled_mask.mif $OUT/template/wmfod_template.mif -nl_warp $IN/subject2template_warp.mif $IN/template2subject_warp.mif -force" 

#11. Compute the template mask (intersection of all subject masks in template space)
step_11a="mrtransform $IN/dwi_unbiased_upsampled_mask.mif -warp $IN/subject2template_warp.mif -interp nearest -datatype bit $IN/dwi_mask_in_template_space.mif -force"

#Compute the template mask as the intersection of all warped masks
step_11b="mrmath analysis/data/*/dwi/fba_dwi_prep/dwi_mask_in_template_space.mif min $OUT/template/template_mask.mif -datatype bit -force"

#12. Compute a white matter template analysis fixel mask
step_12="fod2fixel -mask $OUT/template/template_mask.mif -fmls_peak_value 0.06 $OUT/template/wmfod_template.mif $OUT/template/fixel_mask -force"

#13. Warp FOD images to template space
step_13="mrtransform $IN/wmfod_norm.mif -warp $IN/subject2template_warp.mif -reorient_fod no $IN/fod_in_template_space_NOT_REORIENTED.mif -force" 

#14. Segment FOD images to estimate fixels and their apparent fibre density (FD)
step_14="fod2fixel -mask $OUT/template/template_mask.mif $IN/fod_in_template_space_NOT_REORIENTED.mif $IN/fixel_in_template_space_NOT_REORIENTED -afd fd.mif -force" 

#15. Reorient fixels
step_15="fixelreorient $IN/fixel_in_template_space_NOT_REORIENTED $IN/subject2template_warp.mif $IN/fixel_in_template_space -force"

#16. Assign subject fixels to template fixels
step_16="fixelcorrespondence $IN/fixel_in_template_space/fd.mif $OUT/template/fixel_mask $OUT/template/fd ${subject}.mif -force"

#17. Compute the fibre cross-section (FC) metric
step_17="warp2metric $IN/subject2template_warp.mif -fc $OUT/template/fixel_mask $OUT/template/fc ${subject}.mif -force"

script_name="scripts/dwi_prep/temp_scripts/fba_10-17_${subject}"
echo "$step_10" > $script_name
echo "$step_11a" >> $script_name
echo "$step_11b" >> $script_name
echo "$step_12" >> $script_name
echo "$step_13" >> $script_name
echo "$step_14" >> $script_name
echo "$step_15" >> $script_name
echo "$step_16" >> $script_name
echo "$step_17" >> $script_name
chmod a+rwx $script_name
# Run each script in computing cluster 
echo "Submitting to batch..."
qsub -N "fba_${subject}" -V -l 'nodes=1:ppn=4,walltime=24:00:00,mem=12gb' $script_name # Check resources request for batch processing
done
while true; do qstat; sleep 300; done

rm -r analysis/data/*/dwi/fba_dwi_prep/fixel_in_template_space_NOT_REORIENTED

##For group stats: create a separate fixel directory to store the log(FC) data, copy the fixel index and directions file, and calculate the log(FC) to ensure data are centred around zero and normally distributed
mkdir $OUT/template/log_fc
cp $OUT/template/fc/index.mif $OUT/template/fc/directions.mif $OUT/template/log_fc

for subject in $subjects; do
echo "Processing $subject"
mrcalc $OUT/template/fc/${subject}.mif -log $OUT/template/log_fc/${subject}.mif

#18. Compute a combined measure of fibre density and cross-section (FDC)
mkdir $OUT/template/fdc
cp $OUT/template/fc/index.mif $OUT/template/fdc
cp $OUT/template/fc/directions.mif $OUT/template/fdc

step_18="mrcalc $OUT/template/fd/${subject}.mif $OUT/template/fc/${subject}.mif -mult $OUT/template/fdc/${subject}.mif -force"

#19. Perform whole-brain fibre tractography on the FOD template
step_19="tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 $OUT/template/wmfod_template.mif -seed_image $OUT/template/template_mask.mif -mask $OUT/template/template_mask.mif -select 20000000 -cutoff 0.06 $OUT/tracks_20_million.tck -force" 

#20. Reduce biases in tractogram densities
step_20="tcksift $OUT/tracks_20_million.tck $OUT/template/wmfod_template.mif $OUT/tracks_2_million_sift.tck -term_number 2000000 -out_selection tracks_2_million.txt -force"

#21. Generate fixel-fixel connectivity matrix
step_21="fixelconnectivity $OUT/template/fixel_mask/ $OUT/tracks_2_million_sift.tck $OUT/matrix/" 

#22. Smooth fixel data using fixel-fixel connectivity
step_22a="fixelfilter $OUT/template/fd smooth $OUT/template/fd_smooth -matrix $OUT/matrix/ -force"
step_22b="fixelfilter $OUT/template/log_fc smooth $OUT/template/log_fc_smooth -matrix $OUT/matrix/ -force" 
step_22c="fixelfilter $OUT/template/fdc smooth $OUT/template/fdc_smooth -matrix $OUT/matrix/ -force"

script_name="scripts/dwi_prep/temp_scripts/fba_18-22_${subject}"
echo "$step_18" > $script_name
echo "$step_19" >> $script_name
echo "$step_20" >> $script_name
echo "$step_21" >> $script_name
echo "$step_22a" >> $script_name
echo "$step_22b" >> $script_name
echo "$step_22c" >> $script_name
chmod a+rwx $script_name
# Run each script in computing cluster 
echo "Submitting to batch..."
qsub -N "fba_${subject}" -V -l 'nodes=1:ppn=4,walltime=24:00:00,mem=12gb' $script_name # Check resources request for batch processing
done
while true; do qstat; sleep 300; done

#23. Perform statistical analysis of FD, FC, and FDC
ANALYSIS_DIR=Desktop/3022060.01/analysis

# Create files.txt
fixel_dir="$ANALYSIS_DIR/fixel-based_analysis/template/fdc_smooth"
output_file="$ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/files.txt"

# Make sure the output directory exists
mkdir -p "$(dirname "$output_file")"

# List filenames and write to output file
find "$fixel_dir" -type f -name 'sub-*' -exec basename {} \; > "$output_file"

mkdir -p $ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/fixelcfestats_out/
mkdir -p $ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/scripts/

# Create scripts 
for network in DMN FPN VAN; do
for fba_metric in fdc_smooth fd_smooth log_fc_smooth; do
for design_matrix in design_positive_z design_negative_z; do # z-scored SPSQ dimensions
step_23="fixelcfestats -mask $ANALYSIS_DIR/tractseg/tck2fixel_output/network_mask_${network}_final.mif $ANALYSIS_DIR/fixel-based_analysis/template/${fba_metric} $ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/files.txt $ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/${design_matrix}.csv $ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/contrast.txt $ANALYSIS_DIR/fixel-based_analysis/matrix/ $ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/fixelcfestats_out/${network}_${fba_metric}_${design_matrix}/ -force"
script_name="$ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/scripts/${network}_${fba_metric}_${design_matrix}"
echo "$step_23" > $script_name
chmod a+rwx $script_name
done
done
done

# Run each script in computing cluster 
echo "Submitting to batch..."
for script in $ANALYSIS_DIR/fixel-based_analysis/fixelcfestats/scripts/*; do
  if [[ -f "$script" ]]; then
    qsub -N "$(basename "$script")" -V -l 'nodes=1:ppn=4,walltime=48:00:00,mem=256gb' "$script"
    echo "Submitted: $(basename "$script")"
  else
    echo "Script not found: $script"
  fi
done
while true; do qstat; sleep 300; done
