#! /bin/sh
# FIXEL-BASED ANALYSIS
# https://mrtrix.readthedocs.io/en/latest/fixel_based_analysis/mt_fibre_density_cross-section.html
module load mrtrix
cd Desktop/3022060.01

# Create subject list
subjects=$(find analysis/data -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort)
echo ${subjects}
#subjects="sub-0001_HBU002657EC4DFC5EEB9"

OUT=Desktop/3022060.01/analysis/fixel-based_analysis

for subject in $subjects; do
echo "Processing $subject"
IN=Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep

#4. Computing (average) tissue response functions
dwi2response dhollander $IN/dwi_unbiased.mif $IN/response_wm.txt $IN/response_gm.txt $IN/response_csf.txt 
responsemean analysis/data/*/dwi/fba_dwi_prep/response_wm.txt $OUT/group_average_response_wm.txt
responsemean analysis/data/*/dwi/fba_dwi_prep/response_gm.txt $OUT/group_average_response_gm.txt
responsemean analysis/data/*/dwi/fba_dwi_prep/response_csf.txt $OUT/group_average_response_csf.txt

#5. Upsampling DW images (~8min)
mrgrid $IN/dwi_unbiased.mif regrid -voxel 1.25 $IN/dwi_unbiased_upsampled.mif

#6. Compute upsampled brain mask images (~3min)
dwi2mask $IN/dwi_unbiased_upsampled.mif $IN/dwi_unbiased_upsampled_mask.mif

## Visually check masks; can contain non-brain

#7. Fibre Orientation Distribution estimation (multi-shell multi-tissue spherical deconvolution)
for subject in $subjects; do
echo "dwi2fod msmt_csd Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/dwi_unbiased_upsampled.mif Desktop/3022060.01/analysis/fixel-based_analysis/group_average_response_wm.txt Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/wmfod.mif Desktop/3022060.01/analysis/fixel-based_analysis/group_average_response_gm.txt Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/gm.mif Desktop/3022060.01/analysis/fixel-based_analysis/group_average_response_csf.txt Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/csf.mif -mask Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/dwi_unbiased_upsampled_mask.mif -force" | qsub -V -l 'nodes=1:ppn=4,mem=8gb,walltime=01:30:00' -N ${subject}
done

## Visually check masks; brain masks must contain only brain for next step

#8. Joint bias field correction and intensity normalisation
for subject in $missing_subjects_list; do
echo "mtnormalise Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/wmfod.mif Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/wmfod_norm.mif Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/gm.mif Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/gm_norm.mif Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/csf.mif Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/csf_norm.mif -mask Desktop/3022060.01/analysis/data/${subject}/dwi/fba_dwi_prep/dwi_unbiased_upsampled_mask.mif -force" | qsub -V -l 'mem=8gb,walltime=00:05:00' -N ${subject}
done

#9. Generate a study-specific unbiased FOD template from 40 subjects
https://mrtrix.readthedocs.io/en/latest/reference/commands/population_template.html	 

# Randomly pick a subset of subjects (call external Python script and capture output)
template_subjects=$(python scripts/make_template_subjects_list.py)
echo "Template subjects:"
echo "$template_subjects"

# Create temporary template-building directories
mkdir -p analysis/fixel-based_analysis/template/fod_input analysis/fixel-based_analysis/template/mask_input

# Symbolic link images and masks 
for subject in $template_subjects; do
ln -sr analysis/data/${subject}/dwi/fba_dwi_prep/wmfod_norm.mif analysis/fixel-based_analysis/template/fod_input/${subject}.mif
ln -sr analysis/data/${subject}/dwi/fba_dwi_prep/dwi_unbiased_upsampled_mask.mif analysis/fixel-based_analysis/template/mask_input/${subject}.mif
done

# Build template
# This process will most likely time out of the computing cluster, so need to use -continue option the next time, using the last successfully created file (last shared file between two warp dirs)
#-continue /project/3022028.01/fba_delta/population_template-tmp-8O4C0L/ warps_14/015.mif 
echo "population_template Desktop/3022060.01/analysis/fixel-based_analysis/template/fod_input -mask_dir Desktop/3022060.01/analysis/fixel-based_analysis/template/mask_input Desktop/3022060.01/analysis/fixel-based_analysis/template/wmfod_template.mif -voxel_size 1.25 -scratch Desktop/3022060.01/analysis/fixel-based_analysis/template/ -nocleanup -force > Desktop/3022060.01/analysis/fixel-based_analysis/template/population_template.log 2>&1" | qsub -N 'wmfod_template' -V -l 'nodes=1:ppn=4,walltime=72:00:00,mem=64gb'
while true; do qstat; sleep 300; done

#10. Register all subject FOD images to the FOD template
mkdir $out_dir/${subject}
mrregister $source_dir/${subject}/${subject}_wmfod_norm.mif -mask1 $source_dir/${subject}/${subject}_ses-mri01_acq-diffusion_run-1_dwi_biascorr_mask_upsampled.mif $out_dir/wmfod_template.mif -nl_warp $out_dir/${subject}/${subject}_subject2template_warp.mif $out_dir/${subject}/${subject}_template2subject_warp.mif -force

#11. Compute the template mask (intersection of all subject masks in template space)
mrtransform $source_dir/${subject}/${subject}_ses-mri01_acq-diffusion_run-1_dwi_biascorr_mask_upsampled.mif -warp $out_dir/${subject}/${subject}_subject2template_warp.mif -interp nearest -datatype bit $out_dir/${subject}/${subject}_dwi_mask_in_template_space.mif -force
#Compute the template mask as the intersection of all warped masks
mrmath $out_dir/*/*_dwi_mask_in_template_space.mif min $out_dir/template_mask.mif -datatype bit -force

#12. Compute a white matter template analysis fixel mask
fod2fixel -mask $out_dir/template_mask.mif -fmls_peak_value 0.06 $out_dir/wmfod_template.mif $out_dir/fixel_mask -force

#13. Warp FOD images to template space
mrtransform $source_dir/${subject}/${subject}_wmfod_norm.mif -warp $out_dir/${subject}/${subject}_subject2template_warp.mif -reorient_fod no $out_dir/${subject}/${subject}_fod_in_template_space_NOT_REORIENTED.mif -force

#14. Segment FOD images to estimate fixels and their apparent fibre density (FD)
fod2fixel -mask $out_dir/template_mask.mif $out_dir/${subject}/${subject}_fod_in_template_space_NOT_REORIENTED.mif $out_dir/${subject}/${subject}_fixel_in_template_space_NOT_REORIENTED -afd fd.mif -force

#15. Reorient fixels
fixelreorient $out_dir/${subject}/${subject}_fixel_in_template_space_NOT_REORIENTED $out_dir/${subject}/${subject}_subject2template_warp.mif $out_dir/${subject}/${subject}_fixel_in_template_space -force

#16. Assign subject fixels to template fixels
fixelcorrespondence $out_dir/${subject}/${subject}_fixel_in_template_space/fd.mif $out_dir/fixel_mask $out_dir/fd ${subject}.mif -force

#17. Compute the fibre cross-section (FC) metric
warp2metric $out_dir/${subject}/${subject}_subject2template_warp.mif -fc $out_dir/fixel_mask $out_dir/fc ${subject}.mif -force
##For group stats: create a separate fixel directory to store the log(FC) data, copy the fixel index and directions file, and calculate the log(FC) to ensure data are centred around zero and normally distributed
mkdir $out_dir/log_fc
cp $out_dir/fc/index.mif $out_dir/fc/directions.mif $out_dir/log_fc
mrcalc $out_dir/fc/${subject}.mif -log $out_dir/log_fc/${subject}.mif

#18. Compute a combined measure of fibre density and cross-section (FDC)
mkdir $out_dir/fdc
cp $out_dir/fc/index.mif $out_dir/fdc
cp $out_dir/fc/directions.mif $out_dir/fdc
mrcalc $out_dir/fd/${subject}.mif $out_dir/fc/${subject}.mif -mult $out_dir/fdc/${subject}.mif -force

done

#19. Perform whole-brain fibre tractography on the FOD template
tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 $out_dir/wmfod_template.mif -seed_image $out_dir/template_mask.mif -mask $out_dir/template_mask.mif -select 20000000 -cutoff 0.06 tracks_20_million.tck -force

#20. Reduce biases in tractogram densities
# approximately 1GB per million streamlines
echo "tcksift $out_dir/tracks_20_million.tck $out_dir/wmfod_template.mif $out_dir/tracks_2_million_sift.tck -term_number 2000000 -out_selection tracks_2_million.txt -force" 
#Used resources:	   cput=09:02:01,walltime=09:04:34,mem=49037258752b

#21. Generate fixel-fixel connectivity matrix
echo "fixelconnectivity $out_dir/fixel_mask/ $out_dir/tracks_2_million_sift.tck $out_dir/matrix/" 
#Used resources:	   cput=00:30:51,walltime=00:31:23,mem=47082053632b

#22. Smooth fixel data using fixel-fixel connectivity
echo "fixelfilter $out_dir/fd smooth $out_dir/fd_smooth -matrix $out_dir/matrix/ -force" 
#Used resources:	   cput=00:49:39,walltime=00:49:55,mem=9068638208b
echo "fixelfilter $out_dir/log_fc smooth $out_dir/log_fc_smooth -matrix $out_dir/matrix/ -force" 
#Used resources:	   cput=00:45:33,walltime=00:45:47,mem=423854080b
echo "fixelfilter $out_dir/fdc smooth $out_dir/fdc_smooth -matrix $out_dir/matrix/ -force" 
#Used resources:	   cput=00:41:39,walltime=00:41:53,mem=9054224384b

tckmap tracks_2_million_sift.tck outputTDI.mif -vox 1.25 -dec -force

#23. Perform statistical analysis of FD, FC, and FDC
#need permutations set generated with FSL PALM
#need tract regions-of-interest generated with TractSeg

#create scripts 
DATE=`date +"%Y-%m-%d"`
for F in fdc_smooth fd_smooth log_fc_smooth;do
for D in design_total_demean design_HI_demean;do
for T in CST_left SLF_I_left SLF_II_left SLF_III_left;do
mkdir -p $OUT_DIR/fixelcfestats_out/${DATE}
command_CFE="fixelcfestats -mask /project/3022028.01/fixel_based_analysis/tractseg_w2/tck2fixel_out/${T}.mif -permutations $OUT_DIR/palm/permutations_set_5000.csv $TEMPLATE_DIR/${F} $TEMPLATE_DIR/files.txt $OUT_DIR/${D}.csv $OUT_DIR/contrast.txt $TEMPLATE_DIR/matrix $OUT_DIR/fixelcfestats_out/${DATE}/${F}_${D}_${T}/ -force" 
SCRIPT="$OUT_DIR/fixelcfestats_out/21jul2021/script_${F}_${D}_${T}"
echo "$command_CFE" > $SCRIPT
chmod a+rwx $SCRIPT
done
done
done

#run each script in cluster
for_each * : qsub -N fixelcfestats -l nodes=1:ppn=4,walltime=48:00:00,mem=64gb ${script}