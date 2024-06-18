#! /bin/sh
# Use TractSeg to segment region-of-interest white matter fiber tracts 
# https://github.com/MIC-DKFZ/TractSeg

# Set up environment
module load mrtrix
module load anaconda3/2023.03
conda create --name fba_20240606 -y
source activate fba_20240606
#conda install pytorch torchvision torchaudio cudatoolkit=10.2 -c pytorch -y
#python -c "import torch; print(torch.__version__)"
#pip install TractSeg

# Define user directories
TEMPLATE_DIR=Desktop/3022060.01/analysis/fixel-based_analysis/template
TRACTSEG_DIR=Desktop/3022060.01/analysis/tractseg

# Create directories if they don't exist
mkdir -p $TEMPLATE_DIR
mkdir -p $TRACTSEG_DIR

# Use existing peaks
# Extract peak image from white matter FOD population template: avoids generating the MRtrix CSD peaks every time you run TractSeg by providing them directly by skipping option --raw_diffusion_input; input image must be a peak image (nifti 4D image with dimensions [x,y,z,9])
sh2peaks analysis/fixel-based_analysis/template/wmfod_template.mif -threshold 0.1 analysis/fixel-based_analysis/template/peaks.nii.gz -force

# Segment tracts
# https://github.com/MIC-DKFZ/TractSeg/blob/master/resources/Tutorial.md
# Blobs made up of only a few voxels are removed; deactivate this option by using --no_postprocess
echo "TractSeg -i Desktop/3022060.01/analysis/fixel-based_analysis/template/peaks.nii.gz -o Desktop/3022060.01/analysis/tractseg/tractseg_output" | qsub -N seg_tracts -V -l 'nodes=1:ppn=4,walltime=00:20:00,mem=15gb'
while true; do qstat; sleep 300; done

# Segment bundle start and end regions
echo "TractSeg -i Desktop/3022060.01/analysis/fixel-based_analysis/template/peaks.nii.gz -o Desktop/3022060.01/analysis/tractseg/tractseg_output --output_type endings_segmentation" | qsub -N seg_ends -V -l walltime=00:20:00,mem=20gb

# Create tract orientation maps
echo "TractSeg -i Desktop/3022060.01/analysis/fixel-based_analysis/template/peaks.nii.gz -o Desktop/3022060.01/analysis/tractseg/tractseg_output --output_type TOM" | qsub -N create_TOMs -V -l walltime=00:10:00,mem=15gb

# Track subset of bundles: probabilistic tractography on specific tracts
# Initialize an empty string to hold tract labels
tract_labels=""
# Read tracts index file and process each line
while IFS=, read -r index tract_label tract_name hemisphere network network_label notes; do
    # Skip the header row
    if [[ "$index" == "index" ]]; then
        continue
    fi
    # Check if the network_label matches the desired values
    if [[ "$network_label" == "FPN" || "$network_label" == "DMN" || "$network_label" == "VAN" ]]; then
        # Trim spaces and append the tract_label to the list, separated by commas
        trimmed_tract_label=$(echo "$tract_label" | xargs)  # Trim spaces
        if [[ -z "$tract_labels" ]]; then
            tract_labels="$trimmed_tract_label"
        else
            tract_labels="$tract_labels,$trimmed_tract_label"
        fi
    fi
done < Desktop/3022060.01/analysis/tractseg/index_tracts.csv
# Print list of tract labels
echo "Tract labels: $tract_labels"
# Run  Tracking command for specific tract labels
Tracking -i Desktop/3022060.01/analysis/fixel-based_analysis/template/peaks.nii.gz -o Desktop/3022060.01/analysis/tractseg/tractseg_output --tracking_format tck --bundles $tract_labels

# Convert to a fixel map to segment fixels corresponding to tract
# https://community.mrtrix.org/t/segmenting-white-matter-tracts-for-fba-metric-computation/3746/5
# https://community.mrtrix.org/t/mapping-segmented-tracts-back-to-fixels/4202
for tract in analysis/tractseg/tractseg_output/TOM_trackings/*.tck; do
    # Extract filename without extension
    tract_name=$(basename "$tract" .tck)
    
    # Run tck2fixel command
    tck2fixel "$tract" analysis/fixel-based_analysis/template/fdc_smooth analysis/tractseg/tck2fixel_output "$tract_name".mif 
done


# Below portion needs to be more concise

mv analysis/tractseg/tck2fixel_output/* analysis/tractseg/tck2fixel_output/individual_tracts/

# Create network lists of tracts
DMN_list=$(awk -F ',' '$6 == "DMN" {print $2}' analysis/tractseg/index_tracts.csv)
echo $DMN_list #CC_6 CG_left CG_right T_PREF_left T_PREF_right
FPN_list=$(awk -F ',' '$6 == "FPN" {print $2}' analysis/tractseg/index_tracts.csv)
echo $FPN_list # AF_left AF_right CC_1 CC_2 SLF_I_left SLF_I_right SLF_II_left SLF_II_right SLF_III_left SLF_III_right ST_PREF_left ST_PREF_right
VAN_list=$(awk -F ',' '$6 == "VAN" {print $2}' analysis/tractseg/index_tracts.csv)
echo $VAN_list # UF_left UF_right ST_FO_left ST_FO_right


# Combine masks per network


# Declare an array to hold the tract filenames
tracts=($FPN_list)

# Initialize the network mask
network_mask="analysis/tractseg/tck2fixel_output/network_mask_FPN_final.mif"

# Loop over each tract in the list
for ((i=0; i<${#tracts[@]}; i++)); do
    # Construct the filepath for the current tract
    tract_path="analysis/tractseg/tck2fixel_output/individual_tracts/${tracts[i]}.mif"

    # If it's the first tract, just copy it to the network mask
    if [ $i -eq 0 ]; then
        cp "$tract_path" "$network_mask"
    else
        # Otherwise, combine the current tract with the existing network mask
        mrcalc "$tract_path" "$network_mask" -max "analysis/tractseg/network_mask_FPN_$((i+1)).mif"
        # Update the network mask to the latest one
        network_mask="analysis/tractseg/network_mask_FPN_$((i+1)).mif"
    fi
done

# Rename the final network mask 
mv "$network_mask" "analysis/tractseg/tck2fixel_output/network_mask_FPN_final.mif"

# view fixel map of tract by opening index file in fixel plot overlay and threshold >0.95
# use ${tract}.mif for fixelcfestats option -mask
