#!/bin/bash
# Copy and unzip subject DWI data
# run in cluster: 
# chmod a+rwx scripts/dwi_copy_unzip.sh
# sbatch --time=08:00:00 --mem=10G Desktop/3022060.01/scripts/dwi_copy_unzip.sh

cd Desktop/3022060.01/
data_dir="data_download/data_2024-04-09"
working_dir="analysis/data"

# List all subject directories in working directory
subjects=$(find "$data_dir" -maxdepth 1 -type d -exec basename {} \;)
# subjects="HBU0A6F67FF578F01743"

# Loop through each subject directory
for subject in $subjects; do
    echo "Processing $subject"

    # Construct the path to the subject's directory
    subject_dir=$(find "${working_dir}" -type d -name "*_${subject}")
    if [ -z "$subject_dir" ]; then
        echo "Subject directory not found for $subject."
        continue
    fi

    # Create a directory for the copied files
    mkdir -p "${subject_dir}/dwi"

    # Copy the original file from the subject's data directory into their working directory
    if [ -e "${data_dir}/${subject}/mri_diffusion_1" ]; then
        cp "${data_dir}/${subject}/mri_diffusion_1" "${subject_dir}/dwi/"
    elif [ -e "${data_dir}/${subject}/mri_diffusion_2" ]; then
        cp "${data_dir}/${subject}/mri_diffusion_2" "${subject_dir}/dwi/"
    elif [ -e "${data_dir}/${subject}/mri_diffusion_3" ]; then
        cp "${data_dir}/${subject}/mri_diffusion_3" "${subject_dir}/dwi/"
    else
        echo "Skipping $subject. No DWI files found."
        continue
    fi

    # Change directory to the DWI folder
    cd "${subject_dir}/dwi/"

    # Add .zip to DWI filename 
    files=$(find . -type f -name 'mri_diffusion_*')
    for file in $files; do
        zip_file="${file}.zip"
        mv "$file" "$zip_file"
    done

    # Unzip 
    for zip_file in mri_diffusion_*.zip; do
        unzip -q "$zip_file"
    done

    # Delete .zip file
    rm -f mri_diffusion_*.zip

    # Return to the previous directory
    cd -

done
