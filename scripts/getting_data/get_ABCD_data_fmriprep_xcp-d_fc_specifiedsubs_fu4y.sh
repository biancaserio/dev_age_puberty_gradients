# Script for postprocessing ABCD fMRIprep **4y follow-up** files using XCP-D, which also outputs FC matrices, based on input list from a txt file
# Copied from fu2y version to tweak to make specific to 4y follow-up
# Date creation: 07.03.2025, Author: Bianca Serio
# Date amended: 15.05.2025
# # added check "if sub exists in XCP-D output folder, skip it"

# steps to run script: 
# (0. datalad clone ABCD data -> pt_02667/data/ABCD/ABCD_fmriprep)
# 1. check that directories defined in script match
# 2. ssh on MPI server -> getserver -sL
# 3. activate ssh agent to avoid needing to give pwd for datalad get each time ->  eval "$(ssh-agent -s)" then ssh-add ~/.ssh/id_ed25519_mpi 
# 4. activate miniforge environment to use antsApplyTransform -> . /data/u_serio_software/miniforge3/etc/profile.d/conda.sh;conda activate base
# 4. activate FSL environment to generate BOLD reference using fslroi command -> FSL
# 5. run script: bash get_ABCD_data_fmriprep_xcp-d_fc_specifiedsubs_fu4y.sh /data/p_02667/development/scripts/sub_3_timepoints_qc_val_relabelled.txt
# !!! CAREFUL - format of subject list must be the same as in directories, i.e., “sub-NDARINVXXX” (hence "relabelled") !!!


#!/bin/bash

# Input txt file path - contains a list of subjects (1 column, no header)
subject_list_file="$1"
failed_subjects_file="/data/pt_02667/data/ABCD/XCP-D_output_fu4y/logs/failed_subjects_get_fmriprep.txt"
missing_files_log="/data/pt_02667/data/ABCD/XCP-D_output_fu4y/logs/missing_files_log.txt"

# Read the txt file and create an array of subjects
mapfile -t subs_data_to_download < "$subject_list_file"

# Get the total number of subjects
total_subjects=${#subs_data_to_download[@]}

# Initialize arrays to store failed subject IDs or missing file errors
failed_subjects=()
missing_files_errors=()

# Participant directory
pptDIR="/data/pt_02667/data/ABCD/ABCD_fu4y_fmriprep22.1.1_FSready/fmriprep/"

# Log directory base path
log_base="/data/pt_02667/data/ABCD/XCP-D_output_fu4y/logs"

# Eventname variable
eventname="ses-4YearFollowUpYArm1"

# Counter for tracking the current subject number
counter=0

# Loop over each subject in inputted list
for sub in "${subs_data_to_download[@]}"; do
    # Increment the subject counter
    counter=$((counter + 1))
    
    subject_dir="${pptDIR}${sub}"
    subject_logs="${log_base}/${sub}"
    mkdir -p "$subject_logs"  # Create log directory for the subject

    # Check if the subject directory already exists in the XCP-D output directory
    xcpd_output_dir="/data/pt_02667/data/ABCD/XCP-D_output_fu4y/${sub}"
    if [ -d "$xcpd_output_dir" ]; then
        echo "XCP-D output already exists for $sub at $xcpd_output_dir. Skipping subject." | tee -a "$subject_logs/skipped_subjects.txt"
        continue
    fi
    
    # Check if the subject directory exists
    if [ -d "$subject_dir" ]; then
    	
    	# Change to the subject directory or log failure
        if cd "$subject_dir"; then
            echo "------------------Processing data for $sub------------------"
            echo "Subject number: $counter out of $total_subjects"
            
            # Perform datalad get operations to download the specific files needed for XCP-D pipeline to run (specific anat, all func to be safe)
            datalad get --no-data .
            datalad get -J 10 sub-*.html ${eventname}/func/* ${eventname}/anat/*desc-preproc_T1w* ${eventname}/anat/*desc-brain_mask* ${eventname}/anat/*from-MNI152NLin6Asym_to-T1w_mode-image_xfm.h5 ${eventname}/anat/*from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5
            
            # Get additional files required by XCP-D pipeline that are not available in fmriprep output
            cd ${eventname}/func/ || { 
                echo "Missing functional directory for $sub"; 
                missing_files_errors+=("$sub: Missing functional directory"); 
                continue; 
            }

            # Combined loop for symlink creation and BOLD reference generation
            for bold_file in *task-rest_run-*_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz; do
                # Strip unwanted characters from the path
                bold_file=$(echo "$bold_file" | sed 's/"//g;s/@$//')

                # Check if the file exists
                if [ -f "$bold_file" ]; then
                    # Step 1: Create symlink if it does not exist (desired file: Preprocessed BOLD data in MNI152NLin6Asym space)
                    symlink_file="${bold_file/_desc-smoothAROMAnonaggr/_desc-preproc}"                   
                    if [ ! -f "$symlink_file" ]; then
                    	ln -s "$bold_file" "$symlink_file" && {
                    	    echo "Symlink created: $symlink_file -> $bold_file" | tee -a "$subject_logs/successful_symlinks.txt"
                    	} || {
                    	    echo "Failed to create symlink for $bold_file" | tee -a "$subject_logs/failed_symlinks.txt"
                    	    continue
                    	}
                    else
                    	echo "Symlink already exists: $symlink_file" | tee -a "$subject_logs/successful_symlinks.txt"
                    fi
        
        

                    # Step 2: Generate BOLD reference if it does not exist (desired file: 3D BOLD-reference image in MNI152NLin6Asym space
                    boldref_file="${bold_file/_desc-smoothAROMAnonaggr_bold/_boldref}"
                    if [ ! -f "$boldref_file" ]; then
                    	fslroi "$bold_file" "$boldref_file" 0 1 && {
                    	    echo "running fslroi to create BOLD ref"
                    	    echo "BOLD reference created: $boldref_file" | tee -a "$subject_logs/successful_boldref.txt"
                    	} || {
                    	    echo "Failed to create BOLD reference for $bold_file" | tee -a "$subject_logs/failed_boldref.txt"

                    	    continue
                    	}
                    else
                    	echo "BOLD reference already exists: $boldref_file" | tee -a "$subject_logs/successful_boldref.txt"
                    fi
                else
                    echo "File not found: $bold_file" | tee -a "$subject_logs/missing_files.txt"
                fi
            done


            # 3. Apply transformations for brain mask (desired file: Mask in MNI152NLin6Asym space)
            cd ../anat/ || {
                echo "Missing anatomical directory for $sub";
                missing_files_errors+=("$sub: Missing anatomical directory");
                continue;
            }

            # Dynamically locate the brain mask file
            mask_file=$(ls "${sub}_${eventname}"*desc-brain_mask.nii.gz 2>/dev/null | head -n 1 | sed 's/"//g' | sed 's/@$//')
            if [ -n "$mask_file" ] && [ -f "$mask_file" ]; then
                mask_file=$(realpath "$mask_file") # Resolve symbolic link
                echo "Found brain mask file: $mask_file"
            else
                echo "Brain mask file not found for $sub"
                missing_files_errors+=("$sub: Brain mask file not found")
                continue
            fi

            # Dynamically locate the transform file
            transform_file=$(ls "${sub}_${eventname}"*from-T1w_to-MNI152NLin6Asym_mode-image_xfm.h5 2>/dev/null | head -n 1 | sed 's/"//g' | sed 's/@$//')
            if [ -n "$transform_file" ] && [ -f "$transform_file" ]; then
                transform_file=$(realpath "$transform_file") # Resolve symbolic link
                echo "Found transform file: $transform_file"
            else
                echo "Transform file not found for $sub"
                missing_files_errors+=("$sub: Transform file not found")
                continue
            fi

            output_mask="${sub}_${eventname}_space-MNI152NLin6Asym_desc-brain_mask.nii.gz"
            tpl_file="/data/pt_02667/data/ABCD/bianca/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-01_T1w.nii.gz"

            # Apply transformation if output mask does not already exist
            if [ ! -f "$output_mask" ]; then
                echo "Applying transformation..."
                antsApplyTransforms -d 3 -i "$mask_file" -r "$tpl_file" -n NearestNeighbor \
                    -t "$transform_file" -o "$output_mask" || {
                    echo "Error applying transformation for $mask_file"
                    missing_files_errors+=("$sub: Error applying transformation for $mask_file")
                    continue
                }
                echo "Transformation applied successfully!"
            else
                echo "Mask transformation output already exists for $mask_file"
            fi
            
            # Check if required BOLD data exists in fsLR space
            fsLR_bold_files=$(find "$subject_dir/${eventname}/func/" -name "*_space-fsLR_den-91k_bold.dtseries.nii" 2>/dev/null)
            if [ -z "$fsLR_bold_files" ]; then
                echo "No BOLD data found in fsLR space for $sub. Skipping XCP-D processing."
                missing_files_errors+=("$sub: No BOLD data in fsLR space")
            else
            	# Run XCP-D post-processing pipeline
            	echo "Running XCP-D pipeline"
            
            	singularity run \
            	-B /data/pt_02667/data/ABCD:/data/pt_02667/data/ABCD \
            	--cleanenv /data/u_serio_software/xcp_d-0.10.1.sif \
            	/data/pt_02667/data/ABCD/ABCD_fu4y_fmriprep22.1.1_FSready/fmriprep \
            	/data/pt_02667/data/ABCD/XCP-D_output_fu4y \
            	participant \
            	--mode 'none' \
            	--participant-label "${sub#sub-}" \
            	--bids-filter-file /data/pt_02667/data/ABCD/bianca/bids_filter_file_fu4y.json \
            	--nprocs 36 \
            	--input-type 'fmriprep' \
            	--file-format 'cifti' \
            	--dummy-scans 'auto' \
            	--despike 'y' \
            	--nuisance-regressors /data/pt_02667/data/ABCD/bianca/custom_confounds_24P_csf_wm.yaml \
            	--fd-thresh 0.3 \
            	--output-type 'censored' \
            	--combine-runs 'n' \
            	--smoothing 6 \
            	--motion-filter-type 'none' \
            	--head-radius 50 \
            	--lower-bpf 0.01 \
            	--upper-bpf 0.08 \
            	--bpf-order 2 \
            	--min-time 240 \
            	--atlases '4S456Parcels' \
            	--min-coverage 0.5 \
            	--create-matrices 240 all \
            	--work-dir /data/pt_02667/data/ABCD/work \
            	--warp-surfaces-native2std 'n' \
            	--abcc-qc 'n' \
            	--linc-qc 'y' \
            	1> "${subject_logs}/stdout.log" 2> "${subject_logs}/stderr.log" 
            fi
            
            # Drop all datalad files except .html and confounds files 
            cd "$subject_dir" || { echo "Failed to change directory to $subject_dir"; exit 1; }
            
            # Drop the downloaded data except .html files (in /${eventname}) and confound files (in /func)
            datalad drop --what filecontent --reckless kill -J 10 ${eventname}/anat/* ${eventname}/func/*_AROMA* ${eventname}/func/*_boldref* ${eventname}/func/*_desc-brain* ${eventname}/func/*_desc-MELODIC* ${eventname}/func/*_desc-preproc* ${eventname}/func/*_from* ${eventname}/func/*_hemi* ${eventname}/func/*_space* || { echo "Error dropping files for $sub"; }     
            
        else
            echo "Failed to cd into $subject_dir"
            failed_subjects+=("$sub")
        fi
    else
        echo "Directory $subject_dir does not exist"
        failed_subjects+=("$sub")
    fi
done

# Write the list of failed subjects to a text file
if [ ${#failed_subjects[@]} -ne 0 ]; then
    printf "%s\n" "${failed_subjects[@]}" > "$failed_subjects_file"
    echo "Failed subject IDs saved to $failed_subjects_file"
fi

# Write the missing files errors to a log
if [ ${#missing_files_errors[@]} -ne 0 ]; then
    printf "%s\n" "${missing_files_errors[@]}" > "$missing_files_log"
    echo "Missing files errors logged to $missing_files_log"
fi

