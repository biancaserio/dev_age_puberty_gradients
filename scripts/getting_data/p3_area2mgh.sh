# Script that converts native FreeSurfer surfaces .area to .mgh format (by hemisphere) for total surface area computation - created by Bianca Serio on 11.02.25 
# Date amended: 16.05.2025
# # added check "if sub exists in output dir (and contains the .mgh files), skip it"

# Steps to run script: 
# (0. datalad clone ABCD data freesurfer output)
# (0. create an output folder for surface area files)
# 1. ssh on MPI server -> getserver -sL
# 2. activate ssh agent to avoid needing to give pwd for datalad get each time ->  eval "$(ssh-agent -s)" then ssh-add ~/.ssh/id_ed25519_mpi 
# 3. activate FreeSurfer to use commands (mri_convert) -> FREESURFER

# 4. run script: bash p3_area2mgh.sh /data/p_02667/development/scripts/sub_baseline_qc_val_relabelled.txt /data/pt_02667/data/ABCD/ABCD_freesurfer /data/pt_02667/data/ABCD/ABCD_baseline_SA

# !!! CAREFUL - format of subject list must be the same as in directories, i.e., “sub-NDARINVXXX” (hence "relabelled") !!!

# The script requires 3 arguments (if fewer than 3 arguments are given, the script prints usage instructions and exits)
# 1. `SUBJECT_LIST`: File containing subject IDs (one per line)
# 2. `SUBJECTS_DIR`: Directory containing FreeSurfer subject data
# 3. `OUTPUT_DIR`: Output directory


# FYI FreeSurfer data needed to run this script: /surf/lh.area & /surf/rh.area



#!/bin/bash

print_help() {
    echo "Usage: $0 <SUBJECT_LIST> <SUBJECTS_DIR> <OUTPUT_DIR>"
}

#########################################################################
# Input handling
if [ $# -lt 3 ]; then    
    echo "Not enough arguments"
    print_help
    exit 1
fi

# Variables
SUBJECT_LIST="$1"
SUBJECTS_DIR="$2"
OUTPUT_DIR="$3"
ICO="$4"

FAILED_LOG_FILE="${OUTPUT_DIR}/failed_subjects_log.txt"
SKIPPED_LOG_FILE="${OUTPUT_DIR}/skipped_subjects_log.txt"

# Clear previous failed subjects log if it exists
> "$FAILED_LOG_FILE"
> "$SKIPPED_LOG_FILE"

# Read subjects from input file
mapfile -t SUBJECTS < "$SUBJECT_LIST"
total_subjects=${#SUBJECTS[@]}
failed_subjects=()
skipped_subjects=()
counter=1


# Start time
begin=$(date +%s.%N)

############ Starting script ###########
echo -e "\e[0;44m\n[INFO] Convert .area to .mgh\n \e[0m"

for sub in "${SUBJECTS[@]}"; do
    subject_dir="${SUBJECTS_DIR}/${sub}"

    # Check if the subject directory already exists in the output directory (with .mgh files)
    subject_outdir="${OUTPUT_DIR}/${sub}"
    if [ -e "${subject_outdir}/lh.native.area.mgh" ] && [ -e "${subject_outdir}/rh.native.area.mgh" ]; then
        echo ".mgh output already exists for $sub at $subject_outdir. Skipping subject."
        skipped_subjects+=("$sub")
        echo "$sub" >> "$SKIPPED_LOG_FILE"
        continue
    fi
    

    # Check if subject directory exists
    if [ -d "$subject_dir" ]; then
        echo -e "\e[0;44m\n[PROCESSING] $sub ($counter / $total_subjects)\n\e[0m"

        cd "$subject_dir" || { echo "[ERROR] Failed to enter directory: $subject_dir"; continue; }

        # Perform datalad get operations
        datalad get --no-data .
        if datalad get -J 2 surf/*h.area; then
            subject_outdir="${OUTPUT_DIR}/${sub}"
            mkdir -p "$subject_outdir"

            # Convert FreeSurfer .area to .mgh
            if ! mri_convert -i "surf/lh.area" -o "${subject_outdir}/lh.native.area.mgh"; then
                echo "[ERROR] Conversion failed for: $sub (lh.area)"
                failed_subjects+=("$sub")
                echo "$sub" >> "$FAILED_LOG_FILE"
                continue
            fi

            if ! mri_convert -i "surf/rh.area" -o "${subject_outdir}/rh.native.area.mgh"; then
                echo "[ERROR] Conversion failed for: $sub (rh.area)"
                failed_subjects+=("$sub")
                echo "$sub" >> "$FAILED_LOG_FILE"
                continue
            fi

            # Drop the downloaded data
            datalad drop --what filecontent --reckless kill
        else
            echo "[ERROR] Missing surf/*h.area files for $sub. Skipping."
            failed_subjects+=("$sub")
            echo "$sub" >> "$FAILED_LOG_FILE"
        fi
    else
        echo "[ERROR] Subject directory not found: $subject_dir"
        failed_subjects+=("$sub")
        echo "$sub" >> "$FAILED_LOG_FILE"
    fi
    ((counter++))

done



# End time
end=$(date +%s.%N)
duration=$(echo "$end - $begin" | bc)

echo -e "\033[38;5;220m \n TOTAL running time: ${duration} seconds \n \033[0m"

# Report failed subjects if any
if [ ${#failed_subjects[@]} -ne 0 ]; then
    echo -e "\e[0;41m[WARNING] The following subjects failed:\e[0m"
    printf '%s\n' "${failed_subjects[@]}"
    echo "Failed subjects saved in: $FAILED_LOG_FILE"
fi

# Report skipped subjects if any
if [ ${#skipped_subjects[@]} -ne 0 ]; then
    echo -e "\e[0;41m[WARNING] The following subjects have been skipped:\e[0m"
    printf '%s\n' "${skipped_subjects[@]}"
    echo "Skipped subjects saved in: $SKIPPED_LOG_FILE"
fi
