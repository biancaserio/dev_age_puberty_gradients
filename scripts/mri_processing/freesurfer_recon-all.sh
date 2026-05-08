# Step 1: Autorecon1 with T1w and T2w, skipping skull stripping and pial refinement
recon-all \ -autorecon1 \ -i $INPUT_T1W \ -T2 $INPUT_T2W \ -noskullstrip \ -noT2pial \ -noFLAIRpial \ -openmp 4 \ -subjid $SUBJECT_ID \ -sd $SUBJECTS_DIR

# Step 2: Volume-only autorecon2
recon-all \ -autorecon2-volonly \ -openmp 4 \ -subjid $SUBJECT_ID \ -sd $SUBJECTS_DIR

# Step 3: T2-based pial surface refinement and cortical ribbon generation
recon-all \ -parallel \ -T2pial \ -cortribbon \ -subjid $SUBJECT_ID \ -sd $SUBJECTS_DIR

# Step 4: Selective autorecon3 (with many post-surface steps skipped)
recon-all \ -autorecon3 \ -openmp 4 \ -nosphere \ -nosurfreg \ -nojacobian_white \ -noavgcurv \ -nocortparc \ -nopial \ -noparcstats \ -nocortparc2 \ -noparcstats2 \ -nocortparc3 \ -noparcstats3 \ -nopctsurfcon \ -nocortribbon \ -nobalabels \ -subjid $SUBJECT_ID \ -sd $SUBJECTS_DIR

