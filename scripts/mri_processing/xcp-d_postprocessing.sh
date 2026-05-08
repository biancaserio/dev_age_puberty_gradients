singularity run \
-B $SUBJECTS_DIR:$SUBJECTS_DIR \
--cleanenv $XCPD_DIR/xcp_d-0.10.1.sif \ $SUBJECTS_DIR/ABCD_fMRIprep/fmriprep \ $SUBJECTS_DIR/XCP-D_output \
participant \
--mode 'none' \
--participant-label $SUBJECT_ID \
--bids-filter-file bids_filter_file.json \--nprocs 36 \
--input-type 'fmriprep' \
--file-format 'cifti' \
--dummy-scans 'auto' \
--despike 'y' \
--nuisance-regressors custom_confounds_24P_csf_wm.yaml \
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
--work-dir $SUBJECTS_DIR/work \
--warp-surfaces-native2std 'n' \
--abcc-qc 'n' \
--linc-qc 'y' \
1> stdout.log 2> stderr.log
