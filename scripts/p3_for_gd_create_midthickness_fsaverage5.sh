# Script to create a template midthickness surface in fsaverage space (in order to later compute a template Schaefer-400 GD matrix)

# This script follows what is done by Micapipe in 02_post-structural.sh ("Build the fsLR-32k sphere and midthickness surface"), i.e. using wb_shortcuts -freesurfer-resample-prep (which itself uses wb_command -surface-average, so taking the average of white and pial -NOT equivolumetric surfaces) - see https://github.com/Washington-University/wb_shortcuts/blob/master/wb_shortcuts

# Date creation: 06.11.2025, Author: Bianca Serio


# Before running script: Activate FreeSurfer environment -> FREESURFER (to use mris_convert)



#!/bin/bash

# Path to directory containing fsaverage5 white and pial surfaces
PATH_SURF=/data/pt_02667/data/ABCD/micapipe_geodesic_distances/surfaces/fsaverage5

for hemi in lh rh; do
    echo "=== Processing hemisphere: $hemi ==="

    # Convert FreeSurfer .surf files to GIFTI (if not already done)
    if [ ! -f ${PATH_SURF}/${hemi}.white.surf.gii ]; then
        mris_convert ${PATH_SURF}/${hemi}.white ${PATH_SURF}/${hemi}.white.surf.gii
    fi

    if [ ! -f ${PATH_SURF}/${hemi}.pial.surf.gii ]; then
        mris_convert ${PATH_SURF}/${hemi}.pial ${PATH_SURF}/${hemi}.pial.surf.gii
    fi

    # Compute the midthickness (geometric mean of white and pial)
    wb_command -surface-average \
        ${PATH_SURF}/${hemi}.midthickness.surf.gii \
        -surf ${PATH_SURF}/${hemi}.white.surf.gii \
        -surf ${PATH_SURF}/${hemi}.pial.surf.gii
done
