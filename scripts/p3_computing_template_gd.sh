# Script to create a template GD matrix in Schaefer-400 space using midthickness surface in fsaverage space

# Date creation: 06.11.2025, Author: Bianca Serio

# Before running script: activate conda miniforge gd enviornment in order to be able to use python package pygeodesic (in geoDistMapper.py): camf_gd

#!/bin/bash


Info() { echo "[INFO] $*"; }
Do_cmd() { echo "[CMD] $*"; eval "$*"; }


# Define paths
PATH_SURF="/data/pt_02667/data/ABCD/micapipe_geodesic_distances/surfaces"
PATH_PARC="/data/pt_02667/data/ABCD/micapipe_geodesic_distances/parcellations"
OUT_PATH="/data/pt_02667/data/ABCD/micapipe_geodesic_distances/template_GD"
PATH_SCRIPTS="/data/p_02667/development/scripts/micapipe"

# creates output directory if it doesn't exist
[[ ! -d "$OUT_PATH" ]] && mkdir -p "$OUT_PATH"


#	Timer
aloita=$(date +%s)
Nsteps=0
N=0


# Surfaces and Templates downloaded at:
# - https://github.com/MICA-MNI/micapipe/tree/master/surfaces
# - https://github.com/MICA-MNI/micapipe/tree/master/parcellations



#------------------------------------------------------------------------------#


### Compute geodesic distance on Schaefer 400 parcellation (template (parcellation) space) using template midthickness surface in fsaverage space

# these mithickness surfaces were computed by me as average between white and pial templace surfaces (using p3_forGD_create_midthickness_fsaverage5.sh)
lh_midsurf_fsaverage5="${PATH_SURF}/fsaverage5/lh.midthickness.surf.gii"
rh_midsurf_fsaverage5="${PATH_SURF}/fsaverage5/rh.midthickness.surf.gii"


parc="schaefer-400"
lh_annot="${PATH_PARC}/lh.${parc}_mics.annot"
rh_annot="${PATH_PARC}/rh.${parc}_mics.annot"
outName="${OUT_PATH}/template_atlas-${parc}_GD"

if [ -f "${outName}.shape.gii" ]; then
    Info "Geodesic Distance on $parc, already exists"; ((Nsteps++))
else
    Info "Computing Geodesic Distance from $parc"
    Do_cmd "$PATH_SCRIPTS"/geoDistMapper.py -lh_surf "$lh_midsurf_fsaverage5" -rh_surf "$rh_midsurf_fsaverage5" -outPath "$outName" \
            -lh_annot "$lh_annot" -rh_annot "$rh_annot" -parcel_wise
    if [[ -f "${OUT_PATH}.shape.gii" ]]; then ((Nsteps++)); fi
fi

Do_cmd rm -rf "${OUT_PATH}"/*.func.gii  # this file gets created in the process but we don't need it -> remove


#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)

eri=$((lopuu - aloita))
eri_min=$((eri / 60))
echo "Elapsed time: ${eri_min} minutes"


