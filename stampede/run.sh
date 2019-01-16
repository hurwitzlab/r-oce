#!/bin/bash

#SBATCH -J r-oce
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -A iPlant-Collabs

module load tacc-singularity

IMG=/work/05066/imicrobe/singularity/r-oce-0.0.1.img

set -u

function HELP() {
    printf "Usage:\\n  %s IN_DIR \\n\\n" "$(basename "$0")"
    exit 0
}

[[ $# -ne 1 ]] && HELP

singularity exec $IMG $1

echo "Done."
echo "Comments to Ken Youens-Clark <kyclark@email.arizona.edu>"
