#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -p development
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J rocetest

set -u

./run.sh "$WORK/data/r-oce/x00.input.csv"
