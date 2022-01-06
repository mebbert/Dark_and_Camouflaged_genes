#!/bin/bash

#SBATCH --time=23:00:00
#SBATCH --job-name=Testing_step_03
#SBATCH --ntasks=1
#SBATCH --mem=800G
#SBATCH --partition=jumbo
#SBATCH -A coa_mteb223_uksr
#SBATCH --mail-user=mlps241@uky.edu
#SBATCH --mail-type=END,FAIL


singularity exec /project/mteb223_uksr/singularity_files/rescue_camo_variants_withR.sif bash calc_median_depth.ogs $1 $2 $3
