#!/bin/bash
#SBATCH --time=10:15:00             				# Time limit for the job (REQUIRED).
#SBATCH --job-name=NF_Parent_Rescue_Camo_Vars      # Job name
#SBATCH --ntasks=1                  				# Number of cores for the job. Same as SBATCH -n 1
#SBATCH --mem=5G                    				# Total memory requested
#SBATCH --partition=normal          				# Partition/queue to run the job in. (REQUIRED)
#SBATCH -e slurm-%j.err             				# Error file for this job.
#SBATCH -o slurm-%j.out             				# Output file for this job.
#SBATCH -A coa_mteb223_uksr 	    				# Project allocation account name (REQUIRED)

module load ccs/java/jdk1.8.0_202
nextflow run rescue_camo_variants.nf
