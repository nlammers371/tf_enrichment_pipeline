#!/bin/bash
# Job name:
#SBATCH --job-name=savio_test
#
# Partition:
#SBATCH --partition=savio2
#
# QoS:
#SBATCH --qos=savio_normal
#
# Account:
#SBATCH --account=ac_mhmm
#
# Request one node:
#SBATCH --nodes=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Number of Processors per Node:
#SBATCH --ntasks-per-node=1
#
# Wall clock limit:
#SBATCH --time=72:00:00
#

## NL: This tells Savio to run 10 distinct instances of the job (essentially euivalent to "nBoots")
#SBATCH -a 1-10

## Command(s) to run:
module load matlab

# Make a temporary scratch directory for storing job
# and task information, to coordinate parallelization.
# This directory can then be referenced by assigning it to
# a 'parcluster.JobStorageLocation' property in your script.
mkdir -p /global/scratch/$USER/$SLURM_JOB_ID

## Call the inference function
matlab -nodisplay -nodesktop < ~/repos/tf_enrichment_pipeline/call_inference_savio.m