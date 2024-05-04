#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 8 processor core(s) per node 
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --gres=gpu:a100:1

#SBATCH --job-name="esmfold"
#SBATCH --output="esmfold_%j" # job standard output file (%j replaced by job id)
#SBATCH --mail-user=chuen@iastate.edu   #!!! CHANGE to your next_ID
#SBATCH --mail-type=BEGIN,END,FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# cd to your project directory
cd /work/ratul1/chuen/rfdiffusion/RFdiffusion

# load esmfold
module load esmfold

# run esmfold script
esmfold_inference.py -i esmfold_example.fasta  -o esmfold