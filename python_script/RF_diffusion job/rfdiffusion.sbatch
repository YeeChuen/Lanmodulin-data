#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 8 processor core(s) per node 
#SBATCH --mem=64G   # maximum memory per node
#SBATCH --gres=gpu:a100:1

#SBATCH --job-name="EF2motif_8FNS_gpu"
#SBATCH --output="EF2motif_8FNS_gpu%j" # job standard output file (%j replaced by job id)
#SBATCH --mail-user=curwen@iastate.edu   #!!! CHANGE to your next_ID
#SBATCH --mail-type=BEGIN,END,FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


cd /work/ratul1/curwen/Lanmodulin-data/python_script/RF_diffusion job

# load SE3nv environment for rfdiffusion
module load micromamba
eval "$(micromamba shell hook --shell=bash)"
export MAMBA_ROOT_PREFIX=/work/ratul1/chuen/rfdiffusion/micromamba

micromamba env list

micromamba activate SE3nv

# chain = A
# 29 - 58 (30) <> EF2 = 59 - 70 (12) <> 71 - 133 (63) 
# 29 - 83 (55) <> EF3 = 84 - 95 (12) <> 96 - 133 (38) 
#'contigmap.contigs=[A29-58/A59-70/A71-133]
#'contigmap.contigs=[A29-83/A84-95/A96-133] 

# running EF2motif
./run_inference.py inference.output_prefix=example_outputs/EF2motif_8FNS_design inference.input_pdb=examples/input_pdbs/8FNS.pdb 'contigmap.contigs=[25-35/A59-70/58-68]' inference.num_designs=100 