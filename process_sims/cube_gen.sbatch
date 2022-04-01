#!/bin/bash
#SBATCH -J cube_gen
#SBATCH -N 1 -n 4
#SBATCH -o cube_gen.out
#SBATCH --mail-type=END
#SBATCH --mail-user=panopg@caltech.edu
#SBATCH -t 24:00:00

python3.7 cube_gen.py

