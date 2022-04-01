#!/bin/bash
#SBATCH -J conver2
#SBATCH -N 1 -n 4
#SBATCH -o converter2.out
#SBATCH --mail-type=END
#SBATCH --mail-user=panopg@caltech.edu
#SBATCH -t 72:00:00

python3.7 converter_octree.py

