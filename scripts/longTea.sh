#!/bin/bash

#SBATCH --partition=bch-compute
#SBATCH --job-name=mingyun
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=5G

source /home/ch252274/work2/miniconda/etc/profile.d/conda.sh
conda activate mingyun_te

python longTea.py
