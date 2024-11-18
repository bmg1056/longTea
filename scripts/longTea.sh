#!/bin/bash

#SBATCH --partition=bch-compute
#SBATCH --job-name=mingyun
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=5G
#SBATCH --array=0-17%18

source /home/ch252274/work2/miniconda/etc/profile.d/conda.sh
conda activate mingyun_te

arr=("l1-HG002" "l1-HG005" "l1-HG00438" "l1-HG02257" "l1-HG02486" "l1-HG02622" "alu-HG002" "alu-HG005" "alu-HG00438" "alu-HG02257" "alu-HG02486" "alu-HG02622" "sva-HG002" "sva-HG005" "sva-HG00438" "sva-HG02257" "sva-HG02486" "sva-HG02622")

#arr=("l1-HG002" "alu-HG002" "sva-HG002")
input_string=${arr[$SLURM_ARRAY_TASK_ID]}

part1=$(echo $input_string | cut -d'-' -f1)
part2=$(echo $input_string | cut -d'-' -f2)

python longTea.py ${part1} ${part2}
