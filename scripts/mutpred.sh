#!/bin/bash
#SBATCH --job-name=diamond_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -t 1-00:00:00
#SBATCH --mem=12G


DIR="/sciclone/scr10/abharathan01/PKU/scripts/mutpred2.0"

${DIR}/run_mutpred2.sh -i wildtype.fasta -o mutpred2_results.txt
