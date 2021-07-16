#!/bin/bash
#SBATCH --mem=32000
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --tasks-per-node=1
#SBATCH --partition=broadwl 
#SBATCH --job-name=run_xCell
#SBATCH --output=run_xCell.out

module load R

cd /project2/nobrega/grace/expression/obesity/cortex
Rscript /project2/nobrega/grace/expression/scripts/run_xCell.R /project2/nobrega/grace/expression/obesity/cortex/CMC_GTEx_xCell_tpm.txt
