#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=10000M
#SBATCH --account=rrg-cgreenwo
#SBATCH --mail-user=kaiqiong.zhao@mail.mcgill.ca
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --array=1-500

module load r
echo Started: 
  date
echo $SLURM_ARRAY_TASK_ID
Rscript "/project/6007480/zhaokq/sparse_Simu/sparseSOMNiBUS_Simu_1000SNPs_new_cv/Exp_5_1000_rho0_nknots_10/run_100snps.R" $SLURM_ARRAY_TASK_ID
echo Ended
date