#!/bin/bash -l

# Define settings for a SLURM managed cluster:
#SBATCH --job-name=screen_pan_can_DEMO_genes_one # Job name
#SBATCH --nodes=1 # Nodes requested
#SBATCH --cores=10 # Cores requested
#SBATCH --mem=250G # max memory allocated
#SBATCH --output=Rout/%x.o.%j # Default is where script is run from
#SBATCH --mail-type=ALL # email if failed
#SBATCH --mail-user=ytakemon # email me
#SBATCH --error=Rout/%x.e.%j # Standard error
#SBATCH --time=01-00:00:0 # 13 hours

# Takes $genes are defined by `screen_pan_can_DEMO_genes_all.sh`
/gsc/software/linux-x86_64-centos7/R-4.0.2/lib64/R/bin/Rscript screen_pan_can_small_batch_one.R $genes
