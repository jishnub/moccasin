#!/bin/bash
#SBATCH -p parallel
#SBATCH -o norms_output
#SBATCH --ntasks=250
#SBATCH --time=48:00:00

module purge
source .moduleloader

job="./analyze --compute_norms --instrument HMI --yearnum 1 --nyears 0.2";
ellmax=255

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog $SCRATCH/logs/norms_status.log"

find $SCRATCH/logs -name norms_output.tar -delete

$parallel "srun --exclusive -n1 -N1 $job --ell {1} --ellp {1} &> $SCRATCH/logs/norm_{1}.log && \
tar -C $SCRATCH/logs -uf $SCRATCH/logs/norms_output.tar norm_{1}.log &&  \
rm  $SCRATCH/logs/norm_{1}.log" ::: $(seq $ellmax)