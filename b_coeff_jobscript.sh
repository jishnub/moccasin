#!/bin/bash
#SBATCH -p parallel
#SBATCH -o b_coeff_output
#SBATCH --ntasks=30
#SBATCH --time=48:00:00

module purge
source .moduleloader

find $SCRATCH/logs -name b_output.tar -delete

job_ell_dl() {
	ell=$1
	deltaell=1
	ellp=$(python -c "print($ell+$deltaell)")
	job="./analyze --instrument HMI --yearnum 1 --nyears 0.2";
	srun --exclusive -n1 -N1 $job --ell $ell --ellp $ellp &> $SCRATCH/logs/b_$ell_$ellp.log && \
	tar -C $SCRATCH/logs -uf $SCRATCH/logs/b_output.tar b_$ell_$ellp.log &&  \
	rm  $SCRATCH/logs/b_$ell_$ellp.log
}
export -f job_ell_dl

ellmax=10
parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog $SCRATCH/logs/b_status.log"

$parallel job_ell_dl ::: $(seq $ellmax)