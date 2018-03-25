#!/bin/bash
#
#SBATCH --job-name=AlgoElectroScaling
#SBATCH --output=AlgoElectroScaling.output
#
#SBATCH --time=30:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

VAR=$((echo "### NBR_MPI = $SLURM_NTASKS AND OMP = $SLURM_CPUS_PER_TASK ###"
echo "### input is $INPUTFILE###"
time mpirun ./main -inputfile $INPUTFILE) 2>&1)
flock -x -w 500 AlgoElectroScaling_outfile echo "$VAR" >> AlgoElectroScaling_outfile
