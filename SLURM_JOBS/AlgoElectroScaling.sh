#!/bin/bash
#
#SBATCH --job-name=AlgoElectroScaling
#SBATCH --output=AlgoElectroScaling.output
#
#SBATCH --mem-per-cpu=5000
#SBATCH --time=5:00

export OMP_NUM_THREADS=$OMP_NBR_THRDS

VAR=$((echo "### NBR_MPI = $SLURM_NTASKS AND OMP = $OMP_NUM_THREADS ###"
time mpirun ./main -inputfile $INPUTFILE -whichAlgo ELECTRO) 2>&1)
flock -x -w 500 AlgoElectroScaling_outfile echo "$VAR" >> AlgoElectroScaling_outfile
