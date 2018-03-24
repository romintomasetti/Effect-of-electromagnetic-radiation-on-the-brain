#!/bin/bash
#
#SBATCH --job-name=AlgoElectroSclaing
#SBATCH --output=AlgoElectroScaling.output
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=1:00

export OMP_NUM_THREADS=$OMP_NBR_THRDS

mpirun ./main -inputfile $INPUTFILE -whichAlgo ELECTRO