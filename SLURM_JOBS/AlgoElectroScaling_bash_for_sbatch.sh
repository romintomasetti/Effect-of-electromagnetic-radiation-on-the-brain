#
# This bash file eases the scaling tests.
# The first loop is on the number of MPI processes.
# The second loop is on the number of OPENMP threads.
# The third loop is for statistics.
#
# Maximum simultaneous jobs:
MAXJOBS=5
# Go for it:
# for NBRMPI in {1..2}
# do
	# for NBROMP in {1..2}
	# do
			# for i in {1..2}
			# do
					# NS=$(squeue | grep rtoma)
					# NS=$(echo "$NS" | wc -l)
					# while [ "$NS" -gt "$MAXJOBS" ]
					# do
							# sleep 5
							# NS=$(squeue | grep rtoma)
							# NS=$(echo "$NS" | wc -l)
							# echo "Number of slurm tasks is $NS (max $MAX)"
					# done
					# echo "### NP = $NBRMPI AND K = $NBROMP ###"
					# sbatch --export=OMP_NBR_THRDS=$NBROMP,INPUTFILE="$INPUTFILE" --ntasks $NBRMPI AlgoElectroScaling.sh
			# done
	# done
# done
