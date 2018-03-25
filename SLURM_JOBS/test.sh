rm -f -- AlgoElectroScaling_outfile
for NBRMPI in {1..10}
do
	for NBROMP in {1..16}
        do
                for i in {1..2}
                do
                        NS=$(squeue | grep rtoma)
                        NS=$(echo "$NS" | wc -l)
                        MAX=2
                        while [ "$NS" -gt "$MAX" ]
                        do
                                sleep 30
                                NS=$(squeue | grep rtoma)
                                NS=$(echo "$NS" | wc -l)
                                echo "Number of slurm tasks is $NS (max $MAX)"
                        done
						STR="TESTS/testSourceCenteredInCube.input"
                        echo "### NBR_MPI = $NBRMPI AND NBR_OMP = $NBROMP ### $STR"
						TOTMEM=40000
						# We know the total size of the problem, be smart !
						MEMPERCPU=$((TOTMEM / NBROMP))
						MEMPERCPU=$((MEMPERCPU / NBRMPI))
						echo "    MEM_PER_CPU is $MEMPERCPU (TOT $TOTMEM)"
                        sbatch --export=INPUTFILE=$STR --ntasks $NBRMPI --cpus-per-task $NBROMP --mem-per-cpu MEMPERCPU AlgoElectroScaling.sh
                done
        done

done
