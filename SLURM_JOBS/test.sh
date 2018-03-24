rm -f -- AlgoElectroScaling_outfile
for NBRMPI in {1..10}
do
	for NBROMP in {1..16}
        do
                for i in {1..1}
                do
                        NS=$(squeue | grep rtoma)
                        NS=$(echo "$NS" | wc -l)
                        MAX=8
                        while [ "$NS" -gt "$MAX" ]
                        do
                                sleep 30
                                NS=$(squeue | grep rtoma)
                                NS=$(echo "$NS" | wc -l)
                                echo "Number of slurm tasks is $NS (max $MAX)"
                        done
						STR="TESTS/testSourceCenteredInCube.input"
                        echo "### NP = $NBRMPI AND K = $NBROMP ### $STR"
                        sbatch --export=OMP_NBR_THRDS=$NBROMP,INPUTFILE=$STR --ntasks $NBRMPI AlgoElectroScaling.sh
                done
        done

done
