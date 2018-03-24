# This bash file eases the scaling tests.
# The first loop is on the number of MPI processes.
# The second loop is on the number of OPENMP threads.
# The third loop is for statistics.
#
# Maximum simultaneous jobs:
MAX_JOBS=5
for NBR_MPI in {1..1}
do
        for NBR_OMP in {1..1}
        do
                for i in {1..1}
                do
                        NS=$(squeue | grep rtoma)
                        NS=$(echo "$NS" | wc -l)
                        MAX=8
                        while [ "$NS" -gt "$MAX_JOBS" ]
                        do
                                sleep 5
                                NS=$(squeue | grep rtoma)
                                NS=$(echo "$NS" | wc -l)
                                echo "Number of slurm tasks is $NS (max $MAX)"
                        done
                        echo "### NP = $NP AND K = $K ###"
                        sbatch --export=OMP_NBR_THRDS=$NBR_OMP,INPUTFILE=$INPUTFILE \
                                --ntasks $NBR_MPI \
                                AlgoElectroScaling.sh
                done
        done
done