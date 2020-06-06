#!/bin/bash
#SBATCH --partition=standard
#SBATCH --account=krasny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem-per-cpu=2g
#SBATCH --time=24:00:00
#SBATCH --job-name="c-vs-g"
#SBATCH --output="dataFiles/cpu--%j.%N.out"

TREETYPE=1
THETA=0.8
ORDER=8
CLUSTERSIZE=2000
BATCHSIZE=2000
KAPPA=0.0
SIZECHECKFACTOR=1.0

APPROXIMATIONNAME=lagrange
COMPAREDIRECT=1
VERBOSITY=0
N=1000000
KERNEL=coulomb
SINGULARITYHANDLING=skipping

for THETA in 0.9 0.7 0.5
do
    for ORDER in 2 4 6 8 10 12 14
    do
        echo "LAGRANGE"
        mpirun -np 6 zoltan_example_cpu $N $ORDER $THETA $CLUSTERSIZE $BATCHSIZE $KERNEL $KAPPA $SINGULARITYHANDLING lagrange $TREETYPE 1.0 $COMPAREDIRECT $VERBOSITY
        echo "HERMITE"    
        mpirun -np 6 zoltan_example_cpu $N $ORDER $THETA $CLUSTERSIZE $BATCHSIZE $KERNEL $KAPPA $SINGULARITYHANDLING hermite $TREETYPE 4.0 $COMPAREDIRECT $VERBOSITY
    done
done
