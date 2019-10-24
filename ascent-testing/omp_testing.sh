#!/bin/bash

cd /gpfs/wolf/gen135/proj-shared/BaryTree/profiling/

N=100000
THETA=0.0
ORDER=8
CLUSTERSIZE=4000          
BATCHSIZE=4000
 
KERNEL=0
KAPPA=0.0

NZ=1
NY=1
NX=1
NP=$(($NX * $NY * $NZ))
echo $NP
SOURCES=/gpfs/wolf/proj-shared/gen135/BaryTree/randomPoints/S${N}_${NX}x_${NY}y_${NZ}z.bin
TARGETS=/gpfs/wolf/proj-shared/gen135/BaryTree/randomPoints/T${N}_${NX}x_${NY}y_${NZ}z.bin
OFFSETS=/gpfs/wolf/proj-shared/gen135/BaryTree/randomPoints/offsets${N}_${NX}x_${NY}y_${NZ}z.bin
DIRECT=/gpfs/wolf/proj-shared/gen135/BaryTree/randomPoints/direct_gpu_${N}.bin
OUTPUT=/gpfs/wolf/proj-shared/gen135/BaryTree/randomPoints/njv_tree_gpu_${N}.csv
#jsrun -n ${NP} -a 1 -c 1 -g 1 ~/.local/bin/tree-distributed-cpu $SOURCES $TARGETS $OFFSETS $OFFSETS $DIRECT $OUTPUT $N $N $THETA $ORDER $CLUSTERSIZE $BATCHSIZE $KERNEL $KAPPA
jsrun -n ${NP} -a 1 -c 1 -g 1 ~/.local/bin/tree-distributed-gpu $SOURCES $TARGETS $OFFSETS $OFFSETS $DIRECT $OUTPUT $N $N $THETA $ORDER $CLUSTERSIZE $BATCHSIZE $KERNEL $KAPPA
jsrun -n ${NP} -a 1 -c 1 -g 1 pgprof -o "omp_prof_%p.nvvp" ~/.local/bin/tree-distributed-gpu-omp $SOURCES $TARGETS $OFFSETS $OFFSETS $DIRECT $OUTPUT $N $N $THETA $ORDER $CLUSTERSIZE $BATCHSIZE $KERNEL $KAPPA

