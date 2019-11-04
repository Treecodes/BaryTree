#!/bin/bash


TREETYPE=1
THETA=0.8
ORDER=7
CLUSTERSIZE=100
BATCHSIZE=100
KAPPA=0.5

COMPAREDIRECT=1

for N in 10000 
do
    for NP in 1
    do
        for KERNEL in coulomb #yukawa coulomb_SS yukawa_SS
        #for KERNEL in yukawa_SS coulomb_SS
        do
            echo $NP
            echo $KERNEL
            SOURCES=/Users/nathanvaughn/Desktop/randomPoints/S${N}_${NX}x_${NY}y_${NZ}z.bin
            TARGETS=/Users/nathanvaughn/Desktop/randomPoints/T${N}_${NX}x_${NY}y_${NZ}z.bin
            OFFSETS=/Users/nathanvaughn/Desktop/randomPoints/offsets${N}_${NX}x_${NY}y_${NZ}z.bin
            DIRECT=/Users/nathanvaughn/Desktop/randomPoints/SS_KI_direct_gpu_${N}.bin
            OUTPUT=/Users/nathanvaughn/Desktop/randomPoints/SS_KI_tree_gpu_${N}.csv
            /usr/local/bin/mpirun -n ${NP} zoltan_example_cpu $N $ORDER $THETA $CLUSTERSIZE $BATCHSIZE $KERNEL $KAPPA $TREETYPE $COMPAREDIRECT 
            done
        done   
done    
