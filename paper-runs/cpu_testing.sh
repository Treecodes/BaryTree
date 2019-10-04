#!/bin/bash

 
DIRECT=/Users/nathanvaughn/Desktop/randomPoints/directSum_cpu_Coulomb.csv 



ORDER=8
THETA=0.8
BATCHSIZE=400  
CLUSTERSIZE=400


## COULOMB 
KAPPA=0.0
KERNEL=0
OUTPUT=/Users/nathanvaughn/Desktop/randomPoints/cpu_Coulomb.csv 

NX=$1
NY=$2
NZ=$3
NP=$(($NX * $NY * $NZ))

echo $NX
echo $NY
echo $NZ
echo $NP
for N in 10000 
do
	echo N=$N 
	SOURCES=/Users/nathanvaughn/Desktop/randomPoints/S${N}_${NX}x_${NY}y_${NZ}z.bin    
	TARGETS=/Users/nathanvaughn/Desktop/randomPoints/T${N}_${NX}x_${NY}y_${NZ}z.bin
	OFFSETS=/Users/nathanvaughn/Desktop/randomPoints/offsets${N}_${NX}x_${NY}y_${NZ}z.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/Users/nathanvaughn/Desktop/randomPoints/ex_st_coulomb_$N.bin
	mpirun -np $NP direct-distributed-cpu   $SOURCES $TARGETS $DIRECT $OUTPUT $N $N $KERNEL $KAPPA
	mpirun -np $NP tree-distributed-cpu   $SOURCES $TARGETS $OFFSETS $OFFSETS $DIRECT $OUTPUT $N $N $THETA $ORDER $CLUSTERSIZE $BATCHSIZE $KERNEL $KAPPA
done  


