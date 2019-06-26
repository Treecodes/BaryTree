#!/bin/bash




DS_CSV=../examplesOxygenAtom/local_coulomb.csv  

OUTFILE=/Users/nathanvaughn/Desktop/out.csv
THETA=0.7
ORDER=8
TREETYPE=1
MAXPARNODE=2000 
KAPPA=0.0
POTENTIALTYPE=2
PFLAG=0
SFLAG=1
DFLAG=0
BATCHSIZE=2000 
NUMDEVICES=1 
NUMTHREADS=1


MAXPARNODE=6000
BATCHSIZE=6000
for N in 300500
do 
	#SOURCES=/Users/nathanvaughn/Desktop/S$N.bin
	SOURCES=/scratch/krasny_fluxg/njvaughn/BenzeneRefinement/S$N.bin
	#TARGETS=/Users/nathanvaughn/Desktop/T$N.bin
	TARGETS=/scratch/krasny_fluxg/njvaughn/BenzeneRefinement/T$N.bin
	DIRECTSUM=/Users/nathanvaughn/Desktop/ex_coulombSS_$N.bin
	../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
done


#for N in 195500 300500 486000
#do 
#	SOURCES=/Users/nathanvaughn/Desktop/S$N.bin
#	#SOURCES=/scratch/krasny_fluxg/njvaughn/BenzeneRefinement/S$N.bin
#	TARGETS=/Users/nathanvaughn/Desktop/T$N.bin
#	#TARGETS=/scratch/krasny_fluxg/njvaughn/BenzeneRefinement/T$N.bin
#	DIRECTSUM=/Users/nathanvaughn/Desktop/ex_coulombSS_$N.bin
#	../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
#done


#MAXPARNODE=5000
#BATCHSIZE=5000
#for N in 689000 1455500 2824000 4773500 8294500
#do 
#	SOURCES=/Users/nathanvaughn/Desktop/S$N.bin
#	#SOURCES=/scratch/krasny_fluxg/njvaughn/BenzeneRefinement/S$N.bin
#	TARGETS=/Users/nathanvaughn/Desktop/T$N.bin
#	#TARGETS=/scratch/krasny_fluxg/njvaughn/BenzeneRefinement/T$N.bin
#	DIRECTSUM=/Users/nathanvaughn/Desktop/ex_coulombSS_$N.bin
#	../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
#done

