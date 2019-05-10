#!/bin/bash

#pgititan
#nvidia-smi
#export PGI_ACC_TIME=0




TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

OUTFILE=/home/njvaughn/synchronizedDataFiles/MICDE_Data_2019/gpu_treecode/tc5.csv

KAPPA=0.0
POTENTIALTYPE=0 


 
SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin
NUMSOURCES=21952
NUMTARGETS=21952
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_titan.bin

for BATCHSIZE in 1000 
do
	for MAXPARNODE in 1000
	  do
		for ORDER in 8
		  do 
		     for THETA in 0.8
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
		     done
		 done
	done
done


SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin
NUMSOURCES=79576
NUMTARGETS=79576
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_titan.bin

for BATCHSIZE in 1000 
do
	for MAXPARNODE in 1000
	  do
		for ORDER in 8
		  do 
		     for THETA in 0.8
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
		     done
		 done
	done
done

SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin
NUMSOURCES=348488
NUMTARGETS=348488
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_titan.bin

for BATCHSIZE in 2000 
do
	for MAXPARNODE in 4000
	  do
		for ORDER in 8
		  do 
		     for THETA in 0.8
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
		     done
		 done
	done
done


SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin
NUMSOURCES=636608
NUMTARGETS=636608
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_titan.bin
for BATCHSIZE in 8000 
do
	for MAXPARNODE in 8000
	  do
		for ORDER in 8
		  do 
		     for THETA in 0.8
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
		     done
		 done
	done
done





OUTFILE=/home/njvaughn/synchronizedDataFiles/MICDE_Data_2019/cpu_treecode/tc.csv


SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S3719492.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T3719492.bin
NUMSOURCES=3719492
NUMTARGETS=3719492
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st3719492_coulomb_titan.bin

for BATCHSIZE in 8000 
do
	for MAXPARNODE in 8000
	  do
		for ORDER in 8
		  do 
		     for THETA in 0.8
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
		     done
		 done
	done
done

 





 


