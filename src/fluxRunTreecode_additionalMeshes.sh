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


 
#SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S1328096.bin
#TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T1328096.bin
#NUMSOURCES=1328096
#NUMTARGETS=1328096
#DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st1328096_coulomb_titan.bin
#
#for BATCHSIZE in 8000 
#do
#	for MAXPARNODE in 8000
#	  do
#		for ORDER in 8
#		  do 
#		     for THETA in 0.8
#		     	do
#		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
#		     done
#		 done
#	done
#done
  
  
SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S4727912.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T4727912.bin
NUMSOURCES=47279
NUMTARGETS=47279
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_titan.bin

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
  