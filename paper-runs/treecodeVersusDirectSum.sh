TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0


 

KAPPA=0.5
POTENTIALTYPE=1 

DSFILE=

OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/treecodeVersusDirectSum/yukawa.csv

SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S1328096.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T1328096.bin
NUMSOURCES=1328096
NUMTARGETS=1328096 
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st1328096_yukawa_titan.bin

../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S1328096.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T1328096.bin $DIRECTSUM /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds.tsv 1328096 1328096 0.5 1

for BATCHSIZE in 16000 
do
	for MAXPARNODE in 16000
	  do 
		for ORDER in 6 8 10 12
		  do 
		     for THETA in 0.5 0.6 0.7 0.8 0.9
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
		     done
		 done
	done
done