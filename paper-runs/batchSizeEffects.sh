TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0




## Coulomb Batch Size Study
KAPPA=0.0
POTENTIALTYPE=0 
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/batchSize/coulomb.csv
 
SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S1328096.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T1328096.bin
NUMSOURCES=1328096
NUMTARGETS=1328096
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st1328096_coulomb_titan.bin

for BATCHSIZE in 1000 2000 4000 8000 16000 32000
do
	for MAXPARNODE in 2000 4000 16000
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