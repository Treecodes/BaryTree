TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0




## Coulomb: Hermite vs. Lagrange
KAPPA=0.0
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/initialTests.csv
 
SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S1328096.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T1328096.bin
NUMSOURCES=1328096
NUMTARGETS=1328096
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st1328096_coulomb_titan.bin

for BATCHSIZE in 8000 
do
	for MAXPARNODE in 8000
	  do
		for ORDER in 6
		  do 
		     for THETA in 0.7
		     	do
		     	for POTENTIALTYPE in 0 4
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done