TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

#1328096
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S197000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T197000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st197000_coulomb.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv 197000 197000 0.0 0


## Coulomb: Hermite vs. Lagrange
KAPPA=0.0
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/initialTests.csv
 
SOURCES=/scratch/krasny_fluxg/njvaughn/examplesBenzene/S197000.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesBenzene/T197000.bin
NUMSOURCES=197000
NUMTARGETS=197000
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st197000_coulomb.bin

for BATCHSIZE in 2000 
do
	for MAXPARNODE in 2000
	  do
		for ORDER in 5
		  do 
		     for THETA in 0.7
		     	do
		     	for POTENTIALTYPE in 4
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done