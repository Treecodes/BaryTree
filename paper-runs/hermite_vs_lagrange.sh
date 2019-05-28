TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

#N=821000
N=2365328
#1328096
../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st$N_coulomb.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv $N $N 0.0 0


## Coulomb: Hermite vs. Lagrange
KAPPA=0.0
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/initialTests.csv
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/compareToLagrange.csv


SOURCES=/scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st$N_coulomb.bin

for BATCHSIZE in 10000 
do
	for MAXPARNODE in 10000
	  do
		for ORDER in 4 6 8 
		  do 
		     for THETA in 0.5 0.7 0.9
		     	do
		     	for POTENTIALTYPE in 0 4
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done

for BATCHSIZE in 10000 
do
	for MAXPARNODE in 10000
	  do
		for ORDER in 3 5  
		  do 
		     for THETA in 0.5 0.7 0.9
		     	do
		     	for POTENTIALTYPE in 4
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done


for BATCHSIZE in 10000 
do
	for MAXPARNODE in 10000
	  do
		for ORDER in 7 9  
		  do 
		     for THETA in 0.5 0.7 0.9
		     	do
		     	for POTENTIALTYPE in 0
		     	do
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done
