TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

#N=821000
#N=2365328
N=100000
#1328096
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st$N_coulomb.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv $N $N 0.0 0



## Coulomb: Hermite vs. Lagrange
KAPPA=0.0 
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/initialTests.csv
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/hermite_vs_lagrange_gpu_check4x_smallEnoughLeaf_$N.csv


SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
#SOURCES=/scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin
#TARGETS=/scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
#DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st$N_coulomb.bin     
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st$N_coulomb.bin 
 
#../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv $N $N 0.0 0

for BATCHSIZE in 1000 
do
	for MAXPARNODE in 1000
	  do 
		for ORDER in {4..14}
		  do 
		     for THETA in 0.7
		     	do
		     	for POTENTIALTYPE in 4   
		     	do
		     		#../bin_noACC/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done


