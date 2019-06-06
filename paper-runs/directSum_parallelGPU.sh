TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

N=10000000
#1328096



## Coulomb: Hermite
KAPPA=0.0 
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_directSum_GPU_parallelized.csv


SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
#SOURCES=/scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin
#TARGETS=/scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N

for N in 100000 1000000 10000000
do
	for NUMDEVICES in 1 2
	do
		DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
		../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N 0.0 0 $NUMDEVICES
	done  
done