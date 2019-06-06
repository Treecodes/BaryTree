TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

#N=821000
#N=2365328
N=10000000
#1328096



## Coulomb: Hermite
KAPPA=0.0 
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/K20_hermite_GPU_parallelized_$N.csv
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_hermite_GPU_parallelized_$N.csv
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/bughunt.csv


SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
#SOURCES=/scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin
#TARGETS=/scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin  
 
#../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv $N $N 0.0 0 2

POTENTIALTYPE=4
for ORDER in {8..13}
do
	for THETA in 0.5
	  do    
		for BATCHSIZE in 5000
		  do       
		     for MAXPARNODE in 5000 
		     	do
		     	for NUMDEVICES in 2     
		     	do
		     		#echo Doing Nothing
		     		#../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv $N $N 0.0 0 $NUMDEVICES
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES
		     	done
		     done
		 done
	done
done
 


