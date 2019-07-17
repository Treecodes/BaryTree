TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

N=10000000



## Coulomb: Hermite
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/K20_hermite_GPU_parallelized_$N.csv
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_hermite_GPU_parallelized_nonStatic_$N.csv
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/bughunt.csv

#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/TitanV_hermite_parallel_$N.csv 
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/batchSize_5k_$N.csv 


SOURCES=/oasis/scratch/comet/njvaughn/temp_project/random/S$N.bin    
TARGETS=/oasis/scratch/comet/njvaughn/temp_project/random/T$N.bin
#SOURCES=/scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin
#TARGETS=/scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
#DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulombSS_$N.bin  
#DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_yukawa_$N.bin  
DIRECTSUM=/oasis/scratch/comet/njvaughn/temp_project/random/ex_st_coulomb_$N.bin  
#DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_directSum_GPU_parallelized.csv
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/comet/directSum_Coulomb.csv 


KAPPA=0.0
POTENTIALTYPE=0

for ORDER in 7
do
	for THETA in 0.7
	  do    
		for BATCHSIZE in 5000
		  do       
		     for MAXPARNODE in 5000
		     	do
		     	for NUMDEVICES in 4 3 2 1     
		     		do
				NUMTHREADS=$NUMDEVICES
				#../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     		done
		     	done 
		 done
	done
done 
