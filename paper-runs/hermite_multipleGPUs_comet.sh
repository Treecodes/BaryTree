TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

<<<<<<< HEAD
#N=821000
#N=2365328
N=1000000
#1328096
=======

N=1000000
>>>>>>> 85c431e862710acb5d23b4d83faa3205134dab38



## Coulomb: Hermite
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/K20_hermite_GPU_parallelized_$N.csv
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_hermite_GPU_parallelized_nonStatic_$N.csv
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/bughunt.csv

#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/TitanV_hermite_parallel_$N.csv 
<<<<<<< HEAD
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/comet/testing_$N.csv 
=======
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/batchSize_5k_$N.csv 
>>>>>>> 85c431e862710acb5d23b4d83faa3205134dab38


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
<<<<<<< HEAD
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/comet/directSum_Coulomb.csv 
=======
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/directSum_TitanV_Coulomb.csv  
>>>>>>> 85c431e862710acb5d23b4d83faa3205134dab38


NUMDEVICES=2
NUMTHREADS=2
KAPPA=0.0
<<<<<<< HEAD
POTENTIALTYPE=0
#../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES  
POTENTIALTYPE=0

export OMP_SCHEDULE=guided
for ORDER in 8
do
	for THETA in 0.9
	  do    
		for BATCHSIZE in 5000
		  do       
		     for MAXPARNODE in 5000
		     	do
		     	for NUMDEVICES in 1 2 3 4     
		     	do
			../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES
		     		#../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES
		     	done
		     done 
		 done
	done
done 
=======
POTENTIALTYPE=4
#../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES  

export PGI_ACC_TIME=0
export OMP_SCHEDULE=guided
for ORDER in 8  
do 
	for THETA in 0.7
	  do    
		for BATCHSIZE in 1000
		  do       
		     for MAXPARNODE in 1000
		     	do 
		     	for NUMDEVICES in 1 
		     	do
		     		 
		     		#pgprof --cpu-profiling off --metrics flop_count_dp ../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		     		#../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		     		pgprof --cpu-profiling off ../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     		pgprof --cpu-profiling off  --metrics flop_count_dp ../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     	done
		     done 
		 done
	done
done 


>>>>>>> 85c431e862710acb5d23b4d83faa3205134dab38
