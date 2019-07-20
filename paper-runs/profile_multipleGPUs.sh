TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

N=100000
#1328096


 
## Coulomb: Hermite
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/K20_hermite_GPU_parallelized_$N.csv
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_hermite_GPU_parallelized_nonStatic_$N.csv
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/bughunt.csv

OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/TitanV_Lagrange_yukawaSS_$N.csv 
TCPROF=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/TC2_theta0p9b.prof
DSPROF=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/DS.prof

SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
#SOURCES=/scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin
#TARGETS=/scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulombSS_$N.bin  
#DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_yukawa_$N.bin  
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_directSum_GPU_parallelized.csv



POTENTIALTYPE=1
KAPPA=0.5
#pgprof --export-profile $DSPROF ../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES  
ORDER=8

for POTENTIALTYPE in 0
do
	for THETA in 0.8
	  do    
		for BATCHSIZE in 5000
		  do       
		     for MAXPARNODE in 5000 
		     	do
		     	for NUMDEVICES in 0
		     	do
					NUMTHREADS=$NUMDEVICES 
					NUMTHREADS=4
		     		TCPROF=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/TC_nvprof_MULTISTREAM_$NUMDEVICES.nvvp
		     		#../build/bin/tree-cpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     		#direct-gpu   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		     		tree-cpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     		#pgprof --export-profile $TCPROF -f tree-gpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS

		     		#nvprof --output-profile $TCPROF -f tree-gpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     	done
		     done
		 done
	done
done
 


