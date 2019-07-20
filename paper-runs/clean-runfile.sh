TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0


N=1000000



OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/batchSize_5k_$N.csv 


SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
 
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin  
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/directSum_TitanV_Coulomb.csv  


NUMDEVICES=0
NUMTHREADS=1
KAPPA=0.0
POTENTIALTYPE=0
../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES  

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
		     		## Profiling commands (commented out for now)
		     		#pgprof --cpu-profiling off --metrics flop_count_dp ../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		     		#../bin/direct.exe   $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		     		#pgprof --cpu-profiling off ../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     		#pgprof --cpu-profiling off  --metrics flop_count_dp ../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     		
		     		## Treecode run
		     		../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
		     	done
		     done 
		 done
	done
done 


