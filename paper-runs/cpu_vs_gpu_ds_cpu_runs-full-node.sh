TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



 
#DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_directSum_GPU_parallelized.csv
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/directSum_cpu_Coulomb.csv 


NUMDEVICES=0
NUMTHREADS=6


BATCHSIZE=4000
MAXPARNODE=4000




## COULOMB 
#KAPPA=0.0
#POTENTIALTYPE=4
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/1M_comparison/cpu_coulomb_direct.csv 
#for N in 1000000
#do
#	echo N=$N 
#	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
#	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
#	NUMSOURCES=$N
#	NUMTARGETS=$N
#	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
#	
#	for THETA in 0
#	do
#		for ORDER in 0
#		do
#			direct-cpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
#		done
#	done
#done 


## Yukawa 
KAPPA=0.5
POTENTIALTYPE=5
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/1M_comparison/cpu_yukawa_direct.csv 
for N in 1000000
do
	echo N=$N 
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_yukawa_$N.bin
	
	for THETA in 0
	do
		for ORDER in 0
		do
			direct-cpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		done
	done
done 


