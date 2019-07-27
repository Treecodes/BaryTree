TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



 
#DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_directSum_GPU_parallelized.csv
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/directSum_cpu_Coulomb.csv 


NUMDEVICES=0
NUMTHREADS=1



ORDER=8
THETA=0.8
BATCHSIZE=4000  
MAXPARNODE=4000


## COULOMB 
KAPPA=0.0
POTENTIALTYPE=0
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/cpu_Coulomb.csv 
for N in 100000 1000000 10000000 
do
	echo N=$N 
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
	tree-cpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
							$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
							$NUMDEVICES $NUMTHREADS
done  


## Yukawa 
KAPPA=0.5
POTENTIALTYPE=1
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/cpu_Yukawa.csv 
for N in 100000 1000000 10000000
do
	echo N=$N 
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_yukawa_$N.bin
	tree-cpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
							$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
							$NUMDEVICES $NUMTHREADS
done 



