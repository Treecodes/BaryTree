TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



 
#DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_directSum_GPU_parallelized.csv
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/directSum_gpu_Coulomb.csv 


NUMDEVICES=1
NUMTHREADS=1


BATCHSIZE=4000
MAXPARNODE=4000


## COULOMB 
#KAPPA=0.0
#POTENTIALTYPE=0
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/1M_comparison/gpu_coulomb_lagrange.csv 
#for N in 1000000
#do
#	echo N=$N 
#	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
#	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
#	NUMSOURCES=$N
#	NUMTARGETS=$N
#	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
#	
#	for THETA in 0.3 0.5 0.7 0.9
#	do
#		for ORDER in {1..14}
#		do
#			tree-gpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
#									$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
#									$NUMDEVICES $NUMTHREADS
#		done
#	done
#done 


## Yukawa 
KAPPA=0.5
POTENTIALTYPE=1
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/1M_comparison/gpu_yukawa_lagrange.csv 
for N in 1000000
do
	echo N=$N 
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_yukawa_$N.bin
	
	direct-gpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	
	for THETA in 0.3 0.5 0.7 0.9
	do
		for ORDER in {1..14}
		do
			#tree-gpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
			#						$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
									$NUMDEVICES $NUMTHREADS
		done
	done
done 


## COULOMB 
KAPPA=0.0
POTENTIALTYPE=4
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/1M_comparison/gpu_coulomb_hermite.csv 
for N in 1000000
do
	echo N=$N 
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
	
	for THETA in 0.3 0.5 0.7 0.9
	do
		for ORDER in {1..14}
		do
			tree-gpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
									$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
									$NUMDEVICES $NUMTHREADS
		done
	done
done 


## Yukawa 
KAPPA=0.5
POTENTIALTYPE=5
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/gpu_vs_cpu/1M_comparison/gpu_yukawa_hermite.csv 
for N in 1000000
do
	echo N=$N 
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_yukawa_$N.bin
	
	for THETA in 0.3 0.5 0.7 0.9
	do
		for ORDER in {1..14}
		do
			tree-gpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
									$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
									$NUMDEVICES $NUMTHREADS
		done
	done
done 


