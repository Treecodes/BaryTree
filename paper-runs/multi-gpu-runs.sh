TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

N=100000000

BATCHSIZE=10000
MAXPARNODE=10000  
 

SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin      
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulombSS_$N.bin  
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_directSum_GPU_parallelized.csv


ORDER=7
THETA=0.7 

   
## COULOMB   
KAPPA=0.0
POTENTIALTYPE=0
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/testing/coulomb_lagrange.csv 
for NUMDEVICES in 4 2 1
do   
  	NUMTHREADS=$NUMDEVICES
 	echo $NUMDEVICES
	tree-gpu   	$SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
				$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
				$NUMDEVICES $NUMTHREADS
 done






