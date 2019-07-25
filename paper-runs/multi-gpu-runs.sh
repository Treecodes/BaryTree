TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

N=100000

BATCHSIZE=5000
MAXPARNODE=5000  
 

SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin      
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulombSS_$N.bin  
DS_CSV=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/TitanV_directSum_GPU_parallelized.csv


NUMDEVICES=1
NUMTHREADS=1 

   
## COULOMB   
KAPPA=0.0
POTENTIALTYPE=0
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/testing/coulomb_lagrange.csv 
for ORDER in {1..14}
  do   
     for THETA in 0.4 0.6 0.8 0.9    
     	do
     	echo $THETA
 		tree-gpu   	$SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
 					$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
 					$NUMDEVICES $NUMTHREADS
     done
 done




## Yukawa 
KAPPA=0.5
POTENTIALTYPE=1 
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_yukawa_$N.bin
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/parallelGPU/testing/yukawa_lagrange.csv 
for ORDER in {1..14}
  do   
     for THETA in 0.4 0.6 0.8 0.9    
     	do
     	
 		tree-gpu   	$SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
 					$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
 					$NUMDEVICES $NUMTHREADS
     done
 done


