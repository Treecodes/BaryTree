TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

N=1000000
BATCHSIZE=5000
MAXPARNODE=5000
SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N

NUMDEVICES=4
NUMTHREADS=4

## COULOMB 
KAPPA=0.0
POTENTIALTYPE=0
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
#OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/treecodeVersusDirectSum/coulomb_hermite.csv 
OUTFILE=/home/njvaughn/testing/coulomb_hermite.csv 

for ORDER in 7
  do   
     for THETA in 0.8    
     	do
     	echo $THETA
     	#direct-gpu  $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
 		#tree-gpu   	$SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
 		#			$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
 		#			$NUMDEVICES $NUMTHREADS
     done
 done




## Yukawa 
KAPPA=0.5
POTENTIALTYPE=5 
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_yukawa_$N.bin
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/treecodeVersusDirectSum/yukawa_hermite.csv 
for ORDER in 5
  do   
     for THETA in 0.7    
     	do
     	direct-gpu  $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
 		tree-gpu   	$SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
 					$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
 					$NUMDEVICES $NUMTHREADS
     done
 done
 
 
