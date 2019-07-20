TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

N=10000000
BATCHSIZE=5000
MAXPARNODE=5000
SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N

NUMDEVICES=1
NUMTHREADS=1

## COULOMB 
KAPPA=0.0
POTENTIALTYPE=4
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/treecodeVersusDirectSum/coulomb_hermite.csv 
for ORDER in {1..14}
  do   
     for THETA in 0.4 0.6 0.8 0.9    
     	do
     	echo $THETA
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
for ORDER in {1..14}
  do   
     for THETA in 0.4 0.6 0.8 0.9    
     	do
     	
 		tree-gpu   	$SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
 					$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
 					$NUMDEVICES $NUMTHREADS
     done
 done
 
 
