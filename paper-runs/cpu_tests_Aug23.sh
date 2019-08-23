TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



 
DS_CSV=/Users/nathanvaughn/Desktop/randomPoints/cpuTestingdirectSum_cpu_Coulomb.csv 


NUMDEVICES=0
NUMTHREADS=4


BATCHSIZE=4000
MAXPARNODE=4000


## COULOMB 
KAPPA=0.0
POTENTIALTYPE=0
OUTFILE=/Users/nathanvaughn/Desktop/randomPoints/cpuTesting/cpu_coulomb_lagrange.csv 
for N in 1000 10000
do
	echo N=$N 
	SOURCES=/Users/nathanvaughn/Desktop/randomPoints/S$N.bin    
	TARGETS=/Users/nathanvaughn/Desktop/randomPoints/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/Users/nathanvaughn/Desktop/randomPoints/cpuTesting/ex_st_coulomb_$N.bin
	
	direct-cpu $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	for THETA in 0.7
	do
		for ORDER in 8
		do
			tree-cpu   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER \
									$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
									$NUMDEVICES $NUMTHREADS
		done
	done
done   
