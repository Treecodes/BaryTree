TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



DS_CSV=/Users/nathanvaughn/Desktop/randomPoints/directSum_cpu_Coulomb.csv 


KAPPA=0.0
POTENTIALTYPE=0


ORDER=8
THETA=0.8
BATCHSIZE=1000
MAXPARNODE=1000

NUMDEVICES=0


OUTFILE=/Users/nathanvaughn/Desktop/randomPoints/cpu_Coulomb.csv 
for N in 64000
do
	echo N=$N
	SOURCES=/Users/nathanvaughn/Desktop/randomPoints/S$N.bin    
	TARGETS=/Users/nathanvaughn/Desktop/randomPoints/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/Users/nathanvaughn/Desktop/randomPoints/ex_st_coulomb_$N.bin
	for np in 4
	do
		for NUMTHREADS in 1
		do
			export OMP_NUM_THREADS=$NUMTHREADS
			mpirun -np $np tree-distributed-cpu $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $THETA $ORDER \
							 					$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
							 					$NUMDEVICES $NUMTHREADS
			#mpirun -np $np direct-distributed-cpu $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		done
	done 
done
