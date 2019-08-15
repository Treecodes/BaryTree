TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



DS_CSV=/scratch/krasny_fluxg/njvaughn/random/directSum_cpu_Coulomb.csv 


KAPPA=0.0
POTENTIALTYPE=0


ORDER=5 
THETA=0.8
BATCHSIZE=500 
MAXPARNODE=500

NUMDEVICES=0
NUMTHREADS=1
export OMP_NUM_THREADS=$NUMTHREADS

OUTFILE=/scratch/krasny_fluxg/njvaughn/random/cpu_Coulomb.csv 
for N in 64000
do
	echo N=$N
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
	for np in 4 8
	do
			#mpirun -np $np direct-distributed-cpu $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
			#mpirun -np $np tree-distributed-cpu $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $THETA $ORDER \
			mpirun -np $np tree-distributed-cpu $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $THETA $ORDER \
							 					$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
							 					$NUMDEVICES $NUMTHREADS
	done 
done


