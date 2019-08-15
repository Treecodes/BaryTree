TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



DS_CSV=/scratch/krasny_fluxg/njvaughn/random/directSum_cpu_Coulomb.csv 


KAPPA=0.0
POTENTIALTYPE=0


ORDER=4
THETA=0.8
BATCHSIZE=1
MAXPARNODE=5

NUMDEVICES=0
NUMTHREADS=1
export OMP_NUM_THREADS=$NUMTHREADS

OUTFILE=/scratch/krasny_fluxg/njvaughn/random/cpu_Coulomb.csv 
for N in 64
do
	echo N=$N
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
	for np in 2
	do
			#mpirun -np $np direct-distributed-cpu $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
			mpirun -np $np tree-distributed-cpu $SOURCES $TARGETS $DIRECTSUM $OUTFILE $N $N $THETA $ORDER \
							 					$TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE \
							 					$NUMDEVICES $NUMTHREADS
	done 
done


