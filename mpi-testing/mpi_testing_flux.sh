TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



DS_CSV=/home/njvaughn/synchronizedDataFiles/profiling/ds.csv


KAPPA=0.0
POTENTIALTYPE=0


ORDER=3
THETA=0.8
BATCHSIZE=10 
MAXPARNODE=100  
 
NUMDEVICES=0 


OUTFILE=/home/njvaughn/mpi-testing/test.csv 
for N in 10000
do
	echo N=$N
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N 
	DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st_coulomb_$N.bin
	for np in 24 48
	do
		for NUMTHREADS in 1
		do
			export OMP_NUM_THREADS=$NUMTHREADS
			#export OMP_DISPLAY_ENV=true
			#export OMP_PROC_BIND=spread
			#export OMP_PLACES=cores
			#export OMP_NESTED=true	
			mpirun -np $np --bind-to none direct-distributed-cpu $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
		done
	done 
done

