TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



DS_CSV=/home/njvaughn/synchronizedDataFiles/profiling/ds.csv


NUMDEVICES=1
KAPPA=0.0
POTENTIALTYPE=0


ORDER=3
THETA=0.8
BATCHSIZE=10
MAXPARNODE=100

NUMDEVICES=0
NUMTHREADS=1


OUTFILE=/Users/nathanvaughn/Desktop/randomPoints/cpu_Coulomb.csv 
for N in 100000
do
	echo N=$N
	SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
	TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin  
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=home/njvaughn/synchronizedDataFiles/profiling/_ds_$N.csv
	mpirun -np 8 ../bin/direct-distributed.exe $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	#mpirun -np 4 ../bin/direct-distributed.exe $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	#mpirun -np 2 ../bin/direct-distributed.exe $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	#mpirun -np 1 ../bin/direct-distributed.exe $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	#../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
done 

