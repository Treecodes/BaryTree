TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0



DS_CSV=/Users/nathanvaughn/Desktop/randomPoints/directSum_cpu_Coulomb.csv 


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
for N in 10000
do
	echo N=$N
	SOURCES=/Users/nathanvaughn/Desktop/randomPoints/S$N.bin    
	TARGETS=/Users/nathanvaughn/Desktop/randomPoints/T$N.bin
	NUMSOURCES=$N
	NUMTARGETS=$N
	DIRECTSUM=/Users/nathanvaughn/Desktop/randomPoints/ex_st_coulomb_$N.bin
	mpirun -np 3 ../bin/direct-distributed.exe $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	mpirun -np 2 ../bin/direct-distributed.exe $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	mpirun -np 1 ../bin/direct-distributed.exe $SOURCES $TARGETS $DIRECTSUM $DS_CSV $N $N $KAPPA $POTENTIALTYPE $NUMDEVICES $NUMTHREADS
	#../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES $NUMTHREADS
done 

