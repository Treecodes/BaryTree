TREETYPE=1

SFLAG=1
PFLAG=0
DFLAG=0

#N=821000
#N=2365328
N=1000000
#1328096



## Coulomb: Hermite
KAPPA=0.0 
OUTFILE=/home/njvaughn/synchronizedDataFiles/KITCpaperData/hermiteTesting/coulomb/K40_hermite_batch_node_size_$N.csv


SOURCES=/scratch/krasny_fluxg/njvaughn/random/S$N.bin    
TARGETS=/scratch/krasny_fluxg/njvaughn/random/T$N.bin
#SOURCES=/scratch/krasny_fluxg/njvaughn/examplesBenzene/S$N.bin
#TARGETS=/scratch/krasny_fluxg/njvaughn/examplesBenzene/T$N.bin
NUMSOURCES=$N
NUMTARGETS=$N
#DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st$N_coulomb.bin     
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/random/ex_st$N_coulomb.bin  
 
#../bin35/direct.exe   $SOURCES $TARGETS $DIRECTSUM   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv $N $N 0.0 0

for ORDER in 7 
do
	for THETA in 0.7
	  do    
		for BATCHSIZE in {1600..2400..100}
		  do       
		     for MAXPARNODE in {1600..2400..100} 
		     	do
		     	for POTENTIALTYPE in 4     
		     	do
		     		#../bin_noACC/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     		../bin35/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done
 

for ORDER in 7  
do
	for THETA in 0.7
	  do    
		for BATCHSIZE in {1000..16000..500}
		  do       
		     for MAXPARNODE in {1000..16000..500} 
		     	do
		     	for POTENTIALTYPE in 4     
		     	do
		     		#../bin_noACC/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     		../bin35/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done



for ORDER in 7 
do
	for THETA in 0.7
	  do    
		for BATCHSIZE in {1600..2400..100}
		  do       
		     for MAXPARNODE in {1000..16000..500}   
		     	do
		     	for POTENTIALTYPE in 4     
		     	do
		     		#../bin_noACC/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     		../bin35/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done


for ORDER in 7 
do
	for THETA in 0.7
	  do    
		for BATCHSIZE in {1000..16000..500}
		  do       
		     for MAXPARNODE in {1600..2400..100} 
		     	do
		     	for POTENTIALTYPE in 4     
		     	do
		     		#../bin_noACC/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     		../bin35/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done

for ORDER in 7 
do
	for THETA in 0.7
	  do    
		for BATCHSIZE in {1600..2400..100}
		  do       
		     for MAXPARNODE in {1600..2400..100} 
		     	do
		     	for POTENTIALTYPE in 4     
		     	do
		     		#../bin_noACC/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     		../bin35/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE
		     	done
		     done
		 done
	done
done


