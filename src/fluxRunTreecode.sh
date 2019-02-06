#!/bin/bash
#export OMP_NUM_THREADS=1


#nvidia-smi
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_openACC.bin out.tsv 348488 348488 0.8 8 1 5000 0.0 0 1 0 0 500

#module purge
#module load gcc/5.4.0 openmpi/3.0.0/gcc/5.4.0 mkl/2017.0.1
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb.bin out.tsv 21952 21952 0.8 6 1 5000 0.0 0 1 0 0 500
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb.bin out.tsv 79576 79576 0.8 6 1 5000 0.0 0 1 0 0 3000
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb.bin out.tsv 348488 348488 0.8 6 1 5000 0.0 0 1 0 0 500
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC.bin out.tsv 636608 636608 0.9 6 1 50000 0.0 0 1 0 0 3000



#module purge
#module load cuda cupti pgi openmpi/1.10.2/pgi/16.4 mkl/2017.0.1
#nvidia-smi
#export PGI_ACC_TIME=0

#../bin_openACC/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_openACC.bin  out.csv 21952  21952  0.8 6 1 5000 0.0 0 1 0 0 500
#../bin_openACC/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_openACC.bin  out.csv 79576  79576  0.8 6 1 5000 0.0 0 1 0 0 500
#../bin_openACC/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_openACC.bin out.csv 348488 348488 0.8 6 1 5000 0.0 0 1 0 0 500
#../bin_openACC/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC.bin out.csv 636608 636608 0.9 7 1 5000 0.0 0 1 0 0 2000



SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC.bin
OUTFILE=out636608.csv
NUMSOURCES=636608
NUMTARGETS=636608
TREETYPE=1
KAPPA=0.0
POTENTIALTYPE=0
SFLAG=1
PFLAG=0
DFLAG=0


BATCHSIZE=500
MAXPARNODE=5000
THETA=0.9
ORDER=7

for MAXPARNODE in 500 1000 2000 4000 
  do
	for ORDER in {5..12..1}
	  do 
	     for THETA in 0.5 0.6 0.7 0.8 0.9
	     	do
	     		echo MAXPARNODE $MAXPARNODE ORDER $ORDER THETA $THETA
	     		#../bin_openACC/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
	     done
	 done
done
