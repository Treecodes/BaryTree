!/bin/bash
#export OMP_NUM_THREADS=1

#../bin/tree.exe   ../examplesOxygenAtom/S21952.bin ../examplesOxygenAtom/T21952.bin ../examplesOxygenAtom/ex_st21952_coulomb.bin ../examplesOxygenAtom/out.tsv 21952 21952 0.7 10 1 500 0.0 0 1 0 0 50
#../bin/tree.exe   ../examplesOxygenAtom/S79576.bin ../examplesOxygenAtom/T79576.bin ../examplesOxygenAtom/ex_st79576_coulomb.bin ../examplesOxygenAtom/out.tsv 79576 79576 0.9 8 1 2000 0.0 0 1 0 0 500

#../bin/tree.exe   ../examplesOxygenAtom/S79576.bin ../examplesOxygenAtom/T79576.bin ../examplesOxygenAtom/ex_st79576_coulomb.bin ../examplesOxygenAtom/out.tsv 79576 79576 0.7 7 1 1000 0.0 0 1 0 0 1000


SOURCES=../examplesOxygenAtom/S21952.bin
TARGETS=../examplesOxygenAtom/T21952.bin
DIRECTSUM=../examplesOxygenAtom/ex_st21952_coulomb.bin
OUTFILE=../examplesOxygenAtom/out.tsv
NUMSOURCES=21952
NUMTARGETS=21952
THETA=0.9
ORDER=8
TREETYPE=1
MAXPARNODE=1000
KAPPA=0.0
POTENTIALTYPE=0
PFLAG=0
SFLAG=1
DFLAG=0
BATCHSIZE=1000
NUMDEVICES=0

#export TMPDIR=/tmp
export OMP_NUM_THREADS=4
../bin/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $PFLAG $SFLAG $DFLAG $BATCHSIZE $NUMDEVICES
