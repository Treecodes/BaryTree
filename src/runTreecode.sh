#!/bin/bash
#export OMP_NUM_THREADS=1

#../bin/tree.exe   ../examplesOxygenAtom/S21952.bin ../examplesOxygenAtom/T21952.bin ../examplesOxygenAtom/ex_st21952_coulomb.bin ../examplesOxygenAtom/out.tsv 21952 21952 0.5 5 1 500 0.0 0 1 0 0 50
#../bin/tree.exe   ../examplesOxygenAtom/S348488.bin ../examplesOxygenAtom/T348488.bin ../examplesOxygenAtom/ex_st348488_coulomb.bin ../examplesOxygenAtom/out.tsv 348488 348488 0.5 5 1 2000 0.0 0 1 0 0 50





echo multiple threads
export OMP_NUM_THREADS=6
../bin/tree.exe ../examplesOxygenAtom/S79576.bin ../examplesOxygenAtom/T79576.bin ../examplesOxygenAtom/ex_st79576_coulomb.bin ../examplesOxygenAtom/out.tsv 79576 79576 0.9 8 1 1000 0.0 0 1 0 0 500
#../bin/tree.exe   ../examplesOxygenAtom/S348488.bin ../examplesOxygenAtom/T348488.bin ../examplesOxygenAtom/ex_st348488_coulomb.bin ../examplesOxygenAtom/out.tsv 348488 348488 0.9 9 1 5000 0.0 0 1 0 0 500


#echo single thread
#export OMP_NUM_THREADS=1
#../bin_singleThread/tree.exe ../examplesOxygenAtom/S79576.bin ../examplesOxygenAtom/T79576.bin ../examplesOxygenAtom/ex_st79576_coulomb.bin ../examplesOxygenAtom/out.tsv 79576 79576 0.9 6 1 500 0.0 0 1 0 0 500
        
        
#../bin/tree.exe   ../examplesOxygenAtom/S348488.bin ../examplesOxygenAtom/T348488.bin ../examplesOxygenAtom/ex_st348488_coulomb.bin ../examplesOxygenAtom/out.tsv 348488 348488 0.9 8 1 20000 0.0 0 1 0 0 500
#../bin_noCheck/tree.exe   ../examplesOxygenAtom/S348488.bin ../examplesOxygenAtom/T348488.bin ../examplesOxygenAtom/ex_st348488_coulomb.bin ../examplesOxygenAtom/out.tsv 348488 348488 0.9 8 1 20000 0.0 0 1 0 0 500
