#!/bin/bash
#export OMP_NUM_THREADS=1



#
nvidia-smi
pgprof ../bin/tree.exe  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb.bin out.tsv 79576 79576 0.8 6 1 5000 0.0 0 1 0 0 500

