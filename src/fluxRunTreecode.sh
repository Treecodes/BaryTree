#!/bin/bash
#export OMP_NUM_THREADS=1

../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_openACC.bin  out.csv 348488  348488  0.9 12 1 4000 0.0 0 1 0 0 4000



#module purge
#module load gcc/5.4.0 openmpi/3.0.0/gcc/5.4.0 mkl/2017.0.1
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb.bin out.tsv 21952 21952 0.8 6 1 5000 0.0 0 1 0 0 500
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb.bin out.tsv 79576 79576 0.8 6 1 5000 0.0 0 1 0 0 3000
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb.bin out.tsv 348488 348488 0.8 6 1 5000 0.0 0 1 0 0 500
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC.bin out.tsv 636608 636608 0.9 6 1 50000 0.0 0 1 0 0 3000



#module purge
#module load cuda cupti pgi openmpi/1.10.2/pgi/16.4 mkl/2017.0.1
#nvidia-smi
export PGI_ACC_TIME=0

#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa_SS_0p5_openACC.bin  out.csv 21952  21952  0.7 6 1 500 0.5 3 1 0 0 500
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa_SS_0p5_openACC.bin  out.csv 79576  79576  0.7 6 1 4000 0.5 3 1 0 0 2000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_yukawa_SS_0p5_openACC.bin out.csv 348488 348488 0.7 6 1 5000 0.5 3 1 0 0 2000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa_SS_0p5_openACC.bin out.csv 636608 636608 0.7 6 1 9000 0.5 3 1 0 0 8000
  
  
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_SS_0p5_openACC.bin  out.csv 21952  21952  0.7 6 1 500 0.5 2 1 0 0 500
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_SS_0p5_openACC.bin  out.csv 79576  79576  0.7 6 1 4000 0.5 2 1 0 0 2000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_SS_0p5_openACC.bin out.csv 348488 348488 0.7 6 1 5000 0.5 2 1 0 0 2000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_SS_0p5_openACC.bin out.csv 636608 636608 0.7 6 1 9000 0.5 2 1 0 0 8000
  
  
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa_SS_0p5_openACC.bin  out.csv 79576  79576  0.9 6 1 2000 0.5 3 1 0 0 2000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa_SS_1p0_openACC.bin out.csv 636608 636608 0.7 6 1 8000 1.0 3 1 0 0 2000



#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa_SS_0p5_openACC.bin out.csv 636608 636608 0.7 6 1 8000 0.5 3 1 0 0 8000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa_SS_0p5_openACC.bin out.csv 636608 636608 0.7 7 1 8000 0.5 3 1 0 0 8000
  
   