#!/bin/bash

#../bin_noAcc/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb.bin /home/njvaughn/MICDE_Data_2019/cpu_treecode/ds.tsv 21952 21952 0.0 0
#../bin_noAcc/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb.bin /home/njvaughn/MICDE_Data_2019/cpu_treecode/ds.tsv 79576 79576 0.0 0
#../bin_noAcc/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb.bin /home/njvaughn/MICDE_Data_2019/cpu_treecode/ds.tsv 348488 348488 0.0 0
#../bin_noAcc/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb.bin /home/njvaughn/MICDE_Data_2019/cpu_treecode/ds.tsv 636608 636608 0.0 0
 
 
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC.bin out.tsv 636608 636608 0.0 0


../bin_noAcc/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S3719492.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T3719492.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st3719492_coulomb_titan.bin out.tsv 3719492 3719492 0.0 0
 
  
#nvidia-smi
#export PGI_ACC_TIME=0

#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_SS_0p5_openACC.bin out.tsv 21952 21952 0.5 2
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_SS_0p5_openACC.bin out.tsv 79576 79576 0.5 2
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_SS_0p5_openACC.bin out.tsv 348488 348488 0.5 2
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_SS_0p5_openACC.bin out.tsv 636608 636608 0.5 2



#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa_SS_0p5_openACC.bin out.tsv 21952 21952 0.5 3
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa_SS_0p5_openACC.bin out.tsv 79576 79576 0.5 3
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_yukawa_SS_0p5_openACC.bin out.tsv 348488 348488 0.5 3
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa_SS_0p5_openACC.bin out.tsv 636608 636608 0.5 3


#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_openACC.bin out.tsv 21952 21952 0.0 0
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_openACC.bin out.tsv 79576 79576 0.0 0
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_openACC.bin out.tsv 348488 348488 0.0 0
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC.bin out.tsv 636608 636608 0.0 0
 




#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_SS_0p5_openACC.bin out.tsv 79576 79576 4.0 2
