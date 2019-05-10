#!/bin/bash

## GPU Direct Sum Runs
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds.tsv 21952 21952 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds.tsv 79576 79576 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds.tsv 348488 348488 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds.tsv 636608 636608 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S1328096.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T1328096.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st1328096_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds.tsv 1328096 1328096 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S3719492.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T3719492.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st3719492_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds.tsv 3719492 3719492 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S4727912.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T4727912.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st4727912_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds.tsv 4727912 4727912 0.5 1
 
 
 

SFLAG=1
PFLAG=0
DFLAG=0
TREETYPE=1

THETA=0.8
ORDER=8     
## CPU Treecode Runs

#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 21952 21952 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000 
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 79576 79576 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 348488 348488 $THETA $ORDER $TREETYPE 4000 0.5 1 $SFLAG $PFLAG $DFLAG 4000
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 636608 636608 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S1328096.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T1328096.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st1328096_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 1328096 1328096 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S3719492.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T3719492.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st3719492_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 3719492 3719492 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
#../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S4727912.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T4727912.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st4727912_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 4727912 4727912 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
 

## GPU Treecode Runs

#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 21952 21952 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000 
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 79576 79576 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 348488 348488 $THETA $ORDER $TREETYPE 4000 0.5 1 $SFLAG $PFLAG $DFLAG 4000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 636608 636608 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S1328096.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T1328096.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st1328096_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 1328096 1328096 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S3719492.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T3719492.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st3719492_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 3719492 3719492 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S4727912.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T4727912.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st4727912_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 4727912 4727912 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
 


