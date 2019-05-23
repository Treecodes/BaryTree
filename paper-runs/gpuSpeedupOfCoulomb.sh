#!/bin/bash

## GPU Direct Sum Runs
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S197000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T197000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st197000_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv 197000 197000 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S821000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T821000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st821000_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv 821000 821000 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S1486000.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/T1486000.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st1486000_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv 1486000 1486000 0.5 1
#../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S2365328.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/T2365328.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st2365328_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv 2365328 2365328 0.5 1
../bin/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S2795000.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/T2795000.bin /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st2795000_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds.csv 2795000 2795000 0.5 1
 
 
 

SFLAG=1
PFLAG=0
DFLAG=0
TREETYPE=1

THETA=0.8
ORDER=8     
## CPU Treecode Runs

../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S197000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T197000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st197000_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 197000 197000 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000 
../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S821000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T821000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st821000_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 821000 821000 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000
../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S1486000.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/T1486000.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st1486000_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 1486000 1486000 $THETA $ORDER $TREETYPE 4000 0.5 1 $SFLAG $PFLAG $DFLAG 4000
../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S2365328.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/T2365328.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st2365328_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 2365328 2365328 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
../bin_noAcc/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S2795000.bin /scratch/krasny_fluxg/njvaughn/examplesBenzene/T2795000.bin /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st2795000_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 2795000 2795000 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
 

## GPU Treecode Runs

#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S197000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T197000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st197000_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 197000 197000 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000 
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S821000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/T821000.bin   /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st821000_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 821000 821000 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S1486000.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/T1486000.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st1486000_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 1486000 1486000 $THETA $ORDER $TREETYPE 4000 0.5 1 $SFLAG $PFLAG $DFLAG 4000
#../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S2365328.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/T2365328.bin  /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st2365328_yukawa.bin  /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 2365328 2365328 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesBenzene/S2795000.bin /scratch/krasny_fluxg/njvaughn/examplesBenzene/T2795000.bin /scratch/krasny_fluxg/njvaughn/examplesBenzene/ex_st2795000_yukawa.bin /home/njvaughn/synchronizedDataFiles/KITCpaperData/benzeneData/coulombSpeedup/ds_yukawa.csv 2795000 2795000 $THETA $ORDER $TREETYPE 8000 0.5 1 $SFLAG $PFLAG $DFLAG 8000
 


