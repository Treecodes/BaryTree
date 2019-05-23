
SFLAG=1
PFLAG=0
DFLAG=0
TREETYPE=1

THETA=0.8
ORDER=8    
#pgprof --export-profile timeline.prof ../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 21952 21952 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000 
pgprof ../bin/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa.bin   /home/njvaughn/synchronizedDataFiles/KITCpaperData/treecode-single-convolution/ds_yukawa.tsv 21952 21952 $THETA $ORDER $TREETYPE 1000 0.5 1 $SFLAG $PFLAG $DFLAG 1000 
