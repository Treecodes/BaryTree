####  PBS preamble

#PBS -N treecodeTesting
#PBS -M njvaughn@umich.edu
#PBS -m a

#PBS -A krasny_fluxg
#PBS -l qos=flux
#PBS -q fluxg


#PBS -l nodes=1:gpus=1,mem=8gb
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -V

####  End PBS preamble

if [ -s "$PBS_NODEFILE" ] ; then
    echo "Running on"
    uniq -c $PBS_NODEFILE
fi

if [ -d "$PBS_O_WORKDIR" ] ; then
    cd $PBS_O_WORKDIR
    echo "Running from $PBS_O_WORKDIR"
fi

cd /home/njvaughn/hybrid-gpu-treecode/src


gccmods



../bin_noAcc/direct_with03.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_noAcc_withO3.bin out.tsv 79576 79576 0.0 0
../bin_noAcc/direct_with03.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa0p5_noAcc_withO3.bin out.tsv 79576 79576 0.5 1


../bin_noAcc/direct_without03.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_noAcc.bin out.tsv 79576 79576 0.0 0
../bin_noAcc/direct_without03.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa0p5_noAcc.bin out.tsv 79576 79576 0.5 1
