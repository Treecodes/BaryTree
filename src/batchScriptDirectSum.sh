####  PBS preamble

#PBS -N treecodeTesting
#PBS -M njvaughn@umich.edu
#PBS -m a

#PBS -A krasny_fluxg
#PBS -l qos=flux
#PBS -q fluxg


#PBS -l nodes=1:gpus=1:titanv,mem=8gb
#PBS -l walltime=01:00:00
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


pgititan

../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_titan.bin out.tsv 21952 21952 0.0 0
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_titan.bin out.tsv 79576 79576 0.0 0
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_titan.bin out.tsv 348488 348488 0.0 0
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_titan.bin out.tsv 636608 636608 0.0 0

../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa_0p5_titan.bin out.tsv 21952 21952 0.5 1
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa_0p5_titan.bin out.tsv 79576 79576 0.5 1
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_yukawa_0p5_titan.bin out.tsv 348488 348488 0.5 1
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa_0p5_titan.bin out.tsv 636608 636608 0.5 1

../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_SS_0p5_titan.bin out.tsv 21952 21952 0.5 2
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_SS_0p5_titan.bin out.tsv 79576 79576 0.5 2
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_SS_0p5_titan.bin out.tsv 348488 348488 0.5 2
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_SS_0p5_titan.bin out.tsv 636608 636608 0.5 2

../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_yukawa_SS_0p5_titan.bin out.tsv 21952 21952 0.5 3
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_yukawa_SS_0p5_titan.bin out.tsv 79576 79576 0.5 3
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_yukawa_SS_0p5_titan.bin out.tsv 348488 348488 0.5 3
../binTitan/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_yukawa_SS_0p5_titan.bin out.tsv 636608 636608 0.5 3

nvidia-smi



