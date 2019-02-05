####  PBS preamble

#PBS -N batchGreenIterations
#PBS -M njvaughn@umich.edu
#PBS -m a

#PBS -A krasny_fluxg
#PBS -l qos=flux
#PBS -q fluxg


#PBS -l nodes=1:gpus=1:titanv,mem=8gb
#PBS -l walltime=00:10:00
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


module purge
module load cuda cupti pgi openmpi/1.10.2/pgi/16.4 mkl/2017.0.1
nvidia-smi
export PGI_ACC_TIME=0


../bin_openACC/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_openACC_titanv.bin  out.tsv 21952  21952
../bin_openACC/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_openACC_titanv.bin  out.tsv 79576  79576
../bin_openACC/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_openACC_titanv.bin out.tsv 348488 348488
../bin_openACC/direct.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC_titanv.bin out.tsv 636608 636608


../bin_openACC/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T21952.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st21952_coulomb_openACC_titanv.bin  out.tsv 21952  21952  0.8 6 1 5000 0.0 0 1 0 0 500
../bin_openACC/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T79576.bin  /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st79576_coulomb_openACC_titanv.bin  out.tsv 79576  79576  0.8 6 1 5000 0.0 0 1 0 0 500
../bin_openACC/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T348488.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st348488_coulomb_openACC_titanv.bin out.tsv 348488 348488 0.8 6 1 5000 0.0 0 1 0 0 500
../bin_openACC/tree.exe   /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin /scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC_titanv.bin out.tsv 636608 636608 0.8 6 1 50000 0.0 0 1 0 0 3000
