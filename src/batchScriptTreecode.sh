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


module purge
module load cuda cupti pgi openmpi/1.10.2/pgi/16.4 mkl/2017.0.1
nvidia-smi
export PGI_ACC_TIME=0


SOURCES=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/S636608.bin
TARGETS=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/T636608.bin
DIRECTSUM=/scratch/krasny_fluxg/njvaughn/examplesOxygenAtom/ex_st636608_coulomb_openACC.bin
OUTFILE=out636608.csv
NUMSOURCES=636608
NUMTARGETS=636608
TREETYPE=1
KAPPA=0.0
POTENTIALTYPE=0
SFLAG=1
PFLAG=0
DFLAG=0


BATCHSIZE=500
MAXPARNODE=5000
#THETA=0.9
#ORDER=7

for BATCHSIZE in 2000 4000
do
	for MAXPARNODE in 1000 2000 4000 8000 16000
	  do
		for ORDER in {5..10..1}
		  do 
		     for THETA in 0.6 0.7 0.8
		     	do
		     		../bin_openACC/tree.exe   $SOURCES $TARGETS $DIRECTSUM $OUTFILE $NUMSOURCES $NUMTARGETS $THETA $ORDER $TREETYPE $MAXPARNODE $KAPPA $POTENTIALTYPE $SFLAG $PFLAG $DFLAG $BATCHSIZE
		     done
		 done
	done
done
