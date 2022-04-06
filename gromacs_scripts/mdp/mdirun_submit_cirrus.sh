#!/bin/bash --login

#SBATCH --job-name=[job name]
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err
#SBATCH --partition=gpu-cascade
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=20:00:00
#SBATCH --exclusive
#SBATCH --account=[user account]

# Input files basename
BASE_NAME="1fin-WAT"

# Directories
EM='em'
NVT='nvt'
NPT='npt'
MD='md'

# Load GROMACS and MPI modules
purge_modules () {
    module purge
    module load /scratch/sw/modulefiles/epcc/setup-env
}

gmx_modules () {
    purge_modules
    module load openmpi/4.1.2
    module load cmake/3.17.3
    module load intel-19.5/cmkl
}

gro_gpu_modules () {
    purge_modules
    module load nvidia/nvhpc-nompi/22.2
    module load openmpi/4.1.2-cuda-11.6
    module load cmake/3.17.3
    module load intel-19.5/cmkl
}

GRO_PATH='/work/dc010/dc010/shared/sw/gromacs/2021.5-gpu/'
GMX_CPU_PATH='/work/dc010/dc010/ricci/SOFTWARE/GROMACS/bin/'
export PATH=$PATH:$GRO_PATH/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GRO_PATH/lib64/

export OMP_PLACES=cores
export OMP_NUM_THREADS=10

# * * * * * * * * *
# Minimization
# * * * * * * * * *
mkdir -p $EM
# - Use CPU version
gmx_modules
mpirun  -n 1 \
$GMX_CPU_PATH\gmx_mpi grompp \
    -f minim.mdp \
    -c ${BASE_NAME}.gro \
    -p ${BASE_NAME}.top \
    -po ${EM}/${BASE_NAME}_mdout.mdp \
    -o ${EM}/${EM}.tpr
# - Use GPU version
gro_gpu_modules
srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS \
        mdrun_mpi \
        -ntomp $OMP_NUM_THREADS -nb gpu -bonded cpu \
        -deffnm ${EM}/${EM} \
        -s ${EM}/${EM}.tpr

# * * * * * * * * *
# NVT equilibration
# * * * * * * * * *
mkdir -p $NVT
# - Use CPU version
gmx_modules
mpirun  -n 1 \
$GMX_CPU_PATH\gmx_mpi grompp \
    -f nvt.mdp \
    -c ${EM}/${EM}.gro \
    -r ${EM}/${EM}.gro \
    -p ${BASE_NAME}.top \
    -po ${NVT}/${BASE_NAME}_mdout.mdp \
    -o  ${NVT}/${NVT}.tpr
# - Use GPU version
gro_gpu_modules
srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS \
        mdrun_mpi \
        -ntomp $OMP_NUM_THREADS -nb gpu -bonded cpu \
        -deffnm ${NVT}/${NVT} \
        -s ${NVT}/${NVT}.tpr

# * * * * * * * * *
# NPT equilibration
# * * * * * * * * *
mkdir -p $NPT
# - Use CPU version
gmx_modules
mpirun  -n 1 \
$GMX_CPU_PATH\gmx_mpi grompp \
    -f npt.mdp \
    -c ${NVT}/${NVT}.gro \
    -r ${NVT}/${NVT}.gro \
    -t ${NVT}/${NVT}.cpt \
    -p ${BASE_NAME}.top \
    -po ${NPT}/${BASE_NAME}_mdout.mdp \
    -o  ${NPT}/npt.tpr
# - Use GPU version
gro_gpu_modules
srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS \
        mdrun_mpi \
        -ntomp $OMP_NUM_THREADS -nb gpu -bonded cpu \
        -deffnm ${NPT}/npt \
        -s ${NPT}/npt.tpr

# * * * * * * * * *
# MD production
# * * * * * * * * *
mkdir -p $MD
# - Use CPU version
gmx_modules
mpirun  -n 1 \
$GMX_CPU_PATH\gmx_mpi grompp \
    -f md.mdp \
    -c ${NPT}/${NPT}.gro \
    -r ${NPT}/${NPT}.gro \
    -t ${NPT}/${NPT}.cpt \
    -p ${BASE_NAME}.top \
    -po ${MD}/${BASE_NAME}_mdout.mdp \
    -o  ${MD}/${MD}.tpr
# - Use GPU version
gro_gpu_modules
srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS \
        mdrun_mpi \
        -ntomp $OMP_NUM_THREADS -nb gpu -bonded cpu \
        -deffnm ${MD}/${MD} \
        -s ${MD}/${MD}.tpr