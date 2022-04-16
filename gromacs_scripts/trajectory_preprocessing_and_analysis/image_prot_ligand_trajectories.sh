module load openmpi/4.1.2
module load cmake/3.17.3
module load intel-19.5/cmkl

GMX_CPU_PATH='/work/dc010/dc010/ricci/SOFTWARE/GROMACS/bin/'
EM='em'
SELECTION='Protein_LIG'

# Whole
echo System System | $GMX_CPU_PATH/gmx_mpi \
        trjconv -f md/md.trr \
        -s md/md.tpr \
        -fit none \
        -o md/md_WHOLE.trr \
        -n ${EM}/index.ndx \
        -pbc whole

# No jump
echo System System | $GMX_CPU_PATH/gmx_mpi \
        trjconv -f md/md_WHOLE.trr \
        -s md/md.tpr \
        -fit none \
        -o md/md_NOJUMP.trr \
        -n ${EM}/index.ndx \
        -pbc nojump

# Cluster
SELECTION='Protein_LIG'
echo "$SELECTION" "$SELECTION" System | $GMX_CPU_PATH/gmx_mpi \
        trjconv -f md/md_NOJUMP.trr \
        -s md/md.tpr \
        -fit none \
        -o md/md_CLUSTER.trr \
        -n ${EM}/index.ndx \
        -center -pbc cluster

# Mol and strip water
SELECTION='Protein_LIG'
echo "$SELECTION" "$SELECTION" | $GMX_CPU_PATH/gmx_mpi \
        trjconv -f md/md_CLUSTER.trr \
        -s md/md.tpr \
        -fit none \
        -o md/md_imaged_noWAT.trr \
        -n ${EM}/index.ndx \
        -center -pbc mol -ur compact
