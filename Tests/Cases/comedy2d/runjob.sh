#!/bin/bash

#$ -cwd
#$ -N cement2d_16
#$ -pe openmpi 16
#$ -v OPENMPI=/opt/mpi/openmpi-1.2.5,LD_LIBRARY_PATH=/opt/gcc/gcc-4.3.4/lib64:/opt/mpi/openmpi-1.2.5/lib,SGE_QMASTER_PORT 
# 
# Make sure that OPENMPI and LD_LIBRARY_PATH are set to the same values as in the active comment above.   
OPENMPI=/opt/mpi/openmpi-1.2.5
export LD_LIBRARY_PATH=/opt/gcc/gcc-4.3.4/lib64:$OPENMPI/lib
MPIRUN=$OPENMPI/bin/mpirun

#####################################################################################
MACHINE_FILE=$TMPDIR/machinefile
awk '/^merlin/ {print $1" slots="$2}' $PE_HOSTFILE > $MACHINE_FILE
#####################################################################################
MY_DATE=`date`
echo
echo "Running on HOSTNAME=$HOSTNAME  on $MY_DATE"
echo "PE_HOSTFILE=$PE_HOSTFILE"
cat $PE_HOSTFILE
echo
echo "MACHINE_FILE=$MACHINE_FILE"
cat $MACHINE_FILE
echo
echo "NSLOTS=$NSLOTS"
echo
echo "Running environment:"
env
echo "================================================================"
echo "ls -lA $MPIRUN"
ls -lA $MPIRUN
echo 
echo "ldd $MPIRUN"
ldd $MPIRUN
echo 
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
echo 
echo "ulimit -a"
ulimit -a 
echo 
whoami
pwd
echo "================================================================"

#####################################################################################

echo "#0# #############################################################"
MPICMD="$MPIRUN -machinefile $MACHINE_FILE -np $NSLOTS hostname"
echo "Run test: $MPICMD"
$MPICMD
 
# echo "###############################################################" 
# MPICMD="$MPIRUN --prefix $OPENMPI -x LD_LIBRARY_PATH -machinefile $MACHINE_FILE -np $NSLOTS ./get_ppid.sh"
# echo "A. Running: strace -f -e trace=process $MPICMD"
# strace -f -e trace=process $MPICMD 
# echo "Running: $MPICMD"
# $MPICMD 

echo "###############################################################"
# Set your CMD and ARGS here: 
# CMD=/bin/hostname
CMD=$HOME/OGS-GEM/sources/BUILDMPI/bin/ogs
ARGS='cement2d' 
echo "ldd $CMD"
ldd $CMD
echo 

MPICMD="$MPIRUN --prefix $OPENMPI -x LD_LIBRARY_PATH -machinefile $MACHINE_FILE -np $NSLOTS $CMD $ARGS"
echo 
echo "Run: $MPICMD"
$MPICMD |grep "Time\|failed"

################################################################################
