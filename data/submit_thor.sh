#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -q default
#PBS -N vaspjob
#PBS -r n
source /opt/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64

cd ${PBS_O_WORKDIR}

NPROCS="$(wc -l < ${PBS_NODEFILE} | tr -d '[:blank:]')"

MYMPIPROG="${HOME}/vtst/bin/vasp_std"

mpirun -np ${NPROCS} ${MYMPIPROG} | tee _JOB_OUTPUT.txt
