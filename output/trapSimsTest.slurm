#!/bin/bash

#SBATCH -J mpi_job
#SBATCH -p debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH -t 0-00:04:00
#SBATCH --job-name=UCNt_Trap_Sims
#SBATCH -p debug 
#SBATCH --export=ALL

export SIMS=/N/u/frangonz/BigRed3/UCNtau-Simulations/simulations/sim
export OUTS=/N/slate/frangonz/Sims/Test/

srun ${SIMS} --holdt=20 --dt=0.0005 --ntraj=512 --file=/N/slate/frangonz/Sims/Test/noabs_clean_test.bin > ${OUTS}/tracking_sims_test.txt

