#!/bin/bash
#SBATCH -n 1                   # Total number of mpi tasks requested
#SBATCH -c 4
#SBATCH -t 24:00:00             # Run time (hh:mm:ss)

julia jump-oscillator-4.jl