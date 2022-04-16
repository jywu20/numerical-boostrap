#!/bin/bash
#SBATCH -n 1                   # Total number of mpi tasks requested
#SBATCH -c 4
#SBATCH -t 48:00:00             # Run time (hh:mm:ss)

julia main.jl