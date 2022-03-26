#!/bin/bash
#SBATCH -n 1                   # Total number of mpi tasks requested
#SBATCH -c 4
#SBATCH -t 14:00:00             # Run time (hh:mm:ss)

julia jump-oscillator-3-cosmo-opsrel-1e-10-fix-x2.jl