using Statistics
using LinearAlgebra

include("utils.jl")
include("lattice/abstract-lattice.jl")
include("lattice/square-lattice.jl")
include("abstract-model.jl")

include("flags.jl")

include("gauge-dpi/z2/pure-z2.jl")
include("spin-dpi/tfim.jl")
include("dqmc/z2/z2-fermion.jl")
include("gauge-dpi/z2/z2-ising.jl")