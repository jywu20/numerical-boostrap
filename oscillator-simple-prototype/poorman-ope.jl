using Plots

include("commutation.jl")

##

@show comm(x * p, x^2)
@show comm(x * p, p^2)