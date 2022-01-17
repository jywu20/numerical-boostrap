#region DQMC of fermions coupled to a Z2 gauge field.

# Note that since the fermions in the hopping Hamiltonian H_{hopping} = - t \sum_{⟨i, j⟩} σᶻ_{ij} c†_i c_j
# can be directly integrated out, the contribution of fermions in the discrete path integral is just the 
# determinant of B-matrices, and the "auxiliary field" here is just the Z2 gauge field itself.

"""
`A` is the type of the auxiliary field itself.

`M` is the type of matrices used in DQMC.
"""
abstract type AbstractFermionAuxField{L <: AbstractLattice, 
    A <: AbstractDiscretePathIntegralConfiguration, 
    F <:AbstractFloat} <: AbstractDiscretePathIntegralConfiguration end

"""
type parameters:
- `L`: the lattice.
- `V`: type of field value at each bond, for example `Int`
- `F`: the type of float used, for example `Float64` or `Float32`. 
  Note that `F` must support functions like `sin` or `exp`.
- `S`: a three-dimensional array type used to store auxiliary matrices, for example `Array{Float64, 3}`. 
The indexing convention: [site1, site2, imaginary time]
"""
mutable struct Z2SpinlessFermionSimpleAuxField{
    L <: AbstractLatticeWithPlaquattes, V, F <: AbstractFloat
} <: AbstractFermionAuxField{L, Z2GaugeFieldDPI{L, V}, F}
    σ::Z2GaugeFieldDPI{L, V}
    G::Array{F, 3}
    propagate_count::Int
end

struct Z2SpinlessFermionSimpleDQMC{F <: AbstractFloat}
    t::F
    Δτ::F
    β::F
    n_wrap::Int
end

function Z2SpinlessFermionSimpleDQMC(::Type{F}, σ::Z2GaugeFieldDPI, n_wrap::Int, t::F, Δτ::F) where {
    F <: AbstractFloat
}
    n_τ = time_step_number(σ)

    n_τ = time_step_number(σ)
    β = Δτ * n_τ
    Z2SpinlessFermionSimpleDQMC{F}(t, Δτ, β, n_wrap)
end

include("defs.jl")
include("b-mat.jl")
include("green.jl")
include("bond.jl")
include("update.jl")
include("init.jl")