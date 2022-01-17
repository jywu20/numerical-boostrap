# Initialize auxiliary field configurations

function Z2SpinlessFermionSimpleAuxField(::Type{F}, model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI
    ) where {F <: AbstractFloat}
    lattice = σ.lattice
    n_sites = site_number(lattice)
    n_τ = time_step_number(σ)
    G = zeros(F, n_sites, n_sites, n_τ)
    aux = Z2SpinlessFermionSimpleAuxField(σ, G, 0)
    for τ in 1 : n_τ
        G_τ_τ!(model, aux, τ)
    end
    aux
end