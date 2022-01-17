# Definitions used in DQMC, with poor performance, used for benchmark

"""
B-matrices by definition, without any numerical acceleration.
"""
function B_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    Δτ = model.Δτ
    t = model.t
    n_site = site_number(σ.lattice)
    T = zeros(n_site, n_site)
    for b in bonds(lattice)
        i, j = bond_to_sites(lattice, b)
        T[i, j] = T[j, i] = σ[b, τ]
    end
    exp(Δτ * t * T)
end

function B_τ_0_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    n_site = site_number(σ.lattice)
    B = Matrix(I, n_site, n_site)
    for τ′ in 1 : τ
        B = B_def(model, σ, τ′) * B        
    end
    B
end

function B_β_τ_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    n_site = site_number(σ.lattice)
    n_τ = σ.n_τ
    B = Matrix(I, n_site, n_site)
    for τ′ in n_τ : -1 : τ + 1
        B = B * B_def(model, σ, τ′)
    end
    B
end

function G_τ_τ_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    inv(I + B_τ_0_def(model, σ, τ) * B_β_τ_def(model, σ, τ))
end

function Δ_mat_def(model::Z2SpinlessFermionSimpleDQMC{F}, σ::Z2GaugeFieldDPI{L, V}, b, τ) where {
    L <: AbstractLattice, V, F <: AbstractFloat
}
    t = model.t
    lattice = σ.lattice
    n_site = site_number(lattice)
    i, j = bond_to_sites(lattice, b)
    Δ = zeros(F, n_site, n_site)
    if σ[b, τ] == one(V) 
        Δ[i, i] = Δ[j, j] = cosh(2 * Δτ * t) - 1
        Δ[i, j] = Δ[j, i] = - sinh(2 * Δτ * t)
    else
        Δ[i, i] = Δ[j, j] = cosh(- 2 * Δτ * t) - 1
        Δ[i, j] = Δ[j, i] = - sinh(- 2 * Δτ * t)
    end
    Δ
end

function T_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, τ)
    t = model.t
    n_site = site_number(σ.lattice)
    T = zeros(n_site, n_site)
    for b in bonds(lattice)
        i, j = bond_to_sites(lattice, b)
        T[i, j] = T[j, i] = σ[b, τ]
    end
    - t * T
end

function weight_ratio_def(model::Z2SpinlessFermionSimpleDQMC, σ::Z2GaugeFieldDPI, b, τ)
    σ′ = copy(σ)    
    σ′[b, τ] *= -1
    z = det(I + B_β_τ_def(model, σ, 0))
    z′ = det(I + B_β_τ_def(model, σ′, 0))
    z′ / z
end