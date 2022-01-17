import Base.setindex!, Base.getindex

abstract type AbstractModel end

abstract type AbstractConfigration end

function field_lattice(::C) where {C <: AbstractConfigration}
    error("The lattice of configuration type $C not defined.")
end

abstract type AbstractDiscretePathIntegralConfiguration <: AbstractConfigration end

function timeslice(::C, ::Any) where {C <: AbstractDiscretePathIntegralConfiguration}
    error("The time slice of discrete path integral configuration type $C not defined.")
end

function getindex(::C, ::Any) where {C <: AbstractDiscretePathIntegralConfiguration}
    error("Indexing of discrete path integral configuration type $C not defined.")
end

function setindex!(::C, ::Any, ::Any) where {C <: AbstractDiscretePathIntegralConfiguration}
    error("Indexing of discrete path integral configuration type $C not defined.")
end

function time_steps(::C) where {C <: AbstractDiscretePathIntegralConfiguration}
    error("Time steps of discrete path integral configuration type $C not defined.")
end

function time_step_number(::C) where {C <: AbstractDiscretePathIntegralConfiguration}
    error("Time steps of discrete path integral configuration type $C not defined.")
end
