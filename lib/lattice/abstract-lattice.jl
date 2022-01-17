"""
The abstract type of all lattices.
Note that coordinate types are *not* included in the type parameter of every abstract lattice type, 
as one lattice may have several coordinate systems.
"""
abstract type AbstractLattice{SiteType} end

function sites(::L) where {L <: AbstractLattice}
    error("List of all sites not defined for lattice type $L.")
end

function site_number(::L) where {L <: AbstractLattice}
    error("Total site number not defined for lattice type $L.")
end

function coord_to_site(::L, ::T) where {L <: AbstractLattice, T}
    error("Coordinate of type $T for sites not supported for a lattice with type $L.")
end

function site_to_coord(::L, ::S) where {S, L <: AbstractLattice{S}}
    error("No default site to coordinate rule defined for lattice type $L.")
end

function site_to_coord(::Type{T}, ::L, ::S) where {T, S, L <: AbstractLattice{S}}
    error("No $T coordinates to site rule defined for lattice type $L.")
end

function nearest_neighbors(::L, ::S) where {S, L <: AbstractLattice{S}}
    error("No nearest neighbors defined for lattice type $L.")
end

function nearest_neighbors(::Type{C}, lattice::L, site::S) where {C, S, L <: AbstractLattice{S}}
    convert(C, nearest_neighbors(lattice, site))
end

"""
The abstract type of all lattices with plaquattes. 
"""
abstract type AbstractLatticeWithPlaquattes{SiteType, BondType, PlaquatteType} <: AbstractLattice{SiteType} end

function bonds(::L) where {L <: AbstractLatticeWithPlaquattes}
    error("List of all bonds not defined for lattice type $L.")
end

function bond_number(::L) where {L <: AbstractLatticeWithPlaquattes}
    error("Total bond number not defined for lattice type $L.")
end

function coord_to_bond(::L, ::T) where {L <: AbstractLatticeWithPlaquattes, T}
    error("Coordinate of type $T for bonds not supported for a lattice with type $L.")
end

function bond_to_coord(::L) where {L <: AbstractLatticeWithPlaquattes}
    error("No default bond to coordinate rule defined for lattice type $L.")
end

function bond_to_coord(::Type{T}, ::L) where {L <: AbstractLatticeWithPlaquattes, T}
    error("No $T coordinates to bond rule defined for lattice type $L.")
end

function bond_to_sites(::L, ::B) where {S, P, B, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    error("Bond to sites rule not defined for lattice type $L.")
end

function bond_to_sites(::Type{C}, lattice::L, bond::B) where {S, B, P, C, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    convert(C, bond_to_sites(lattice, bond))
end

function sites_to_bond(::L, ::S, ::S)::B where {S, B, P, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    error("Sites to bond rule not defined for lattice type $L.")
end

function plaquattes(::L) where {L <: AbstractLatticeWithPlaquattes}
     error("Plaquattes not defined for lattice type $L.")
end

function plaquatte_to_bonds(::L, ::P) where {B, P, S, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    error("Plaquatte to bonds rule not defined for lattice type $L.")
end

function plaquatte_to_bonds(::Type{C}, lattice::L, plaquatte::P) where {B, P, S, C, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    convert(C, plaquatte_to_bonds(lattice, plaquatte))
end

function bonds_around_site(::L, ::S) where {B, P, S, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    error("Bonds around site not defined for lattice type $L.")
end

function bonds_around_site(::Type{C}, lattice::L, site::S) where {B, P, C, S, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    convert(C, bonds_around_site(lattice, site))
end

function plaquatte_containing_bond(::L, ::B) where {B, P, S, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    error("Plaquattes with bond not defined for lattice type $L.")
end

function plaquatte_containing_bond(::Type{C}, lattice::L, bond::B) where {B, P, C, S, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    convert(C, plaquatte_containing_bond(lattice, bond))
end

function plaquatte_to_corners(::L, ::P) where {B, P, S, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    error("Corners of plaquatte not defined for lattice type $L.")
end

function plaquatte_to_corners(::Type{C}, lattice::L, plaquatte::P) where {C, B, P, S, L <: AbstractLatticeWithPlaquattes{S, B, P}}
    convert(C, plaquatte_to_corners(lattice, plaquatte))
end