struct SquareLattice2DPeriodic <: AbstractLatticeWithPlaquattes{Int, Int, Int}
    n_side::Int
    site_list::Vector{Tuple{Int, Int}}
    inverse_site_list::Matrix{Int} 
    nn_list::Vector{Tuple{Int, Int, Int, Int}}
    bond_list::Matrix{Int}
    inverse_bond_list::Vector{Tuple{Int, Int}}
    bonds_around_site_list::Vector{Tuple{Int, Int, Int, Int}}
    plaquatte_corners_list::Vector{Tuple{Int, Int, Int, Int}}
    plaquatte_bonds_list::Vector{Tuple{Int, Int, Int, Int}}
    bond_sharing_plaquattes_list::Vector{Tuple{Int, Int}}
end

function site_to_coord(n_side::Int, i::Int)::Tuple{Int, Int}
    x = i % n_side
    if x == 0
        x = n_side
    end
    y = 1 + Int((i - x) / n_side)
    
    (x, y)
end

function coord_to_site(n_side::Int, x::Int, y::Int)::Int
    (y - 1) * n_side + x 
end

function nearest_neighbors(n_side::Int, i::Int)::Tuple{Int, Int, Int, Int}
    x, y = site_to_coord(n_side, i)
    zero_to_max(i) = i > 0 ? i : i + n_side
    (
        coord_to_site(n_side, zero_to_max(mod(x + 1, n_side)), y),
        coord_to_site(n_side, zero_to_max(mod(x - 1, n_side)), y),
        coord_to_site(n_side, x, zero_to_max(mod(y + 1, n_side))),
        coord_to_site(n_side, x, zero_to_max(mod(y - 1, n_side)))
    )
end

function bond_to_sites(n_side::Int, bond::Int)
    if bond > n_side^2
        bond = bond - n_side^2
        x, y = site_to_coord(n_side, bond)
        return (bond, coord_to_site(n_side, back_into_range(x + 1, n_side), y))
    end
    x, y = site_to_coord(n_side, bond)
    (bond, coord_to_site(n_side, x, back_into_range(y + 1, n_side)))
end

function SquareLattice2DPeriodic(n_side::Int)::SquareLattice2DPeriodic 
    site_list = Vector{Tuple{Int, Int}}(undef, n_side * n_side)
    inverse_list = Matrix{Int}(undef, n_side, n_side)
    for x in 1 : n_side
        for y in 1 : n_side 
            site = coord_to_site(n_side, x, y)
            site_list[site] = (x, y)
            inverse_list[x, y] = site
        end
    end

    inverse_bond_list = Vector{Tuple{Int, Int}}(undef, 2 * n_side * n_side)
    bond_list = zeros(Int, n_side^2, n_side^2)
    for bond in 1 : 2 * n_side * n_side
        x, y = bond_to_sites(n_side, bond)
        inverse_bond_list[bond] = (x, y)
        bond_list[x, y] = bond
        bond_list[y, x] = bond
    end

    nn_list = Vector{Tuple{Int, Int, Int, Int}}(undef, n_side * n_side)
    bonds_around_site_list = Vector{Tuple{Int, Int, Int, Int}}(undef, n_side * n_side)
    for site in 1 : n_side * n_side
        nn_list[site] = nearest_neighbors(n_side, site)
        bonds_around_site_list[site] = map(x -> bond_list[x, site], nn_list[site])
    end

    plaquatte_corner_list = Vector{Tuple{Int, Int, Int, Int}}(undef, n_side * n_side)
    plaquatte_bonds_list = Vector{Tuple{Int, Int, Int, Int}}(undef, n_side * n_side)
    for plaquatte in 1 : n_side^2
        x, y = site_to_coord(n_side, plaquatte)
        i, j, k, l = plaquatte_corner_list[plaquatte] = (
            coord_to_site(n_side, x, y),
            coord_to_site(n_side, back_into_range(x + 1, n_side), y),
            coord_to_site(n_side, x, back_into_range(y + 1, n_side)),
            coord_to_site(n_side, back_into_range(x + 1, n_side), back_into_range(y + 1, n_side))
        )
        plaquatte_bonds_list[plaquatte] = (bond_list[i, j], bond_list[i, k], bond_list[j, l], bond_list[k, l])
    end

    bond_sharing_plaquatte_list = Vector{Tuple{Int, Int}}(undef, 2 * n_side * n_side)
    for bond in 1 : n_side^2
        x, y = site_to_coord(n_side, bond)
        x_moved_site = coord_to_site(n_side, back_into_range(x - 1, n_side), y)
        y_moved_site = coord_to_site(n_side, x, back_into_range(y - 1, n_side))
        bond_sharing_plaquatte_list[bond] = (x_moved_site, bond)
        bond_sharing_plaquatte_list[bond + n_side^2] = (y_moved_site, bond)
    end

    SquareLattice2DPeriodic(n_side, 
        site_list, inverse_list, 
        nn_list, bond_list, inverse_bond_list, 
        bonds_around_site_list, 
        plaquatte_corner_list, plaquatte_bonds_list, bond_sharing_plaquatte_list)
end

sites(lattice::SquareLattice2DPeriodic) = 1 : lattice.n_side^2 
site_number(lattice::SquareLattice2DPeriodic) = lattice.n_side^2
coord_to_site(lattice::SquareLattice2DPeriodic, coord) = lattice.inverse_site_list[coord...]
site_to_coord(lattice::SquareLattice2DPeriodic, site::Int) = lattice.site_list[site]

bond_number(lattice::SquareLattice2DPeriodic) = 2 * lattice.n_side^2
bonds(lattice::SquareLattice2DPeriodic) = 1 : 2 * lattice.n_side^2 

nearest_neighbors(lattice::SquareLattice2DPeriodic, site::Int) = lattice.nn_list[site]
bonds_around_site(lattice::SquareLattice2DPeriodic, site::Int) = lattice.bonds_around_site_list[site]
bond_to_sites(lattice::SquareLattice2DPeriodic, bond::Int) = lattice.inverse_bond_list[bond]
sites_to_bond(lattice::SquareLattice2DPeriodic, site1::Int, site2::Int) = lattice.bond_list[site1, site2]

plaquattes(lattice::SquareLattice2DPeriodic) = sites(lattice)
plaquatte_to_corners(lattice::SquareLattice2DPeriodic, plaquatte::Int) = lattice.plaquatte_corners_list[plaquatte]
plaquatte_to_bonds(lattice::SquareLattice2DPeriodic, plaquatte::Int) = lattice.plaquatte_bonds_list[plaquatte]
plaquatte_containing_bond(lattice::SquareLattice2DPeriodic, bond::Int) = lattice.bond_sharing_plaquattes_list[bond]