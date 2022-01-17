"""
Fermionic opertor string with spin-1/2. `False` means ↑, `True` means ↓. Note that in Julia 1 is euivalent to 
`True` (which can be checked by `BitArray([1, 0, 1])[1] == True`), and usually we denote 0 to ↑.
"""
struct ElectronString{S} where {S}
    creation_positions::S[]
    annihilation_positions::S[]
    create_spins::Bool[]
    annihilation_spins::Bool[]
end

## Since we are going to work with an infinity square lattice, we are NOT going to use any finite lattice.

