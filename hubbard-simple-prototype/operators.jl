# Functionalities concerning operator algebra

#region Purely symbolic calculation of operator multiplication

"""
Note that this function may involve operators outside of `hubbard_opstr_basis`. 
If the final result of the normal ordered multiplication of `opstr_index_1` and `opstr_index_2`
involves operators outside of `hubbard_opstr_basis`, they are *not* truncated here.
They are truncated when constructing `hubbard_opstr_normal_order`.
"""
function hubbard_opstr_mul(opstr_index_1::HubbardOperatorString, opstr_index_2::HubbardOperatorString)
    
end

#endregion

#region Construct the M matrix 

hubbard_opstr_zero = zeros(ComplexF64, length(hubbard_opstr_basis))

hubbard_opstr_normal_order = Matrix{Vector{ComplexF64}}(undef, 
    hubbard_opstr_basis_length, hubbard_opstr_basis_length)

for opstr_index_1 in 1 : hubbard_opstr_basis_length
    for opstr_index_2 in 1 : hubbard_opstr_basis_length
        opstr_1 = hubbard_opstr_basis[opstr_index_1]
        opstr_2 = hubbard_opstr_basis[opstr_index_2]

        
    end
end

#endregion

#region 

#endregion