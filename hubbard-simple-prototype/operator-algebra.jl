# Functionalities concerning operator algebra
using QuantumAlgebra

#region Construct all operator involved (see labbook.md#2022.4.7 for some discussion), and do purely 
# symbolic calculation of operator multiplication

@fermion_ops c

↑ = 1
↓ = -1

begin
    local HubbardOperator = Tuple{Symbol, Int, Int}
    local HubbardOperatorString = Vector{HubbardOperator}
    local qualified_opstr_create_half = HubbardOperatorString[]
    local qualified_opstr_annihilate_half = HubbardOperatorString[]

    #region Truning spin configurations at each site into operator string
    for current_op_str in qualified_opstr_site_configuration
        current_create_half = HubbardOperator[]
        current_annihilate_half = HubbardOperator[]
        for i in 1 : length(site_list)
            @match current_op_str[i] begin
                :no => Nothing
                :up => begin
                    # QuantumAlgebra recognize :↑ and :↓ as indices, so δ(↑↓) may appear in normal ordered 
                    # operators. We need to avoid this case. Using string instead of symbol to label spins
                    # creates the same problem. Nor can we label spins as 1/2 or 1//2, since QuantumAlgebra
                    # don't support float or rational. Nor can we use c↑.
                    # The only choice is to use 1 to denote ↑ and -1 ↓.
                    push!(current_create_half, (:cdag, i, ↑));
                    push!(current_annihilate_half, (:c, i, ↑))
                end
                :dn => begin
                    push!(current_create_half, (:cdag, i, ↓));
                    push!(current_annihilate_half, (:c, i, ↓))
                end
                :both => begin
                    # Later we need to multiply all operators together, and QuantumAlgebra places ↓ before ↑,
                    # so we also need to do so
                    push!(current_create_half, (:cdag, i, ↓));
                    push!(current_create_half, (:cdag, i, ↑));
                    push!(current_annihilate_half, (:c, i, ↓));
                    push!(current_annihilate_half, (:c, i, ↑));
                end
            end
        end
        push!(qualified_opstr_create_half, current_create_half)
        push!(qualified_opstr_annihilate_half, current_annihilate_half)
    end
    #endregion

    #region Imposing the l(O) ≤ K constrint
    hubbard_opstr_basis = reshape(map(seq_pair -> begin
            cdag_seq, c_seq = seq_pair
            [cdag_seq..., c_seq...]
        end, 
        Iterators.product(qualified_opstr_create_half, qualified_opstr_annihilate_half) |> collect), 
        length(qualified_opstr_create_half)^2)

    local hubbard_opstr_size(opstr) = begin
        if opstr == []
            return 0.0
        end
        sum(map(op -> site_norm_1_from_center[op[2]], opstr)) + length(opstr)
    end

    filter!(opstr -> hubbard_opstr_size(opstr) ≤ K, hubbard_opstr_basis)
    hubbard_opstr_basis = convert(Vector{Vector{Tuple{Symbol, Int, Int}}}, hubbard_opstr_basis)

    #endregion

    #region output

    hubbard_opstr_basis_length = length(hubbard_opstr_basis)

    hubbard_opstr_basis_sites = map(opstr -> map(op -> op[2], opstr), hubbard_opstr_basis)

    hubbard_opstr_basis_size = map(hubbard_opstr_size, hubbard_opstr_basis)

    # Turning operator labels into actual QuantumAlgebra objects. Redefine `hubbard_opstr_basis`.
    # 
    # Note that multiplication of two elements in `hubbard_opstr_basis` may involve operators 
    # outside of `hubbard_opstr_basis`. 
    # If the final result of the normal ordered multiplication of `opstr_index_1` and `opstr_index_2`
    # involves operators outside of `hubbard_opstr_basis`, they are *not* truncated here.
    # They are truncated when constructing `hubbard_opstr_normal_order`.
    hubbard_opstr_basis = map(hubbard_opstr_basis) do opstr_labels
        if length(opstr_labels) == 0
            return one(QuExpr)
        end
        map(opstr_labels) do op_label
            if op_label[1] == :c
                return c(op_label[2], op_label[3])
            end
            return cdag(op_label[2], op_label[3])
        end |> prod
    end

    hubbard_opstr_index = Dict{QuExpr, Int}()
    for opstr_index in 1 : length(hubbard_opstr_basis)
        hubbard_opstr_index[hubbard_opstr_basis[opstr_index]] = opstr_index
    end

    if show_hubbard_opstr_basis
        open(working_path * output_name, "a") do file
            for opstr in hubbard_opstr_basis
                println(file, opstr)
            end
            println()
        end
    end

    hubbard_opstr_zero = zeros(ComplexF64, length(hubbard_opstr_basis))

    function hubbard_opstr_coefficients(opstr)
        res = copy(hubbard_opstr_zero)
        for basis_coefficient_pair in opstr.terms
            # Note: `opstr.terms`'s key type is QuantumAlgebra.QuTerm, a private type, 
            # so here we need to get the QuTerm object corresponding to `basis`
            basis = QuExpr(Dict(basis_coefficient_pair[1] => 1))

            coefficient = basis_coefficient_pair[2]
            res[hubbard_opstr_index[basis]] = coefficient
        end
        res
    end
    
    #endregion
end

#endregion

#region Construct the M matrix 

hubbard_opstr_normal_order = Matrix{Vector{ComplexF64}}(undef, 
    hubbard_opstr_basis_length, hubbard_opstr_basis_length)

# Indices of operators qualified to span the M matrix in `hubbard_opstr_basis`
# The condition: O_i is qualified if ⟨ O†_i O_j ⟩ is in `hubbard_opstr_basis` for any O_j
M_mat_spanning_opstr_indices = filter(collect(1 : hubbard_opstr_basis_length)) do opstr_index
    hubbard_opstr_basis_size[opstr_index] ≤ K / 2
end

M_coefficient = Matrix{Vector{Float64}}(undef, 
    length(M_mat_spanning_opstr_indices), length(M_mat_spanning_opstr_indices))

for (i, opstr_index_1) in enumerate(M_mat_spanning_opstr_indices)
    for (j, opstr_index_2) in enumerate(M_mat_spanning_opstr_indices)
        opstr_1 = hubbard_opstr_basis[opstr_index_1]
        opstr_2 = hubbard_opstr_basis[opstr_index_2]

        Mij = normal_form(opstr_1' * opstr_2)
        M_coefficient[i, j] = hubbard_opstr_coefficients(Mij)
    end
end

#endregion

#region Construct equational constraints



#endregion