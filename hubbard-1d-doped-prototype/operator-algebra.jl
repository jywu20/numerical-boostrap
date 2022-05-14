# Functionalities concerning operator algebra
using QuantumAlgebra
using Match
using SparseArrays

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
                :no => nothing
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

    hubbard_opstr_basis_sites = map(opstr -> map(op -> op[2], opstr) |> Set |> collect, hubbard_opstr_basis)

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
            for (i, opstr) in enumerate(hubbard_opstr_basis)
                println(file, "hubbard_opstr_basis[$i]  =  ", opstr)
            end
            println(file)
        end
    end

    hubbard_opstr_zero = spzeros(length(hubbard_opstr_basis))

    function hubbard_opstr_coefficients(opstr)
        res = copy(hubbard_opstr_zero)
        for basis_coefficient_pair in opstr.terms
            # Note: `opstr.terms`'s key type is QuantumAlgebra.QuTerm, a private type, 
            # so here we need to get the QuTerm object corresponding to `basis`
            basis = QuExpr(Dict(basis_coefficient_pair[1] => 1))

            coefficient = basis_coefficient_pair[2]
            if ! haskey(hubbard_opstr_index, basis)
                return nothing
            end
            res[hubbard_opstr_index[basis]] = coefficient
        end
        res
    end
    
    #endregion
end

#endregion

#region Construct the M matrix 

open(working_path * output_name, "a") do file
    println(file, "Start to build M matrix.")
end

#hubbard_opstr_normal_order = Matrix{SparseVector{Float64, Int64}}(undef, 
#    hubbard_opstr_basis_length, hubbard_opstr_basis_length)

# Indices of operators qualified to span the M matrix in `hubbard_opstr_basis`
# The condition: O_i is qualified if ⟨ O†_i O_j ⟩ is in `hubbard_opstr_basis` for any O_j
M_mat_spanning_opstr_indices = filter(collect(1 : hubbard_opstr_basis_length)) do opstr_index
    hubbard_opstr_basis_size[opstr_index] ≤ K / 2
end

open(working_path * output_name, "a") do file
    println(file, "Spanning operators chosen. Total number = $(length(M_mat_spanning_opstr_indices))")
end

M_coefficient = Matrix{SparseVector{Float64, Int64}}(undef, 
    length(M_mat_spanning_opstr_indices), length(M_mat_spanning_opstr_indices)
)

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

begin 
    H_hubbard = zero(QuExpr)
    local displacements = [1, -1]
    for (i, i_pos) in enumerate(site_list)
        # println(H_hubbard)
        # The U-term: U n_{i↑} n_{i↓}
        global H_hubbard += U * normal_form(cdag(i, ↑) * c(i, ↑) * cdag(i, ↓) * c(i, ↓))
        # The t-term
        for d in displacements
            j_pos = i_pos + d
            if haskey(inverse_list, j_pos)
                j = inverse_list[j_pos]
                global H_hubbard += - t * (cdag(j, ↑) * c(i, ↑) + cdag(j, ↓) * c(i, ↓))
            end
        end
    end
end

if show_hamiltonian
    open(full_output_name, "a") do file
        println(file, "H = ", H_hubbard)
        println(file)
    end
end

open(full_output_name, "a") do file
    println(file, "Constraints caused by H:")
    println(file)
end

H_constraints_coefficients = SparseVector{Float64, Int64}[]
for opstr_basis_index in 1 : hubbard_opstr_basis_length
    constraint_op = comm(H_hubbard, hubbard_opstr_basis[opstr_basis_index]) |> normal_form
    coefficients = constraint_op |> hubbard_opstr_coefficients
    if coefficients !== nothing
        push!(H_constraints_coefficients, coefficients)

        if show_constraints
            open(full_output_name, "a") do file
                # Show the constraint and what operator causes it
                println(file, hubbard_opstr_basis[opstr_basis_index], "  =>  ", constraint_op, " = 0")
            end
        end
    end
end

#endregion

#region Finding operators related with symmetry operations 

N = sum(i -> cdag(i, 1) * c(i, 1) + cdag(i, -1) * c(i, -1), 1 : length(site_list)) |> normal_form
Sz = sum(i -> cdag(i, 1) * c(i, 1) - cdag(i, -1) * c(i, -1), 1 : length(site_list)) |> normal_form

particle_number_constraint_ops = map(hubbard_opstr_basis) do op
    comm_res = comm(op, N) |> normal_form
    if length(comm_res.terms) == 1
        op_val_pair = collect(comm_res.terms)[1]
        return QuExpr(Dict(op_val_pair[1] => 1))
    end
end

filter!(x -> ! isnothing(x), particle_number_constraint_ops)
particle_number_constraint_ops = convert(Vector{QuExpr}, particle_number_constraint_ops)

spin_constraint_ops = map(hubbard_opstr_basis) do op
    comm_res = comm(op, Sz) |> normal_form
    if length(comm_res.terms) == 1
        op_val_pair = collect(comm_res.terms)[1]
        return QuExpr(Dict(op_val_pair[1] => 1))
    end
end

filter!(x -> ! isnothing(x), spin_constraint_ops)
spin_constraint_ops = convert(Vector{QuExpr}, spin_constraint_ops)

function hubbard_opstr_translate(op, Δx)
    opstr_string = string(op)
    opstr_string = replace(opstr_string, "†" => "dag")
    opstr_each_site_strings = split(opstr_string)

    opstr_factors = QuExpr[]
    for str in opstr_each_site_strings
        str = str[1 : end - 2]
        spin = 1

        if str[end] == '-'
            str = str[1 : end - 1]
            head, site_idx = split(str, "(")
            spin = -1 
        else 
            head, site_idx = split(str, "(")
            spin = 1
        end

        if ! haskey(inverse_list, site_list[parse(Int, site_idx)] + Δx)
            return nothing
        end
        i = inverse_list[site_list[parse(Int, site_idx)] + Δx]
        
        if head == "c"
            push!(opstr_factors, c(i, spin))
        else 
            push!(opstr_factors, cdag(i, spin))
        end
    end

    prod(opstr_factors)
end

function hubbard_opstr_reflect(op)
    opstr_string = string(op)
    opstr_string = replace(opstr_string, "†" => "dag")
    opstr_each_site_strings = split(opstr_string)

    opstr_factors = QuExpr[]
    for str in opstr_each_site_strings
        str = str[1 : end - 2]
        spin = 1

        if str[end] == '-'
            str = str[1 : end - 1]
            head, site_idx = split(str, "(")
            spin = -1 
        else 
            head, site_idx = split(str, "(")
            spin = 1
        end

        if ! haskey(inverse_list, - site_list[parse(Int, site_idx)])
            return nothing
        end
        i = inverse_list[- site_list[parse(Int, site_idx)]]
        
        if head == "c"
            push!(opstr_factors, c(i, spin))
        else 
            push!(opstr_factors, cdag(i, spin))
        end
    end

    prod(opstr_factors)
end

if show_constraints
    open(full_output_name, "a") do file
        println(file)
        println(file)
        println(file, "Constraints caused by translational symmetry:")
    end
end

translational_constraint_coefficients = map(hubbard_opstr_basis[2:end]) do op
    translated_op = hubbard_opstr_translate(op, 1)
    if translated_op === nothing
        return
    end
    cons_ops = op - translated_op
    if show_constraints
        open(full_output_name, "a") do file
            println(file, op, "  =>  ", cons_ops, " = 0")
        end
    end
    cons_ops |> hubbard_opstr_coefficients
end 

filter!(translational_constraint_coefficients) do c
    c !== nothing
end

if show_constraints
    open(full_output_name, "a") do file
        println(file)
        println(file)
        println(file, "Constraints caused by reflectional symmetry:")
    end
end

reflectional_constraints_coefficients = map(hubbard_opstr_basis[2:end]) do op
    cons_ops = op - hubbard_opstr_reflect(op)
    if show_constraints
        open(full_output_name, "a") do file
            println(file, op, "  =>  ", cons_ops, " = 0")
        end
    end
    cons_ops |> hubbard_opstr_coefficients
end 

filter!(reflectional_constraints_coefficients) do c
    c !== nothing
end

#endregion