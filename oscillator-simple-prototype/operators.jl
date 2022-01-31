using SparseArrays
import Base.setindex!, Base.getindex, Base.*, Base.+, Base.^, Base.show

struct OneDimXOperator 
end

x = OneDimXOperator 

struct OneDimPOperator 
end

p = OneDimPOperator

"""
An operator string made of x and p. x is placed in front of p.
"""
struct XPOpString 
    # Maximal power + 1 (for example if K = 10, x^0, x^1, ..., x^9 are allowed)
    K::Int
    coefficient::SparseMatrixCSC{ComplexF64}
end

const default_maximal_power = 10

function XPOpString(K::Int)
    XPOpString(K, spzeros(K, K))
end

index_to_number(x) = x - 1
number_to_index(x) = x + 1

function XPOpString(coefficient::AbstractMatrix)
    K1, K2 = size(coefficient)
    if K1 != K2
        error("The coefficient matrix must be a square matrix.")
    end
    K = K1 
    XPOpString(K, sparse(coefficient))
end

function XPOpString(K::Int, x_order, p_order)
    res = XPOpString(K)
    res.coefficient[x_order, p_order] = 1.0
    res
end

function XPOpString(K::Int, x_order, p_order, val)
    res = XPOpString(K)
    res.coefficient[x_order, p_order] = val
    res
end

function getindex(ops::XPOpString, x_order, p_order)
    ops.coefficient[x_order, p_order]
end

function setindex!(ops::XPOpString, val, x_order, p_order)
    ops.coefficient[x_order, p_order] = val
end

function change_maximal_power(ops::XPOpString, K_new::Integer)
    K = ops.K
    if K == K_new
        return ops
    end
    if K > K_new
        return XPOpString(K_new, ops.coefficient[1:K_new, 1:K_new])
    end

    return XPOpString(K_new, [
        ops.coefficient          spzeros(K, K_new - K) ;
        spzeros(K_new - K, K)    spzeros(K_new - K, K_new - K)
    ])
end

function plus(o1::XPOpString, o2::XPOpString; K::Integer = max(o1.K, o2.K))
    o1 = change_maximal_power(o1, K)
    o2 = change_maximal_power(o2, K)
    XPOpString(K, o1.coefficient + o2.coefficient)
end

+(o1::XPOpString, o2::XPOpString) = plus(o1, o2)
+(o1::XPOpString, ::Type{OneDimXOperator}) = begin
    new_coefficient = copy(o1.coefficient)
    new_coefficient[number_to_index(1), number_to_index(0)] += 1
    XPOpString(o1.K, new_coefficient)
end

mul(::Type{OneDimXOperator}, ::Type{OneDimXOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(2), number_to_index(0)] = 1.0
    XPOpString(default_maximal_power, mat)
end

mul(::Type{OneDimPOperator}, ::Type{OneDimPOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(0), number_to_index(2)] = 1.0
    XPOpString(default_maximal_power, mat)
end

mul(::Type{OneDimXOperator}, ::Type{OneDimPOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(1), number_to_index(1)] = 1.0
    XPOpString(default_maximal_power, mat)
end

mul(::Type{OneDimPOperator}, ::Type{OneDimXOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(0), number_to_index(0)] = - im
    mat[number_to_index(1), number_to_index(1)] = 1.0
    XPOpString(default_maximal_power, mat)
end

function mul(::Type{OneDimXOperator}, o2::XPOpString)
    K = o2.K
    XPOpString(K, [
        spzeros(1, K) ;
        o2.coefficient[1 : end - 1, :] 
    ])
end

function mul(o1::XPOpString, ::Type{OneDimPOperator})
    K = o1.K
    XPOpString(K, [
        spzeros(K, 1)   o1.coefficient[:, 1 : end - 1] 
    ])
end

*(o1::Union{XPOpString, Type{OneDimXOperator}, Type{OneDimPOperator}}, 
    o2::Union{XPOpString, Type{OneDimXOperator}, Type{OneDimPOperator}}) = mul(o1, o2)

*(x::Number, o::XPOpString) = XPOpString(o.K, x * o.coefficient)
*(o::XPOpString, x::Number) = XPOpString(o.K, x * o.coefficient)

*(x::Number, ::Type{OneDimXOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(1), number_to_index(0)] = x 
    XPOpString(default_maximal_power, mat)
end

*(x::Number, ::Type{OneDimPOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(0), number_to_index(1)] = x 
    XPOpString(default_maximal_power, mat)
end

^(::Type{OneDimXOperator}, x::Number) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(x), number_to_index(0)] = 1
    XPOpString(default_maximal_power, mat)
end

^(::Type{OneDimPOperator}, x::Number) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(0), number_to_index(x)] = 1
    XPOpString(default_maximal_power, mat)
end

function power_str(var, pow)
    if pow == 0
        return ""
    end
    if pow == 1
        return var 
    end
    return "$var^$pow"
end

function show(io::IO, ::MIME"text/plain", ops::XPOpString)
    terms = String[]
    row_idxs, col_idxs, entries = findnz(ops.coefficient)

    for i in 1 : length(row_idxs)
        row_idx = row_idxs[i]
        col_idx = col_idxs[i]
        entry = entries[i]
        push!(terms, "$entry $(power_str("x", index_to_number(row_idx))) $(power_str("p", index_to_number(col_idx)))")
    end
    print(io, join(terms, " + "))
end

comm(::Type{OneDimXOperator}, ::Type{OneDimPOperator}) = im * XPOpString(default_maximal_power, number_to_index(0), number_to_index(0))
comm(::Type{OneDimPOperator}, ::Type{OneDimXOperator}) = - im * XPOpString(default_maximal_power, number_to_index(0), number_to_index(0))
comm(::Type{OneDimXOperator}, ::Type{OneDimXOperator}) = 0 * XPOpString(default_maximal_power, number_to_index(0), number_to_index(0))
comm(::Type{OneDimPOperator}, ::Type{OneDimPOperator}) = 0 * XPOpString(default_maximal_power, number_to_index(0), number_to_index(0))

op_sequence(x_number, p_number) = [repeat([x], x_number) ; repeat([p], p_number) ]

function comm(op_string_1::AbstractVector{O}, op_string_2::AbstractVector{O}; times = *, term_type = O) where {O}
    len1 = length(op_string_1)
    len2 = length(op_string_2)

    terms = term_type[]

    for i in 1 : len1, j in 1 : len2
        head1 = op_string_1[1 : i - 1]
        head2 = op_string_2[1 : j - 1]
        tail1 = op_string_1[i + 1 : end]
        tail2 = op_string_2[j + 1 : end] 
        o1 = op_string_1[i]
        o2 = op_string_2[j]
        term_ij = reduce(times, [head1 ; head2 ; [comm(o1, o2)] ; tail2 ; tail1])
        push!(terms, term_ij)
    end

    sum(terms)
end

function comm(ops1::XPOpString, ops2::XPOpString)
    row_idxs_1, col_idxs_1, coefficients_1 = findnz(ops1.coefficient)
    row_idxs_2, col_idxs_2, coefficients_2 = findnz(ops2.coefficient)
    len1 = length(row_idxs_1)
    len2 = length(row_idxs_2)

    terms = []

    for i in 1 : len1, j in 1 : len2
        opseq1 = op_sequence(index_to_number(row_idxs_1[i]), index_to_number(col_idxs_1[i]))
        opseq2 = op_sequence(index_to_number(row_idxs_2[j]), index_to_number(col_idxs_2[j]))
        push!(terms, coefficients_1[i] * coefficients_2[j] * comm(opseq1, opseq2, term_type = XPOpString))
    end

    sum(terms)
end