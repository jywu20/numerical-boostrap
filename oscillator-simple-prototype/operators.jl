using SparseArrays
import Base.setindex!, Base.getindex, Base.*, Base.+, Base.-, Base.^, Base.show

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

const default_maximal_power = 20

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

-(o1::XPOpString, o2::XPOpString) = o1 + (-1) * o2

*(::Type{OneDimXOperator}, ::Type{OneDimXOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(2), number_to_index(0)] = 1.0
    XPOpString(default_maximal_power, mat)
end

*(::Type{OneDimPOperator}, ::Type{OneDimPOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(0), number_to_index(2)] = 1.0
    XPOpString(default_maximal_power, mat)
end

*(::Type{OneDimXOperator}, ::Type{OneDimPOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(1), number_to_index(1)] = 1.0
    XPOpString(default_maximal_power, mat)
end

*(::Type{OneDimPOperator}, ::Type{OneDimXOperator}) = begin
    mat = spzeros(ComplexF64, default_maximal_power, default_maximal_power)
    mat[number_to_index(0), number_to_index(0)] = - im
    mat[number_to_index(1), number_to_index(1)] = 1.0
    XPOpString(default_maximal_power, mat)
end

function *(::Type{OneDimXOperator}, o2::XPOpString)
    K = o2.K
    XPOpString(K, [
        spzeros(1, K) ;
        o2.coefficient[1 : end - 1, :] 
    ])
end

function *(o1::XPOpString, ::Type{OneDimPOperator})
    K = o1.K
    XPOpString(K, [
        spzeros(K, 1)   o1.coefficient[:, 1 : end - 1] 
    ])
end

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

op_sequence(x_number, p_number) = [repeat([x], x_number) ; repeat([p], p_number) ]

function *(ops1::XPOpString, ops2::XPOpString)
    row_idxs1, col_idxs1, entries1 = findnz(ops1.coefficient)
    row_idxs2, col_idxs2, entries2 = findnz(ops2.coefficient)

    K = max(ops1.K, ops2.K)
    result = XPOpString(K, number_to_index(0), number_to_index(0))

    for i in 1 : length(row_idxs1), j in 1 : length(row_idxs2)
        x_num_1 = index_to_number(row_idxs1[i])
        p_num_1 = index_to_number(col_idxs1[i])
        x_num_2 = index_to_number(row_idxs2[j])
        p_num_2 = index_to_number(col_idxs2[j])

        
    end
end
