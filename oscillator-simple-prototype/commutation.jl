function comm(op_string_1::AbstractVector{O}, op_string_2::AbstractVector{O}; times = *) where {O}
    len1 = length(op_string_1)
    len2 = length(op_string_2)

    terms = O[]

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

