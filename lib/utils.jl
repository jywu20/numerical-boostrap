function back_into_range(idx, upper)
    if idx > upper
        return idx % upper
    end
    (idx - upper) % upper + upper
end

relative_err(m1, m2) = norm(m1 - m2) / norm(m2)