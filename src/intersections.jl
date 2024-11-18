
"""
    compare_for_overlap(i1::Interval, i2::Interval)

Compare two intervals for overlap. Returns -1 if i1 is before i2, 
1 if i1 is after i2, and 0 if they overlap.
"""
function compare_for_overlap(i1::Interval{T1,Closed,Closed},
    i2::Interval{T2,Closed,Closed}) where {T1,T2}
    if i1.last < i2.first      # i1 is before i2
        return -1
    elseif i2.last < i1.first  # i1 is after i2
        return 1
    else
        return 0               # i1 and i2 overlap
    end
end


"""
    compare_for_inclusion_1in2(i1::Interval, i2::Interval)

Compare two intervals for inclusion. Returns -2 if i1 is before i2,
2 if i1 is after i2, -1 if i1 partially overlaps i2 at the beginning,
1 if i1 partially overlaps i2 at the end, and 0 if i1 is inside i2.
"""
function compare_for_inclusion_1in2(i1::Interval{T1,Closed,Closed},
    i2::Interval{T2,Closed,Closed}) where {T1,T2}
    if i1.last < i2.first      # i1 is before i2
        return -2
    elseif i2.last < i1.first  # i1 is after i2
        return 2
    elseif i1.first < i2.first # partial overlap
        return -1
    elseif i1.last > i2.last   # partial overlap
        return 1
    else
        return 0               # i1 is inside i2 
    end
end

"""
    compare_for_inclusion_2in1(i1::Interval, i2::Interval)

Compare two intervals for inclusion. Returns -2 if i1 is before i2,
2 if i1 is after i2, -1 if i2 partially overlaps i1 at the beginning,
1 if i2 partially overlaps i1 at the end, and 0 if i2 is inside i1.
"""
function compare_for_inclusion_2in1(i1::Interval{T1,Closed,Closed},
    i2::Interval{T2,Closed,Closed}) where {T1,T2}
    if i1.last < i2.first      # i1 is before i2
        return -2
    elseif i2.last < i1.first  # i1 is after i2
        return 2
    elseif i2.first < i1.first # partial overlap
        return 1
    elseif i2.last > i1.last   # partial overlap
        return -1
    else
        return 0               # i2 is inside i1
    end
end


"""
    find_intersections(x, y, compare=compare_for_overlap)

Find intersections between two arrays of intervals. The `compare` function
is used to compare two intervals. The default is `compare_for_overlap`.

Returns an array of arrays, where the i-th element contains the indices of
intervals in `y` that intersect with the i-th interval in `x`.
"""
function find_intersections(x, y, compare=compare_for_overlap)
    sx = sortperm(x)
    sy = sortperm(y)
    nx = 1
    qf = 1
    ql = 0
    qi = 1
    rx = 0
    ry = 0

    results = [Int[] for i in 1:length(x)]
    
    while true
        if qi > ql  # queue is empty
            if qi <= length(sy) # add interval in y
                ql += 1
            else # no more intervals in y, advance x
                nx += 1
                nx > length(sx) && break
                qi = qf
            end
        else
            rx = sx[nx]
            ry = sy[qi]
            c = compare(x[rx], y[ry])
            qi += 1
            if c < 0 # advance x:  [---IX---]  [---IY---]
                nx += 1  
                nx > length(sx) && break
                qi = qf
            elseif c == 0 # intersection found
                push!(results[rx], ry)
            else  # c > 0   [---IY---]  [---IX---]
                if qi == qf + 1 # noting else can intersect front of queue
                    qf += 1
                end
            end
        end
    end
    results
end

