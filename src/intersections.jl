
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

The intervals in `x` and `y` do not need to be sorted. However the function
will sort them internally, so for repeated calls it is more efficient to sort 
them before calling this function.
"""
function find_intersections(x, y, compare=compare_for_overlap)
    sortedIndices_x = sortperm(x)
    sortedIndices_y = sortperm(y)
    pos_x = 1
    queue_first = 1
    queue_last = 0
    queue_i = 1
    index_x = 0
    index_y = 0

    results = [Int[] for i in 1:length(x)]
    
    while true
        if queue_i > queue_last  # queue is empty
            if queue_i <= length(y) # add interval in y
                queue_last += 1
            else # no more intervals in y, advance x
                pos_x += 1
                pos_x > length(x) && break
                queue_i = queue_first
            end
        else
            index_x = sortedIndices_x[pos_x]
            index_y = sortedIndices_y[queue_i]
            c = compare(x[index_x], y[index_y])
            queue_i += 1
            if c < 0 # advance x:  [---IX---]  [---IY---]
                pos_x += 1  
                pos_x > length(sortedIndices_x) && break
                queue_i = queue_first
            elseif c == 0 # intersection found
                push!(results[index_x], index_y)
            else  # c > 0   [---IY---]  [---IX---]
                if queue_i == queue_first + 1 # noting else can intersect front of queue
                    queue_first += 1
                end
            end
        end
    end
    results
end

