
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


# `aggregator(results, index_x, index_y)` records one intersection and returns
# `true` to keep scanning the current interval in `x`, or `false` to move on to
# the next one. Aggregators that only need to know *whether* an interval in `x`
# is hit return `false` and skip the remaining comparisons for it. Because each
# aggregator returns a literal, the branch below folds away when specialized.
function _find_intersections(results, x, y, aggregator,
        compare=compare_for_overlap)
    sortedIndices_x = sortperm(x)
    sortedIndices_y = sortperm(y)

    # The live set holds positions into sortedIndices_y whose intervals may still
    # overlap the current or a later interval in x. It is an explicit list rather
    # than a contiguous window so that a dead interval can be dropped from
    # anywhere in it: y is ordered by where intervals start, but an interval dies
    # by where it ends, and those two orders are unrelated. With a window, only a
    # prefix can be dropped, so one wide interval in y pins the front and forces
    # every interval in x to rescan everything behind it.
    queue = Int[]
    sizehint!(queue, 64)
    next_y = 1

    for pos_x in eachindex(sortedIndices_x)
        index_x = sortedIndices_x[pos_x]
        ix = x[index_x]

        # admit every interval in y that starts at or before the end of ix
        while next_y <= length(sortedIndices_y) &&
                !(ix.last < y[sortedIndices_y[next_y]].first)
            push!(queue, next_y)
            next_y += 1
        end

        # Walk the live set, compacting it in place: `r` reads, `w` writes back
        # only what is still live. Costs nothing, since we walk it to compare
        # anyway. Deadness is tested geometrically rather than taken from
        # `compare`, whose positive return means "after" for one comparator but
        # "partially overlapping" -- still live -- for another.
        w = 1
        keep_scanning = true
        for r in eachindex(queue)
            qy = queue[r]
            index_y = sortedIndices_y[qy]
            iy = y[index_y]
            # behind ix, and so behind every later interval in x: drop it
            iy.last < ix.first && continue
            queue[w] = qy
            w += 1
            # nothing more is needed for this ix, but keep compacting
            keep_scanning || continue
            if compare(ix, iy) == 0
                keep_scanning = aggregator(results, index_x, index_y)
            end
        end
        resize!(queue, w - 1)
    end
    results
end




"""
    find_intersections(::Type{Vector{Vector{T}}}, x, y) where T <: Integer

Like [`find_intersections`](@ref), but return the indices as `Vector{Vector{T}}`
instead of `Vector{Vector{Int}}`.
"""
function find_intersections(::Type{Vector{Vector{T}}}, x, y) where T <: Integer
    results = [T[] for i in 1:length(x)]
    aggregator = (r, x, y) -> (push!(r[x], y); true)
    _find_intersections(results, x, y, aggregator)
end


"""
    find_intersections(::Type{Vector{T}}, x, y) where T <: Integer

Like [`find_intersections`](@ref), but return a `Vector{T}` where the i-th
element is the number of intervals in `y` that intersect with the i-th interval in `x`.
"""
function find_intersections(::Type{Vector{T}}, x, y) where T <: Integer
    results = zeros(T, length(x))
    aggregator = (r, x, y) -> (r[x] += 1; true)
    _find_intersections(results, x, y, aggregator)
end

"""
    find_intersections(::Type{Vector{Bool}}, x, y)

Like [`find_intersections`](@ref), but return a `Vector{Bool}` where the i-th
element indicates whether the i-th interval in `x` intersects with any interval in `y`.
"""
function find_intersections(::Type{Vector{Bool}}, x, y)
    results = fill(false, length(x))
    aggregator = (r, x, y) -> (r[x] = true; false)  # one hit settles it
    _find_intersections(results, x, y, aggregator)
end



"""
    find_intersections(x, y)

Find intersections between two arrays of intervals. 

Returns an array of arrays, where the i-th element contains the indices of
intervals in `y` that intersect with the i-th interval in `x`.

The intervals in `x` and `y` do not need to be sorted. However the function
will sort them internally, so for repeated calls it is more efficient to sort 
them before calling this function.
"""
find_intersections(x, y) = find_intersections(Vector{Vector{Int}}, x, y)
