"""
    GenomicCoordinates

A package for working with genomic coordinates. It provides:

- [`GenomicPosition`](@ref), a single position on a chromosome.
- [`GenomicInterval`](@ref), a closed interval on a chromosome, built on top of
  [Intervals.jl](https://github.com/invenia/Intervals.jl).
- [`chr2int`](@ref), for converting chromosome names to `Int` for fast comparison and sorting.
- [`find_intersections`](@ref), for efficiently finding overlaps between two collections of intervals.
"""
module GenomicCoordinates

using Intervals

export GenomicPosition, GenomicInterval, 
    chr2int

include("types.jl")
include("intersections.jl")



end
