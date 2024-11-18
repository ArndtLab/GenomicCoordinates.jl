module GenomicCoordinates

using Intervals

export GenomicPosition, GenomicInterval, 
    chr2int

include("types.jl")
include("intersections.jl")



end
