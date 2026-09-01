
"""
    GenomicPosition{C,P}

A genomic position.
"""
struct GenomicPosition{C,P<:Integer}
    chr::C
    pos::P
end



Base.show(io::IO, p::GenomicPosition) = print(io, "GenomicPosition(", p.chr, ", ", p.pos, ")")

"""
    chr(p::GenomicPosition)

Return the chromosome of the genomic position `p`.
"""
chr(p::GenomicPosition) = p.chr

"""
    pos(p::GenomicPosition)

Return the position (offset on the chromosome) of the genomic position `p`.
"""
pos(p::GenomicPosition) = p.pos

function Base.isless(a::GenomicPosition, b::GenomicPosition)
    (a.chr, a.pos) < (b.chr, b.pos)
end


function Base.zero(::Type{GenomicPosition{C, P}}) where {C, P}
    GenomicPosition(zero(C), zero(P))
end


"""
    GenomicInterval{C,P}

A genomic interval.
"""
const GenomicInterval{C,P} = Intervals.Interval{GenomicCoordinates.GenomicPosition{C,P},Closed,Closed}


"""
    GenomicInterval(chr, first::P, last::P)

Create a `GenomicInterval` with chromosome `chr`, start position `first`, and end position `last`.
"""
GenomicInterval(chr, first::P, last::P) where {P<:Integer} =
    Interval{GenomicPosition{typeof(chr),typeof(first)},Closed,Closed}(
        GenomicPosition(chr, first),
        GenomicPosition(chr, last)
    )

"""
    GenomicInterval(chr, pos::P)

Create a single-position `GenomicInterval` on chromosome `chr` at position `pos`.
"""
GenomicInterval(chr, pos::P) where {P<:Integer} =
    Interval{GenomicPosition{typeof(chr),typeof(pos)},Closed,Closed}(
        GenomicPosition(chr, pos),
        GenomicPosition(chr, pos)
    )

Base.show(io::IO, iv::GenomicInterval) =
    print(io, "GenomicInterval(", chr(iv.first), ": ",
        pos(iv.first), " - ", pos(iv.last), ")")




function Base.isless(a::GenomicInterval{C, P}, b::GenomicInterval{C, P}) where {C, P}
    (a.first, a.last) < (b.first, b.last)
end


"""
    chr2int(chr)::Int

Convert a chromosome name to an integer.
"""
function chr2int(chr::AbstractString)::Int64
    # SubString avoids copying the name just to strip the "chr" prefix
    name = startswith(chr, "chr") ? SubString(chr, 4) : SubString(chr, 1)
    if !isempty(name) && isdigit(@inbounds name[1])
        n = tryparse(Int64, name)
        n === nothing || return n
    elseif name == "X"
        return 23
    elseif name == "Y"
        return 24
    elseif name == "M" || name == "MT"
        return 25
    end
    2^8 + hash(name) % 2^(64 - 8 - 1)  # Use a hash for non-numeric chromosomes
    # to ensure uniqueness and avoid collisions with numeric chromosomes
    # This is a simple way to convert non-numeric chromosome names to integers.
end

function chr2int(chr::Int)::Int
    return chr
end




"""
    segment_length(a::GenomicInterval{C,P})

Return the length of the genomic segment `a`.
"""
function segment_length(a::Interval{GenomicPosition{C,P},L,R}) where {C,P,L,R}
    chr(a.first) == chr(a.last) ?
    pos(a.last) - pos(a.first) - 1 + (L === Closed ? 1 : 0) + (R === Closed ? 1 : 0) :
    zero(P)
end
