```@meta
CurrentModule = GenomicCoordinates
```

# GenomicCoordinates

Documentation for [GenomicCoordinates](https://github.com/ArndtLab/GenomicCoordinates.jl).

`GenomicCoordinates` is a Julia package for working with genomic coordinates. It provides:

- [`chr2int`](@ref) to convert chromosome names to `Int` for fast comparison and sorting.
- [`GenomicPosition`](@ref) to represent a single position on a chromosome.
- [`GenomicInterval`](@ref) to represent a closed interval on a chromosome, built on top of
  [Intervals.jl](https://github.com/invenia/Intervals.jl).
- [`find_intersections`](@ref) to efficiently find overlaps between two collections of intervals.

## Installation

```julia
using Pkg
Pkg.add("GenomicCoordinates")
```

## Quick start

```@example quickstart
using GenomicCoordinates

# convert chromosome names to Int, for fast comparison and sorting
chr2int("chrX")
```

```@example quickstart
# define genomic intervals
i1 = GenomicInterval(1,   0, 200)
i2 = GenomicInterval(1, 201, 400)
i3 = GenomicInterval(1, 401, 600)

gene1 = GenomicInterval(1,  50,  70)
gene2 = GenomicInterval(1, 150, 200)
gene3 = GenomicInterval(1, 350, 400)

# find, for each interval in the first vector, the indices of overlapping
# intervals in the second vector
GenomicCoordinates.find_intersections([i1, i2, i3], [gene1, gene2, gene3])
```

Intervals also support the standard comparisons and set operations from
[Intervals.jl](https://github.com/invenia/Intervals.jl), such as `<`, `==`, `in`, and `intersect`:

```@example quickstart
using Intervals

i1 < i2
```

```@example quickstart
GenomicPosition(1, 100) in i1
```

## API reference

```@index
```

```@autodocs
Modules = [GenomicCoordinates]
```
