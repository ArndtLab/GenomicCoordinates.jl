# GenomicCoordinates

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ArndtLab.github.io/GenomicCoordinates.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ArndtLab.github.io/GenomicCoordinates.jl/dev/)
[![Build Status](https://github.com/ArndtLab/GenomicCoordinates.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArndtLab/GenomicCoordinates.jl/actions/workflows/CI.yml?query=branch%3Amain)




## Overview

`GenomicCoordinates` is a Julia package designed to handle genomic coordinates efficiently. It allows to:


- convert chromosome names to `Int` for higher performace during comparison and sorting operations
- define genomic intervals
- effeciently find overlaps in two vectors of genomic intervals.


## Usage

```julia
using GenomicCoordinates

# convert chromosome names to Int
chr2int("chrX")

# define intervals
i1 = GenomicInterval(1,   0, 200)
i2 = GenomicInterval(1, 201, 400)
i3 = GenomicInterval(1, 401, 600)

gene1 = GenomicInterval(1,  50,  70)
gene2 = GenomicInterval(1, 150, 200)
gene3 = GenomicInterval(1, 350, 400)

# find overlapping intervals
GenomicCoordinates.find_intersections([i1, i2, i3], [gene1, gene2])
# [ [1,2], [3], [] ]
```