using TestItems
using TestItemRunner
@run_package_tests verbose=true

  

@testitem "chr2int" begin
    @test chr2int("X") == 23
    @test chr2int("chrX") == 23
    @test chr2int("Y") == 24
    @test chr2int("1") == 1
    @test chr2int("chr1") == 1
    @test chr2int("chr01") == 1
    @test chr2int("22") == 22
    @test chr2int("M") == 25
    @test chr2int("MT") == 25
    
    @test chr2int("someplasmid") > 25
    @test chr2int("xxxMTxxx"[4:5]) == 25
    
    @test chr2int(1) == 1
    @test chr2int(23) == 23
    @test chr2int("ML679985.1") != chr2int("VOCO01011341.1")
end


@testitem "GenomicInterval" begin
    gene1 = GenomicInterval(1, 100, 200)
    @show gene1

    gene2 = GenomicInterval(1, 150, 150)
    @show gene2

    gene3 = GenomicInterval(1, 150)
    @show gene3
end

@testitem "Test Chr Int" begin
    using Intervals
    gene1 = GenomicInterval(1, 100, 200)
    gene2 = GenomicInterval(1, 150, 250)
    gene3 = GenomicInterval(1, 100, 200)
    gene4 = GenomicInterval(2, 100, 200)

    @test gene1 < gene2
    @test gene1 == gene3
    @test !(gene1 == gene2)
    @test gene1 < gene4
    @test gene2 < gene4
    @test isclosed(gene1)
    @test !isopen(gene1)
    @test Intervals.isbounded(gene1)

    @test !in(GenomicPosition(1,99), gene1)
    @test GenomicPosition(1,100) in gene1
    @test GenomicPosition(1,133) in gene1
    @test GenomicPosition(1,200) in gene1
    @test !in(GenomicPosition(1,201), gene1)

    
    @test GenomicCoordinates.segment_length(intersect(gene1, gene2)) == 51
    @show intersect(gene1, gene4)
    @test GenomicCoordinates.segment_length(intersect(gene1, gene4)) == -1
end





@testitem "Test Chr String" begin
    using Intervals
    gene1 = GenomicInterval("chr1", 100, 200)
    gene2 = GenomicInterval("chr1", 150, 250)
    gene3 = GenomicInterval("chr1", 100, 200)
    gene4 = GenomicInterval("chr2", 100, 200)

    @test gene1 < gene2
    @test gene1 == gene3
    @test !(gene1 == gene2)
    @test gene1 < gene4
    @test gene2 < gene4
end



@testitem "Compare Intervals Overlap" begin
    using Intervals
    @test GenomicCoordinates.compare_for_overlap(Interval(1, 10), Interval(11, 20)) == -1
    @test GenomicCoordinates.compare_for_overlap(Interval(1, 10), Interval(0, 20)) == 0
    @test GenomicCoordinates.compare_for_overlap(Interval(1, 10), Interval(10, 20)) == 0
    @test GenomicCoordinates.compare_for_overlap(Interval(101, 110), Interval(11, 20)) == 1
end


@testitem "Compare Intervals Inclusion" begin
    using Intervals

    @test GenomicCoordinates.compare_for_inclusion_1in2(Interval(1, 10), Interval(1, 20)) == 0
    @test GenomicCoordinates.compare_for_inclusion_1in2(Interval(1, 10), Interval(1, 10)) == 0
    @test GenomicCoordinates.compare_for_inclusion_1in2(Interval(0, 10), Interval(1, 20)) == -1
    @test GenomicCoordinates.compare_for_inclusion_1in2(Interval(2, 10), Interval(1, 20)) == 0

    r = map(1:20) do i
        GenomicCoordinates.compare_for_inclusion_1in2(Interval(i, i+2), Interval(10, 15))
    end
    @test issorted(r)
    @test sum(r .== 0) == 4
    @test length(unique(r)) == 5
    @test extrema(r) == (-2, 2)


    @test GenomicCoordinates.compare_for_inclusion_2in1(Interval(1, 10), Interval(1, 20)) == -1
    @test GenomicCoordinates.compare_for_inclusion_2in1(Interval(0, 10), Interval(1, 20)) == -1
    @test GenomicCoordinates.compare_for_inclusion_2in1(Interval(1, 10), Interval(1, 10)) == 0
    @test GenomicCoordinates.compare_for_inclusion_2in1(Interval(1, 10), Interval(1, 9)) == 0
    @test GenomicCoordinates.compare_for_inclusion_2in1(Interval(1, 10), Interval(2, 10)) == 0

    r = map(1:20) do i
        GenomicCoordinates.compare_for_inclusion_2in1(Interval(i, i+2), Interval(10, 15))
    end
    @test issorted(r)
    @test sum(r .== 0) == 0
    @test length(unique(r)) == 4
    @test extrema(r) == (-2, 2)
end



@testitem "Compare and Sort Intervals" begin
    i0 = GenomicInterval(2, 100, 200)
    i1 = GenomicInterval(1, 100, 200)
    i2 = GenomicInterval(1, 150, 250)
    i3 = GenomicInterval(1, 100, 150)

    @test i1 < i2
    @test i3 < i1
    @test i3 < i0
    @test i2 < i0

    v = [i0, i1, i2, i3]
    @test sort(v) == [i3, i1, i2, i0]
    sort!(v)
    @test v == [i3, i1, i2, i0]
end



@testitem "Overlaps" begin
    using Intervals
    
    v1 = [ 
        GenomicInterval(1, 350, 450),
        GenomicInterval(1, 100, 200), 
        GenomicInterval(1, 150, 250), 
        GenomicInterval(1, 300, 400),
        GenomicInterval(1, 300, 350),
        ]
    v2 = [
        GenomicInterval(1, 50, 150),
        GenomicInterval(1, 200, 300),
        GenomicInterval(2, 200, 300),
        GenomicInterval(1, 350, 450)
        ]




    r1 = Intervals.find_intersections(v1, v2)
    r2 = GenomicCoordinates.find_intersections(v1, v2)
    @test r1 == r2
    r3 = GenomicCoordinates.find_intersections(Vector{Int}, v1, v2)
    @test r3 == map(length, r2)
    r4 = GenomicCoordinates.find_intersections(Vector{Bool}, v1, v2)
    @test r4 == map(v -> length(v) > 0, r2)
    @test r4 == map(>(0), r3)


    v1 = reverse(v1)

    r1 = Intervals.find_intersections(v1, v2)
    r2 = GenomicCoordinates.find_intersections(v1, v2)
    @test r1 == r2
    r3 = GenomicCoordinates.find_intersections(Vector{Int}, v1, v2)
    @test r3 == map(length, r2)
    r4 = GenomicCoordinates.find_intersections(Vector{Bool}, v1, v2)
    @test r4 == map(v -> length(v) > 0, r2)
    @test r4 == map(>(0), r3)


    
end

include("eachoverlap.jl")