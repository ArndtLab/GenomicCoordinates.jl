



@testitem "overlap 1" begin
    using DataFrames


    function eachoverlap(intervals1, intervals2;
        chr1=:chr, first1=:first, last1=:last,
        chr2=:chr, first2=:first, last2=:last
    )

        gi1 = GenomicInterval.(intervals1[:, chr1], intervals1[:, first1], intervals1[:, last1])
        gi2 = GenomicInterval.(intervals2[:, chr2], intervals2[:, first2], intervals2[:, last2])
        r = GenomicCoordinates.find_intersections(gi1, gi2)
        mapreduce(vcat, enumerate(r)) do (i, x)
            map(x) do k
                (i, k)
            end
        end
    end



    df = DataFrame(chr=["1", "1", "2", "2", "2"],
        firstbp=[1, 10, 3, 4, 9],
        lastbp=[3, 13, 10, 4, 18])

    i = 0

    @show eachoverlap(df, df, first1=:firstbp, first2=:firstbp, last1=:lastbp, last2=:lastbp)

    for (a, b) in eachoverlap(df, df,
        first1=:firstbp, first2=:firstbp,
        last1=:lastbp, last2=:lastbp)

        global i += 1
    end
    @test i == 9

    i = 0
    for (a, b) in eachoverlap(df, df,
        first1=2, first2=2,
        last1=3, last2=3)

        global i += 1
    end
    @test i == 9

    # ----------------


    L = 100
    N = 1000
    intervals = DataFrame(chr=String[], first=Int[], last=Int[])
    for i in 1:N
        push!(intervals, ("A", i * L + 1, i * L + L))
    end

    for Le in (1, 2, 10, 40, 100)
        events = DataFrame(chr=String[], first=Int[], last=Int[])
        for i in 1:(N+3)*L
            push!(events, ("A", i, i + Le - 1))
        end

        intervals[!, :nevents] .= 0
        @time for (a, b) in eachoverlap(intervals, events)
            intervals[a, :nevents] += 1
        end

        @test all(intervals.nevents .== L + Le - 1)

        intervals[!, :nevents] .= 0
        @time for (a, b) in eachoverlap(events, intervals)
            intervals[b, :nevents] += 1
        end

        @test all(intervals.nevents .== L + Le - 1)
    end

    # ----------------



    Na = 1000
    dfa = DataFrame(chr=String[], first=Int[], last=Int[])
    for i in 1:Na
        c = rand(["1", "2", "3"])
        f = rand(1:10)
        l = f + rand(0:10)
        push!(dfa, (c, f, l))
    end

    Nb = 1000
    dfb = DataFrame(chr=String[], first=Int[], last=Int[])
    for i in 1:Nb
        c = rand(["1", "2", "3"])
        f = rand(1:10)
        l = f + rand(0:10)
        push!(dfb, (c, f, l))
    end

    sort!(dfa, [:chr, :first])
    sort!(dfb, [:chr, :first])

    c = length([a for (a, b) in eachoverlap(dfb, dfa)])
    @test c == length([a for (a, b) in eachoverlap(dfb, dfa)])

    sort!(dfa, [:chr, :first, :last])
    @test c == length([a for (a, b) in eachoverlap(dfb, dfa)])
    @test c == length([a for (a, b) in eachoverlap(dfa, dfb)])

    sort!(dfb, [:chr, :first, :last])
    @test c == length([a for (a, b) in eachoverlap(dfb, dfa)])
    @test c == length([a for (a, b) in eachoverlap(dfa, dfb)])


end

