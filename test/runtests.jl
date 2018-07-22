using NMR
using Base.Test

ϵ = 0.001

# Phenyl boronic acid in CDCl3 (300 MHz)
dataset_path = joinpath(@__DIR__, "data", "PhB(OH)2", "1")

# write your own tests here
@testset "Read Bruker" begin
    s = Spectrum(dataset_path, 1)
    @test s.default_proc == 1
    @test length(s) == length(s[1])  == s["SI"] == 32768
    @test s["TD"] == 18486
    @test typeof(s[:]) == typeof(s[1][:]) <: Array{T,1} where T <: Number
    @test limits(s)[1] ≈ -0.4997 atol=ϵ
    @test limits(s)[2] ≈ 13.4996 atol=ϵ
    @test all( [1.2, 3.5, 7.9, 12.9] .∈ Ref(s) )
    @test all( [14.2, 30.9, -2.0, -0.7] .∉ Ref(s) )
end

@testset "Unit conversions" begin
    s = Spectrum(dataset_path, 1)
    sf = s["SF"]
    bf = s["BF1"]
    @test NMR.sr(sf, bf) ≈ 0.0
    @test -2 < ppmtoindex(s, 6.4502) - 16502 < 2
    let rng = ppmtoindex(s, limits(s))
        @test 1 ∈ rng
        @test s["SI"] - 1 ∈ rng
    end
end