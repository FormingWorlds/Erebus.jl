using Test
using Erebus

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_darcy_flow.jl"))

const SDF = SphericalDarcyFlowValidation

@testset "Spherical Darcy Flow Validation" begin
    runs = Dict(
        (n, sf) => SDF.run_case(n; spherical_flow=sf) for n in (33, 65), sf in (false, true)
    )

    @testset "excess pressure profile against the exact solution (spherical_flow = $sf)" for sf in
                                                                                             (
        false, true
    )
        e33 = runs[(33, sf)].l2_inner
        e65 = runs[(65, sf)].l2_inner
        @test e65 < 0.08
        @test log2(e33 / e65) > 0.8
    end

    @testset "vented flux equals the discrete source (spherical_flow = $sf)" for sf in
                                                                                 (false, true)
        m = runs[(65, sf)]
        @test isapprox(m.vent_flux, m.source_flux; rtol=1.0e-8)
    end

    @testset "centre excess pressure on/off is 2/3" begin
        ratio = runs[(65, true)].dp_centre / runs[(65, false)].dp_centre
        @test isapprox(ratio, 2 / 3; atol=0.01)
    end

    @testset "matrix stays rigid" for sf in (false, true)
        m = runs[(65, sf)]
        @test m.max_solid_speed < 1.0e-12
    end
end
