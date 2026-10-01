# Pinning tests for canonical physical constants in Erebus.jl.
# Pinned to CODATA 2018 exact definitions and IUPAC 2013 standard atomic weights.

using Test
using Erebus

@testset "Canonical Constants Pinning (CODATA 2018 / IUPAC 2013)" begin
    @testset "Universal Gas Constant R_GAS" begin
        # CODATA 2018 exact definition: N_A * k_B = 8.31446261815324 J/(mol*K)
        codata_r_gas = 8.31446261815324
        @test isapprox(Erebus.R_GAS, codata_r_gas; atol=0.0, rtol=1e-15)
        @test isapprox(
            Erebus.R_GAS,
            Erebus.AVOGADRO_CONSTANT * Erebus.BOLTZMANN_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test Erebus.RG === Erebus.R_GAS
        @test Erebus.K_BOLTZMANN === Erebus.BOLTZMANN_CONSTANT
    end

    @testset "Julian Year Length SECONDS_PER_YEAR" begin
        # Standard Julian year: 365.25 days * 86,400 s/day = 31,557,600.0 s
        iau_seconds_per_year = 31_557_600.0
        @test isapprox(Erebus.SECONDS_PER_YEAR, iau_seconds_per_year; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SECONDS_PER_YEAR, 365.25 * 86_400.0; atol=0.0, rtol=1e-15)
        @test Erebus.SEC_PER_YEAR === Erebus.SECONDS_PER_YEAR
    end

    @testset "Astronomical and Particle Constants" begin
        @test isapprox(Erebus.AU_METERS, 1.495978707e11; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_SUN_KG, 1.98847e30; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.BOLTZMANN_CONSTANT, 1.380649e-23; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.AVOGADRO_CONSTANT, 6.02214076e23; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.ATOMIC_MASS_UNIT, 1.66053906660e-27; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_PROTON, 1.67262192e-27; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.GRAVITATIONAL_CONSTANT, 6.67430e-11; atol=0.0, rtol=1e-15)
        @test Erebus.G_GRAV === Erebus.GRAVITATIONAL_CONSTANT
    end

    @testset "Elemental Atomic Weights [kg/mol] (IUPAC 2013)" begin
        @test isapprox(Erebus.M_H, 0.001008; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_C, 0.012011; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_N, 0.014007; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_O, 0.015999; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_S, 0.03206; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_Fe, 0.055845; atol=0.0, rtol=1e-15)
    end

    @testset "Compound Species Stoichiometric Invariants" begin
        # Invariant relations: M_compound == sum(nu_i * M_i)
        @test isapprox(Erebus.M_H2O, 2.0 * Erebus.M_H + Erebus.M_O; atol=0.0, rtol=1e-15)
        @test Erebus.MH₂O === Erebus.M_H2O
        @test Erebus.MH2O === Erebus.M_H2O

        @test isapprox(Erebus.M_CO, Erebus.M_C + Erebus.M_O; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_CO2, Erebus.M_C + 2.0 * Erebus.M_O; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_CH4, Erebus.M_C + 4.0 * Erebus.M_H; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_N2, 2.0 * Erebus.M_N; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_NH3, Erebus.M_N + 3.0 * Erebus.M_H; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_H2S, 2.0 * Erebus.M_H + Erebus.M_S; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_S2, 2.0 * Erebus.M_S; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_SO2, Erebus.M_S + 2.0 * Erebus.M_O; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.M_FeS, Erebus.M_Fe + Erebus.M_S; atol=0.0, rtol=1e-15)
    end

    @testset "Species AMU Dictionary Consistency" begin
        # Derived from M_* * 1000.0 to prevent dual tables
        @test isapprox(Erebus.SPECIES_AMU[:H], Erebus.M_H * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:C], Erebus.M_C * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:N], Erebus.M_N * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:O], Erebus.M_O * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:S], Erebus.M_S * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:Fe], Erebus.M_Fe * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:H2O], Erebus.M_H2O * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:CO], Erebus.M_CO * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:CO2], Erebus.M_CO2 * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:CH4], Erebus.M_CH4 * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:N2], Erebus.M_N2 * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:NH3], Erebus.M_NH3 * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:H2S], Erebus.M_H2S * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:SO2], Erebus.M_SO2 * 1000.0; atol=0.0, rtol=1e-15)
        @test isapprox(Erebus.SPECIES_AMU[:S2], Erebus.M_S2 * 1000.0; atol=0.0, rtol=1e-15)
    end

    @testset "Atmospheric Escape Species Masses" begin
        @test isapprox(
            Erebus.MASS_H2O_KG,
            Erebus.M_H2O / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_H2_KG,
            (2.0 * Erebus.M_H) / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_N2_KG,
            Erebus.M_N2 / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_NH3_KG,
            Erebus.M_NH3 / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_CO_KG,
            Erebus.M_CO / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_CO2_KG,
            Erebus.M_CO2 / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_CH4_KG,
            Erebus.M_CH4 / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_H2S_KG,
            Erebus.M_H2S / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_S2_KG,
            Erebus.M_S2 / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
        @test isapprox(
            Erebus.MASS_SO2_KG,
            Erebus.M_SO2 / Erebus.AVOGADRO_CONSTANT;
            atol=0.0,
            rtol=1e-15,
        )
    end
end
