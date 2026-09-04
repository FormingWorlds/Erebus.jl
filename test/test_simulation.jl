@testset "Simulation" begin
    @testset "setup_dynamic_simulation_parameters()" begin
        # set up dynamic simulation parameters
        (timestep, dt, timesum, marknum, hrsolidm, hrfluidm, YERRNOD) = Erebus.setup_dynamic_simulation_parameters()
        # verification & test
        @test timestep == start_step
        @test dt == dt_longest
        @test timesum == start_time
        @test marknum == start_marknum
        @test hrsolidm == start_hrsolidm
        @test hrfluidm == start_hrfluidm
        @test YERRNOD == zeros(nplast)
    end # testset "setup_dynamic_simulation_parameters()"

    @testset "s_to_Ma()" begin
        for _ in 1:1:10
            s = rand(rgen, Int)
            Ma = Erebus.s_to_Ma(s)
            @test Ma ≈ s * 1e-6 / yearlength rtol=1e-9
        end
    end # testset "s_to_Ma()"
end
