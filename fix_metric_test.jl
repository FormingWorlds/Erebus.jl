content = read("test/test_dehydration_darcy_coupling.jl", String)
content = replace(content, 
    "@testset \"ReactionConfig Schema & Serialization\" begin" => 
    """
    @testset "Metric L3D Computation" begin
        # Metric testing: L3D = (4/3)*R
        @test isapprox(Erebus.Physics.compute_l3d_metric(3000.0), 4000.0, rtol=1e-10)
        @test_throws DomainError Erebus.Physics.compute_l3d_metric(0.0)
        @test_throws DomainError Erebus.Physics.compute_l3d_metric(-100.0)
    end

    @testset "ReactionConfig Schema & Serialization" begin""")
write("test/test_dehydration_darcy_coupling.jl", content)
