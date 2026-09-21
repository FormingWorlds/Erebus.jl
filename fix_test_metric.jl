content = read("test/test_dehydration_darcy_coupling.jl", String)
content = replace(content, "Erebus.Physics.compute_l3d_metric" => "Erebus.compute_l3d_metric")
write("test/test_dehydration_darcy_coupling.jl", content)
