content = read("src/simulation/loop.jl", String)
content = replace(content, "L_3D_equiv = 2.0 * rplanet_val" => "L_3D_equiv = compute_l3d_metric(rplanet_val)")
write("src/simulation/loop.jl", content)

content2 = read("src/physics/magma_ocean_degassing.jl", String)
content2 = replace(content2, "L_3D_equiv = 2.0 * rplanet_val" => "L_3D_equiv = compute_l3d_metric(rplanet_val)")
write("src/physics/magma_ocean_degassing.jl", content2)
