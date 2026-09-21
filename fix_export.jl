content = read("src/Erebus.jl", String)
content = replace(content, 
    "compute_rhofluid," => 
    "compute_rhofluid,\n        compute_l3d_metric,")
write("src/Erebus.jl", content)
