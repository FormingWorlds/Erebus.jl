content = read("src/physics/magma_ocean_degassing.jl", String)
content = replace(content, 
    "is_degassing_zone = (r_sq >= r_degas_sq) && (F_curr >= F_thresh || F_curr > 0.01)" => 
    "is_degassing_zone = (r_sq >= r_degas_sq) && (F_curr >= F_thresh || F_curr > 0.01) && F_curr > 0.0")
write("src/physics/magma_ocean_degassing.jl", content)
