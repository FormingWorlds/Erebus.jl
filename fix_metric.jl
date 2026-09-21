content = read("src/physics/materials_metric.jl", String)
content = replace(content, 
    "function compute_l3d_metric(rplanet::Real)\n    return (4.0 / 3.0) * rplanet\nend" => 
    "function compute_l3d_metric(rplanet::Real)\n    rplanet <= 0.0 && throw(DomainError(rplanet, \"Planet radius must be positive\"))\n    return (4.0 / 3.0) * rplanet\nend")
write("src/physics/materials_metric.jl", content)
