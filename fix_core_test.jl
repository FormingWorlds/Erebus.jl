content = read("test/test_core_volatile_partitioning.jl", String)
content = replace(content, 
    "println(\"Xfe_C_m[1] = \", Xfe_C_m[1])" => 
    "@test isapprox(Xfe_C_m[1], 282.78307775976026, rtol=1e-8)")
write("test/test_core_volatile_partitioning.jl", content)
