content = read("test/runtests.jl", String)
old_check = """
all_test_files = filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(@__DIR__))
listed_tests = Set(vcat(unit_tests, integration_tests))
for f in all_test_files
    if f ∉ listed_tests && f != "test_helpers.jl" && f != "test_constants.jl"
        error("Test file \$f is not listed in runtests.jl (neither unit nor integration group).")
    end
end
"""
new_check = """
all_test_files = filter(f -> endswith(f, ".jl") && f != "runtests.jl", readdir(@__DIR__))
listed_tests = Set(vcat(unit_tests, integration_tests))
for f in all_test_files
    if f ∉ listed_tests && f != "test_helpers.jl" && f != "test_constants.jl" && f != "mpi_worker_tests.jl"
        error("Test file \$f is not listed in runtests.jl (neither unit nor integration group).")
    end
end
"""
content = replace(content, old_check => new_check)
write("test/runtests.jl", content)
