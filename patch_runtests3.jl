content = read("test/runtests.jl", String)
idx = findfirst("files_to_run = String[]", content)
insert_point = idx[1]

file_check = """
all_test_files = filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(@__DIR__))
listed_tests = Set(vcat(unit_tests, integration_tests))
for f in all_test_files
    if f ∉ listed_tests && f != "test_helpers.jl" && f != "test_constants.jl"
        error("Test file \$f is not listed in runtests.jl (neither unit nor integration group).")
    end
end

"""

content = content[1:(insert_point - 1)] * file_check * content[insert_point:end]
write("test/runtests.jl", content)
