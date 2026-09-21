content = read("test/runtests.jl", String)
old_code = """
if test_group == "all" || test_group == "unit"
    append!(files_to_run, unit_tests)
end
if test_group == "all" || test_group == "integration"
    append!(files_to_run, integration_tests)
end
"""
new_code = """
if test_group == "all"
    append!(files_to_run, unit_tests)
    append!(files_to_run, integration_tests)
elseif test_group == "unit"
    append!(files_to_run, unit_tests)
elseif test_group == "integration"
    append!(files_to_run, integration_tests)
else
    error("Unknown EREBUS_TEST_GROUP: \$test_group")
end

if isempty(files_to_run)
    error("No tests found to run!")
end
"""
content = replace(content, old_code => new_code)
write("test/runtests.jl", content)
