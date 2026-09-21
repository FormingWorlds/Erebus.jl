content = read("test/runtests.jl", String)

# Remove the previous file_check if it exists
if occursin("all_test_files", content)
    lines = split(content, '\n')
    new_lines = []
    skip = false
    for l in lines
        if startswith(l, "all_test_files =")
            skip = true
        end
        if skip && startswith(l, "end") && length(l) == 3
            skip = false
            continue
        end
        if !skip
            push!(new_lines, l)
        end
    end
    # Wait, the previous block has an 'end' for the if and for loop.
    # It's easier to just recreate the file from scratch or use regex.
end
