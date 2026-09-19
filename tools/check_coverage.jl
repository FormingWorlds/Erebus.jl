using Coverage

cov = process_folder()
covered, total = get_summary(cov)

if total == 0
    println(stderr, "Error: No coverage data collected! (total lines = 0)")
    exit(1)
end

ratio = covered / total
println("Test coverage: ", round(ratio * 100, digits=2), "% (", covered, "/", total, " lines)")

if ratio < 0.90 || isnan(ratio)
    println(stderr, "Error: Coverage is below 90.0% floor!")
    exit(1)
end
