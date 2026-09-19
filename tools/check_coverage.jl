using Pkg
Pkg.add("Coverage")
using Coverage
cov = process_folder()
covered, total = get_summary(cov)
ratio = covered / total
println("Test coverage: ", round(ratio * 100, digits=2), "%")
if ratio < 0.90
    println(stderr, "Error: Coverage is below 90.0% floor!")
    exit(1)
end
