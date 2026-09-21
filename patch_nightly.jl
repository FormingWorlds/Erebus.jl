content = read(".github/workflows/nightly.yml", String)
content = replace(
    content,
    "run: julia tools/check_coverage.jl" => "run: julia --project=tools tools/check_coverage.jl",
)
write(".github/workflows/nightly.yml", content)
