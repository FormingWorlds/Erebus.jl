content = read(".github/workflows/CI.yml", String)
content = replace(
    content,
    "julia --project=tools tools/check_coverage.jl" => "julia --project=tools -e 'using Pkg; Pkg.instantiate()'\n          julia --project=tools tools/check_coverage.jl",
)
content = replace(
    content,
    "runs-on: ubuntu-latest\n    steps:\n      - uses: actions/checkout@v4" => "runs-on: ubuntu-latest\n    timeout-minutes: 45\n    steps:\n      - uses: actions/checkout@v4",
)
write(".github/workflows/CI.yml", content)

content2 = read(".github/workflows/nightly.yml", String)
content2 = replace(
    content2,
    "julia --project=tools tools/check_coverage.jl" => "julia --project=tools -e 'using Pkg; Pkg.instantiate()'\n        julia --project=tools tools/check_coverage.jl",
)
content2 = replace(
    content2,
    "runs-on: \${{ matrix.os }}\n    strategy:" => "runs-on: \${{ matrix.os }}\n    timeout-minutes: 60\n    strategy:",
)
write(".github/workflows/nightly.yml", content2)
