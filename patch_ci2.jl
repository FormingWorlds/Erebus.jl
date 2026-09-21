content = read(".github/workflows/CI.yml", String)

test_job_code = """
      - uses: julia-actions/julia-runtest@v1
        with:
          coverage: true
        env:
          EREBUS_TEST_GROUP: unit
      - name: Check Coverage
        run: julia tools/check_coverage.jl
"""

test_job_replacement = """
      - uses: julia-actions/julia-runtest@v1
        with:
          coverage: false
        env:
          EREBUS_TEST_GROUP: unit

  coverage:
    name: Coverage
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v1
        with:
          version: '1.10'
      - uses: julia-actions/cache@v1
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-runtest@v1
        with:
          coverage: true
        env:
          EREBUS_TEST_GROUP: all
      - name: Check Coverage
        run: |
          julia --project=tools tools/check_coverage.jl
"""
content = replace(content, test_job_code => test_job_replacement)
write(".github/workflows/CI.yml", content)
