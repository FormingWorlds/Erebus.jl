content = read(".github/workflows/CI.yml", String)

test_job_code = """
  test:
    name: Tests - Julia \${{ matrix.version }} - \${{ matrix.os }} - \${{ matrix.arch }} - \${{ github.event_name }}
    runs-on: \${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        version:
          - '1.10'
          - '1.11'
        os:
          - ubuntu-latest
        arch:
          - x64
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v1
        with:
          version: \${{ matrix.version }}
          arch: \${{ matrix.arch }}
      - uses: julia-actions/cache@v1
        with:
          cache-registries: "true"
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-runtest@v1
        with:
          coverage: false
"""

test_job_replacement = """
  test:
    name: Tests - Julia \${{ matrix.version }} - \${{ matrix.os }} - \${{ matrix.arch }} - \${{ github.event_name }}
    runs-on: \${{ matrix.os }}
    timeout-minutes: 10
    strategy:
      fail-fast: false
      matrix:
        version:
          - '1.10'
          - '1.11'
        os:
          - ubuntu-latest
        arch:
          - x64
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v1
        with:
          version: \${{ matrix.version }}
          arch: \${{ matrix.arch }}
      - uses: julia-actions/cache@v1
        with:
          cache-registries: "true"
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-runtest@v1
        with:
          coverage: true
        env:
          EREBUS_TEST_GROUP: unit
      - name: Check Coverage
        run: julia tools/check_coverage.jl
"""
content = replace(content, test_job_code => test_job_replacement)
write(".github/workflows/CI.yml", content)
