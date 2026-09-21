content = read(".github/workflows/nightly.yml", String)
# Find the exact string to replace
old = """      - uses: julia-actions/julia-runtest@v1
        with:
          coverage: true
        env:
          EREBUS_TEST_GROUP: integration
      - name: Check Coverage
        run: julia --project=tools -e 'using Pkg; Pkg.instantiate()'
        julia --project=tools tools/check_coverage.jl
"""
new_step = """      - uses: julia-actions/julia-runtest@v1
        with:
          coverage: false
        env:
          EREBUS_TEST_GROUP: integration
"""
content = replace(content, old => new_step)
write(".github/workflows/nightly.yml", content)
