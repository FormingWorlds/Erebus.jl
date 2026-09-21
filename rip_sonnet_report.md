=== Sonnet RIP Agent 1 (Correctness) ===
All three workflow YAML files parse cleanly after the diff. Here is the review.

## Findings

**1. README.md now makes a false claim about CI coverage (confirmed, verified against the actual workflow file)**

`README.md` line 157 was changed to:

> "Continuous Integration verifies test execution on Julia 1.12 and 1.13, and enforces formatting and documentation builds on Julia 1.12 and 1.13."

I read the resulting `.github/workflows/CI.yml`. The `docs`, `format`, and `test-quality` jobs each pin a single fixed version:

```yaml
- uses: julia-actions/setup-julia@v1
  with:
    version: '1.12'
```

None of these three jobs use a matrix; they never run on 1.13. Before this diff the README correctly said "enforces formatting and documentation builds on Julia 1.10" (singular, matching the single-version jobs of that time). The diff changed the version number but also changed "Julia 1.10" to "Julia 1.12 **and 1.13**", turning a previously accurate sentence into a false one. The "test execution on Julia 1.12 and 1.13" half of the sentence is accurate (the `test` job matrix really does include both), so this is a targeted inaccuracy in the second clause only.

Failure scenario: a contributor reads the README, believes formatting/doc-build failures are caught on 1.13, and doesn't test doc builds under 1.13 before merging — the actual CI can't catch that.

**2. `nightly.yml`'s coverage step was deleted rather than fixed — plausible but not a functional regression**

I compared against the pre-diff file (`git show 5eb612d:.github/workflows/nightly.yml`) and confirmed with `yaml.safe_load` that the old "Check Coverage" step was already invalid YAML (`run:` scalar followed by a second unindented command line with no key — "could not find expected ':'"). So the nightly workflow's coverage check was already broken before this diff; it likely never ran.

This diff removes the step and sets `coverage: false` instead of fixing the syntax. I checked `test/runtests.jl`: `EREBUS_TEST_GROUP=all` (used by `CI.yml`'s dedicated `coverage` job on every push/PR to `main`) already runs both `unit_tests` and `integration_tests` together with `coverage: true` and `tools/check_coverage.jl`. So integration-test coverage is still enforced elsewhere, and removing the broken, redundant step doesn't leave a coverage gap. I'm flagging this only because the diff silently drops behavior instead of repairing it, with no note that it's already covered by a different job — worth a one-line justification if this is reviewed by someone unfamiliar with `CI.yml`.

## What I checked vs. did not check

- Verified: all three post-diff workflow YAML files parse without error; the pre-diff `nightly.yml` "Check Coverage" step does not parse; `test/runtests.jl` group semantics (`all` = `unit` + `integration`); the exact jobs/versions in the current `CI.yml`; no leftover `1.10`/`1.11` references in README.md or `docs/src/howto/installation.md`.
- Not checked: whether the package actually builds and its test suite passes under Julia 1.12/1.13 (no Julia toolchain run) — I have not verified the compat floor bump (`julia = "1"` → `"1.12"` in `Project.toml`) is technically justified by any version-gated code, only that the syntax is valid semver-caret notation.

=== Sonnet RIP Agent 2 (Completeness) ===
## Adversarial review report

**Scope reviewed**: the diff moving CI/compat from Julia 1.10/1.11 to 1.12/1.13 (`.github/workflows/CI.yml`, `Coverage.yml`, `nightly.yml`, `Project.toml`, `README.md`, `docs/src/howto/installation.md`). Cross-checked against git history (`5eb612d` → `e92eeb5` → `e36cdac`) and the current working tree.

### Finding 1 — Nightly coverage enforcement silently removed, not fixed (confirmed by execution)

The diff changes `coverage: true` → `false` in `nightly.yml` and deletes the whole "Check Coverage" step (the one calling `tools/check_coverage.jl`, which enforces the 90% floor). Nothing in the diff or the commit context explains this; it rides along inside a commit titled "Update to Julia 1.13 and remove older compatibility."

I checked why: the pre-diff "Check Coverage" step was invalid YAML. I extracted it and the full pre-diff file and parsed both with PyYAML:
```
yaml.scanner.ScannerError: while scanning a simple key ... could not find expected ':'
```
It fails on both the isolated step and the full file (`git show 5eb612d:.github/workflows/nightly.yml`). The `run:` value was written on two lines without a `run: |` block scalar, so the second line (`julia --project=tools tools/check_coverage.jl`) is not attached to `run:` at all — the file cannot be parsed as YAML. I confirmed the current `nightly.yml`, `CI.yml`, and `Coverage.yml` all parse cleanly now.

So the pre-existing nightly workflow was non-functional (GitHub Actions would reject the whole file, not just fail one step) since `5eb612d` (#75, "enforce strict coverage floor"). The diff under review does make the file valid again — but only by deleting the check outright, rather than fixing the syntax (`run: |`) and keeping coverage enforcement on the nightly integration run. The correct fix (`CI.yml`'s own "coverage" job already shows the right pattern, `run: |` + two lines) was sitting right there in the same repo.

Net effect: the nightly schedule (Julia 1.12/1.13, `EREBUS_TEST_GROUP: integration`) no longer collects or checks coverage at all. Project-wide 90% floor enforcement still exists via `CI.yml`'s separate `coverage` job (`EREBUS_TEST_GROUP: all`, runs on every push/PR), so this is not a total loss of the floor — but it is an undisclosed regression of a guarantee a very recent PR (#75) explicitly added, bundled into an unrelated version-bump commit.

### Finding 2 — Julia 1.13 is untested by most non-matrix jobs (confirmed by reading)

`CI.yml`'s `test` job matrix now covers `1.12` and `1.13`, but `docs`, `format`/lint, `test-quality`, and `Coverage.yml`'s `coverage` job are all pinned to `'1.12'` only. So JuliaFormatter checks, the test-quality linter (`check_test_quality.jl`), doc builds, and the coverage floor are never exercised under 1.13 — only the plain test run is. If 1.13 support is a real goal, this leaves a gap; if it's just "make sure tests pass," it should probably not be advertised as full 1.13 CI coverage in the README ("CI verifies test execution on Julia 1.12 and 1.13, and enforces formatting and documentation builds on Julia 1.12 and 1.13" — the README's own new wording overstates this, since formatting/docs are 1.12-only).

### Finding 3 — No demonstrated reason for the 1.12 compat floor (informational)

`Project.toml`'s `julia` compat moved from `"1"` to `"1.12"`, and 1.10/1.11 were dropped from the test matrix entirely (not just left untested-but-allowed). I grepped `src/` for `VERSION` checks, `@static`, `pkgversion` — none exist. Nothing in this diff shows a technical reason 1.10/1.11 no longer work; it reads as a policy choice ("remove older compatibility") with no changelog entry or justification, meaning downstream users pinned to 1.10/1.11 LTS get a hard `Pkg` resolution failure with no explanation of why.

### What I did and did not verify

- Verified by execution: the pre-diff nightly.yml step (and full file) fails YAML parsing; the post-diff `CI.yml`/`Coverage.yml`/`nightly.yml` all parse.
- Verified by reading: version pins across all four workflow files, `Project.toml` compat line, README/docs wording, absence of version-gated code in `src/`.
- Not verified: I did not actually run the GitHub Actions workflows (no CI execution available here), so I cannot confirm the jobs succeed end-to-end on 1.12/1.13, only that the YAML is now syntactically valid and internally consistent with the stated matrix.
- Out of scope: I did not review the many untracked scratch files in the working tree (`patch_*.jl`, `rip_*` files, etc.) — those aren't part of the diff under review.

=== Sonnet RIP Agent 3 (Design) ===
**Scope of what I checked**: I read the post-diff contents of `.github/workflows/{CI,Coverage,nightly}.yml`, `README.md`, `Project.toml`, and `test/runtests.jl`'s `EREBUS_TEST_GROUP` branches, and diffed against commit `5eb612d` (which introduced the nightly coverage step) to see its original form. I did not run any GitHub Actions workflow or a YAML parser against the pre-diff nightly.yml, so the "invalid YAML" claim on that step is based on visual inspection of indentation, not an executed parse — worth confirming with a linter before treating it as certain. I did not check whether Erebus's source actually uses any Julia-1.12-only language feature that would justify the `Project.toml` compat bump.
