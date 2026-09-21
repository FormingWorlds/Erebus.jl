=== Opus RIP Agent 1 (Correctness) ===
I have verified the substantive points directly. Here is my adversarial review.

## Adversarial review: Julia 1.10/1.11 → 1.12/1.13 version bump

**Scope.** The diff changes CI matrices, pinned tool versions, `Project.toml` compat, and two doc files. It also removes the nightly coverage step. This is a mechanical version bump plus one behavioral change (nightly coverage). No `src/` logic changes.

### What I verified (facts, not reasoning)

1. All three current workflow files parse: `nightly.yml`, `CI.yml`, `Coverage.yml` all load without error (Python `yaml.safe_load`).
2. The removed nightly "Check Coverage" step was invalid YAML. I reconstructed it from the diff and the parser failed with `could not find expected ':'` at the bare `julia --project=tools tools/check_coverage.jl` line. That line sat at the same indent as `run:`, so it was a scalar where a mapping key was expected.
3. No leftover `1.10`/`1.11` references exist in `README.md` or `docs/src/howto/installation.md`. The remaining grep hits are scratch files (`patch_*.jl`, `rip_*.txt`, `pr_diff.txt`) and numeric coincidences (solubility constant `1.1008`, mass `1.1124e21`).
4. `[compat] julia = "1.12"` resolves to `>=1.12, <2` under Julia semver, so the `1.13` matrix entry is allowed. Compat and matrix are consistent.
5. The 90% coverage floor still exists in `CI.yml`'s `coverage` job (`EREBUS_TEST_GROUP: all`, then `check_coverage.jl`). `all` = unit + integration per `test/runtests.jl:81-90`. The diff does not touch that job.

### Findings, most severe first

**F1. Julia 1.13 availability is unverified and is the highest-impact risk (must verify).**
Every 1.13 matrix entry in `CI.yml` and `nightly.yml` calls `setup-julia@v1` with `version: '1.13'`. If 1.13 is not a released stable version that the action can resolve, every 1.13 job fails at setup. The pinned `1.12` jobs (coverage, docs, format, test-quality) have the same dependency on 1.12 being resolvable. I cannot confirm the Julia release state from this sandbox, and my knowledge cutoff is January 2026, so I do not know whether 1.13 reached general availability by 2026-09-20. Action: confirm 1.13 is GA (not only a release candidate) before you trust this CI, for example run `setup-julia` once or check the Julia release list. Confidence: this is a real failure mode; I did not verify whether it triggers.

**F2. Dependency support on 1.12/1.13 is not checked by the diff (verify by resolving).**
The bump raises the floor to 1.12 but does not prove the `[deps]` resolve on 1.12 or 1.13. Several deps carry tight caps: `AMD = "0.5.0 - 0.5.3"`, `AppleAccelerate = "0.4, 0.7"`, `Pardiso = "0.5, 1"`, `ExtendableSparse = "1, 2"`. If any of these caps its own `julia` compat below 1.13, the resolver drops it or fails. Action: run `Pkg.resolve()` / `Pkg.instantiate()` on both 1.12 and 1.13 before merge. I did not run this; the sandbox has no matching Julia toolchain set up for it. Confidence: unverified risk, not a confirmed bug.

**F3. Nightly no longer produces any coverage signal (low severity; the "regression" reading is wrong).**
The diff sets `coverage: false` and deletes the check step, so the nightly schedule collects no coverage. This is not a real loss of a working guarantee, for two verified reasons: (a) the old step was invalid YAML (F-fact 2), so the nightly workflow could not run successfully before this change; (b) the old check ran on `EREBUS_TEST_GROUP: integration` (4 files) but `check_coverage.jl` measures all of `src/` and demands 90%, so it could never pass. The floor stays enforced in `CI.yml`'s `coverage` job. Net: this is a cleanup. The one true consequence is that nightly now gives no coverage number, which is the correct outcome given the above.

**F4. Consistency checks pass; no bug found in the mechanical bump.**
Docs, README, matrices, pinned versions, and compat are internally consistent. Dropping 1.10 (an LTS line) is a policy choice, not a defect; the matrix still exercises the compat lower bound (1.12). I found no correctness or logic error in the changed lines themselves.

### Out of scope but worth noting
The working tree holds many untracked scratch files (`patch_*.jl`, `fix_*.jl`, `rip_*.txt`, `pr_*.txt`, `run_*.sh`). They are not part of this diff. Confirm they stay out of any commit.

### Boundary of what I checked
I verified YAML validity, text consistency, and compat/matrix logic directly. I did NOT verify: that Julia 1.12 and 1.13 are resolvable by the CI action (F1), that all package dependencies resolve on 1.12/1.13 (F2), or that the test suite passes on these versions. Those need a real toolchain run, which this read-only sandbox cannot do.

=== Opus RIP Agent 2 (Completeness) ===
I have enough to report. Network verification of 1.13 was denied, so I flag that as an open item rather than assert it.

## Adversarial review report: Julia 1.12/1.13 version bump

**Scope reviewed**: the diff only. It moves CI, docs, and compat from Julia 1.10/1.11 to 1.12/1.13 across `.github/workflows/CI.yml`, `Coverage.yml`, `nightly.yml`, `Project.toml`, `README.md`, and `docs/src/howto/installation.md`. I ignored the many untracked `patch_*.jl` and `rip_*.txt` files in the tree; they are scratch artifacts from earlier runs, not part of this diff.

### What I verified directly
- Local Julia is 1.12.7, so 1.12 is a real, resolvable version.
- `test/runtests.jl:14` defaults `EREBUS_TEST_GROUP` to `"all"`, so `Coverage.yml` (no group set) still runs the full suite for coverage.
- No stale `1.10` or `1.11` string remains in any git-tracked CI or docs file. The only `1.1*` hits in tracked files are physics constants and untracked scratch files.
- Version strings are internally consistent: `Project.toml` floor `1.12`; README states tests on 1.12 and 1.13 and docs/format on 1.12; installation says 1.12+. These match the workflows.

### Findings

**1. [Medium] The compat floor `julia = "1.12"` drops the LTS and 1.11 with no stated reason.**
`Project.toml:50` changes `julia = "1"` to `julia = "1.12"`. This excludes Julia 1.10 (the current LTS) and 1.11 from install resolution, not just from CI. The diff shows no 1.12-only language feature that forces this floor. If the code does not need a 1.12 feature, this floor is over-restrictive and blocks every LTS-channel user from installing the package. This is a compatibility policy decision. Confirm it is intended, and if so, record the reason. The package `version` is still `1.0.0-DEV`, so the break is pre-release, but the exclusion is real.

**Verify by**: run the suite on Julia 1.10/1.11, or grep `src/` for 1.12-only syntax. If nothing needs 1.12, lower the floor to keep LTS support.

**2. [Verify, high impact if wrong] Julia 1.13 is unverified as a released, resolvable version.**
`CI.yml:23`, `nightly.yml:16`, and the README claim Julia 1.13. If `julia-actions/setup-julia@v1` cannot resolve `'1.13'`, every `test` matrix job in both workflows fails at the setup step, and the README points users at a version they cannot install. I confirmed 1.12 exists (installed locally). I could NOT confirm 1.13: the network fetch of the Julia versions list was denied by the sandbox. Given today is 2026-09-20 and 1.12 is already released, 1.13 is plausible but unproven.

**Verify by**: check that `1.13.x` exists in `https://julialang-s3.julialang.org/bin/versions.json`, or watch the first CI run resolve the matrix.

**3. [Low, scope creep] `nightly.yml` disables coverage and deletes the "Check Coverage" step inside a version-bump commit.**
The diff sets `coverage: true` to `coverage: false` and removes the whole Check Coverage step. Two points:
- The removed step was malformed YAML in the original (a two-line `run:` value without a `run: |` block scalar), so it likely never ran the second command correctly. Removal is defensible cleanup, not a loss.
- No coverage gap results: `CI.yml:44-63` keeps a `coverage` job that runs `EREBUS_TEST_GROUP: all` with `check_coverage.jl`, and `Coverage.yml` still uploads to Codecov. So the floor still exists elsewhere.

The problem is bundling: a behavior change to the nightly rides inside a commit titled for a Julia version bump, with no note. Confirm the removal is intended and separate it in the commit message.

### What is missing
- No test or CI check guards the compat-floor claim. If the code does not actually need 1.12, nothing catches the over-restriction in Finding 1.
- The change does not touch the ageing action pins (`julia-actions/setup-julia@v1`, `julia-actions/cache@v1`). These are outside the diff scope, so not a defect of this change, but an older `setup-julia@v1` is a candidate reason a new version like 1.13 might not resolve. I did not verify this.

### Bottom line
This is a clean, complete string-level version bump with consistent docs. The two substantive risks are not typos: the LTS drop (Finding 1) is a policy choice that needs a decision, and the existence of Julia 1.13 (Finding 2) is a reachability fact I could not confirm and that will fail the whole matrix if wrong. Verify both before merge.
