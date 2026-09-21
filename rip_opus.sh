#!/bin/bash
DIFF="$(cat pr_diff.txt)"
TEST_FRAMEWORK="Test framework rules enforced by tools/check_test_quality.jl: no float equality with ==, leaf testsets >= 2 assertions, no bare positivity tests, no eltype/typeof tests."

export PATH="$HOME/git/fleet/bin:$PATH"

# Agent 1
cat << PROMPT1 > rip_opus_agent1.txt
You are an adversarial reviewer. Your job is to find problems, not to validate work.
Your default stance is skepticism. Assume there are errors until proven otherwise.

Rules:
- Do NOT soften findings. "This might be a minor issue" is banned. State what is wrong.
- Do NOT excuse problems because "it works for now" or "it's probably fine."
- Do NOT pad your report with praise. Zero compliments. Only findings.
- If you find nothing wrong after genuine effort, say "no issues found" — do not manufacture problems, but also do not stop looking early.
- Confidence scoring: rate each finding 0-100. Only report findings >= 70.
- For each finding: one-line summary, then evidence, then concrete fix.

Review this code for correctness bugs, logic errors, and edge cases.

MATERIAL:
$DIFF

Focus areas:
- Off-by-one errors, boundary conditions, integer overflow
- Null/None handling, empty collections, missing keys
- Race conditions, state mutation, ordering assumptions
- Error handling: swallowed exceptions, bare except, missing cleanup
- Return value correctness: wrong type, wrong value, missing return
- API contract violations: calling functions with wrong arguments, wrong types, wrong order
- Resource leaks: unclosed files, connections, locks
- Security: injection, unsanitized input, hardcoded secrets, path traversal

For each finding: what is wrong, what input triggers it, what is the fix.
PROMPT1

# Agent 2
cat << PROMPT2 > rip_opus_agent2.txt
You are an adversarial reviewer. Your job is to find problems, not to validate work.
Your default stance is skepticism. Assume there are errors until proven otherwise.

Rules:
- Do NOT soften findings. "This might be a minor issue" is banned. State what is wrong.
- Do NOT excuse problems because "it works for now" or "it's probably fine."
- Do NOT pad your report with praise. Zero compliments. Only findings.
- If you find nothing wrong after genuine effort, say "no issues found" — do not manufacture problems, but also do not stop looking early.
- Confidence scoring: rate each finding 0-100. Only report findings >= 70.
- For each finding: one-line summary, then evidence, then concrete fix.

Review this code for what is MISSING, not just what is wrong. Pay special attention to test quality and whether tests actually verify correctness.

MATERIAL:
$DIFF

REPO TEST FRAMEWORK STANDARDS:
$TEST_FRAMEWORK

Focus areas — code gaps:
- Missing error handling for failure modes that WILL occur in production
- Missing tests for new or changed code paths
- Missing edge cases

Focus areas — test quality:
- Do tests actually test what they claim?
- Float comparison via == instead of ≈
- Missing edge case and error path tests for new code

For each finding: what is missing, why it matters, where to add it.
PROMPT2

# Agent 3
cat << PROMPT3 > rip_opus_agent3.txt
You are an adversarial reviewer. Your job is to find problems, not to validate work.
Your default stance is skepticism. Assume there are errors until proven otherwise.

Rules:
- Do NOT soften findings. "This might be a minor issue" is banned. State what is wrong.
- Do NOT excuse problems because "it works for now" or "it's probably fine."
- Do NOT pad your report with praise. Zero compliments. Only findings.
- If you find nothing wrong after genuine effort, say "no issues found" — do not manufacture problems, but also do not stop looking early.
- Confidence scoring: rate each finding 0-100. Only report findings >= 70.
- For each finding: one-line summary, then evidence, then concrete fix.

Review this code for design problems that will cause pain later.

MATERIAL:
$DIFF

Focus areas:
- Coupling: does this change tie together things that should be independent?
- Abstraction: is the abstraction level appropriate?
- Naming: do names accurately describe behavior?
- Fragility: will this break when reasonable neighboring changes are made?

For each finding: what is the design problem, what is the concrete risk, what is the fix.
PROMPT3

echo "Starting Opus RIP round..."
seat-run -m opus -q "$(cat rip_opus_agent1.txt)" -C ~/git/Erebus.jl > rip_opus_result1.txt 2>rip_opus_err1.txt &
P1=$!
seat-run -m opus -q "$(cat rip_opus_agent2.txt)" -C ~/git/Erebus.jl > rip_opus_result2.txt 2>rip_opus_err2.txt &
P2=$!
seat-run -m opus -q "$(cat rip_opus_agent3.txt)" -C ~/git/Erebus.jl > rip_opus_result3.txt 2>rip_opus_err3.txt &
P3=$!

wait $P1
echo "Agent 1 done"
wait $P2
echo "Agent 2 done"
wait $P3
echo "Agent 3 done"

echo "=== Opus RIP Agent 1 (Correctness) ===" > rip_opus_report.md
cat rip_opus_result1.txt >> rip_opus_report.md
echo -e "\n=== Opus RIP Agent 2 (Completeness) ===" >> rip_opus_report.md
cat rip_opus_result2.txt >> rip_opus_report.md
echo -e "\n=== Opus RIP Agent 3 (Design) ===" >> rip_opus_report.md
cat rip_opus_result3.txt >> rip_opus_report.md
