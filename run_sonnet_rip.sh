#!/bin/bash
export PATH="$HOME/git/fleet/bin:$PATH"
RIPDIR=$(mktemp -d "${TMPDIR:-/tmp}/rip.XXXXXX")
REPO="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"

git -C "$REPO" diff main...HEAD > "$RIPDIR/wip.patch"

# Prepare Lenses
cat << 'PROMPT1' > "$RIPDIR/lens_1.txt"
You are an adversarial reviewer. Your job is to find problems, not to validate work.
Your default stance is skepticism. Assume there are errors until proven otherwise.
Review this code for correctness bugs, logic errors, and edge cases.
MATERIAL:
PROMPT1
cat "$RIPDIR/wip.patch" >> "$RIPDIR/lens_1.txt"

cat << 'PROMPT2' > "$RIPDIR/lens_2.txt"
You are an adversarial reviewer. Your job is to find problems, not to validate work.
Your default stance is skepticism. Assume there are errors until proven otherwise.
Review this code for what is MISSING, not just what is wrong. Pay special attention to test quality and whether tests actually verify correctness.
MATERIAL:
PROMPT2
cat "$RIPDIR/wip.patch" >> "$RIPDIR/lens_2.txt"

cat << 'PROMPT3' > "$RIPDIR/lens_3.txt"
You are an adversarial reviewer. Your job is to find problems, not to validate work.
Your default stance is skepticism. Assume there are errors until proven otherwise.
Review this code for design problems that will cause pain later.
MATERIAL:
PROMPT3
cat "$RIPDIR/wip.patch" >> "$RIPDIR/lens_3.txt"

echo "Starting Sonnet RIP round..."
seat-run -m sonnet -q "$(cat "$RIPDIR/lens_1.txt")" -C ~/git/Erebus.jl > rip_sonnet_result1.txt 2>rip_sonnet_err1.txt &
P1=$!
seat-run -m sonnet -q "$(cat "$RIPDIR/lens_2.txt")" -C ~/git/Erebus.jl > rip_sonnet_result2.txt 2>rip_sonnet_err2.txt &
P2=$!
seat-run -m sonnet -q "$(cat "$RIPDIR/lens_3.txt")" -C ~/git/Erebus.jl > rip_sonnet_result3.txt 2>rip_sonnet_err3.txt &
P3=$!

wait $P1
wait $P2
wait $P3

echo "=== Sonnet RIP Agent 1 (Correctness) ===" > rip_sonnet_report.md
cat rip_sonnet_result1.txt >> rip_sonnet_report.md
echo -e "\n=== Sonnet RIP Agent 2 (Completeness) ===" >> rip_sonnet_report.md
cat rip_sonnet_result2.txt >> rip_sonnet_report.md
echo -e "\n=== Sonnet RIP Agent 3 (Design) ===" >> rip_sonnet_report.md
cat rip_sonnet_result3.txt >> rip_sonnet_report.md
