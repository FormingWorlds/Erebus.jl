#!/bin/bash
while true; do
  OUTPUT=$(gh pr checks 76)
  if echo "$OUTPUT" | grep -q "Some checks are still pending"; then
    sleep 30
  elif echo "$OUTPUT" | grep -q "Some checks were not successful"; then
    echo "FAILED"
    exit 1
  elif echo "$OUTPUT" | grep -q "All checks were successful"; then
    gh pr merge 76 --squash --admin --delete-branch
    echo "MERGED"
    exit 0
  else
    echo "UNKNOWN STATE"
    echo "$OUTPUT"
    sleep 30
  fi
done
