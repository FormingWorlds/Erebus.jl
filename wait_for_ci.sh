#!/bin/bash
while true; do
  STATUS=$(gh pr checks 76 --json state 2>/dev/null)
  # Actually, `gh pr checks` doesn't have a reliable JSON output for aggregate status that's simple. 
  # Better to use `gh pr status` or just `gh pr checks` and parse it.
  
  if gh pr checks 76 | grep -q "pending"; then
    sleep 30
  elif gh pr checks 76 | grep -q "failing\|cancelled"; then
    echo "CI failed or was cancelled."
    exit 1
  else
    echo "CI passed!"
    exit 0
  fi
done
