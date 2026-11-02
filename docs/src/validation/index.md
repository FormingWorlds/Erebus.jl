# Validation

## What this section checks

Validation means that we run the code on a case where we know the right answer from another source, and we check that the code gives that answer. The other source can be an exact mathematical solution, a second independent solver, or a published measurement.

The pages in [Physics Function Reference](../reference/benchmarks/index.md) mostly check one function at a time. The pages here check how the parts work together, in most cases with the complete `simulation_loop`, where a fault in any step of the time loop changes the result.

Each case on a page has the same layout:

- why the check exists, which names the fault it catches;
- what it tests, with the setup and the expected result;
- how we know the right answer, which names the source of the reference;
- why it matters, which names the model results that depend on it;
- a figure, a note on how to read it, and a box with the test details.

## Pages

| Page | Part of the model | Reference |
|---|---|---|
| [Spherical Heat Conduction](thermal_conduction.md) | Heat conduction, the spherical geometry term, radiogenic heating in the time loop | Exact solution, backward-Euler series, 1-D finite-volume solver |

## How the results are stored

Each validation page has a record file in `benchmarks/records/`. The record holds the measured values, the tolerances, a hash of the configuration, and the last commit that changed `src/`. The table on the page comes from the record, and `tools/check_doc_numbers.jl` fails if the two differ.

`Pkg.test` runs each case again and compares it with its record. To write a new record after an intended physics change, run the export tool for that record. The tool also writes the table on the page. It refuses to run with uncommitted changes in `src/` or `benchmarks/`. It also refuses to write a value that leaves its regression band (the allowed range around the recorded value), unless you pass `--accept-physics-change`, and a value that fails its reference tolerance in all cases.
