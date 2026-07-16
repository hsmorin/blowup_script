# Q1 double-octic audit

The main deliverable is `q1_revised_report.pdf`; its editable source is
`q1_revised_report.tex`. Statements shown in red are mistakes in either the
expository PDF or the earlier Rust audit report.

## Reproducible checks

- `q1_exact_checks.sing` checks the seven transverse Milnor numbers, the
  isolated point, its Hessian, and the three proposed small-resolution
  divisors. Run it with `Singular -q q1_exact_checks.sing`.
- `q1_stage33_curve_checks.sh` verifies that the displayed ideals in terminal
  Rust charts 212 and 213 define one-dimensional singular subschemes. Run it
  from this directory with `./q1_stage33_curve_checks.sh`.
- The corresponding captured outputs are `q1_exact_checks.txt` and
  `q1_stage33_curve_checks.txt`.
- `q1_revised_findings.json` records the corrected conclusions in a compact,
  machine-readable form.

## Source and prior artifacts

- `q1_expository_report.pdf` is the newly supplied expository source.
- `q1_branch_polynomial.txt` is the expanded branch octic used by the Rust
  program and agrees exactly with the discriminant in the source PDF.
- `seven_line_*.json` are the prior Rust artifacts. They are retained as audit
  evidence, but they do not document a completed resolution. In particular,
  `seven_line_initial.json` omits the isolated singular point and
  `seven_line_stage33_raw.json` ends with its terminal audit disabled.
