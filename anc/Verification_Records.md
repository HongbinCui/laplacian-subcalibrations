# Verification records

The supplied native records were checked against the author's uploads of
2 October 2026. Preparing this compact package did not execute the programs.

| Unit | Assertions | Supplied evidence |
|---|---:|---|
| Main 13-unit suite | 550 | Mathematica 12.0.0, Windows-x86-64; 2 October 2026, 15:48:47–15:48:48; all passed |
| General cubic algebra | 31 | Evaluated notebook and native Mathematica 12.0.0 TXT/CSV; all passed; supplied 2 October |
| Exceptional cubics | 28 | Mathematica 12.0.0; 2 October, 17:37:56–17:37:57; all passed |
| Cubic update | 38 | Mathematica 12.0.0; 2 October, 17:38:31; all passed |
| Pullback transfer | 13 | Mathematica 12.0.0; 2 October, 17:39:14; all passed; separate exact-Python results also included |
| Jordan pencil | 65 | Passing Python/SymPy results |
| Intrinsic five-generator table | 62 | Passing Python/SymPy coordinate results |

Times are copied from the native logs, which do not state a time-zone offset.
The 17 canonical native CSV tables contain **660 rows, all True**: 13 main
unit tables, one general-cubic table and three notebook-run tables. Together
with the two Python-only units, the inventory contains **787 program
assertions including overlaps**, not 787 distinct mathematical results.
Repeated implementations of the same assertions are not counted again.

## Evidence locations

Paths below are relative to this directory.

- Main outputs: `mathematica/chapters/verification_results_v6/`;
  evaluated entry notebook: `START_HERE_Verification_v5.nb`.
- General-cubic outputs:
  `mathematica/supplemental/verification_results_v6_supplemental/`;
  evaluated notebook in its parent directory.
- Exceptional, cubic-update and pullback native outputs: respectively
  `mathematica/supplemental/exceptional_cubics/`,
  `mathematica/supplemental/cubic_update_20260927/` and
  `mathematica/supplemental/pullback_transfer_20260930/`.
  Each contains the evaluated `.nb`, `notebook_run_N_output.txt`,
  `notebook_run_N_results.csv` and `notebook_run_N_summary.json`.
- Python-only outputs: `verification_results.json` in
  `mathematica/supplemental/jordan_pencil_20260928/` and
  `mathematica/supplemental/intrinsic_five_generators_20260929/`.

The three notebook-run summaries report their expected counts (28, 38 and
13), successful mathematical checks and evaluation, and no evaluation
messages. Their source SHA-256 hashes match the supplied `.wl` files.
The evaluated notebooks record successful exports and
`COMPLETE AND PASSED: True`. They use `Notebook_Run_Helper.wl` r2.

The cubic-update Python log and metadata are an earlier execution of the
same assertions before local variable renaming; their recorded source
hashes refer to that execution, not to a new run of the current files.
The pullback `exact_results.json` is the separate Python rational-arithmetic
record. Neither is counted as an additional set of assertions.

## Stable identifiers and scope

The retained executable files, evaluated notebooks and result records are
unchanged. Superseded compatibility logs and notation-audit archives are
omitted. No execution record or source hash has been relabeled as a new run.
Dates and v5/v6 in program names and output paths remain stable identifiers;
current manuscript mapping is in `mathematica/chapters/VERIFICATION_GUIDE.md`.

The programs check the specified identities, finite tests and coordinate
models. They do not replace the manuscript's general structural theorems,
actual-image inequalities, singular-set analysis or geometric comparison.
