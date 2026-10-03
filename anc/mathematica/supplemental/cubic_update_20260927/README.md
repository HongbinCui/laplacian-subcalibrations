# Exact checks for cubic certificates

This supplement contains **38 exact assertions** for the cubic chapters
and nonharmonic appendix, separate from the other verification units.

The two implementations recompute derivatives from the displayed differential
tables and compare exact rational expressions. Each check clears denominators
and asserts that the resulting polynomial is identically zero. A failed check
stops execution. No floating-point arithmetic or sampling is used.

## Run

Python 3 with SymPy:

```sh
python3 verify_cubic_update.py > python_execution.log 2>&1
```

Install SymPy if it is not already available. The supplied execution used
Python 3.12.14 and SymPy 1.14.0. A successful run prints
`SUCCESS: 38 exact assertions passed.` and writes
`python_execution_metadata.json`, including the assertion names, engine
versions, UTC execution time, and SHA-256 hashes of the two source programs.

Wolfram Mathematica 12, in a fresh kernel:

```sh
wolframscript -file Verify_Cubic_Update_20260927.wl > mathematica_execution.log 2>&1
```

Alternatively, evaluate `Get["/full/path/Verify_Cubic_Update_20260927.wl"]`
in a fresh Mathematica kernel and save the output. The script clears the
`Global` context at entry and exits the kernel with status 1 upon failure.
The final line records the engine version; retain that information with the
author's execution record.

**Native Mathematica status — 2 October 2026: all 38 checks passed.**
The author ran Mathematica 12.0.0, Windows-x86-64, at 17:38:31 (kernel-recorded
time). The evaluated `Verify_Cubic_Update_20260927.nb` and
`notebook_run_38_output.txt`, `notebook_run_38_results.csv`, and
`notebook_run_38_summary.json` preserve the run. Mathematical checks,
evaluation and record exports succeeded, no evaluation messages were
recorded, and the source SHA-256 matches the supplied `.wl`.

To reproduce, open the notebook in a fresh kernel and choose
Evaluation > Evaluate Notebook; require `COMPLETE AND PASSED: True`.
The shared helper is `../Notebook_Run_Helper.wl` version r2.

The supplied Python log and metadata preserve an earlier execution of the
same assertions before local variable renaming. Their source hashes refer
to that execution; the native notebook-run summary identifies the current
Wolfram source. The additional implementation is not counted a second time.

## Correspondence with the revised paper

| Checks | Mathematical content | Source labels |
|---|---|---|
| 1–4 | Gradient, Laplacian, and H/J defect identity for the unified sextic potential; trace-corollary remainder | General cubic section / nonharmonic appendix: `eq:cubic-unified-gradient`, `eq:cubic-unified-laplacian`, `eq:cubic-unified-defect`, `eq:cubic-unified-J`, `cor:cubic-unified-trace` |
| 5–12 | Five normalized Clifford gradient pairings, normalized Laplacian, H44 decomposition, and parameter increment | Clifford chapter: `eq:Clifford-normalized-data`, `eq:Clifford-unified-table`, `eq:Clifford-H44`, `prop:Clifford-common-potential` |
| 13–14 | Endpoint and interior values used by the shortened all-powers obstruction | Clifford chapter: `eq:Clifford-boundary-coefficient`, `eq:obstruction-lambda0`, `eq:obstruction-lambda23` |
| 15–19 | General radial multiplication formula and its degree-five, dimension-26 specialization | Appendix A: `eq:cubic-radial-lift`; applied in `thm:nonradial-minimizing`; also used in the nonharmonic appendix |
| 20–27 | Five-invariant nonradial closure consequences, derivatives of H, and exact N defect | Nonradial chapter: `eq:nonradial-gradient-square`, `eq:nonradial-minimal-equation`, `eq:nonradial-T`, `eq:nonradial-H-derivatives`, `eq:nonradial-certificate`, `eq:nonradial-N` |
| 28–32 | Seven-summand polynomial identity and arithmetic supporting its bounds and the H Delta(H) lower bound | Nonradial chapter: `eq:nonradial-positive-certificate`, `eq:nonradial-bounds`, `prop:cubic-radial-lift` |
| 33–38 | General nonharmonic H derivatives, B/C identities with the full Bochner correction, and the harmonic radial specialization | General cubic section / nonharmonic appendix: `eq:cubic-nonharmonic-defect`, `eq:cubic-nonharmonic-leading-H`; Nonharmonic appendix: `eq:cubic-nonharmonic-B0`–`eq:cubic-nonharmonic-B2`, `eq:cubic-nonharmonic-coefficients`, `eq:cubic-nonharmonic-H-defect` |

Stable source labels are given because theorem and equation numbers can change
during editing. Older directories named `Chapter_07` belong to the earlier
organization and now correspond to material in Chapters 7 and 8. An older
`Appendix_F` verification entry concerns historical homogeneity checks, whereas
the current Appendix F contains the nonharmonic certificate. The former Jordan-recovery appendix was replaced by elementary arguments in the intrinsic chapter (current Chapter 10). Those
older program names have not been repurposed in this supplement.

## Scope and normalizations

- In the Clifford block, `ff` is the original manuscript cubic; the criterion
  is checked for `f = ff/2` and `s = |grad f|^2 = S/4`. The obstruction
  variable `jl` denotes current `J = |grad S|^2` on its normalized edge.
  The normalized Laplacian remains `k = Delta s` in this program.
- In the nonradial block, `w = 4 u* Z^2 u` is a quartic generator, while `p`
  is the squared gradient norm. The closure is specialized to quaternionic
  dimension `d = 4`, hence ambient dimension 26.
- In the nonharmonic block, the contractions `gft`, `gfq`, etc. are formal
  symbols. The substitution `gft = k/2 - h2` retains the full Bochner identity;
  here `h2 = H_2 = |Hess F|^2`, while `tt = T = Gamma(S,S)`;
  it does not impose harmonicity. The polynomial positive-semidefiniteness
  assumptions themselves are not asserted to hold for an unknown model.

The computations verify the displayed **identities conditional on the input
differential tables**. They do not prove the tables from matrix coordinates,
spectral classification or Jordan recovery, the analytic inequalities on the
actual invariant image, singular-set removability, or the geometric comparison
theorem. Those points are proved in the manuscript. In particular, an exact
seven-summand decomposition does not by itself establish nonnegativity without
the four inequalities proved in the text. No nonharmonic irreducible minimal
cubic is constructed or certified here.
