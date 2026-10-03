# Verification guide

All 660 Wolfram assertions have passing native Mathematica 12.0.0 records.
The two Python-only units add 127 passing assertions. See
`../../Verification_Records.md` for the evidence inventory and
`../../Verification_Details.pdf` for the detailed supplement.

Program names, contexts and result directories are stable identifiers;
they do not track manuscript chapter numbers or publication versions.
Keep each `.nb` wrapper, its `.wl` source and the relative directory layout.
Fully extract the archive before running any notebook.

## Main suite: 550 assertions

Start a fresh kernel, open `../../START_HERE_Verification_v5.nb`, and choose
**Evaluation > Evaluate Notebook**. The entry reports 13 units and 550
assertions. Every unit must pass and export both TXT and CSV records.
The results are written to `verification_results_v6/`. Save the evaluated
notebook explicitly. Its displayed V6 release identifier is intentional.

To run a single unit, open its same-named `.nb`, or evaluate in a saved
notebook in this directory:

```wl
Block[{$IterationLimit = Infinity,
       $ContextPath = {"System`"}, $Context = "Global`"},
 Get[FileNameJoin[{NotebookDirectory[], "Verify_Chapter_08.wl"}]]]
```

| Stable unit | Current manuscript location | Scope | Assertions |
|---|---|---|---:|
| Chapter04 | Chapter 4 | General multipliers, spherical formulas, Hardy square and Simons potential | 6 |
| Chapter05 | Chapter 5 | Isoparametric corrections, Lawson and Lin potentials | 16 |
| Chapter06 | Chapter 6 | Bihomogeneous quartics, Hopf constructions and Hardy coefficients | 49 |
| Chapter07 | Chapters 7–8 | Specialized SF identities, normalized closure, spectrum and quaternionic model | 18 |
| Chapter08 | Chapter 9 | Clifford cubic calculations | 193 |
| Chapter09 | Chapter 11 | Riedler, Cartan and Lawson pullbacks | 50 |
| Chapter10 | Chapter 12 | Finite-order determinant recurrences and scalar inequalities | 134 |
| AppendixA | Appendix A | Gradient, radial and combined multipliers; boundary coefficient | 8 |
| AppendixB | Appendix B | Bernstein coefficients and Sturm endpoint signs | 15 |
| AppendixC | Appendix C | Repeated cutoff-exponent check | 1 |
| AppendixD | Appendix D | Odd-degree rigidity and moment elimination | 15 |
| AppendixE | Ancillary Hessian note | Formal Clifford-algebra Hessian contractions | 18 |
| AppendixF | Retained homogeneity checks | Homogeneity and endpoint denominator cancellation | 27 |
| **Total** | | **13 units** | **550** |

Stable AppendixF does not refer to current Appendix F, which contains the
nonharmonic certificate checked in the separate 38-assertion unit.
The Hessian derivation is in `../../notes/Clifford_Hessian_Identities.pdf`.
Assertion labels, scopes and TeX correspondence are in
`verification_manifest.json`. The total includes the stated repeated
cutoff-exponent check.

After running Chapter08, selected outputs are available as:

```wl
MatrixForm /@ MergedVerificationV5`Chapter08`BernsteinMatrices
MergedVerificationV5`Chapter08`JacobiMatrix // MatrixForm
MergedVerificationV5`Chapter08`NegativeJacobiValue
```

The last value is `-9952/969969`.

## Supplemental units

Paths below are relative to `../supplemental/`. Use a fresh kernel for each
Wolfram notebook and save the evaluated notebook. The 28-, 38- and 13-check
notebooks require `Notebook_Run_Helper.wl`; their final gate must be
`COMPLETE AND PASSED: True`.

| Entry | Assertions | Current scope |
|---|---:|---|
| `Verify_Chapter_07_General_Cubic_Algebra.nb` | 31 | Chapter 7 cubic dictionary, Jordan/Euler identities, general-q jets and variable-K SF identities |
| `exceptional_cubics/Verify_Exceptional_Cubics.nb` | 28 | Chapter 8 exceptional cubic calculations |
| `cubic_update_20260927/Verify_Cubic_Update_20260927.nb` | 38 | Chapters 7–9 and Appendix F identities and certificate decompositions |
| `pullback_transfer_20260930/verify_pullback_transfer.nb` | 13 | Chapter 11 radial and Clifford transfer identities |
| `jordan_pencil_20260928/verify_jordan_pencil.py` | 65 | Chapter 8 pencil and real-coordinate checks |
| `intrinsic_five_generators_20260929/check_intrinsic_five_generators.py` | 62 | Chapter 10 intrinsic table in real and complex coordinate models |

The two Python-only units require Python 3 and SymPy; local README files
state commands and dependencies. The pullback `verify_exact.py` supplies
another implementation of its same 13 checks, using only Python's standard
library. Additional implementations are not counted again.

## Mathematical scope and notation

The checks use exact integers, rational or algebraic numbers, and symbolic
expressions. `PowerExpand` is restricted to the positive bases specified in
the manuscript; it does not assume F positive or justify extension across
S=0 or Q=0. The general-cubic 31-check unit uses polynomial operations.

The determinant unit samples orders 2–10 and scalar parameters 2–100;
these finite checks do not replace the all-orders proof. Chapter08 includes
110 coordinate-contraction checks at specified points in two Clifford
systems; AppendixE separately verifies formal Clifford-algebra contractions.
AppendixC checks an exponent, not the local geometric construction.
The general structural, spectral, positivity and geometric arguments remain
in the manuscript.

In the pullback supplement, the weighted-Laplacian assertion checks the
Clifford specialization of the entered radial formula. Independent
second-order differentiation is supplied by `CheckRadialMaster` in
`Verify_Appendix_A.wl`. Other local variable conventions and normalizations
are recorded in `../../NOTATION_GUIDE.md`.

The total is 550+31+28+38+13+65+62=787 program assertions, including stated
overlaps. Preparing the compact package did not change or rerun programs
or relabel execution records.
