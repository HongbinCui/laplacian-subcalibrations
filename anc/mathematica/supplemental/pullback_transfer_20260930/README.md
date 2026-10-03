# Pullback transfer supplement

This supplement checks differential formulas used in the
Riedler pullback chapter (Chapter 11). Its 13 checks form a separate unit.

## Execution status

- **Executed:** `verify_exact.py`, Python 3.12.14, exact arithmetic using the
  standard-library `fractions.Fraction`; all **13 identities passed**.
- **Native Mathematica executed on 2 October 2026:** all **13 checks passed**
  in Mathematica 12.0.0, Windows-x86-64, at 17:39:14 (kernel-recorded time).
  The evaluated `verify_pullback_transfer.nb`, `notebook_run_13_output.txt`,
  `notebook_run_13_results.csv`, `notebook_run_13_summary.json`, and
  `mathematica_native_results.json` preserve this run. Checks, evaluation
  and exports succeeded with no evaluation messages. The summary's source
  SHA-256 matches the supplied `.wl` file.
- `exact_results.json` remains the earlier Python exact-arithmetic report;
  it is separate from the new native Mathematica records.
- The two scripts perform the same identity checks with separate
  implementations. The Python implementation uses sparse Laurent
  polynomials and exact coefficient comparison. The Wolfram Language
  implementation uses differentiation and `FullSimplify`.

To reproduce in Mathematica, open `verify_pullback_transfer.nb` in a fresh
kernel and choose Evaluation > Evaluate Notebook. Require
`COMPLETE AND PASSED: True`; retain `../Notebook_Run_Helper.wl` r2 and save
the evaluated notebook. The helper exports TXT, CSV and JSON records.

Alternatively, run the original native script:

```sh
wolframscript -file verify_pullback_transfer.wl
```

A successful run prints thirteen `PASS` lines and writes
`mathematica_native_results.json`, including the actual kernel version.
To reproduce the executed exact check:

```sh
python3 verify_exact.py
```

This regenerates `exact_results.json`. No third-party Python packages are
required.

## Notation and checks

Write `R=|x|^2` for the source squared radius, `q=|y|^2` for the target
squared radius, `h=H(P)`, `W_H=|grad H|^2`, `D_H=Delta H`, and
`L_H=L(H)`. Target quantities in source formulas are evaluated at `P`.
The source gradient square is `|grad_x(H o P)|^2 = 4 R W_H(P)`;
`W_H(P)` itself is the target gradient square evaluated at the image.

The operator convention is

```text
L(v) = |grad v|^2 Delta v - (1/2) Gamma(v,|grad v|^2).
```

The calculations first derive the radial multiplication identity for an
arbitrary homogeneous source function, using Euler's identity and the
definition of `L`. They separately check the Clifford chain-rule identity

```text
L(h) = 16 R^2 L_H(P) - 8 ell h W_H(P).
```

They then check the gradient, defect, and preserved ratio identities for
`G=R^a h`, and the consistency of the weighted Laplacian under Clifford
specialization, together with the finite difference of
`G L(G)` when the source dimension changes at fixed `(R,P)`. Powers `R^a`, `R^(2a)`, and `R^(3a)` are factored
out before calculation. Consequently the general real exponent `a`
remains a free indeterminate without requiring branch manipulations.

The weighted-Laplacian assertion specializes the entered radial formula;
it is not a new independent differentiation. An independent second-order
chain-rule check is provided by `CheckRadialMaster` in
`../../chapters/Verify_Appendix_A.wl` (the `radial Laplacian` assertion).

The target Simons example is calculated independently from its quotient
relations

```text
Gamma(q,q)=4q, Gamma(q,f)=4f, Gamma(f,f)=4q,
Delta q=4k, Delta f=0.
```

For `H=q^(3/4) f` the result is

```text
W_H = 4 q^(5/2) + (33/4) q^(1/2) f^2,
Delta H = (3/4)(4k+7) q^(-1/4) f,
L(H) = q^(1/4) f [(12k-47)q^2 + (33/8)(6k+1)f^2],
H Delta H = (3/4)(4k+7)q^(1/2)f^2.
```

Therefore `H L(H)>=0` for `k>=4`, and `H Delta H>=0`, with
`ell=7/2` and ratio constant `M=0`. Taking `a=1` in the general dimension
bound gives `N>=41/2` (equivalently integer `N>=21`). This is compatible
with the admissible Clifford family dimensions `N>=32` used in the
chapter; it does not assert that a Clifford representation exists at
every integer dimension satisfying the analytic bound.

These are checks of differential identities. They do not independently
establish the singular-set cutoff argument, area minimization as an
integral current, the smallest-algebra identity, or literature priority.
