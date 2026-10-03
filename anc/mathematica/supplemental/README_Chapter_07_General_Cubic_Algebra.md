# General cubic algebra: supplemental verification

Release: `20260923-v6`. Run `Verify_Chapter_07_General_Cubic_Algebra.nb`, or evaluate
`Get["/absolute/path/Verify_Chapter_07_General_Cubic_Algebra.wl"]` in a fresh kernel.
The source targets Mathematica 12 and Mathics and uses exact polynomial arithmetic.
Its context is ``MergedVerificationV6`Chapter07General` ``.

This unit has **31 assertions** and is independent of the original 13 units and
550 assertions. The retained 13-unit entry point does not invoke this supplemental
unit. A successful execution has 29 exact-zero residuals and two
`True` results. Outputs are written beside the source in
`verification_results_v6_supplemental`.

The final source return value is `True` only if all mathematical assertions pass,
the count is exactly 31, no group reports evaluation messages, and both CSV and
TXT exports succeed. The TXT records mathematical status; the notebook output and source return value also record export status.

## Mathematical scope

Assertions 01–18 independently check coordinate algebra with general coefficients
in three dimensions. They support the displayed identities and do not replace
the dimension-independent proofs. The quadratic weight in 09 and 13–18 is general:
no assumption that it is a multiple of the Euclidean norm is made. Assertion 09
checks the `2/3` coefficient in the multiplication form of the general weight
equation. Assertions 16 and 18 are residual identities conditional on the Hessian
polynomial relation, rather than assertions that an arbitrary symmetric matrix
satisfies that relation.

Assertions 19–27 use the stated differential table with `K = Delta S` represented
by an independent symbol. They neither set `K` to a scalar multiple of `R` nor
require derivatives of `K`; only assertion 27 deliberately tests that specialization.
In particular, freedom of `K` does not make the cubic weight nonradial: assertion
19 checks that the assumed differential identities imply `L(F) = -2 R F`.

Assertions 28–29 substitute the two positive-`S` stationary branch values. They
do not prove existence of a maximum or exhaust the possible branches; those
arguments remain in the manuscript. Assertions 30–31 check nonnegative polynomial
coefficients on the explicitly assumed domain
`R >= 0`, `S >= 0`, `t = F^2 >= 0`, `K - 16 R >= 0`,
`R^2/3 - S >= 0`. At `F != 0`, strict positivity additionally uses the manuscript's
Euler implication `R > 0` and `S > 0` and the summand `150 R^2 S t`.
No sign of `F` is assumed. No `PowerExpand`, negative real power, or radical occurs,
so no unproved positive-base simplification is used.

The program does **not** prove geometric comparison, singular-set measure or
removability, global positivity hypotheses, existence/classification of algebras,
or existence of a realization with nonradial `K`.

## Manuscript mapping

| Assertions | TeX labels / location |
|---|---|
| 01–08 | `eq:cubic-polarized-product`, `eq:cubic-multiplication-dictionary`, `eq:cubic-variable-trace` and adjacent dictionary/annihilator discussion |
| 09 | `eq:cubic-general-weight` and the equivalent algebraic multiplication identity |
| 10–12 | `thm:cubic-one-dimensional-degree-three`, Jordan and Euler operator calculations in its proof |
| 13–16 | `eq:cubic-general-weight`, `eq:cubic-general-hessian-polynomial` |
| 17–18 | `eq:cubic-positive-multiplier-leading-sign` |
| 19 | `eq:cubic-gradient-identities-only` |
| 20–23 | `eq:cubic-variable-trace-gradient` and the following contraction/derivative table |
| 24–26 | `eq:cubic-variable-trace-positive-decomposition` |
| 27 | `rem:cubic-algebraic-certificate-scope` |
| 28–29 | `eq:cubic-gradient-only-S-bound`, stationary equations in its proof |
| 30–31 | `eq:cubic-variable-trace-gradient`, `eq:cubic-variable-trace-positive-decomposition`, under the stated sign domain |

Exact assertion labels and artifact hashes are recorded in `supplemental_manifest.json`.

## Supplied verification record

All **31/31** assertions passed in the author's native Mathematica 12.0.0,
Windows-x86-64 run: **29 exact-zero residuals and two `True` results**, with
no evaluation messages. The groups contain 9, 3, 6 and 13 assertions.
The evaluated notebook and TXT/CSV records are supplied; the TXT records
mathematical status, and the final notebook return additionally requires
both exports to succeed. This documentation update did not rerun the unit.
