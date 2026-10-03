# Intrinsic five-generator coefficient checks

This supplement checks the coefficients in `nonradial_five_generators.tex` by
ordinary Euclidean differentiation of the explicit real and complex models.
It is a regression check of those models, **not** a proof of the universal
structural theorem or its classification-free scope.

Let `Z` be a trace-free real symmetric or complex Hermitian 3-by-3 matrix, with
metric `Re tr(ZZ')/2`, and let `u` belong to `R^3` or `C^3` with its usual real
Euclidean metric. The script constructs orthonormal real coordinates and uses

- `r = tr(Z^2)/2`, `s = u* u`;
- `a = (3 sqrt(3)/2) det Z`, `p = (u* Z u)/(2 sqrt(3))`;
- `h = (u* Z^2 u)/12`, `F = (a+3p)/6`.

All derivatives are computed from these coordinate polynomials, independently
of the differential table being checked. For each model the script checks:

1. the two coordinate metric normalizations;
2. all 15 independent gradient products;
3. the five Laplacians;
4. harmonicity and the normalized Hsiang equation;
5. six generator-recovery identities;
6. one exact nonzero rank-five Jacobian minor.

Run with Python 3 and SymPy installed:

```sh
python check_intrinsic_five_generators.py
```

The script writes `verification_results.json` in its own directory. Every
identity is checked exactly; there is no floating-point sampling. The rank
certificate uses `Z = diag(1,-1,0)` and `u = (1,1,1)`.

The checked identities are coordinate verifications for `d = 1,2`. The abstract
proof for every admissible dimension is given in the manuscript. No native
Mathematica execution is claimed by this supplement.
