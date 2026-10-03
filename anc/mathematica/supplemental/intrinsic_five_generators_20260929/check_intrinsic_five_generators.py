#!/usr/bin/env python3
"""Exact coordinate checks of the real and complex five-generator models.

This script does not assume the displayed differential table in computing
derivatives.  It constructs Hermitian matrices in orthonormal real coordinates,
then compares their ordinary Euclidean derivatives with the claimed formulas.
It verifies these explicit models, not the universal structural theorem.
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path

import sympy as sp


def verify_model(d: int) -> dict:
    started = time.monotonic()
    sqrt3 = sp.sqrt(3)
    x, z = sp.symbols("x z", real=True)
    real_off = sp.symbols("b12 b13 b23", real=True)
    imaginary_off = sp.symbols("c12 c13 c23", real=True) if d == 2 else ()
    off = [real_off[i] + (sp.I * imaginary_off[i] if d == 2 else 0)
           for i in range(3)]
    Z = sp.Matrix([
        [x + z / sqrt3, off[0], off[1]],
        [sp.conjugate(off[0]), -x + z / sqrt3, off[2]],
        [sp.conjugate(off[1]), sp.conjugate(off[2]), -2 * z / sqrt3],
    ])
    ure = sp.symbols("u1 u2 u3", real=True)
    uim = sp.symbols("v1 v2 v3", real=True) if d == 2 else ()
    u = sp.Matrix([ure[i] + (sp.I * uim[i] if d == 2 else 0)
                   for i in range(3)])
    b_coords = (x, z) + real_off + imaginary_off
    u_coords = ure + uim
    coords = b_coords + u_coords
    r = sp.expand(sp.trace(Z * Z) / 2)
    s = sp.expand((u.conjugate().T * u)[0])
    a = sp.expand(3 * sqrt3 * Z.det() / 2)
    p = sp.expand((u.conjugate().T * Z * u)[0] / (2 * sqrt3))
    h = sp.expand((u.conjugate().T * Z * Z * u)[0] / 12)
    functions = {"r": r, "s": s, "a": a, "p": p, "h": h}
    checks = []

    def check(name, lhs, rhs):
        residual = sp.expand(lhs - rhs)
        if residual != 0:
            residual = sp.simplify(residual)
        passed = residual == 0
        checks.append({"name": name, "passed": bool(passed)})
        if not passed:
            raise AssertionError(f"d={d}, {name}: {residual}")

    def grad(f):
        return tuple(sp.diff(f, c) for c in coords)

    def gamma(f, g):
        return sp.expand(sum(v * w for v, w in zip(grad(f), grad(g))))

    def laplacian(f):
        return sp.expand(sum(sp.diff(f, c, 2) for c in coords))

    # This first check ensures that the displayed coordinates have exactly the
    # metric used by the manuscript, so no inverse metric factors are omitted.
    check("orthonormal_B_coordinates", r, sum(c**2 for c in b_coords))
    check("orthonormal_U_coordinates", s, sum(c**2 for c in u_coords))

    # Independent Euclidean derivatives, compared with all 15 table entries.
    expected_gamma = {
        ("r", "r"): 4*r,
        ("r", "s"): 0,
        ("r", "a"): 6*a,
        ("r", "p"): 2*p,
        ("r", "h"): 4*h,
        ("s", "s"): 4*s,
        ("s", "a"): 0,
        ("s", "p"): 4*p,
        ("s", "h"): 4*h,
        ("a", "a"): 9*r**2,
        ("a", "p"): 18*h-r*s,
        ("a", "h"): r*p+s*a/3,
        ("p", "p"): s**2/9+4*h,
        ("p", "h"): r*p/3+2*s*p/9+s*a/27,
        ("h", "h"): (r+s)*h/3+p**2/9+a*p/27,
    }
    for (first, second), value in expected_gamma.items():
        check(f"Gamma({first},{second})", gamma(functions[first], functions[second]), value)

    m = len(u_coords)
    expected_delta = {"r": 2*(m+2), "s": 2*m, "a": 0, "p": 0,
                      "h": (m*r+(m+2)*s)/9}
    for name, value in expected_delta.items():
        check(f"Delta({name})", laplacian(functions[name]), value)

    R = r+s
    F = sp.expand((a+3*p)/6)
    S = gamma(F, F)
    q = sp.Rational(3, 2)*(r+sp.Rational(7, 27)*s)
    L = sp.expand(S*laplacian(F)-gamma(F, S)/2)
    check("F_is_harmonic", laplacian(F), 0)
    check("normalized_Hsiang_equation", L, -q*F)
    D = gamma(r, F)
    check("D_equals_a_plus_p", D, a+p)
    check("recover_r", (18*q-7*R)/20, r)
    check("recover_s", (27*R-18*q)/20, s)
    check("recover_p", 3*F-D/2, p)
    check("recover_a", 3*D/2-3*F, a)
    check("recover_h", (gamma(a, p)+r*s)/18, h)

    # A concrete rank-five point: Z = diag(1,-1,0), u = (1,1,1).
    # These same five real coordinate directions suffice also in d=2.
    point = {c: 0 for c in coords}
    point[x] = 1
    point.update({c: 1 for c in ure})
    columns = (x, z) + ure
    jacobian = sp.Matrix([[sp.diff(f, c) for c in columns]
                          for f in (r, a, s, p, h)])
    minor = sp.simplify(jacobian.subs(point).det())
    check("rank_five_exact_minor", minor, 4*sqrt3/3)

    return {
        "division_algebra": "R" if d == 1 else "C",
        "d": d,
        "ambient_real_dimension": len(coords),
        "B_dimension": len(b_coords),
        "U_dimension": m,
        "checks_passed": len(checks),
        "checks": checks,
        "rank_certificate": {
            "Z": "diag(1,-1,0)", "u": "(1,1,1)",
            "invariant_row_order": ["r", "a", "s", "p", "h"],
            "coordinate_column_order": [str(c) for c in columns],
            "jacobian_minor": str(minor),
        },
        "elapsed_seconds": round(time.monotonic()-started, 3),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("verification_results.json"))
    args = parser.parse_args()
    models = [verify_model(d) for d in (1, 2)]
    report = {
        "scope": "Exact ordinary Euclidean differentiation of the explicit real and complex Hermitian models; not a proof of the universal structural theorem.",
        "software": {"python": platform.python_version(), "sympy": sp.__version__},
        "normalization": {
            "r": "tr(Z^2)/2", "s": "u* u", "a": "(3 sqrt(3)/2) det Z",
            "p": "(u* Z u)/(2 sqrt(3))", "h": "(u* Z^2 u)/12",
            "F": "(a+3p)/6", "q": "(3/2)(r+7s/27)",
        },
        "models": models,
        "total_checks_passed": sum(m["checks_passed"] for m in models),
        "all_passed": True,
        "native_mathematica_executed": False,
    }
    args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"all_passed": True, "total_checks_passed": report["total_checks_passed"],
                      "models": [{k: m[k] for k in ("d", "ambient_real_dimension", "checks_passed", "elapsed_seconds")}
                                 for m in models], "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
