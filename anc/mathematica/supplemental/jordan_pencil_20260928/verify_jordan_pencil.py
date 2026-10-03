#!/usr/bin/env python3
"""Exact Jordan-pencil checks, with an independent real-coordinate check.

Requires Python 3 and SymPy. Run `python verify_jordan_pencil.py`; the JSON
report is written next to this script. All identities use rational
coefficients and formal variables, never floating-point sampling.
"""
from datetime import datetime, timezone
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import sympy as sp

Q = sp.Rational
CHECKS = []


def zero(name, expression, category):
    residual = sp.expand(sp.fraction(sp.together(expression))[0])
    passed = residual == 0
    row = {"name": name, "category": category, "passed": passed}
    if not passed:
        row["residual"] = str(residual)
    CHECKS.append(row)
    return passed


def condition(name, value, category):
    passed = bool(value)
    CHECKS.append({"name": name, "category": category, "passed": passed})
    return passed


def grad(expression, variables):
    return sp.Matrix([sp.diff(expression, x) for x in variables])


def gamma(a, b, variables, metric_inverse):
    return (grad(a, variables).T * metric_inverse * grad(b, variables))[0]


def lap(expression, variables, table, generator_laplacians):
    return sum(sp.diff(expression, x) * lx
               for x, lx in zip(variables, generator_laplacians)) + sum(
        table[i, j] * sp.diff(expression, x, y)
        for i, x in enumerate(variables)
        for j, y in enumerate(variables))


def formal_closure():
    rho, nu, A, B, W = variables = sp.symbols("rho nu A B W")
    t, d = sp.symbols("t d", real=True)
    # Transcribed from nonradial_models.tex, equations
    # nonradial-laplacians, nonradial-quadratic-pairings, nonradial-closure.
    table = sp.Matrix([
        [4*rho, 0, 6*A, 2*B, 4*W],
        [0, 4*nu, 0, 4*B, 4*W],
        [6*A, 0, Q(4, 3)*rho**2, Q(4, 3)*rho*nu-W/2,
         16*nu*A-Q(16, 3)*rho*B],
        [2*B, 4*B, Q(4, 3)*rho*nu-W/2, Q(4, 3)*nu**2+W,
         16*rho*B-16*nu*A+Q(32, 3)*nu*B],
        [4*W, 4*W, 16*nu*A-Q(16, 3)*rho*B,
         16*rho*B-16*nu*A+Q(32, 3)*nu*B,
         16*(rho+nu)*W+Q(64, 3)*B**2-64*A*B]])
    laps = sp.Matrix([6*d+4, 6*d, 0, 0,
                      16*d*rho+Q(16, 3)*(3*d+2)*nu])
    f = A+t*B
    s = Q(4, 3)*(rho+t*nu)**2+t*(t-1)*W
    lf = lap(f, variables, table, laps)
    actual_s = gamma(f, f, variables, table)
    actual_L = sp.expand(actual_s*lf-gamma(f, actual_s, variables, table)/2)
    expected_L = (-8*(rho+t**2*(2-t)*nu)*f
                  -Q(8, 3)*t*(t-1)*(3*t+1)*(rho+t*nu)*B)
    zero("Delta(F_t)=0 for formal d", lf, "formal closure")
    zero("S_t full polynomial", actual_s-s, "formal closure")
    zero("L(F_t) full polynomial", actual_L-expected_L, "formal closure")
    for k in range(3):
        zero(f"S_t coefficient of t^{k}",
             sp.expand(actual_s-s).coeff(t, k), "coefficient checks")
    # The displayed form has t^4 terms which must cancel.
    for k in range(5):
        zero(f"L(F_t) coefficient of t^{k}",
             sp.expand(actual_L-expected_L).coeff(t, k), "coefficient checks")
    zero("L(F_t) displayed t^4 cancellation",
         sp.expand(expected_L).coeff(t, 4), "coefficient checks")
    special_weights = {sp.Integer(0): 8*rho,
                       sp.Integer(1): 8*(rho+nu),
                       -Q(1, 3): 8*(rho+Q(7, 27)*nu)}
    for value, weight in special_weights.items():
        zero(f"minimal equation at t={value}",
             actual_L.subs(t, value)+weight*f.subs(t, value),
             "special parameters")
    zero("existing nonradial S at t=-1/3",
         s.subs(t, -Q(1, 3))
         -Q(4, 3)*(rho-nu/3)**2-Q(4, 9)*W, "special parameters")
    zero("existing nonradial Gamma(F,S)=2qF",
         gamma(f, s, variables, table).subs(t, -Q(1, 3))
         -16*(rho+Q(7, 27)*nu)*(A-B/3), "special parameters")
    return variables, t, d, table, laps, f, s, actual_L


def independent_real_coordinates(closure):
    variables, t, d, table, laps, f, s, formal_L = closure
    x, y, a, b, c, u1, u2, u3 = coords = sp.symbols(
        "x y a b c u1 u2 u3", real=True)
    Z = sp.Matrix([[x, c, b], [c, y, a], [b, a, -x-y]])
    u = sp.Matrix([u1, u2, u3])
    # In these rational coordinates 1/2 tr(dZ^2)+|du|^2 has diagonal
    # 2x2 block [[1,1/2],[1/2,1]] and an identity block on the rest.
    inverse_metric = sp.eye(8)
    inverse_metric[:2, :2] = sp.Matrix([[Q(4, 3), -Q(2, 3)],
                                      [-Q(2, 3), Q(4, 3)]])
    generators = [sp.expand(sp.trace(Z*Z)/2),
                  (u.T*u)[0], sp.expand(Z.det()),
                  sp.expand(-(u.T*Z*u)[0]),
                  sp.expand(4*(u.T*Z*Z*u)[0])]
    substitute = dict(zip(variables, generators))
    for i, vi in enumerate(variables):
        zero(f"real-coordinate Delta({vi})",
             lap(generators[i], coords, inverse_metric, [0]*8)
             -laps[i].subs(d, 1).subs(substitute),
             "independent d=1 coordinates")
        for j in range(i, len(variables)):
            zero(f"real-coordinate Gamma({vi},{variables[j]})",
                 gamma(generators[i], generators[j], coords, inverse_metric)
                 -table[i, j].subs(substitute),
                 "independent d=1 coordinates")
    actual_f = generators[2]+t*generators[3]
    actual_s = sp.expand(gamma(actual_f, actual_f, coords, inverse_metric))
    actual_lap = lap(actual_f, coords, inverse_metric, [0]*8)
    # Direct Hessian contraction checks the convention for L independently
    # of the formal Gamma(F,S)/2 calculation.
    gradient_vector = inverse_metric*grad(actual_f, coords)
    actual_L = sp.expand(actual_s*actual_lap
                        -(gradient_vector.T*sp.hessian(actual_f, coords)
                          *gradient_vector)[0])
    zero("real-coordinate Delta(F_t)", actual_lap,
         "independent d=1 coordinates")
    zero("real-coordinate S_t", actual_s-s.subs(substitute),
         "independent d=1 coordinates")
    zero("real-coordinate Hessian definition of L(F_t)",
         actual_L-formal_L.subs(substitute),
         "independent d=1 coordinates")
    return Z, u, coords, actual_f, actual_s, actual_L


def regular_zero_obstruction(closure, coordinate_data):
    variables, t, d, table, laps, f, s, formal_L = closure
    rho, nu, A, B, W = variables
    Z, u, coords, actual_f, actual_s, actual_L = coordinate_data
    x, y, a, b, c, u1, u2, u3 = coords
    tp = sp.Symbol("t_positive", positive=True)
    tm = sp.Symbol("t_negative", negative=True)
    branches = [
        ("t>0", tp, sp.Matrix([0, 0, 1/sp.sqrt(tp)]),
         1/tp, 16/tp, 4, 16*tp, -Q(64, 3)*(tp-1)*(3*tp+1)),
        ("t<0", tm, sp.Matrix([sp.sqrt(-2/tm), 0, 0]),
         -2/tm, -8/tm, 1, -8*tm, -Q(16, 3)*(tm-1)*(3*tm+1))]
    witnesses = []
    Z0 = sp.diag(1, 1, -2)
    for label, tb, ub, nub, wb, rb, norm_u_sq, expected_L in branches:
        witness = {rho: 3, nu: nub, A: -2, B: 2/tb, W: wb, t: tb}
        coord_witness = dict(zip(coords, [1, 1, 0, 0, 0, *ub]))
        coord_witness[t] = tb
        zero(f"{label} witness nu", (ub.T*ub)[0]-nub, "regular zeros")
        zero(f"{label} witness A", Z0.det()+2, "regular zeros")
        zero(f"{label} witness B", -(ub.T*Z0*ub)[0]-2/tb,
             "regular zeros")
        zero(f"{label} witness W", 4*(ub.T*Z0*Z0*ub)[0]-wb,
             "regular zeros")
        zero(f"{label} witness F_t=0", f.subs(witness), "regular zeros")
        zero(f"{label} witness rho+t*nu", (rho+t*nu).subs(witness)-rb,
             "regular zeros")
        derivative_u = -2*tb*Z0*ub
        zero(f"{label} witness |D_u F_t|^2",
             (derivative_u.T*derivative_u)[0]-norm_u_sq, "regular zeros")
        condition(f"{label} witness is regular",
                  norm_u_sq.is_positive is True, "regular zeros")
        zero(f"{label} L(F_t) at witness",
             formal_L.subs(witness)-expected_L, "regular zeros")
        zero(f"{label} direct coordinate F_t=0",
             actual_f.subs(coord_witness), "regular zeros")
        zero(f"{label} direct coordinate L(F_t)",
             actual_L.subs(coord_witness)-expected_L, "regular zeros")
        witnesses.append({
            "parameter_branch": label, "Z": "diag(1,1,-2)",
            "u": str(ub.T), "A": "-2", "B": str(2/tb),
            "rho": "3", "nu": str(nub), "W": str(wb),
            "rho_plus_t_nu": str(rb), "F_t": "0",
            "D_u_F_t_squared_norm": str(norm_u_sq),
            "L_at_regular_zero": str(expected_L)})
    condition("exceptional nonzero roots are exactly -1/3 and 1",
              set(sp.solve((t-1)*(3*t+1), t)) == {-Q(1, 3), sp.Integer(1)},
              "regular zeros")
    # This regular zero works for every real t, including t=0.
    common_zero = {x: 1, y: -1, a: 0, b: 0, c: 0,
                   u1: 0, u2: 0, u3: 0}
    zero("all parameters have common zero diag(1,-1,0),u=0",
         actual_f.subs(common_zero), "regular zeros")
    zero("common zero gradient squared norm equals 4/3",
         actual_s.subs(common_zero)-Q(4, 3), "regular zeros")
    return witnesses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("verification_results.json"))
    args = parser.parse_args()
    started = time.perf_counter()
    closure = formal_closure()
    coordinate_data = independent_real_coordinates(closure)
    witnesses = regular_zero_obstruction(closure, coordinate_data)
    variables, t, d, table, laps, f, s, actual_L = closure
    report = {
        "description": "Exact verification of the Jordan pencil F_t=A+tB",
        "utc_execution_time": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "sympy": sp.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "arithmetic": "exact rational polynomial and symbolic radical identities",
        "elapsed_seconds": round(time.perf_counter()-started, 3),
        "status": "PASS" if all(row["passed"] for row in CHECKS) else "FAIL",
        "passed": sum(row["passed"] for row in CHECKS), "total": len(CHECKS),
        "formulas": {
            "F_t": str(f), "S_t": str(s), "Delta_F_t": "0",
            "L_F_t": "-8*(rho+t**2*(2-t)*nu)*F_t - (8/3)*t*(t-1)*(3*t+1)*(rho+t*nu)*B",
            "expanded_L_coefficients_in_t": {
                str(k): str(sp.factor(actual_L.coeff(t, k))) for k in range(4)},
            "minimal_parameters": ["0", "1", "-1/3"],
            "q_at_minus_one_third": "8*(rho+7*nu/27)"},
        "regular_zero_witnesses": witnesses,
        "scope": [
            "The formal calculation assumes the stated invariant closure for d=1,2,4.",
            "The full closure and Laplacian table are independently verified in real coordinates for d=1.",
            "The diagonal real witnesses embed unchanged into K=R,C,H.",
            "The necessity argument uses L nonzero at a regular zero for every excluded real t.",
            "The geometric equivalence of L=0 on regular zeros with minimality is used, not formally proved by this script.",
            "The Cayley-Dickson identification N=A+B, complex/quaternionic coordinate closure, irreducibility, stability, and minimization are outside this check."],
        "checks": CHECKS}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(f"{report['status']}: {report['passed']}/{report['total']} exact checks")
    print(f"Report: {args.output.resolve()}")
    if report["status"] != "PASS":
        for row in CHECKS:
            if not row["passed"]:
                print(json.dumps(row))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
