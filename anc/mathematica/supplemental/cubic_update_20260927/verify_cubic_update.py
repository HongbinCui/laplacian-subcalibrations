#!/usr/bin/env python3
"""Exact coefficient checks for the September 27 cubic revision.

Run with Python 3 and SymPy. Every reported check tests a rational
identity, after clearing denominators. The geometric inequalities,
removability arguments, and spectral classification require the proofs
in the paper; this script does not substitute numerical sampling for them.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import sympy as sp

Q = sp.Rational
CHECKS = []


def zero(name, expression):
    numerator = sp.fraction(sp.cancel(expression))[0]
    residual = sp.expand(numerator)
    if residual != 0:
        raise AssertionError(f"FAIL {name}: {residual}")
    CHECKS.append(name)
    print(f"PASS {len(CHECKS):02d} {name}", flush=True)


def grad(expression, variables):
    return sp.Matrix([sp.diff(expression, x) for x in variables])


def gamma(a, b, variables, table):
    return (grad(a, variables).T * table * grad(b, variables))[0]


def lap(expression, variables, table, laps):
    return sum(sp.diff(expression, x) * d for x, d in zip(variables, laps)) + sum(
        table[i, j] * sp.diff(expression, x, y)
        for i, x in enumerate(variables) for j, y in enumerate(variables))


def unified_radial():
    r, f, s, d, n, k, lapd = sp.symbols("r f s d n k lapd")
    variables = [r, f, s, d]
    table = sp.Matrix([[4*r, 6*f, 8*s, 12*d],
                       [6*f, s, 4*r*f, 2*s*f],
                       [8*s, 4*r*f, 4*r*s+12*d, 8*r*d],
                       [12*d, 2*s*f, 8*r*d, 4*s*d]])
    # l = grad(log G), G^2 = r*s*f^2. These formulas avoid radicals.
    ell = sp.Matrix([1/(2*r), 1/f, 1/(2*s), 0])
    q = r*s**2 + f**2*(3*d*r/s + 5*r**2 + 11*s)
    dl = (n+9)/r + k/(2*s) + 3*r/s - 3*d/s**2
    q_from_chain = r*s*f**2*(ell.T*table*ell)[0]
    dl_from_chain = (ell.T*sp.Matrix([2*n, 0, k, lapd]))[0] + sum(
        table[i, j]*(sp.diff(ell[i], variables[j])+ell[i]*ell[j])
        for i in range(4) for j in range(4))
    zero("unified gradient norm", q_from_chain-q)
    zero("unified Laplacian ratio", dl_from_chain-dl)
    hh = (n-12)*s**2 + r*(k-16*r)*s/2 - 12*r*d
    jj = (r*d*(Q(3, 2)*r*(k-16*r)+15*(r**2-3*s)+3*(n-8)*s)
          + Q(5, 2)*(k-16*r)*r**3*s + 45*r**4*s
          + Q(11, 2)*(k-16*r)*r*s**2
          + (5*n+33)*r**2*s**2 + (11*n+44)*s**3)
    m_over_g = q*dl - (ell.T*table*grad(q, variables))[0]/2
    zero("unified H/J defect identity", r*s**2*m_over_g-r*s**2*hh-f**2*jj)
    zero("trace corollary remainder", hh-(n-16)*s**2-r*(k-16*r)*s/2-4*(s**2-3*r*d))


def clifford():
    u, v, tt, ff, m, q = sp.symbols("u v tt ff m q")
    variables = [u, v, tt, ff]
    # ff is the manuscript cubic F, not its normalized half.
    table = sp.Matrix([[4*u, 0, 8*tt, 4*ff], [0, 4*v, 0, 2*ff],
                       [8*tt, 0, 16*u*tt, 8*u*ff],
                       [4*ff, 2*ff, 8*u*ff, 4*u*v+tt]])
    r, f, s, d = u+v, ff/2, u*v+tt/4, v*tt/4
    k = 4*m*v+4*(q+1)*u
    for a, b, expected, name in [
        (f, s, 4*r*f, "f,s"), (s, s, 4*r*s+12*d, "s,s"),
        (f, d, 2*s*f, "f,D"), (s, d, 8*r*d, "s,D"),
        (d, d, 4*s*d, "D,D")]:
        zero("Clifford Gamma("+name+")", gamma(a, b, variables, table)-expected)
    zero("Clifford normalized Laplacian", lap(s, variables, table, [4*m, 2*(q+1), 8*(q+1)*u, 0])-k)
    hh = (2*m+q-11)*s**2+r*(k-16*r)*s/2-12*r*d
    h44 = tt*(8*u**2+tt)/16+(u**2-tt)*v*(2*u+3*v)
    zero("Clifford H44 decomposition", hh.subs({m: 4, q: 4})-h44)
    zero("Clifford parameter increment", hh-h44-(2*(m-4)+q-4)*s**2-2*r*((m-4)*v+(q-4)*u)*s)
    # The shortened all-power obstruction, with the original S normalization.
    lam, al, ga = sp.symbols("lam al ga")
    sl = lam*(4-3*lam)
    # jl denotes current J = |grad S|^2 on the normalized edge.
    jl = 16*(4*lam*(1-lam)+(4-3*lam)*lam**2)
    hl = 32*lam*(1-lam)**2+32*lam**2*(1-lam)+8*lam**3-16*(1-lam)*lam**2
    boundary = 64*ga-8-ga*(ga+2)*jl/sl-ga*hl/sl+2*al*(5-8*ga-2*al)*sl
    zero("Clifford power obstruction endpoint", sp.limit(boundary, lam, 0)+8*(ga-1)*(2*ga-1))
    zero("Clifford power obstruction interior", boundary.subs(lam, Q(2, 3))+Q(8, 3)*(2*al**2+(8*ga-5)*al+8*ga**2-6*ga+3))


def radial_lift():
    r, h, w, ell, z, d, n = sp.symbols("r h w ell z d n")
    variables = [r, h, w]
    table = sp.Matrix([[4*r, 2*d*h, 4*(d-1)*w],
                       [2*d*h, w, z], [4*(d-1)*w, z, 0]])
    lg = sp.Matrix([1/(2*r), 1/h, 0])
    qg = r*w+(2*d+1)*h**2
    zero("radial lift gradient identity", r*h**2*(lg.T*table*lg)[0]-qg)
    dl = ell/h+(2*d+n-1)/r
    actual = qg*dl-(lg.T*table*grad(qg, variables))[0]/2
    expected = r/h*(w*ell-z/2+(2*d+1)*h**2*ell/r
                         +(n-2*d-2)*h*w/r+(2*d+1)*(d+n-1)*h**3/r**2)
    zero("radial lift defect identity", actual-expected)
    zero("dimension 26 radial lift coefficient 11", (2*d+1).subs(d, 5)-11)
    zero("dimension 26 radial lift coefficient 14", (n-2*d-2).subs({n: 26, d: 5})-14)
    zero("dimension 26 radial lift coefficient 330", ((2*d+1)*(d+n-1)).subs({n: 26, d: 5})-330)


def nonradial_26():
    r, nu, a, b, w = sp.symbols("rho nu a b w")
    variables = [r, nu, a, b, w]
    table = sp.Matrix([
        [4*r, 0, 6*a, 2*b, 4*w],
        [0, 4*nu, 0, 4*b, 4*w],
        [6*a, 0, Q(4, 3)*r**2, Q(4, 3)*r*nu-w/2, 16*nu*a-Q(16, 3)*r*b],
        [2*b, 4*b, Q(4, 3)*r*nu-w/2, Q(4, 3)*nu**2+w, 16*r*b-16*nu*a+Q(32, 3)*nu*b],
        [4*w, 4*w, 16*nu*a-Q(16, 3)*r*b, 16*r*b-16*nu*a+Q(32, 3)*nu*b,
         16*(r+nu)*w+Q(64, 3)*b**2-64*a*b]])
    laps = [28, 24, 0, 0, Q(16, 3)*(12*r+14*nu)]
    gm = lambda x, y: gamma(x, y, variables, table)
    f = a-b/3
    p = Q(4, 3)*(r-nu/3)**2+Q(4, 9)*w
    q = 8*(r+Q(7, 27)*nu)
    k = Q(32, 27)*(78*r+14*nu)
    tt = Q(64, 27)*(9*r+nu)*p-Q(1024, 81)*f*b
    zero("nonradial Gamma(F,F)", gm(f, f)-p)
    zero("nonradial harmonicity", lap(f, variables, table, laps))
    zero("nonradial minimality", gm(f, p)-2*q*f)
    zero("nonradial Delta(P)", lap(p, variables, table, laps)-k)
    zero("nonradial Gamma(P,P)", gm(p, p)-tt)
    logh = grad(f, variables)/f+grad(p, variables)/(2*p)
    qh = p**2+2*q*f**2+f**2*tt/(4*p)
    ll = Q(512, 9)*r+Q(320, 27)*nu+Q(256, 81)*f*b/p
    zero("nonradial H gradient", f**2*p*(logh.T*table*logh)[0]-qh)
    zero("nonradial H Laplacian", 2*q+k/2-tt/(4*p)-ll)
    nn = ((405*r+81*nu)*p**3+144*f*b*p**2
          +f**2*p*(46692*r**2+20040*r*nu+2212*nu**2-4239*p)
          -(4992*r+896*nu)*f**3*b-256*nu*f**4)
    mh_over_h = qh*ll/p-(logh.T*table*grad(qh, variables))[0]/2
    zero("nonradial exact N defect", 729*p**2*mh_over_h-16*nn)
    # The positive decomposition is a polynomial identity in formal P,F,B.
    pp, ff = sp.symbols("pp ff")
    nformal = ((405*r+81*nu)*pp**3+144*ff*b*pp**2
               +ff**2*pp*(46692*r**2+20040*r*nu+2212*nu**2-4239*pp)
               -(4992*r+896*nu)*ff**3*b-256*nu*ff**4)
    pmax = Q(4, 3)*r**2+Q(40, 27)*r*nu+Q(4, 27)*nu**2
    h0 = Q(8, 9)*(36639*r**2+12634*r*nu+1555*nu**2)
    k0 = 4992*r+896*nu
    terms = [405*r*pp**3, 36*pp**2*(ff+2*b)**2,
             pp**2*(81*nu*pp-144*b**2), ff**2*pp*h0,
             k0*ff**2*(Q(3, 16)*(9*r+nu)*pp-ff*b),
             256*nu*ff**2*((r+nu)*pp/9-ff**2),
             4275*ff**2*pp*(pmax-pp)]
    zero("nonradial seven-summand decomposition", nformal-sum(terms))
    zero("nonradial rho=0 strictness reduction", nformal.subs({r: 0, ff: 0, b: 0})-81*nu*pp**3)
    zero("nonradial P maximum remainder", pmax-p-Q(4, 9)*(Q(16, 3)*r*nu-w))
    zero("nonradial square bound for H Delta(H)", (r+2*nu)**2-4*nu*(r+nu)-r**2)
    lower = Q(512, 9)*r+Q(320, 27)*nu-Q(32, 81)*(r+2*nu)
    zero("nonradial H Delta(H) lower coefficients", lower-(4576*r+896*nu)/81)


def nonharmonic():
    # h2 = H_2 = |Hess F|^2, k = K = Delta S, tt = T = Gamma(S,S).
    r, f, s, t, q, k, tt, h2, n = sp.symbols("r f s t q k tt h2 n")
    gft, gfq, gftt, gst, gsq, gstt = sp.symbols("gft gfq gftt gst gsq gstt")
    variables = [f, s, t, q, tt]
    fp = 2*s*t+2*q*f
    gamma_f = lambda expr: sum(sp.diff(expr, x)*v for x, v in zip(variables, [s, fp, gft, gfq, gftt]))
    gamma_s = lambda expr: sum(sp.diff(expr, x)*v for x, v in zip(variables, [fp, tt, gst, gsq, gstt]))
    qh = s**2+2*f*s*t+f**2*(2*q+tt/(4*s))
    dl = 3*t/f+(2*q+k/2-tt/(4*s))/s
    zero("nonharmonic H gradient chain rule", s*(s+f*fp/s+f**2*tt/(4*s**2))-qh)
    zero("nonharmonic H Laplacian chain rule", t/f+k/(2*s)+fp/(f*s)-tt/(4*s**2)-dl)
    b0 = 16*s**2*(s*(h2-2*q+3*t**2)-tt)
    b1 = -16*s**2*gfq-8*s**2*gst+16*k*s**2*t+80*s**2*q*t-4*s*tt*t-2*s*gftt
    b2 = -8*s*gsq+16*k*q*s+2*k*tt+32*q**2*s-gstt
    bb = b0+f*b1+f**2*b2
    mh_over_h = qh*dl-gamma_f(qh)/(2*f)-gamma_s(qh)/(4*s)
    zero("nonharmonic B identity with Bochner correction", (16*s**2*mh_over_h-bb).subs(gft, k/2-h2))
    c0 = r**2*b0+16*(n-12)*r*s**4
    c1 = r**2*b1+16*(2*n+9)*r*s**3*t
    c2 = (r**2*b2+r*(88*k*s**2+32*(n-1)*q*s**2+4*(n-23)*s*tt)
          +176*(n+4)*s**3)
    from_lift = (r**2*bb+176*r*(3*f*s**3*t+f**2*s**2*(2*q+k/2-tt/(4*s)))
                 +16*r*(n-12)*s**2*qh+176*(n+4)*f**2*s**3)
    zero("nonharmonic C identity", from_lift-c0-f*c1-f**2*c2)
    hgen = (n-12)*s**2+r*s*(h2-2*q+3*t**2)-r*tt
    zero("nonharmonic leading coefficient", c0-16*r*s**2*hgen)
    d = sp.Symbol("d")
    hrad = (n-12)*s**2+r*(k-16*r)*s/2-12*r*d
    zero("nonharmonic to radial leading term", hgen.subs({t: 0, q: 2*r, tt: 4*r*s+12*d, h2: k/2})-hrad)


def main():
    print("Exact cubic verification revision 20260927", flush=True)
    print("Python", platform.python_version(), "SymPy", sp.__version__, flush=True)
    unified_radial()
    clifford()
    radial_lift()
    nonradial_26()
    nonharmonic()
    print(f"SUCCESS: {len(CHECKS)} exact assertions passed.", flush=True)
    here = Path(__file__).resolve().parent
    metadata = {
        "status": "passed", "exact_assertions": len(CHECKS), "checks": CHECKS,
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(), "sympy": sp.__version__,
        "native_mathematica_run": False,
        "scope": "Exact identities conditional on the displayed differential tables; not a proof of the tables, inequalities, spectral classification, or GMT comparison.",
        "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in [Path(__file__).resolve(), here/"Verify_Cubic_Update_20260927.wl"] if p.exists()}}
    (here/"python_execution_metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")


if __name__ == "__main__":
    main()
