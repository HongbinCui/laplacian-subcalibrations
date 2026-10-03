#!/usr/bin/env python3
"""Exact independent coefficient checks; standard library only, no CAS.

Sparse Laurent polynomials have Fraction exponents and coefficients.
The general source powers R^a are factored before calculation so that a
remains an indeterminate. This file does not emulate or run Mathematica.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
import platform
import sys

NAMES = 'r h U D L a n k ell wh dh lh M q f'.split()
DIM = len(NAMES)
ZERO = (Q(0),) * DIM


class Poly:
    def __init__(self, data=0):
        if isinstance(data, Poly):
            self.terms = dict(data.terms)
        elif isinstance(data, dict):
            self.terms = {t: Q(c) for t, c in data.items() if c}
        else:
            self.terms = {} if not data else {ZERO: Q(data)}

    def __add__(self, other):
        out = dict(self.terms)
        for t, c in Poly(other).terms.items():
            out[t] = out.get(t, Q(0)) + c
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({t: -c for t, c in self.terms.items()})

    def __sub__(self, other):
        return self + -Poly(other)

    def __rsub__(self, other):
        return Poly(other) + -self

    def __mul__(self, other):
        out = {}
        for t, c in self.terms.items():
            for s, d in Poly(other).terms.items():
                v = tuple(x+y for x, y in zip(t, s))
                out[v] = out.get(v, Q(0)) + c*d
        return Poly(out)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        exponent = Q(exponent)
        if exponent.denominator == 1 and exponent >= 0:
            out = Poly(1)
            for _ in range(int(exponent)):
                out *= self
            return out
        if len(self.terms) != 1:
            raise ValueError('Nonintegral or negative power requires monomial')
        (t, c), = self.terms.items()
        if c != 1:
            raise ValueError('Monomial coefficient must be 1')
        return Poly({tuple(x*exponent for x in t): 1})

    def diff(self, name):
        i = NAMES.index(name)
        out = {}
        for t, c in self.terms.items():
            if t[i]:
                v = list(t)
                v[i] -= 1
                out[tuple(v)] = c*t[i]
        return Poly(out)

    def subs(self, replacements):
        out = Poly()
        for t, c in self.terms.items():
            monomial = Poly(c)
            for name, exponent in zip(NAMES, t):
                if exponent:
                    monomial *= replacements.get(name, VARS[name]) ** exponent
            out += monomial
        return out


VARS = {}
for index, name in enumerate(NAMES):
    exponent = list(ZERO)
    exponent[index] = Q(1)
    VARS[name] = Poly({tuple(exponent): 1})
globals().update(VARS)
results = []


def check(label, actual, expected):
    remainder = actual - expected
    passed = not remainder.terms
    results.append({'identity': label, 'passed': passed,
                    'nonzero_remainder_terms': len(remainder.terms)})
    if not passed:
        raise AssertionError(label + ': nonzero exact remainder')


# Generic homogeneous source u=h of degree k, U=|grad h|^2,
# D=Delta h, L=L(h). These are only the differential data needed.
metric = [[4*r, 2*k*h, 4*(k-1)*U],
          [2*k*h, U, 2*(U*D-L)],
          [4*(k-1)*U, 2*(U*D-L), Poly(0)]]


def pair(left, right, g):
    return sum((left[i]*g[i][j]*right[j]
                for i in range(len(g)) for j in range(len(g))), Poly())


# grad(R^a h)/R^a in coordinates (R,h,U).
grad_g = [a*h*r**-1, Poly(1), Poly(0)]
norm_bar = pair(grad_g, grad_g, metric)
norm_expected = U + 4*a*(a+k)*r**-1*h**2
check('radial master: gradient square / R^(2a)', norm_bar, norm_expected)

lap_bar = D + 2*a*(n+2*a+2*k-2)*r**-1*h
grad_norm_bar = [norm_bar.diff('r') + 2*a*r**-1*norm_bar,
                 norm_bar.diff('h'), norm_bar.diff('U')]
defect_bar = norm_bar*lap_bar-Q(1, 2)*pair(grad_g, grad_norm_bar, metric)
master_expected = (L + 2*a*(n-2*k-2*a-1)*r**-1*h*U
                   +4*a*(k+a)*r**-1*h**2*D
                   +4*a*(k+a)*(2*a*(n-1)+k)*r**-2*h**3)
check('radial master: defect / R^(3a)', defect_bar, master_expected)

# Independent pullback chain rule, using ell-degree H and W_H.
pullback_metric = [[4*r, 4*ell*h, 4*(2*ell-2)*wh],
                   [4*ell*h, 4*r*wh, 8*r*(wh*dh-lh)],
                   [4*(2*ell-2)*wh, 8*r*(wh*dh-lh), Poly(0)]]
pullback_gradient = [Poly(0), Poly(1), Poly(0)]
pullback_norm = 4*r*wh
pullback_lap = 4*r*dh
norm_diff = [pullback_norm.diff(v) for v in ['r', 'h', 'wh']]
pullback_defect = pullback_norm*pullback_lap-Q(1, 2)*pair(
    pullback_gradient, norm_diff, pullback_metric)
check('unweighted Clifford pullback: defect', pullback_defect,
      16*r**2*lh-8*ell*h*wh)

pullback = {'k': 2*ell, 'U': 4*r*wh, 'D': 4*r*dh,
            'L': 16*r**2*lh-8*ell*h*wh}
b = a*(a+2*ell)
A = a*(n-1)-ell-2*b
source_norm = norm_bar.subs(pullback)
source_lap = lap_bar.subs(pullback)
source_defect = defect_bar.subs(pullback)
check('weighted Clifford pullback: gradient / R^(2a)', source_norm,
      4*r*wh+4*b*r**-1*h**2)
check('weighted Clifford pullback: Laplacian / R^a', source_lap,
      4*r*dh+2*a*(n+2*a+4*ell-2)*r**-1*h)
check('weighted Clifford pullback: defect / R^(3a)', source_defect,
      16*r**2*lh+8*A*h*wh+16*b*h**2*dh
      +8*b*(a*(n-1)+ell)*r**-2*h**3)
check('weighted Clifford pullback: preserved M-ratio / R^(2a)',
      h*source_lap+M*source_norm,
      4*r*(h*dh+M*wh)
      +(2*a*(n+2*a+4*ell-2)+4*M*a*(a+2*ell))*r**-1*h**2)

# Here k denotes a second source dimension N'; all target data are fixed.
check('dimension monotonicity: finite difference / R^(4a)',
      h*(source_defect.subs({'n': k})-source_defect),
      2*a*(k-n)*r**-1*h**2*source_norm)

# Independent target calculation in the exact two-variable quotient.
# q=|u|^2+|v|^2, f=|u|^2-|v|^2 on R^k+R^k.
target_metric = [[4*q, 4*f], [4*f, 4*q]]


def target_gamma(v, w):
    return pair([v.diff('q'), v.diff('f')],
                [w.diff('q'), w.diff('f')], target_metric)


def target_lap(v):
    return (4*q*(v.diff('q').diff('q')+v.diff('f').diff('f'))
            +8*f*v.diff('q').diff('f')+4*k*v.diff('q'))


target_h = q**Q(3, 4)*f
target_wh = target_gamma(target_h, target_h)
target_dh = target_lap(target_h)
target_lh = target_wh*target_dh-Q(1, 2)*target_gamma(target_h, target_wh)
check('Simons fractional target: W_H', target_wh,
      4*q**Q(5, 2)+Q(33, 4)*q**Q(1, 2)*f**2)
check('Simons fractional target: Delta H', target_dh,
      Q(3, 4)*(4*k+7)*q**Q(-1, 4)*f)
check('Simons fractional target: L(H)', target_lh,
      q**Q(1, 4)*f*((12*k-47)*q**2+Q(33, 8)*(6*k+1)*f**2))
check('Simons fractional target: H Delta H', target_h*target_dh,
      Q(3, 4)*(4*k+7)*q**Q(1, 2)*f**2)
dimension_bound = 1+4*ell*(1+M)+2*a*(1+M)+ell*a**-1
check('Simons fractional target: a=1, ell=7/2, M=0 dimension bound',
      dimension_bound.subs({'ell': Poly(Q(7, 2)), 'M': Poly(0), 'a': Poly(1)}),
      Poly(Q(41, 2)))

report = {
    'method': 'Exact sparse Laurent polynomial arithmetic over fractions.Fraction',
    'engine': 'Python standard library; no SymPy, Mathics, or Mathematica used',
    'python_version': platform.python_version(),
    'native_mathematica_executed': False,
    'passed': all(item['passed'] for item in results),
    'identity_count': len(results),
    'checks': results,
    'scope': 'Differential identities only; does not verify geometric measure theory arguments.',
}
output = Path(__file__).with_name('exact_results.json')
output.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
sys.exit(0 if report['passed'] else 1)
