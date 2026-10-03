"""Exact checks for four-generator exceptional cubic subcalibrations.

Requires Python 3 and SymPy.  Run this file directly.  The computations check
the stated differential identities; the global inequalities, critical-set
geometry away from the displayed rank-one point, and removable-singularity
argument are proved in the manuscript, not inferred from these calculations.
"""
import datetime
import itertools
import platform
import sympy as sp

print("UTC:", datetime.datetime.now(datetime.timezone.utc).isoformat())
print("Python:", platform.python_version(), "SymPy:", sp.__version__)
passed = []


def check(label, residual):
    entries = list(residual) if isinstance(residual, sp.MatrixBase) else [residual]
    values = [sp.factor(sp.cancel(sp.expand(value))) for value in entries]
    if any(value != 0 for value in values):
        raise AssertionError((label, values))
    passed.append(label)
    print("PASS:", label)


# The free variables have degrees 2, 4, 6, 3, respectively.
R, S = sp.symbols("R S", positive=True)
D, F, c, d, dr, dd = sp.symbols("D F c d dr dd", real=True)
variables = (R, S, D, F)
metric = sp.Matrix([
    [4*R, 8*S, 12*D, 6*F],
    [8*S, 4*R*S+12*D, 8*R*D, 4*R*F],
    [12*D, 8*R*D, 4*D*S, 2*S*F],
    [6*F, 4*R*F, 2*S*F, S],
])


def gradient(expression, coordinates=variables):
    return sp.Matrix([sp.diff(expression, x) for x in coordinates])


def gamma(left, right):
    return (gradient(left).T*metric*gradient(right))[0]


def lap(expression, drift):
    return (gradient(expression).T*sp.Matrix(drift))[0] + sum(
        metric[i, j]*sp.diff(expression, variables[i], variables[j])
        for i in range(4) for j in range(4))


def defect(expression, drift):
    q = sp.expand(gamma(expression, expression))
    return q*lap(expression, drift)-gamma(expression, q)/2


G = S*F
Q = S**3+12*F**2*(R*S+D)
P = ((c-16)*R*S**3-30*D*S**2
     + F**2*((12*c-72)*R*D+12*c*R**2*S-96*S**2))
positive = ((c-sp.Rational(58, 3))*R*S**3
            + sp.Rational(10, 3)*S**2*(R*S-9*D)
            + F**2*((12*c-72)*R*D+12*c*S*(R**2-3*S)
                    +(36*c-96)*S**2))
check("SF gradient norm", gamma(G, G)-Q)
check("SF Laplacian", lap(G, (dr, c*R, dd, 0))-(c+8)*R*F)
check("SF defect, with arbitrary Delta R and Delta D",
      defect(G, (dr, c*R, dd, 0))-F*P)
check("SF positive decomposition", P-positive)
check("SF threshold at balanced spectrum and F=0",
      P.subs({R: 3, S: 3, D: 1, F: 0})-81*(c-sp.Rational(58, 3)))

Groot = F*sp.sqrt(R*S)
Qroot = R*S**2+F**2*(3*D*R/S+5*R**2+11*S)
Pbase = (6*R*S**2*(S**2-2*R*D)+15*D*F**2*R*(R**2-S)
         +45*F**2*R**4*S+123*F**2*R**2*S**2+242*F**2*S**3)
check("sqrt(RS)F gradient norm", gamma(Groot, Groot)-Qroot)
check("sqrt(RS)F certificate for Jordan d=2",
      defect(Groot, (36, 16*R, 4*S, 0))*R*S**2-Groot*Pbase)
transfer = (d-2)*Qroot*(12*sp.diff(Groot, R)+4*R*sp.diff(Groot, S))
check("Jordan drift transfer",
      defect(Groot, (12*(d+1), 4*(d+2)*R, 4*S, 0))
      -defect(Groot, (36, 16*R, 4*S, 0))-transfer)
check("Sign of the transfer bracket",
      12*sp.diff(Groot, R)+4*R*sp.diff(Groot, S)-Groot*(6/R+2*R/S))

# Independent Gram-coordinate computation for an n by 3 real matrix X.
# B_ij = <X_i,X_j>; Delta B_ii=2n; Delta B_ij=0 for i != j.
n = sp.symbols("n", integer=True, positive=True)
gramvars = sp.symbols("b11 b22 b33 b12 b13 b23", real=True)
u, v, w, p, q, t = gramvars
B = sp.Matrix([[u, p, q], [p, v, t], [q, t, w]])
pairs = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
grammetric = sp.Matrix(6, 6, lambda alpha, beta: sum([
    int(pairs[alpha][0] == pairs[beta][0])*B[pairs[alpha][1], pairs[beta][1]],
    int(pairs[alpha][0] == pairs[beta][1])*B[pairs[alpha][1], pairs[beta][0]],
    int(pairs[alpha][1] == pairs[beta][0])*B[pairs[alpha][0], pairs[beta][1]],
    int(pairs[alpha][1] == pairs[beta][1])*B[pairs[alpha][0], pairs[beta][0]],
]))
gramdrift = sp.Matrix([2*n, 2*n, 2*n, 0, 0, 0])
Rg, Sg, Dg = sp.trace(B), u*v+u*w+v*w-p*p-q*q-t*t, B.det()


def gramgamma(left, right):
    return (gradient(left, gramvars).T*grammetric*gradient(right, gramvars))[0]


def gramlap(expression):
    return (gradient(expression, gramvars).T*gramdrift)[0]+sum(
        grammetric[i, j]*sp.diff(expression, gramvars[i], gramvars[j])
        for i in range(6) for j in range(6))


check("Gram Delta R", gramlap(Rg)-6*n)
check("Gram Delta S", gramlap(Sg)-4*(n-1)*Rg)
check("Gram Delta D", gramlap(Dg)-2*(n-2)*Sg)
for label, left, right, expected in [
    ("RR", Rg, Rg, 4*Rg), ("RS", Rg, Sg, 8*Sg),
    ("RD", Rg, Dg, 12*Dg), ("SS", Sg, Sg, 4*Rg*Sg+12*Dg),
    ("SD", Sg, Dg, 8*Rg*Dg), ("DD", Dg, Dg, 4*Sg*Dg),
]:
    check("Gram Gamma "+label, gramgamma(left, right)-expected)
check("G2 Delta D = 10S", gramlap(Dg).subs(n, 7)-10*Sg)
# For any alternating trilinear F, Gamma(F,B_ij)=2F delta_ij.
check("Alternating F: Gamma(F,S)=4RF",
      2*F*sum(sp.diff(Sg, variable) for variable in gramvars[:3])-4*Rg*F)
check("Alternating F: Gamma(F,D)=2SF",
      2*F*sum(sp.diff(Dg, variable) for variable in gramvars[:3])-2*Sg*F)

# Compact G2 form: 123+145+167+246-257-347-356.
triples = [((0, 1, 2), 1), ((0, 3, 4), 1), ((0, 5, 6), 1),
           ((1, 3, 5), 1), ((1, 4, 6), -1), ((2, 3, 6), -1),
           ((2, 4, 5), -1)]
a = sp.Matrix(sp.symbols("a1:8", real=True))
b = sp.Matrix(sp.symbols("b1:8", real=True))
z = sp.Matrix(sp.symbols("z1:8", real=True))
coordinates = tuple(a)+tuple(b)+tuple(z)
FG = 0
for triple, sign in triples:
    for permutation in itertools.permutations(range(3)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(3) for j in range(i+1, 3))
        i, j, k = (triple[index] for index in permutation)
        FG += sign*(-1)**inversions*a[i]*b[j]*z[k]
SG = ((a.dot(a)*b.dot(b)-a.dot(b)**2)
      +(a.dot(a)*z.dot(z)-a.dot(z)**2)
      +(b.dot(b)*z.dot(z)-b.dot(z)**2))
gradF = gradient(FG, coordinates)
check("G2 coordinate gradient norm = Gram e2", gradF.dot(gradF)-SG)
check("G2 harmonicity", sum(sp.diff(FG, x, x) for x in coordinates))

# A concrete four-parameter slice: a=r e1, b=s e2, z=u e3+v e4.
rr, ss, uu, vv = sp.symbols("r s u v", real=True)
restricted = sp.Matrix([
    rr**2+ss**2+uu**2+vv**2,
    rr**2*ss**2+(rr**2+ss**2)*(uu**2+vv**2),
    rr**2*ss**2*(uu**2+vv**2), rr*ss*uu])
jac = restricted.jacobian((rr, ss, uu, vv)).det()
sample = jac.subs({rr: 1, ss: 2, uu: 3, vv: 4})
if sample == 0:
    raise AssertionError("G2 algebraic independence witness vanishes")
print("G2 nonzero Jacobian witness:", sample)
check("G2 independence: explicit Jacobian formula",
      jac+8*rr**2*ss**2*vv*(rr**2-ss**2)
      *(rr**2-uu**2-vv**2)*(ss**2-uu**2-vv**2))

# At (r e1,0,0), the Hessian has eigenvalues +r and -r, each six times.
rankone = {x: 0 for x in coordinates}
rankone[a[0]] = rr
H = sp.hessian(FG, coordinates).subs(rankone)
check("G2 rank-one Hessian cubic identity", H**3-rr**2*H)
check("G2 rank-one Hessian trace", sp.trace(H))
check("G2 rank-one Hessian squared trace", sp.trace(H**2)-12*rr**2)
lam = sp.symbols("lambda")
check("G2 rank-one Hessian characteristic polynomial",
      H.charpoly(lam).as_expr()-lam**9*(lam**2-rr**2)**6)
print("All", len(passed), "exact checks passed.")
print("These computations do not constitute a Wolfram Mathematica execution.")
