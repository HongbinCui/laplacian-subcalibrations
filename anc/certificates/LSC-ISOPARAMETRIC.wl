(* ::Package:: *)
(* LSC-ISOPARAMETRIC -- Isoparametric.
   Compatible with Mathematica 12. No external package or file dependencies.
   Get returns checks; only the central runner writes execution records.

   The positivity tables accompanying these identities are checked separately in LSC-INTERVAL-POSITIVITY.
*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`Isoparametric`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`Isoparametric`Private`*"];

VerificationLog = {};
EvaluationMessages = {};
EvaluationOK = True;
VerificationVerbose = False;

(* Real powers occur only on the positive-base domains stated in the paper.
   In particular, the defining polynomial F is not assumed positive.
   Matrix equalities count as one assertion, with every entry tested. *)
AlgebraicResidual[expr_] := Module[{residual},
  residual = Factor[Together[PowerExpand[Expand[expr]]]];
  If[residual === 0, 0, Simplify[residual]]
];
CheckEqual[lhs_, rhs_, label_String] := Module[{residuals, passed},
  residuals = AlgebraicResidual /@ Flatten[{lhs - rhs}];
  passed = And @@ (TrueQ[# === 0] & /@ residuals);
  AppendTo[VerificationLog,
    {label, passed, If[passed, "0", ToString[residuals, InputForm]]}];
  If[passed,
    If[TrueQ[VerificationVerbose], Print["PASS: ", label]],
    Print["FAIL: ", label, "; residual = ", InputForm[residuals]]];
  passed
];
CheckTrue[condition_, label_String] := Module[{passed = TrueQ[condition]},
  AppendTo[VerificationLog,
    {label, passed, If[passed, "True", ToString[condition, InputForm]]}];
  If[passed,
    If[TrueQ[VerificationVerbose], Print["PASS: ", label]],
    Print["FAIL: ", label, "; condition = ", InputForm[condition]]];
  passed
];
RunVerificationGroup[label_String, fun_Symbol] := Module[{status, evaluationResult, messages},
  status = Block[{System`$MessageList = {}},
    evaluationResult = CheckAbort[Catch[fun[]; True, _, ($Failed &)], $Aborted];
    messages = ToString[#, InputForm] & /@ System`$MessageList;
    evaluationResult === True && Length[messages] == 0
  ];
  If[!TrueQ[status],
    EvaluationOK = False;
    AppendTo[EvaluationMessages, <|"Group" -> label, "Messages" -> messages,
      "Status" -> ToString[evaluationResult, InputForm]|>]];
  If[evaluationResult === $Aborted, Abort[]];
];

(* Exact algorithm helpers. No execution or global clearing. *)

GammaExpr[p_, q_, vars_List, matrix_List] := Module[{dp, dq},
  dp = Table[D[p, vars[[i]]], {i, Length[vars]}];
  dq = Table[D[q, vars[[i]]], {i, Length[vars]}];
  dp . matrix . dq
];

LaplaceExpr[p_, vars_List, matrix_List, drift_List] :=
  Total[Table[drift[[i]] D[p, vars[[i]]], {i, Length[vars]}]] +
  Total[Flatten[Table[matrix[[i, j]] D[p, vars[[i]], vars[[j]]],
    {i, Length[vars]}, {j, Length[vars]}]]];

CurvatureExpr[p_, vars_List, matrix_List, drift_List] := Module[{v},
  v = GammaExpr[p, p, vars, matrix];
  Expand[v LaplaceExpr[p, vars, matrix, drift] -
    GammaExpr[p, v, vars, matrix]/2]
];

IsoNumerator[g_, d_, t_, ell_, h_, z_] := Module[{aa, bb, energy},
  aa = g^2 (1 - (z + t)^2);
  bb = -g^2 t - g (g + d) z;
  energy = ell^2 h^2 + aa D[h, z]^2;
  Expand[(aa D[h, {z, 2}] + bb D[h, z] + (d + 1) ell h) energy -
    aa D[h, z] D[energy, z]/2]
];

BernsteinCoefficientsExact[polynomial_, variable_, lo_, hi_, degree_Integer] :=
  Module[{tau, transformed, actualDegree},
    transformed = Expand[polynomial /. variable -> lo + (hi - lo) tau];
    actualDegree = Exponent[transformed, tau];
    If[degree < actualDegree,
      Print["FAIL: Bernstein target degree is too small"];
      Return[{}]];
    Table[Total[Table[Coefficient[transformed, tau, a]
          Binomial[i, a]/Binomial[degree, a],
        {a, 0, Min[i, actualDegree]}]], {i, 0, degree}]
  ];

PolynomialRemainderExact[dividend_, divisor_, variable_] :=
  Module[{remainder, divisorDegree, degree, leading},
    remainder = Expand[dividend];
    divisorDegree = Exponent[divisor, variable];
    leading = Coefficient[divisor, variable, divisorDegree];
    While[remainder =!= 0 && Exponent[remainder, variable] >= divisorDegree,
      degree = Exponent[remainder, variable];
      remainder = Expand[remainder -
        Coefficient[remainder, variable, degree]/leading
          variable^(degree - divisorDegree) divisor];
    ];
    remainder
  ];

SturmChainExact[polynomial_, variable_] := Module[{chain, remainder},
  chain = {Expand[polynomial], D[polynomial, variable]};
  While[True,
    remainder = -PolynomialRemainderExact[chain[[-2]], chain[[-1]], variable];
    remainder = Together[Expand[remainder]];
    If[TrueQ[remainder === 0], Break[]];
    AppendTo[chain, remainder];
  ];
  chain
];

ElementarySymmetricValues[xs_List] := Module[{es, count},
  count = Length[xs];
  es = Join[{1}, ConstantArray[0, count]];
  Do[Do[es[[k + 1]] = es[[k + 1]] + x es[[k]],
      {k, count, 1, -1}], {x, xs}];
  es
];

FormalDeterminantC[m_Integer, es_List] := Module[{e, h, prod, val},
  (* es[[k+1]] represents e_k. Values above degree m are set to zero,
     even when es was evaluated in more than m variables. *)
  e[k_Integer] := If[0 <= k <= m, es[[k + 1]], 0];
  h[i_Integer, j_Integer] := Total[Table[
    (i + j + 1 - 2 b) e[i + j + 1 - b] e[b], {b, 0, Min[i, j]}]];
  prod = Times @@ Table[e[j], {j, 1, m - 1}];
  val = 2 Total[Table[j (j + 1) e[j + 1]/e[j], {j, 1, m - 1}]] -
    2 Total[Table[h[j, j]/e[j]^2, {j, 1, m - 1}]] -
    Total[Flatten[Table[h[i, j]/(e[i] e[j]),
      {i, 1, m - 1}, {j, 1, m - 1}]]] +
    Total[Table[(e[1] e[j + 1] - (j + 2) e[j + 2])/e[j],
      {j, 1, m - 1}]]/e[1] - 2 e[2]/e[1];
  e[1] prod^2 val
];

DirectDeterminantC[xs_List] := Module[
  {m, es, prod, aa, vv, yy, omitted, ratios, p, y},
  m = Length[xs];
  es = ElementarySymmetricValues[xs];
  prod = Times @@ Table[es[[j + 1]], {j, 1, m - 1}];
  aa = 0; vv = 0; yy = 0;
  Do[
    y = xs[[i]];
    omitted = ElementarySymmetricValues[Delete[xs, i]];
    ratios = Table[omitted[[j + 1]]/es[[j + 1]], {j, 1, m - 1}];
    p = Total[ratios];
    aa = aa + y (Total[Table[j ratios[[j]], {j, 1, m - 1}]] -
      Total[ratios^2]);
    vv = vv + y p^2;
    yy = yy + y^2 p,
    {i, 1, m}];
  es[[2]] prod^2 (2 aa + yy/es[[2]] - vv - 2 es[[3]]/es[[2]])
];

DeterminantRemainder[r_Integer, es_List] := Module[{val, pref},
  val = r (r + 1);
  Do[val = val - (r + 1) (2 + If[p == 1, 1, 0] + If[2 p == r, 1, 0])
      es[[r + 1]]/(es[[p + 1]] es[[r - p + 1]]),
    {p, 1, Floor[r/2]}];
  Do[Do[val = val - (r - k + 1) (2 + If[k == r - 2 d, 1, 0])
        es[[k + 1]] es[[r + 1]]/(es[[k + d + 1]] es[[r - d + 1]]),
      {k, 1, r - 2 d}], {d, 1, Floor[(r - 1)/2]}];
  pref = es[[2]]^3 es[[r + 1]] es[[r + 2]]
    (Times @@ Table[es[[j + 1]]^2, {j, 2, r - 1}]);
  pref val
];

(* Isoparametric identities.
   Definitions only; run through the accompanying chapter harness.
   Exact calculations retain their original mathematical expressions. *)

CheckIsoparametricIdentities[] := Module[
  {z, g, d, t, ell, mu, nu, radialDefect, rows, gi, di, ti, li, c,
   expected, factor, coefficients, p, intervals, lo, hi, minimum,
   sturm, endpoint, signs, endpoints},
mu = (d + 1) ell - ell^2 - g d;
  nu = (d + 1) ell - g d;
  radialDefect = g^2 (1 - t^2) mu +
    g^2 t (ell^2 - 2 nu) z + (ell^2 - g^2) nu z^2;
  CheckEqual[IsoNumerator[g, d, t, ell, z, z], z radialDefect,
    "rank-one radial defect"];
  rows = {
    {2, 6, 1/3, 4, 3/10, (648/125) z^4 (11 z + 30)},
    {4, 14, 1/7, 8, 7/54,
      -(64/19683) z^4 (2548 z^2 - 132237 z - 910926)},
    {4, 14, 3/7, 8, 7/15,
      -(128/1125) z^4 (588 z^2 - 8015 z - 16550)}
  };

Do[
    {gi, di, ti, li, c, expected} = row;
    CheckEqual[z IsoNumerator[gi, di, ti, li, z (1 + c z), z],
      expected, "isoparametric correction g=" <> ToString[gi] <>
        ", t=" <> ToString[ti, InputForm]],
 {row, rows}];
  p = 189 z^5 + 1778 z^4 + 5492 z^3 + 4848 z^2 - 763 z + 42;
  CheckEqual[z IsoNumerator[4, 16, 3/4, 10, z (1 + z/2), z],
    z^2 p, "(4,1,7) correction"];
];

CheckIsoparametricExtended[] := Module[
  {g, t, lambda, z, d, ell, h, signPolynomial, cubicCoefficient,
   degree, u, v, vars, matrix, drift, radiusSquared, f, potential,
   gradientSquared, laplacian, expectedDefect, expectedCount = 7},
(* Critical dimension N=4g and endpoint potential degree ell=2g. *)
  d = 4 g - 2;
  ell = 2 g;
  h = z (1 + lambda z);
  signPolynomial = Expand[z IsoNumerator[g, d, t, ell, h, z]];
  Do[
    CheckEqual[Coefficient[signPolynomial, z, degree], 0,
      "critical isoparametric coefficient z^" <> ToString[degree] <>
        " vanishes"],
    {degree, 0, 2}];
  cubicCoefficient = 2 g^3 ((2 g + 1) (1 - t^2) lambda - 2 g t);
  CheckEqual[Coefficient[signPolynomial, z, 3], cubicCoefficient,
    "forced critical isoparametric cubic coefficient"];

  (* u=|x|^2 for x in R^2, v=|y|^2 for y in R^6. *)
  vars = {u, v};
  matrix = {{4 u, 0}, {0, 4 v}};
  drift = {4, 12};
  radiusSquared = u + v;
  f = v - 5 u;
  potential = radiusSquared f;
  gradientSquared = GammaExpr[potential, potential, vars, matrix];
  laplacian = LaplaceExpr[potential, vars, matrix, drift];
  CheckEqual[gradientSquared, 16 (25 u^3 + 24 u^2 v + v^3),
    "Lin one-sided gradient square"];
  CheckEqual[laplacian, 16 (v - 8 u),
    "Lin one-sided Laplacian"];
  expectedDefect = 64 radiusSquared f^2 (v - 17 u);
  CheckEqual[CurvatureExpr[potential, vars, matrix, drift], expectedDefect,
    "Lin one-sided defect 64 R F^2 (v-17u)"];

  expectedCount
];

CheckLawsonMotivatingPotential[] := Module[
  {u, v, f, p, q, xx, yy, vars, matrix, drift, changed,
   rad, split, quad, expected, expectedLap, i, j, bigU, bigV,
   defining, radius},
bigU = xx/p; bigV = yy/q;
  defining = bigU - bigV; radius = xx + yy;
  CheckEqual[(bigU + bigV)/4,
    (2 radius + (q - p) defining)/(4 (p + q)), "Lawson multiplier"];
  CheckEqual[(bigU^2 - bigV^2)/4,
    (bigU + bigV) defining/4, "Lawson quartic potential"];
  CheckEqual[p (radius + q defining)/(p + q), xx,
    "inverse squared x-radius"];
  CheckEqual[q (radius - p defining)/(p + q), yy,
    "inverse squared y-radius"];
];



Block[{$IterationLimit = Infinity, $Assumptions = True},
  RunVerificationGroup["CheckIsoparametricIdentities", CheckIsoparametricIdentities];
  RunVerificationGroup["CheckIsoparametricExtended", CheckIsoparametricExtended];
  RunVerificationGroup["CheckLawsonMotivatingPotential", CheckLawsonMotivatingPotential];
];

CertificateResult = <|"CertificateID" -> "LSC-ISOPARAMETRIC",
  "ExpectedChecks" -> 16,
  "Checks" -> (<|"Name" -> #[[1]], "Passed" -> TrueQ[#[[2]]],
      "Detail" -> If[StringQ[#[[3]]], #[[3]], ToString[#[[3]], InputForm]]|> &) /@ VerificationLog,
  "EvaluationOK" -> TrueQ[EvaluationOK],
  "EvaluationMessages" -> EvaluationMessages|>;
End[];
EndPackage[];
LaplacianCertificates`Isoparametric`Private`CertificateResult
