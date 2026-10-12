(* ::Package:: *)
(* LSC-POWER-MULTIPLIERS -- Power Multipliers.
   Compatible with Mathematica 12. No external package or file dependencies.
   Get returns checks; only the central runner writes execution records.

   Symbolic real powers use the original positive-base convention for radial and gradient-square variables.
*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`PowerMultipliers`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`PowerMultipliers`Private`*"];

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

(* Power Multipliers identities.
   Definitions only; run through the accompanying chapter harness.
   Exact calculations retain their original mathematical expressions. *)

CheckGradientPower[] := Module[
  {f, s, u, v, uf, us, vf, vs, delta, b, t,
   vars, matrix, h, norm, lap, lambda, defect, direct},
(* The unused bottom-right block may be set to zero: the first argument
     of every contraction here depends only on f and s. *)
  vars = {f, s, u, v};
  matrix = {{s, u, uf, vf}, {u, v, us, vs},
            {uf, us, 0, 0}, {vf, vs, 0, 0}};
  h = s^t f;
  norm = s + 2 t f u/s + t^2 f^2 v/s^2;
  lap = delta + 2 t u/s + t f b/s + t (t - 1) f v/s^2;
  lambda = (s delta - u/2)/f;
  defect = lambda + t b +
    (2 t u delta - t uf - (t^2 + 3 t/2) v)/s +
    t (t + 1) u^2/s^2 +
    t^2 f/s^2 (2 u b + v delta - us - vf/2) +
    t^3 f^2/s^3 (v b - vs/2);
  CheckEqual[GammaExpr[h, h, vars, matrix], s^(2 t) norm,
    "gradient-power gradient"];
  direct = s^(3 t) norm lap -
    GammaExpr[h, s^(2 t) norm, vars, matrix]/2;
  CheckEqual[direct/s^(3 t), f defect,
    "gradient-power defect for symbolic real exponent"];
];

CheckRadialMaster[] := Module[
  {r, h, s, u, delta, m, a, n, vars, matrix, drift,
   c, potential, expected},
vars = {r, h, s};
  matrix = {{4 r, 2 m h, 4 (m - 1) s},
            {2 m h, s, u}, {4 (m - 1) s, u, 0}};
  drift = {2 n, delta, 0};
  c = 4 a (m + a);
  potential = r^a h;
  expected = s delta - u/2 + c h^2 delta/r +
    (2 a (n - 1) - c) h s/r +
    c (m + 2 a (n - 1)) h^3/r^2;
  CheckEqual[GammaExpr[potential, potential, vars, matrix]/r^(2 a),
    s + c h^2/r, "radial gradient"];
  CheckEqual[CurvatureExpr[potential, vars, matrix, drift]/r^(3 a),
    expected, "radial curvature"];
  CheckEqual[LaplaceExpr[potential, vars, matrix, drift]/r^a,
    delta + 2 a (2 m + n + 2 a - 2) h/r, "radial Laplacian"];
];

CheckCombinedPowers[] := Module[
  {f, a, p, q, p1, p2, p3, q1, q2, q3,
   f11, f12, f13, f22, f23, f33,
   a11, a12, a13, a22, a23, a33,
   hf, ha, trf, tra, ssJet, ttJet, vvJet,
   lf, la, gg, hg, direct, expected, d1, d2, composedGrad,
   composedHess, i, j,
   radius, z, ell, ambientDim, angular, metric, drift, vars,
   angularMetric, angularDrift, potential, norm, w, numerator,
   rr, ff, ss, uu, vv, uf, us, vf, vs, delta, bb,
   exponent, alpha, gamma, m, cc, bc, vgamma, dgamma, qgamma,
   vcombined, qcombined, coneDim, newEll},
(* Combined powers are checked directly, not only by substituting the
     separate radial and gradient-power identities. The unused lower
     block of the formal gradient matrix is zero: the first argument
     in every required contraction depends only on rr, ff, and ss. *)
  vars = {rr, ff, ss, uu, vv};
  metric = {
    {4 rr, 2 exponent ff, 4 (exponent - 1) ss,
      2 (3 exponent - 4) uu, 2 (4 exponent - 6) vv},
    {2 exponent ff, ss, uu, uf, vf},
    {4 (exponent - 1) ss, uu, vv, us, vs},
    {2 (3 exponent - 4) uu, uf, us, 0, 0},
    {2 (4 exponent - 6) vv, vf, vs, 0, 0}};
  drift = {2 ambientDim, delta, bb, 0, 0};
  potential = rr^alpha ss^gamma ff;
  norm = GammaExpr[potential, potential, vars, metric];
  m = exponent + 2 gamma (exponent - 1);
  cc = 4 alpha (m + alpha);
  bc = 2 alpha (ambientDim - 1) - cc;
  vgamma = ss + 2 gamma ff uu/ss + gamma^2 ff^2 vv/ss^2;
  dgamma = delta + 2 gamma uu/ss + gamma ff bb/ss +
    gamma (gamma - 1) ff vv/ss^2;
  qgamma = (ss delta - uu/2)/ff + gamma bb +
    (2 gamma uu delta - gamma uf - (gamma^2 + 3 gamma/2) vv)/ss +
    gamma (gamma + 1) uu^2/ss^2 +
    gamma^2 ff/ss^2 (2 uu bb + vv delta - us - vf/2) +
    gamma^3 ff^2/ss^3 (vv bb - vs/2);
  vcombined = vgamma + cc ff^2/rr;
  qcombined = qgamma + cc ff dgamma/rr + bc vgamma/rr +
    cc (m + 2 alpha (ambientDim - 1)) ff^2/rr^2;
  CheckEqual[norm/(rr^(2 alpha) ss^(2 gamma)), vcombined,
    "combined powers: direct gradient identity"];
  direct = norm LaplaceExpr[potential, vars, metric, drift] -
    GammaExpr[potential, norm, vars, metric]/2;
  CheckEqual[direct/(rr^(3 alpha) ss^(3 gamma)), ff qcombined,
    "combined powers: direct curvature identity"];
];

CheckRadialBoundaryPerturbation[] := Module[
  {f, a, p, q, p1, p2, p3, q1, q2, q3,
   f11, f12, f13, f22, f23, f33,
   a11, a12, a13, a22, a23, a33,
   hf, ha, trf, tra, ssJet, ttJet, vvJet,
   lf, la, gg, hg, direct, expected, d1, d2, composedGrad,
   composedHess, i, j,
   radius, z, ell, ambientDim, angular, metric, drift, vars,
   angularMetric, angularDrift, potential, norm, w, numerator,
   rr, ff, ss, uu, vv, uf, us, vf, vs, delta, bb,
   exponent, alpha, gamma, m, cc, bc, vgamma, dgamma, qgamma,
   vcombined, qcombined, coneDim, newEll},
CheckEqual[(1 - ell) (coneDim - 1 - ell) -
      (1 - newEll) (coneDim - 1 - newEll),
    (ell - coneDim/2)^2 - (newEll - coneDim/2)^2,
    "radial perturbation: exact boundary coefficient change"];
];



Block[{$IterationLimit = Infinity, $Assumptions = True},
  RunVerificationGroup["CheckGradientPower", CheckGradientPower];
  RunVerificationGroup["CheckRadialMaster", CheckRadialMaster];
  RunVerificationGroup["CheckCombinedPowers", CheckCombinedPowers];
  RunVerificationGroup["CheckRadialBoundaryPerturbation", CheckRadialBoundaryPerturbation];
];

CertificateResult = <|"CertificateID" -> "LSC-POWER-MULTIPLIERS",
  "ExpectedChecks" -> 8,
  "Checks" -> (<|"Name" -> #[[1]], "Passed" -> TrueQ[#[[2]]],
      "Detail" -> If[StringQ[#[[3]]], #[[3]], ToString[#[[3]], InputForm]]|> &) /@ VerificationLog,
  "EvaluationOK" -> TrueQ[EvaluationOK],
  "EvaluationMessages" -> EvaluationMessages|>;
End[];
EndPackage[];
LaplacianCertificates`PowerMultipliers`Private`CertificateResult
