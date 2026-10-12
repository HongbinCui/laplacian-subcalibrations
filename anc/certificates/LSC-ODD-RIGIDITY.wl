(* ::Package:: *)
(* LSC-ODD-RIGIDITY -- Odd Rigidity.
   Compatible with Mathematica 12. No external package or file dependencies.
   Get returns checks; only the central runner writes execution records.

*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`OddRigidity`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`OddRigidity`Private`*"];

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

(* Odd Rigidity identities.
   Definitions only; run through the accompanying chapter harness.
   Exact calculations retain their original mathematical expressions. *)

CheckOddDegreeRigidityAndMoments[] := Module[
  {f, s, e, n, alpha, k, p, x, r, q, y, z, vars, matrix,
   momentResidual, phi, jphi, expected,
   discriminantRoot, smallRoot, eigenvalue, potential, norm,
   quarticLaplacian, eigenvalues, powerTwo, powerThree, powerFour,
   modelF, modelS, gradientF, gradientS, projection, diagonalHessian,
   offDiagonalTrace, offDiagonalSquare, angularJacobi, radialJacobi,
   certificateBound, i, j, a, aSquared, maxS, t, b, c, testF, testS,
   hessS, residual, cubicMatrix, cubicDrift},
(* Second-order rigidity and cubic spherical moments: 15 checks. *)
  aSquared = 2/(e^3 (e - 1));
  maxS = 2/(e (e - 1));
  CheckEqual[e^2 (4 (e - 1)^2 maxS - 4), 4 e (e - 2),
    "odd-degree coefficient normalization at a maximum"];
  matrix = {{s - e^2 f^2, f (4 - 2 e (e - 1) s)},
    {f (4 - 2 e (e - 1) s),
     4 s + 4 e (e - 2) f^2 - 4 (e - 1)^2 s^2}};
  CheckTrue[And @@ (TrueQ[AlgebraicResidual[#] === 0] & /@
    Flatten[Expand[matrix /. {f -> a, s -> maxS}] /. a^2 -> aSquared]),
    "spherical gradient matrix at a maximum"];

  testF = a + b t^2/2;
  testS = maxS + c t^2/2;
  residual = Expand[(b t)^2 - (testS - e^2 testF^2)];
  CheckEqual[Coefficient[residual, t, 2],
    b^2 + e^2 a b - c/2, "first spherical Hessian equation"];
  residual = Expand[b c t^2 - testF (4 - 2 e (e - 1) testS)];
  CheckEqual[Coefficient[residual, t, 2],
    b c + e (e - 1) a c, "mixed spherical Hessian equation"];
  residual = Expand[c^2 t^2 -
    (4 testS + 4 e (e - 2) testF^2 - 4 (e - 1)^2 testS^2)];
  CheckEqual[Coefficient[residual, t, 2],
    c^2 - (8/e - 6) c - 4 e (e - 2) a b,
    "last spherical Hessian equation"];
  hessS = 2 b (b + e^2 a);
  CheckEqual[(b + e (e - 1) a) hessS,
    2 b (b + e^2 a) (b + e (e - 1) a),
    "factorization of the Hessian eigenvalue restriction"];
  CheckEqual[Expand[(c^2 - (8/e - 6) c - 4 e (e - 2) a b) /.
    {b -> -e^2 a, c -> 0}] /. a^2 -> aSquared,
    8 (e - 2)/(e - 1), "excluded Hessian eigenvalue residual"];
  CheckEqual[Expand[hessS /. b -> -e (e - 1) a] /. a^2 -> aSquared,
    -4/e, "remaining Hessian eigenvalue of the gradient square"];
  CheckEqual[Expand[(c^2 - (8/e - 6) c - 4 e (e - 2) a b) /.
    {b -> -e (e - 1) a, c -> -4/e}] /. a^2 -> aSquared,
    8 (e - 2) (e - 3)/e^2, "second-order degree rigidity obstruction"];

  vars = {f, s};
  cubicMatrix = matrix /. e -> 3;
  cubicDrift = {-3 (n + 1) f, alpha - 4 (n + 2) s};
  CheckEqual[LaplaceExpr[f^2, vars, cubicMatrix, cubicDrift],
    2 (s - 3 (n + 4) f^2),
    "first cubic spherical moment differential identity"];
  CheckEqual[GammaExpr[f^2, s, vars, cubicMatrix] +
    f^2 LaplaceExpr[s, vars, cubicMatrix, cubicDrift],
    (alpha + 8) f^2 - 4 (n + 8) f^2 s,
    "mixed cubic spherical moment differential identity"];
  CheckEqual[GammaExpr[f, f s, vars, cubicMatrix] +
    f s LaplaceExpr[f, vars, cubicMatrix, cubicDrift],
    s^2 + 4 f^2 - 3 (n + 8) f^2 s,
    "third cubic spherical moment differential identity"];
  CheckEqual[LaplaceExpr[s^2, vars, cubicMatrix, cubicDrift],
    2 ((alpha + 4) s + 12 f^2 - 4 (n + 6) s^2),
    "fourth cubic spherical moment differential identity"];
  momentResidual = (n + 6) (3 alpha + 8) -
    3 (n + 4) (alpha + 4) - 12;
  CheckEqual[momentResidual, 6 alpha - 4 n - 12,
    "cubic spherical moment elimination"];
  CheckEqual[momentResidual /. alpha -> 2 (n + 3)/3, 0,
    "normalized cubic Bochner coefficient"];
];



Block[{$IterationLimit = Infinity, $Assumptions = True},
  RunVerificationGroup["CheckOddDegreeRigidityAndMoments", CheckOddDegreeRigidityAndMoments];
];

CertificateResult = <|"CertificateID" -> "LSC-ODD-RIGIDITY",
  "ExpectedChecks" -> 15,
  "Checks" -> (<|"Name" -> #[[1]], "Passed" -> TrueQ[#[[2]]],
      "Detail" -> If[StringQ[#[[3]]], #[[3]], ToString[#[[3]], InputForm]]|> &) /@ VerificationLog,
  "EvaluationOK" -> TrueQ[EvaluationOK],
  "EvaluationMessages" -> EvaluationMessages|>;
End[];
EndPackage[];
LaplacianCertificates`OddRigidity`Private`CertificateResult
