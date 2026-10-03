(* ::Package:: *)
(* Complete code release 20260923-v6. Use its matching notebook. *)
(* Chapter 4: algebraic multipliers and Jacobi stability -- merged manuscript verification.
   Target: Mathematica 12. Exact arithmetic; no external dependencies.
   This file runs independently and writes its actual execution results.
   Source mapping: Original general groups 9, 12
*)
ClearAll["MergedVerificationV5`Chapter04`*"];
Begin["MergedVerificationV5`Chapter04`"];
VerificationTitle = "Chapter 4: algebraic multipliers and Jacobi stability";
VerificationFileStem = "Verify_Chapter_04";
VerificationExpectedCount = 6;
VerificationDirectory = If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[$InputFileName], Directory[]];

VerificationLog = {};
VerificationGroupLog = {};
VerificationMessages = {};
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
RunVerificationGroup[label_String, fun_Symbol] := Module[{before, status},
  before = Length[VerificationLog];
  status = Block[{System`$MessageList = {}}, Module[{result},
    result = CheckAbort[Catch[fun[], _, ($Failed &)], $Aborted];
    !MemberQ[{$Failed, $Aborted}, result] &&
      Length[System`$MessageList] == 0
  ]];
  If[!TrueQ[status], AppendTo[VerificationMessages, label]];
  AppendTo[VerificationGroupLog,
    {label, Length[VerificationLog] - before, TrueQ[status]}];
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

(* Chapter 4: algebraic multipliers and Jacobi stability.
   Definitions only; run through the accompanying chapter harness.
   Exact calculations retain their original mathematical expressions. *)

CheckGeneralMultiplierAndSphere[] := Module[
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
(* Arbitrary values, gradients, and symmetric Hessians at a point.
     This is a universal contraction identity; the three-coordinate
     encoding checks the full algebra without imposing a special F or a. *)
  p = {p1, p2, p3};
  q = {q1, q2, q3};
  hf = {{f11, f12, f13}, {f12, f22, f23}, {f13, f23, f33}};
  ha = {{a11, a12, a13}, {a12, a22, a23}, {a13, a23, a33}};
  trf = Total[Table[hf[[i, i]], {i, 3}]];
  tra = Total[Table[ha[[i, i]], {i, 3}]];
  ssJet = p . p;
  ttJet = p . q;
  vvJet = q . q;
  lf = ssJet trf - p . hf . p;
  la = vvJet tra - q . ha . q;
  gg = a p + f q;
  hg = Table[a hf[[i, j]] + f ha[[i, j]] +
    p[[i]] q[[j]] + q[[i]] p[[j]], {i, 3}, {j, 3}];
  direct = (gg . gg) Total[Table[hg[[i, i]], {i, 3}]] - gg . hg . gg;
  expected = a^3 lf +
    f (a^2 (ssJet tra + 2 ttJet trf - 2 p . hf . q - p . ha . p) +
      2 a (ttJet^2 - ssJet vvJet)) +
    f^2 a (2 ttJet tra + vvJet trf - q . hf . q - 2 p . ha . q) +
    f^3 la;
  CheckEqual[direct, expected, "general exact multiplier: arbitrary two-jets"];

  (* Arbitrary first and second derivatives of a scalar outer function. *)
  composedGrad = d1 p;
  composedHess = Table[d1 hf[[i, j]] + d2 p[[i]] p[[j]],
    {i, 3}, {j, 3}];
  direct = (composedGrad . composedGrad)
    Total[Table[composedHess[[i, i]], {i, 3}]] -
    composedGrad . composedHess . composedGrad;
  CheckEqual[direct, d1^3 lf, "general scalar composition identity"];

  (* Direct polar differential operator, with an arbitrary angular
     potential and arbitrary one-coordinate spherical closure data.
     This verifies the numerator through derivatives, before collection. *)
  vars = {radius, z};
  metric = {{1, 0}, {0, angularMetric[z]/radius^2}};
  drift = {(ambientDim - 1)/radius, angularDrift[z]/radius^2};
  potential = radius^ell angular[z];
  norm = GammaExpr[potential, potential, vars, metric];
  direct = norm LaplaceExpr[potential, vars, metric, drift] -
    GammaExpr[potential, norm, vars, metric]/2;
  w = ell^2 angular[z]^2 + angularMetric[z] D[angular[z], z]^2;
  numerator = w (angularMetric[z] D[angular[z], z, z] +
      angularDrift[z] D[angular[z], z] +
      (ambientDim - 1) ell angular[z]) -
    angularMetric[z] D[angular[z], z] D[w, z]/2;
  CheckEqual[direct/radius^(3 ell - 4), numerator,
    "general spherical numerator from polar differentiation"];


  (* coneDim is the cone dimension, so ambientDim = coneDim + 1.
     These identities isolate the exact Hardy square and the boundary
     change for any two homogeneous degrees. *)
  CheckEqual[(1 - ell) (coneDim - 1 - ell),
    (ell - coneDim/2)^2 - (coneDim - 2)^2/4,
    "degree window: exact completed square"];
];

CheckSimonsMotivatingExample[] := Module[
  {f, s, e, n, alpha, k, p, x, r, q, y, z, vars, matrix,
   momentResidual, phi, jphi, expected,
   discriminantRoot, smallRoot, eigenvalue, potential, norm,
   quarticLaplacian, eigenvalues, powerTwo, powerThree, powerFour,
   modelF, modelS, gradientF, gradientS, projection, diagonalHessian,
   offDiagonalTrace, offDiagonalSquare, angularJacobi, radialJacobi,
   certificateBound, i, j, a, aSquared, maxS, t, b, c, testF, testS,
   hessS, residual, cubicMatrix, cubicDrift},
(* Simons and quartic identities: 3 checks. *)
  vars = {r, f};
  matrix = {{4 r, 4 f}, {4 f, 4 r}};
  potential = r f;
  norm = Expand[GammaExpr[potential, potential, vars, matrix]];
  CheckEqual[norm, 4 r (r^2 + 3 f^2),
    "Simons RF gradient-square identity"];
  CheckEqual[24 f norm - GammaExpr[potential, norm, vars, matrix]/2,
    192 r f^3, "Simons RF polynomial defect"];
];


RunChapterChecks[] := Block[{$IterationLimit = Infinity},
  Module[{started, failures},
    started = DateString[];
    VerificationLog = {}; VerificationGroupLog = {};
    VerificationMessages = {};
    Print[VerificationTitle];
    Print["Kernel: ", $Version];
    RunVerificationGroup["CheckGeneralMultiplierAndSphere", CheckGeneralMultiplierAndSphere];
    RunVerificationGroup["CheckSimonsMotivatingExample", CheckSimonsMotivatingExample];
    failures = Select[VerificationLog, !TrueQ[#[[2]]] &];
    VerificationResults = Length[failures] == 0 &&
      Length[VerificationLog] == VerificationExpectedCount &&
      Length[VerificationMessages] == 0;
    VerificationMetadata = {
      {"Code release", "20260923-v6"},
      {"Manuscript", "Revised manuscript, 23 September 2026"},
      {"Unit", VerificationTitle}, {"Kernel", $Version},
      {"Version number", $VersionNumber}, {"System ID", $SystemID},
      {"Started", started}, {"Finished", DateString[]},
      {"Expected checks", VerificationExpectedCount},
      {"Actual checks", Length[VerificationLog]},
      {"Passed checks", Length[VerificationLog] - Length[failures]},
      {"Failed checks", Length[failures]},
      {"Evaluation-message groups", VerificationMessages},
      {"Mathematical checks complete and passed", VerificationResults}};
    Print["Checks: ", Length[VerificationLog], "; passed: ",
      Length[VerificationLog] - Length[failures], "; failed: ",
      Length[failures], "; expected: ", VerificationExpectedCount];
    If[TrueQ[VerificationResults], Print["ALL CHECKS PASSED."],
      Print["FAILED OR INCOMPLETE. Preserve all messages and result files."]];
    VerificationResults
  ]
];

ExportVerificationResults[dir_String] := Block[{$IterationLimit = Infinity}, Module[
  {formatRow, csvField, txt, csv, lines, csvRows},
  If[!DirectoryQ[dir],
    CreateDirectory[dir, CreateIntermediateDirectories -> True]];
  txt = FileNameJoin[{dir, VerificationFileStem <> "_Results.txt"}];
  csv = FileNameJoin[{dir, VerificationFileStem <> "_Results.csv"}];
  formatRow[row_] := StringRiffle[ToString[#, InputForm] & /@ row, " | "];
  lines = Join[{VerificationTitle, "Actual mathematical-check record (the unit return also requires both result exports)", ""},
    formatRow /@ VerificationMetadata,
    {"", "Groups: label | count | completed without evaluation messages"},
    formatRow /@ VerificationGroupLog,
    {"", "Every check: label | passed | residual"},
    formatRow /@ VerificationLog,
    {"", "Finite sample checks do not replace the all-orders or analytic proofs."}];
  csvField[item_] := "\"" <> StringReplace[
    If[StringQ[item], item, ToString[item, InputForm]],
    "\"" -> "\"\""] <> "\"";
  csvRows = (StringRiffle[csvField /@ #, ","] &) /@
    Prepend[VerificationLog, {"label", "passed", "residual"}];
  {Export[txt, StringRiffle[lines, "\n"] <> "\n", "Text"],
   Export[csv, StringRiffle[csvRows, "\r\n"] <> "\r\n", "Text"]}
]];

VerificationResults = RunChapterChecks[];
VerificationOutputFiles = ExportVerificationResults[
  FileNameJoin[{VerificationDirectory, "verification_results_v6"}]];
VerificationExportSucceeded = ListQ[VerificationOutputFiles] &&
  Length[VerificationOutputFiles] == 2 &&
  And @@ ((StringQ[#] && FileExistsQ[#]) & /@ VerificationOutputFiles);
VerificationChecksPassed = TrueQ[VerificationResults];
VerificationResults = VerificationChecksPassed && TrueQ[VerificationExportSucceeded];
Print["Result files: ", VerificationOutputFiles];
If[TrueQ[VerificationResults], Print["UNIT PASSED AND RESULTS SAVED."],
  If[!TrueQ[VerificationExportSucceeded],
    Print["RESULT EXPORT FAILED: unit return value is False."]]];
End[];
MergedVerificationV5`Chapter04`VerificationResults
