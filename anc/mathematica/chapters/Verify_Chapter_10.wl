(* ::Package:: *)
(* Complete code release 20260923-v6. Use its matching notebook. *)
(* Manuscript Chapter 12 (stable unit Chapter10): complex determinant cones -- merged manuscript verification.
   Target: Mathematica 12. Exact arithmetic; no external dependencies.
   This file runs independently and writes its actual execution results.
   Source mapping: Original general groups 6
*)
ClearAll["MergedVerificationV5`Chapter10`*"];
Begin["MergedVerificationV5`Chapter10`"];
VerificationTitle = "Manuscript Chapter 12 (stable unit Chapter10): complex determinant cones";
VerificationFileStem = "Verify_Chapter_10";
VerificationExpectedCount = 134;
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

(* Manuscript Chapter 12 (stable unit Chapter10): complex determinant cones.
   Definitions only; run through the accompanying chapter harness.
   Exact calculations retain their original mathematical expressions. *)

CheckDeterminants[] := Module[
  {u, v, d, r, denominator, b, rhs, expected, xs, es, scalarI, scalarJ},
d = u + 1;
  r = v + 2 d + 1;
  denominator = (r - 1) (d + 1) (d + 2) (r - d)^2;
  b = 2 (r - d + 1) (r + 2 d^2 + 3 d + 2)/((d + 1) (d + 2)) + 2 d + 1;
  rhs = 2 r^3/(r - 1) (1/(d + 1) -
    (r - 2 d) (r - 2 d - 1)/((d + 2) (r - d)^2));
  expected = 2 u v^4 + (4 u^3 + 25 u^2 + 29 u) v^3 +
    (16 u^4 + 102 u^3 + 206 u^2 + 140 u + 14) v^2 +
    (18 u^5 + 137 u^4 + 385 u^3 + 484 u^2 + 258 u + 44) v +
    4 u^6 + 38 u^5 + 140 u^4 + 252 u^3 + 230 u^2 + 104 u + 24;
  CheckEqual[(rhs - b) denominator, expected,
    "symbolic identity for the telescoping numerator"];
  Do[
    xs = Range[m];
    es = ElementarySymmetricValues[xs];
    CheckEqual[DirectDeterminantC[xs], FormalDeterminantC[m, es],
      "direct/formal determinant C, order " <> ToString[m]];
    CheckTrue[FormalDeterminantC[m, es] >= 0,
      "nonnegative determinant sample, order " <> ToString[m]];
    Do[
      es = ElementarySymmetricValues[Range[n]];
      If[m > 2,
        CheckEqual[FormalDeterminantC[m, es],
          es[[m]]^2 FormalDeterminantC[m - 1, es] +
            DeterminantRemainder[m - 1, es],
          "determinant recurrence, order " <> ToString[m] <>
            ", variables " <> ToString[n]]],
      {n, {m, m + 2}}],
    {m, 2, 10}];
  Do[
    scalarI = Total[Table[
      (rr + 1) (2 + If[p == 1, 1, 0] + If[2 p == rr, 1, 0])/
        Binomial[rr, p], {p, 1, Floor[rr/2]}]];
    scalarJ = Total[Flatten[Table[
      (rr - k + 1) (2 + If[k == rr - 2 dd, 1, 0])
        Binomial[k + dd, dd]/Binomial[rr, dd],
      {dd, 1, Floor[(rr - 1)/2]}, {k, 1, rr - 2 dd}]]];
    CheckTrue[scalarI + scalarJ <= rr (rr + 1),
      "exact scalar budget r=" <> ToString[rr]],
    {rr, 2, 100}];
];


RunChapterChecks[] := Block[{$IterationLimit = Infinity},
  Module[{started, failures},
    started = DateString[];
    VerificationLog = {}; VerificationGroupLog = {};
    VerificationMessages = {};
    Print[VerificationTitle];
    Print["Kernel: ", $Version];
    RunVerificationGroup["CheckDeterminants", CheckDeterminants];
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
MergedVerificationV5`Chapter10`VerificationResults
