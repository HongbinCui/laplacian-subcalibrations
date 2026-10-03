(* ::Package:: *)
(* Complete code release 20260923-v6. Use its matching notebook. *)
(* Manuscript Chapter 8 (stable unit Chapter07): gradient-square closed singular cubics -- merged manuscript verification.
   Target: Mathematica 12. Exact arithmetic; no external dependencies.
   This file runs independently and writes its actual execution results.
   Source mapping: Original general groups 3, 9
*)
ClearAll["MergedVerificationV5`Chapter07`*"];
Begin["MergedVerificationV5`Chapter07`"];
VerificationTitle = "Manuscript Chapter 8 (stable unit Chapter07): gradient-square closed singular cubics";
VerificationFileStem = "Verify_Chapter_07";
VerificationExpectedCount = 18;
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

(* Manuscript Chapter 8 (stable unit Chapter07): gradient-square closed singular cubics.
   Definitions only; run through the accompanying chapter harness.
   Exact calculations retain their original mathematical expressions. *)

CheckCubic[] := Module[
  {r, s, f, k, n, vars, matrix, drift, potential, norm, defect},
vars = {r, f, s};
  matrix = {{4 r, 6 f, 8 s}, {6 f, s, 4 r f},
            {8 s, 4 r f, 4 r s + 12 f^2}};
  n = 6 k - 3;
  drift = {2 n, 0, 4 k r};
  potential = s f;
  norm = s^3 + 12 r s f^2 + 12 f^4;
  defect = (4 k - 16) r s^3 +
    f^2 s (48 k r^2 - 126 s) + 12 (4 k - 6) r f^4;
  CheckEqual[GammaExpr[potential, potential, vars, matrix], norm,
    "cubic SF gradient"];
  CheckEqual[LaplaceExpr[potential, vars, matrix, drift],
    (4 k + 8) r f, "cubic SF Laplacian"];
  CheckEqual[CurvatureExpr[potential, vars, matrix, drift], f defect,
    "cubic SF defect"];
  CheckEqual[defect - (defect /. k -> 4), 4 (k - 4) r norm,
    "cubic parameter propagation"];
];

CheckCubicSpectrum[] := Module[
  {f, s, e, n, alpha, k, p, x, r, q, y, z, vars, matrix,
   momentResidual, phi, jphi, expected,
   discriminantRoot, smallRoot, eigenvalue, potential, norm,
   quarticLaplacian, eigenvalues, powerTwo, powerThree, powerFour,
   modelF, modelS, gradientF, gradientS, projection, diagonalHessian,
   offDiagonalTrace, offDiagonalSquare, angularJacobi, radialJacobi,
   certificateBound, i, j, a, aSquared, maxS, t, b, c, testF, testS,
   hessS, residual, cubicMatrix, cubicDrift},
(* Angular eigenfunction, cutoff exponent, sharp Hardy value: 5 checks. *)
  phi = x^(-p);
  jphi = (4 x - 16 x^2) D[phi, {x, 2}] +
    (4 k - 2 - (24 k - 8) x) D[phi, x] + 2 (k - 1) phi/x;
  expected = (4 p^2 + (6 - 4 k) p + 2 k - 2)/x +
    4 p (6 k - 6 - 4 p);
  CheckEqual[jphi/phi, expected, "cubic angular power Jacobi identity"];
  discriminantRoot = Sqrt[4 k^2 - 20 k + 17];
  smallRoot = (2 k - 3 - discriminantRoot)/4;
  eigenvalue = 4 smallRoot (6 k - 6 - 4 smallRoot);
  CheckEqual[4 smallRoot^2 + (6 - 4 k) smallRoot + 2 k - 2, 0,
    "cubic smaller eigenfunction root"];
  CheckEqual[2 k - 3 - 4 smallRoot, discriminantRoot,
    "cubic singular cutoff exponent"];
  CheckEqual[9 (k - 1)^2 - eigenvalue, (k + discriminantRoot)^2,
    "cubic sharp Hardy square"];
  CheckEqual[(k + discriminantRoot)^2 /. k -> 5, 42 + 10 Sqrt[17],
    "cubic sharp Hardy value in dimension 27"];
];

CheckQuaternionicCubicModel[] := Module[
  {f, s, e, n, alpha, k, p, x, r, q, y, z, vars, matrix,
   momentResidual, phi, jphi, expected,
   discriminantRoot, smallRoot, eigenvalue, potential, norm,
   quarticLaplacian, eigenvalues, powerTwo, powerThree, powerFour,
   modelF, modelS, gradientF, gradientS, projection, diagonalHessian,
   offDiagonalTrace, offDiagonalSquare, angularJacobi, radialJacobi,
   certificateBound, i, j, a, aSquared, maxS, t, b, c, testF, testS,
   hessS, residual, cubicMatrix, cubicDrift},
(* Quaternionic Hermitian tracefree 4 by 4 model: 6 checks.
     Conjugation invariance reduces the scalar identities to a real
     diagonal matrix. Each off-diagonal block has multiplicity four. *)
  eigenvalues = {x, y, z, -x - y - z};
  powerTwo = Total[Table[eigenvalues[[i]]^2, {i, 4}]];
  powerThree = Total[Table[eigenvalues[[i]]^3, {i, 4}]];
  powerFour = Total[Table[eigenvalues[[i]]^4, {i, 4}]];
  modelF = powerThree/3;
  modelS = powerFour - powerTwo^2/4;
  gradientF = Table[eigenvalues[[i]]^2 - powerTwo/4, {i, 4}];
  gradientS = Table[4 eigenvalues[[i]]^3 - powerThree -
    powerTwo eigenvalues[[i]], {i, 4}];
  CheckEqual[gradientF.gradientF, modelS,
    "quaternionic cubic normalized gradient square"];
  CheckEqual[gradientF.gradientS, 4 powerTwo modelF,
    "quaternionic cubic mixed gradient closure"];
  CheckEqual[gradientS.gradientS, 4 powerTwo modelS + 12 modelF^2,
    "quaternionic cubic quartic gradient closure"];
  CheckTrue[(gradientF /. {x -> 1, y -> 1, z -> -1}) ===
    ConstantArray[0, 4], "quaternionic cubic nonzero critical point"];
  projection = IdentityMatrix[4] - ConstantArray[1, {4, 4}]/4;
  diagonalHessian = projection.DiagonalMatrix[2 eigenvalues].projection;
  offDiagonalTrace = 4 Total[Flatten[Table[
    eigenvalues[[i]] + eigenvalues[[j]], {i, 4}, {j, i + 1, 4}]]];
  offDiagonalSquare = 4 Total[Flatten[Table[
    (eigenvalues[[i]] + eigenvalues[[j]])^2,
    {i, 4}, {j, i + 1, 4}]]];
  CheckEqual[Tr[diagonalHessian] + offDiagonalTrace, 0,
    "quaternionic cubic harmonic Hessian trace"];
  CheckEqual[Tr[diagonalHessian.diagonalHessian] + offDiagonalSquare,
    10 powerTwo, "quaternionic cubic Hessian norm and Bochner coefficient"];
];

CheckCubicCertificateHardyBound[] := Module[
  {f, s, e, n, alpha, k, p, x, r, q, y, z, vars, matrix,
   momentResidual, phi, jphi, expected,
   discriminantRoot, smallRoot, eigenvalue, potential, norm,
   quarticLaplacian, eigenvalues, powerTwo, powerThree, powerFour,
   modelF, modelS, gradientF, gradientS, projection, diagonalHessian,
   offDiagonalTrace, offDiagonalSquare, angularJacobi, radialJacobi,
   certificateBound, i, j, a, aSquared, maxS, t, b, c, testF, testS,
   hessS, residual, cubicMatrix, cubicDrift},
(* The SF certificate bound is distinct from the sharp value: 3 checks. *)
  p = 3/2;
  angularJacobi = (4 p^2 + (6 - 4 k) p + 2 k - 2)/s +
    4 p (6 k - 6 - 4 p);
  radialJacobi = (-6) (-6 + (6 k - 4) - 2);
  CheckEqual[-(angularJacobi + radialJacobi), 4 (k - 4)/s,
    "cubic SF boundary Jacobi coefficient"];
  certificateBound = 9 (k - 3)^2 + 16 (k - 4);
  CheckEqual[(7 - (6 k - 4)/2)^2 + 16 (k - 4), certificateBound,
    "cubic SF certificate Hardy bound"];
  CheckEqual[certificateBound /. k -> 5, 52,
    "cubic SF certificate value in dimension 27"];
];


RunChapterChecks[] := Block[{$IterationLimit = Infinity},
  Module[{started, failures},
    started = DateString[];
    VerificationLog = {}; VerificationGroupLog = {};
    VerificationMessages = {};
    Print[VerificationTitle];
    Print["Kernel: ", $Version];
    RunVerificationGroup["CheckCubic", CheckCubic];
    RunVerificationGroup["CheckCubicSpectrum", CheckCubicSpectrum];
    RunVerificationGroup["CheckQuaternionicCubicModel", CheckQuaternionicCubicModel];
    RunVerificationGroup["CheckCubicCertificateHardyBound", CheckCubicCertificateHardyBound];
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
MergedVerificationV5`Chapter07`VerificationResults
