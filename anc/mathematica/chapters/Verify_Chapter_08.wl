(* ::Package:: *)
(* Base code release 20260923-v6; notation revision 20260928-v3.
   Use its matching notebook. Assertions and exact formulas are unchanged. *)
(* Manuscript Chapter 9 (stable unit Chapter08): Area-minimizing Clifford cubic cones -- merged manuscript verification.
   Target: Mathematica 12. Exact arithmetic; no external dependencies.
   This file runs independently and writes its actual execution results.
   Source mapping: Original Clifford assertions
*)
ClearAll["MergedVerificationV5`Chapter08`*"];
Begin["MergedVerificationV5`Chapter08`"];
VerificationTitle = "Manuscript Chapter 9 (stable unit Chapter08): Area-minimizing Clifford cubic cones";
VerificationFileStem = "Verify_Chapter_08";
VerificationExpectedCount = 193;
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


(* Manuscript Chapter 9 (stable unit Chapter08): all 193 assertions from the former Paper II program.
   CheckEqual and CheckTrue are supplied by the independent-file wrapper.
   Matrix checks and all coordinate-loop assertions are retained unchanged.
   The historical Mathematica 12 run belongs to the original program only. *)
CheckCliffordClassification[] := Module[{a, actual, al, b, b0, b1, bernstein, bg1, block, c, cc, condition, cubic, d, dd, deg, dhat, dim, direction, ds, edge, ee, ei, ej, expected, expectedJacobi, exponents, expr, f, f2, fval, ga, gf, gfv, gradient, gs, gsv, hessian, hf, hfv, hh, hhat, hopf, hs, hssUU, hssUZ, hssZZ, hsv, i, id4, integrand, j, jacobi, jacobiRec, k, jLink, jGradS, jValue, l, label, lam, m, m0, m1, matrices, matrixCheck, mm, moment, momentGamma, momentRecurrence, must, n0, n083, norm2, normalized, nx, ny, omega, p, p43, pg1, pg2, poly, positive, q, qhat, qq, qv, rebuilt, report, residuals, rr, rules, rv, sj, sphereGrad, spoly, ss, sv, t, tau, theta, tm, tv, u, um, uv, v, vars, vhat, vv, vv0, weighted, x, x1, x2, x3, x4, x5, x6, x7, x8, xs, y, yExamples, ys, yval, z, zExamples, zero, zs, zval, zz4},
 report[parts___] := Null;
 Catch[Block[{$IterationLimit = 1000000},
(* Block 1. Exact checks and conversion to the Bernstein basis. *)
must[condition_, label_] := If[!TrueQ[CheckTrue[condition, label]],
  Throw[$Failed, "CliffordClassificationFailure"]];
zero[expr_, label_] := If[!TrueQ[CheckEqual[expr, 0, label]],
  Throw[$Failed, "CliffordClassificationFailure"]];
bernstein[poly_, x_, y_, nx_, ny_] := Module[{p, a, b, rebuilt},
  p = Expand[Cancel[poly]];
  must[PolynomialQ[p, {x, y}] && Exponent[p, x] <= nx &&
    Exponent[p, y] <= ny, "Bernstein bidegree"];
  a = Table[Coefficient[Coefficient[p, x, k], y, l],
    {k, 0, nx}, {l, 0, ny}];
  b = Table[Total[Flatten[Table[
    a[[k + 1, l + 1]] Binomial[i, k]/Binomial[nx, k]
    Binomial[j, l]/Binomial[ny, l], {k, 0, i}, {l, 0, j}]]],
    {i, 0, nx}, {j, 0, ny}];
  rebuilt = Total[Flatten[Table[
    b[[i + 1, j + 1]] Binomial[nx, i] x^i (1 - x)^(nx - i)
    Binomial[ny, j] y^j (1 - y)^(ny - j), {i, 0, nx}, {j, 0, ny}]]];
  zero[p - rebuilt, "Bernstein reconstruction"]; b];
matrixCheck[actual_, expected_, label_, positive_:False] := (
  must[Dimensions[actual] === Dimensions[expected], label];
  zero[actual - expected, label];
  If[positive, must[And @@ (TrueQ[# >= 0] & /@ Flatten[actual]), label]]);

(* Block 2. The multiplier identity; f2 denotes F^2.
   Current notation: jGradS = J = Gamma(S,S); ds = K = Delta S.
   The earlier arXiv program called jGradS kk; no formula is changed. *)
rr = u + v; ss = 4 u v + t;
jGradS = 16 (4 u v rr + (u + 4 v) t);
mm = 4 u^2 + 12 u v + 4 v^2 + t;
ds = 16 m v + 16 (q + 1) u;
dim = 2 m + q + 1; deg = 3 + 4 ga;
cc = ga/ss; ee = ga (ga - 1)/ss^2;
vv = ss + 32 cc rr f2 + cc^2 jGradS f2;
dd = ga (ds + 32 rr)/ss + ga (ga - 1) jGradS/ss^2;
hssUU = 32 u v^2 + 32 u^2 v + 8 u t - 16 v t + 96 f2;
hssUZ = 128 u^2 + 640 u v + 128 v^2 + 96 t;
hssZZ = 512 u v^3 + 2048 u^2 v^2 + 512 u^3 v
  + 1280 v^2 t + 1920 u v t + 128 u^2 t + 64 t^2;
qq = dd vv - (8 rr + cc jGradS + 32 cc^2 f2 mm)
  - 2 cc (16 rr + cc jGradS) (ss + 16 cc rr f2)
  - cc (hssUU + 2 cc f2 hssUZ + cc^2 f2 hssZZ)
  - ee f2 (16 rr + cc jGradS)^2;
vhat = vv + 4 al (deg + al) f2/rr;
dhat = dd + 2 al (2 deg + dim + 2 al - 2)/rr;
hhat = vv dd - qq + 4 al (deg - 1) vv/rr
  + 4 al^2 deg (deg - 1) f2/rr^2
  + 4 al (deg + 2 al) (vv + 2 al deg f2/rr)/rr
  + 2 al vhat/rr + 4 al (al - 1) (deg + 2 al)^2 f2/rr^2;
qhat = qq + 4 al (deg + al) f2 dd/rr
  + 2 al (dim - 2 deg - 2 al - 1) vv/rr
  + 4 al (deg + al) (deg + 2 al (dim - 1)) f2/rr^2;
zero[vhat dhat - hhat - qhat, "Radial identity"];
zero[D[qhat, {f2, 2}], "Affine dependence on F^2"];
zero[D[qhat, m] - vhat (16 ga v/ss + 4 al/rr), "Derivative in m"];
zero[D[qhat, q] - vhat (16 ga u/ss + 2 al/rr), "Derivative in q"];
report["Master quotient, affine dependence, and propagation: PASS"];

(* Block 3. The three Bernstein matrices and the power obstruction. *)
normalized = {u -> lam, v -> 1 - lam, t -> theta lam^2,
  f2 -> omega theta lam^2 (1 - lam)};
hh = 4 (1 - lam) + theta lam;
pg1 = Together[qhat /. {al -> 0, ga -> 3/4} /. normalized];
pg2 = Together[qhat /. {al -> 1/2, ga -> 1/2} /. normalized];
n0 = 93 lam^2 theta - 48 lam^2 - 128 lam theta + 236 lam - 188
  + m (-12 lam^2 theta + 48 lam^2 + 12 lam theta - 96 lam + 48)
  + q (12 lam^2 theta - 48 lam^2 + 48 lam);
zero[(pg1 /. omega -> 0) hh - n0, "G1 endpoint numerator"];
bg1 = bernstein[(pg1 /. {m -> 8, q -> 3, omega -> 1}) hh^3,
  lam, theta, 4, 3];
matrixCheck[bg1, {{3136, 7448, 11760, 16072}, {16, 1432, 3293, 5599},
  {0, 1108/3, 2855/3, 1747}, {0, 0, 277/2, 1405/4}, {0, 0, 0, 1}},
  "G1 Bernstein matrix", True];
n083 = Expand[n0 /. {m -> 8, q -> 3}];
zero[n083 - lam (33 lam - 32) theta - 192 lam^2 + 388 lam - 196,
  "G1 base numerator"];
zero[(n083 /. theta -> 1) - (15 lam - 14)^2, "G1 square"];
zero[(n083 /. theta -> 0) - 4 (lam - 1) (48 lam - 49), "G1 boundary"];
zero[(pg1 /. {m -> 4, q -> 4, omega -> 0, lam -> 1/3, theta -> 1}) + 1,
  "G1 failure at (4,4)"];
m0 = Cancel[(pg2 /. {m -> 4, q -> 4, omega -> 0}) hh];
zero[m0 - lam (lam^2 theta^2 - 8 lam^2 theta + 16 lam^2
  + 64 lam theta - 64 lam - 48 theta + 48), "G2 endpoint numerator"];
b0 = bernstein[m0, lam, theta, 3, 2];
matrixCheck[b0, {{0, 0, 0}, {16, 8, 0}, {32/3, 16/3, 0}, {0, 4, 9}},
  "G2 omega=0 Bernstein matrix", True];
m1 = Cancel[(pg2 /. {m -> 4, q -> 4, omega -> 1}) hh^3];
b1 = bernstein[m1, lam, theta, 6, 4];
matrixCheck[b1, {{0, 720, 1440, 2160, 2880},
  {128, 864, 4880/3, 2416, 3232},
  {256/3, 4112/5, 14648/9, 12504/5, 51632/15},
  {128/5, 1464/5, 10976/15, 6726/5, 2136},
  {0, 976/15, 1064/5, 7379/15, 4756/5},
  {0, 0, 244/9, 302/3, 814/3}, {0, 0, 0, 2, 9}},
  "G2 omega=1 Bernstein matrix", True];
p43 = Together[qhat /. {m -> 4, q -> 3} /. normalized];
edge = Cancel[p43 /. {theta -> 1, omega -> 0}];
zero[Limit[edge, lam -> 0, Direction -> -1] + 8 (ga - 1) (2 ga - 1),
  "Power obstruction at lambda=0"];
zero[(p43 /. {lam -> 2/3, theta -> 1, omega -> 0})
  + 8/3 (2 al^2 + (8 ga - 5) al + 8 ga^2 - 6 ga + 3),
  "Power obstruction at lambda=2/3"];
report["Three Bernstein matrices and power obstruction: PASS"];

(* Block 4. Jacobi integrals: Gamma values and rational recurrences.
   jLink is the restriction of J = Gamma(S,S), not K = Delta S.
   JacobiMatrix is the bilinear-form matrix M_J, not the scalar J. *)
sj = 4 u (1 - u) + u^2 theta;
jLink = 16 (4 u (1 - u) + (u + 4 (1 - u)) u^2 theta);
integrand[a_, b_] := Expand[(a b + 1/2) jLink sj^(a + b - 3/2)
  - 32 sj^(a + b - 1/2) + (81/4 - 16 a b) sj^(a + b + 1/2)];
momentGamma[i_, j_] := Gamma[i + 3] Gamma[3/2]/Gamma[i + 9/2]
  Gamma[j + 3/2] Gamma[1/2]/(Pi Gamma[j + 2]);
um[0] = 16/105; tm[0] = 1/2;
um[i_Integer] /; i > 0 := um[i] = um[i - 1] 2 (i + 2)/(2 i + 7);
tm[j_Integer] /; j > 0 := tm[j] = tm[j - 1] (2 j + 1)/(2 j + 2);
momentRecurrence[i_, j_] := um[i] tm[j];
weighted[poly_, moment_] := Module[{a = CoefficientList[poly, {u, theta}]},
  Total[Flatten[Table[a[[i, j]] moment[i - 1, j - 1],
    {i, Length[a]}, {j, Length[a[[1]]]}]]]];
exponents = {3/4, 7/4};
jacobi = Table[weighted[integrand[a, b], momentGamma],
  {a, exponents}, {b, exponents}];
jacobiRec = Table[weighted[integrand[a, b], momentRecurrence],
  {a, exponents}, {b, exponents}];
expectedJacobi = {{96/715, 471936/1616615},
  {471936/1616615, 571904/969969}};
matrixCheck[jacobi, expectedJacobi, "Jacobi matrix"];
matrixCheck[jacobiRec, expectedJacobi, "Independent Jacobi moments"];
direction = {1, -1/2};
zero[direction.jacobi.direction + 9952/969969, "Negative Jacobi form"];
report["Jacobi matrix and independent rational moments: PASS"];
report["Negative Jacobi value: ", direction.jacobi.direction];

(* Block 5. The quaternionic Clifford model and coordinate contractions. *)
id4 = IdentityMatrix[4]; zz4 = ConstantArray[0, {4, 4}];
ei = {{0, -1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, -1}, {0, 0, 1, 0}};
ej = {{0, 0, -1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}, {0, -1, 0, 0}};
block[a_, b_, c_, d_] := Join[MapThread[Join, {a, b}],
  MapThread[Join, {c, d}]];
matrices = Join[{block[id4, zz4, zz4, -id4],
  block[zz4, id4, id4, zz4]}, block[zz4, #, -#, zz4] & /@ {ei, ej, ei.ej}];
Do[zero[matrices[[i]] - Transpose[matrices[[i]]], "Symmetry"];
  zero[Tr[matrices[[i]]], "Trace"];
  Do[zero[matrices[[i]].matrices[[j]] + matrices[[j]].matrices[[i]]
    - 2 If[i == j, 1, 0] IdentityMatrix[8], "Clifford relation"],
    {j, 5}], {i, 5}];
report["Quaternionic Clifford relations: PASS"];
yExamples = {{1, 2, -1, 0, 2, 0, 1, -2}, {1, 0, 0, 0, 1, 0, 0, 0},
  {1, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 1, 1, 0, 0, 0}};
zExamples = {{1, -2, 1, 3, -1}, {1, 1, 0, 0, 0}, {1, 0, 0, 0, 0},
  {1, 2, 0, -1, 1}, {0, 0, 0, 0, 0}};
gradient[f_, vars_] := D[f, #] & /@ vars;
hessian[f_, vars_] := Table[D[f, a, b], {a, vars}, {b, vars}];
Do[Module[{ys = Array[y, 8], zs = Array[z, qv + 1], vars, tau,
  cubic, spoly, gf, gs, hf, hs, yval, zval, rules, fval, uv, vv0,
  tv, rv, sv, jValue, gfv, gsv, hfv, hsv, residuals},
  vars = Join[ys, zs];
  tau = Table[ys.matrices[[i]].ys, {i, qv + 1}];
  cubic = Expand[tau.zs]; gf = gradient[cubic, vars];
  spoly = Expand[gf.gf]; gs = gradient[spoly, vars];
  hf = hessian[cubic, vars]; hs = hessian[spoly, vars];
  Do[yval = yExamples[[j]]; zval = Take[zExamples[[j]], qv + 1];
    rules = Thread[vars -> Join[yval, zval]];
    fval = cubic /. rules; uv = yval.yval; vv0 = zval.zval;
    tv = (tau.tau) /. rules; rv = uv + vv0; sv = 4 uv vv0 + tv;
    jValue = 16 (4 uv vv0 rv + (uv + 4 vv0) tv);
    {gfv, gsv, hfv, hsv} = {gf, gs, hf, hs} /. rules;
    residuals = {gfv.gfv - sv, gsv.gsv - jValue, gfv.gsv - 16 rv fval,
      Tr[hfv], Tr[hsv] - 64 vv0 - 16 (qv + 1) uv,
      gfv.hfv.gfv - 8 rv fval, gfv.hfv.gsv - jValue/2,
      gsv.hfv.gsv - 32 fval (4 uv^2 + 12 uv vv0 + 4 vv0^2 + tv),
      gfv.hsv.gfv - (32 uv vv0^2 + 32 uv^2 vv0 + 8 uv tv - 16 vv0 tv
        + 96 fval^2),
      gfv.hsv.gsv - fval (128 uv^2 + 640 uv vv0 + 128 vv0^2 + 96 tv),
      gsv.hsv.gsv - (512 uv vv0^3 + 2048 uv^2 vv0^2 + 512 uv^3 vv0
        + 1280 vv0^2 tv + 1920 uv vv0 tv + 128 uv^2 tv + 64 tv^2)};
    Do[zero[residuals[[k]], StringJoin["Ambient q=", ToString[qv],
      ", point=", ToString[j], ", identity=", ToString[k]]],
      {k, Length[residuals]}], {j, 5}];
  report["55 ambient contractions at (m,q)=(4,", qv, "): PASS"]],
  {qv, {3, 4}}];

(* Block 6. Symbolic Hopf identities on the unit sphere. *)
xs = {x1, x2, x3, x4, x5, x6, x7, x8}; norm2 = xs.xs;
hopf = Table[xs.matrices[[i]].xs, {i, 5}];
zero[hopf.hopf - norm2^2, "Hopf norm"];
sphereGrad = Table[2 (matrices[[i]].xs - hopf[[i]] xs), {i, 5}];
Do[zero[sphereGrad[[i]].sphereGrad[[j]]
  - 4 (If[i == j, 1, 0] - hopf[[i]] hopf[[j]])
  - 4 (norm2 - 1) (If[i == j, 1, 0] + hopf[[i]] hopf[[j]]),
  "Hopf spherical gradient identity"], {i, 5}, {j, i, 5}];
report["Symbolic Hopf norm and spherical gradient identities: PASS"];

 BernsteinMatrices={bg1,b0,b1};
 JacobiMatrix=jacobi;
 NegativeJacobiValue=direction.jacobi.direction;
 True], "CliffordClassificationFailure"]
];


RunChapterChecks[] := Block[{$IterationLimit = Infinity},
  Module[{started, failures},
    started = DateString[];
    VerificationLog = {}; VerificationGroupLog = {};
    VerificationMessages = {};
    Print[VerificationTitle];
    Print["Kernel: ", $Version];
    RunVerificationGroup["CheckCliffordClassification", CheckCliffordClassification];
    failures = Select[VerificationLog, !TrueQ[#[[2]]] &];
    VerificationResults = Length[failures] == 0 &&
      Length[VerificationLog] == VerificationExpectedCount &&
      Length[VerificationMessages] == 0;
    VerificationMetadata = {
      {"Code release", "20260923-v6"},
      {"Notation revision", "20260928-v3"},
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
MergedVerificationV5`Chapter08`VerificationResults
