(* ::Package:: *)
(* Independent supplemental unit, release 20260923-v6.
   New Chapter 7: cubic algebras, general quadratic weight, and variable-K SF.
   Mathematica 12 / Mathics; exact polynomial arithmetic; no PowerExpand.
   NOT part of the historical 13-unit / 550-assertion entry point.
   Read README_Chapter_07_General_Cubic_Algebra.md for scope and hypotheses. *)
ClearAll["MergedVerificationV6`Chapter07General`*"];
Begin["MergedVerificationV6`Chapter07General`"];
VerificationTitle = "Chapter 7 supplemental: general cubic algebra and variable-K SF";
VerificationFileStem = "Verify_Chapter_07_General_Cubic_Algebra";
VerificationExpectedCount = 31;
VerificationDirectory = If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
  DirectoryName[$InputFileName], Directory[]];
VerificationLog = {}; VerificationGroupLog = {}; VerificationMessages = {};

CheckEqual[lhs_, rhs_, label_String] := Module[{residuals, passed},
  residuals = (Factor[Together[Expand[#]]] &) /@ Flatten[{lhs - rhs}];
  passed = And @@ ((# === 0) & /@ residuals);
  AppendTo[VerificationLog,
    {label, passed, If[passed, "0", ToString[residuals, InputForm]]}];
  If[!passed, Print["FAIL: ", label, " : ", InputForm[residuals]]];
  passed
];
CheckTrue[condition_, label_String] := Module[{passed = TrueQ[condition]},
  AppendTo[VerificationLog,
    {label, passed, If[passed, "True", ToString[condition, InputForm]]}];
  If[!passed, Print["FAIL: ", label, " : ", InputForm[condition]]];
  passed
];
RunGroup[label_String, fun_Symbol] := Module[{before, status},
  before = Length[VerificationLog];
  status = Block[{System`$MessageList = {}}, Module[{result},
    result = CheckAbort[Catch[fun[]; True, _, ($Failed &)], $Aborted];
    TrueQ[result] && Length[System`$MessageList] == 0
  ]];
  If[!TrueQ[status], AppendTo[VerificationMessages, label]];
  AppendTo[VerificationGroupLog, {label, Length[VerificationLog] - before, TrueQ[status]}];
  Print["Group completed: ", label, "; checks: ", Length[VerificationLog] - before]
];
GradExpr[p_, xs_List] := D[p, #] & /@ xs;
HessExpr[p_, xs_List] := Table[D[p, xs[[i]], xs[[j]]], {i, Length[xs]}, {j, Length[xs]}];
LapExpr[p_, xs_List] := Tr[HessExpr[p, xs]];
CurvExpr[p_, xs_List] := Module[{g = GradExpr[p, xs], h = HessExpr[p, xs]}, (g.g) Tr[h] - g.h.g];
GammaExpr[p_, q_, xs_List, gram_List] := GradExpr[p, xs].gram.GradExpr[q, xs];
ReducedLapExpr[p_, xs_List, gram_List, drift_List] := drift.GradExpr[p, xs] +
  Total[Flatten[Table[gram[[i, j]] D[p, xs[[i]], xs[[j]]],
    {i, Length[xs]}, {j, Length[xs]}]]];
NonnegativeCoefficients[p_, xs_List] :=
  And @@ (TrueQ[# >= 0] & /@ Flatten[CoefficientList[Expand[p], xs]]);

(* Generic coefficient checks in three coordinates support, but do not replace,
   the dimension-independent proofs. All ten cubic coefficients are free. *)
CheckCubicDictionary[] := Module[
  {xs = {x, y, z}, p, tensor, mul, xx, xxx, lx, vv = {v1, v2, v3}, lv, ss, qform, i, j, k},
  p = c300 x^3 + c030 y^3 + c003 z^3 + c210 x^2 y + c201 x^2 z +
    c120 x y^2 + c021 y^2 z + c102 x z^2 + c012 y z^2 + c111 x y z;
  tensor = Table[D[p, xs[[i]], xs[[j]], xs[[k]]], {i, 3}, {j, 3}, {k, 3}];
  mul[u_List, v_List] := Table[Total[Flatten[Table[tensor[[i, j, k]] u[[j]] v[[k]],
    {j, 3}, {k, 3}]]], {i, 3}];
  xx = mul[xs, xs]; xxx = mul[xs, xx];
  lx = Table[Total[Table[tensor[[i, j, k]] xs[[k]], {k, 3}]], {i, 3}, {j, 3}];
  lv = Table[Total[Table[tensor[[i, j, k]] vv[[k]], {k, 3}]], {i, 3}, {j, 3}];
  ss = GradExpr[p, xs].GradExpr[p, xs];
  CheckEqual[p, xx.xs/6, "C07G-01 cubic dictionary F = <x^2,x>/6"];
  CheckEqual[GradExpr[p, xs], xx/2, "C07G-02 cubic dictionary gradient = x^2/2"];
  CheckEqual[HessExpr[p, xs], lx, "C07G-03 cubic dictionary Hessian = L_x"];
  CheckEqual[LapExpr[p, xs], Tr[lx], "C07G-04 cubic dictionary Laplacian = trace L_x"];
  CheckEqual[ss, xx.xx/4, "C07G-05 gradient square = |x^2|^2/4"];
  CheckEqual[GradExpr[ss, xs], xxx, "C07G-06 gradient S = x^3"];
  CheckEqual[LapExpr[ss, xs], 2 Tr[lx.lx] + 2 GradExpr[p, xs].GradExpr[LapExpr[p, xs], xs],
    "C07G-07 general Bochner identity including nonharmonic term"];
  CheckEqual[vv.GradExpr[p, xs], xs.lv.xs/2,
    "C07G-08 translation derivative and annihilator identity"];
  qform = q11 x^2 + q22 y^2 + q33 z^2 + q12 x y + q13 x z + q23 y z;
  CheckEqual[4 (CurvExpr[p, xs] + qform p),
    (xx.xx) Tr[lx] - xx.xxx + (2/3) qform (xx.xs),
    "C07G-09 general quadratic-weight cubic-algebra residual"]
];

CheckQuadraticOperators[] := Module[
  {xs = {x, y, z}, bb, cc, dd, qb, qc, p, eb, ec, bracket, i, j},
  bb = {{b11, b12, b13}, {b12, b22, b23}, {b13, b23, b33}};
  cc = {{c11, c12, c13}, {c12, c22, c23}, {c13, c23, c33}};
  qb = xs.bb.xs; qc = xs.cc.xs;
  CheckEqual[GradExpr[qb, xs].GradExpr[qc, xs], 4 xs.((bb.cc + cc.bb)/2).xs,
    "C07G-10 quadratic Gamma and Jordan product"];
  p = p300 x^3 + p030 y^3 + p003 z^3 + p210 x^2 y + p201 x^2 z +
    p120 x y^2 + p021 y^2 z + p102 x z^2 + p012 y z^2 + p111 x y z;
  eb[h_] := (bb.xs).GradExpr[h, xs]; ec[h_] := (cc.xs).GradExpr[h, xs];
  bracket = bb.cc - cc.bb;
  CheckEqual[eb[ec[p]] - ec[eb[p]], -(bracket.xs).GradExpr[p, xs],
    "C07G-11 Euler-operator commutator sign"];
  dd = DiagonalMatrix[{d1, d2, d3}];
  CheckEqual[dd.(dd.cc - cc.dd) - (dd.cc - cc.dd).dd,
    Table[({d1, d2, d3}[[i]] - {d1, d2, d3}[[j]])^2 cc[[i, j]], {i, 3}, {j, 3}],
    "C07G-12 diagonal double-commutator entries"]
];

(* General q: q0, q_i, q_ij are independent; A is a general symmetric matrix.
   The Hessian polynomial is a conditional residual identity. *)
CheckCriticalJets[] := Module[
  {xs = {x, y, z}, aa, h, cubic, gg, hh, lh, qjet, residualMatrix,
   multiplier, g2, gradientJet, hessianJet, leading, tt, a0, q0},
  aa = {{a11, a12, a13}, {a12, a22, a23}, {a13, a23, a33}};
  h = xs.aa.xs/2;
  cubic = t300 x^3 + t030 y^3 + t003 z^3 + t210 x^2 y + t201 x^2 z +
    t120 x y^2 + t021 y^2 z + t102 x z^2 + t012 y z^2 + t111 x y z;
  gg = tt (aa.xs) + tt^2 GradExpr[cubic, xs];
  hh = aa + tt HessExpr[cubic, xs]; lh = CurvExpr[h, xs];
  CheckEqual[Coefficient[Expand[(gg.gg) Tr[hh] - gg.hh.gg], tt, 2], lh,
    "C07G-13 critical cubic curvature quadratic Taylor coefficient"];
  CheckEqual[lh, Tr[aa] xs.(aa.aa).xs - xs.(aa.aa.aa).xs,
    "C07G-14 quadratic Hessian curvature identity"];
  qjet = q0 + tt (q1 x + q2 y + q3 z) + tt^2
    (q11 x^2 + q22 y^2 + q33 z^2 + q12 x y + q13 x z + q23 y z);
  CheckEqual[Coefficient[Expand[qjet (tt^2 h + tt^3 cubic)], tt, 2], q0 h,
    "C07G-15 general quadratic q contributes q(c) h at order two"];
  residualMatrix = aa.aa.aa - Tr[aa] aa.aa - q0 aa/2;
  CheckEqual[lh + q0 h, -xs.residualMatrix.xs,
    "C07G-16 general-q Hessian polynomial residual"];
  multiplier = a0 + tt (m1 x + m2 y + m3 z) + tt^2
    (m11 x^2 + m22 y^2 + m33 z^2 + m12 x y + m13 x z + m23 y z);
  g2 = Coefficient[Expand[multiplier (tt^2 h + tt^3 cubic)], tt, 2];
  CheckEqual[g2, a0 h, "C07G-17 smooth-multiplier quadratic Taylor jet"];
  (* Only the displayed leading jets affect order four; higher jets are free. *)
  gradientJet = tt a0 (aa.xs) + tt^2 {u1, u2, u3};
  hessianJet = a0 aa + tt {{v11, v12, v13}, {v12, v22, v23}, {v13, v23, v33}};
  leading = Coefficient[Expand[(tt^2 h + tt^3 cubic)
    ((gradientJet.gradientJet) Tr[hessianJet] - gradientJet.hessianJet.gradientJet)], tt, 4];
  CheckEqual[leading + a0^3 q0 h^2, -a0^3 h xs.residualMatrix.xs,
    "C07G-18 positive-multiplier order-four sign under Hessian relation"]
];

CheckVariableTraceSF[] := Module[
  {xs = {r, f, s}, gram, drift, g, q, lapg, curvature, defect, positive,
   kgap, sgap, t, k, n, kval, qPositive, positiveSlack, fsquared},
  gram = {{4 r, 6 f, 8 s}, {6 f, s, 4 r f}, {8 s, 4 r f, 4 r s + 12 f^2}};
  (* kval denotes the value of the quadratic K=Delta S at a point.
     It is independent of R. No Gamma(K,.) or Delta K is needed. *)
  drift = {2 n, 0, kval}; g = s f; q = s^3 + 12 r s f^2 + 12 f^4;
  lapg = ReducedLapExpr[g, xs, gram, drift];
  curvature = GammaExpr[g, g, xs, gram] lapg - GammaExpr[g, GammaExpr[g, g, xs, gram], xs, gram]/2;
  defect = (kval - 16 r) s^3 + f^2 s (12 r kval - 126 s) + 12 f^4 (kval - 6 r);
  positive = (kval - 16 r) q + 150 r^2 s f^2 + 126 f^2 s (r^2/3 - s) + 120 r f^4;
  CheckEqual[GammaExpr[f, f, xs, gram] ReducedLapExpr[f, xs, gram, drift] -
    GammaExpr[f, GammaExpr[f, f, xs, gram], xs, gram]/2, -2 r f,
    "C07G-19 gradient identities imply radial weight q = 2R"];
  CheckEqual[GammaExpr[g, g, xs, gram], q, "C07G-20 variable-K SF gradient square"];
  CheckEqual[lapg, (kval + 8 r) f, "C07G-21 variable-K SF Laplacian"];
  CheckEqual[Table[GammaExpr[g, xs[[i]], xs, gram], {i, 3}],
    {14 s f, s^2 + 4 r f^2, 8 r s f + 12 f^3}, "C07G-22 SF contractions with R F S"];
  CheckEqual[GradExpr[q, xs], {12 s f^2, 24 r s f + 48 f^3, 3 s^2 + 12 r f^2},
    "C07G-23 SF gradient-square partial derivatives"];
  CheckEqual[curvature, f defect, "C07G-24 variable-K SF curvature polynomial"];
  CheckEqual[defect, positive, "C07G-25 variable-K SF positive decomposition"];
  CheckEqual[defect /. f -> 0, (kval - 16 r) s^3,
    "C07G-26 quotient polynomial extension at F = 0"];
  CheckEqual[(defect /. kval -> 4 k r) - (defect /. kval -> 16 r), 4 (k - 4) r q,
    "C07G-27 radial specialization K = 4kR"];
  (* For F!=0 divide the first stationary equation by F; fsquared denotes F^2.
     Neither sign of F is assumed. Branch existence/exhaustion is a text proof. *)
  CheckEqual[{12 s - 4, 16 s^2 - 4 s - 12 fsquared} /.
    {s -> 1/3, fsquared -> 1/27}, {0, 0},
    "C07G-28 positive-S spherical maximum branch F nonzero"];
  CheckEqual[{12 s f - 4 f, 16 s^2 - 4 s - 12 f^2} /. {s -> 1/4, f -> 0}, {0, 0},
    "C07G-29 positive-S spherical maximum branch F zero"];
  (* Sign domain: r>=0, s>=0, t=f^2>=0, kgap=K-16r>=0,
     sgap=r^2/3-s>=0. r>0 and s>0 when F!=0 are geometric/Euler inputs.
     These assertions test coefficients only, not the global hypotheses.
     The strictly positive summand when F!=0 is 150 r^2 s t. *)
  qPositive = s^3 + 12 r s t + 12 t^2;
  positiveSlack = kgap qPositive + 150 r^2 s t + 126 t s sgap + 120 r t^2;
  CheckTrue[NonnegativeCoefficients[qPositive, {r, s, t}],
    "C07G-30 Q nonnegative coefficients on stated nonnegative domain"];
  CheckTrue[NonnegativeCoefficients[positiveSlack, {r, s, t, kgap, sgap}],
    "C07G-31 SF decomposition nonnegative coefficients on stated slack domain"]
];

RunChecks[] := Block[{$IterationLimit = Infinity}, Module[{failures},
  VerificationLog = {}; VerificationGroupLog = {}; VerificationMessages = {};
  Print[VerificationTitle]; Print["Kernel: ", $Version];
  RunGroup["Cubic dictionary", CheckCubicDictionary];
  RunGroup["Quadratic operators", CheckQuadraticOperators];
  RunGroup["General-q critical jets", CheckCriticalJets];
  RunGroup["Variable-K SF", CheckVariableTraceSF];
  failures = Select[VerificationLog, !TrueQ[#[[2]]] &];
  VerificationMathPassed = Length[failures] == 0 &&
    Length[VerificationLog] == VerificationExpectedCount && Length[VerificationMessages] == 0;
  VerificationMetadata = {{"Code release", "20260923-v6"}, {"Kernel", $Version},
    {"Expected checks", VerificationExpectedCount}, {"Actual checks", Length[VerificationLog]},
    {"Passed checks", Length[VerificationLog] - Length[failures]}, {"Failed checks", Length[failures]},
    {"Evaluation-message groups", VerificationMessages}, {"Mathematical checks passed", VerificationMathPassed}};
  Print["Checks: ", Length[VerificationLog], "; passed: ", Length[VerificationLog] - Length[failures],
    "; failed: ", Length[failures], "; expected: ", VerificationExpectedCount];
  VerificationMathPassed
]];

EnsureOutputDirectory[dir_String] := Module[{ancestor = ExpandFileName[dir]},
  (* Avoid Return inside While: its scope differs in Mathics. *)
  While[!DirectoryQ[ancestor] && !FileExistsQ[ancestor] &&
      ancestor =!= DirectoryName[ancestor], ancestor = DirectoryName[ancestor]];
  If[!DirectoryQ[ancestor], False,
    If[DirectoryQ[dir], True,
      StringQ[CheckAbort[CreateDirectory[dir, CreateIntermediateDirectories -> True], $Aborted]] && DirectoryQ[dir]]]
];

ExportResults[dir_String] := Module[{txt, csv, field, row, txtBody, csvBody, written},
  If[!TrueQ[EnsureOutputDirectory[dir]], Return[False]];
  txt = FileNameJoin[{dir, VerificationFileStem <> "_Results.txt"}];
  csv = FileNameJoin[{dir, VerificationFileStem <> "_Results.csv"}];
  If[DirectoryQ[txt] || DirectoryQ[csv], Return[False]];
  field[item_] := "\"" <> StringReplace[If[StringQ[item], item, ToString[item, InputForm]], "\"" -> "\"\""] <> "\"";
  row[items_] := StringRiffle[ToString[#, InputForm] & /@ items, " | "];
  txtBody = StringRiffle[Join[{VerificationTitle, "Actual supplemental execution record", ""},
    row /@ VerificationMetadata, {"", "Groups: label | count | no evaluation messages"},
    row /@ VerificationGroupLog, {"", "Checks: label | passed | residual"}, row /@ VerificationLog,
    {"", "The returned source value additionally requires BOTH exports to succeed.",
     "Geometry, singular-set removability, and global positivity hypotheses are not CAS conclusions."}], "\n"] <> "\n";
  csvBody = StringRiffle[(StringRiffle[field /@ #, ","] &) /@
    Prepend[VerificationLog, {"label", "passed", "residual"}], "\r\n"] <> "\r\n";
  written = {CheckAbort[Export[txt, txtBody, "Text"], $Aborted], CheckAbort[Export[csv, csvBody, "Text"], $Aborted]};
  VerificationOutputFiles = written;
  written === {txt, csv} && FileExistsQ[txt] && FileExistsQ[csv]
];

RunChecks[];
VerificationExportSucceeded = TrueQ[ExportResults[
  FileNameJoin[{VerificationDirectory, "verification_results_v6_supplemental"}]]];
VerificationResults = TrueQ[VerificationMathPassed] && TrueQ[VerificationExportSucceeded];
Print["Exports succeeded: ", VerificationExportSucceeded];
Print["Complete and passed: ", VerificationResults];
End[];
MergedVerificationV6`Chapter07General`VerificationResults
