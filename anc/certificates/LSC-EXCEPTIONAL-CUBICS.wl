(* ::Package:: *)
(* LSC-EXCEPTIONAL-CUBICS -- Exceptional Cubics.
   Compatible with Mathematica 12. No external package or file dependencies.
   Get returns checks; only the central runner writes execution records.

   Formal drift/certificate identities and explicit G2 coordinate calculations; rank-stratum and removability arguments remain analytic.
*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`ExceptionalCubics`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`ExceptionalCubics`Private`*"];

VerificationLog = {};
EvaluationMessages = {};
EvaluationOK = Block[{$IterationLimit = Infinity, $Assumptions = True, System`$MessageList = {}},
  ExecutionStatus = CheckAbort[Catch[(
residual[expression_] := Factor[Together[PowerExpand[Expand[expression]]]];
check[label_String, expression_] := Module[{values, success},
  values = residual /@ Flatten[{expression}];
  success = And @@ (TrueQ[# === 0] & /@ values);
  AppendTo[VerificationLog, {label, success, If[success, "0", ToString[values, InputForm]]}];
  
  If[!success, Print[InputForm[values]]]; success
];

(* Free generators, of degrees 2, 4, 6, 3. *)
vars = {rr, ss, dd, ff};
metric = {{4 rr, 8 ss, 12 dd, 6 ff},
  {8 ss, 4 rr ss + 12 dd, 8 rr dd, 4 rr ff},
  {12 dd, 8 rr dd, 4 dd ss, 2 ss ff},
  {6 ff, 4 rr ff, 2 ss ff, ss}};
gradient[expression_, coordinates_List] := D[expression, #] & /@ coordinates;
gamma[left_, right_] := gradient[left, vars].metric.gradient[right, vars];
lap[expression_, drift_List] := gradient[expression, vars].drift +
  Total[Flatten[Table[With[{xi = vars[[i]], xj = vars[[j]]},
    metric[[i, j]] D[expression, xi, xj]], {i, 4}, {j, 4}]]];
defect[expression_, drift_List] := Module[{norm = Expand[gamma[expression, expression]]},
  norm lap[expression, drift] - gamma[expression, norm]/2];

gg = ss ff;
normSF = ss^3 + 12 ff^2 (rr ss + dd);
polySF = (cc - 16) rr ss^3 - 30 dd ss^2 +
  ff^2 ((12 cc - 72) rr dd + 12 cc rr^2 ss - 96 ss^2);
positiveSF = (cc - 58/3) rr ss^3 + (10/3) ss^2 (rr ss - 9 dd) +
  ff^2 ((12 cc - 72) rr dd + 12 cc ss (rr^2 - 3 ss) + (36 cc - 96) ss^2);
check["SF gradient norm", gamma[gg, gg] - normSF];
check["SF Laplacian", lap[gg, {dr, cc rr, dt, 0}] - (cc + 8) rr ff];
check["SF defect, with arbitrary Delta R and Delta D",
  defect[gg, {dr, cc rr, dt, 0}] - ff polySF];
check["SF positive decomposition", polySF - positiveSF];
check["SF threshold at balanced spectrum and F=0",
  (polySF /. {rr -> 3, ss -> 3, dd -> 1, ff -> 0}) - 81 (cc - 58/3)];

ggRoot = ff Sqrt[rr ss];
normRoot = rr ss^2 + ff^2 (3 dd rr/ss + 5 rr^2 + 11 ss);
polyRoot = 6 rr ss^2 (ss^2 - 2 rr dd) + 15 dd ff^2 rr (rr^2 - ss) +
  45 ff^2 rr^4 ss + 123 ff^2 rr^2 ss^2 + 242 ff^2 ss^3;
check["sqrt(RS)F gradient norm", gamma[ggRoot, ggRoot] - normRoot];
check["sqrt(RS)F certificate for Jordan d=2",
  defect[ggRoot, {36, 16 rr, 4 ss, 0}] rr ss^2 - ggRoot polyRoot];
transfer = (dim - 2) normRoot (12 D[ggRoot, rr] + 4 rr D[ggRoot, ss]);
check["Jordan drift transfer",
  defect[ggRoot, {12 (dim + 1), 4 (dim + 2) rr, 4 ss, 0}] -
  defect[ggRoot, {36, 16 rr, 4 ss, 0}] - transfer];
check["Sign of the transfer bracket",
  12 D[ggRoot, rr] + 4 rr D[ggRoot, ss] - ggRoot (6/rr + 2 rr/ss)];

(* Independent Gram-coordinate calculation for an n by 3 real matrix. *)
gramVars = {b11, b22, b33, b12, b13, b23};
gram = {{b11, b12, b13}, {b12, b22, b23}, {b13, b23, b33}};
pairs = {{1, 1}, {2, 2}, {3, 3}, {1, 2}, {1, 3}, {2, 3}};
gramMetric = Table[With[{i = pairs[[alpha, 1]], j = pairs[[alpha, 2]],
    k = pairs[[beta, 1]], l = pairs[[beta, 2]]},
  If[i == k, 1, 0] gram[[j, l]] + If[i == l, 1, 0] gram[[j, k]] +
  If[j == k, 1, 0] gram[[i, l]] + If[j == l, 1, 0] gram[[i, k]]],
  {alpha, 6}, {beta, 6}];
gramDrift = {2 nn, 2 nn, 2 nn, 0, 0, 0};
rGram = Tr[gram];
sGram = b11 b22 + b11 b33 + b22 b33 - b12^2 - b13^2 - b23^2;
dGram = Det[gram];
gramGamma[left_, right_] := gradient[left, gramVars].gramMetric.gradient[right, gramVars];
gramLap[expression_] := gradient[expression, gramVars].gramDrift +
  Total[Flatten[Table[With[{xi = gramVars[[i]], xj = gramVars[[j]]},
    gramMetric[[i, j]] D[expression, xi, xj]], {i, 6}, {j, 6}]]];
check["Gram Delta R", gramLap[rGram] - 6 nn];
check["Gram Delta S", gramLap[sGram] - 4 (nn - 1) rGram];
check["Gram Delta D", gramLap[dGram] - 2 (nn - 2) sGram];
gramChecks = {{"RR", rGram, rGram, 4 rGram}, {"RS", rGram, sGram, 8 sGram},
  {"RD", rGram, dGram, 12 dGram}, {"SS", sGram, sGram, 4 rGram sGram + 12 dGram},
  {"SD", sGram, dGram, 8 rGram dGram}, {"DD", dGram, dGram, 4 sGram dGram}};
Do[check["Gram Gamma " <> item[[1]], gramGamma[item[[2]], item[[3]]] - item[[4]]],
  {item, gramChecks}];
check["G2 Delta D = 10S", (gramLap[dGram] /. nn -> 7) - 10 sGram];
check["Alternating F: Gamma(F,S)=4RF",
  2 ff Total[D[sGram, #] & /@ Take[gramVars, 3]] - 4 rGram ff];
check["Alternating F: Gamma(F,D)=2SF",
  2 ff Total[D[dGram, #] & /@ Take[gramVars, 3]] - 2 sGram ff];

(* Compact G2 three-form: 123+145+167+246-257-347-356. *)
triples = {{{1, 2, 3}, 1}, {{1, 4, 5}, 1}, {{1, 6, 7}, 1},
  {{2, 4, 6}, 1}, {{2, 5, 7}, -1}, {{3, 4, 7}, -1}, {{3, 5, 6}, -1}};
aa = Array[a, 7]; bb = Array[b, 7]; zz = Array[z, 7];
coordinates = Join[aa, bb, zz];
signedPermutations = {{1, 2, 3, 1}, {1, 3, 2, -1}, {2, 1, 3, -1},
  {2, 3, 1, 1}, {3, 1, 2, 1}, {3, 2, 1, -1}};
fG2 = Total[Flatten[Table[With[{triple = item[[1]], sign = item[[2]]},
  sign perm[[4]] aa[[triple[[perm[[1]]]]]] bb[[triple[[perm[[2]]]]]]
    zz[[triple[[perm[[3]]]]]]], {item, triples}, {perm, signedPermutations}]]];
sG2 = (aa.aa) (bb.bb) - (aa.bb)^2 + (aa.aa) (zz.zz) - (aa.zz)^2 +
  (bb.bb) (zz.zz) - (bb.zz)^2;
gradFG2 = gradient[fG2, coordinates];
check["G2 coordinate gradient norm = Gram e2", gradFG2.gradFG2 - sG2];
check["G2 harmonicity", Total[D[fG2, {#, 2}] & /@ coordinates]];

(* Four-parameter slice a=r e1, b=s e2, z=u e3+v e4. *)
slice = {r^2 + s^2 + u^2 + v^2,
  r^2 s^2 + (r^2 + s^2) (u^2 + v^2), r^2 s^2 (u^2 + v^2), r s u};
jacobian = Det[gradient[#, {r, s, u, v}] & /@ slice];

check["G2 independence: explicit Jacobian formula",
  jacobian + 8 r^2 s^2 v (r^2 - s^2) (r^2 - u^2 - v^2) (s^2 - u^2 - v^2)];

rankOneRules = Join[{a[1] -> r}, Thread[Rest[coordinates] -> 0]];
hessian = (gradient[#, coordinates] & /@ gradFG2) /. rankOneRules;
check["G2 rank-one Hessian cubic identity", hessian.hessian.hessian - r^2 hessian];
check["G2 rank-one Hessian trace", Tr[hessian]];
check["G2 rank-one Hessian squared trace", Tr[hessian.hessian] - 12 r^2];
check["G2 rank-one Hessian characteristic polynomial",
  Det[lambda IdentityMatrix[21] - hessian] - lambda^9 (lambda^2 - r^2)^6];


    ); True, _, ($Failed &)], $Aborted];
  EvaluationMessages = ToString[#, InputForm] & /@ System`$MessageList;
  ExecutionStatus === True && Length[EvaluationMessages] == 0
];

CertificateResult = <|"CertificateID" -> "LSC-EXCEPTIONAL-CUBICS",
  "ExpectedChecks" -> 28,
  "Checks" -> (<|"Name" -> #[[1]], "Passed" -> TrueQ[#[[2]]],
      "Detail" -> If[StringQ[#[[3]]], #[[3]], ToString[#[[3]], InputForm]]|> &) /@ VerificationLog,
  "EvaluationOK" -> TrueQ[EvaluationOK],
  "EvaluationMessages" -> EvaluationMessages|>;
End[];
EndPackage[];
LaplacianCertificates`ExceptionalCubics`Private`CertificateResult
