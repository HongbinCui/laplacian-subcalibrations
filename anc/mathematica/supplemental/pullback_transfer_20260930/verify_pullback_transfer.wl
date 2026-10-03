(* Supplement: exact differential identities for the 2026-09-30 pullback
   revision. Mathematica / Wolfram Language; no external packages required.
   Run: wolframscript -file verify_pullback_transfer.wl
   A native run creates mathematica_native_results.json in this directory.
   This supplied script has NOT been executed with a native Wolfram kernel.

   r is source R. wH is target W_H = |grad_y H|^2 evaluated at P.
   The source square |grad_x(H o P)|^2 equals 4 R W_H(P).
   targetq is the target squared radius.
*)
ClearAll["Global`*"];
$Assumptions = r > 0 && targetq > 0 &&
  Element[{aa, nn, nnPrime, degreeK, ell, hh, uu, lapU, defU,
    wH, dH, lH, mu, targetf, targetK}, Reals];
checks = {};
check[label_, actual_, expected_] := Module[{ok},
  ok = TrueQ[FullSimplify[actual - expected == 0]];
  AppendTo[checks, <|"identity" -> label, "passed" -> ok|>];
  Print[If[ok, "PASS: ", "FAIL: "], label];
];

(* Homogeneous source h of degree degreeK. Factor r^aa before taking
   derivatives, so all remaining expressions are rational functions. *)
sourceGram = {
  {4 r, 2 degreeK hh, 4 (degreeK - 1) uu},
  {2 degreeK hh, uu, 2 (uu lapU - defU)},
  {4 (degreeK - 1) uu, 2 (uu lapU - defU), 0}
};
gradGBar = {aa hh/r, 1, 0};
normGBar = Expand[gradGBar . sourceGram . gradGBar];
check["radial master: gradient square / R^(2a)", normGBar,
  uu + 4 aa (aa + degreeK) hh^2/r];
lapGBar = lapU + 2 aa (nn + 2 aa + 2 degreeK - 2) hh/r;
gradNormBar = {
  D[normGBar, r] + 2 aa normGBar/r,
  D[normGBar, hh], D[normGBar, uu]
};
defectGBar = normGBar lapGBar - (gradGBar . sourceGram . gradNormBar)/2;
masterExpected = defU + 2 aa (nn - 2 degreeK - 2 aa - 1) hh uu/r +
  4 aa (degreeK + aa) hh^2 lapU/r +
  4 aa (degreeK + aa) (2 aa (nn - 1) + degreeK) hh^3/r^2;
check["radial master: defect / R^(3a)", defectGBar, masterExpected];

(* Clifford chain rule, independently from the radial master identity. *)
pullbackGram = {
  {4 r, 4 ell hh, 4 (2 ell - 2) wH},
  {4 ell hh, 4 r wH, 8 r (wH dH - lH)},
  {4 (2 ell - 2) wH, 8 r (wH dH - lH), 0}
};
pullbackNorm = 4 r wH;
pullbackLaplacian = 4 r dH;
pullbackDefect = pullbackNorm pullbackLaplacian -
  ({0, 1, 0} . pullbackGram .
    {D[pullbackNorm, r], D[pullbackNorm, hh], D[pullbackNorm, wH]})/2;
check["unweighted Clifford pullback: defect", pullbackDefect,
  16 r^2 lH - 8 ell hh wH];

pullbackRules = {degreeK -> 2 ell, uu -> 4 r wH,
  lapU -> 4 r dH, defU -> 16 r^2 lH - 8 ell hh wH};
bb = aa (aa + 2 ell);
bigA = aa (nn - 1) - ell - 2 bb;
weightedNormBar = normGBar /. pullbackRules;
weightedLaplacianBar = lapGBar /. pullbackRules;
weightedDefectBar = defectGBar /. pullbackRules;
check["weighted Clifford pullback: gradient / R^(2a)",
  weightedNormBar, 4 r wH + 4 bb hh^2/r];
(* Specialization consistency for the entered radial Laplacian.
   Independent second-order differentiation: Verify_Appendix_A.wl,
   CheckRadialMaster, assertion "radial Laplacian". *)
(* Specialization consistency for the entered radial Laplacian.
   Independent second-order differentiation: Verify_Appendix_A.wl,
   CheckRadialMaster, assertion "radial Laplacian". *)
check["weighted Clifford pullback: Laplacian / R^a",
  weightedLaplacianBar, 4 r dH + 2 aa (nn + 2 aa + 4 ell - 2) hh/r];
check["weighted Clifford pullback: defect / R^(3a)",
  weightedDefectBar, 16 r^2 lH + 8 bigA hh wH + 16 bb hh^2 dH +
    8 bb (aa (nn - 1) + ell) hh^3/r^2];
check["weighted Clifford pullback: preserved M-ratio / R^(2a)",
  hh weightedLaplacianBar + mu weightedNormBar,
  4 r (hh dH + mu wH) +
    (2 aa (nn + 2 aa + 4 ell - 2) + 4 mu aa (aa + 2 ell)) hh^2/r];

(* Finite difference at fixed (R,P); nnPrime is a second source dimension. *)
check["dimension monotonicity: finite difference / R^(4a)",
  hh ((weightedDefectBar /. nn -> nnPrime) - weightedDefectBar),
  2 aa (nnPrime - nn) hh^2 weightedNormBar/r];

(* Simons quotient on R^targetK + R^targetK; all differentiation is
   performed independently in (targetq,targetf). *)
targetGamma[v_, w_] :=
  4 targetq (D[v, targetq] D[w, targetq] + D[v, targetf] D[w, targetf]) +
  4 targetf (D[v, targetq] D[w, targetf] + D[v, targetf] D[w, targetq]);
targetDelta[v_] := 4 targetq (D[v, {targetq, 2}] + D[v, {targetf, 2}]) +
  8 targetf D[v, targetq, targetf] + 4 targetK D[v, targetq];
targetH = targetq^(3/4) targetf;
targetWH = targetGamma[targetH, targetH];
targetDH = targetDelta[targetH];
targetLH = targetWH targetDH - targetGamma[targetH, targetWH]/2;
check["Simons fractional target: W_H", targetWH,
  4 targetq^(5/2) + (33/4) targetq^(1/2) targetf^2];
check["Simons fractional target: Delta H", targetDH,
  (3/4) (4 targetK + 7) targetq^(-1/4) targetf];
check["Simons fractional target: L(H)", targetLH,
  targetq^(1/4) targetf ((12 targetK - 47) targetq^2 +
    (33/8) (6 targetK + 1) targetf^2)];
check["Simons fractional target: H Delta H", targetH targetDH,
  (3/4) (4 targetK + 7) targetq^(1/2) targetf^2];
dimensionBound = 1 + 4 ell (1 + mu) + 2 aa (1 + mu) + ell/aa;
check["Simons fractional target: a=1, ell=7/2, M=0 dimension bound",
  dimensionBound /. {aa -> 1, ell -> 7/2, mu -> 0}, 41/2];

passed = And @@ (# ["passed"] & /@ checks);
report = <|"engine" -> "Wolfram Language", "version" -> $Version,
  "native_mathematica_executed" -> True, "passed" -> passed,
  "identity_count" -> Length[checks], "checks" -> checks,
  "scope" -> "Differential identities only; no geometric measure theory verification."|>;
outputDirectory = If[StringLength[$InputFileName] > 0,
  DirectoryName[ExpandFileName[$InputFileName]], Directory[]];
Export[FileNameJoin[{outputDirectory, "mathematica_native_results.json"}], report, "RawJSON"];
Print["All identities passed: ", passed];
Exit[If[passed, 0, 1]];
