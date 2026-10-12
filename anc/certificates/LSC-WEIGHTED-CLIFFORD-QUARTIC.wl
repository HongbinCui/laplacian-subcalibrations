(* Mathematica 12 certificate: LSC-WEIGHTED-CLIFFORD-QUARTIC.
   Stable manuscript references: eq:weighted-Clifford-table,
   eq:weighted-Clifford-certificate, cor:weighted-Clifford-Hardy,
   eq:homogeneous-strict-stability. Current locations are in the crosswalk.

   The differential table is INPUT. Sign conclusions are conditional on
   R >= 0, W >= 0, K = R^2 - 4 W >= 0, and T = F^2 >= 0.
   The boundary quotient additionally requires R > 0 and W > 0.
   This file does not prove the coordinate table, the actual invariant
   image, regularity, removability, comparison, or the Hardy theorem.
   Loading this file returns an association and writes no files. *)

(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`WeightedCliffordQuartic`"];
Begin["`Private`"];

certificateReport = Module[
  {checks = {}, zero, positiveCoefficients, radius, weight, definingQuartic,
   slack, square, variables, gammaTable, generatorLaplacians, gradient,
   gamma, laplacian, potential, gradientSquare, deltaG, curvature,
   bracket, positiveBracket, lOverF, boundaryBeta, hardyGap},

  zero[name_String, expression_] := Module[{residual},
    residual = Expand[Numerator[Together[expression]]];
    AppendTo[checks, <|"Name" -> name,
      "Passed" -> TrueQ[residual === 0],
      "Detail" -> If[TrueQ[residual === 0], "Exact residual: 0",
        "Residual: " <> ToString[residual, InputForm]]|>]
  ];
  positiveCoefficients[name_String, expression_, indeterminates_List] :=
    Module[{coefficients, passed},
      coefficients = Last /@ CoefficientRules[Expand[expression], indeterminates];
      passed = Length[coefficients] > 0 &&
        And @@ (TrueQ[# > 0] & /@ coefficients);
      AppendTo[checks, <|"Name" -> name, "Passed" -> TrueQ[passed],
        "Detail" -> "Nonzero coefficients: " <>
          ToString[coefficients, InputForm] <>
          "; all supplied indeterminates are assumed nonnegative."|>]
    ];

  variables = {radius, weight, definingQuartic};
  gammaTable = {{4 radius, 8 weight, 8 definingQuartic},
    {8 weight, 4 radius weight, 4 radius definingQuartic},
    {8 definingQuartic, 4 radius definingQuartic, 4 radius weight}};
  generatorLaplacians = {64, 28 radius, 0};
  gradient[expression_] := (D[expression, #] & /@ variables);
  gamma[left_, right_] := Expand[gradient[left] . gammaTable . gradient[right]];
  laplacian[expression_] := Expand[
    Total[MapThread[D[expression, #1] #2 &, {variables, generatorLaplacians}]] +
    Total[Flatten[Table[
      gammaTable[[i, j]] D[expression, variables[[i]], variables[[j]]],
      {i, 1, Length[variables]}, {j, 1, Length[variables]}]]]];

  potential = weight definingQuartic;
  gradientSquare = gamma[potential, potential];
  deltaG = laplacian[potential];
  curvature = Expand[gradientSquare deltaG - gamma[potential, gradientSquare]/2];
  bracket = (3 radius^2 - 2 weight) weight^2 +
    (21 radius^2 - 6 weight) definingQuartic^2;
  zero["weighted_clifford_gradient_square",
    gradientSquare - 4 radius weight (weight^2 + 3 definingQuartic^2)];
  zero["weighted_clifford_Laplacian", deltaG - 36 radius definingQuartic];
  zero["weighted_clifford_curvature_defect", curvature - 16 weight definingQuartic bracket];
  zero["weighted_clifford_degree_eight",
    Total[MapThread[#1 #2 D[potential, #2] &, {{2, 4, 4}, variables}]] -
      8 potential];

  (* Slack K = R^2 - 4 W is assumed nonnegative; its sign is not proved. *)
  positiveBracket = 3 slack weight^2 + 10 weight^3 +
    21 slack square + 78 weight square;
  zero["weighted_clifford_slack_decomposition",
    bracket - (positiveBracket /. {slack -> radius^2 - 4 weight,
      square -> definingQuartic^2})];
  positiveCoefficients["weighted_clifford_positive_bracket_coefficients",
    positiveBracket, {slack, weight, square}];
  zero["weighted_clifford_signed_defect_certificate",
    potential curvature - (16 weight^2 square positiveBracket /.
      {slack -> radius^2 - 4 weight, square -> definingQuartic^2})];
  positiveCoefficients["weighted_clifford_Laplacian_product_coefficients",
    Expand[potential deltaG] /. definingQuartic^2 -> square, {radius, weight, square}];

  (* Cancel F before restricting to F = 0. The original denominator
     W^3 Gamma(F,F) requires R,W > 0 at the regular boundary. *)
  lOverF = Cancel[curvature/definingQuartic];
  boundaryBeta = Cancel[(radius lOverF/(weight^3 gamma[definingQuartic, definingQuartic])) /.
    definingQuartic -> 0];
  zero["weighted_clifford_boundary_beta",
    boundaryBeta - (12 radius^2/weight - 8)];
  zero["weighted_clifford_beta_lower_bound_remainder",
    weight (boundaryBeta - 40) - 12 (radius^2 - 4 weight)];
  positiveCoefficients["weighted_clifford_beta_slack_coefficient", 12 slack, {slack}];

  (* Arithmetic substitution in the manuscript Hardy formula only. *)
  hardyGap = (8 - 31/2)^2 + 40;
  zero["weighted_clifford_Hardy_gap_385_over_4", hardyGap - 385/4];

  <|"CertificateID" -> "LSC-WEIGHTED-CLIFFORD-QUARTIC",
    "ExpectedChecks" -> 12, "Checks" -> checks,
    "Scope" -> "Formal identities conditional on the supplied differential table and slack inequalities; boundary beta >= 40 and Hardy arithmetic 385/4."|>
];

End[];
EndPackage[];
LaplacianCertificates`WeightedCliffordQuartic`Private`certificateReport
