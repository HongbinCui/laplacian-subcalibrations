(* Mathematica 12 certificate: LSC-REVIEW-STRENGTHENINGS.
   Stable manuscript references: eq:nr-five-radius-recovery,
   eq:nr-cartan-reduction, eq:nr-cartan-reconstruction,
   cor:nr-cartan-dimensions, cor:quartic-nine-weights, bh:potentials,
   eq:quartic-nine-certificate, eq:quartic-nine-discriminant,
   eq:spectral-radial-Newton, eq:spectral-Hardy-gap,
   eq:octonionic-product-Hardy-gap, eq:pfaffian-Hardy-gap.
   Current locations are in the accompanying certificate crosswalk.

   Intrinsic cubic and quartic differential tables are INPUT. Checks do
   not prove those tables, Cartan classification, quartic realization,
   removability, Newton inequalities, or geometric comparison. The
   excluded low-weight cases use the exact original 3 by 7 enumeration.
   Loading this file returns an association and writes no files. *)

(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`ReviewStrengthenings`"];
Begin["`Private`"];

certificateReport = Module[
  {checks = {}, zero, positiveCoefficients, r, s, a, p, h, m,
   variables, matrix, deltas, gamma, laplacian, f, cartan, radius,
   gradientSquare, k, r0, c, t, z, alpha, u, aa, bb, p10, coeffA,
   coeffB, coeffC, p9, discriminant, lower, discriminant8, cliffordDelta,
   cases, pair, unusedWeight, sorted, index, ambientDimension, degree},

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

  variables = {r, s, a, p, h};
  matrix = {
    {4 r, 0, 6 a, 2 p, 4 h},
    {0, 4 s, 0, 4 p, 4 h},
    {6 a, 0, 9 r^2, 18 h - r s, r p + s a/3},
    {2 p, 4 p, 18 h - r s, s^2/9 + 4 h,
      r p/3 + 2 s p/9 + s a/27},
    {4 h, 4 h, r p + s a/3, r p/3 + 2 s p/9 + s a/27,
      (r + s) h/3 + p^2/9 + a p/27}
  };
  deltas = {2 (m + 2), 2 m, 0, 0, (m r + (m + 2) s)/9};
  gamma[left_, right_] := Expand[
    (D[left, #] & /@ variables) . matrix . (D[right, #] & /@ variables)];
  laplacian[expression_] := Expand[
    Total[MapThread[#1 D[expression, #2] &, {deltas, variables}]] +
    Total[Flatten[Table[
      matrix[[i, j]] D[expression, variables[[i]], variables[[j]]],
      {i, 1, Length[variables]}, {j, 1, Length[variables]}]]]];

  f = (a + 3 p)/6;
  cartan = a - 9 p;
  radius = r + s;
  gradientSquare = gamma[f, f];
  zero["intrinsic_gradient_square", gradientSquare - (r - s/3)^2/4 - 4 h];
  zero["intrinsic_laplacian_gradient_square",
    laplacian[gradientSquare] - ((10 m + 36) r + (2 m + 4) s)/9];
  zero["intrinsic_radius_recovery",
    (9 laplacian[gradientSquare] - (2 m + 4) radius)/(8 (m + 4)) - r];
  zero["Cartan_reduction_from_R_F_algebra", cartan - 6 (gamma[r, f] - 5 f)];
  zero["Cartan_reduction_eikonal_identity", gamma[cartan, cartan] - 9 radius^2];
  zero["Cartan_reduction_harmonicity", laplacian[cartan]];
  zero["Cartan_reconstruction", (4 a/6 - cartan/6)/3 - f];

  aa = 1 - 2 c;
  bb = c (1 - c);
  p10 = bb (k - 8 + (2 k - 23) t) +
    aa (k - 4 + (2 k - 14) t) z +
    (3 k + (15 k - 33) t + (18 k + 27) t^2) z^2 +
    2 (r0 - 1) t (bb + aa z + 3 (1 + 3 t) z^2);
  coeffA = 27 + 144 t + 315 t^2;
  coeffB = 5 + 18 t;
  coeffC = 1 + 9 t;
  p9 = bb coeffC + aa coeffB z + coeffA z^2;
  zero["quartic_h_nine_r0_eight", (p10 /. {k -> 9, r0 -> 8}) - p9];
  zero["quartic_h_nine_weight_remainder",
    (p10 /. k -> 9) - p9 -
      2 (r0 - 8) t (bb + aa z + 3 (1 + 3 t) z^2)];
  discriminant = (1 - alpha) coeffA coeffC - alpha coeffB^2;
  lower = (23 + 4788 t + 30564 t^2 + 68040 t^3)/49;
  zero["quartic_h_nine_negative_discriminant",
    -Discriminant[p9, z] - (discriminant /. alpha -> aa^2)];
  zero["quartic_h_nine_discriminant_bound",
    (discriminant /. alpha -> 25/49) - lower];
  zero["quartic_h_nine_parameter_remainder",
    discriminant - lower - (25/49 - alpha) (coeffA coeffC + coeffB^2)];
  positiveCoefficients["quartic_h_nine_positive_discriminant_bound", lower, {t}];
  positiveCoefficients["quartic_h_nine_positive_leading_coefficient", coeffA, {t}];
  discriminant8 = 12 (k - 1) (k - 9) bb - (k - 5)^2 aa^2;
  zero["degree_eight_h_eleven", (discriminant8 /. {k -> 11, c -> 1/9}) - 52/27];
  zero["degree_eight_h_ten_active_weights_at_least_three",
    (discriminant8 /. {k -> 10, c -> 1/4}) - 14];

  cliffordDelta[j_Integer] :=
    {1, 2, 4, 4, 8, 8, 8, 8}[[Mod[j - 1, 8] + 1]] 16^Quotient[j - 1, 8];
  cases = Flatten[Table[
    sorted = Sort[Prepend[pair, unusedWeight]];
    index = cliffordDelta[sorted[[1]] + sorted[[2]] - 1];
    <|"Weights" -> Prepend[pair, unusedWeight], "Ell" -> unusedWeight + 8,
      "Delta" -> index, "Excluded" -> (Mod[unusedWeight + 8, index] != 0)|>,
    {pair, {{2, 7}, {3, 6}, {4, 5}}}, {unusedWeight, 1, 7}], 1];
  AppendTo[checks, <|"Name" -> "quartic_h_nine_low_unused_weights_excluded",
    "Passed" -> TrueQ[Length[cases] == 21 && And @@ (#["Excluded"] & /@ cases)],
    "Detail" -> "Exact 21-case integer enumeration: " <> ToString[cases, InputForm]|>];
  positiveCoefficients["quartic_h_nine_Hardy_positive",
    (2 (8 + u + 9 - 1) - 17)^2/4, {u}];

  zero["spectral_necessary_dimension_arithmetic",
    (2 degree - 2) (ambientDimension + 2 degree - 4) -
      (24 - 16/degree) degree (degree - 1)/2 -
      2 (degree - 1) (ambientDimension - 4 degree)];
  zero["spectral_Hardy_homogeneity",
    (2 degree - 1 - (ambientDimension - 1)/2)^2 -
      (ambientDimension - 4 degree + 1)^2/4];
  zero["octonionic_Hardy_coefficient",
    ((ambientDimension - 4 degree + 1)^2/4 /. ambientDimension -> 8 degree) -
      (4 degree + 1)^2/4];
  zero["Pfaffian_Hardy_coefficient",
    ((ambientDimension - 4 degree + 1)^2/4 /.
      ambientDimension -> 2 degree (2 degree - 1)) -
      (4 degree^2 - 6 degree + 1)^2/4];
  positiveCoefficients["spectral_Hardy_positive_at_necessary_dimension",
    (ambientDimension - 4 degree + 1)^2/4 /.
      ambientDimension -> 4 degree + u, {u}];

  <|"CertificateID" -> "LSC-REVIEW-STRENGTHENINGS",
    "ExpectedChecks" -> 23, "Checks" -> checks,
    "Scope" -> "Conditional intrinsic cubic identities; quartic discriminant and endpoint identities; exact 21-case exclusions; spectral and Hardy arithmetic. Classification and analytic proofs remain in the manuscript."|>
];

End[];
EndPackage[];
LaplacianCertificates`ReviewStrengthenings`Private`certificateReport
