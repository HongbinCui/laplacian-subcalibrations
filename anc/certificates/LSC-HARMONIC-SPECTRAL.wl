(* Mathematica 12 certificate: LSC-HARMONIC-SPECTRAL.
   Stable manuscript references: eq:harmonic-log-defect,
   eq:harmonic-log-Jacobi, eq:spectral-polynomial-pairings,
   eq:spectral-Q, eq:spectral-threshold, eq:spectral-explicit-certificate.
   Current locations are in the accompanying certificate crosswalk.

   Scope: nine symbolic parameter identities and 52 finite exact checks
   of elementary-symmetric identities in orders 2 through 5. The square
   root calculations assume s > 0. These computations do not replace
   the arbitrary-order inequalities, actual-image arguments, removal of
   singular sets, or the geometric comparison theorem.
   Loading this file returns an association and writes no files. *)

(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`HarmonicSpectral`"];
Begin["`Private`"];

certificateReport = Module[
  {checks = {}, zero, elementary, f, s, tau, hh, lap, q, hs, fq, hq,
   aux, norm, gammaFV, gammaHV, raw, d0, polynomialCoefficient, kappa,
   kk, tt, lapT, normalHessT, tangentialGradTSquared, curvatureSquared,
   jacobiRatio, special, b, c, d, aa, jj, delta, weight, lapG, pairings,
   curvature, polynomial, certificate, z, u, v, w, qStar, numerator,
   positiveDecomposition, order, balancedS, balancedB, balancedC,
   balanced, beta, orderInt, roots, es, ii, index, ell, left, right,
   pairs, rootSymbols, iRoot, jRoot},

  zero[name_String, expression_] := Module[{residual},
    residual = FullSimplify[Factor[Together[expression]],
      Assumptions -> s > 0];
    AppendTo[checks, <|"Name" -> name,
      "Passed" -> TrueQ[residual === 0],
      "Detail" -> If[TrueQ[residual === 0], "Exact residual: 0",
        "Residual: " <> ToString[residual, InputForm]]|>]
  ];
  elementary[values_List, degree_Integer] :=
    If[degree < 0 || degree > Length[values], 0,
      Total[(Times @@ # &) /@ Subsets[values, {degree}]]];

  aux = hh + 2 tau;
  norm = s + f^2 aux;
  gammaFV = 2 q f + 2 f s aux + f^2 fq;
  gammaHV = hs + 2 tau f^2 aux + f^2 hq;
  raw = f norm (lap + aux) -
    (2 f norm (tau + hh) + gammaFV + f gammaHV)/2;
  d0 = s (lap - tau - hh) - q - hs/2;
  polynomialCoefficient = aux lap - hq/2;
  zero["general_log_multiplier_identity",
    raw - f (d0 + f^2 polynomialCoefficient - f fq/2)];
  zero["divisible_two_coefficient_identity",
    (raw /. fq -> kappa f) -
      f (d0 + f^2 (polynomialCoefficient - kappa/2))];

  (* tLog = h + Log[s]/2 and w = Exp[-tLog] on the regular cone. *)
  lapT = lap + kk/(2 s) - tt/(2 s^2);
  normalHessT = tau + q/s - hs/(2 s) - tt/(4 s^2);
  tangentialGradTSquared = hh + hs/s + tt/(4 s^2);
  curvatureSquared = kk/(2 s) - tt/(2 s^2);
  jacobiRatio = -(lapT - normalHessT) + tangentialGradTSquared +
    curvatureSquared;
  zero["Jacobi_boundary_coefficient", jacobiRatio + d0/s];
  special = d0 /. {lap -> kk/(2 s) - tt/(2 s^2), tau -> q/s,
    hh -> tt/(4 s^2), hs -> tt/(2 s)};
  zero["gradient_norm_multiplier_D0", special - (kk/2 - 2 q - tt/s)];

  weight = s^2 + 5 b f^2 + 3 c d f^2/s;
  lapG = f/s^(3/2) ((2 delta + 3) b s - 3 c d);
  pairings = {
    {f, s^(3/2) + 2 b f^2/Sqrt[s]},
    {s, 6 b f Sqrt[s] + 6 c d f/Sqrt[s]},
    {b, 10 c f Sqrt[s] + 8 aa d f/Sqrt[s]},
    {c, 14 aa f Sqrt[s] + 10 jj d f/Sqrt[s]},
    {d, 2 f s^(3/2) + 4 b d f/Sqrt[s]}
  };
  curvature = weight lapG -
    Total[(D[weight, #[[1]]] #[[2]] &) /@ pairings]/2;
  polynomial = (20 delta + 10) b^2 s + (12 delta - 18) b c d -
    56 c s^2 - 82 aa s d - 30 jj d^2;
  certificate = f/(2 s^(3/2))
    (4 s^2 ((delta - 4) b s - 6 c d) + f^2 polynomial);
  zero["spectral_explicit_certificate_all_orders", curvature - certificate];

  qStar = 5 (z u + v)/(z + u) + 3 (z v + w)/(z + u)^2;
  numerator = z (5 u^2 - 8 v) + 5 u^3 - 2 u v - 6 w;
  zero["reciprocal_root_derivative", D[qStar, z] - numerator/(z + u)^3];
  positiveDecomposition = z u^2 + 3 u^3 + 4 z (u^2 - 2 v) +
    u (u^2 - 2 v) + (u^3 - 6 w);
  zero["antitone_numerator_decomposition", numerator - positiveDecomposition];

  balancedS = order;
  balancedB = order (order - 1)/2;
  balancedC = order (order - 1) (order - 2)/6;
  balanced = 2 (delta - 4) balancedB - 12 balancedC/balancedS;
  zero["balanced_threshold_all_orders",
    balanced - order (order - 1) (delta - 6 + 4/order)];
  zero["threshold_normalization",
    ((delta - 6 + 4/order) /. delta -> beta/4) -
      (beta - 24 + 16/order)/4];

  (* Keep exactly the original finite range: orders 2, 3, 4, and 5. *)
  rootSymbols = Table[Unique["spectralRoot"], {5}];
  Do[
    roots = Take[rootSymbols, orderInt];
    es[index_Integer] := elementary[roots, index];
    Do[
      left = Total[(# D[es[ii], #] D[es[index], #] &) /@ roots];
      right = Total[Table[(index - ii + 1 + 2 ell) es[ii - 1 - ell]
        es[index + ell], {ell, 0, ii - 1}]];
      zero["pairing_order_" <> ToString[orderInt] <> "_indices_" <>
        ToString[ii] <> "_" <> ToString[index], left - right],
      {ii, 1, orderInt}, {index, ii, orderInt}];
    Do[
      zero["root_derivation_order_" <> ToString[orderInt] <> "_index_" <>
        ToString[index], Total[(D[es[index], #] &) /@ roots] -
          (orderInt - index + 1) es[index - 1]],
      {index, 1, orderInt}];
    pairs = Total[Flatten[Table[
      If[iRoot != jRoot, roots[[iRoot]]^2 roots[[jRoot]], 0],
      {iRoot, 1, orderInt}, {jRoot, 1, orderInt}]]];
    zero["Newton_expansion_order_" <> ToString[orderInt],
      es[1] es[2] - 3 es[3] - pairs],
    {orderInt, 2, 5}];

  <|"CertificateID" -> "LSC-HARMONIC-SPECTRAL",
    "ExpectedChecks" -> 61, "Checks" -> checks,
    "Scope" -> "Nine symbolic identities; 52 exact finite checks in orders 2 through 5. Analytic and arbitrary-order inequality proofs remain in the manuscript."|>
];

End[];
EndPackage[];
LaplacianCertificates`HarmonicSpectral`Private`certificateReport
