(* ::Package:: *)
(*
  LSC-INTRINSIC-FIVE-GENERATORS
  Mathematica 12 compatible; exact arithmetic throughout.

  Scope: ordinary Euclidean differentiation of the explicit real and complex
  Hermitian coordinate models (ambient real dimensions 8 and 14).  The table
  being checked is NOT used to compute any derivative.  The nonzero rank-five
  minor is an exact algebraic-independence certificate, not numerical sampling.
  This file does not prove the universal structural theorem or the quaternionic
  model, nor does it certify area minimization or a singular-set argument.

  Manuscript references (semantic labels are authoritative):
    thm:nonradial-five-generators
    eq:nr-five-invariants
    eq:nr-five-gamma
    eq:nr-five-laplacians
    eq:nr-five-gradient-square
    eq:nr-five-radius-recovery
  See the current certificate crosswalk for compiled section/equation numbers.

  All scalar coordinate variables are real.  Z is traceless Hermitian, and the
  conjugates below are written explicitly in these real coordinates.  With
  r=Tr[Z.Z]/2, s=u* u, a=(3 Sqrt[3]/2) Det[Z],
  p=(u* Z u)/(2 Sqrt[3]), h=(u* Z.Z u)/12,
  the normalization is F=(a+3 p)/6 and q=(3/2)(r+7 s/27).

  The 31 named checks per model are retained.  The two radius
  recovery checks now also verify the current manuscript's Delta[S] formula;
  the older q-based equalities remain checked in the same records.

  Get[this file] returns the result Association.  It performs no file writes,
  uses no external source or saved output, and makes no native-run claim.
*)

(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

Begin["LaplacianCertificates`IntrinsicFiveGenerators`"];

Module[{verifyModel, allChecks, result},
  verifyModel[divisionDimension_Integer] := Module[
    {x, z, b12, b13, b23, c12, c13, c23, u1, u2, u3, v1, v2, v3,
     realOff, imaginaryOff, off, offBar, zMatrix, uReal, uImaginary,
     uVector, uBar, bCoordinates, uCoordinates, coordinates,
     r, s, a, p, h, functions, functionNames, checks = {},
     check, grad, gamma, laplacian, expectedGamma, expectedDelta,
     entry, j, m, radius, cubic, gradientSquare, weight, hsiang,
     radialDerivative, deltaGradientSquare, recoveredRadius,
     point, columns, jacobian, minor, prefix},

    prefix = "d=" <> ToString[divisionDimension] <> ": ";
    realOff = {b12, b13, b23};
    imaginaryOff = If[divisionDimension == 2, {c12, c13, c23}, {}];
    off = realOff + If[divisionDimension == 2,
      I imaginaryOff, ConstantArray[0, 3]];
    offBar = realOff - If[divisionDimension == 2,
      I imaginaryOff, ConstantArray[0, 3]];
    zMatrix = {
      {x + z/Sqrt[3], off[[1]], off[[2]]},
      {offBar[[1]], -x + z/Sqrt[3], off[[3]]},
      {offBar[[2]], offBar[[3]], -2 z/Sqrt[3]}
    };
    uReal = {u1, u2, u3};
    uImaginary = If[divisionDimension == 2, {v1, v2, v3}, {}];
    uVector = uReal + If[divisionDimension == 2,
      I uImaginary, ConstantArray[0, 3]];
    uBar = uReal - If[divisionDimension == 2,
      I uImaginary, ConstantArray[0, 3]];
    bCoordinates = Join[{x, z}, realOff, imaginaryOff];
    uCoordinates = Join[uReal, uImaginary];
    coordinates = Join[bCoordinates, uCoordinates];

    r = Expand[Tr[zMatrix.zMatrix]/2];
    s = Expand[uBar.uVector];
    a = Expand[3 Sqrt[3] Det[zMatrix]/2];
    p = Expand[(uBar.zMatrix.uVector)/(2 Sqrt[3])];
    h = Expand[(uBar.zMatrix.zMatrix.uVector)/12];
    functions = {r, s, a, p, h};
    functionNames = {"r", "s", "a", "p", "h"};

    check[name_String, lhs_, rhs_, description_String] := Module[
      {residual, passed},
      residual = Flatten[{Expand[lhs - rhs]}];
      If[!And @@ (SameQ[#, 0] & /@ residual),
        residual = FullSimplify[residual]];
      passed = TrueQ[And @@ (SameQ[#, 0] & /@ residual)];
      AppendTo[checks, <|
        "Name" -> prefix <> name,
        "Passed" -> passed,
        "Detail" -> If[passed, description,
          description <> " Residual: " <> ToString[residual, InputForm]]
      |>];
      passed
    ];
    grad[polynomial_] := D[polynomial, #] & /@ coordinates;
    gamma[first_, second_] := Expand[grad[first].grad[second]];
    laplacian[polynomial_] := Expand[
      Total[(D[polynomial, {#, 2}] &) /@ coordinates]];

    check["orthonormal_B_coordinates", r, Total[bCoordinates^2],
      "Tr[Z.Z]/2 equals the sum of squares of all B coordinates."];
    check["orthonormal_U_coordinates", s, Total[uCoordinates^2],
      "u* u equals the sum of squares of all U coordinates."];

    (* Every unordered pair is differentiated independently in real coordinates. *)
    expectedGamma = {
      {1, 1, 4 r}, {1, 2, 0}, {1, 3, 6 a}, {1, 4, 2 p}, {1, 5, 4 h},
      {2, 2, 4 s}, {2, 3, 0}, {2, 4, 4 p}, {2, 5, 4 h},
      {3, 3, 9 r^2}, {3, 4, 18 h - r s}, {3, 5, r p + s a/3},
      {4, 4, s^2/9 + 4 h}, {4, 5, r p/3 + 2 s p/9 + s a/27},
      {5, 5, (r + s) h/3 + p^2/9 + a p/27}
    };
    Do[
      check["Gamma(" <> functionNames[[entry[[1]]]] <> "," <>
          functionNames[[entry[[2]]]] <> ")",
        gamma[functions[[entry[[1]]]], functions[[entry[[2]]]]], entry[[3]],
        "Exact polynomial equality from ordinary Euclidean first derivatives."],
      {entry, expectedGamma}
    ];

    m = Length[uCoordinates];
    expectedDelta = {2 (m + 2), 2 m, 0, 0, (m r + (m + 2) s)/9};
    Do[
      check["Delta(" <> functionNames[[j]] <> ")",
        laplacian[functions[[j]]], expectedDelta[[j]],
        "Exact polynomial equality from ordinary Euclidean second derivatives."],
      {j, 1, 5}
    ];

    radius = r + s;
    cubic = Expand[(a + 3 p)/6];
    gradientSquare = gamma[cubic, cubic];
    weight = (3/2) (r + 7 s/27);
    hsiang = Expand[gradientSquare laplacian[cubic] -
      gamma[cubic, gradientSquare]/2];
    check["F_is_harmonic", laplacian[cubic], 0,
      "The normalized cubic has zero Euclidean Laplacian."];
    check["normalized_Hsiang_equation", hsiang, -weight cubic,
      "The exact normalized Hsiang identity is L[F]=-q F."];
    radialDerivative = gamma[r, cubic];
    check["D_equals_a_plus_p", radialDerivative, a + p,
      "D=Gamma(r,F)=a+p by independent differentiation."];
    deltaGradientSquare = laplacian[gradientSquare];
    recoveredRadius = (9 deltaGradientSquare - (2 m + 4) radius)/(8 (m + 4));
    check["recover_r",
      {(18 weight - 7 radius)/20, recoveredRadius}, {r, r},
      "Both the q formula and Delta[S] formula recover r."];
    check["recover_s",
      {(27 radius - 18 weight)/20, radius - recoveredRadius}, {s, s},
      "Both the q formula and Delta[S] formula recover s."];
    check["recover_p", 3 cubic - radialDerivative/2, p,
      "The cubic and D recover p=3 F-D/2."];
    check["recover_a", 3 radialDerivative/2 - 3 cubic, a,
      "The cubic and D recover a=3 D/2-3 F."];
    check["recover_h", (gamma[a, p] + r s)/18, h,
      "Independent differentiation recovers h=(Gamma(a,p)+r s)/18."];

    (* One nonzero exact minor proves generic rank five.  This is the same
       real point and the same five coordinate directions in both models. *)
    point = Join[{x -> 1}, Thread[uReal -> ConstantArray[1, 3]],
      Thread[Complement[coordinates, Join[{x}, uReal]] ->
        ConstantArray[0, Length[coordinates] - 4]]];
    columns = Join[{x, z}, uReal];
    jacobian = Table[D[functions[[j]], columns[[entry]]],
      {j, {1, 3, 2, 4, 5}}, {entry, 1, 5}];
    minor = FullSimplify[Det[jacobian /. point]];
    check["rank_five_exact_minor", minor, 4 Sqrt[3]/3,
      "At Z=diag(1,-1,0), u=(1,1,1), the rows (r,a,s,p,h) and " <>
      "columns (x,z,u1,u2,u3) give the nonzero exact minor 4 Sqrt[3]/3."];
    checks
  ];

  allChecks = Join[verifyModel[1], verifyModel[2]];
  result = <|
    "CertificateID" -> "LSC-INTRINSIC-FIVE-GENERATORS",
    "ExpectedChecks" -> 62,
    "Checks" -> allChecks
  |>;
  End[];
  result
]
