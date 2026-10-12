(* ::Package:: *)
(*
  LSC-JORDAN-PENCIL -- exact Jordan-pencil certificate, Mathematica 12.

  Current manuscript: nonradial Jordan pencils; stable references
    eq:nonradial-pencil,
    eq:nonradial-closure,
    eq:nonradial-pencil-gradient,
    eq:nonradial-pencil-minimality.
  Stable labels, rather than page numbers, identify the statements.

  This unit contains 65 exact checks. Only integer/rational
  arithmetic and symbolic square roots occur.  No random sampling occurs.

  Scope: the formal calculation assumes the displayed invariant closure
  for d = 1, 2, 4.  Its full closure/Laplacian table is independently
  checked in real coordinates for d = 1.  The diagonal real witnesses
  embed unchanged in R, C, and H.  Their nonzero mean-curvature numerator
  excludes every real parameter other than 0, 1, -1/3.  The equivalence
  between this numerator's vanishing at regular zeros and geometric
  minimality is used, not proved by this program.  The identification
  N = A + B, full complex/quaternionic coordinate closure, irreducibility,
  stability, and minimization are outside the certificate's scope.
  No independent rank test is included.

  Get[thisFile] returns an Association.  This file never writes files,
  clears Global`, exits the kernel, or presents a previous run as current.
*)

(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

Begin["LaplacianCertificates`JordanPencil`"];

certificateResult = Module[
  {checks = {}, zero, condition, grad, gamma, lap, rho, nu, aa, bb, ww,
   tt, dim, variables, table, laps, ff, ss, lf, actualS, actualL,
   expectedL, weights, pair, value, weight, xx, yy, za, zb, zc, u1, u2, u3,
   coords, zMatrix, uVector, inverseMetric, generators, substitute,
   actualF, coordinateS, actualLap, gradientVector, coordinateL,
   tPositive, tNegative, branches, branch, label, tb, ub, nub, wb,
   rb, normUSquared, expectedWitnessL, assumption, witness,
   coordinateWitness, z0, derivativeU, witnesses = {}, roots,
   commonZero, i, j, k, category, scope},

  zero[name_String, expression_, group_String, assumptions_: True] :=
    Module[{residual, passed},
      residual = FullSimplify[
        Expand[Numerator[Together[expression]]], assumptions];
      passed = TrueQ[residual === 0];
      AppendTo[checks, <|
        "Name" -> name, "Passed" -> passed,
        "Detail" -> If[passed, group <> ": exact zero",
          group <> ": residual " <> ToString[residual, InputForm]]
      |>];
      passed
    ];
  condition[name_String, expression_, group_String] :=
    Module[{passed = TrueQ[expression]},
      AppendTo[checks, <|
        "Name" -> name, "Passed" -> passed,
        "Detail" -> If[passed, group <> ": exact condition holds",
          group <> ": condition failed or remained unevaluated"]
      |>];
      passed
    ];
  grad[expression_, vars_List] := (D[expression, #] & /@ vars);
  gamma[first_, second_, vars_List, inverse_] :=
    grad[first, vars] . inverse . grad[second, vars];
  lap[expression_, vars_List, inverse_, generatorLaps_List] :=
    Sum[D[expression, vars[[i]]] generatorLaps[[i]],
      {i, Length[vars]}] +
    Sum[inverse[[i, j]] D[expression, vars[[i]], vars[[j]]],
      {i, Length[vars]}, {j, Length[vars]}];

  (* Seventeen formal polynomial/coefficient/special-parameter checks. *)
  variables = {rho, nu, aa, bb, ww};
  table = {
    {4 rho, 0, 6 aa, 2 bb, 4 ww},
    {0, 4 nu, 0, 4 bb, 4 ww},
    {6 aa, 0, (4/3) rho^2, (4/3) rho nu - ww/2,
      16 nu aa - (16/3) rho bb},
    {2 bb, 4 bb, (4/3) rho nu - ww/2, (4/3) nu^2 + ww,
      16 rho bb - 16 nu aa + (32/3) nu bb},
    {4 ww, 4 ww, 16 nu aa - (16/3) rho bb,
      16 rho bb - 16 nu aa + (32/3) nu bb,
      16 (rho + nu) ww + (64/3) bb^2 - 64 aa bb}
  };
  laps = {6 dim + 4, 6 dim, 0, 0,
    16 dim rho + (16/3) (3 dim + 2) nu};
  ff = aa + tt bb;
  ss = (4/3) (rho + tt nu)^2 + tt (tt - 1) ww;
  lf = lap[ff, variables, table, laps];
  actualS = gamma[ff, ff, variables, table];
  actualL = Expand[actualS lf -
    gamma[ff, actualS, variables, table]/2];
  expectedL = -8 (rho + tt^2 (2 - tt) nu) ff -
    (8/3) tt (tt - 1) (3 tt + 1) (rho + tt nu) bb;
  zero["Delta(F_t)=0 for formal d", lf, "formal closure"];
  zero["S_t full polynomial", actualS - ss, "formal closure"];
  zero["L(F_t) full polynomial", actualL - expectedL, "formal closure"];
  Do[zero["S_t coefficient of t^" <> ToString[k],
    Coefficient[Expand[actualS - ss], tt, k], "coefficient checks"],
    {k, 0, 2}];
  Do[zero["L(F_t) coefficient of t^" <> ToString[k],
    Coefficient[Expand[actualL - expectedL], tt, k], "coefficient checks"],
    {k, 0, 4}];
  zero["L(F_t) displayed t^4 cancellation",
    Coefficient[Expand[expectedL], tt, 4], "coefficient checks"];
  weights = {{0, 8 rho}, {1, 8 (rho + nu)},
    {-1/3, 8 (rho + (7/27) nu)}};
  Do[
    {value, weight} = pair;
    zero["minimal equation at t=" <> ToString[value, InputForm],
      (actualL /. tt -> value) + weight (ff /. tt -> value),
      "special parameters"],
    {pair, weights}];
  zero["existing nonradial S at t=-1/3",
    (ss /. tt -> -1/3) - (4/3) (rho - nu/3)^2 - (4/9) ww,
    "special parameters"];
  zero["existing nonradial Gamma(F,S)=2qF",
    (gamma[ff, ss, variables, table] /. tt -> -1/3) -
      16 (rho + (7/27) nu) (aa - bb/3), "special parameters"];

  (* Twenty-three independent d=1 real-coordinate checks. *)
  coords = {xx, yy, za, zb, zc, u1, u2, u3};
  zMatrix = {{xx, zc, zb}, {zc, yy, za}, {zb, za, -xx - yy}};
  uVector = {u1, u2, u3};
  (* Metric: (1/2) tr(dZ^2) + |du|^2. *)
  inverseMetric = IdentityMatrix[8];
  inverseMetric[[1 ;; 2, 1 ;; 2]] = {{4/3, -2/3}, {-2/3, 4/3}};
  generators = Expand /@ {
    Tr[zMatrix . zMatrix]/2, uVector . uVector, Det[zMatrix],
    -uVector . zMatrix . uVector,
    4 uVector . zMatrix . zMatrix . uVector};
  substitute = Thread[variables -> generators];
  Do[
    zero["real-coordinate Delta(" <>
      {"rho", "nu", "A", "B", "W"}[[i]] <> ")",
      lap[generators[[i]], coords, inverseMetric, ConstantArray[0, 8]] -
        (laps[[i]] /. dim -> 1 /. substitute),
      "independent d=1 coordinates"];
    Do[
      zero["real-coordinate Gamma(" <>
        {"rho", "nu", "A", "B", "W"}[[i]] <> "," <>
        {"rho", "nu", "A", "B", "W"}[[j]] <> ")",
        gamma[generators[[i]], generators[[j]], coords, inverseMetric] -
          (table[[i, j]] /. substitute),
        "independent d=1 coordinates"],
      {j, i, Length[variables]}],
    {i, Length[variables]}];
  actualF = generators[[3]] + tt generators[[4]];
  coordinateS = Expand[gamma[actualF, actualF, coords, inverseMetric]];
  actualLap = lap[actualF, coords, inverseMetric, ConstantArray[0, 8]];
  gradientVector = inverseMetric . grad[actualF, coords];
  (* Independent Hessian contraction fixes the convention for L. *)
  coordinateL = Expand[coordinateS actualLap -
    gradientVector .
      Table[D[actualF, coords[[i]], coords[[j]]], {i, 8}, {j, 8}] .
      gradientVector];
  zero["real-coordinate Delta(F_t)", actualLap,
    "independent d=1 coordinates"];
  zero["real-coordinate S_t", coordinateS - (ss /. substitute),
    "independent d=1 coordinates"];
  zero["real-coordinate Hessian definition of L(F_t)",
    coordinateL - (actualL /. substitute), "independent d=1 coordinates"];

  (* Twenty-five exact regular-zero/obstruction checks. *)
  branches = {
    {"t>0", tPositive, {0, 0, 1/Sqrt[tPositive]},
      1/tPositive, 16/tPositive, 4, 16 tPositive,
      -(64/3) (tPositive - 1) (3 tPositive + 1), tPositive > 0},
    {"t<0", tNegative, {Sqrt[-2/tNegative], 0, 0},
      -2/tNegative, -8/tNegative, 1, -8 tNegative,
      -(16/3) (tNegative - 1) (3 tNegative + 1), tNegative < 0}
  };
  z0 = DiagonalMatrix[{1, 1, -2}];
  category = "regular zeros";
  Do[
    {label, tb, ub, nub, wb, rb, normUSquared, expectedWitnessL,
      assumption} = branch;
    witness = {rho -> 3, nu -> nub, aa -> -2, bb -> 2/tb,
      ww -> wb, tt -> tb};
    coordinateWitness = Join[
      Thread[coords -> Join[{1, 1, 0, 0, 0}, ub]], {tt -> tb}];
    zero[label <> " witness nu", ub . ub - nub, category, assumption];
    zero[label <> " witness A", Det[z0] + 2, category, assumption];
    zero[label <> " witness B", -ub . z0 . ub - 2/tb,
      category, assumption];
    zero[label <> " witness W", 4 ub . z0 . z0 . ub - wb,
      category, assumption];
    zero[label <> " witness F_t=0", ff /. witness, category, assumption];
    zero[label <> " witness rho+t*nu", ((rho + tt nu) /. witness) - rb,
      category, assumption];
    derivativeU = -2 tb z0 . ub;
    zero[label <> " witness |D_u F_t|^2",
      derivativeU . derivativeU - normUSquared, category, assumption];
    condition[label <> " witness is regular",
      FullSimplify[normUSquared > 0, assumption], category];
    zero[label <> " L(F_t) at witness",
      (actualL /. witness) - expectedWitnessL, category, assumption];
    zero[label <> " direct coordinate F_t=0", actualF /. coordinateWitness,
      category, assumption];
    zero[label <> " direct coordinate L(F_t)",
      (coordinateL /. coordinateWitness) - expectedWitnessL,
      category, assumption];
    AppendTo[witnesses, <|
      "ParameterBranch" -> label, "Z" -> z0, "u" -> ub,
      "rho" -> 3, "nu" -> nub, "A" -> -2, "B" -> 2/tb,
      "W" -> wb, "rhoPlusTNu" -> rb, "F_t" -> 0,
      "SquaredUGradientNorm" -> normUSquared,
      "LAtRegularZero" -> expectedWitnessL
    |>],
    {branch, branches}];
  roots = Sort[tt /. Solve[(tt - 1) (3 tt + 1) == 0, tt]];
  condition["exceptional nonzero roots are exactly -1/3 and 1",
    roots === {-1/3, 1}, category];
  commonZero = {xx -> 1, yy -> -1, za -> 0, zb -> 0, zc -> 0,
    u1 -> 0, u2 -> 0, u3 -> 0};
  zero["all parameters have common zero diag(1,-1,0),u=0",
    actualF /. commonZero, category];
  zero["common zero gradient squared norm equals 4/3",
    (coordinateS /. commonZero) - 4/3, category];

  scope = {
    "Formal differential closure is assumed for d=1,2,4.",
    "The complete table is independently verified in real coordinates for d=1.",
    "The diagonal real witnesses embed unchanged into R,C,H.",
    "The regular-zero obstruction excludes all real t outside {0,1,-1/3}.",
    "The geometric minimality equivalence is assumed, not formally proved.",
    "Full complex/quaternionic coordinate closure, irreducibility, rank, stability and minimization are outside this certificate."
  };
  <|
    "CertificateID" -> "LSC-JORDAN-PENCIL",
    "ExpectedChecks" -> 65, "Checks" -> checks,
    "Scope" -> scope, "RegularZeroWitnesses" -> witnesses,
    "References" -> {"eq:nonradial-pencil", "eq:nonradial-closure",
      "eq:nonradial-pencil-gradient", "eq:nonradial-pencil-minimality"}
  |>
];

End[];
LaplacianCertificates`JordanPencil`certificateResult
