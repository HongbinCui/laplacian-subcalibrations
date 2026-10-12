(* ::Package:: *)
(* LSC-CLIFFORD-PULLBACK -- Clifford Pullback.
   Compatible with Mathematica 12. No external package or file dependencies.
   Get returns checks; only the central runner writes execution records.

*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`CliffordPullback`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`CliffordPullback`Private`*"];

VerificationLog = {};
EvaluationMessages = {};
EvaluationOK = True;
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
RunVerificationGroup[label_String, fun_Symbol] := Module[{status, evaluationResult, messages},
  status = Block[{System`$MessageList = {}},
    evaluationResult = CheckAbort[Catch[fun[]; True, _, ($Failed &)], $Aborted];
    messages = ToString[#, InputForm] & /@ System`$MessageList;
    evaluationResult === True && Length[messages] == 0
  ];
  If[!TrueQ[status],
    EvaluationOK = False;
    AppendTo[EvaluationMessages, <|"Group" -> label, "Messages" -> messages,
      "Status" -> ToString[evaluationResult, InputForm]|>]];
  If[evaluationResult === $Aborted, Abort[]];
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

(* Clifford Pullback identities.
   Definitions only; run through the accompanying chapter harness.
   Exact calculations retain their original mathematical expressions. *)

CheckQuartic[] := Module[
  {r, q, f, n, k, vars, matrix, drift, potential, norm, defect,
   weight, beta, fs, z, sphericalMatrix, sphericalDrift, h, w, numerator},
vars = {r, q, f};
  matrix = {{4 r, 8 q, 8 f}, {8 q, 16 r q, 16 r f},
            {8 f, 16 r f, 16 r q}};
  drift = {2 n, 16 k r, 0};
  potential = r q^(3/4) f;
  norm = r Sqrt[q] (16 r^2 q^2 + 33 r^2 f^2 + 32 q f^2);
  defect = r q^(1/4) (16 (12 k - 47) r^4 q^2 +
    16 (2 n - 41) r^2 q^3 +
    f^2 ((396 k + 66) r^4 + (384 k + 66 n - 681) r^2 q +
      (64 n + 160) q^2));
  CheckEqual[GammaExpr[potential, potential, vars, matrix], norm,
    "Riedler corrected gradient"];
  CheckEqual[CurvatureExpr[potential, vars, matrix, drift], f defect,
    "Riedler corrected defect"];
  CheckEqual[Det[matrix], 1024 r (r^2 - q) (q^2 - f^2),
    "quartic cometric determinant"];
  weight = r q^(3/4);
  beta = r (defect /. f -> 0)/(weight^3 16 r q);
  CheckEqual[beta, 2 n - 41 + (12 k - 47) r^2/q,
    "quartic Jacobi coefficient"];
  sphericalMatrix = {{16 z - 16 fs^2, 16 fs (1 - z)},
                     {16 fs (1 - z), 16 z (1 - z)}};
  sphericalDrift = {-4 (n + 2) fs, 16 k - 4 (n + 2) z};
  h = z^(3/4) fs;
  w = 81 h^2 + GammaExpr[h, h, {fs, z}, sphericalMatrix];
  numerator = w (LaplaceExpr[h, {fs, z}, sphericalMatrix,
       sphericalDrift] + 9 (n - 1) h) -
    GammaExpr[h, w, {fs, z}, sphericalMatrix]/2;
  CheckEqual[w, norm /. {r -> 1, q -> z, f -> fs},
    "quartic spherical gradient"];
  CheckEqual[numerator, (f defect) /. {r -> 1, q -> z, f -> fs},
    "quartic spherical defect"];
];

CheckPullbackExtended[] := Module[
  {t, r, q, f, ff, g, m, n, c, d, b, z, targetVars, targetMatrix,
   targetDrift, h, targetSquare, targetDefect, targetTransfer, pullRules,
   sourceVars, sourceMatrix, sourceDrift, potential, sourceSquare,
   sourceDefect, cartanMatrix, expectedCartanMatrix, cartanDrift,
   dimension, lawsonMatrix, lawsonDrift, hb, hbSquare, coefficientA,
   coefficientB, coefficientC, sourceLawsonMatrix, sourceLawsonDrift,
   boundaryDefect, beta},
(* General harmonic Cartan--Muenzer target and weighted transfer. *)
  targetVars = {t, f};
  targetMatrix = {{4 t, 2 g f}, {2 g f, g^2 t^(g - 1)}};
  targetDrift = {2 m, 0};
  h = t f;
  targetSquare = GammaExpr[h, h, targetVars, targetMatrix];
  targetDefect = CurvatureExpr[h, targetVars, targetMatrix, targetDrift];
  CheckEqual[targetSquare, g^2 t^(g + 1) + 4 (g + 1) t f^2,
    "harmonic CM target gradient"];
  CheckEqual[LaplaceExpr[h, targetVars, targetMatrix, targetDrift],
    2 (m + 2 g) f, "harmonic CM target Laplacian"];
  targetTransfer = 2 t targetDefect - (g + 2) h targetSquare;
  CheckEqual[targetTransfer,
    t^2 f (g^2 t^g (4 m - 2 g^2 - 7 g - 14) +
      4 (g + 1) f^2 (4 m + g - 6)),
    "harmonic CM weighted transfer polynomial"];

  (* Source closure derived from horizontal conformality and Euler. *)
  pullRules = {t -> q, f -> ff};
  sourceVars = {r, q, ff};
  sourceMatrix = {
    {4 r, 8 q, 4 g ff},
    {8 q, 4 r (GammaExpr[t, t, targetVars, targetMatrix] /. pullRules),
      4 r (GammaExpr[t, f, targetVars, targetMatrix] /. pullRules)},
    {4 g ff, 4 r (GammaExpr[f, t, targetVars, targetMatrix] /. pullRules),
      4 r (GammaExpr[f, f, targetVars, targetMatrix] /. pullRules)}};
  sourceDrift = {2 n,
    4 r (LaplaceExpr[t, targetVars, targetMatrix, targetDrift] /. pullRules),
    4 r (LaplaceExpr[f, targetVars, targetMatrix, targetDrift] /. pullRules)};
  CheckEqual[Det[sourceMatrix],
    256 g^2 r (r^2 - q) (q^g - ff^2),
    "harmonic CM source cometric determinant"];
  potential = q ff;
  sourceSquare = GammaExpr[potential, potential, sourceVars, sourceMatrix];
  sourceDefect = CurvatureExpr[potential, sourceVars, sourceMatrix, sourceDrift];
  CheckEqual[sourceSquare, 4 r q (g^2 q^g + 4 (g + 1) ff^2),
    "harmonic CM source gradient"];
  CheckEqual[LaplaceExpr[potential, sourceVars, sourceMatrix, sourceDrift],
    8 (m + 2 g) r ff, "harmonic CM source Laplacian"];
  CheckEqual[sourceDefect,
    16 r^2 (targetDefect /. pullRules) -
      8 (g + 2) (h targetSquare /. pullRules),
    "harmonic CM pullback defect identity"];
  CheckEqual[sourceDefect,
    8 q ff (g^2 q^g ((4 m - 2 g^2 - 6 g - 12) r^2 - (g + 2) q) +
      4 (g + 1) ff^2 ((4 m + 2 g - 4) r^2 - (g + 2) q)),
    "harmonic CM explicit source certificate"];

  (* Cubic Cartan targets and degree-ten source potential. *)
  cartanMatrix = sourceMatrix /. g -> 3;
  expectedCartanMatrix = {{4 r, 8 q, 12 ff},
    {8 q, 16 r q, 24 r ff}, {12 ff, 24 r ff, 36 r q^2}};
  CheckTrue[cartanMatrix === expectedCartanMatrix,
    "Cartan source closure specialization"];
  cartanDrift = {2 n, 8 m r, 0};
  CheckEqual[GammaExpr[ff, ff, sourceVars, cartanMatrix], 36 r q^2,
    "Cartan sextic gradient"];
  CheckEqual[CurvatureExpr[ff, sourceVars, cartanMatrix, cartanDrift],
    -216 ff q (q + 4 r^2), "Cartan sextic defining defect"];
  CheckEqual[CurvatureExpr[potential, sourceVars, cartanMatrix, cartanDrift],
    8 q ff (9 q^3 (4 (m - 12) r^2 - 5 q) +
      16 ff^2 ((4 m + 2) r^2 - 5 q)),
    "Cartan degree-ten source defect"];
  CheckEqual[Det[cartanMatrix], 2304 r (r^2 - q) (q^3 - ff^2),
    "Cartan source determinant"];
  Do[
    CheckTrue[4 dimension - 53 > 0,
      "Cartan first transfer coefficient m=" <> ToString[dimension]];
    CheckTrue[4 dimension - 3 > 0,
      "Cartan second transfer coefficient m=" <> ToString[dimension]],
    {dimension, {14, 26}}];
  CheckEqual[2 r D[potential, r] + 4 q D[potential, q] +
      6 ff D[potential, ff], 10 potential,
    "Cartan source potential homogeneous degree ten"];

  (* General unequal Lawson target and symbolic radial power. *)
  lawsonMatrix = {{4 t, 4 f}, {4 f, 4 (d t + c f)}};
  lawsonDrift = {2 m, 2 c};
  CheckEqual[CurvatureExpr[f, targetVars, lawsonMatrix, lawsonDrift],
    -8 d f, "unequal Lawson target defining defect"];
  hb = t^b f;
  hbSquare = GammaExpr[hb, hb, targetVars, lawsonMatrix];
  coefficientA = -4 b^2 + (2 m - 11) b - 3;
  coefficientC = -2 b^2 + (2 m - 7) b - 1;
  coefficientB = b (b + 2) ((2 m - 3) b + 1);
  CheckEqual[(2 t CurvatureExpr[hb, targetVars, lawsonMatrix, lawsonDrift] -
       (2 b + 2) hb hbSquare)/(8 t^(3 b + 1) f),
    d coefficientA + c coefficientC f/t + coefficientB f^2/t^2,
    "unequal Lawson symbolic radial-power transfer"];
  CheckEqual[2 t CurvatureExpr[h, targetVars, lawsonMatrix, lawsonDrift] -
      4 h GammaExpr[h, h, targetVars, lawsonMatrix],
    16 t^2 f (d (m - 9) t^2 + c (m - 5) t f + 3 (m - 1) f^2),
    "unequal Lawson linear radial-power transfer"];

  (* Target dimensions p=5, q=6; all square completions are exact. *)
  CheckEqual[40 + 6 z + 30 z^2, 30 (z + 1/10)^2 + 397/10,
    "unequal Lawson transfer square completion"];
  sourceLawsonMatrix = {{4 r, 8 q, 8 ff},
    {8 q, 16 r q, 16 r ff},
    {8 ff, 16 r ff, 16 r (20 q + ff)}};
  sourceLawsonDrift = {2 n, 88 r, 8 r};
  CheckEqual[GammaExpr[potential, potential, sourceVars, sourceLawsonMatrix],
    16 r q (20 q^2 + q ff + 3 ff^2),
    "unequal Lawson source gradient"];
  CheckEqual[CurvatureExpr[ff, sourceVars, sourceLawsonMatrix, sourceLawsonDrift],
    -64 ff (40 r^2 + 20 q + ff), "unequal Lawson source defining defect"];
  CheckEqual[CurvatureExpr[potential, sourceVars, sourceLawsonMatrix,
      sourceLawsonDrift],
    128 q ff ((r^2 - q) (60 q^2 + 7 q ff + 33 ff^2) +
      q (40 q^2 + 6 q ff + 30 ff^2)),
    "unequal Lawson explicit source certificate"];
  CheckEqual[60 + 7 z + 33 z^2, 33 (z + 7/66)^2 + 7871/132,
    "unequal Lawson target defect square completion"];
  CheckEqual[20 + z + 3 z^2, 3 (z + 1/6)^2 + 239/12,
    "unequal Lawson target gradient square completion"];
  CheckEqual[2 r D[potential, r] + 4 q D[potential, q] +
      4 ff D[potential, ff], 8 potential,
    "unequal Lawson source potential homogeneous degree eight"];

  (* Cartan and Lawson Jacobi coefficients and Hardy remainders. *)
  boundaryDefect = Cancel[CurvatureExpr[potential, sourceVars, cartanMatrix,
    cartanDrift]/ff];
  beta = Cancel[(r boundaryDefect/
      (q^3 GammaExpr[ff, ff, sourceVars, cartanMatrix])) /. ff -> 0];
  CheckEqual[beta, (8 m - 96) r^2/q - 10,
    "Cartan boundary Jacobi coefficient"];
  Do[
    CheckTrue[8 dimension - 106 > 0,
      "Cartan positive Hardy remainder m=" <> ToString[dimension]],
    {dimension, {14, 26}}];
  boundaryDefect = Cancel[CurvatureExpr[potential, sourceVars,
    sourceLawsonMatrix, {128, 88 r, 8 r}]/ff];
  beta = Cancel[(r boundaryDefect/
      (q^3 GammaExpr[ff, ff, sourceVars, sourceLawsonMatrix])) /. ff -> 0];
  CheckEqual[beta, 24 r^2/q - 8, "unequal Lawson boundary Jacobi coefficient"];
  CheckEqual[(8 - 63/2)^2 + 16, 2273/4,
    "unequal Lawson Hardy bound 2273/4"];
];

CheckQuarticGradientSquareLaplacian[] := Module[
  {f, s, e, n, alpha, k, p, x, r, q, y, z, vars, matrix,
   momentResidual, phi, jphi, expected,
   discriminantRoot, smallRoot, eigenvalue, potential, norm,
   quarticLaplacian, eigenvalues, powerTwo, powerThree, powerFour,
   modelF, modelS, gradientF, gradientS, projection, diagonalHessian,
   offDiagonalTrace, offDiagonalSquare, angularJacobi, radialJacobi,
   certificateBound, i, j, a, aSquared, maxS, t, b, c, testF, testS,
   hessS, residual, cubicMatrix, cubicDrift},
  (* Compute Delta(|grad F|^2) from the closure table itself. *)
  quarticLaplacian = LaplaceExpr[16 r q, {r, q, f},
    {{4 r, 8 q, 8 f}, {8 q, 16 r q, 16 r f},
     {8 f, 16 r f, 16 r q}}, {2 n, 16 k r, 0}];
  CheckEqual[quarticLaplacian, 32 (n + 8) q + 256 k r^2,
    "quartic gradient-square Laplacian from closure"];
];

CheckQuarticStructureAndObstruction[] := Module[
  {rr,e0,e1,e2,r0,r1,r2,nDim,a0,a1,a2,vars,gram,drift,ff,
   sigma,second,third,gamma,lap,minOp,q,f,w,hnorm,gnorm,hfn,
   hq,dcq,linkGrad,linkDrift,phi,curv,u,v,d},
vars={rr,e0,e1,e2};
  gram={{4rr,8e0,8e1,8e2},{8e0,16rr e0,0,0},
        {8e1,0,16rr e1,0},{8e2,0,0,16rr e2}};
  drift={2nDim,8r0 rr,8r1 rr,8r2 rr};
  gamma[x_,y_] := GammaExpr[x,y,vars,gram];
  lap[x_] := LaplaceExpr[x,vars,gram,drift];
  minOp[x_] := CurvatureExpr[x,vars,gram,drift];
  ff=a0 e0+a1 e1+a2 e2; sigma=r0 a0+r1 a1+r2 a2;
  second=a0^2 e0+a1^2 e1+a2^2 e2;
  third=a0^3 e0+a1^3 e1+a2^3 e2;
  CheckEqual[gamma[ff,ff],16rr second,"orthogonal quartic gradient norm"];
  CheckEqual[lap[ff],8rr sigma,"orthogonal quartic Laplacian"];
  CheckEqual[minOp[ff],128rr^2(sigma second-third)-64ff second,
    "orthogonal quartic minimal-level numerator"];
  vars={u,v,q};
  gram={{4u,0,4q},{0,4v,4q},{4q,4q,4(u+v)q}};
  CheckEqual[gamma[(u-v)^2,4q]^2+
    gamma[(u-v)^2,4(u v-q)]^2+gamma[4q,4(u v-q)]^2,
    0,"orthogonality of the three splitting quartics"];
  vars={rr,q,f};
  gram={{4rr,8q,8f},{8q,16rr q,16rr f},{8f,16rr f,16rr q}};
  drift={32,64rr,0};
  w=gamma[f,f];
  hnorm=lap[w]/2-gamma[f,lap[f]];
  gnorm=gamma[w,w]/(4w); hfn=gamma[f,w]/(2w);
  curv=(hnorm-2gnorm+hfn^2)/w;
  CheckEqual[curv/.f->0,24rr/q+14/rr,"Riedler sixteen-dimensional curvature"];
  hq=(gamma[f,gamma[f,q]]-gamma[q,w]/2)/w;
  dcq=lap[q]-hq;
  CheckEqual[dcq/.f->0,56rr+4q/rr,"Riedler induced Laplacian of Q"];
  linkGrad=(gamma[q,q]-gamma[q,f]^2/w)/. {rr->1,f->0};
  CheckEqual[linkGrad-16q^2,16q(1-q),"Riedler link gradient norm"];
  linkDrift=(dcq/.{rr->1,f->0})-4(4+15-2)q;
  CheckEqual[linkDrift,56-64q,"Riedler link Laplacian"];
  phi=1/q;
  CheckEqual[16q(1-q)D[phi,{q,2}]+(56-64q)D[phi,q]+(24/q+14)phi,
    46phi,"Riedler angular eigenfunction"];
  CheckEqual[(13/2)^2-46,-15/4,"Riedler negative Hardy quotient"];
];



Block[{$IterationLimit = Infinity, $Assumptions = True},
  RunVerificationGroup["CheckQuartic", CheckQuartic];
  RunVerificationGroup["CheckPullbackExtended", CheckPullbackExtended];
  RunVerificationGroup["CheckQuarticGradientSquareLaplacian", CheckQuarticGradientSquareLaplacian];
  RunVerificationGroup["CheckQuarticStructureAndObstruction", CheckQuarticStructureAndObstruction];
];

CertificateResult = <|"CertificateID" -> "LSC-CLIFFORD-PULLBACK",
  "ExpectedChecks" -> 50,
  "Checks" -> (<|"Name" -> #[[1]], "Passed" -> TrueQ[#[[2]]],
      "Detail" -> If[StringQ[#[[3]]], #[[3]], ToString[#[[3]], InputForm]]|> &) /@ VerificationLog,
  "EvaluationOK" -> TrueQ[EvaluationOK],
  "EvaluationMessages" -> EvaluationMessages|>;
End[];
EndPackage[];
LaplacianCertificates`CliffordPullback`Private`CertificateResult
