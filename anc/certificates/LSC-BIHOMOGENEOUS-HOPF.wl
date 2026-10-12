(* ::Package:: *)
(* LSC-BIHOMOGENEOUS-HOPF -- Bihomogeneous Hopf.
   Compatible with Mathematica 12. No external package or file dependencies.
   Get returns checks; only the central runner writes execution records.

   Four formerly separate Hardy residual calculations are now checked and logged with labels. No external script is needed.
*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`BihomogeneousHopf`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`BihomogeneousHopf`Private`*"];

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

(* Bihomogeneous Hopf identities.
   Definitions only; run through the accompanying chapter harness.
   Exact calculations retain their original mathematical expressions. *)

CheckBihomogeneous[] := Module[
  {u,v,q,rr,tt,ell,mm,cc,aa,bb,ff,gg,ww,vars,gram,lapVector,
   gamma,lap,minOp,assertZero,minimalRule,defect8,p10,ss,zz,defect10,
   p8,p8bal,disc8,cstar,disc10,wf,hessSquared,normalDerivativeSquared,
   hessNormal,curvature,hessTNormal,linkLap,tvar,linkPotential,
   data,degree,kappa,level,pvar,sourceF},
assertZero[name_String, expression_] := CheckEqual[expression,0,name];
vars = {u, v, q};
rr = u + v; tt = u v;
gram = {{4 u, 0, 4 q}, {0, 4 v, 4 q}, {4 q, 4 q, 4 (u + v) q}};
lapVector = {2 ell, 2 ell, 2 mm (u + v)};
gamma[arg1_, arg2_] := Total[Flatten[Table[
  D[arg1, vars[[i]]] gram[[i,j]] D[arg2, vars[[j]]],
  {i, 1, 3}, {j, 1, 3}]]];
lap[arg_] := Total[Table[lapVector[[i]] D[arg, vars[[i]]], {i, 1, 3}]] +
  Total[Flatten[Table[gram[[i,j]] D[arg, vars[[i]], vars[[j]]],
    {i, 1, 3}, {j, 1, 3}]]];
minOp[arg_] := gamma[arg, arg] lap[arg] - gamma[arg, gamma[arg, arg]]/2;
minimalRule = mm -> 1 + cc (ell - 2);
assertZero["splitting quartic gradient determinant", Det[gram] - 64 rr q (tt-q)];
ff = q - cc tt; aa = 1 - 2 cc; bb = cc (1 - cc);
assertZero["norm of defining quartic", gamma[ff, ff] - 4 rr (aa ff + bb tt)];
assertZero["Laplacian of defining quartic", (lap[ff] /. minimalRule) - 2 aa rr];
assertZero["minimal quartic identity", (minOp[ff] /. minimalRule) +
  8 ff (bb (rr^2 + 2 tt) + 2 aa ff)];
assertZero["two complementary quartics are orthogonal", gamma[q, tt-q]];
assertZero["complementary quartic norm", gamma[tt-q, tt-q] - 4 rr (tt-q)];

gg = tt ff;
ww = 4 rr tt (bb tt^2 + aa tt ff + 3 ff^2);
assertZero["degree-eight gradient norm", gamma[gg, gg] - ww];
assertZero["degree-eight Laplacian", (lap[gg] /. minimalRule) -
  rr (2 aa tt + 2 (ell+4) ff)];
defect8 = 8 tt ff (bb tt^2 ((ell-8) rr^2-4 tt) +
  aa tt ff ((ell-4) rr^2-4 tt) + 3 ff^2 (ell rr^2-4 tt));
assertZero["degree-eight divergence numerator", (minOp[gg] /. minimalRule) - defect8];
p10 = bb (ell-8+(2 ell-23) ss) +
  aa (ell-4+(2 ell-14) ss) zz +
  (3 ell+(15 ell-33) ss+(18 ell+27) ss^2) zz^2;
defect10 = 8 rr^5 tt^3 ff (p10 /. {ss -> tt/rr^2, zz -> ff/tt});
assertZero["degree-ten divergence numerator", (minOp[rr gg] /. minimalRule) - defect10];

p8 = bb (ell-8-4 ss) + aa (ell-4-4 ss) zz + 3 (ell-4 ss) zz^2;
assertZero["monotonicity in radius ratio", D[p8, ss] + 4 (bb+aa zz+3 zz^2)];
p8bal = p8 /. ss -> 1/4;
assertZero["endpoint for degree-eight positivity", (p8bal /. zz -> -cc) -
  4 cc ((ell-1) cc-1)];
disc8 = 12 (ell-1)(ell-9) cc (1-cc) - (ell-5)^2 (1-2 cc)^2;
cstar = (ell-5)/(8 (ell-2));
assertZero["entry discriminant", (disc8 /. cc -> cstar) -
  3 (ell-1)(ell-5)(ell^2-14 ell+21)/(4 (ell-2)^2)];
assertZero["dimension twelve endpoint discriminant", (disc8 /. {ell -> 12, cc -> 1/10}) - 107/25];
disc10 = 4 bb (2-3 ss)(30+117 ss+207 ss^2) - aa^2 (6+6 ss)^2;
assertZero["dimension ten worst imbalance", (disc10 /. cc -> 1/8) -
  (6 + 45 ss/2 + 117 ss^2/16 - 4347 ss^3/16)];

wf = gamma[ff, ff];
hessSquared = lap[wf]/2 - gamma[ff, lap[ff]];
normalDerivativeSquared = gamma[wf, wf]/(4 wf);
hessNormal = gamma[ff, wf]/(2 wf);
curvature = (hessSquared - 2 normalDerivativeSquared + hessNormal^2)/wf;
assertZero["curvature on the quartic cone", ((curvature /. minimalRule) /. q -> cc u v) -
  ((ell-2) rr/tt + 2 (ell-1)/rr)];
hessTNormal = (gamma[ff, gamma[ff, tt]] - gamma[tt, wf]/2)/wf;
assertZero["closure formula for normal Hessian of T",
  (hessTNormal /. q -> cc u v) - (2 rr-4 tt/rr)];
linkLap[h_] := (2 (ell-1)-8 ell tvar) D[h, tvar] +
  4 tvar (1-4 tvar) D[h, {tvar, 2}];
linkPotential = (ell-2)/tvar + 2 (ell-1);
assertZero["dimension sixteen angular eigenfunction", ((linkLap[1/tvar] +
  linkPotential/tvar) /. ell -> 8) - 46/tvar];
assertZero["dimension sixteen negative Hardy margin", (13/2)^2 - 46 + 15/4];

Do[
  data = item; degree = data[[1]]; kappa = data[[2]]; level = data[[3]];
  vars = {u, v, pvar};
  gram = {{4u, 0, 2 degree pvar}, {0, 4v, 2 degree pvar},
    {2 degree pvar, 2 degree pvar, degree^2 rr tt^(degree-1)}};
  lapVector = {16, 16, If[kappa == 0, 0, kappa rr tt^(degree/2-1)]};
  sourceF = pvar - If[level == 0, 0, level tt^(degree/2)];
  assertZero["octonionic minimal level, g=" <> ToString[degree],
    minOp[sourceF] /. pvar -> If[level == 0, 0, level tt^(degree/2)]];
, {item, {{2,-4,-1/3}, {3,0,0}, {4,8,1/3}, {6,0,0}}}];
];

CheckHopfPairing[] := Module[
  {u, v, f, r, vars, matrix, drift, rad, prod, potential,
   expected, fkm, ivars, imatrix, idrift, intrinsic},
vars = {u, v, f};
  matrix = {{4 u, 0, 4 f}, {0, 4 v, 4 f},
    {4 f, 4 f, 4 u v (u + v)}};
  drift = {32, 32, 0};
  rad = u + v; prod = u v; potential = prod f;
  CheckEqual[CurvatureExpr[f, vars, matrix, drift],
    -8 (rad^2 + 2 prod) f, "Hopf quartic minimality"];
  CheckEqual[LaplaceExpr[potential, vars, matrix, drift],
    40 rad f, "Hopf ambient Laplacian"];
  CheckEqual[GammaExpr[potential, potential, vars, matrix],
    4 rad prod (prod^2 + 3 f^2), "Hopf ambient gradient"];
  expected = 32 prod f (prod^2 (2 rad^2 - prod) +
    3 f^2 (4 rad^2 - prod));
  CheckEqual[CurvatureExpr[potential, vars, matrix, drift],
    expected, "Hopf ambient defect"];
  CheckEqual[(rad Cancel[expected/f]/
      (prod^3 GammaExpr[f, f, vars, matrix])) /. f -> 0,
    16 rad^2/prod - 8, "Hopf Jacobi boundary weight"];
  fkm = 4 f - (u - v)^2;
  CheckEqual[GammaExpr[fkm, fkm, vars, matrix],
    16 rad^3, "Hopf FKM gradient"];
  CheckEqual[LaplaceExpr[fkm, vars, matrix, drift],
    -8 rad, "Hopf FKM Laplacian"];
  ivars = {r, f};
  imatrix = {{4 r, 8 f}, {8 f, r^3}};
  idrift = {62, 4 f/r}; intrinsic = r f;
  CheckEqual[LaplaceExpr[intrinsic, ivars, imatrix, idrift],
    82 f, "Hopf intrinsic Laplacian"];
  CheckEqual[GammaExpr[intrinsic, intrinsic, ivars, imatrix],
    r^5 + 20 r f^2, "Hopf intrinsic gradient"];
  CheckEqual[CurvatureExpr[intrinsic, ivars, imatrix, idrift],
    16 r f (2 r^4 + 85 f^2), "Hopf intrinsic defect"];
];

CheckChangedHopfGenerators[] := Module[
  {u, v, f, p, q, xx, yy, vars, matrix, drift, changed,
   rad, split, quad, expected, expectedLap, i, j, bigU, bigV,
   defining, radius},
vars = {u, v, f};
  matrix = {{4 u, 0, 4 f}, {0, 4 v, 4 f},
    {4 f, 4 f, 4 u v (u + v)}};
  drift = {32, 32, 0};
  rad = u + v; split = u - v; quad = u^2 + v^2 - 2 f;
  changed = {rad, split, quad};
  expected = {{4 rad, 4 split, 8 quad},
    {4 split, 4 rad, 8 rad split},
    {8 quad, 8 rad split, 16 rad quad}};
  expectedLap = {64, 0, 72 rad};
  Do[CheckEqual[
    GammaExpr[changed[[i]], changed[[j]], vars, matrix],
    expected[[i, j]], "changed Gamma " <> ToString[{i, j}]],
    {i, 3}, {j, i, 3}];
  Do[CheckEqual[LaplaceExpr[changed[[i]], vars, matrix, drift],
    expectedLap[[i]], "changed Laplacian " <> ToString[i]], {i, 3}];
];

CheckHopfIncidence[] := Module[
  {rr,e0,e1,e2,r0,r1,r2,nDim,a0,a1,a2,vars,gram,drift,ff,
   sigma,second,third,gamma,lap,minOp,q,f,w,hnorm,gnorm,hfn,
   hq,dcq,linkGrad,linkDrift,phi,curv,u,v,d},
gamma[x_,y_] := GammaExpr[x,y,vars,gram];
lap[x_] := LaplaceExpr[x,vars,gram,drift];
minOp[x_] := CurvatureExpr[x,vars,gram,drift];
  vars={u,v,f};
  gram={{4u,0,4f},{0,4v,2f},{4f,2f,u^2+4u v}};
  drift={4d,2(d+1),0};
  CheckEqual[minOp[f],-8(u+v)f,"Hopf incidence minimal cubic"];
  CheckEqual[Det[gram],16(u+4v)(u^2 v-f^2),"Hopf incidence gradient determinant"];
  CheckEqual[lap[u+v],2(3d+1),"Hopf incidence ambient dimension"];
];

CheckBihomogeneousHardyCoefficients[] := Module[
 {r,t,f,c,ell,a,b,vars,matrix,drift,mult,expected,defect,beta,data},
 vars={r,t,f}; a=1-2c; b=c(1-c);
 matrix={{4r,8t,8f},{8t,4r t,4r f},{8f,4r f,4r(a f+b t)}};
 drift={4ell,2ell r,2a r};
 Do[
  mult=data[[1]]; expected=data[[2]];
  defect=Factor[Together[CurvatureExpr[mult f,vars,matrix,drift]/f]];
  beta=Together[(r defect/(mult^3 GammaExpr[f,f,vars,matrix]))/.f->0];
  CheckEqual[beta,expected,"bihomogeneous Hardy boundary coefficient: "<>data[[3]]],
  {data,{{t,2(ell-8)r^2/t-8,"TF"},{r t,2(ell-8)r^2/t+4ell-46,"RTF"}}}];
 CheckEqual[(8-(2ell-1)/2)^2+8(ell-8)-8,ell^2-9ell+1/4,
  "bihomogeneous Hardy lower bound from TF"];
 CheckEqual[(10-(2ell-1)/2)^2+8(ell-8)+4ell-46,ell^2-9ell+1/4,
  "bihomogeneous Hardy lower bound from RTF"];
];



Block[{$IterationLimit = Infinity, $Assumptions = True},
  RunVerificationGroup["CheckBihomogeneous", CheckBihomogeneous];
  RunVerificationGroup["CheckHopfPairing", CheckHopfPairing];
  RunVerificationGroup["CheckChangedHopfGenerators", CheckChangedHopfGenerators];
  RunVerificationGroup["CheckHopfIncidence", CheckHopfIncidence];
  RunVerificationGroup["CheckBihomogeneousHardyCoefficients", CheckBihomogeneousHardyCoefficients];
];

CertificateResult = <|"CertificateID" -> "LSC-BIHOMOGENEOUS-HOPF",
  "ExpectedChecks" -> 49,
  "Checks" -> (<|"Name" -> #[[1]], "Passed" -> TrueQ[#[[2]]],
      "Detail" -> If[StringQ[#[[3]]], #[[3]], ToString[#[[3]], InputForm]]|> &) /@ VerificationLog,
  "EvaluationOK" -> TrueQ[EvaluationOK],
  "EvaluationMessages" -> EvaluationMessages|>;
End[];
EndPackage[];
LaplacianCertificates`BihomogeneousHopf`Private`CertificateResult
