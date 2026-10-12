(* ::Package:: *)
(* LSC-CLIFFORD-HOMOGENEITY -- Clifford Homogeneity.
   Compatible with Mathematica 12. No external package or file dependencies.
   Get returns checks; only the central runner writes execution records.
   Supplemental homogeneity derivation; no current numbered statement is assigned.
   Weighted Euler identities check homogeneity of the full rational master quotient.
   All four endpoint denominators are cleared symbolically and the resulting polynomial expressions are finite at H=0.
   No value is assigned to the normalized gradient at H=0, and these algebra checks do not replace the analytic comparison argument.
*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`CliffordHomogeneity`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`CliffordHomogeneity`Private`*"];

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


(* Supplemental homogeneity and denominator clearing, independent of
   LSC-CLIFFORD-CUBIC. jGradS = J = Gamma(S,S), ds = K = Delta S.
   Nonharmonic certificates are covered by LSC-CUBIC-CERTIFICATES.
   This supplemental derivation has no current numbered manuscript counterpart.
   Algebraic checks do not prove the analytic singular-set
   comparison or the continuity argument; those remain in the manuscript. *)
CheckCliffordHomogeneityAppendix[] := Module[{al, cc, dd, deg, dhat, dim, ds, ee, endpointPolynomials, eq, f2, ga, hh, hhat, hssUU, hssUZ, hssZZ, j, jGradS, lam, m, mm, normalized, omega, p, pg1, pg2, q, qhat, qq, radius, rr, ss, t, theta, truth, u, v, vhat, vv, weightedEuler, lhs, rhs, label, cond, allPassed},
 allPassed=True;
 eq[lhs_,rhs_,label_] := If[!TrueQ[CheckEqual[lhs,rhs,label]],allPassed=False];
 truth[cond_,label_] := If[!TrueQ[CheckTrue[cond,label]],allPassed=False];
rr = u + v; ss = 4 u v + t;
jGradS = 16 (4 u v rr + (u + 4 v) t);
mm = 4 u^2 + 12 u v + 4 v^2 + t;
ds = 16 m v + 16 (q + 1) u;
dim = 2 m + q + 1; deg = 3 + 4 ga;
cc = ga/ss; ee = ga (ga - 1)/ss^2;
vv = ss + 32 cc rr f2 + cc^2 jGradS f2;
dd = ga (ds + 32 rr)/ss + ga (ga - 1) jGradS/ss^2;
hssUU = 32 u v^2 + 32 u^2 v + 8 u t - 16 v t + 96 f2;
hssUZ = 128 u^2 + 640 u v + 128 v^2 + 96 t;
hssZZ = 512 u v^3 + 2048 u^2 v^2 + 512 u^3 v
  + 1280 v^2 t + 1920 u v t + 128 u^2 t + 64 t^2;
qq = dd vv - (8 rr + cc jGradS + 32 cc^2 f2 mm)
  - 2 cc (16 rr + cc jGradS) (ss + 16 cc rr f2)
  - cc (hssUU + 2 cc f2 hssUZ + cc^2 f2 hssZZ)
  - ee f2 (16 rr + cc jGradS)^2;
vhat = vv + 4 al (deg + al) f2/rr;
dhat = dd + 2 al (2 deg + dim + 2 al - 2)/rr;
hhat = vv dd - qq + 4 al (deg - 1) vv/rr
  + 4 al^2 deg (deg - 1) f2/rr^2
  + 4 al (deg + 2 al) (vv + 2 al deg f2/rr)/rr
  + 2 al vhat/rr + 4 al (al - 1) (deg + 2 al)^2 f2/rr^2;
qhat = qq + 4 al (deg + al) f2 dd/rr
  + 2 al (dim - 2 deg - 2 al - 1) vv/rr
  + 4 al (deg + al) (deg + 2 al (dim - 1)) f2/rr^2;

 weightedEuler[p_] := 2 u D[p,u]+2 v D[p,v]+4 t D[p,t]+6 f2 D[p,f2];
 eq[3 (2 al+4 ga+3)-4,6 al+12 ga+5,"F: degree of L(G)"];
 eq[2(3 al+1)+4(3 ga)+3,6 al+12 ga+5,"F: degree of prefactor"];
 eq[(2 al+4 ga+3)/.{al->0,ga->3/4},6,"F: degree of G1"];
 eq[(2 al+4 ga+3)/.{al->1/2,ga->1/2},6,"F: degree of G2"];
 eq[weightedEuler[vv],4 vv,"F: V_gamma has weighted degree four"];
 eq[weightedEuler[dd],-2 dd,"F: D_gamma has weighted degree minus two"];
 eq[weightedEuler[qq],2 qq,"F: Q_gamma has weighted degree two"];
 eq[weightedEuler[qhat],2 qhat,"F: Qhat has weighted degree two"];
 eq[weightedEuler[qhat/rr],0,"F: normalized quotient has degree zero"];
 normalized={u->lam,v->1-lam,t->theta lam^2,
   f2->omega theta lam^2 (1-lam)};
 hh=4(1-lam)+theta lam;
 eq[ss/.normalized,lam hh,"F: S/R squared equals lambda H"];
 eq[(ss/.{u->radius lam,v->radius(1-lam),t->radius^2 theta lam^2}),
   radius^2 lam hh,"F: S=R squared lambda H before normalization"];
 eq[hh/.lam->0,4,"F: H on lambda=0"];
 eq[hh/.lam->1,theta,"F: H on lambda=1"];
 eq[hh/.theta->0,4(1-lam),"F: H on theta=0"];
 eq[hh/.theta->1,4-3lam,"F: H on theta=1"];
 pg1=Together[qhat/.{al->0,ga->3/4}/.normalized];
 pg2=Together[qhat/.{al->1/2,ga->1/2}/.normalized];
 endpointPolynomials={Cancel[(pg1/.omega->0)hh],
   Cancel[(pg1/.omega->1)hh^3],
   Cancel[(pg2/.omega->0)hh],
   Cancel[(pg2/.omega->1)hh^3]};
 Do[truth[PolynomialQ[endpointPolynomials[[j]],{lam,theta}],
   "F: cleared endpoint is polynomial, case "<>ToString[j]],{j,4}];
 Do[eq[Together[endpointPolynomials[[j]]/
   If[MemberQ[{1,3},j],hh,hh^3]],
   {pg1/.omega->0,pg1/.omega->1,pg2/.omega->0,pg2/.omega->1}[[j]],
   "F: denominator clearing identity, case "<>ToString[j]],{j,4}];
 (* Endpoint polynomials remain finite at H=0.  This does not assign
    values to the rational quotients or the normalized gradient there. *)
 Do[truth[FreeQ[endpointPolynomials[[j]]/.{lam->1,theta->0},
   Indeterminate|ComplexInfinity|DirectedInfinity],
   "F: polynomial value at H=0 is finite, case "<>ToString[j]],{j,4}];

 allPassed
];



Block[{$IterationLimit = Infinity, $Assumptions = True},
  RunVerificationGroup["CheckCliffordHomogeneityAppendix", CheckCliffordHomogeneityAppendix];
];

CertificateResult = <|"CertificateID" -> "LSC-CLIFFORD-HOMOGENEITY",
  "ExpectedChecks" -> 27,
  "Checks" -> (<|"Name" -> #[[1]], "Passed" -> TrueQ[#[[2]]],
      "Detail" -> If[StringQ[#[[3]]], #[[3]], ToString[#[[3]], InputForm]]|> &) /@ VerificationLog,
  "EvaluationOK" -> TrueQ[EvaluationOK],
  "EvaluationMessages" -> EvaluationMessages|>;
End[];
EndPackage[];
LaplacianCertificates`CliffordHomogeneity`Private`CertificateResult
