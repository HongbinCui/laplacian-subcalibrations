(* ::Package:: *)
(* LSC-GENERAL-QUARTIC: conditional exact orthogonal-quartic certificates.
   Current manuscript: orthogonal quartic data and exact certificates; see CertificateIndex.md.
   Mathematica 12 compatible. The input differential table and invariant-image
   inequalities are hypotheses. Geometry, realizability and removability are
   proved in the paper. The h=9 extension is checked separately in
   LSC-REVIEW-STRENGTHENINGS. No previous execution record is imported.
   Domain: R>0, 0<T<=R^2/4, -c<=F/T<=1-c. The positive range is
   integer h>=10, r0>=1, 1/(h-2)<=c<=1-1/(h-2), or h=8,c=1/2,r0>=5.
*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`GeneralQuartic`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`GeneralQuartic`Private`*"];
checks = {};
expectedChecks = 46;
zero[name_, expression_] := Module[{residual},
  residual = Expand[Numerator[Together[expression]]];
  AppendTo[checks, <|"name" -> name, "passed" -> TrueQ[residual === 0],
    "residual" -> ToString[residual, InputForm]|>]];
positiveCoefficients[name_, expression_, variables_] := Module[{coefficients},
  coefficients = Last /@ CoefficientRules[Expand[expression], variables];
  AppendTo[checks, <|"name" -> name,
    "passed" -> (Length[coefficients] > 0 && And @@ (TrueQ[# > 0] & /@ coefficients)),
    "polynomial" -> ToString[Expand[expression], InputForm],
    "variables_assumed_nonnegative" -> (ToString /@ variables),
    "nonzero_coefficients" -> (ToString[#, InputForm] & /@ coefficients)|>]];

variables = {rr, tt, ff};
aa = 1-2 cc; bb = cc (1-cc); ambient = 2 (r0+hh-1);
gammaTable = {{4 rr,8 tt,8 ff},{8 tt,4 rr tt,4 rr ff},
  {8 ff,4 rr ff,4 rr (aa ff+bb tt)}};
laplacians = {2 ambient,2 hh rr,2 aa rr};
gradient[expression_] := (D[expression,#] & /@ variables);
gamma[left_,right_] := Expand[gradient[left].gammaTable.gradient[right]];
laplacian[expression_] := Expand[gradient[expression].laplacians +
  Sum[gammaTable[[i,j]] D[expression,variables[[i]],variables[[j]]],
    {i,3},{j,3}]];
defect[expression_] := Module[{square=gamma[expression,expression]},
  Expand[square laplacian[expression]-gamma[expression,square]/2]];

zero["minimal_quartic_defect", defect[ff]+8 ff (bb (rr^2+2 tt)+2 aa ff)];
zero["image_linear_positive_decomposition",
  bb+aa zz-((1-cc)^2 (zz+cc)+cc^2 (1-cc-zz))];
p8 = bb (hh-8-4 ss)+aa (hh-4-4 ss) zz+3 (hh-4 ss) zz^2;
base10 = bb (hh-8+(2 hh-23) ss)+aa (hh-4+(2 hh-14) ss) zz+
  (3 hh+(15 hh-33) ss+(18 hh+27) ss^2) zz^2;
correction = 2 (r0-1) ss (bb+aa zz+3 (1+3 ss) zz^2);
p10 = base10+correction;
g8 = tt ff; g10 = rr tt ff; l8 = defect[g8]; l10 = defect[g10];
zero["degree_eight_gradient_square",gamma[g8,g8]-4 rr tt (bb tt^2+aa tt ff+3 ff^2)];
zero["degree_eight_Laplacian",laplacian[g8]-2 rr (aa tt+(hh+4) ff)];
zero["degree_eight_normalized_defect",l8-8 rr^2 tt^3 ff (p8/.{ss->tt/rr^2,zz->ff/tt})];
zero["degree_ten_gradient_square",
  gamma[g10,g10]-4 rr tt (rr^2 (bb tt^2+aa tt ff+3 ff^2)+9 tt ff^2)];
zero["degree_ten_Laplacian",laplacian[g10]-2 aa rr^2 tt-
  2 (hh+4) rr^2 ff-(2 ambient+32) tt ff];
zero["degree_ten_normalized_defect",l10-8 rr^5 tt^3 ff (p10/.{ss->tt/rr^2,zz->ff/tt})];
Do[zero["potential_homogeneous_degree_"<>ToString[item[[1]]],
  Total[MapThread[#1 #2 #3 &,{ {2,4,4},variables,gradient[item[[2]]]}]]-
  item[[1]] item[[2]]],{item,{{8,g8},{10,g10}}}];
zero["degree_ten_dimension_correction",l10-(l10/.r0->1)-
  4 (r0-1) tt ff gamma[g10,g10]];

(* Cancel ff before restricting to ff=0; rr,tt must be positive. *)
beta8 = Cancel[(rr Cancel[l8/ff]/(tt^3 gamma[ff,ff]))/.ff->0];
beta10 = Cancel[(rr Cancel[l10/ff]/((rr tt)^3 gamma[ff,ff]))/.ff->0];
zero["degree_eight_boundary_beta",beta8-(2 (hh-8) rr^2/tt-8)];
zero["degree_ten_boundary_beta",beta10-(2 (hh-8) rr^2/tt+4 (r0+hh)-50)];
beta08 = 8 hh-72; beta010 = 12 hh+4 r0-114;
zero["degree_eight_beta_lower_bound_remainder",
  beta8-beta08-2 (hh-8) (rr^2-4 tt)/tt];
zero["degree_ten_beta_lower_bound_remainder",
  beta10-beta010-2 (hh-8) (rr^2-4 tt)/tt];
hardy = (ambient-17)^2/4+8 hh-72;
zero["degree_eight_Hardy_coefficient",(8-(ambient-1)/2)^2+beta08-hardy];
zero["degree_ten_Hardy_coefficient",(10-(ambient-1)/2)^2+beta010-hardy];
positiveCoefficients["Hardy_positive_h_at_least_ten",hardy/.{hh->10+uu,r0->1+vv},{uu,vv}];
positiveCoefficients["Hardy_positive_balanced_h_eight",hardy/.{hh->8,r0->5+vv},{vv}];

(* Interval proof for integer hh>=12, cc in [1/(hh-2),1-1/(hh-2)]. *)
qq = p8/.ss->1/4;
zero["degree_eight_monotonicity_in_s",D[p8,ss]+4 (bb+aa zz+3 zz^2)];
zero["degree_eight_reflection_symmetry",p8-(p8/.{cc->1-cc,zz->-zz})];
zero["degree_eight_left_endpoint",(qq/.zz->-cc)-4 cc ((hh-1) cc-1)];
cstar = (hh-5)/(8 (hh-2)); zstar = -aa (hh-5)/(6 (hh-1));
zero["degree_eight_vertex",D[qq,zz]/.zz->zstar];
zero["degree_eight_vertex_enters_interval",zstar+cc-4 (hh-2) (cc-cstar)/(3 (hh-1))];
d8 = 12 (hh-1) (hh-9) bb-(hh-5)^2 aa^2;
zero["degree_eight_negative_discriminant",-Discriminant[qq,zz]-d8];
zero["degree_eight_discriminant_monotonicity",D[d8,cc]-16 (1-2 cc) (hh^2-10 hh+13)];
positiveCoefficients["degree_eight_discriminant_derivative_factor",(hh^2-10 hh+13)/.hh->12+uu,{uu}];
zero["degree_eight_discriminant_at_vertex_threshold",(d8/.cc->cstar)-
  3 (hh-1) (hh-5) (hh^2-14 hh+21)/(4 (hh-2)^2)];
positiveCoefficients["degree_eight_threshold_factor_h_at_least_thirteen",
  (hh^2-14 hh+21)/.hh->13+uu,{uu}];
zero["degree_eight_h_twelve_endpoint_discriminant",(d8/.{hh->12,cc->1/10})-107/25];

(* At hh=10 or 11, bound alpha=aa^2 and use 0<=ss<=1/4. *)
Do[Module[{kk=item[[1]],maxA2=item[[2]],bound=item[[3]],ap,bp,cp,disc,lower,
    cubicCoefficient,slackForm,prefix},
  prefix = "degree_ten_h_"<>ToString[kk]<>"_";
  ap = kk-8+(2 kk-23) ss; bp = kk-4+(2 kk-14) ss;
  cp = 3 kk+(15 kk-33) ss+(18 kk+27) ss^2;
  disc = Expand[(1-alpha) ap cp-alpha bp^2]; lower = Expand[disc/.alpha->maxA2];
  zero[prefix<>"negative_discriminant",-Discriminant[base10/.hh->kk,zz]-(disc/.alpha->aa^2)];
  zero[prefix<>"discriminant_parameter_remainder",disc-lower-(maxA2-alpha) (ap cp+bp^2)];
  positiveCoefficients[prefix<>"constant_factor",ap/.ss->1/4-ww,{ww}];
  positiveCoefficients[prefix<>"leading_factor",cp,{ss}];
  cubicCoefficient = Coefficient[lower,ss,3];
  slackForm = Expand[lower-cubicCoefficient ss^3+cubicCoefficient (1/64-uu)];
  zero[prefix<>"positive_lower_bound_decomposition",lower-(slackForm/.uu->1/64-ss^3)];
  zero[prefix<>"lower_bound_constant",(slackForm/.{ss->0,uu->0})-bound];
  positiveCoefficients[prefix<>"positive_discriminant_bound",slackForm,{ss,uu}]],
  {item,{{10,9/16,1797/1024},{11,49/81,1309/162}}}];
balanced = (2 r0-9) ss/4+(24+(6 r0+81) ss+(18 r0+153) ss^2) zz^2;
zero["balanced_h_eight_normalized_defect",(p10/.{hh->8,cc->1/2})-balanced];
positiveCoefficients["balanced_h_eight_r0_at_least_five_positive_coefficients",
  balanced/.{r0->5+uu,zz^2->vv},{ss,uu,vv}];

CertificateResult = <|
 "CertificateID" -> "LSC-GENERAL-QUARTIC", "ExpectedChecks" -> 46,
 "Checks" -> (Function[row, <|"Name" -> row["name"],
    "Passed" -> TrueQ[row["passed"]],
    "Detail" -> ToString[KeyDrop[row, {"name", "passed"}], InputForm]|>] /@ checks)
|>;
End[];
EndPackage[];
LaplacianCertificates`GeneralQuartic`Private`CertificateResult
