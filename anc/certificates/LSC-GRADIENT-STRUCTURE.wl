(* Mathematica 12 certificate: LSC-GRADIENT-STRUCTURE.
   Exact identities conditional on the manuscript differential table.
   It does not prove geometric existence or global subcalibration positivity. *)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`GradientStructure`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`GradientStructure`Private`*"];
checks = {};
zero[name_String, expression_] := Module[{residual},
  residual = FullSimplify[Together[expression],
    Assumptions -> Element[{rr,ff,ss,dd,cc,lam,mu,nn,sigma,alpha,beta,eta,kk,hh,bb,uu,aa,gamma,rho,ee,xx,gg,eps}, Reals] &&
      ee>=3 && gg>0 && nn>=3 && rr>0 && ss>0 && lam<0 && kk>0 && hh>=8 && bb>=1 && 0<uu<1 && 0<kk rho<Pi/2];
  AppendTo[checks, <|"Name" -> name, "Passed" -> TrueQ[residual === 0],
    "Detail" -> "Exact residual: " <> ToString[residual, InputForm]|>]];

zero["radial_potential_gradient_square",
  dd*ff^2*sigma - ff^2*sigma*(2*dd + sigma) + (1/2)*ff*sigma*(2*dd*ff + 2*ff*sigma)];
zero["radial_potential_Laplacian",
  2*dd*ff*sigma + ff*nn*sigma + ff*sigma*(sigma - 2) - ff*sigma*(2*dd + nn + sigma - 2)];
zero["radial_potential_normalized_defect",
  -ff^2*sigma*(dd + sigma*(nn - 1))*(2*dd + sigma) - ff*(cc*sigma*(2*dd + sigma) + mu) - lam - sigma*ss*(-2*dd + nn - sigma - 1) + (-cc*ss + ff^2*mu + ff*lam - 1/4*ff*sigma*ss*(4*dd - 4) - 1/2*(2*dd*ff + 2*ff*sigma)*((1/2)*dd*ff^2*sigma^2 - dd*ff^2*sigma + (1/4)*ff*sigma^2*(2*dd*ff + 2*ff*sigma) - 1/2*ff*sigma*(2*dd*ff + 2*ff*sigma) + (1/2)*ff*sigma*(dd*ff*sigma + ff*sigma^2) + (1/2)*sigma*ss + (1/2)*sigma*(dd*ff^2*sigma + ss)) - 1/2*(dd*ff^2*sigma + ss)*(2*dd*ff*sigma + (1/2)*ff*sigma*(2*dd + 2*sigma) + (1/2)*sigma*(2*dd*ff + 2*ff*sigma)) + (dd*ff^2*sigma + (1/2)*ff*sigma*(2*dd*ff + 2*ff*sigma) + ss)*(cc + 2*dd*ff*sigma + ff*nn*sigma + ff*sigma*(sigma - 2)))/ff];
zero["Bochner_Hessian_norm",
  (1/2)*beta*ff - 2*cc*dd*ff*((1/2)*dd - 1) - ff*((1/2)*beta - cc*dd*(dd - 2))];
zero["critical_quartic_reduction",
  2*cc*eta*lam + 4*eta^4 - eta^2*(4*cc^2 - 2*lam) - (2*cc + 2*eta)*(-2*cc*eta^2 + 2*eta^3 + eta*lam)];
zero["critical_nonzero_eigenvalue",
  -2*cc*((1/2)*cc + (1/2)*(cc^2 - 2*lam)^(1/2)) + lam + 2*((1/2)*cc + (1/2)*(cc^2 - 2*lam)^(1/2))^2];
zero["cone_gradient_square",
  -4*cc^2*ss + 2*lam*ss + ss*(4*cc^2 - 2*lam)];
zero["cone_Laplacian",
  -hh*lam - lam*(1 - hh) + lam];
zero["cone_second_fundamental_norm",
  (1/2)*lam*(hh - 2)/ss + (-1/2*hh*lam + lam)/ss];
zero["normalized_link_gradient",
  4*kk^4*(-lam^2*uu^2/kk^2 + lam^2*uu/kk^2)/lam^2 - 4*kk^2*uu*(1 - uu)];
zero["normalized_link_Laplacian",
  -2*kk^2*(hh - 1) - 2*kk^2*(lam*(1 - hh) + lam*uu*(2*kk + nn - 3)/kk)/lam + 2*kk*uu*(2*kk + nn - 3)];
zero["maximum_rank_dimension_identity",
  2*bb*kk^2 + 2*kk^2*(hh - 1) - 2*kk*(kk*(bb + hh - 3) + 2*kk)];
zero["distance_Laplacian",
  -kk*uu^(1/2)*(1 - uu)^(1/2)*((uu - 1)^(-1) + uu^(-1)) - kk*(-uu*(bb - 1) + (1 - uu)*(hh - 2))/(uu^(1/2)*(1 - uu)^(1/2)) + kk*(hh - uu*(bb + hh - 1) - 1)/(uu^(1/2)*(1 - uu)^(1/2))];
zero["normal_flow_volume_density",
  -kk*(bb - 1)*Sin[kk*rho]/Cos[kk*rho] + kk*(hh - 2)*Cos[kk*rho]/Sin[kk*rho] - kk*(-(bb - 1)*Tan[kk*rho] + (hh - 2)*Cot[kk*rho])];
zero["intrinsic_power_Laplacian",
  4*aa*kk^2*(1 - uu)*(aa - 1)/uu + 2*aa*kk^2*(hh - uu*(bb + hh - 1) - 1)/uu - kk^2*(-4*aa^2 - aa*(2*bb + 2*hh - 6) + (4*aa^2 + aa*(2*hh - 6))/uu)];
zero["cone_gradient_power_Jacobi",
  -2*aa*lam*(aa - 1)/ss - aa*lam*(hh - 1)/ss - 1/2*lam*(hh - 2)/ss + lam*(2*aa^2 + aa*(hh - 3) + (1/2)*hh - 1)/ss];
zero["cone_separated_Jacobi",
  4*aa*kk^2*(1 - uu)*(aa - 1)/uu + 2*aa*kk^2*(hh - uu*(bb + hh - 1) - 1)/uu - kk^2*(-4*aa^2 - aa*(2*bb + 2*hh - 6) + (4*aa^2 + aa*(2*hh - 6) + hh - 2)/uu) + kk^2*(hh - 2)/uu];
zero["optimal_indicial_root",
  hh + (2*hh - 6)*(-1/4*hh + (1/4)*(hh^2 - 10*hh + 17)^(1/2) + 3/4) + 4*(-1/4*hh + (1/4)*(hh^2 - 10*hh + 17)^(1/2) + 3/4)^2 - 2];
zero["optimal_multiplier_Jacobi_gap",
  (1/4)*kk^2*(bb + (hh^2 - 10*hh + 17)^(1/2))^2 - 1/4*kk^2*(bb + hh - 3)^2 + 4*kk^2*(1 - uu)*(-1/4*hh + (1/4)*(hh^2 - 10*hh + 17)^(1/2) - 1/4)*(-1/4*hh + (1/4)*(hh^2 - 10*hh + 17)^(1/2) + 3/4)/uu + kk^2*(hh - 2)/uu + 2*kk^2*(-1/4*hh + (1/4)*(hh^2 - 10*hh + 17)^(1/2) + 3/4)*(hh - uu*(bb + hh - 1) - 1)/uu];
zero["optimal_multiplier_degree",
  (1/2)*kk*(bb + (hh^2 - 10*hh + 17)^(1/2)) - 1/2*kk*(bb + hh - 3) + 2*kk*((1/4)*hh - 1/4*(hh^2 - 10*hh + 17)^(1/2) - 5/4) + kk];
zero["optimal_normal_speed_degree",
  -1/2*kk*(bb + (hh^2 - 10*hh + 17)^(1/2)) + (1/2)*kk*(bb + hh - 3) - 2*kk*((1/4)*hh - 1/4*(hh^2 - 10*hh + 17)^(1/2) - 3/4)];
zero["optimal_positive_exponent_identity",
  -hh^2 + 10*hh + (hh - 5)^2 - 25];
zero["odd_flow_first_integral",
  2*ee^2*ff*(ff^2 - gg/ee - (1/2)/(ee^3*(ee - 1))) + 2*ff*gg*(ee - 1) + 2*ff*(-ee^2*ff^2 + gg + (1/2)/(ee*(ee - 1)))];
zero["odd_minimum_Q_elimination",
  2*ee^3*(ee - 1)*((1/2)*eps*xx/(ee*(ee - 1)) - 1/2*eps/(ee*(ee - 1)) - xx^2/ee^2 + ee^(-2) - 1/2*(-ee*eps + 4*ee - 6)/(ee^3*(ee - 1))) - 2*ee^2 - ee*eps*(ee*xx - ee + 1) + 2*ee*xx^2*(ee - 1) + 6*ee - 6];
zero["odd_minimum_Hessian_elimination",
  -ee*eps*(ee*xx - ee + 1) + 2*ee*(ee - 1)*((1 - xx)^2 - 1/2*(1 - xx)*(eps - 4*xx*(ee - 1)/ee) - 1/2*xx*(-ee*eps + 4*ee - 6)/(ee*(ee - 1))) - 2*ee*(ee - 1) + xx^2*(ee - 2)*(2*ee - 2) - xx*(10 - 8*ee)];
zero["odd_elimination_factor",
  -2*ee^2 + 2*ee*xx^2*(ee - 1) + 2*ee*(ee - 1) + 6*ee - xx^2*(ee - 2)*(2*ee - 2) + xx*(10 - 8*ee) - (2*xx - 2)*(-2*ee + xx*(2*ee - 2) + 3) - 6];
zero["odd_maximum_Hessian_delta",
  ((1/2)*eps - 2*(ee - 1)/ee)/ee + (1/2)*(-ee*eps + 4*ee - 6)/ee^2 + ee^(-2)];
zero["odd_forced_epsilon",
  2*ee^2 - 2*ee*(ee - 1)*(2*ee - 3)^2/(2*ee - 2)^2 - 6*ee + 6 + (4*ee^2 - 15*ee + 12)*(ee*(2*ee - 3)/(2*ee - 2) - ee + 1)/(ee - 2)];
zero["odd_forced_delta",
  4*ee - ee/(ee - 2) - 6 - (4*ee^2 - 15*ee + 12)/(ee - 2)];
zero["odd_forced_Q_derivative",
  (ee - 2)^(-1) - 4*(ee - 1)*(2*ee - 3)/(ee*(2*ee - 2)) + (4*ee^2 - 15*ee + 12)/(ee*(ee - 2))];
zero["odd_two_rank_separation",
  -(ee + nn)/(2*ee - 2) - (ee + nn - 4)/(2*ee - 2) + (ee + nn - 2)/(ee - 1)];
zero["stable_dimension_slack",
  -kk*(bb + hh - 9) + kk*(bb + hh - 3) - 6*kk];

CertificateResult = <|"CertificateID" -> "LSC-GRADIENT-STRUCTURE",
  "ExpectedChecks" -> 32, "Checks" -> checks|>;
End[];
EndPackage[];
LaplacianCertificates`GradientStructure`Private`CertificateResult
