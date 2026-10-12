(* ::Package:: *)
(* LSC-QUADRATIC-COMPOSITIONS -- release 20261011_memo_v7.
   Mathematica 12.0 compatible, exact arithmetic, no external dependencies.
   Get returns one certificate association; only the central runner writes logs.
   Formal table identities are conditional on the differential tables proved
   in the manuscript. Coordinate identities in H^2 and C^2 are exact for those
   dimensions; finite integer matrix checks are explicitly finite tests.
   Irreducibility, realizability, removability, global minimization,
   noncongruence and tangent-current convergence are not computer-certified.
   Extra Hopf and initial-form identities are algebraic cross-checks only.
*)
BeginPackage["LaplacianCertificates`QuadraticCompositions`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`QuadraticCompositions`Private`*"];
checks = {};
expectedChecks = 158;
EvaluationMessages = {};
EvaluationOK = True;
Block[{System`$MessageList = {}, $IterationLimit = Infinity},
 evaluationResult = CheckAbort[Catch[
  zero[name_, expression_, group_] := Module[{residual},
    residual = Expand[Numerator[Together[expression]]];
    AppendTo[checks, <|"Name" -> name, "Group" -> group,
      "Passed" -> TrueQ[residual === 0],
      "Method" -> "exact polynomial/rational identity",
      "Detail" -> "Exact residual: " <> ToString[residual, InputForm]|>]];
  positiveCoefficients[name_, expression_, variables_, group_] := Module[{rules, coefficients, constant},
    rules = CoefficientRules[Expand[expression], variables];
    coefficients = Last /@ rules;
    constant = Expand[expression] /. Thread[variables -> 0];
    AppendTo[checks, <|"Name" -> name, "Group" -> group,
      "Passed" -> TrueQ[Length[coefficients] > 0 && And @@ (TrueQ[# > 0] & /@ coefficients) && constant > 0],
      "Method" -> "positive nonzero coefficients and positive constant",
      "VariablesAssumedNonnegative" -> (ToString[#, InputForm] & /@ variables),
      "Detail" -> "Polynomial with nonnegative indeterminates: " <> ToString[Expand[expression], InputForm]|>]];
  gamma[f_, g_] := Expand[(D[f, #] & /@ variables).metric.(D[g, #] & /@ variables)];
  laplace[f_] := Expand[Sum[D[f, variables[[i]]] laplacians[[i]], {i, Length[variables]}] +
    Sum[metric[[i,j]] D[f, variables[[i]], variables[[j]]], {i, Length[variables]}, {j, Length[variables]}]];
  defect[f_] := Module[{square = gamma[f,f]}, Expand[square laplace[f]-gamma[f,square]/2]];

  (* Start with the underlying R,W,Q table. *)
  variables = {rr, ww, qq};
  metric = {{4 rr,8 ww,8 qq},{8 ww,4 rr ww,4 rr qq},{8 qq,4 rr qq,4 rr qq}};
  laplacians = {4 ell,2 (ell-1) rr,2 mm rr};
  aa = ell-mm-2; bb = mm-1; dd = aa-bb;
  ff = (aa+bb) qq-bb ww; bigD = aa^2 qq+bb^2 (ww-qq);
  group = "orthogonal_quartic_universal_table";
  zero["quartic_delta_F",laplace[ff]-2 dd rr,group];
  zero["quartic_gamma_F_W",gamma[ff,ww]-4 rr ff,group];
  zero["quartic_gradient_square",gamma[ff,ff]-4 rr bigD,group];
  zero["quartic_minimality_defect",defect[ff]+8 ff (2 bigD+aa bb rr^2),group];
  zero["quartic_recover_W_without_polynomial_division",laplace[gamma[ff,ff]]-
    16 (ell+4) bigD-8 (dd^2+aa bb (ell-1)) rr^2,group];
  blocks = {rr^2-4 ww,4 qq,4 (ww-qq)}; weights = {2,mm,ell-mm-1};
  zero["quartic_block_sum",Total[blocks]-rr^2,group];
  Do[
    zero["quartic_block_laplacian_"<>ToString[i-1],laplace[blocks[[i]]]-8 weights[[i]] rr,group];
    Do[zero["quartic_block_gamma_"<>ToString[i-1]<>"_"<>ToString[j-1],
      gamma[blocks[[i]],blocks[[j]]]-If[i==j,16 rr blocks[[i]],0],group],{j,i,3}],{i,3}];
  gg = ww ff;
  polynomial = (3 (ell-1)-12 tt) zz^2+dd (ell-5-4 tt) zz+aa bb (ell-9-4 tt);
  zero["quartic_potential_laplacian",laplace[gg]-2 rr ((ell+3) ff+dd ww),group];
  zero["quartic_potential_gradient_square",gamma[gg,gg]-4 rr ww (3 ff^2+dd ff ww+aa bb ww^2),group];
  zero["quartic_potential_defect",defect[gg]-8 rr^2 ff ww^3 (polynomial/.{zz->ff/ww,tt->ww/rr^2}),group];
  zero["quartic_P_monotonicity",D[polynomial,tt]+4 (3 zz^2+dd zz+aa bb),group];
  zero["quartic_P_boundary",(polynomial/.tt->1/4)-(3 (ell-2) zz^2+dd (ell-6) zz+aa bb (ell-10)),group];
  zero["quartic_gradient_positive_decomposition",((dd zz+aa bb)/.zz->(aa+bb) uu-bb)-
    (aa^2 uu+bb^2 (1-uu)),group];
  beta = Cancel[rr Cancel[defect[gg]/ff]/(ww^3 gamma[ff,ff])];
  zero["quartic_boundary_Jacobi_potential",Cancel[beta/.qq->bb ww/(aa+bb)]-
    (2 (ell-9) rr^2/ww-8),group];
  zero["quartic_Hardy_coefficient",((2 ell-17)/2)^2+8 ell-80-(ell^2-9 ell-31/4),group];

  group = "quaternionic_family_strict_positivity";
  boundary = polynomial/.tt->1/4; p2 = Expand[boundary/.mm->2]; p3 = Expand[boundary/.mm->3];
  zero["m2_ell12_discriminant",Discriminant[p2/.ell->12,zz]+156,group];
  zero["m2_left_endpoint",(p2/.zz->-1)-4,group];
  zero["m2_left_derivative",(D[p2,zz]/.zz->-1)-(ell-3) (ell-14),group];
  positiveCoefficients["m2_ell_at_least16_increasing",D[p2,zz]/.{ell->16+vv,zz->-1+ss},{vv,ss},group];
  Do[zero["m3_ell"<>ToString[item[[1]]]<>"_discriminant",
    Discriminant[p3/.ell->item[[1]],zz]-item[[2]],group],{item,{{12,-2460},{16,-14076},{20,-31676}}}];
  zero["m3_shifted_positive_form",(p3/.zz->ss-2)-(3 (ell-2) ss^2+(ell-3) (ell-22) ss+8 (ell-1)),group];
  positiveCoefficients["m3_ell_at_least24_positive",p3/.{ell->24+vv,zz->ss-2},{vv,ss},group];
  zero["m3_ell12_completed_square",(p3/.ell->12)-(30 (zz+1/2)^2+41/2),group];
  zero["k3_Hardy_constant",((ell^2-9 ell-31/4)/.ell->12)-113/4,group];
  zero["k2_transverse_Hardy_gap",(6-2)^2/4-5+1,group];

  (* Independently differentiate x,y in H^2, using right multiplication. *)
  group = "quaternionic_coordinate_table_H2";
  xx = Array[x,8]; yy = Array[y,8]; coordinates = Join[xx,yy];
  cg[f_,g_] := Expand[Total[(D[f,#] D[g,#]) & /@ coordinates]];
  cl[f_] := Expand[Total[D[f,{#,2}] & /@ coordinates]];
  rightUnits[list_] := Transpose[(Function[block,With[
    {a0=block[[1]],b0=block[[2]],c0=block[[3]],d0=block[[4]]},
    {{-b0,a0,d0,-c0},{-c0,-d0,a0,b0},{-d0,c0,-b0,a0}}]]) /@ Partition[list,4]];
  units = Flatten /@ rightUnits[xx];
  rc = xx.xx+yy.yy; wc = (xx.xx) (yy.yy)-(xx.yy)^2; qis = (#.yy & /@ units);
  Do[qc = Total[Take[qis,j]^2];
    zero["coordinate_m"<>ToString[j]<>"_delta_Q",cl[qc]-2 j rc,group];
    zero["coordinate_m"<>ToString[j]<>"_gamma_QQ",cg[qc,qc]-4 rc qc,group];
    zero["coordinate_m"<>ToString[j]<>"_gamma_WQ",cg[wc,qc]-4 rc qc,group],{j,{2,3}}];
  zero["coordinate_delta_W",cl[wc]-14 rc,group];
  zero["coordinate_gamma_WW",cg[wc,wc]-4 rc wc,group];

  (* Universal complex-bilinear Hermitian Gram differential table. *)
  group = "complex_Gram_cycle_sextic";
  h = Array[hh,{3,3}]; variables = Flatten[h]; pairs = Tuples[Range[3],2];
  metric = Table[2 (If[pairs[[i,1]]==pairs[[j,2]],h[[pairs[[j,1]],pairs[[i,2]]]],0]+
    If[pairs[[i,2]]==pairs[[j,1]],h[[pairs[[i,1]],pairs[[j,2]]]],0]),{i,9},{j,9}];
  laplacians = (If[#[[1]]==#[[2]],4 kk,0]& /@ pairs);
  cycle = h[[1,2]] h[[2,3]] h[[3,1]]; cycleBar = h[[2,1]] h[[3,2]] h[[1,3]];
  fs = (cycle-cycleBar)/(2 I); ts = (cycle+cycleBar)/2; rs = Diagonal[h];
  edges = {h[[2,3]] h[[3,2]],h[[3,1]] h[[1,3]],h[[1,2]] h[[2,1]]};
  rsum = Total[rs]; esum = Total[edges];
  squareS = Sum[rs[[i]] edges[[i]] (esum-edges[[i]]),{i,3}]-2 ts esum;
  bigB = Sum[edges[[i]] (rs[[i]]+rs[[Mod[i,3]+1]]) (rs[[i]]+rs[[Mod[i+1,3]+1]])-edges[[i]]^2,{i,3}]-2 rsum ts;
  zero["sextic_harmonic_all_row_counts",laplace[fs],group];
  zero["sextic_gradient_square",gamma[fs,fs]-squareS,group];
  zero["sextic_gamma_F_S",gamma[fs,squareS]-4 bigB fs,group];
  zero["sextic_minimality_defect",defect[fs]+2 bigB fs,group];
  zero["sextic_delta_complex_cycle",laplace[cycle]-4 esum,group];
  zero["sextic_gamma_F_real_cycle",gamma[fs,ts]-2 esum fs,group];
  Do[
    zero["sextic_gamma_F_radius_"<>ToString[i],gamma[fs,rs[[i]]]-4 fs,group];
    zero["sextic_gamma_F_edge_"<>ToString[i],gamma[fs,edges[[i]]]-2 (rsum-rs[[i]]) fs,group],{i,3}];
  zero["sextic_regular_point_gradient_square",(squareS/.Thread[Flatten[h]->Flatten[{{1,1,1},{1,2,1},{1,1,1}}]])-2,group];

  group = "complex_Gram_coordinates_and_tangents";
  rv = Array[ar,6]; iv = Array[ai,6]; vectors = Partition[rv+I iv,2];
  hv = Table[ComplexExpand[Conjugate[vectors[[i]]].vectors[[j]]],{i,3},{j,3}];
  cycleIm = Expand[ComplexExpand[Im[hv[[1,2]] hv[[2,3]] hv[[3,1]]]]];
  hopf[list_] := ComplexExpand[{2 Re[Conjugate[list[[1]]] list[[2]]],
    2 Im[Conjugate[list[[1]]] list[[2]]],list[[1]] Conjugate[list[[1]]]-list[[2]] Conjugate[list[[2]]]}];
  zero["sextic_Hopf_determinant_identity",cycleIm-Det[Transpose[hopf /@ vectors]]/4,group];
  scalarVector = Array[cr,3]+I Array[ci,3];
  oneRow = Table[ComplexExpand[Conjugate[scalarVector[[i]]] scalarVector[[j]]],{i,3},{j,3}];
  zero["sextic_one_row_identically_zero",fs/.Thread[Flatten[h]->Flatten[oneRow]],group];
  tangentVectors = {{1,0},{1,1},{eps (zr1+I zi1),eps (zr2+I zi2)}};
  ht = Table[ComplexExpand[Conjugate[tangentVectors[[i]]].tangentVectors[[j]]],{i,3},{j,3}];
  zero["sextic_nonvertex_quadratic_initial_term",(fs/.Thread[Flatten[h]->Flatten[ht]])-
    eps^2 (zi2 zr1-zr2 zi1),group];

  group = "determinant_nonvertex_cubic_initial_form";
  matrixY = Array[dy,{4,4}]; det4 = Det[DiagonalMatrix[{1,0,0,0}]+eps matrixY];
  det3 = Det[matrixY[[2;;4,2;;4]]];
  zero["det4_rank1_exact_expansion",det4-eps^3 det3-eps^4 Det[matrixY],group];
  zero["det4_rank1_lower_orders_vanish",Sum[Coefficient[Expand[det4],eps,j] eps^j,{j,0,2}],group];
  zero["det4_rank1_cubic_initial_form",Coefficient[Expand[det4],eps,3]-det3,group];

  (* Rectangular M_{p,3}(R) family: formal trace reductions for all p.
     The trace differential table is an input, proved in the manuscript. *)
  group = "rectangular_quartic_formal_trace_identities";
  rectE2 = (rectR^2-rectP2)/2;
  rectE3 = (rectP3-rectR rectP2+rectR rectE2)/3;
  rectP4 = rectR rectP3-rectE2 rectP2+rectE3 rectR;
  rectP5 = rectR rectP4-rectE2 rectP3+rectE3 rectP2;
  rectF = 2 rectP2-rectR^2;
  rectS = 64 rectP3-64 rectR rectP2+16 rectR^3;
  rectRFSrules = {rectR -> pr,rectP2 -> (pf+pr^2)/2,
    rectP3 -> (ps+32 pr pf+16 pr^3)/64};
  zero["rectangular_trace_P4",(rectP4/.rectRFSrules)-
    (pr ps/48+5 pr^2 pf/12+pr^4/8+pf^2/8),group];
  zero["rectangular_trace_P5",(rectP5/.rectRFSrules)-
    (5 pr^2 ps/256+5 pr^3 pf/16+pr^5/16+5 pr pf^2/24+5 pf ps/768),group];
  variables = {rectR,rectP2,rectP3};
  metric = {{4 rectR,8 rectP2,12 rectP3},
    {8 rectP2,16 rectP3,24 rectP4},{12 rectP3,24 rectP4,36 rectP5}};
  laplacians = {6 rows,4 (rows+4) rectR,6 ((rows+5) rectP2+rectR^2)};
  zero["rectangular_delta_R",laplace[rectR]-6 rows,group];
  zero["rectangular_delta_F",laplace[rectF]-4 (6-rows) rectR,group];
  zero["rectangular_delta_S",laplace[rectS]-32 (rows+6) rectR^2-448 rectF,group];
  zero["rectangular_gamma_RR",gamma[rectR,rectR]-4 rectR,group];
  zero["rectangular_gamma_RF",gamma[rectR,rectF]-8 rectF,group];
  zero["rectangular_gamma_RS",gamma[rectR,rectS]-12 rectS,group];
  zero["rectangular_gamma_FF",gamma[rectF,rectF]-rectS,group];
  zero["rectangular_gamma_FS",gamma[rectF,rectS]-
    8 rectR rectS-128 (rectR^2+rectF) rectF,group];
  zero["rectangular_gamma_SS",gamma[rectS,rectS]-
    192 (rectR^2+rectF) rectS-1024 rectR^3 rectF-2048 rectR rectF^2,group];
  rectDefect = defect[rectF];
  zero["rectangular_general_minimality_defect",rectDefect-
    4 (5-rows) rectR rectS+64 (rectR^2+rectF) rectF,group];
  zero["rectangular_p5_minimality_defect",(rectDefect/.rows->5)+
    64 (rectR^2+rectF) rectF,group];
  zero["rectangular_regular_point_F",rectF/.{rectR->6,rectP2->18,rectP3->66},group];
  zero["rectangular_regular_point_S",(rectS/.{rectR->6,rectP2->18,rectP3->66})-768,group];
  rectEigenvalues = {t1,t2,t3};
  rectEigenProfiles = {rectR,rectF,rectS}/.{rectR->Total[rectEigenvalues],
    rectP2->Total[rectEigenvalues^2],rectP3->Total[rectEigenvalues^3]};
  rectJacobian = Table[D[rectEigenProfiles[[i]],rectEigenvalues[[j]]],{i,3},{j,3}];
  zero["rectangular_RFS_eigenvalue_Jacobian_squared",Det[rectJacobian]^2-
    768^2 (t1-t2)^2 (t1-t3)^2 (t2-t3)^2,group];

  (* Independent finite coordinate tests: three explicitly fixed integer
     matrices in each M_{p,3}(R), p=3,4,5,6, plus one critical matrix.
     These tests are not a general coordinate proof in p. Each Hessian
     action is evaluated in matrix coordinates, without closure reduction. *)
  coordinateData[inputMatrix_List] := Module[
    {p=Length[inputMatrix],b,r,f,grad,hessianTimes,basis,columns,h,hg,
     gradS,s,deltaF,cc0,deltaS},
    b=Transpose[inputMatrix].inputMatrix;
    r=Total[Flatten[inputMatrix^2]];
    f=2 Tr[b.b]-r^2;
    grad=8 inputMatrix.b-4 r inputMatrix;
    hessianTimes[v_List] := 8 v.b+
      8 inputMatrix.(Transpose[v].inputMatrix+Transpose[inputMatrix].v)-
      8 Total[Flatten[inputMatrix v]] inputMatrix-4 r v;
    columns=Table[
      basis=Partition[UnitVector[3 p,j],3];
      Flatten[hessianTimes[basis]],{j,3 p}];
    h=Transpose[columns];
    hg=hessianTimes[grad];
    gradS=2 hg;
    s=Total[Flatten[grad^2]];
    deltaF=Tr[h]; cc0=4 (6-p);
    (* Bochner with Delta F = cc0 R; that identity is independently
       checked against the coordinate Hessian trace in every sample. *)
    deltaS=2 Total[Flatten[h^2]]+16 cc0 f;
    <|"R"->r,"F"->f,"S"->s,"DeltaF"->deltaF,"DeltaS"->deltaS,
      "GammaFS"->Total[Flatten[grad gradS]],
      "GammaSS"->Total[Flatten[gradS^2]],
      "Defect"->s deltaF-Total[Flatten[grad hg]],"Hessian"->h,"B"->b|>
  ];
  group="rectangular_quartic_finite_integer_coordinate_checks";
  Do[
    sampleMatrix=Table[Mod[17 i+11 j+3 i j+7 sample+sample i^2,5]-2,{i,p},{j,3}];
    sampleData=coordinateData[sampleMatrix];
    cr=sampleData["R"]; cf=sampleData["F"]; cs=sampleData["S"];
    sampleName="rectangular_coordinate_p"<>ToString[p]<>"_sample"<>ToString[sample]<>"_";
    zero[sampleName<>"delta_F",sampleData["DeltaF"]-4 (6-p) cr,group];
    zero[sampleName<>"delta_S",sampleData["DeltaS"]-32 (p+6) cr^2-448 cf,group];
    zero[sampleName<>"gamma_FS",sampleData["GammaFS"]-8 cr cs-128 (cr^2+cf) cf,group];
    zero[sampleName<>"gamma_SS",sampleData["GammaSS"]-
      192 (cr^2+cf) cs-1024 cr^3 cf-2048 cr cf^2,group];
    zero[sampleName<>"minimality_defect",sampleData["Defect"]-
      4 (5-p) cr cs+64 (cr^2+cf) cf,group];
    zero[sampleName<>"gradient_trace_identity",cs-
      64 Tr[MatrixPower[sampleData["B"],3]]+32 cr cf+16 cr^3,group],
    {p,{3,4,5,6}},{sample,1,3}];
  criticalMatrix=Table[If[i==j && i<=2,1,0],{i,5},{j,3}];
  criticalData=coordinateData[criticalMatrix];
  criticalR=criticalData["R"]; criticalH=criticalData["Hessian"];
  zero["rectangular_critical_RFS",(criticalR-2)^2+
    criticalData["F"]^2+criticalData["S"]^2,group];
  zero["rectangular_critical_Hessian_polynomial",
    Total[Flatten[(MatrixPower[criticalH,3]-
      4 criticalR MatrixPower[criticalH,2]-32 criticalR^2 criticalH)^2]],group];
  zero["rectangular_critical_Hessian_trace",Tr[criticalH]-4 criticalR,group];
  zero["rectangular_critical_Hessian_squared_trace",
    Tr[MatrixPower[criticalH,2]]-176 criticalR^2,group];

  (* Radial Clifford profile identities for two explicit homogeneous
     profiles. These are checks in fixed profile degree, not a proof of
     the general homogeneous-profile statement. *)
  group="radial_Clifford_profile_fixed_degree_checks";
  variables={rad,y1,y2};
  metric={{4 rad,4 y1,4 y2},{4 y1,4 rad,0},{4 y2,0,4 rad}};
  laplacians={2 ambientN,0,0};
  cubicProfile=rad^3+2 rad y1^2+y1 y2^2;
  profileK=D[cubicProfile,y1]^2+D[cubicProfile,y2]^2-D[cubicProfile,rad]^2;
  profileWave=D[cubicProfile,{y1,2}]+D[cubicProfile,{y2,2}]-D[cubicProfile,{rad,2}];
  profileLorentz=D[cubicProfile,y1] D[profileK,y1]+
    D[cubicProfile,y2] D[profileK,y2]-D[cubicProfile,rad] D[profileK,rad];
  zero["radial_profile_cubic_gradient_identity",gamma[cubicProfile,cubicProfile]-
    4 rad profileK-24 cubicProfile D[cubicProfile,rad],group];
  zero["radial_profile_cubic_laplacian_identity",laplace[cubicProfile]-
    4 rad profileWave-(2 ambientN+16) D[cubicProfile,rad],group];
  profileRHS=16 rad^2 (profileK profileWave-profileLorentz/2)+
    8 (ambientN-2) rad profileK D[cubicProfile,rad];
  zero["radial_profile_cubic_defect_modulo_profile",
    PolynomialRemainder[Expand[defect[cubicProfile]-profileRHS],cubicProfile,rad],group];
  quadraticProfile=y1^2+y2^2-2 rad^2/(ambientN-2);
  zero["radial_profile_quadratic_selected_level",
    PolynomialRemainder[Expand[defect[quadraticProfile]],quadraticProfile,y1],group];

  True, _, ($Failed &) ], $Aborted];
 EvaluationMessages=ToString[#,InputForm]& /@ System`$MessageList;
 EvaluationOK=TrueQ[evaluationResult===True && Length[EvaluationMessages]==0];
];
If[evaluationResult===$Aborted,Abort[]];
CertificateResult=<|"CertificateID"->"LSC-QUADRATIC-COMPOSITIONS",
  "Release"->"20261011_memo_v7","ExpectedChecks"->expectedChecks,
  "Checks"->checks,"EvaluationOK"->EvaluationOK,
  "EvaluationMessages"->EvaluationMessages,
  "Scope"->{
    "Formal orthogonal-quartic, rectangular-trace and complexified Hermitian-Gram identities are conditional on the differential tables proved in the manuscript.",
    "Quaternionic H^2 and complex C^2 coordinate identities are exact for the stated finite dimensions; the rectangular coordinate tests use exactly twelve integer matrices and one critical matrix.",
    "The radial Clifford profile checks use one degree-three and one degree-two profile; the general proof is in the manuscript.",
    "Positivity checks use the stated nonnegative shifted variables and the parameter ranges justified in the manuscript.",
    "The Hopf interpretation and nonvertex initial-form checks are supplementary algebraic identities; they do not assert additional geometric theorems in the manuscript.",
    "No check proves irreducibility, realizability, singular-set removal, global minimization, noncongruence or tangent-current convergence."
  }|>;
End[];
EndPackage[];
LaplacianCertificates`QuadraticCompositions`Private`CertificateResult
