(* ::Package:: *)
(* LSC-CUBIC-CERTIFICATES -- Cubic Certificates.
   Compatible with Mathematica 12. No external package or file dependencies.
   Get returns checks; only the central runner writes execution records.
   Scope: the nonharmonic checks verify the displayed algebraic defect.
   Square-freeness, singular-set removability, and geometric comparison
   are proved in the manuscript, not by these algebraic assertions.

   Exact identities conditional on displayed differential tables; not proofs of the tables, actual-image inequalities, spectral classification, or geometric comparison.
*)
(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`CubicCertificates`"];
Begin["`Private`"];
ClearAll["LaplacianCertificates`CubicCertificates`Private`*"];

VerificationLog = {};
EvaluationMessages = {};
EvaluationOK = Block[{$IterationLimit = Infinity, $Assumptions = True, System`$MessageList = {}},
  ExecutionStatus = CheckAbort[Catch[(
zero[name_String, expression_] := Module[{residual, passed},
  residual = Expand[Numerator[Together[expression]]];
  passed = TrueQ[residual === 0];
  AppendTo[VerificationLog,
    {name, passed, If[passed, "0", ToString[residual, InputForm]]}];
  passed
];
grad[expression_, variables_List] := (D[expression, #] & /@ variables);
gam[a_, b_, variables_List, table_List] :=
  grad[a, variables].table.grad[b, variables];
lap[expression_, variables_List, table_List, laps_List] :=
  grad[expression, variables].laps +
  Total[Flatten[Table[table[[i, j]] D[expression, variables[[i]], variables[[j]]],
      {i, Length[variables]}, {j, Length[variables]}]]];

(* Cubic certificate: G = F Sqrt[R S], with the stated sextic Gamma table. *)
Module[{r, f, s, dd, nn, k, lapd, variables, table, ell, qq, dl,
        qChain, dlChain, hh, jj, mOverG},
  variables = {r, f, s, dd};
  table = {{4 r, 6 f, 8 s, 12 dd}, {6 f, s, 4 r f, 2 s f},
           {8 s, 4 r f, 4 r s + 12 dd, 8 r dd},
           {12 dd, 2 s f, 8 r dd, 4 s dd}};
  ell = {1/(2 r), 1/f, 1/(2 s), 0};
  qq = r s^2 + f^2 (3 dd r/s + 5 r^2 + 11 s);
  dl = (nn + 9)/r + k/(2 s) + 3 r/s - 3 dd/s^2;
  qChain = r s f^2 (ell.table.ell);
  dlChain = ell.{2 nn, 0, k, lapd} +
    Total[Flatten[Table[table[[i, j]] (D[ell[[i]], variables[[j]]] + ell[[i]] ell[[j]]),
        {i, 4}, {j, 4}]]];
  zero["unified gradient norm", qChain - qq];
  zero["unified Laplacian ratio", dlChain - dl];
  hh = (nn - 12) s^2 + r (k - 16 r) s/2 - 12 r dd;
  jj = r dd (3 r (k - 16 r)/2 + 15 (r^2 - 3 s) + 3 (nn - 8) s) +
       5 (k - 16 r) r^3 s/2 + 45 r^4 s + 11 (k - 16 r) r s^2/2 +
       (5 nn + 33) r^2 s^2 + (11 nn + 44) s^3;
  mOverG = qq dl - ell.table.grad[qq, variables]/2;
  zero["unified H/J defect identity", r s^2 mOverG - r s^2 hh - f^2 jj];
  zero["trace corollary remainder",
       hh - (nn - 16) s^2 - r (k - 16 r) s/2 - 4 (s^2 - 3 r dd)]
];

(* Cubic Certificates identities. *)
Module[{u, v, tt, ff, m, q, variables, table, r, f, s, dd, k,
        hh, h44, lam, al, ga, sl, jl, hl, boundary},
  variables = {u, v, tt, ff};
  table = {{4 u, 0, 8 tt, 4 ff}, {0, 4 v, 0, 2 ff},
           {8 tt, 0, 16 u tt, 8 u ff},
           {4 ff, 2 ff, 8 u ff, 4 u v + tt}};
  r = u + v; f = ff/2; s = u v + tt/4; dd = v tt/4;
  k = 4 m v + 4 (q + 1) u;
  zero["Clifford Gamma(f,s)", gam[f, s, variables, table] - 4 r f];
  zero["Clifford Gamma(s,s)", gam[s, s, variables, table] - 4 r s - 12 dd];
  zero["Clifford Gamma(f,D)", gam[f, dd, variables, table] - 2 s f];
  zero["Clifford Gamma(s,D)", gam[s, dd, variables, table] - 8 r dd];
  zero["Clifford Gamma(D,D)", gam[dd, dd, variables, table] - 4 s dd];
  zero["Clifford normalized Laplacian",
       lap[s, variables, table, {4 m, 2 (q + 1), 8 (q + 1) u, 0}] - k];
  hh = (2 m + q - 11) s^2 + r (k - 16 r) s/2 - 12 r dd;
  h44 = tt (8 u^2 + tt)/16 + (u^2 - tt) v (2 u + 3 v);
  zero["Clifford H44 decomposition", (hh /. {m -> 4, q -> 4}) - h44];
  zero["Clifford parameter increment",
       hh - h44 - (2 (m - 4) + q - 4) s^2 - 2 r ((m - 4) v + (q - 4) u) s];
  sl = lam (4 - 3 lam);
  jl = 16 (4 lam (1 - lam) + (4 - 3 lam) lam^2);
  hl = 32 lam (1 - lam)^2 + 32 lam^2 (1 - lam) + 8 lam^3 - 16 (1 - lam) lam^2;
  boundary = 64 ga - 8 - ga (ga + 2) jl/sl - ga hl/sl + 2 al (5 - 8 ga - 2 al) sl;
  zero["Clifford power obstruction endpoint", Limit[boundary, lam -> 0] + 8 (ga - 1) (2 ga - 1)];
  zero["Clifford power obstruction interior",
       (boundary /. lam -> 2/3) + 8 (2 al^2 + (8 ga - 5) al + 8 ga^2 - 6 ga + 3)/3]
];

(* General radial multiplication identity; h has homogeneous degree dd. *)
Module[{r, h, w, ell, z, dd, nn, variables, table, lg, qg, dl, actual, expected},
  variables = {r, h, w};
  table = {{4 r, 2 dd h, 4 (dd - 1) w},
           {2 dd h, w, z}, {4 (dd - 1) w, z, 0}};
  lg = {1/(2 r), 1/h, 0};
  qg = r w + (2 dd + 1) h^2;
  zero["radial lift gradient identity", r h^2 (lg.table.lg) - qg];
  dl = ell/h + (2 dd + nn - 1)/r;
  actual = qg dl - lg.table.grad[qg, variables]/2;
  expected = r/h (w ell - z/2 + (2 dd + 1) h^2 ell/r +
             (nn - 2 dd - 2) h w/r + (2 dd + 1) (dd + nn - 1) h^3/r^2);
  zero["radial lift defect identity", actual - expected];
  zero["dimension 26 radial lift coefficient 11", (2 dd + 1 /. dd -> 5) - 11];
  zero["dimension 26 radial lift coefficient 14", (nn - 2 dd - 2 /. {nn -> 26, dd -> 5}) - 14];
  zero["dimension 26 radial lift coefficient 330", ((2 dd + 1) (dd + nn - 1) /. {nn -> 26, dd -> 5}) - 330]
];

(* Clifford cubic certificate: the five-invariant closure at d = 4, ambient dimension 26.
   w denotes 4 u* Z^2 u; p is the gradient norm square. *)
Module[{rho, nu, a, b, w, variables, table, laps, f, p, q, k, tt, logh,
        qh, ll, polynomialN, mOverH, pp, ff, nFormal, pmax, h0, k0, terms, lower},
  variables = {rho, nu, a, b, w};
  table = {{4 rho, 0, 6 a, 2 b, 4 w},
           {0, 4 nu, 0, 4 b, 4 w},
           {6 a, 0, 4 rho^2/3, 4 rho nu/3 - w/2, 16 nu a - 16 rho b/3},
           {2 b, 4 b, 4 rho nu/3 - w/2, 4 nu^2/3 + w, 16 rho b - 16 nu a + 32 nu b/3},
           {4 w, 4 w, 16 nu a - 16 rho b/3, 16 rho b - 16 nu a + 32 nu b/3,
            16 (rho + nu) w + 64 b^2/3 - 64 a b}};
  laps = {28, 24, 0, 0, 16 (12 rho + 14 nu)/3};
  f = a - b/3; p = 4 (rho - nu/3)^2/3 + 4 w/9;
  q = 8 (rho + 7 nu/27); k = 32 (78 rho + 14 nu)/27;
  tt = 64 (9 rho + nu) p/27 - 1024 f b/81;
  zero["nonradial Gamma(F,F)", gam[f, f, variables, table] - p];
  zero["nonradial harmonicity", lap[f, variables, table, laps]];
  zero["nonradial minimality", gam[f, p, variables, table] - 2 q f];
  zero["nonradial Delta(P)", lap[p, variables, table, laps] - k];
  zero["nonradial Gamma(P,P)", gam[p, p, variables, table] - tt];
  logh = grad[f, variables]/f + grad[p, variables]/(2 p);
  qh = p^2 + 2 q f^2 + f^2 tt/(4 p);
  ll = 512 rho/9 + 320 nu/27 + 256 f b/(81 p);
  zero["nonradial H gradient", f^2 p (logh.table.logh) - qh];
  zero["nonradial H Laplacian", 2 q + k/2 - tt/(4 p) - ll];
  polynomialN = (405 rho + 81 nu) p^3 + 144 f b p^2 +
    f^2 p (46692 rho^2 + 20040 rho nu + 2212 nu^2 - 4239 p) -
    (4992 rho + 896 nu) f^3 b - 256 nu f^4;
  mOverH = qh ll/p - logh.table.grad[qh, variables]/2;
  zero["nonradial exact N defect", 729 p^2 mOverH - 16 polynomialN];
  nFormal = (405 rho + 81 nu) pp^3 + 144 ff b pp^2 +
    ff^2 pp (46692 rho^2 + 20040 rho nu + 2212 nu^2 - 4239 pp) -
    (4992 rho + 896 nu) ff^3 b - 256 nu ff^4;
  pmax = 4 rho^2/3 + 40 rho nu/27 + 4 nu^2/27;
  h0 = 8 (36639 rho^2 + 12634 rho nu + 1555 nu^2)/9;
  k0 = 4992 rho + 896 nu;
  terms = {405 rho pp^3, 36 pp^2 (ff + 2 b)^2,
           pp^2 (81 nu pp - 144 b^2), ff^2 pp h0,
           k0 ff^2 (3 (9 rho + nu) pp/16 - ff b),
           256 nu ff^2 ((rho + nu) pp/9 - ff^2), 4275 ff^2 pp (pmax - pp)};
  zero["nonradial seven-summand decomposition", nFormal - Total[terms]];
  zero["nonradial rho=0 strictness reduction", (nFormal /. {rho -> 0, ff -> 0, b -> 0}) - 81 nu pp^3];
  zero["nonradial P maximum remainder", pmax - p - 4 (16 rho nu/3 - w)/9];
  zero["nonradial square bound for H Delta(H)", (rho + 2 nu)^2 - 4 nu (rho + nu) - rho^2];
  lower = 512 rho/9 + 320 nu/27 - 32 (rho + 2 nu)/81;
  zero["nonradial H Delta(H) lower coefficients", lower - (4576 rho + 896 nu)/81]
];

(* Cubic Certificates identities.
   H2 = |Hess F|^2, K = Delta S, T = Gamma(S,S).
   Gamma(F,t) = K/2 - H2 is the full Bochner identity, not harmonicity. *)
Module[{r, f, s, t, q, k, tt, h2, nn, gft, gfq, gftt, gst, gsq, gstt,
        variables, fp, qh, dl, b0, b1, b2, bb, mOverH, c0, c1, c2,
        fromLift, hgen, dd, hrad, gammaFQ, gammaSQ},
  variables = {f, s, t, q, tt}; fp = 2 s t + 2 q f;
  qh = s^2 + 2 f s t + f^2 (2 q + tt/(4 s));
  dl = 3 t/f + (2 q + k/2 - tt/(4 s))/s;
  zero["nonharmonic H gradient chain rule", s (s + f fp/s + f^2 tt/(4 s^2)) - qh];
  zero["nonharmonic H Laplacian chain rule", t/f + k/(2 s) + fp/(f s) - tt/(4 s^2) - dl];
  b0 = 16 s^2 (s (h2 - 2 q + 3 t^2) - tt);
  b1 = -16 s^2 gfq - 8 s^2 gst + 16 k s^2 t + 80 s^2 q t - 4 s tt t - 2 s gftt;
  b2 = -8 s gsq + 16 k q s + 2 k tt + 32 q^2 s - gstt;
  bb = b0 + f b1 + f^2 b2;
  gammaFQ = grad[qh, variables].{s, fp, gft, gfq, gftt};
  gammaSQ = grad[qh, variables].{fp, tt, gst, gsq, gstt};
  mOverH = qh dl - gammaFQ/(2 f) - gammaSQ/(4 s);
  zero["nonharmonic B identity with Bochner correction", (16 s^2 mOverH - bb) /. gft -> k/2 - h2];
  c0 = r^2 b0 + 16 (nn - 12) r s^4;
  c1 = r^2 b1 + 16 (2 nn + 9) r s^3 t;
  c2 = r^2 b2 + r (88 k s^2 + 32 (nn - 1) q s^2 + 4 (nn - 23) s tt) + 176 (nn + 4) s^3;
  fromLift = r^2 bb + 176 r (3 f s^3 t + f^2 s^2 (2 q + k/2 - tt/(4 s))) +
             16 r (nn - 12) s^2 qh + 176 (nn + 4) f^2 s^3;
  zero["nonharmonic C identity", fromLift - c0 - f c1 - f^2 c2];
  hgen = (nn - 12) s^2 + r s (h2 - 2 q + 3 t^2) - r tt;
  zero["nonharmonic leading coefficient", c0 - 16 r s^2 hgen];
  hrad = (nn - 12) s^2 + r (k - 16 r) s/2 - 12 r dd;
  zero["nonharmonic to radial leading term", (hgen /. {t -> 0, q -> 2 r, tt -> 4 r s + 12 dd, h2 -> k/2}) - hrad]
];


    ); True, _, ($Failed &)], $Aborted];
  EvaluationMessages = ToString[#, InputForm] & /@ System`$MessageList;
  ExecutionStatus === True && Length[EvaluationMessages] == 0
];

CertificateResult = <|"CertificateID" -> "LSC-CUBIC-CERTIFICATES",
  "ExpectedChecks" -> 38,
  "Checks" -> (<|"Name" -> #[[1]], "Passed" -> TrueQ[#[[2]]],
      "Detail" -> If[StringQ[#[[3]]], #[[3]], ToString[#[[3]], InputForm]]|> &) /@ VerificationLog,
  "EvaluationOK" -> TrueQ[EvaluationOK],
  "EvaluationMessages" -> EvaluationMessages|>;
End[];
EndPackage[];
LaplacianCertificates`CubicCertificates`Private`CertificateResult
