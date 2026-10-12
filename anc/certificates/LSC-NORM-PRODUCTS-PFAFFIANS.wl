(* ::Package:: *)
(*
   LSC-NORM-PRODUCTS-PFAFFIANS -- Mathematica 12 compatible.

   Coverage: 52 complete symbolic elementary-symmetric identities;
   5880 exact quaternion/octonion product and block derivative checks;
   83 exact complex Pfaffian checks, including three independent symbolic
   real-coordinate normalizations. Total: 6015 assertions.

   The fixed integer fixtures below have coordinate values in -2..2.
   Embedding fixtures avoids dependence on a particular random generator.

   Scope: the manuscript proves the arbitrary-degree statements. Finite
   exact coordinate checks support, and do not replace, those proofs.
   This source contains no claimed execution result. Get[thisFile] runs
   the checks and returns an association for the common runner.
*)

(* Release: 20261011_memo_v7. Prepared: 2026-10-11.
   Current statement numbers and scope are recorded in CertificateIndex.md.
   Execute RunAll.wl to generate fresh results. *)

BeginPackage["LaplacianCertificates`NormProductsPfaffians`"];
ClearAll["LaplacianCertificates`NormProductsPfaffians`*"];

result = Module[
  {check, elementary, gammaClosed, quaternion, conjugateVector, octonion,
   bracketedProduct, normedProductFixtures, pfaffianFixtures, records},

  check[name_String, actual_, expected_: 0] := Module[{residual, passed},
    residual = Expand[actual - expected];
    passed = TrueQ[residual === 0];
    Sow[<|"Name" -> name, "Passed" -> passed,
      "Detail" -> If[passed, "Exact residual is zero.",
        "Nonzero or unevaluated residual: " <> ToString[residual, InputForm]]|>,
      "CertificateCheck"];
    Null
  ];

  elementary[values_List, k_Integer] := If[k < 0 || k > Length[values], 0,
    Total[Times @@@ Subsets[values, {k}]]];

  gammaClosed[e_, i0_Integer, j0_Integer] := Module[
    {i = Min[i0, j0], j = Max[i0, j0]},
    4 Total[Table[(j - i + 1 + 2 l) e[i - 1 - l] e[j + l], {l, 0, i - 1}]]
  ];

  quaternion[a_List, b_List] := {
    a[[1]] b[[1]] - a[[2]] b[[2]] - a[[3]] b[[3]] - a[[4]] b[[4]],
    a[[1]] b[[2]] + a[[2]] b[[1]] + a[[3]] b[[4]] - a[[4]] b[[3]],
    a[[1]] b[[3]] - a[[2]] b[[4]] + a[[3]] b[[1]] + a[[4]] b[[2]],
    a[[1]] b[[4]] + a[[2]] b[[3]] - a[[3]] b[[2]] + a[[4]] b[[1]]
  };

  conjugateVector[a_List] := Prepend[-Rest[a], First[a]];

  octonion[x_List, y_List] := Module[{a, b, c, d},
    {a, b} = Partition[x, 4]; {c, d} = Partition[y, 4];
    Join[quaternion[a, c] - quaternion[conjugateVector[d], b],
      quaternion[d, a] + quaternion[b, conjugateVector[c]]]
  ];

  bracketedProduct[xs_List, mode_String, multiply_] := Module[{cut},
    If[Length[xs] == 1, Return[First[xs]]];
    cut = Switch[mode, "left", Length[xs] - 1, "right", 1,
      "balanced", Quotient[Length[xs], 2]];
    multiply[bracketedProduct[Take[xs, cut], mode, multiply],
      bracketedProduct[Drop[xs, cut], mode, multiply]]
  ];

  (* Each entry is {algebra dimension, degree, bracketing, sample index, xs}.
     The sample indices agree with the old record: 0 and 1. *)
  normedProductFixtures = {
    {4, 2, "left", 0, {{-1, -2, 0, 1}, {0, 2, -1, -1}}},
    {4, 2, "left", 1, {{1, -1, 1, 2}, {0, -1, -1, -2}}},
    {4, 2, "right", 0, {{-1, -1, 0, 1}, {0, -1, -1, -1}}},
    {4, 2, "right", 1, {{-2, -1, 1, 1}, {-1, 0, -2, -1}}},
    {4, 2, "balanced", 0, {{0, 1, 0, 0}, {0, 2, 0, -2}}},
    {4, 2, "balanced", 1, {{-2, 2, 0, 2}, {-1, -1, 0, 2}}},
    {4, 3, "left", 0, {{2, 1, -2, -1}, {-1, 0, 0, 1}, {2, 0, -2, -1}}},
    {4, 3, "left", 1, {{-2, 2, 2, 2}, {0, 2, 0, 2}, {0, 2, 1, 1}}},
    {4, 3, "right", 0, {{1, 2, 2, 0}, {1, -1, 2, 0}, {1, -2, 2, -2}}},
    {4, 3, "right", 1, {{1, 2, 0, 1}, {0, 2, -2, -1}, {-1, 2, 2, 0}}},
    {4, 3, "balanced", 0, {{2, 2, 2, 0}, {0, 0, 1, 0}, {-1, -1, 0, 1}}},
    {4, 3, "balanced", 1, {{1, 1, 2, -1}, {-1, -2, 1, -1}, {-2, -1, -2, 2}}},
    {4, 4, "left", 0, {{-1, 0, -2, -1}, {-2, 0, -1, -2}, {0, 1, 1, -1}, {-1, 0, 0, 2}}},
    {4, 4, "left", 1, {{-2, 0, -2, 0}, {-2, 2, 2, 2}, {-1, 2, -1, -2}, {1, -2, -2, 1}}},
    {4, 4, "right", 0, {{0, 0, 2, 0}, {2, 0, -1, 0}, {1, 2, 1, 1}, {1, 0, 1, -2}}},
    {4, 4, "right", 1, {{1, -2, 1, -2}, {-1, 0, -1, -1}, {-1, 2, -1, 0}, {1, 2, -2, 0}}},
    {4, 4, "balanced", 0, {{-2, 1, 1, 0}, {1, -1, -2, -1}, {-2, -2, 0, 0}, {-2, 2, 2, 2}}},
    {4, 4, "balanced", 1, {{0, 0, 1, -2}, {2, 0, 2, 2}, {-1, 1, 0, 1}, {0, -1, -1, -2}}},
    {4, 5, "left", 0, {{-2, 1, 1, -1}, {0, 1, -1, -2}, {2, 1, 0, -2}, {0, 1, 2, 2}, {0, -2, 0, -1}}},
    {4, 5, "left", 1, {{-1, -2, 2, 1}, {-1, -1, 2, -2}, {-2, -1, -2, 1}, {-2, 2, -1, 2}, {2, 0, 1, 0}}},
    {4, 5, "right", 0, {{2, 0, -1, -1}, {0, 0, -1, -1}, {0, 0, 2, 0}, {2, 1, -1, 1}, {-2, 2, 0, 2}}},
    {4, 5, "right", 1, {{2, 0, -2, 2}, {-2, 1, -1, 2}, {1, 1, -2, -2}, {1, 1, 2, 0}, {-2, 0, 0, 2}}},
    {4, 5, "balanced", 0, {{0, -2, 2, 0}, {0, -2, 0, 2}, {-1, -1, 1, 1}, {0, 1, -1, 0}, {-2, 2, -2, -1}}},
    {4, 5, "balanced", 1, {{-1, -2, 1, 0}, {-2, -2, 0, -1}, {0, 0, -1, 1}, {1, 2, 2, -1}, {2, 1, -1, 2}}},
    {4, 6, "left", 0, {{1, 0, -2, 2}, {2, 0, -1, 2}, {-1, -1, 0, 2}, {0, -1, -2, 0}, {-1, -1, 1, -2}, {-1, 2, 0, 1}}},
    {4, 6, "left", 1, {{1, 2, 1, 0}, {2, 2, -2, -2}, {0, 0, -2, 1}, {-2, -2, 1, 1}, {1, 2, 2, 2}, {1, 1, 1, 2}}},
    {4, 6, "right", 0, {{0, 1, -1, 2}, {2, 0, 2, -2}, {-1, -2, -2, -2}, {-1, 0, 1, -1}, {-1, 0, 1, -1}, {1, 1, 0, -1}}},
    {4, 6, "right", 1, {{0, 1, -2, 1}, {0, 0, -1, -2}, {-1, -1, 2, 1}, {1, 1, 2, -1}, {2, 2, -1, 0}, {0, 1, 2, -2}}},
    {4, 6, "balanced", 0, {{1, 1, 0, -1}, {1, 2, -2, 2}, {-1, -2, -2, -1}, {-1, 2, -2, 2}, {2, 1, 2, -1}, {2, 1, 0, 0}}},
    {4, 6, "balanced", 1, {{0, 0, 2, 0}, {0, -1, -1, 0}, {2, 0, -2, 1}, {2, -1, 2, 2}, {1, -2, 1, 1}, {-1, 2, 2, 1}}},
    {8, 2, "left", 0, {{0, 0, -1, 0, 2, 2, 1, -2}, {-1, 2, 0, 0, 0, -2, 0, 2}}},
    {8, 2, "left", 1, {{0, 1, 2, -2, -1, 2, 0, 2}, {1, 1, 0, 2, -2, -2, 1, -2}}},
    {8, 2, "right", 0, {{2, 1, 2, -1, 1, -1, 1, -2}, {2, -1, 1, 2, 0, -2, 1, -2}}},
    {8, 2, "right", 1, {{0, -2, 2, -2, 1, -2, -1, -2}, {2, -2, -1, 2, 1, 1, 1, 0}}},
    {8, 2, "balanced", 0, {{1, -2, -1, 2, 1, 0, -1, 2}, {2, 2, -2, 1, -1, 2, -1, 0}}},
    {8, 2, "balanced", 1, {{1, 0, 0, 2, -1, 2, -2, -2}, {-2, -1, 1, 2, -2, -2, -1, -1}}},
    {8, 3, "left", 0, {{-1, 2, 1, 2, 1, 1, 2, 0}, {-1, 2, -2, 2, 1, -1, 1, 2}, {2, 0, 1, -1, 1, 1, 0, -1}}},
    {8, 3, "left", 1, {{0, -2, 1, 1, -2, 1, 2, 1}, {-1, -2, 1, 1, 2, 0, -2, -1}, {-1, 0, -2, -1, -2, -1, 1, 0}}},
    {8, 3, "right", 0, {{2, 1, 1, 2, 2, 1, -2, 0}, {-1, -2, 1, 0, 2, 1, 1, 1}, {1, -2, 1, -2, 0, 1, -2, 0}}},
    {8, 3, "right", 1, {{2, -2, -2, 1, -1, -1, 0, -2}, {1, 2, 1, 2, 1, -1, 0, 2}, {1, 2, 1, 1, 2, -2, -2, -1}}},
    {8, 3, "balanced", 0, {{1, 0, 0, -2, -1, 1, 2, -2}, {0, -1, -1, -1, 2, 2, 1, -1}, {2, -2, -1, 1, -2, -1, 1, -1}}},
    {8, 3, "balanced", 1, {{2, 2, 1, -1, -2, 1, 2, -1}, {0, 0, -1, 2, -2, 2, 2, 0}, {0, 0, 0, 2, -2, -2, 2, 1}}},
    {8, 4, "left", 0, {{-1, -1, 1, -1, 2, 0, 2, 0}, {-2, -1, 0, 0, 2, 0, 2, -2}, {0, -2, 0, -2, -2, 2, 2, -1}, {0, 0, -1, -1, 2, 2, -1, -1}}},
    {8, 4, "left", 1, {{-2, -1, 2, -2, 0, -1, 2, 0}, {0, -1, 0, 2, 1, 2, -1, 2}, {-2, -2, 0, 0, 1, 2, 0, 1}, {0, 2, -2, 0, 1, 0, 0, -1}}},
    {8, 4, "right", 0, {{1, 1, 1, 0, -2, -1, -2, 0}, {-1, 0, -1, 2, 1, 1, -1, 0}, {-1, -1, 0, -1, 1, -2, -2, -1}, {-2, -2, -1, 2, 2, 2, -2, 1}}},
    {8, 4, "right", 1, {{1, -1, 0, 2, 2, 1, 2, 1}, {0, 0, 2, 2, 2, 0, -1, 1}, {-2, -2, -2, 2, 1, 2, -2, 1}, {0, -1, 1, -1, 1, -2, 2, -2}}},
    {8, 4, "balanced", 0, {{0, -2, 2, 0, 0, -1, 2, 0}, {2, 0, 1, -2, -2, 1, 2, 0}, {2, 1, 2, 2, -1, 0, 0, 2}, {1, -2, -1, 0, -1, 0, 0, 2}}},
    {8, 4, "balanced", 1, {{-2, 2, -2, 0, 2, -1, 1, 1}, {-2, -1, 0, -1, 0, -1, 2, 2}, {1, -1, 0, -1, 2, 0, 1, 2}, {1, 0, 2, 0, 2, 1, -1, -2}}},
    {8, 5, "left", 0, {{2, -1, -2, -2, -2, -2, -1, 1}, {2, 0, 2, 1, 0, -1, -2, 2}, {1, -2, 2, -2, -2, -2, 0, 1}, {-1, -2, -2, 2, -2, -1, 0, -1}, {-1, -2, 1, 1, -2, 0, 0, 2}}},
    {8, 5, "left", 1, {{2, -2, 1, 0, -1, -2, -2, -1}, {1, 2, 2, -2, 1, -1, -1, 2}, {0, 0, 1, 0, 0, -2, 1, -1}, {-2, -1, -2, -2, 2, 0, -2, -1}, {0, 2, -2, 1, 0, 1, 1, 1}}},
    {8, 5, "right", 0, {{1, -1, 2, -2, 2, 1, -1, -2}, {1, 2, 2, 0, 1, 0, 2, 2}, {2, 0, -2, 1, 2, -1, 1, 1}, {-2, -2, 0, 2, -1, -2, 1, 1}, {-1, 1, 1, 2, 2, -1, 1, 2}}},
    {8, 5, "right", 1, {{1, 1, 2, -1, -2, 1, 1, -2}, {-2, -1, 0, -2, 0, 1, 1, 2}, {-1, -1, 2, -1, 0, 0, -1, 2}, {1, 2, -2, -2, 2, -2, -2, -1}, {0, 0, 2, -2, -1, -2, -2, 0}}},
    {8, 5, "balanced", 0, {{0, 2, -1, -1, 2, 0, 0, -2}, {1, -2, -2, 2, -1, -2, -1, -2}, {0, 0, 2, -2, 2, 2, 1, 0}, {-1, 2, 1, 1, 1, 1, -2, -1}, {-1, -1, 0, -2, 0, -1, -2, 1}}},
    {8, 5, "balanced", 1, {{-2, -1, -1, 0, -1, -1, 1, 0}, {0, 2, 1, 1, -1, -2, 1, 0}, {0, -1, 0, -1, 2, 1, 1, 1}, {0, -2, 2, -2, 1, -2, 0, 1}, {1, 2, -1, 2, 1, 1, 2, -2}}},
    {8, 6, "left", 0, {{1, -1, 2, 0, -1, 0, -2, 0}, {-1, 2, 2, -1, 2, 2, 1, 1}, {1, -1, 1, 1, -1, 2, -1, 1}, {0, -1, 0, 0, 1, -2, 0, -2}, {2, 2, -2, 1, -2, 1, 1, 1}, {-1, 2, 0, 0, 2, -1, 1, -1}}},
    {8, 6, "left", 1, {{-2, -2, -2, 1, 0, 2, 1, 1}, {-1, 2, 2, 2, 2, 0, 2, 2}, {1, 2, -1, 0, 0, -1, 1, -1}, {0, 0, 2, 2, 0, 1, 2, -1}, {0, -1, -2, 2, -1, 1, 0, 2}, {2, 2, 2, 2, 1, 0, 1, -1}}},
    {8, 6, "right", 0, {{2, 2, 1, -1, 0, -2, 1, 1}, {-1, 0, 0, -1, 2, -1, 1, 1}, {1, 1, 2, 0, -2, 2, 0, -1}, {-1, -1, 0, 0, 0, 2, 2, 2}, {-2, 0, 2, 0, -2, 1, 2, 0}, {1, -1, 0, 1, 2, -1, -1, -2}}},
    {8, 6, "right", 1, {{2, -2, -2, -1, -1, 2, -2, 2}, {-1, -2, -2, -1, -2, -2, 0, 2}, {-2, 1, 0, -2, -1, -2, 1, 2}, {1, 1, 2, -2, -2, 0, -1, 2}, {1, -1, 0, 1, -1, -2, -1, -1}, {-1, 1, -2, 1, -1, -2, 1, 2}}},
    {8, 6, "balanced", 0, {{1, -1, 2, -1, -1, 2, -1, -1}, {2, 0, 0, 1, -2, -1, 2, 1}, {1, -2, -1, 2, 0, 1, -2, 0}, {-1, 0, -2, -2, 2, -1, 1, -2}, {-1, 1, -2, 0, -2, 2, 0, -1}, {-1, 2, 2, 0, 1, 1, -2, -2}}},
    {8, 6, "balanced", 1, {{0, 1, 2, 0, 2, 1, 2, 0}, {2, 2, -1, 0, 0, 2, 0, 0}, {1, -1, 2, -1, 1, 2, 0, -2}, {-2, -2, 0, -1, 1, 0, -1, -2}, {0, 2, 1, -2, -2, -2, 1, 0}, {2, -1, -2, 1, 1, 1, 1, -1}}}
  };

  (* Each entry is {Pfaffian degree, sample index, upper-triangular entries}.
     Coordinates are ordered by Subsets[Range[2 degree], {2}]. *)
  pfaffianFixtures = {
    {2, 0, {(0 + (2) I), (0 + (1) I), (-1 + (-2) I), (2 + (2) I), (-2 + (-2) I), (-2 + (0) I)}},
    {2, 1, {(2 + (0) I), (1 + (1) I), (2 + (1) I), (-1 + (-2) I), (-1 + (0) I), (1 + (0) I)}},
    {3, 0, {(0 + (1) I), (0 + (0) I), (-1 + (-2) I), (-1 + (-1) I), (0 + (0) I), (2 + (2) I), (1 + (-1) I), (-2 + (-1) I), (2 + (1) I), (2 + (-1) I), (1 + (2) I), (0 + (0) I), (-1 + (0) I), (-2 + (2) I), (-2 + (-2) I)}},
    {3, 1, {(1 + (-2) I), (-2 + (1) I), (-2 + (-1) I), (2 + (1) I), (2 + (-2) I), (-1 + (-2) I), (2 + (0) I), (0 + (2) I), (0 + (0) I), (0 + (-1) I), (-2 + (-2) I), (0 + (2) I), (0 + (-2) I), (0 + (-1) I), (0 + (-2) I)}},
    {4, 0, {(-1 + (1) I), (-1 + (-1) I), (0 + (0) I), (0 + (2) I), (0 + (2) I), (2 + (2) I), (-2 + (2) I), (-2 + (1) I), (-1 + (0) I), (0 + (0) I), (2 + (-1) I), (-1 + (2) I), (2 + (2) I), (1 + (1) I), (2 + (-2) I), (1 + (0) I), (-1 + (-2) I), (2 + (0) I), (-1 + (1) I), (1 + (0) I), (1 + (0) I), (0 + (-2) I), (-1 + (2) I), (0 + (-1) I), (-2 + (0) I), (-1 + (0) I), (1 + (0) I), (-2 + (1) I)}},
    {4, 1, {(1 + (0) I), (2 + (1) I), (-2 + (2) I), (-2 + (-2) I), (2 + (1) I), (-2 + (2) I), (0 + (1) I), (-2 + (1) I), (-2 + (-2) I), (0 + (1) I), (1 + (0) I), (0 + (-2) I), (1 + (-1) I), (2 + (1) I), (1 + (1) I), (2 + (1) I), (2 + (-1) I), (2 + (1) I), (1 + (2) I), (0 + (-2) I), (-2 + (-1) I), (-2 + (2) I), (0 + (-2) I), (-1 + (2) I), (1 + (0) I), (1 + (1) I), (-2 + (0) I), (0 + (0) I)}}
  };

  records = Last[Reap[

    (* Complete polynomial identities in formal squared-norm variables. *)
    Do[Module[{r, rr, es, e, direct},
      r = Array[rr, d];
      es = Table[elementary[r, k], {k, 0, d}];
      e[k_Integer] := If[0 <= k <= d, es[[k + 1]], 0];
      Do[
        check["symmetric Euler d=" <> ToString[d] <> ", k=" <> ToString[k],
          Total[(D[e[k], #] &) /@ r], (d - k + 1) e[k - 1]],
        {k, 1, d}];
      Do[
        direct = 4 Total[(# D[e[i], #] D[e[j], #] &) /@ r];
        check["symmetric Gamma d=" <> ToString[d] <> ", i=" <> ToString[i] <>
          ", j=" <> ToString[j], direct, gammaClosed[e, i, j]],
        {i, 1, d}, {j, i, d}];
      check["top generator recovery d=" <> ToString[d],
        gammaClosed[e, 2, d - 1],
        4 ((d - 2) e[1] e[d - 1] + d e[d])]
    ], {d, 2, 5}];

    (* Product norm, every upper-triangular entry of every block Gram
       matrix, each block Euler identity, and the scalar gradient norm. *)
    Do[Module[
      {delta, d, mode, sample, xs, multiply, units, rs, p, label,
       gradientNorm, columns, scale, perturbed},
      {delta, d, mode, sample, xs} = fixture;
      multiply = If[delta == 4, quaternion, octonion];
      units = IdentityMatrix[delta];
      rs = (Total[#^2] &) /@ xs;
      p = bracketedProduct[xs, mode, multiply];
      label = "K" <> ToString[delta] <> ",d" <> ToString[d] <> "," <>
        mode <> ",sample" <> ToString[sample];
      check["norm " <> label, Total[p^2], Times @@ rs];
      gradientNorm = 0;
      Do[
        columns = Table[
          perturbed = ReplacePart[xs, a -> units[[b]]];
          bracketedProduct[perturbed, mode, multiply], {b, 1, delta}];
        scale = Times @@ Delete[rs, a];
        Do[
          check["block conformality " <> label <> ",a" <> ToString[a - 1] <>
            "," <> ToString[i - 1] <> "," <> ToString[j - 1],
            columns[[i]].columns[[j]], If[i == j, scale, 0]],
          {i, 1, delta}, {j, i, delta}];
        gradientNorm += Total[columns[[All, 1]]^2];
        check["block Euler " <> label <> ",a" <> ToString[a - 1],
          xs[[a]].columns[[All, 1]], p[[1]]],
        {a, 1, d}];
      check["scalar gradient " <> label, gradientNorm, elementary[rs, d - 1]]
    ], {fixture, normedProductFixtures}];

    (* Full skew-matrix coordinate models of orders 4, 6 and 8. Pfaffian
       polynomials are recursively expanded before differentiating.
       For a holomorphic minor p, grad[k] stores sum dp Conjugate[p];
       real-coordinate Gamma is 4 Re[grad[k].Conjugate[grad[ell]]]. *)
    Do[Module[
      {n, pairs, zz, z, coordinateMatrix, pf, minors, derivatives, full,
       fixtures, at, evaluate, values, gradient, laps, e, fvalue, ff,
       df, v, dv, j, direct, sample, xx, yy, x, y, vars, realF},
      n = 2 d;
      pairs = Subsets[Range[n], {2}];
      zz = Array[z, Length[pairs]];
      coordinateMatrix = ConstantArray[0, {n, n}];
      Do[
        coordinateMatrix[[pairs[[a, 1]], pairs[[a, 2]]]] = zz[[a]];
        coordinateMatrix[[pairs[[a, 2]], pairs[[a, 1]]]] = -zz[[a]],
        {a, Length[pairs]}];
      pf[{}] = 1;
      pf[indices_List] := pf[indices] = Expand[Total[Table[
        (-1)^a coordinateMatrix[[First[indices], indices[[a]]]]
          pf[Delete[Rest[indices], a - 1]], {a, 2, Length[indices]}]]];
      minors = Table[pf /@ Subsets[Range[n], {2 k}], {k, 0, d}];
      derivatives = Table[
        Table[D[minor, #] & /@ zz, {minor, minors[[k + 1]]}], {k, 0, d}];
      full = pf[Range[n]];
      fixtures = Select[pfaffianFixtures, First[#] == d &];
      Do[
        sample = fixture[[2]];
        at = Thread[zz -> fixture[[3]]];
        evaluate[expr_] := Expand[expr /. at];
        values = ConstantArray[0, d + 1];
        gradient = ConstantArray[0, {d + 1, Length[zz]}];
        laps = ConstantArray[0, d + 1];
        Do[
          Do[
            v = evaluate[minors[[k + 1, a]]];
            dv = evaluate /@ derivatives[[k + 1, a]];
            values[[k + 1]] += v Conjugate[v];
            gradient[[k + 1]] += dv Conjugate[v];
            laps[[k + 1]] += 4 Total[dv Conjugate[dv]],
            {a, Length[minors[[k + 1]]]}];
          values[[k + 1]] = Expand[values[[k + 1]]];
          gradient[[k + 1]] = Expand /@ gradient[[k + 1]];
          laps[[k + 1]] = Expand[laps[[k + 1]]],
          {k, 0, d}];
        e[k_Integer] := If[0 <= k <= d, values[[k + 1]], 0];
        fvalue = evaluate[full]; ff = Re[fvalue];
        df = evaluate /@ (D[full, #] & /@ zz);
        check["Pf gradient normalization d=" <> ToString[d] <>
          ",sample=" <> ToString[sample], Total[df Conjugate[df]], e[d - 1]];
        Do[
          j = d - k + 1;
          check["Pf Laplacian d=" <> ToString[d] <> ",k=" <> ToString[k] <>
            ",sample=" <> ToString[sample], laps[[k + 1]],
            4 j (2 j - 1) e[k - 1]];
          check["Pf mixed Gamma d=" <> ToString[d] <> ",k=" <> ToString[k] <>
            ",sample=" <> ToString[sample],
            2 Re[df.Conjugate[gradient[[k + 1]]]], 2 j e[k - 1] ff];
          Do[
            direct = 4 Re[gradient[[k + 1]].Conjugate[gradient[[ell + 1]]]];
            check["Pf spectral Gamma d=" <> ToString[d] <> ",k=" <>
              ToString[k] <> ",ell=" <> ToString[ell] <> ",sample=" <>
              ToString[sample], direct, gammaClosed[e, k, ell]],
            {ell, k, d}],
          {k, 1, d}],
        {fixture, fixtures}];

      (* Independent real-coordinate check fixes the metric convention. *)
      If[d == 2,
        xx = Table[Unique["xr"], {Length[zz]}];
        yy = Table[Unique["yi"], {Length[zz]}];
        vars = Join[xx, yy];
        realF = Expand[ComplexExpand[Re[full /. Thread[zz -> xx + I yy]]]];
        check["Pf d2 real-coordinate gradient norm",
          Total[(D[realF, #]^2 &) /@ vars], Total[vars^2]];
        check["Pf d2 real-coordinate harmonicity",
          Total[(D[realF, {#, 2}] &) /@ vars]];
        check["Pf d2 real-coordinate ambient Laplacian",
          Total[(D[Total[vars^2], {#, 2}] &) /@ vars], 24]
      ]
    ], {d, 2, 4}];

    Null,
    "CertificateCheck"
  ]];
  records = If[records === {}, {}, First[records]];

  <|"CertificateID" -> "LSC-NORM-PRODUCTS-PFAFFIANS",
    "ExpectedChecks" -> 6015, "Checks" -> records,
    "Scope" -> <|
      "SymbolicClosureChecks" -> 52,
      "NormedProductChecks" -> 5880,
      "PfaffianChecks" -> 83,
      "ProductAlgebras" -> {"quaternions", "octonions"},
      "ProductDegrees" -> Range[2, 6],
      "Bracketings" -> {"left", "right", "balanced"},
      "SamplesPerCase" -> 2,
      "ComplexMatrixOrders" -> {4, 6, 8},
      "Limitation" -> "Finite exact checks do not replace the arbitrary-degree proofs."
    |>|>
];

EndPackage[];
LaplacianCertificates`NormProductsPfaffians`result
