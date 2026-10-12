(* Release 20261011_memo_v7. Prepared 2026-10-11. *)
(* Evaluate Get[".../RunAll.wl"] in a fresh Mathematica 12.0 kernel. *)
Module[{root = DirectoryName[$InputFileName], runner},
  runner = FileNameJoin[{root, "Kernel", "Runner.wl"}];
  If[!FileExistsQ[runner],
    Return[Failure["RunnerMissing", <|"File" -> runner|>]]];
  If[Check[Get[runner], $Failed] === $Failed,
    Return[Failure["RunnerLoadFailed", <|"File" -> runner|>]]];
  LaplacianCertificates`Runner`RunSuite[]
]
