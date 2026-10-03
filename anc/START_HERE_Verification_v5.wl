(* Verification v6: stable v5-named root entry point; filenames and contexts are stable unit IDs. *)
Block[{$IterationLimit = Infinity, $ContextPath = {"System`"}, $Context = "Global`"},
 Get[FileNameJoin[{DirectoryName[$InputFileName], "mathematica", "chapters",
   "Run_All_Verifications.wl"}]]
]
