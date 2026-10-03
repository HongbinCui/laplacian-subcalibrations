(* Complete verification v6, 20260923. Stable v5 filenames/contexts are unit IDs.
   One lexical scope for the runner.
   Per-chapter mathematics is unchanged. No package contexts are exported. *)
Block[{$IterationLimit = Infinity},
Module[{baseDirectory, units, summary = {}, started, release = "20260923-v6",
        unit, overall, totalCount, outputDirectory, lines, outputFile},
  ClearAll["MergedVerificationV5`RunAll`*"];
  baseDirectory = If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    DirectoryName[$InputFileName], Directory[]];
  units = {
    {"Verify_Chapter_04.wl", "MergedVerificationV5`Chapter04`", 6},
    {"Verify_Chapter_05.wl", "MergedVerificationV5`Chapter05`", 16},
    {"Verify_Chapter_06.wl", "MergedVerificationV5`Chapter06`", 49},
    {"Verify_Chapter_07.wl", "MergedVerificationV5`Chapter07`", 18},
    {"Verify_Chapter_08.wl", "MergedVerificationV5`Chapter08`", 193},
    {"Verify_Chapter_09.wl", "MergedVerificationV5`Chapter09`", 50},
    {"Verify_Chapter_10.wl", "MergedVerificationV5`Chapter10`", 134},
    {"Verify_Appendix_A.wl", "MergedVerificationV5`AppendixA`", 8},
    {"Verify_Appendix_B.wl", "MergedVerificationV5`AppendixB`", 15},
    {"Verify_Appendix_C.wl", "MergedVerificationV5`AppendixC`", 1},
    {"Verify_Appendix_D.wl", "MergedVerificationV5`AppendixD`", 15},
    {"Verify_Appendix_E.wl", "MergedVerificationV5`AppendixE`", 18},
    {"Verify_Appendix_F.wl", "MergedVerificationV5`AppendixF`", 27}
  };
  started = DateString[];
  Print["COMPLETE VERIFICATION V6 (20260923) -- ", $Version];
  Print["Source directory: ", baseDirectory];
  Block[{$IterationLimit = Infinity}, Do[
    Module[{file = unit[[1]], contextName = unit[[2]], expected = unit[[3]],
            value, log, passed},
      Print["\n========== ", file, " =========="];
      value = Block[{$ContextPath = {"System`"}, $Context = "Global`"},
        CheckAbort[Get[FileNameJoin[{baseDirectory, file}]], $Aborted]];
      If[value === $Aborted, Abort[]];
      log = ToExpression[contextName <> "VerificationLog"];
      passed = TrueQ[value] && ListQ[log] && Length[log] == expected;
      AppendTo[summary, {file, expected, If[ListQ[log], Length[log], 0], passed}]
    ], {unit, units}]];
  overall = Length[summary] == Length[units] && And@@(TrueQ[#[[4]]] & /@ summary);
  totalCount = Total[summary[[All,3]]];
  outputDirectory = FileNameJoin[{baseDirectory, "verification_results_v6"}];
  If[!DirectoryQ[outputDirectory],
    CreateDirectory[outputDirectory, CreateIntermediateDirectories -> True]];
  lines = Join[
    {"Revised manuscript: actual v6 execution of the retained 13-unit suite",
     "Code release: " <> release,
     "Kernel: " <> $Version, "System ID: " <> $SystemID,
     "Started: " <> started, "Finished: " <> DateString[],
     "Overall passed: " <> ToString[overall, InputForm],
     "Total retained-suite assertions, including stated appendix repetitions: " <>
       ToString[totalCount, InputForm], "",
     "File | expected | actual | passed"},
    (StringRiffle[ToString[#, InputForm] & /@ #, " | "] &) /@ summary];
  outputFile = Export[FileNameJoin[{outputDirectory, "Merged_Verification_Summary.txt"}],
    StringRiffle[lines,"\n"] <> "\n", "Text"];
  If[!StringQ[outputFile] || !FileExistsQ[outputFile], overall = False];
  MergedVerificationV5`RunAll`VerificationSummary = summary;
  MergedVerificationV5`RunAll`VerificationTotalCount = totalCount;
  MergedVerificationV5`RunAll`VerificationResults = overall;
  MergedVerificationV5`RunAll`VerificationOutputDirectory = outputDirectory;
  Print["\nTotal assertions: ", totalCount];
  If[overall, Print["ALL CHAPTERS AND APPENDICES PASSED."],
    Print["AT LEAST ONE UNIT FAILED, WAS INCOMPLETE, OR COULD NOT BE SAVED."]];
  overall
]
]
