(* Notebook interface 2026-10-02 r2. Original mathematical .wl files are unchanged.
   Get supplies the read-only $InputFileName; never assign or localize it. *)
BeginPackage["SupplementNotebook20261002`"];
RunSupplementNotebook::usage = "RunSupplementNotebook[directory, source, kind, expected] runs one supplemental unit and exports its execution record.";
Begin["`Private`"];

RunSupplementNotebook[directory_String, source_String, kind_String, expected_Integer] :=
Module[{path, prefix, stream, started, finished, sourceText, outcome,
        rows = {}, actual = 0, mathPassed = False, noMessages = True,
        exports, exportsOK, report, sourceHash, result, executionSucceeded,
        executionPath, temporarySource = False, sourceWritten,
        tag = "SupplementNotebook20261002-stop"},
  path = FileNameJoin[{directory, source}];
  prefix = FileNameJoin[{directory, "notebook_run_" <> ToString[expected]}];
  If[!FileExistsQ[path], Print["SOURCE FILE NOT FOUND: ", path]; Return[False]];
  If[!MemberQ[{"exceptional", "cubic", "pullback"}, kind], Return[False]];
  sourceText = Import[path, "Text"];
  If[!StringQ[sourceText], Print["Cannot read source."]; Return[False]];
  sourceHash = IntegerString[FileHash[path, "SHA256"], 16, 64];

  (* Only command-line termination is adapted. All assertion expressions stay intact. *)
  If[kind === "cubic",
    If[StringCount[sourceText, "Quit[1]"] =!= 1,
      Print["Unexpected source version: Quit adapter not applied."]; Return[False]];
    sourceText = StringReplace[sourceText,
      "Quit[1]" -> "Throw[$Failed, \"SupplementNotebook20261002-stop\"]"]
  ];
  If[kind === "pullback",
    If[StringCount[sourceText, "Exit[If[passed, 0, 1]];"] =!= 1,
      Print["Unexpected source version: Exit adapter not applied."]; Return[False]];
    sourceText = StringReplace[sourceText, "Exit[If[passed, 0, 1]];" -> "Null;"]
  ];
  executionPath = path;
  If[kind =!= "exceptional",
    (* Keep the temporary adapter beside the original so the source's own
       DirectoryName[$InputFileName] still selects the correct output folder. *)
    executionPath = FileNameJoin[{directory, "notebook_adapter_" <> CreateUUID[] <> ".wl"}];
    sourceWritten = Export[executionPath, sourceText, "Text"];
    If[!StringQ[sourceWritten] || !FileExistsQ[executionPath],
      Print["Cannot create notebook adapter source."]; Return[False]];
    temporarySource = True
  ];
  stream = OpenWrite[prefix <> "_output.txt", CharacterEncoding -> "UTF-8", PageWidth -> Infinity];
  If[Head[stream] =!= OutputStream,
    If[temporarySource, DeleteFile[executionPath]];
    Print["Cannot create output log."]; Return[False]];
  started = DateString[];
  Block[{$Output = Append[$Output, stream], $Messages = Append[$Messages, stream]},
    Print["Notebook interface: 20261002-notebook-r2; expected assertions: ", expected];
    Print["Kernel: ", $Version, "; system: ", $SystemID];
    Print["Started: ", started]; Print["Source: ", path];
    Print["Source SHA256: ", sourceHash];
    outcome = CheckAbort[
      Check[
        Catch[
          Block[{$Context = "Global`", $ContextPath = {"System`", "Global`"},
                 $IterationLimit = Infinity},
            Get[executionPath]
          ], tag],
        noMessages = False; $Failed],
      $Aborted];
    Switch[kind,
      "exceptional",
        If[ListQ[ExceptionalCubicVerification`log], rows = ExceptionalCubicVerification`log],
      "cubic",
        If[ListQ[Global`checkNames], rows = ({#, True, "0"} & /@ Global`checkNames)],
      "pullback",
        If[ListQ[Global`checks],
          rows = ({#["identity"], TrueQ[#["passed"]],
                    If[TrueQ[#["passed"]], "0", "See output"]} & /@ Global`checks)]
    ];
    actual = Length[rows];
    mathPassed = actual === expected && And @@ (TrueQ[#[[2]]] & /@ rows);
    executionSucceeded = outcome =!= $Failed && outcome =!= $Aborted;
    finished = DateString[];
    Print["Finished: ", finished];
    Print["Expected checks: ", expected, "; actual checks: ", actual];
    Print["Mathematical checks passed: ", mathPassed];
    Print["Evaluation succeeded: ", executionSucceeded];
    Print["No evaluation messages: ", noMessages];
  ];
  Close[stream];
  If[temporarySource, DeleteFile[executionPath]];
  report = <|"interface_version" -> "20261002-notebook-r2",
    "source_file" -> source, "source_sha256" -> sourceHash,
    "kernel_version" -> $Version, "system_id" -> $SystemID,
    "started" -> started, "finished" -> finished,
    "expected_checks" -> expected, "actual_checks" -> actual,
    "mathematical_checks_passed" -> mathPassed, "evaluation_succeeded" -> executionSucceeded,
    "no_evaluation_messages" -> noMessages,
    "scope" -> "Existing symbolic assertions; no new geometric proof."|>;
  exports = {
    Export[prefix <> "_results.csv", Prepend[rows, {"Check", "Passed", "Residual"}], "CSV"],
    Export[prefix <> "_summary.json", report, "RawJSON"]};
  exportsOK = And @@ (StringQ[#] && FileExistsQ[#] & /@ exports) && FileExistsQ[prefix <> "_output.txt"];
  result = TrueQ[mathPassed] && TrueQ[executionSucceeded] && TrueQ[noMessages] && TrueQ[exportsOK];
  Print["Record exports succeeded: ", exportsOK];
  Print["COMPLETE AND PASSED: ", result];
  Print["Records: ", prefix, "_output.txt / _results.csv / _summary.json"];
  result
];
End[];
EndPackage[];
