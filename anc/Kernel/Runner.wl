(* ::Package:: *)
(* Release 20261011_memo_v7. Prepared 2026-10-11.
   Native Mathematica 12.0 compatible. Each execution creates a new run folder.
   Every reported result is evaluated in this run; no prior result is read. *)
BeginPackage["LaplacianCertificates`Runner`"];
RunSuite::usage = "RunSuite[] evaluates the full manifest. Options: CertificateIDs -> All and ResultsDirectory -> Automatic (both option names are strings).";
Begin["`Private`"];

packageRoot = DirectoryName[DirectoryName[$InputFileName]];
Options[RunSuite] = {"CertificateIDs" -> All, "ResultsDirectory" -> Automatic};

utcNow[] := DateString[Now,
  {"ISODate", "T", "Hour", ":", "Minute", ":", "Second", "Z"}, TimeZone -> 0];
localNow[] := DateString[Now, {"ISODate", " ", "Time"}];
sourcePath[path_String] := FileNameJoin[Prepend[StringSplit[path, "/"], packageRoot]];
sha256[file_String] := If[FileExistsQ[file],
  Quiet[Check[IntegerString[FileHash[file, "SHA256"], 16, 64], "HASH_ERROR"]],
  "MISSING"];
snapshot[paths_List] := Association @@ ((# -> sha256[sourcePath[#]]) & /@ paths);

invalidateOutput[directory_String, failed_String] := (
  Quiet[Check[Export[FileNameJoin[{directory, "status.json"}],
    <|"Status" -> "INCOMPLETE", "FailedOutput" -> failed,
      "FinishedUTC" -> utcNow[]|>, "RawJSON"], Null]];
  Throw[Failure["OutputWriteFailed", <|"File" -> failed|>], "RunFailure"]
);
safeExport[file_String, value_, format_String] := Module[{answer},
  answer = Quiet[Check[Export[file, value, format], $Failed]];
  If[!StringQ[answer] || !FileExistsQ[file],
    invalidateOutput[DirectoryName[file], FileNameTake[file]]];
  answer
];
safePut[file_String, value_] := Module[{answer},
  answer = Quiet[Check[Put[value, file]; True, False]];
  If[!TrueQ[answer] || !FileExistsQ[file],
    invalidateOutput[DirectoryName[file], FileNameTake[file]]];
  file
];
appendText[file_String, text_String] := Module[{stream, answer},
  stream = Quiet[Check[OpenAppend[file, CharacterEncoding -> "UTF8"], $Failed]];
  If[Head[stream] =!= OutputStream,
    invalidateOutput[DirectoryName[file], FileNameTake[file]]];
  answer = Quiet[Check[WriteString[stream, text]; Close[stream]; True, False]];
  If[!TrueQ[answer],
    Quiet[Close[stream]];
    invalidateOutput[DirectoryName[file], FileNameTake[file]]];
  Null
];
progress[directory_String, message_String] := (
  Print[message];
  appendText[FileNameJoin[{directory, "progress.log"}], utcNow[] <> "  " <> message <> "\n"]
);
validCheck[row_] := AssociationQ[row] &&
  StringQ[Lookup[row, "Name", Null]] && KeyExistsQ[row, "Passed"] &&
  StringQ[Lookup[row, "Detail", Null]];
normalizeCheck[row_Association, id_String, index_Integer] := <|
  "CheckID" -> id <> "::" <> IntegerString[index, 10, 5],
  "Name" -> row["Name"], "Passed" -> TrueQ[row["Passed"]],
  "Detail" -> row["Detail"]|>;

RunSuite[OptionsPattern[]] := Catch[Module[
  {manifest, fullUnits, units, chosen, base, output, runName, paths, before, after,
   begun, begunLocal, begunClock, reports = {}, rows, id, file, raw, elapsed,
   messages, evaluationOK, countOK, catalogOK, complete, selectedPassed,
   summary, totals, report, status, expectedHashes, actualIDs, diagnostic,
   requiredPaths, aborted = False, kernelInfo, progressStatus, sourceMismatch,
   allRowsValid, checkLines, release, finished},
  If[!NumberQ[$VersionNumber] || $VersionNumber < 12 || !StringQ[$Version] ||
    StringContainsQ[$Version, "Mathics"],
    Return[Failure["NativeKernelRequired", <|"Message" -> "Use native Mathematica 12.0 or later."|>]]];
  manifest = Quiet[Check[Import[FileNameJoin[{packageRoot, "manifest.json"}], "RawJSON"], $Failed]];
  If[!AssociationQ[manifest], Return[Failure["ManifestUnreadable", <||>]]];
  fullUnits = Lookup[manifest, "Certificates", {}];
  If[!ListQ[fullUnits] || fullUnits === {} || !And @@ (AssociationQ /@ fullUnits) ||
    !And @@ (Function[u,
      StringQ[Lookup[u, "CertificateID", Null]] &&
      StringQ[Lookup[u, "File", Null]] &&
      IntegerQ[Lookup[u, "ExpectedChecks", Null]] && Lookup[u, "ExpectedChecks", 0] > 0 &&
      ListQ[Lookup[u, "CheckNames", Null]] &&
      Length[Lookup[u, "CheckNames", {}]] === Lookup[u, "ExpectedChecks", -1] &&
      And @@ (StringQ /@ Lookup[u, "CheckNames", {}]) &&
      Lookup[u, "File", ""] === "certificates/" <> Lookup[u, "CertificateID", ""] <> ".wl"
    ] /@ fullUnits) ||
    !DuplicateFreeQ[Lookup[fullUnits, "CertificateID"]] ||
    Length[fullUnits] =!= Lookup[manifest, "ExpectedCertificateCount", -1] ||
    Total[Lookup[fullUnits, "ExpectedChecks"]] =!= Lookup[manifest, "ExpectedChecks", -1],
    Return[Failure["InvalidCertificateManifest", <||>]]];
  expectedHashes = Lookup[manifest, "SourceSHA256", <||>];
  requiredPaths = Join[{"RunAll.wl", "Kernel/Runner.wl"}, Lookup[fullUnits, "File"]];
  If[!AssociationQ[expectedHashes] ||
    !And @@ (KeyExistsQ[expectedHashes, #] & /@ requiredPaths) ||
    !And @@ (StringQ[#] && StringMatchQ[#, RegularExpression["[0-9a-f]{64}"]] & /@ Values[expectedHashes]),
    Return[Failure["InvalidSourceManifest", <||>]]];
  paths = Keys[expectedHashes];
  before = snapshot[paths];
  sourceMismatch = Select[paths, before[#] =!= expectedHashes[#] &];
  If[sourceMismatch =!= {},
    Return[Failure["SourceHashMismatch", <|"Files" -> sourceMismatch,
      "Message" -> "Use the unchanged prepared files. Source edits require a new release manifest."|>]]];
  chosen = OptionValue["CertificateIDs"];
  units = fullUnits;
  If[chosen =!= All,
    If[!ListQ[chosen] || chosen === {} || !DuplicateFreeQ[chosen] ||
      !And @@ (MemberQ[Lookup[fullUnits, "CertificateID"], #] & /@ chosen),
      Return[Failure["UnknownCertificateIDs", <|"Requested" -> ToString[chosen, InputForm]|>]]];
    units = Select[fullUnits, MemberQ[chosen, #["CertificateID"]] &]];
  base = Replace[OptionValue["ResultsDirectory"],
    Automatic :> FileNameJoin[{packageRoot, "results"}]];
  If[!StringQ[base], Return[Failure["InvalidResultsDirectory", <||>]]];
  runName = "run-" <> DateString[Now,
    {"Year", "Month", "Day", "T", "Hour", "Minute", "Second", "Z"}, TimeZone -> 0] <>
    "-" <> StringTake[CreateUUID[], 8];
  output = FileNameJoin[{ExpandFileName[base], runName}];
  If[DirectoryQ[output], Return[Failure["RunDirectoryAlreadyExists", <|"Directory" -> output|>]]];
  If[Quiet[Check[CreateDirectory[output, CreateIntermediateDirectories -> True], $Failed]] === $Failed,
    Return[Failure["RunDirectoryCreationFailed", <|"Directory" -> output|>]]];
  begun = utcNow[]; begunLocal = localNow[]; begunClock = AbsoluteTime[];
  release = Lookup[manifest, "Release", "20261011_memo_v7"];
  kernelInfo = <|"WolframVersion" -> $Version, "VersionNumber" -> $VersionNumber,
    "SystemID" -> $SystemID, "OperatingSystem" -> $OperatingSystem,
    "MachineType" -> $MachineType, "ProcessorType" -> $ProcessorType,
    "SystemWordLength" -> $SystemWordLength, "MachinePrecision" -> $MachinePrecision,
    "TimeZoneHours" -> $TimeZone|>;
  progressStatus = <|"Status" -> "RUNNING", "Release" -> release,
    "StartedUTC" -> begun, "StartedLocal" -> begunLocal,
    "ExpectedUnits" -> Length[units], "CompletedUnits" -> 0,
    "ExpectedChecks" -> Total[Lookup[units, "ExpectedChecks"]], "Kernel" -> kernelInfo|>;
  safeExport[FileNameJoin[{output, "status.json"}], progressStatus, "RawJSON"];
  safeExport[FileNameJoin[{output, "source_snapshot.json"}], before, "RawJSON"];
  safeExport[FileNameJoin[{output, "manifest_snapshot.json"}], manifest, "RawJSON"];
  safeExport[FileNameJoin[{output, "checks.log"}],
    "CheckID\tResult\tName\tDetail\n", "Text"];
  progress[output, "Release " <> release <> ": " <> ToString[Length[units]] <>
    " units; " <> ToString[Total[Lookup[units, "ExpectedChecks"]]] <> " expected checks."];
  progress[output, "Kernel: " <> $Version];
  progress[output, "Run directory: " <> output];
  Do[
    id = unit["CertificateID"]; file = sourcePath[unit["File"]];
    progress[output, "Running " <> id <> " (" <> ToString[unit["ExpectedChecks"]] <> " checks) ..."];
    safeExport[FileNameJoin[{output, "status.json"}], Join[progressStatus,
      <|"CurrentUnit" -> id, "CompletedUnits" -> Length[reports]|>], "RawJSON"];
    messages = {}; evaluationOK = True;
    {elapsed, raw} = AbsoluteTiming[Block[
      {$Context = "Global`", $ContextPath = {"System`"}, $MessageList = {},
       $Assumptions = True, $IterationLimit = Infinity, $RecursionLimit = 4096, $HistoryLength = 0},
      diagnostic = Catch[
        CheckAbort[Check[Get[file], evaluationOK = False; $Failed], $Aborted],
        _, Function[{value, tag}, evaluationOK = False;
          <|"UnexpectedThrow" -> ToString[{tag, value}, InputForm]|>]];
      messages = ToString[#, InputForm] & /@ $MessageList;
      diagnostic
    ]];
    If[raw === $Aborted, aborted = True];
    rows = If[AssociationQ[raw], Lookup[raw, "Checks", {}], {}];
    allRowsValid = ListQ[rows] && And @@ (validCheck /@ rows);
    evaluationOK = TrueQ[evaluationOK] && AssociationQ[raw] &&
      TrueQ[Lookup[raw, "EvaluationOK", True]] && messages === {} &&
      Lookup[raw, "CertificateID", ""] === id &&
      Lookup[raw, "ExpectedChecks", -1] === unit["ExpectedChecks"] && allRowsValid;
    If[!TrueQ[allRowsValid], rows = {}];
    countOK = Length[rows] === unit["ExpectedChecks"];
    (* Repeated display names are permitted; exact order and counts define coverage. *)
    catalogOK = Lookup[rows, "Name", {}] === unit["CheckNames"];
    rows = MapIndexed[normalizeCheck[#1, id, First[#2]] &, rows];
    selectedPassed = TrueQ[evaluationOK && countOK && catalogOK &&
      And @@ Lookup[rows, "Passed", {}]];
    report = <|"CertificateID" -> id,
      "Status" -> If[aborted, "ABORTED", If[selectedPassed, "PASS", "FAIL"]],
      "ExpectedChecks" -> unit["ExpectedChecks"], "ActualChecks" -> Length[rows],
      "PassedChecks" -> Count[Lookup[rows, "Passed", {}], True],
      "EvaluationOK" -> evaluationOK, "CountOK" -> countOK, "CheckCatalogOK" -> catalogOK,
      "Messages" -> messages, "UnitMessages" -> If[AssociationQ[raw],
        ToString[Lookup[raw, "EvaluationMessages", {}], InputForm], ToString[raw, InputForm]],
      "ElapsedSeconds" -> elapsed, "SourceSHA256" -> sha256[file], "Checks" -> rows|>;
    safeExport[FileNameJoin[{output, id <> ".json"}], report, "RawJSON"];
    checkLines = (StringRiffle[{#["CheckID"], If[TrueQ[#["Passed"]], "PASS", "FAIL"],
      StringReplace[#["Name"], {"\t" -> " ", "\n" -> " "}],
      StringReplace[#["Detail"], {"\t" -> " ", "\n" -> " "}]}, "\t"] &) /@ rows;
    appendText[FileNameJoin[{output, "checks.log"}], StringRiffle[checkLines, "\n"] <> "\n"];
    AppendTo[reports, KeyDrop[report, "Checks"]];
    safePut[FileNameJoin[{output, "progress.wl"}], reports];
    progress[output, report["Status"] <> ": " <> ToString[report["PassedChecks"]] <> "/" <>
      ToString[report["ExpectedChecks"]] <> "  " <> id <> "  (" <> ToString[elapsed, InputForm] <> " s)"];
    If[aborted, Break[]],
    {unit, units}];
  after = snapshot[paths];
  safeExport[FileNameJoin[{output, "source_snapshot_after.json"}], after, "RawJSON"];
  actualIDs = Lookup[reports, "CertificateID", {}];
  complete = actualIDs === Lookup[units, "CertificateID"];
  selectedPassed = !aborted && complete && before === after &&
    And @@ ((Lookup[#, "Status", "FAIL"] === "PASS") & /@ reports);
  totals = {Total[Lookup[units, "ExpectedChecks"]], Total[Lookup[reports, "ActualChecks", {}]],
    Total[Lookup[reports, "PassedChecks", {}]]};
  status = If[aborted, "ABORTED", If[selectedPassed,
    If[chosen === All, "PASS", "PARTIAL_PASS"], "FAIL"]];
  finished = utcNow[];
  summary = <|"Status" -> status, "Release" -> release,
    "CompleteSuite" -> (chosen === All && complete && !aborted),
    "AllSelectedUnitsPassed" -> selectedPassed, "StartedUTC" -> begun,
    "StartedLocal" -> begunLocal, "FinishedUTC" -> finished, "FinishedLocal" -> localNow[],
    "ElapsedSeconds" -> (AbsoluteTime[] - begunClock), "Kernel" -> kernelInfo,
    "ExpectedUnits" -> Length[units], "CompletedUnits" -> Length[reports],
    "ExpectedChecks" -> totals[[1]], "ActualChecks" -> totals[[2]], "PassedChecks" -> totals[[3]],
    "SourcesUnchangedDuringRun" -> (before === after), "SourceSHA256" -> before,
    "Manuscript" -> Lookup[manifest, "Manuscript", <||>], "Units" -> reports|>;
  safePut[FileNameJoin[{output, "summary.wl"}], summary];
  safeExport[FileNameJoin[{output, "summary.txt"}], StringRiffle[{
    "Computational Verification -- " <> release, "Status: " <> status,
    "Units: " <> ToString[Length[reports]] <> "/" <> ToString[Length[units]],
    "Checks evaluated: " <> ToString[totals[[2]]] <> "/" <> ToString[totals[[1]]],
    "Checks passed: " <> ToString[totals[[3]]] <> "/" <> ToString[totals[[1]]],
    "Kernel: " <> $Version, "System: " <> $SystemID <> " (" <> $OperatingSystem <> ")",
    "Started UTC: " <> begun, "Finished UTC: " <> finished,
    "Elapsed seconds: " <> ToString[summary["ElapsedSeconds"], InputForm],
    "Sources unchanged during run: " <> ToString[before === after, InputForm],
    "Detailed checks: checks.log and LSC-*.json.",
    "A complete successful run requires PASS in both summary.json and status.json.",
    "PARTIAL_PASS, RUNNING, INCOMPLETE, ABORTED and FAIL are not complete verification."}, "\n"] <> "\n", "Text"];
  safeExport[FileNameJoin[{output, "summary.json"}], summary, "RawJSON"];
  progress[output, "Finished: " <> status <> ". Results: " <> output];
  (* Written last: this is the authoritative successful-completion marker. *)
  safeExport[FileNameJoin[{output, "status.json"}], KeyTake[summary,
    {"Status", "Release", "CompleteSuite", "ExpectedUnits", "CompletedUnits", "ExpectedChecks",
      "ActualChecks", "PassedChecks", "StartedUTC", "FinishedUTC", "ElapsedSeconds"}], "RawJSON"];
  summary
], "RunFailure"];

End[];
EndPackage[];
