# Computational Verification

**Release:** `20261011_memo_v7`  
**Prepared:** 2026-10-11  
**Manuscript:** *Laplacian algebras and subcalibrations for minimal hypercones*  
**Target environment:** Mathematica 12.0, Windows 64-bit.

This package contains Wolfram Language source code for the current manuscript. It is prepared for a new native Mathematica run. It contains no prior execution logs or inherited PASS results. The preparation date is not an execution date.

## Run the complete suite

1. Extract the complete archive into a writable folder. Keep its directory structure.
2. Open `StartHere.nb` in Mathematica. Use a fresh kernel: **Evaluation > Quit Kernel > Local**.
3. Evaluate the notebook's input cell with **Shift+Enter**:

   ```wolfram
   Get[FileNameJoin[{NotebookDirectory[], "RunAll.wl"}]]
   ```

Alternatively, evaluate `Get["C:/full/path/to/RunAll.wl"]` in a fresh notebook. No external package, network access or additional programming language is required.

The manifest contains **26 units**. Its `ExpectedChecks` field gives the exact assertion count, which is also printed at the start of the run. All supplementary identities are included in the full run. All units are executed sequentially from `manifest.json`; no manual supplementary step is needed. Some symbolic calculations take substantially longer than others. The runner prints and records each unit's start and completion, and it flushes each progress entry to disk.

## Results created by this run

Every invocation creates a distinct directory:

```text
results/run-<actual UTC date and time>-<unique suffix>/
```

It contains:

- `status.json`: final completion marker; it says `RUNNING` during evaluation.
- `summary.txt`, `summary.json`, `summary.wl`: totals, actual Mathematica version, system/OS information, real execution times, elapsed seconds and source SHA-256 hashes.
- `LSC-*.json`: all individual assertion results and diagnostics for each unit.
- `checks.log`: a readable tab-separated record of the assertions actually evaluated.
- `progress.log`, `progress.wl`: unit-level progress, saved as the run proceeds.
- `manifest_snapshot.json`, `source_snapshot.json`, `source_snapshot_after.json`: the manifest and source hashes bound to this run.

A successful complete execution has `Status: PASS` in **both** `summary.json` and `status.json`, with `CompleteSuite: true`, **all 26 units completed, with `ActualChecks` and `PassedChecks` equal to `ExpectedChecks`**. Each assertion must evaluate to literal `True`. Messages, unexpected throws, a wrong check count or ordered name catalog, and changed source bytes prevent a PASS. Each unit has its own package context; the runner also resets assumptions and evaluation settings around each unit.

If you abort evaluation, an intercepted abort is recorded as `ABORTED`. If the kernel or computer is terminated before the final marker is written, the result remains incomplete. Do not reinterpret `RUNNING`, `INCOMPLETE`, `ABORTED`, `FAIL` or `PARTIAL_PASS` as success. Re-running creates a fresh directory and does not overwrite the interrupted record.

## Review and publication

Preserve the complete new run directory. Send its `summary.txt` and the full results folder for review before publishing the execution record. The source package and its manifest should accompany the actual output on GitHub. Date fields in the new results are generated at execution time; they are never copied from the release date.

`CertificateIndex.md` and `manifest.json` state the mathematical scope and manuscript references. Finite coordinate checks retain that scope and do not replace a general proof. Code execution verifies the listed algebraic certificates; analytic arguments remain in the manuscript.

## Optional rerun of selected units

For diagnosis only:

```wolfram
Get["C:/full/path/to/Kernel/Runner.wl"];
LaplacianCertificates`Runner`RunSuite[
  "CertificateIDs" -> {"LSC-MULTIPLIER-JACOBI"}
]
```

A selected-unit run can return only `PARTIAL_PASS`, never a complete `PASS`. The optional `"ResultsDirectory" -> "C:/some/writable/folder"` chooses the parent directory under which a new timestamped run folder is created. Run the unchanged full package again to obtain a complete execution record.
