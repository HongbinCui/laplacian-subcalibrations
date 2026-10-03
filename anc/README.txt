Ancillary material for
Laplacian algebras and subcalibrations for minimal hypercones
Hongbin Cui

START HERE
Fully extract the package before opening the notebooks. In Mathematica 12,
start a fresh kernel, open START_HERE_Verification_v5.nb, and choose
Evaluation > Evaluate Notebook. The main entry runs 13 units, with 550
assertions; all units must pass and export their TXT and CSV records.
Save the evaluated notebook explicitly.

SEPARATE ENTRIES (paths relative to this directory)
1. mathematica/supplemental/Verify_Chapter_07_General_Cubic_Algebra.nb: 31 checks.
2. mathematica/supplemental/exceptional_cubics/Verify_Exceptional_Cubics.nb: 28 checks.
3. mathematica/supplemental/cubic_update_20260927/Verify_Cubic_Update_20260927.nb: 38 checks.
4. mathematica/supplemental/pullback_transfer_20260930/verify_pullback_transfer.nb: 13 checks.
5. mathematica/supplemental/jordan_pencil_20260928/verify_jordan_pencil.py: 65 checks.
6. mathematica/supplemental/intrinsic_five_generators_20260929/check_intrinsic_five_generators.py: 62 checks.

Use a fresh kernel for each separate Wolfram entry. Entries 2--4 load
Notebook_Run_Helper.wl from their parent directory; require the final
COMPLETE AND PASSED: True. The helper exports notebook_run_N_output.txt,
_results.csv and _summary.json. Keep the .wl files alongside their .nb
wrappers. Entries 5--6 require Python/SymPy; see their local README files.

CURRENT EVIDENCE
All 660 Wolfram assertions (550+31+28+38+13) have passing author-supplied
native Mathematica 12.0.0, Windows-x86-64 records. The two Python-only units
retain passing records for 127 assertions (65+62). The total of 787 counts
program assertions, including stated overlaps, not distinct mathematical
results. Verification_Records.md identifies the evidence and its scope.
The programs and retained results were not changed or rerun in preparing
this compact package. Older compatibility and audit archives are omitted.

DOCUMENTATION
The manuscript Appendix G gives a concise overview. Verification_Details.pdf
provides the detailed electronic supplement. NOTATION_GUIDE.md maps local
code variables to the paper. mathematica/chapters/VERIFICATION_GUIDE.md and
the manifests give the current chapter correspondence and running details.

Program filenames, contexts and output directories containing v5, v6 or
dates are stable program identifiers, not manuscript versions. They are
retained so that notebook dependencies and execution records remain valid.
Ancillary material is not needed to compile manuscript.tex.

SCOPE
The programs verify the stated algebraic identities, finite consistency
tests and coordinate models. General structural theorems, spectral rigidity,
actual-image inequalities, singular-set removability and geometric area
comparison are proved in the manuscript.
