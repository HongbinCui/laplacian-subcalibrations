Exceptional cubic verification supplement
========================================

verify_exceptional_cubics.py and Verify_Exceptional_Cubics.wl contain the
same 28 exact checks.  The first needs Python 3 and SymPy.  The second uses
syntax available in Wolfram Mathematica 12 and needs no external packages.

To run:
  python verify_exceptional_cubics.py
or, in Mathematica:
  Get["/full/path/to/Verify_Exceptional_Cubics.wl"]

Native status — 2 October 2026:
All 28 checks passed in Mathematica 12.0.0, Windows-x86-64, at
17:37:56–17:37:57 (time as recorded by the author's kernel). The evaluated
Verify_Exceptional_Cubics.nb, notebook_run_28_output.txt, _results.csv and
_summary.json preserve this run. Checks, evaluation and exports succeeded;
no evaluation messages were recorded, and the source SHA-256 matches.

To reproduce through the notebook, use a fresh kernel and choose
Evaluation > Evaluate Notebook; require COMPLETE AND PASSED: True.

The checks cover:
  * the SF gradient norm, Laplacian, defect, and positive decomposition;
  * the balanced-spectrum threshold for this particular multiplier;
  * the sqrt(RS)F certificate at Jordan parameter d=2 and drift transfer;
  * all Gram-coordinate gradient pairings and Laplacians, for an n by 3
    real matrix, including Delta(det Gram)=2(n-2)e2(Gram);
  * the compact G2 form's Euclidean gradient identity and harmonicity;
  * a nonzero four-generator Jacobian and the rank-one Hessian spectrum.

For the square-root formulas R>0 and S>0.  F may have either sign.  The
Wolfram program's PowerExpand use is restricted to these positive bases.
The abstract variables are ordered (R,S,D,F), with degrees (2,4,6,3).
Delta S=cR; in the Jordan family c=4(d+2), whereas for G2 c=24.
The G2 coefficient Delta D=10S differs from the Jordan coefficient 4S.

These are coefficient and coordinate checks.  The manuscript separately
proves the global inequalities, describes the critical sets, establishes
removability, and performs the area comparison.  Symbolic verification
does not replace those arguments or establish priority of the results.
