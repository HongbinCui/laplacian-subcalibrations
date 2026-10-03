# Laplacian algebras and subcalibrations: computational verification

Verification code accompanying **Laplacian algebras and subcalibrations for minimal hypercones**, by **Hongbin Cui**.

The programs check differential identities, exact positivity certificates, and specified finite or coordinate-model calculations from the paper. The manuscript contains the mathematical proofs; the code provides reproducible algebraic checks.

## Getting started

Download the complete repository using **Code > Download ZIP**, then extract it before opening a notebook. Keep the directory structure intact: the `.nb` notebooks load accompanying `.wl` files, and some supplemental notebooks also load a shared helper.

In Mathematica 12, start a fresh kernel, open [START_HERE_Verification.nb](anc/START_HERE_Verification_v5.nb), and choose **Evaluation > Evaluate Notebook**. This entry runs the main thirteen units, containing **550 assertions**. All units must pass and successfully export their TXT and CSV records. Save the evaluated notebook explicitly.

The additional units below run separately; the main entry does not run them. Use a fresh kernel for each separate Mathematica entry.

## Verification units

| Unit | Assertions | Entry |
|---|---:|---|
| Main thirteen-unit suite | 550 | [Mathematica notebook](anc/START_HERE_Verification_v5.nb) |
| General cubic algebra | 31 | [Mathematica notebook](anc/mathematica/supplemental/Verify_Chapter_07_General_Cubic_Algebra.nb) |
| Exceptional cubics | 28 | [Mathematica notebook](anc/mathematica/supplemental/exceptional_cubics/Verify_Exceptional_Cubics.nb) |
| Cubic certificates | 38 | [Mathematica notebook](anc/mathematica/supplemental/cubic_update_20260927/Verify_Cubic_Update_20260927.nb) |
| Pullback transfer | 13 | [Mathematica notebook](anc/mathematica/supplemental/pullback_transfer_20260930/verify_pullback_transfer.nb) |
| Cartan/Jordan pencil | 65 | [Python program and instructions](anc/mathematica/supplemental/jordan_pencil_20260928/README.md) |
| Intrinsic five-generator table | 62 | [Python program and instructions](anc/mathematica/supplemental/intrinsic_five_generators_20260929/README.md) |

The 28-, 38-, and 13-assertion notebooks should report **`COMPLETE AND PASSED: True`**. They export a text transcript, a CSV table, and a JSON summary. The two Python-only units use Python 3 and SymPy; their linked instructions give commands and dependencies.

## Supplied results

The author's native **Mathematica 12.0.0 (Windows 64-bit)** runs cover all **660 Wolfram Language assertions**, with all checks passing. The two Python-only units have passing records for **127 assertions**. The resulting total of **787 program assertions includes stated overlaps**; it does not count 787 distinct mathematical results.

Evaluated notebooks and execution records are included. [Verification Records](anc/Verification_Records.md) identifies their locations, dates, environments, and source-hash scope. Preparing this repository did not change or rerun the programs or retained results. Dates and `v5`/`v6` in program names are stable identifiers, not manuscript versions.

## Verification scope and documentation

Exact identities and interval positivity certificates are distinguished from finite consistency tests. In particular, bounded-order determinant checks do not prove the all-orders assertions, and coordinate checks of particular nonradial models do not replace the intrinsic five-generator theorem. Those proofs, together with the geometric comparison and singular-set arguments, are in the manuscript.

- [Verification Details (PDF)](anc/Verification_Details.pdf): formulas, methods, and the correspondence with the paper.
- [Running guide and current section mapping](anc/mathematica/chapters/VERIFICATION_GUIDE.md).
- [Notation guide](anc/NOTATION_GUIDE.md): local code variables and manuscript notation.
- [Ancillary instructions](anc/README.txt): additional entry and export details.

When referring to these computations, cite the accompanying manuscript and identify the repository commit or release used.
