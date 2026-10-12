# Laplacian Algebras and Subcalibrations for Minimal Hypercones

Hongbin Cui — manuscript and computational verification materials.

**Distribution updated:** 2026-10-13 · **Manuscript:** arXiv v2 submission, 166 pages  
**Recorded verification:** PASS — 26/26 units, 7134/7134 checks  
**Environment:** Wolfram Mathematica 12.0.0, Microsoft Windows 64-bit

## Start here

- [Current manuscript (PDF)](paper/Laplacian_algebras_and_subcalibrations_20261012_v15.pdf)
- [LaTeX source](paper/source/manuscript.tex)
- [Verification notebook](anc/StartHere.nb)
- [Completed run summary](anc/results/run-20261012T140950Z-cdffc097/summary.txt)
- [Machine-readable summary](anc/results/run-20261012T140950Z-cdffc097/summary.json)
- [Final completion marker](anc/results/run-20261012T140950Z-cdffc097/status.json)
- [Certificate index (frozen reference)](anc/CertificateIndex.md)

The paper studies variational properties of algebraic minimal hypercones,
their associated Laplacian algebras, and constructions through polynomial
realizations of algebraic closure relations. The code checks the algebraic
identities and finite assertions specified in its certificate index.
These checks supplement the mathematical proofs; finite examples do not
establish arbitrary-order statements or replace analytic arguments.

## Reproduce the calculation

Download and extract the complete repository or release archive into a writable
folder. Open `anc/StartHere.nb` in Mathematica 12.0, quit the existing local
kernel, and evaluate its input cell with Shift+Enter:

```wolfram
Get[FileNameJoin[{NotebookDirectory[], "RunAll.wl"}]]
```

Alternatively, use an absolute path in a fresh Mathematica notebook:

```wolfram
Get["C:/full/path/to/laplacian-subcalibrations/anc/RunAll.wl"]
```

All verification code is Wolfram Language. No Python, external package, or
network access is needed to run the suite. Each invocation creates a new
run directory under `anc/results/`; the supplied record is not overwritten.
Full success requires `PASS` in both `summary.json` and `status.json`,
`CompleteSuite: true`, and all 26 units and 7134 checks completed and passed.

## Supplied execution record

The author ran the supplied verification code in Wolfram Mathematica 12.0.0
(Windows 64-bit). The archived run began on 2026-10-12 at 14:09:51 UTC and
finished at 14:09:53 UTC (22:09:51–22:09:53 in China).
Its complete, unmodified output is in
`anc/results/run-20261012T140950Z-cdffc097/`.
All 34 recorded source SHA-256 hashes match the supplied verification files
and the before/after snapshots. No new execution is claimed by this packaging.

## Manuscript and verification versions

The current paper is the **arXiv v2 submission prepared on 2026-10-13**.
This update restores the author's original full title and includes the
anonymous-reader acknowledgment and AI-assisted figure-generation disclosure.
The PDF, LaTeX source, and this README cite the same release,
`verification-20261013-v2`.
The PDF filename is retained to preserve existing links. No mathematical
statements, proofs, or verification files have changed since V15.
The unchanged verification suite retains its
original identifier **20261011_memo_v7**. These identify different artifacts.
Subsequent manuscript revisions did not change the tested identities or code.

The files inside `anc/` are preserved byte for byte, including their historical
preparation documents. Their statements such as “prepared for a new run” and
`PREPARED_NOT_YET_RUN` describe the state before the supplied execution;
the dated PASS records above establish the completed run. This root README
provides the current publication status.

`anc/manuscript/` contains the V7 PDF and label snapshot used by the frozen
certificate index. It is a provenance reference, not the current paper.
The index's old page and equation/theorem locations should be resolved by
semantic LaTeX labels in `paper/source/` when consulting the current paper. Do not delete
or edit the frozen reference files: they belong to the recorded source hashes.

## Compile the paper

With a standard TeX Live installation and latexmk:

```sh
cd paper/source
latexmk -pdf manuscript.tex
```

The paper uses ordinary `amsart`, 11pt, A4, with the contents immediately after
the abstract. Local figures are included. This repository layout is intended
for GitHub; an arXiv submission should be packaged separately.

## Integrity and citation

`SHA256SUMS` lists every distributed file other than itself, relative to the
repository root. On a system with sha256sum, run `sha256sum -c SHA256SUMS`
before generating new output or modifying the files.

For citation, use the paper's arXiv record and the release cited in the paper:
[verification-20261013-v2](https://github.com/HongbinCui/laplacian-subcalibrations/releases/tag/verification-20261013-v2).
Download the attached `Laplacian_GitHub_complete_20261013.zip` for the current
paper and the complete, unchanged verification suite and execution record.
This corrected attachment supersedes the earlier attachment of the same name.
GitHub's automatic `Source code` downloads follow the tagged commit and are
separate from this maintained attachment. This package does not assert that
arXiv v2 is already public. See [release notes](RELEASE_NOTES.md).
