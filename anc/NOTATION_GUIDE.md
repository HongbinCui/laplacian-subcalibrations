# Code-to-paper notation

Stable program filenames refer to the verification units, not current chapter
numbers. See `mathematica/chapters/VERIFICATION_GUIDE.md` for the current map.
The following local variable names are intentional; they do not require global
renaming of the executable programs.

| Scope | Paper notation | Code notation or interpretation |
|---|---|---|
| Radial calculations, especially `CheckRadialMaster` in `Verify_Appendix_A.wl` | $R=\lvert x\rvert^2$, ambient dimension $N$ | `r=R`, `n=N`; code `r` is a squared radius |
| Clifford cubics, Chapter 9 | $S=\lvert\nabla F\rvert^2$, $J=\Gamma(S,S)=\lvert\nabla S\rvert^2$, $K=\Delta S$ | `ss` or the locally specified gradient-square variable; `jGradS/jLink/jValue` for $J$, `ds` for $K$ |
| Clifford Jacobi test | $\mathsf M_J$ | `JacobiMatrix`; it is a matrix, not the scalar $J$ |
| Nonharmonic certificate, Appendix F | $H_2=\lvert\operatorname{Hess}F\rvert^2$, $T=\Gamma(S,S)$ | `h2`, `tt` |
| Complex determinants, Chapter 12 | $P=a^{-2}\lvert\nabla G\rvert^2=Q+F^2K$, $Q=E_{m-1}$, $K=2c+t$ | This local $K$ is a coefficient and is not the Clifford $K$ or $J$ |
| Retained Cartan/Lawson pullback unit, Chapter 11 | Target $q=\lvert y\rvert^2$, source $Q=q(P)$ | Local code `t` is target $q$; local code `q` is source $Q$ |
| Pullback-transfer supplement | $R=\lvert x\rvert^2$, $h=H\circ P$, $W_H=\lvert\nabla_yH\rvert^2$ | `r`, `hh`, `wH`; `wH` means $W_H(P)$ and $\lvert\nabla_xh\rvert^2=4R\,W_H(P)$ |
| Pullback-transfer supplement | Target $q=\lvert y\rvert^2$, exponent $a$, source dimension $N$ | `targetq`, `aa`, `nn` |

## Cubic normalizations

For the Clifford polynomial of Chapter 9, its normalization for the common
radial criterion in Chapter 7 is

$$
f=\frac F2,\qquad s=\frac S4,\qquad
\kappa=\frac K4,\qquad D=\frac{vT}{4}.
$$

Consequently $f\sqrt{Rs}=\frac14F\sqrt{RS}$; the two potentials have the same
normalized gradient field wherever their gradients do not vanish.

The explicit nonradial models in Chapter 8 use
$(\rho,\nu,A,\mathsf B,W)$ and $F_8=A-\mathsf B/3$.
Their realization in the intrinsic normalization of Chapter 10 is

$$
(r,s,a,p,h)
 =\left(\rho,\nu,\frac{3\sqrt3}{2}A,
          -\frac{\mathsf B}{2\sqrt3},\frac{W}{48}\right).
$$

Here the intrinsic cubic is $F_{10}=(a+3p)/6$, so that

$$
F_{10}=\frac{\sqrt3}{4}F_8,\qquad
q_{10}=\frac3{16}q_8
\quad\text{when }\mathcal L(F_i)=-q_iF_i.
$$

The real and complex coordinate realizations of this dictionary occur in
`mathematica/supplemental/intrinsic_five_generators_20260929/`.
The lowercase $r,s$ in this intrinsic table denote the squared norms of the
two orthogonal components, not the source radius and gradient square of
another chapter.

## What the Laplacian assertion checks

In `mathematica/supplemental/pullback_transfer_20260930/`, the assertion
`weighted Clifford pullback: Laplacian / R^a` checks the specialization
of the entered radial Laplacian formula under the Clifford chain rule.
An independent second-order differentiation of the radial formula is
performed by `CheckRadialMaster` in
`mathematica/chapters/Verify_Appendix_A.wl`, using `LaplaceExpr`.
The new assertion is a consistency check and does not increase the number
of independently derived formulas.

## Additional identifications

- `Verify_Appendix_D.wl`: local `n` is the ambient dimension $N$.
- In Chapter 8, the Cartan normalization is $F=\Phi/(3\sqrt3)$ when
  $\Phi$ is the Cartan--Münzner cubic of Chapter 5. It gives
  $S=R^2/3$, $D=R^3/27$, and $K=4(N+2)R/3$.
- Chapter 7's radial certificate has $q=2R$; Chapter 10's maximal-idempotent
  radial normalization has $q=3R/2$. They are related by
  $F_7=(2/\sqrt3)F_{10}$, and generate the same smallest algebra.
- In the cubic pullback, $T=S\circ P$ is the target gradient square
  composed with the Clifford map; the source gradient square is
  $|\nabla_x(f\circ P)|^2=4RT$.
- In `DirectDeterminantC[xs]` and `FormalDeterminantC[m,es]`, `xs` are the
  reciprocal singular-value variables $y_i=x_i^{-1}$, and `es[[j+1]]`
  represents $e_j(y)$. Their local `n` counts formal variables and `r`
  is an induction index; neither is a geometric radius or cone dimension.
- For minimal defining equations, $\lambda_F=-q$ when
  $\mathcal L(F)=-qF$. This $q$ is unrelated to the nonnegative divergence
  weight $q_G=q_{X_G}$ in the comparison theorem.
- Radial exponents satisfy $r^\sigma=R^{\sigma/2}$; for
  $G=r^\sigma S^\gamma F$ with $\deg F=e$, the potential degree is
  $\ell=e+\sigma+2\gamma(e-1)$.

The numerical coefficient field may contain exact algebraic numbers such as
$\sqrt3$; exact arithmetic here is not limited to rational arithmetic.
Program identifiers and local variable names are retained to preserve the
correspondence with the supplied execution records.
