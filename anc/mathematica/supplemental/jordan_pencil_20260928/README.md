# Exact checks for the Jordan pencil

This 65-assertion supplement uses Python 3 and SymPy 1.14.0; it is separate
from the Wolfram verification units.

## Reproduce

From this directory, install the pinned dependencies in an appropriate
Python environment and run:

```bash
python -m pip install -r requirements.txt
python verify_jordan_pencil.py
```

The script writes `verification_results.json` beside itself and exits with
status 0 exactly when every check passes. Use `--output path/to/report.json`
to choose a different report location. The report records the Python and
SymPy versions, execution time, script SHA-256, every check, the resulting
coefficients, and the regular-zero witnesses. The supplied run passed
**65 of 65 exact checks**. There is no numerical sampling.

## Identities checked

Use the manuscript normalization

$$
\rho=\tfrac12\operatorname{tr}Z^2,\qquad \nu=|u|^2,\qquad
A=\det Z,\qquad B=-u^*Zu,\qquad W=4u^*Z^2u,
$$

with the metric $\tfrac12\operatorname{Re}\operatorname{tr}(dZ^2)+|du|^2$.
For $F_t=A+tB$, the script checks

$$
\Delta F_t=0,\qquad
S_t=|\nabla F_t|^2=\tfrac43(\rho+t\nu)^2+t(t-1)W,
$$

and, with
$\mathcal L(F)=|\nabla F|^2\Delta F-\operatorname{Hess}F(\nabla F,\nabla F)$,

$$
\mathcal L(F_t)=-8[\rho+t^2(2-t)\nu]F_t
-\tfrac83t(t-1)(3t+1)(\rho+t\nu)B.
$$

The formal calculation transcribes the five-generator closure and
Laplacian table from `nonradial_models.tex`. It checks both full polynomial
identities and their coefficients in $t$, including cancellation of
the apparent fourth-degree term in the displayed formula for
$\mathcal L(F_t)$. At $t=-1/3$ it recovers the existing gradient-square
formula and $q=8(\rho+7\nu/27)$.

A second calculation independently checks **every pairing and Laplacian
in that table for $d=1$**, directly from real coordinates

$$
Z=\begin{pmatrix}x&c&b\\c&y&a\\b&a&-x-y\end{pmatrix},
\qquad u=(u_1,u_2,u_3).
$$

The coordinate metric has diagonal-coordinate block
$\left(\begin{smallmatrix}1&1/2\\1/2&1\end{smallmatrix}\right)$
and an identity block on the other six coordinates. Its inverse therefore
has block
$\left(\begin{smallmatrix}4/3&-2/3\\-2/3&4/3\end{smallmatrix}\right)$.
This avoids irrational orthonormal coordinates. The second calculation
also obtains $S_t$, $\Delta F_t$, and $\mathcal L(F_t)$ directly,
using the Hessian contraction in the definition of $\mathcal L$.

## Why the minimal parameters are exactly $0,1,-1/3$

At these three parameters the remainder term vanishes identically, so
$\mathcal L(F_t)=-q_tF_t$, with

| Parameter | Weight $q_t$ |
| --- | --- |
| $0$ | $8\rho$ |
| $1$ | $8(\rho+\nu)$ |
| $-1/3$ | $8(\rho+7\nu/27)$ |

For necessity, take $Z=\operatorname{diag}(1,1,-2)$. The following
witnesses are real and thus lie in all three spaces
$\mathbb K=\mathbb R,\mathbb C,\mathbb H$.

| Range | Vector $u$ | $\nu$ | $W$ | $\rho+t\nu$ | $\lvert D_uF_t\rvert^2$ | $\mathcal L(F_t)$ at the zero |
| --- | --- | --- | --- | --- | --- | --- |
| $t>0$ | $e_3/\sqrt t$ | $1/t$ | $16/t$ | $4$ | $16t$ | $-\tfrac{64}{3}(t-1)(3t+1)$ |
| $t<0$ | $\sqrt{-2/t}\,e_1$ | $-2/t$ | $-8/t$ | $1$ | $-8t$ | $-\tfrac{16}{3}(t-1)(3t+1)$ |

In both rows, $\rho=3$, $A=-2$, $B=2/t$, and $F_t=0$.
The strictly positive value of $|D_uF_t|^2$, where
$D_uF_t=-2tZu$, proves regularity. The last column is nonzero for
every excluded parameter: the positive branch permits only $t=1$
and the negative branch permits only $t=-1/3$. Consequently each
other nonzero real parameter has a regular zero where mean curvature
is nonzero. The parameter $t=0$ is already covered by its minimal
equation. For every real $t$,
$(Z,u)=(\operatorname{diag}(1,-1,0),0)$ is a regular zero with
$S_t=4/3$, so the regular set is nonempty.

This argument proves that **the regular zero hypersurface is minimal
if and only if $t\in\{0,1,-1/3\}$**. It does not require an
irreducibility claim for the other members of the pencil.

## Scope and limitations

- Exact formal calculations verify the pencil formulas conditional on
  the manuscript closure for $d=1,2,4$; they are not independent
  coordinate proofs of that closure for $d=2,4$.
- Independent direct-coordinate checks establish the complete closure,
  Laplacian table, and pencil formulas for $d=1$.
- Symbolic witness calculations establish the algebra and regularity
  needed for the necessity proof. The standard geometric identity
  $\operatorname{div}(\nabla F/|\nabla F|) =\mathcal L(F)/|\nabla F|^3$ on the regular set supplies the
  connection with minimality; the script does not formalize that
  differential-geometric theorem.
- The Cayley--Dickson identification $N=A+B$, irreducibility,
  absence of cylindrical directions, stability, and area minimization
  are outside these checks.

The script encodes the formulas explicitly; it does not parse the TeX
source or attest that every later editorial version matches them.
