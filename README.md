# py_cubic_polynomial_roots

Calculates the real and complex roots of a cubic equation using Cardano's method.

```python
from py_cubic_polynomial_roots import CubicRoots

cr = CubicRoots(a=1, b=-6, c=11, d=-6)   # (x-1)(x-2)(x-3)
for r in cr.roots:
    print(r)   # (1+0j), (2+0j), (3+0j)
```

`CubicRoots(a, b, c, d)` solves `a·x³ + b·x² + c·x + d = 0` and exposes
`.roots`, a list of three `complex` numbers sorted by `(real, imag)`.
The leading coefficient `a` must be non-zero; passing `a == 0` raises
`ZeroDivisionError`.

## Why this exists

The standard library ships no cubic root finder. For most code the right
answer is to reach for `numpy.roots`, but `numpy` is a heavy dependency
when the only thing you need is three roots of one polynomial. This
library is pure standard-library Python (only `math` and `cmath`) so it
works in any CPython 3.8+ environment with no install step.

The trade-off is precision and speed: Cardano's formula is exact in
infinite precision but numerically noisy for near-degenerate cubics
(e.g. a triple root, or three roots clustered within machine epsilon).
Internally tiny real/imaginary residuals are snapped to zero at a
`1e-12` threshold so that a real root does not come back as
`1.0000000000001 + 0.0000000000002j`. If your coefficients span more
than ~10 orders of magnitude you will want something more numerically
robust.

## Edge you will hit

When the discriminant is near zero (repeated roots), the three candidate
roots that Cardano produces can differ only at the level of float
round-off. The ordering by `(real, imag)` is deterministic, but the
*values* in that regime carry noise on the order of `sqrt(eps) ≈ 1e-8`.
Compare with a tolerance rather than `==`.
