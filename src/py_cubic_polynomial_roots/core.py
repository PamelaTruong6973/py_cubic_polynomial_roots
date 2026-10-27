import cmath
import math


def _sorted_by_real_imag(roots):
    """Stable ordering by real part, then imaginary part.

    Without a fixed ordering, float round-off would make the position of
    a particular root non-deterministic across platforms and test runs.
    """
    return sorted(roots, key=lambda z: (z.real, z.imag))


def _cube_root_cplx(z):
    """Return all three cube roots of a complex *z*.

    We always take the complex route even for negative real discriminants
    so that the three branches are consistent. ``cmath.pow`` chooses a
    single branch, which is fine for one root but does not give us all
    three; here we rotate explicitly by the two non-trivial cube roots of
    unity (±sqrt(3)/2 i halves).
    """
    r, theta = cmath.polar(z)
    if r == 0.0:
        return [complex(0.0, 0.0)] * 3
    r = r ** (1.0 / 3.0)
    c0 = cmath.rect(r, theta / 3.0)
    c1 = cmath.rect(r, (theta + 2.0 * math.pi) / 3.0)
    c2 = cmath.rect(r, (theta + 4.0 * math.pi) / 3.0)
    return [c0, c1, c2]


def _solve_cubic_cardano(A2, A1, A0):
    """Solve x^3 + A2 x^2 + A1 x + A0 = 0 via Cardano's formula.

    Returns the three roots (complex) in a fixed canonical order.
    """
    p = A1 - A2 * A2 / 3.0
    q = (2.0 * A2 * A2 * A2) / 27.0 - (A2 * A1) / 3.0 + A0

    D = (q * q) / 4.0 + (p * p * p) / 27.0

    u_roots = _cube_root_cplx(-q / 2.0 + cmath.sqrt(D))
    v_roots = _cube_root_cplx(-q / 2.0 - cmath.sqrt(D))

    # Product constraint: u * v = -p / 3 (from Cardano).
    # Each u must be paired with the v whose product is closest to -p/3.
    target = -p / 3.0
    candidates = []
    used_v = set()
    for u in u_roots:
        best_err = None
        best_j = None
        for j, v in enumerate(v_roots):
            if j in used_v:
                continue
            err = abs(u * v - target)
            if best_err is None or err < best_err:
                best_err = err
                best_j = j
        used_v.add(best_j)
        candidates.append(u + v_roots[best_j])

    shift = -A2 / 3.0
    candidates = [c + shift for c in candidates]

    return _sorted_by_real_imag(candidates)


def _clean(z, eps=1e-12):
    """Snap tiny real/imag parts to zero to avoid spurious complex noise."""
    re = 0.0 if abs(z.real) < eps else z.real
    im = 0.0 if abs(z.imag) < eps else z.imag
    return complex(re, im)


class CubicRoots:
    """Roots of the cubic a*x**3 + b*x**2 + c*x + d.

    The leading coefficient ``a`` must be non-zero; a ``ZeroDivisionError``
    is raised otherwise because the polynomial is no longer cubic.

    Cardano's method is used internally. Results are returned as Python
    ``complex`` numbers; roots that are real in exact arithmetic come out
    with an imaginary part of exactly 0.0 (see ``_clean``).
    """

    def __init__(self, a, b, c, d):
        if a == 0:
            raise ZeroDivisionError("leading coefficient 'a' must be non-zero")
        # Normalise to monic: x^3 + A2 x^2 + A1 x + A0 = 0.
        A2 = b / a
        A1 = c / a
        A0 = d / a
        raw = _solve_cubic_cardano(A2, A1, A0)
        self._roots = [_clean(r) for r in raw]

    @property
    def roots(self):
        """List of three roots as complex numbers, sorted by (real, imag)."""
        return list(self._roots)
