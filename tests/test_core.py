import math
import unittest

from py_cubic_polynomial_roots import CubicRoots


class TestCubicRoots(unittest.TestCase):
    def _approx(self, got, want, tol=1e-9):
        self.assertEqual(len(got), len(want))
        for g, w in zip(got, want):
            self.assertLess(abs(g.real - w.real), tol, msg=f"real: {g} vs {w}")
            self.assertLess(abs(g.imag - w.imag), tol, msg=f"imag: {g} vs {w}")

    def test_three_real_distinct(self):
        # (x-1)(x-2)(x-3) = x^3 - 6 x^2 + 11 x - 6
        cr = CubicRoots(1, -6, 11, -6)
        self._approx(cr.roots, [complex(1, 0), complex(2, 0), complex(3, 0)])

    def test_one_real_two_complex(self):
        # (x-1)(x^2 + 1) = x^3 - x^2 + x - 1
        cr = CubicRoots(1, -1, 1, -1)
        self._approx(cr.roots, [complex(0, -1), complex(0, 1), complex(1, 0)])

    def test_triple_root(self):
        # (x-2)^3 = x^3 - 6 x^2 + 12 x - 8
        cr = CubicRoots(1, -6, 12, -8)
        self._approx(cr.roots, [complex(2, 0)] * 3)

    def test_double_real_and_single(self):
        # (x-1)^2 (x+2) = x^3 - 3 x + 2  (note: x^2 coefficient is 0)
        cr = CubicRoots(1, 0, -3, 2)
        self._approx(cr.roots, [complex(-2, 0), complex(1, 0), complex(1, 0)])

    def test_negative_leading(self):
        # -(x-1)(x-2)(x-3) = -x^3 + 6 x^2 - 11 x + 6
        cr = CubicRoots(-1, 6, -11, 6)
        self._approx(cr.roots, [complex(1, 0), complex(2, 0), complex(3, 0)])

    def test_complex_conjugate_pair_with_offset(self):
        # (x-5)((x-1)^2 + 4) = (x-5)(x^2 - 2x + 5)
        #                    = x^3 - 7 x^2 + 15 x - 25
        cr = CubicRoots(1, -7, 15, -25)
        self._approx(
            cr.roots,
            [complex(1, -2), complex(1, 2), complex(5, 0)],
        )

    def test_all_complex_no_real_root(self):
        # x^3 + 1 has roots -1, (1 ± i sqrt(3))/2
        cr = CubicRoots(1, 0, 0, 1)
        expected = [
            complex(-1.0, 0.0),
            complex(0.5, -math.sqrt(3) / 2),
            complex(0.5, math.sqrt(3) / 2),
        ]
        self._approx(cr.roots, expected)

    def test_leading_zero_raises(self):
        with self.assertRaises(ZeroDivisionError):
            CubicRoots(0, 1, 2, 3)

    def test_fractional_coefficients(self):
        # x^3 - 0.5 x^2 - 0.25 x + 0.125
        # Roots at 0.5, 0.5, -0.5  -> (x-0.5)^2 (x+0.5)
        # = (x^2 - x + 0.25)(x + 0.5)
        # = x^3 - 0.5 x^2 - 0.25 x + 0.125
        cr = CubicRoots(1, -0.5, -0.25, 0.125)
        self._approx(
            cr.roots,
            [complex(-0.5, 0), complex(0.5, 0), complex(0.5, 0)],
        )

    def test_large_coefficients(self):
        # 2 x^3 - 4 x^2 - 2 x + 4 = 2 (x^3 - 2 x^2 - x + 2)
        # = 2 (x-2)(x-1)(x+1)
        cr = CubicRoots(2, -4, -2, 4)
        self._approx(cr.roots, [complex(-1, 0), complex(1, 0), complex(2, 0)])

    def test_returns_complex_type(self):
        cr = CubicRoots(1, -6, 11, -6)
        for r in cr.roots:
            self.assertIsInstance(r, complex)

    def test_roots_property_returns_copy(self):
        cr = CubicRoots(1, -6, 11, -6)
        first = cr.roots
        first.clear()
        second = cr.roots
        self.assertEqual(len(second), 3)


if __name__ == "__main__":
    unittest.main()
