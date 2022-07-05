import unittest
from padepy import directAlgorithm as da
import sympy as sp


x = sp.Symbol("x")


class decimal_precision(unittest.TestCase):
    def test_pi_65(self):
        func = sp.pi 
        d = 65
        da_pade00 = da.pade(0, 0, func, d)
        self.assertEqual(str(da_pade00), "3.1415926535897932384626433832795028841971693993751058209749445923")

    def test_sqrtpi_65(self):
        func = sp.sqrt(sp.pi) 
        d = 65
        da_pade00 = da.pade(0, 0, func, d)
        self.assertEqual(str(da_pade00), "1.7724538509055160272981674833411451827975494561223871282138077899")

    def test_e_65(self):
        func = sp.E 
        d = 65
        da_pade00 = da.pade(0, 0, func, d)
        self.assertEqual(str(da_pade00), "2.7182818284590452353602874713526624977572470936999595749669676277")

    def test_sqrt2_65(self):
        func = sp.sqrt(2) 
        d = 65
        da_pade00 = da.pade(0, 0, func, d)
        self.assertEqual(str(da_pade00), "1.4142135623730950488016887242096980785696718753769480731766797380")



if __name__ == "__main__":
    unittest.main()
