import unittest
from padepy import directAlgorithm as da
import sympy as sp


x = sp.Symbol("x")


class da_wrong_input(unittest.TestCase):
    def test_da_inputObj(self):
        pade12 = da.pade(1, 2, 2)
        self.assertEqual(pade12, None)
        func = [1, 2, "f", 3, 4]
        pade55 = da.pade(5, 5, func)
        self.assertEqual(pade55, None)
        func = (1, 2, 3, 3, 4)
        pade22 = da.pade(2, 2, func)
        self.assertEqual(pade22, None)

    def test_da_ComplexNumber(self):
        func = [1, 2, sp.I, 3, 4]
        pade22 = da.pade(2, 2, func)
        self.assertEqual(pade22, None)

    def test_da_nbCoeff(self):
        func = [1, 2, 3, 4, 5]
        pade33 = da.pade(3, 3, func)
        self.assertEqual(pade33, None)

    def test_da_inputString(self):
        pade8S = da.pade(8, "q", sp.exp(x))
        self.assertEqual(pade8S, None)
        pade8S = da.pade(8, 1.2, sp.exp(x))
        self.assertEqual(pade8S, None)
        padeS9 = da.pade("p", 9, sp.exp(x))
        self.assertEqual(pade8S, None)
        padeS9 = da.pade(1.5, 9, sp.exp(x))
        self.assertEqual(padeS9, None)
        pade11 = da.pade(1, 1, sp.exp(x), "d")
        self.assertEqual(pade11, None)
        pade11 = da.pade(1, 1, sp.exp(x), "32")
        self.assertEqual(pade11, None)


class da_bad_output(unittest.TestCase):
    def test_da_exp(self):
        func = sp.cos(x)
        pade33 = da.pade(3, 3, func)
        self.assertEqual(pade33, None)

    def test_da_exp(self):
        func = sp.sin(x)
        pade33 = da.pade(0, 2, func)
        self.assertEqual(pade33, None)


class da_correct_output(unittest.TestCase):
    def test_da_exp_if(self):
        pade00 = da.pade(0, 0, sp.exp(x))
        self.assertEqual(pade00, 1)
        pade44 = da.pade(4, 4, sp.exp(x))
        self.assertEqual(f"{pade44}", "(x**4/1680 + x**3/84 + 3*x**2/28 + x/2 + 1)/(x**4/1680 - x**3/84 + 3*x**2/28 - x/2 + 1)")
        pade43 = da.pade(4, 3, sp.exp(x))
        self.assertEqual(f"{pade43}", "(x**4/840 + 2*x**3/105 + x**2/7 + 4*x/7 + 1)/(-x**3/210 + x**2/14 - 3*x/7 + 1)")
        pade34 = da.pade(3, 4, sp.exp(x))
        self.assertEqual(f"{pade34}", "(x**3/210 + x**2/14 + 3*x/7 + 1)/(x**4/840 - 2*x**3/105 + x**2/7 - 4*x/7 + 1)")   

    def test_da_exp_fp_8(self):
        pade03 = da.pade(0, 3, sp.exp(x), 16)
        self.assertEqual(f"{pade03}", "1.0/(-0.1666666666666667*x**3 + 0.5*x**2 - 1.0*x + 1.0)")
        pade33 = da.pade(3, 3, sp.exp(x), 16)
        self.assertEqual(f"{pade33}", "(0.00833333333333333*x**3 + 0.09999999999999996*x**2 + 0.4999999999999999*x + 1.0)/(-0.00833333333333334*x**3 + 0.1*x**2 - 0.5000000000000001*x + 1.0)")
        pade32 = da.pade(3, 2, sp.exp(x), 16)
        self.assertEqual(f"{pade32}", "(0.01666666666666667*x**3 + 0.15*x**2 + 0.6*x + 1.0)/(0.05*x**2 - 0.4*x + 1.0)")

    def test_da_exp_fp_16(self):
        pade03 = da.pade(0, 3, sp.exp(x), 16)
        self.assertEqual(f"{pade03}", "1.0/(-0.1666666666666667*x**3 + 0.5*x**2 - 1.0*x + 1.0)")
        pade33 = da.pade(3, 3, sp.exp(x), 16)
        self.assertEqual(f"{pade33}", "(0.00833333333333333*x**3 + 0.09999999999999996*x**2 + 0.4999999999999999*x + 1.0)/(-0.00833333333333334*x**3 + 0.1*x**2 - 0.5000000000000001*x + 1.0)")
        pade32 = da.pade(3, 2, sp.exp(x), 16)
        self.assertEqual(f"{pade32}", "(0.01666666666666667*x**3 + 0.15*x**2 + 0.6*x + 1.0)/(0.05*x**2 - 0.4*x + 1.0)")

    def test_da_exp_fp_32(self):
        pade03 = da.pade(0, 3, sp.exp(x), 16)
        self.assertEqual(f"{pade03}", "1.0/(-0.1666666666666667*x**3 + 0.5*x**2 - 1.0*x + 1.0)")
        pade33 = da.pade(3, 3, sp.exp(x), 16)
        self.assertEqual(f"{pade33}", "(0.00833333333333333*x**3 + 0.09999999999999996*x**2 + 0.4999999999999999*x + 1.0)/(-0.00833333333333334*x**3 + 0.1*x**2 - 0.5000000000000001*x + 1.0)")
        pade32 = da.pade(3, 2, sp.exp(x), 16)
        self.assertEqual(f"{pade32}", "(0.01666666666666667*x**3 + 0.15*x**2 + 0.6*x + 1.0)/(0.05*x**2 - 0.4*x + 1.0)")


if __name__ == "__main__":
    unittest.main()
