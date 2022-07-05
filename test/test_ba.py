import unittest
from padepy import bakerAlgorithm as ba
import sympy as sp


x = sp.Symbol("x")


class ba_wrong_input(unittest.TestCase):
    def test_input_obj(self):
        ba_pade12, ba_path12 = ba.pade(1, 2, "2")
        self.assertEqual(ba_pade12, None)
        self.assertEqual(ba_path12, None)
        func = [1, 2, "f", 3, 4]
        ba_pade55 = ba.pade(5, 5, func)
        self.assertEqual(ba_pade55, None)
        func = (1, 2, "f", 3, 4)
        ba_pade55 = ba.pade(5, 5, func)
        self.assertEqual(ba_pade55, None)

    def test_complex_number(self):
        func = [1, 2, sp.I, 3, 4]
        ba_pade22 = ba.pade(2, 2, func)
        self.assertEqual(ba_pade22, None)

    def test_nn_coeff(self):
        func = [1, 2, 3, 4, 5]
        ba_pade55 = ba.pade(3, 3, func)
        self.assertEqual(ba_pade55, None)

    def test_string_inputs(self):
        ba_pade8S, ba_path8S = ba.pade(8, "q", sp.exp(x))
        self.assertEqual(ba_pade8S, None)
        self.assertEqual(ba_path8S, None)
        ba_padeS10, ba_pathS10 = ba.pade("p", 10, sp.exp(x))
        self.assertEqual(ba_padeS10, None)
        self.assertEqual(ba_pathS10, None)
        ba_pade11, ba_path11 = ba.pade(1, 1, sp.exp(x), "d")
        self.assertEqual(ba_pade11, None)
        ba_pade11, ba_path11 = ba.pade(1, 1, sp.exp(x), "16")
        self.assertEqual(ba_pade11, None)


class ba_bad_output(unittest.TestCase):
    def test_cos(self):
        pade14, path = ba.pade(1, 4, sp.cos(x))
        self.assertEqual(pade14, None)

    def test_epx_1cos_16(self):
        func = sp.exp(x) + (1/sp.cos(x))
        pade019, path = ba.pade(0, 19, func, 16)
        self.assertEqual(pade019, None)        
    
    def test_epx_1cos_16(self):
        func = sp.exp(x) + (1/sp.cos(x))
        pade021, path = ba.pade(0, 21, func, 16)
        self.assertEqual(pade021, None)   


class ba_correct_output(unittest.TestCase):
    def test_exp(self):
        pade00, path00 = ba.pade(0, 0, sp.exp(x))
        self.assertEqual(pade00, 1)
        self.assertEqual(path00, 1)
        pade05, path05 = ba.pade(0, 5, sp.exp(x))
        self.assertEqual(f"{pade05}", "1/(-x**5/120 + x**4/24 - x**3/6 + x**2/2 - x + 1)")
        self.assertEqual(f"{path05}", "Matrix([[x**5/120 + x**4/24 + x**3/6 + x**2/2 + x + 1], [x**4/24 + x**3/6 + x**2/2 + x + 1], [(x**4/120 + x**3/15 + 3*x**2/10 + 4*x/5 + 1)/(1 - x/5)], [(x**3/24 + x**2/4 + 3*x/4 + 1)/(1 - x/4)], [(x**3/60 + 3*x**2/20 + 3*x/5 + 1)/(x**2/20 - 2*x/5 + 1)], [(x**2/12 + x/2 + 1)/(x**2/12 - x/2 + 1)], [(x**2/20 + 2*x/5 + 1)/(-x**3/60 + 3*x**2/20 - 3*x/5 + 1)], [(x/4 + 1)/(-x**3/24 + x**2/4 - 3*x/4 + 1)], [(x/5 + 1)/(x**4/120 - x**3/15 + 3*x**2/10 - 4*x/5 + 1)], [1/(x**4/24 - x**3/6 + x**2/2 - x + 1)], [1/(-x**5/120 + x**4/24 - x**3/6 + x**2/2 - x + 1)]])")
    
    def test_Expression(self):
        pade12, path12 = ba.pade(1, 2)
        self.assertEqual(f"{pade12}", "(a0 + x*(-a0**2*a3 + 2*a0*a1*a2 - a1**3)/(a0*a2 - a1**2))/(x**2*(a1*a3 - a2**2)/(a0*a2 - a1**2) + x*(-a0*a3 + a1*a2)/(a0*a2 - a1**2) + 1)")        
        self.assertEqual(f"{path12}", "Matrix([[a0 + a1*x + a2*x**2 + a3*x**3], [a0 + a1*x + a2*x**2], [(a0 + x**2*(-a1*a3 + a2**2)/a2 + x*(-a0*a3 + a1*a2)/a2)/(1 - a3*x/a2)], [(a0 + x*(-a0*a2 + a1**2)/a1)/(1 - a2*x/a1)], [(a0 + x*(-a0**2*a3 + 2*a0*a1*a2 - a1**3)/(a0*a2 - a1**2))/(x**2*(a1*a3 - a2**2)/(a0*a2 - a1**2) + x*(-a0*a3 + a1*a2)/(a0*a2 - a1**2) + 1)]])")
    
    def test_exp_fp_16(self):
        pade03, path03 = ba.pade(0, 3, sp.exp(x), 16)
        self.assertEqual(f"{pade03}", "1.0/(-0.1666666666666667*x**3 + 0.5*x**2 - 1.0*x + 1.0)")
        pade33, path33 = ba.pade(3, 3, sp.exp(x), 16)
        self.assertEqual(f"{pade33}", "(0.008333333333333348*x**3 + 0.1000000000000001*x**2 + 0.5000000000000001*x + 1.0)/(-0.008333333333333327*x**3 + 0.09999999999999996*x**2 - 0.4999999999999999*x + 1.0)")
        pade32, path32 = ba.pade(3, 2, sp.exp(x), 16)
        self.assertEqual(f"{pade32}", "(0.01666666666666667*x**3 + 0.15*x**2 + 0.6*x + 1.0)/(0.05*x**2 - 0.4*x + 1.0)")    


if __name__ == "__main__":
    unittest.main()
