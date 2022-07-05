import unittest
import sympy as sp
from padepy import directAlgorithm as da
from padepy import bakerAlgorithm as ba
from padepy.directAlgorithm import r

x = sp.Symbol("x")

class func_finite_precision(unittest.TestCase):
    def test_exp_11(self):
        func = sp.exp(x)
        d = 16
        ba_pade11, path = ba.pade(1, 1, func, d)
        da_pade11 = da.pade(1, 1, func, d)
        self.assertEqual(str(da_pade11), str(ba_pade11))

    def test_sin_cos_23(self):
        func = sp.sin(x) + sp.cos(x) 
        d = 16
        ba_pade23, path = ba.pade(1, 2, func, d)
        da_pade23 = da.pade(1, 2, func, d)
        self.assertEqual(str(da_pade23), str(ba_pade23))        

    def test_sin_cos_23(self):
        func = sp.sin(x) + sp.cos(x) 
        d = 16
        ba_pade33, path = ba.pade(3, 3, func, d)
        da_pade33 = da.pade(3, 3, func, d)
        self.assertEqual(str(da_pade33), str(ba_pade33))


class func_infinite_precision(unittest.TestCase):
    def test_exp_26(self):
        func = sp.exp(x)
        ba_pade26, path = ba.pade(2, 6, func)
        da_pade26 = da.pade(2, 6, func)
        self.assertEqual(str(da_pade26), str(ba_pade26))

    def test_sin_cos_43(self):
        func = sp.sin(x) + sp.cos(x) 
        ba_pade43, path = ba.pade(4, 3, func)
        da_pade43 = da.pade(4, 3, func)
        self.assertEqual(str(da_pade43), str(ba_pade43))    

    def test_sin_cos_95(self):
        func = sp.sin(x) + sp.cos(x) 
        ba_pade95, path = ba.pade(9, 5, func)
        da_pade95 = da.pade(9, 5, func)
        self.assertEqual(str(da_pade95), str(ba_pade95))    


class list_infinite_precision(unittest.TestCase):
    def test_Euler_22(self):
        func = [1, -1, 2, -6, 24, -120, 720, -5040, 40320]
        d = 0
        ba_pade22, path = ba.pade(2, 2, func, d)
        da_pade22 = da.pade(2, 2, func, d)
        self.assertEqual(str(da_pade22), str(ba_pade22))

    def test_Baker_23(self):
        func = [1, r( 1, 2), -r(5 , 8), r(13 ,16 ), -r(141 , 128), r(399 ,256 ), -r(2353 , 1024), r( 7205,2048 ), ]
        d = 0
        ba_pade23, path = ba.pade(2, 3, func, d)
        da_pade23 = da.pade(2, 3, func, d)
        self.assertEqual(str(da_pade23), str(ba_pade23))    
    
    def test_Baker_60(self):
        func = [1, r( 1, 2), -r(5 , 8), r(13 ,16 ), -r(141 , 128), r(399 ,256 ), -r(2353 , 1024), r( 7205,2048 ), ]
        d = 0
        ba_pade60, path = ba.pade(6, 0, func, d)
        da_pade60 = da.pade(6, 0, func, d)
        self.assertEqual(str(da_pade60), str(ba_pade60))             


class list_finite_precision(unittest.TestCase):
    def test_Euler_22(self):
        func = [1, -1, 2, -6, 24, -120, 720, -5040, 40320]
        d = 16
        ba_pade22, path = ba.pade(2, 2, func, d)
        da_pade22 = da.pade(2, 2, func, d)
        self.assertEqual(str(da_pade22), str(ba_pade22))

    def test_Baker_22(self):
        func = [1, r( 1, 2), -r(5 , 8), r(13 ,16 ), -r(141 , 128), r(399 ,256 ), -r(2353 , 1024), r( 7205,2048 ), ]
        d = 5
        ba_pade22, path = ba.pade(2, 2, func, d)
        da_pade22 = da.pade(2, 2, func, d)
        self.assertEqual(str(da_pade22), str(ba_pade22))   
    
    def test_Baker_70(self):
        func = [1, r( 1, 2), -r(5 , 8), r(13 ,16 ), -r(141 , 128), r(399 ,256 ), -r(2353 , 1024), r( 7205,2048 ), ]
        d = 32
        ba_pade70, path = ba.pade(7, 0, func, d)
        da_pade70 = da.pade(7, 0, func, d)
        self.assertEqual(str(da_pade70), str(ba_pade70))             
  


if __name__ == "__main__":
    unittest.main()