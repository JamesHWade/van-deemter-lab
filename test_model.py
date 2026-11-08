import math
import unittest
from dataclasses import replace
import numpy as np
from model import Conditions, calculate, chromatogram

class PhysicsTests(unittest.TestCase):
    def test_default_units(self):
        m=calculate(Conditions())
        self.assertEqual(m['H'],20)
        self.assertEqual(m['N'],7500)
        self.assertEqual(m['t2'],654)
    def test_length_scaling(self):
        a=calculate(Conditions());b=calculate(Conditions(L=300))
        self.assertAlmostEqual(b['Rs']/a['Rs'],math.sqrt(2))
        self.assertEqual(b['t2']/a['t2'],2)
    def test_no_selectivity(self):
        self.assertEqual(calculate(Conditions(alpha=1))['Rs'],0)
    def test_minimum(self):
        p=Conditions(B=12);u=math.sqrt(p.B/p.C)
        low=calculate(replace(p,u=u))['H']
        for v in [u-.01,u+.01]:self.assertGreater(calculate(replace(p,u=v))['H'],low)
    def test_peak_areas_and_normalization(self):
        for norm in [False,True]:
            x,a,b=chromatogram(Conditions(),norm)
            self.assertAlmostEqual(float(np.trapezoid(a,x)),1,places=5)
            self.assertAlmostEqual(float(np.trapezoid(b,x)),1,places=5)
    def test_resolution_matches_widths(self):
        m=calculate(Conditions())
        self.assertAlmostEqual(m['Rs'],2*(m['t2']-m['t1'])/(4*m['s1']+4*m['s2']))
    def test_boundaries(self):
        for u in [.1,5]:
            for alpha in [1,1.3]:
                m=calculate(Conditions(u=u,alpha=alpha))
                self.assertTrue(all(math.isfinite(v) for v in m.values()))
        with self.assertRaises(ValueError):calculate(Conditions(u=0))

if __name__=='__main__':unittest.main()
