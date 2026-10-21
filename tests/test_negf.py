"""Physics and independent numerical checks; run with unittest discovery."""
from pathlib import Path
import sys
import unittest
from dataclasses import replace
import numpy as np
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from negf import Device, transport, self_energy, supply, integrate_current, KB, Q
from reference import direct_banded, continuum, rectangular_exact
from experiments import first_resonance


class TransportPhysics(unittest.TestCase):
    def setUp(self): self.d=Device()

    def test_01_contact_causality_and_quadratic(self):
        t=self.d.hopping
        e=np.array([-3.,.01,.08,2*t,4*t+3.])
        s=self_energy(e,0.,t)
        self.assertTrue(np.all(s.imag<=0))
        np.testing.assert_allclose(s*s-(e-2*t)*s+t*t,0,atol=2e-12)
        self.assertEqual(s[0].imag,0); self.assertEqual(s[-1].imag,0)

    def test_02_uniform_chain_is_transparent(self):
        d=replace(self.d,height_eV=0.)
        r=transport(d,np.linspace(.001,.8,101))
        np.testing.assert_allclose(r['T'],1,atol=5e-13)
        self.assertLess(max(r['R']),1e-22)

    def test_03_unitarity_without_clipping(self):
        for v in [0.,.115,.2,.35]:
            r=transport(self.d,np.linspace(.0001,.6,1501),v)
            np.testing.assert_allclose(r['R']+r['T'],1,atol=5e-12)
            self.assertGreaterEqual(min(r['T']),0.)
            self.assertLessEqual(max(r['T']),1+5e-12)

    def test_04_closed_right_contact_blocks_current(self):
        r=transport(self.d,[.01,.03],-.08)
        np.testing.assert_array_equal(r['T'],0.)
        np.testing.assert_allclose(r['R'],1.,atol=2e-12)

    def test_05_gauge_invariance(self):
        e=np.linspace(.001,.5,100)
        a=transport(self.d,e,.1)['T']
        b=transport(self.d,e+.7,.1,shift=.7)['T']
        np.testing.assert_allclose(a,b,atol=2e-11)

    def test_06_contact_spectral_identity(self):
        r=transport(self.d,np.linspace(.01,.25,101),.115,local=True)
        self.assertLess(r['spectral_identity_error'],2e-9)
        self.assertTrue(np.all(r['LDOS']>=0))

    def test_07_bond_current_conservation_both_injections(self):
        r=transport(self.d,[.03,.05,.08,.12],.115,local=True)
        np.testing.assert_allclose(r['bond_left'],np.broadcast_to(r['T'][:,None],r['bond_left'].shape),atol=2e-11)
        np.testing.assert_allclose(r['bond_right'],-np.broadcast_to(r['T'][:,None],r['bond_right'].shape),atol=2e-11)

    def test_08_direct_banded_agrees_with_recursion(self):
        for e in [.0001,.025,.05,.083,.11,.3,.5]:
            for v in [0.,.1,.2]:
                self.assertAlmostEqual(direct_banded(self.d,e,v),transport(self.d,[e],v)['T'][0],places=10)

    def test_09_continuum_single_barrier_exact(self):
        for e in [.005,.04,.1,.25,.279,.281,.35,.7]:
            self.assertAlmostEqual(continuum(e,[(3.,.28)]),rectangular_exact(e,.28,3.),places=12)
        self.assertAlmostEqual(rectangular_exact(.28,.28,3.),rectangular_exact(.28+1e-9,.28,3.),places=7)

    def test_10_continuum_transparency_and_adding_spacers(self):
        self.assertAlmostEqual(continuum(.07,[(12.,0.)]),1.,places=13)
        a=continuum(.086,self.d.layers)
        b=continuum(.086,[(11.,0.),*self.d.layers,(7.,0.)])
        self.assertAlmostEqual(a,b,places=12)

    def test_11_second_order_mesh_approach(self):
        es=np.linspace(.06,.11,201)
        reference=np.array([continuum(e,self.d.layers) for e in es])
        errs=[]
        for h in [.2,.1,.05]:
            tr=transport(replace(self.d,spacing_nm=h),es)['T']; errs.append(max(abs(tr-reference)))
        self.assertTrue(3.7<errs[0]/errs[1]<4.3)
        self.assertTrue(3.7<errs[1]/errs[2]<4.3)

    def test_12_transverse_supply_independent_integral(self):
        vt=KB*300/Q
        for e in [.01,.08,.2]:
            # Integral over transverse kinetic energy, with spin in external current prefactor.
            f=lambda et: 1/(1+np.exp((e+et-.08)/vt))-1/(1+np.exp((e+et-.08+.1)/vt))
            q,_=quad(f,0.,1.,epsabs=1e-13)
            self.assertAlmostEqual(q,supply(e,.1,300.),places=12)

    def test_13_equilibrium_current_zero(self):
        self.assertEqual(integrate_current(self.d,0.)['J_A_m2'],0.)

    def test_14_bias_reversal_in_symmetric_device(self):
        for v in [.1,.2]:
            a=integrate_current(self.d,v)['J_A_m2']; b=integrate_current(self.d,-v)['J_A_m2']
            self.assertLess(abs(a+b)/a,2e-7)

    def test_15_energy_refinement_and_tail_bound(self):
        a=integrate_current(self.d,.1,rtol=1e-8,panel_eV=.008)
        b=integrate_current(self.d,.1,rtol=1e-8,panel_eV=.002,tail_kT=35.)
        self.assertLess(abs(a['J_A_m2']/b['J_A_m2']-1),3e-8)
        self.assertLess(a['tail_bound_A_m2']/a['J_A_m2'],1e-10)

    def test_16_thicker_barrier_narrows_ground_resonance(self):
        a=first_resonance(replace(self.d,barrier_nm=2.))
        b=first_resonance(replace(self.d,barrier_nm=5.))
        self.assertGreater(a[2]/b[2],20)
        self.assertAlmostEqual(b[1],1.,places=8)

    def test_17_ndr_is_an_interior_current_drop(self):
        peak=integrate_current(self.d,.115)['J_A_m2']
        valley=integrate_current(self.d,.2)['J_A_m2']
        self.assertGreater(peak/valley,4.)
        self.assertGreater(integrate_current(self.d,.35)['J_A_m2'],valley)

    def test_18_geometry_validation(self):
        with self.assertRaises(ValueError): Device(spacing_nm=.13)
        with self.assertRaises(ValueError): Device(mass=-.1)
        with self.assertRaises(ValueError): Device(barrier_nm=0.)
        with self.assertRaises(ValueError): integrate_current(self.d,.1,temperature=0.)

if __name__=='__main__': unittest.main()
