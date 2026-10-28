#!/usr/bin/env python3
"""Numerical figure checks; run with .venv/bin/python tests/unit/test_workflow_spectral.py."""
import sys
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm
PROJECT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PROJECT))
from tools.manuscript.render_workflow_spectral import (
    MODELS, VERTICES, MODE_BAR_TIME, MODE_BAR_WEIGHTS, spectral_model, probability, mode_weights,
)


def main():
    checks=0
    for name,spec,_ in MODELS:
        m=spectral_model(spec); q,pi=m['q'],m['pi']
        assert_allclose(q.sum(1),0,atol=1e-13)
        assert_allclose(pi@q,0,atol=1e-13)
        assert_allclose(-pi@np.diag(q),1,atol=1e-13)
        assert_allclose(q.T@m['modes'],m['modes']*m['values'][1:],atol=1e-13)
        assert_allclose(m['modes'].sum(0),0,atol=1e-13)
        # The manuscript's h_i and the displayed v_i share the same left basis.
        assert_allclose(m['modes'].T@np.diag(1/pi)@m['modes'],np.eye(3),atol=1e-13)
        for p0 in [*np.eye(4),np.array([.1,.2,.3,.4]),pi]:
            for t in [0,.001,.5,2,8,200]:
                p=probability(m,p0,t)
                assert_allclose(p,expm(q.T*t)@p0,atol=2e-12)
                assert_allclose(p.sum(),1,atol=1e-13)
                assert np.min(p)>-1e-13
                if name=='JC':
                    assert_allclose(p,.25+(p0-.25)*np.exp(-4*t/3),atol=1e-13)
                checks+=1
            assert_allclose(probability(m,p0,200),pi,atol=1e-13)
            assert_allclose(probability(m,probability(m,p0,.5),2),probability(m,p0,2.5),atol=1e-13)
        expected={'JC':[-4/3]*3,'K80':[-.8,-1.6,-1.6]}
        if name in expected:
            assert_allclose(m['values'][1:],expected[name],atol=1e-13)
        if name=='JC':
            assert_allclose(mode_weights(m,MODE_BAR_TIME),np.ones(3),atol=1e-13)
        elif name=='K80':
            assert_allclose(mode_weights(m,MODE_BAR_TIME),
                            [1,np.exp(-.8*MODE_BAR_TIME),np.exp(-.8*MODE_BAR_TIME)],atol=1e-13)
        else:
            # Recover modal decay from an independent matrix exponential.
            start=np.array([1.,0.,0.,0.])
            initial=np.linalg.solve(m['basis'],start)[1:]
            evolved=np.linalg.solve(m['basis'],expm(q.T*MODE_BAR_TIME)@start)[1:]
            relative=(evolved/initial)/(evolved[0]/initial[0])
            assert_allclose(MODE_BAR_WEIGHTS,relative,atol=1e-13)
        # Orthogonal basis rotation within a repeated eigenspace must not change p(t).
        if name in ('JC','K80'):
            other={**m,'basis':m['basis'].copy()}
            angle=.57; rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
            other['basis'][:,2:4]=m['basis'][:,2:4]@rotation
            assert_allclose(probability(other,np.eye(4)[0],.8),probability(m,np.eye(4)[0],.8),atol=1e-13)
    assert np.linalg.matrix_rank(VERTICES[1:]-VERTICES[0])==3
    print(f'PASS: {checks} probability cases against scipy.linalg.expm; stationarity, normalization, left-eigenvector metric, JC/K80 spectra, model-specific weights, limits, semigroup and repeated-eigenspace invariance.')


if __name__=='__main__':
    main()
