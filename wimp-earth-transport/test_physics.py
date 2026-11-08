#!/usr/bin/env python3
"""Numerical physics tests distinct from the file-level independent verifier."""
import math
import numpy as np
from transport_core import (Medium,AMU_G,AMU_GEV,helm_sq,nuclear_sigma_zero,
    recoil_endpoint_keV,scatter_elastic,load_json,load_csv,validate_inputs)


def test_reference_first_collision():
    halo=load_json('halo_parameters.json');inter=load_json('interaction_parameters.json')
    det=load_json('detector_config.json');sim=load_json('simulation_config.json')
    atmosphere=load_csv('atmosphere_layers.csv');crust=load_csv('crust_composition.csv')
    validate_inputs(halo,inter,det,sim,atmosphere,crust)
    med=Medium(halo,inter,sim,1e-34,atmosphere,crust)
    # Isolate pure Si-28 homogeneous slab and exact unity form factor.
    A=28; density=2.7/(A*AMU_G)
    med.layers=[(0.,1e5,((A,density),),'crust')]
    med.top=0.;med.depth=1e5
    med.tables[A]=(np.ones_like(med.vgrid),np.tile(med.xgrid,(len(med.vgrid),1)))
    med.sigma0[A]=nuclear_sigma_zero(100.,float(inter['proton_mass_GeV']),A,1e-34)
    lamb=1/(density*med.sigma0[A]);tau=1e5/lamb
    rng=np.random.default_rng(20261102)
    n=10000;without=0;first=[]
    for _ in range(n):
        r=med.propagate(np.array([0.,0.,220.]),rng)
        if r.collisions==0:without+=1
        if r.first_collision_path_cm is not None:first.append(r.first_collision_path_cm)
        assert r.max_energy_relative_error<1e-8
        assert r.max_momentum_relative_error<1e-8
    predicted=math.exp(-tau)
    stderr=math.sqrt(predicted*(1-predicted)/n)
    assert abs(without/n-predicted)<=4*stderr+1/n,(without/n,predicted)
    # Distances *conditional* on a collision within slab have truncated exponential CDF.
    for frac in (.25,.5,.75):
        p=(1-math.exp(-frac*tau))/(1-math.exp(-tau))
        measured=np.mean(np.asarray(first)<=frac*1e5)
        sem=math.sqrt(p*(1-p)/len(first))
        assert abs(measured-p)<4*sem+1/len(first),(frac,measured,p)
    print('PASS slab: P0=',without/n,'expected=',predicted,'first collision count=',len(first))


def test_isotropic_and_kinematics():
    rng=np.random.default_rng(7373)
    n=40000
    x=rng.random(n)
    assert abs(x.mean()-.5)<4/math.sqrt(12*n)
    assert abs(np.mean(x<=.5)-.5)<4*.5/math.sqrt(n)
    assert abs(float(helm_sq(0.,131))-1)<1e-13
    assert np.all(helm_sq(np.linspace(0,100,300),131)<=1.00001)
    for _ in range(2000):
        v=rng.normal(0,200,3);v[2]=abs(v[2])+10
        A=int(rng.choice([14,16,28,40,56,131]));xx=rng.random()
        after,recoil,eerr,perr=scatter_elastic(v,100.,A,xx,rng.random()*2*math.pi)
        assert 0<=recoil<=recoil_endpoint_keV(np.linalg.norm(v),100.,A)*(1+1e-12)
        assert eerr<1e-8 and perr<1e-8,(eerr,perr)
        assert np.linalg.norm(after)<=np.linalg.norm(v)*(1+1e-12)
    print('PASS recoil endpoint, Helm F(0), isotropic recoil fraction and two-body conservation')


def test_input_rejection():
    h=load_json('halo_parameters.json');i=load_json('interaction_parameters.json');d=load_json('detector_config.json');s=load_json('simulation_config.json')
    a=load_csv('atmosphere_layers.csv');c=load_csv('crust_composition.csv')
    altered=dict(d);altered['recoil_bin_edges_keV']=[0,2,1]
    try:validate_inputs(h,i,altered,s,a,c)
    except ValueError:pass
    else:raise AssertionError('Non-monotone recoil bins not rejected')
    altered=dict(d);altered['xenon_isotope_number_fractions']={'131':.8}
    try:validate_inputs(h,i,altered,s,a,c)
    except ValueError:pass
    else:raise AssertionError('Non-normalized xenon isotope fractions not rejected')
    print('PASS input schema rejects invalid energy bins and isotope fractions')


if __name__=='__main__':
    test_input_rejection()
    test_isotropic_and_kinematics()
    test_reference_first_collision()
    print('ALL PHYSICS TESTS PASS')
