#!/usr/bin/env python3
"""Independent, file-based benchmark verifier. Does not import or rerun solution.py.

Uses independent nuclear formulas, deterministic truncated-Maxwell integration,
analytic homogeneous-slab results, and the frozen acceptance rules.
"""
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

from scipy.integrate import quad

P=Path(__file__).resolve().parent
I=P/'inputs'
checks=[]

def record(label,condition,detail=''):
    success=bool(condition)
    checks.append({'name':label,'passed':success,'detail':str(detail)})
    print(('PASS' if success else 'FAIL')+': '+label+((' -- '+str(detail)) if detail else ''))


def load_json(name):return json.loads((I/name).read_text())
def load_csv(name):
    with (I/name).open(newline='') as file:return list(csv.DictReader(file))
def is_finite_positive(x,allow_zero=True):
    try:v=float(x)
    except (ValueError,TypeError):return False
    return math.isfinite(v) and (v>=0 if allow_zero else v>0)


def independent_flux(halo):
    """Integrate 1-D truncated-Maxwell marginal, independent from halo MC."""
    v0=float(halo['halo_speed_scale_km_s'])
    sd=v0/math.sqrt(2)
    esc=float(halo['escape_speed_km_s'])
    lab=halo['lab_velocity_galactic_km_s']
    V=math.sqrt(sum(float(x)*float(x) for x in lab))
    x=esc/(math.sqrt(2)*sd)
    Z=math.erf(x)-math.sqrt(2/math.pi)*(esc/sd)*math.exp(-x*x)
    const=math.exp(-esc*esc/(2*sd*sd))
    def marg(u):
        if abs(u)>=esc:return 0.
        return (math.exp(-u*u/(2*sd*sd))-const)/(Z*math.sqrt(2*math.pi)*sd)
    mean=quad(lambda u:max(0.,V+u)*marg(u),-esc,esc,epsabs=1e-9,epsrel=1e-11)[0]
    return float(halo['dark_matter_density_GeV_cm3'])/float(halo['wimp_mass_GeV'])*mean*1e5


def main():
    mandatory=['halo_parameters.json','interaction_parameters.json','detector_config.json','simulation_config.json','atmosphere_layers.csv','crust_composition.csv']
    record('six required input fixtures exist',all((I/f).is_file() for f in mandatory))
    if not all((I/f).is_file() for f in mandatory):return 1
    halo=load_json('halo_parameters.json');inter=load_json('interaction_parameters.json')
    det=load_json('detector_config.json');sim=load_json('simulation_config.json')
    atm=load_csv('atmosphere_layers.csv');crust=load_csv('crust_composition.csv')
    results=list(csv.DictReader((P/'results.csv').open(newline=''))) if (P/'results.csv').exists() else []
    summary=json.loads((P/'summary.json').read_text()) if (P/'summary.json').exists() else {}
    record('nonempty benchmark outputs',bool(results and summary))
    if not results:return 1
    atom_abun=[float(v) for v in det['xenon_isotope_number_fractions'].values()]
    record('natural xenon isotope numbers normalized',len(atom_abun)>=7 and abs(sum(atom_abun)-1)<1e-7)
    rock_mass=[float(r['mass_fraction']) for r in crust]
    record('mixed crust mass fractions normalized',len(crust)>1 and abs(sum(rock_mass)-1)<1e-7)
    record('atmospheric molecular fractions normalized', all(abs(sum(float(r[x]) for x in ('molecular_fraction_N2','molecular_fraction_O2','molecular_fraction_Ar'))-1)<1e-3 for r in atm))
    atm_coords=sorted((float(r['altitude_bottom_km']),float(r['altitude_top_km'])) for r in atm)
    record('contiguous atmospheric layers',atm_coords[0][0]==0 and all(abs(b[0]-a[1])<1e-10 for a,b in zip(atm_coords[:-1],atm_coords[1:])))
    edges=list(map(float,det['recoil_bin_edges_keV']))
    record('physical detector energy bins',len(edges)>=2 and edges[0]>=0 and all(b>a for a,b in zip(edges,edges[1:])))
    expected={(float(m),float(s),edges[i],edges[i+1]) for m in halo.get('wimp_mass_grid_GeV',[halo['wimp_mass_GeV']]) for s in inter['sigma_p_cm2_grid'] for i in range(len(edges)-1)}
    try:seen={(float(r['mchi_GeV']),float(r['sigma_p_cm2']),float(r['ER_low_keV']),float(r['ER_high_keV'])) for r in results}
    except (ValueError,KeyError):seen=set()
    record('exact configured mass-cross-section-bin coverage',len(seen)==len(results)==len(expected) and seen==expected, f'{len(results)} rows')
    columns=['rate_collisionless','rate_ballistic','rate_transport','rate_scattered']
    valid=True
    for row in results:
        if row.get('rate_unit')!='events / kg / day / keV': valid=False
        if not all(is_finite_positive(row.get(k),allow_zero=True) for k in columns+[k+'_stderr' for k in columns]):valid=False
        if row.get('collisionless_bin_resolved') not in ('True','False'):valid=False
        if row.get('ratio_transport_to_collisionless') and not is_finite_positive(row['ratio_transport_to_collisionless']):valid=False
    record('correct rates, uncertainties and declared units',valid)
    record('scattered-only component bounded by full transport', all(float(row['rate_scattered'])<=float(row['rate_transport'])*(1+1e-10) for row in results))
    record('analytic ballistic spectrum never exceeds collisionless',all(float(r['rate_ballistic'])<=float(r['rate_collisionless'])*(1+1e-10)+1e-14 for r in results))
    # Independent 100 GeV pure silicon, F=1 analytic reference.
    mchi=100.;A=28.;mp=float(inter['proton_mass_GeV']);amu=float(inter['atomic_mass_unit_GeV']);amug=float(inter['atomic_mass_unit_g'])
    sigma_p=1e-34;v=220.;light=float(inter['speed_of_light_km_s'])
    mn=A*amu
    mup=mchi*mp/(mchi+mp);mun=mchi*mn/(mchi+mn)
    nsil=2.7/(A*amug)
    sigma_si=sigma_p*A*A*(mun/mup)**2
    lambda_cm=1/(nsil*sigma_si)
    tau=1e5/lambda_cm
    endpoint=2*mun*mun/mn*(v/light)**2*1e6
    record('independent pure silicon reference tau and recoil endpoint',abs(tau-0.225)<0.002 and abs(endpoint-17.7)<0.25,f'tau={tau:.6g}; Emax={endpoint:.5g} keV; mean-free-path={lambda_cm/1e5:.4g} km')
    # Standard test fixture: Gaussian *radial* truncation, integrate 1D marginal.
    bymass={float(d['mass_GeV']):d for d in summary['flights']}
    good_flux=True
    for mass in halo.get('wimp_mass_grid_GeV',[halo['wimp_mass_GeV']]):
        mass_h=dict(halo);mass_h['wimp_mass_GeV']=float(mass)
        exact=independent_flux(mass_h)
        data=bymass[float(mass)]
        measured=float(data['incoming_flux_cm2_s']);sem=float(data['incoming_flux_stderr_cm2_s'])
        good_flux &= abs(measured-exact) <= 5*sem+.005*exact
    record('independent normalized incident halo flux',good_flux)
    fingerprint={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(I.glob('*')) if f.suffix in ('.json','.csv')}
    fingerprint_hash=hashlib.sha256(json.dumps(fingerprint,sort_keys=True).encode()).hexdigest()
    record('frozen input file SHA256 provenance',fingerprint_hash==summary.get('configuration_sha256') and fingerprint==summary.get('input_sha256'))
    expected_count=int(sim['benchmark_histories_per_seed'])*len(sim['benchmark_seed_list'])
    record('independent seed bookkeeping',expected_count==summary['histories_per_mass_cross_section'] and len(set(summary['seeds']))>=2)
    groups=defaultdict(list)
    for row in results:groups[(float(row['mchi_GeV']),float(row['sigma_p_cm2']))].append(row)
    details=[]
    for info in summary['results']:
        mass=float(info['mchi_GeV']); sig=float(info['sigma_p_cm2']); name=f'm={mass:g}, sigma={sig:.3g}'
        group=groups[mass,sig]
        record(f'statuses partition ({name})',sum(info['statuses'].values())==expected_count)
        record(f'first crossing / zero and scattered classification ({name})',info['statuses'].get('reached',0)==info['zero_scatter_reached']+info['scattered_reached'])
        record(f'finite physically valid scattering diagnostics ({name})',float(info['max_relative_energy_error'])<float(sim['benchmark_tolerance_rel_energy']) and float(info['max_relative_momentum_error'])<float(sim['benchmark_tolerance_rel_momentum']))
        expected0=float(info['zero_scatter_expected_count']); expectedvar=float(info['zero_scatter_expected_variance'])
        z=(info['zero_scatter_reached']-expected0)/math.sqrt(expectedvar) if expectedvar>0 else 0
        record(f'analytic zero-scatter test ({name})',abs(z)<=float(sim['benchmark_significance_z'])+1/math.sqrt(max(1.,expectedvar)),f'z={z:.2f}')
        rec_errors=[]
        for kind in columns:
            direct=sum(float(r[kind])*(float(r['ER_high_keV'])-float(r['ER_low_keV'])) for r in group)
            saved=info['integrated_rates_events_per_kg_day'][kind.removeprefix('rate_')]
            rec_errors.append(abs(direct-saved)/max(1,abs(saved)))
        record(f'independent integrated-bin-width recovery ({name})',all(x<1e-9 for x in rec_errors))
        for kind in ('collisionless','ballistic','transport','scattered'):
            if not is_finite_positive(info['integrated_stderr_events_per_kg_day'][kind]):record(f'stderr ({kind}, {name})',False)
        # Multiple seed comparison on *independent* history groups.
        seeds=info['per_seed_integrated_rates_events_per_kg_day']
        seed_se=info['per_seed_integrated_statistical_stderr_events_per_kg_day']
        for kind in ('collisionless','ballistic','transport'):
            diff=abs(seeds[0][kind]-seeds[1][kind]);se=math.hypot(seed_se[0][kind],seed_se[1][kind])
            record(f'seed stability {kind} ({name})',diff<=float(sim['benchmark_significance_z'])*se+1e-10,f'z={diff/se if se>0 else 0:.2f}')
        details.append({'mchi_GeV':mass,'sigma_p_cm2':sig,'zero_scatter_z':z,'rate_precision_met':info['rate_precision_met'],'truncated':info['statuses'].get('truncated',0)})
    record('collisionless rate linear in sigma_p', all(
        all(abs((sum(float(row[k])*(float(row['ER_high_keV'])-float(row['ER_low_keV'])) for row in groups[m,s])/s)-(
            sum(float(row[k])*(float(row['ER_high_keV'])-float(row['ER_low_keV'])) for row in groups[m,min(sg for mm,sg in groups if mm==m)])/min(sg for mm,sg in groups if mm==m)))
            <= 1e-8*max(1,abs(sum(float(row[k])*(float(row['ER_high_keV'])-float(row['ER_low_keV'])) for row in groups[m,s])/s))
            for k in ('rate_collisionless',))
        for m,s in groups
    ))
    conclusion={'overall_pass':all(c['passed'] for c in checks),'checks':checks,'diagnostics':details,
                'scope':'Independent input, analytic, halo-normalisation and file-level checker; supplemental transport tests in test_physics.py',
                'precision_warning':'An overall PASS is a set of frozen checks, not an experimental validation or proof of statistical precision.'}
    (P/'verification_report.json').write_text(json.dumps(conclusion,indent=2)+'\n')
    print(f'\nOVERALL: {"PASS" if conclusion["overall_pass"] else "FAIL"} ({sum(c["passed"] for c in checks)}/{len(checks)} checks)')
    return 0 if conclusion['overall_pass'] else 1

if __name__=='__main__':sys.exit(main())
