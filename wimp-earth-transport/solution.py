#!/usr/bin/env python3
"""WIMP Earth-shielding: benchmark reference implementation, sections 1-9.

Run from any directory: python solution.py
Produces results.csv, summary.json and diagnostics in the project directory.
This is an illustrative *hypothetical model*, not an experimental prediction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from collections import Counter

import numpy as np
from transport_core import (
    ROOT, INPUTS, SEC_PER_DAY, Medium, XenonResponse, incoming_halo_velocity,
    estimate_incoming_flux, load_json, load_csv, validate_inputs, row_stderr
)

KINDS = ('collisionless','ballistic','transport','scattered')


def config_fingerprint():
    files=sorted(INPUTS.glob('*.json'))+sorted(INPUTS.glob('*.csv'))
    hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    combined=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode('utf-8')).hexdigest()
    return combined,hashes


def build_summary_metrics(sum_x, sum_x2, n, factor, flux_rel_se):
    mean=sum_x/n
    stderr=row_stderr(sum_x,sum_x2,n)
    rates=mean*factor
    statistical_err=stderr*factor
    combined_err=np.hypot(statistical_err,rates*flux_rel_se)
    return rates,combined_err,statistical_err


def compute_mass(mchi, halo_base, interactions, detector, sim, atmosphere, crust):
    halo=dict(halo_base)
    halo['wimp_mass_GeV']=float(mchi)
    edges=np.asarray(detector['recoil_bin_edges_keV'],float)
    widths=np.diff(edges)
    isotopic=XenonResponse(halo,interactions,detector)
    flux,flux_err=estimate_incoming_flux(halo,int(sim['benchmark_flux_sampling_seed']),int(sim['benchmark_flux_estimation_samples']))
    flux_rel_se=flux_err/flux
    factor=flux*SEC_PER_DAY*1000/isotopic.mean_atomic_mass_g
    seed_list=list(map(int,sim['benchmark_seed_list']))
    perseed=int(sim['benchmark_histories_per_seed'])
    rng_init={seed:np.random.default_rng(seed) for seed in seed_list}
    # Common random incident histories across cross sections, independent across seeds.
    initials={seed:[incoming_halo_velocity(rng_init[seed],halo) for _ in range(perseed)] for seed in seed_list}
    out_rows=[];per_sigma=[]
    for sigma_idx,sigma_p in enumerate(interactions['sigma_p_cm2_grid']):
        sigma_p=float(sigma_p)
        medium=Medium(halo,interactions,sim,sigma_p,atmosphere,crust)
        nbins=len(widths)
        sums=np.zeros((4,nbins)); sums2=np.zeros_like(sums)
        integrated_sums=np.zeros(4);integrated_squares=np.zeros(4)
        perseed_integrated=[]
        perseed_sum=[]; perseed_sumsq=[]
        status=Counter(); multiplicity=Counter(); targets=Counter()
        n_scattered_reached=0;n_zero_reached=0
        expected_zero=expected_zero_var=0.
        total_max_e=total_max_p=0.
        grazing_weight_max=0.
        incoming_speed_sum=arrival_speed_sum=arrival_cos_sum=0.
        arrival_count=0
        for seed in seed_list:
            rng=np.random.default_rng(seed + 1000*(sigma_idx+1) + int(round(mchi)))
            batch_sum=np.zeros_like(sums);batch_sq=np.zeros_like(sums)
            batch_integrated=np.zeros(4);batch_integrated_sq=np.zeros(4)
            for initial in initials[seed]:
                speed=np.linalg.norm(initial); dz=float(initial[2]/speed)
                incoming_speed_sum+=speed
                tau=medium.straight_tau(initial)
                prob0=math.exp(-tau)
                expected_zero+=prob0
                expected_zero_var+=prob0*(1-prob0)
                base=isotopic.cross_section_bins(initial,sigma_p)/dz
                contributions=np.zeros((4,nbins))
                contributions[0]=base
                contributions[1]=base*prob0
                history=medium.propagate(initial,rng)
                status[history.status]+=1
                multiplicity[history.collisions]+=1
                targets['atmosphere']+=history.atmosphere_collisions
                targets['crust']+=history.crust_collisions
                total_max_e=max(total_max_e,history.max_energy_relative_error)
                total_max_p=max(total_max_p,history.max_momentum_relative_error)
                if history.status=='reached':
                    arrival_count+=1
                    arr_speed=float(np.linalg.norm(history.velocity))
                    arrival_speed_sum+=arr_speed
                    cosine=float(history.velocity[2]/arr_speed)
                    arrival_cos_sum+=cosine
                    angular_factor=1/cosine
                    grazing_weight_max=max(grazing_weight_max,angular_factor)
                    contributions[2]=isotopic.cross_section_bins(history.velocity,sigma_p)*angular_factor
                    if history.collisions==0:
                        n_zero_reached+=1
                    else:
                        n_scattered_reached+=1
                        contributions[3]=contributions[2]
                batch_sum+=contributions
                batch_sq+=contributions*contributions
                totals=contributions.sum(axis=1)
                batch_integrated+=totals
                batch_integrated_sq+=totals*totals
            perseed_sum.append(batch_sum)
            perseed_sumsq.append(batch_sq)
            perseed_integrated.append((batch_integrated,batch_integrated_sq))
            integrated_sums+=batch_integrated;integrated_squares+=batch_integrated_sq
            sums+=batch_sum
            sums2+=batch_sq
        N=len(seed_list)*perseed
        rates,errs,stat_errs=build_summary_metrics(sums,sums2,N,factor,flux_rel_se)
        integrated_rates= rates.sum(axis=1)
        integrated_stat_err=row_stderr(integrated_sums,integrated_squares,N)*factor
        integrated_err=np.hypot(integrated_stat_err,integrated_rates*flux_rel_se)
        ess=np.divide(integrated_sums**2, integrated_squares, out=np.zeros_like(integrated_sums), where=integrated_squares>0)
        z_zero=(n_zero_reached-expected_zero)/math.sqrt(expected_zero_var) if expected_zero_var>0 else None
        low_stat_threshold=float(sim['benchmark_significance_z'])
        bin_resolution=(rates[0]>0)&np.isfinite(stat_errs[0])&(rates[0]>=low_stat_threshold*stat_errs[0])
        sigma_info={
            'mchi_GeV':float(mchi),'sigma_p_cm2':sigma_p,'seeds':seed_list,
            'histories_per_seed':perseed,'histories':N,
            'statuses':dict(status),'zero_scatter_reached':n_zero_reached,
            'scattered_reached':n_scattered_reached,
            'zero_scatter_expected_count':expected_zero,
            'zero_scatter_expected_variance':expected_zero_var,
            'zero_scatter_z':z_zero,
            'zero_scatter_pass':z_zero is None or abs(z_zero)<low_stat_threshold,
            'collision_multiplicity':dict(sorted(multiplicity.items())),
            'collision_target_counts':dict(targets),
            'mean_incoming_speed_km_s':incoming_speed_sum/N,
            'mean_arrival_speed_km_s':arrival_speed_sum/arrival_count if arrival_count else None,
            'mean_arrival_cosine':arrival_cos_sum/arrival_count if arrival_count else None,
            'max_inverse_arrival_cosine':grazing_weight_max,
            'max_relative_energy_error':total_max_e,
            'max_relative_momentum_error':total_max_p,
            'integrated_rates_events_per_kg_day':dict(zip(KINDS,integrated_rates.tolist())),
            'integrated_stderr_events_per_kg_day':dict(zip(KINDS,integrated_err.tolist())),
            'integrated_statistical_stderr_events_per_kg_day':dict(zip(KINDS,integrated_stat_err.tolist())),
            'effective_history_count_by_rate':dict(zip(KINDS,ess.tolist())),
            'precision_target_relative_stderr':sim['benchmark_relative_stderr_target'],
            'rate_precision_met':{label:bool(err<=sim['benchmark_relative_stderr_target']*val) if val>0 else False for label,val,err in zip(KINDS,integrated_rates,integrated_err)},
            'scattered_component_well_sampled':bool(ess[3]>=25),
            'scattered_contribution_rare_statistics_warning':bool(ess[3]<25),
            'recoil_bin_edges_keV':edges.tolist(),
            'rate_differential_units':'events / kg / day / keV',
            'unresolved_collisionless_bins':int(np.sum(~bin_resolution)),
            'scattered_component_isolated':True,
            'rate_stderr_includes_flux_estimation':True,
            'per_seed_integrated_rates_events_per_kg_day':[
                dict(zip(KINDS,(part/perseed*factor).sum(axis=1).tolist())) for part in perseed_sum],
            'per_seed_integrated_statistical_stderr_events_per_kg_day':[
                dict(zip(KINDS,(row_stderr(mom,second,perseed)*factor).tolist()))
                for mom,second in perseed_integrated
            ],
            'convergence_first_seed_vs_combined':{
                kind:{'N_first_seed':perseed,'N_combined':N,
                      'first_seed_standard_error':float(row_stderr(perseed_integrated[0][0],perseed_integrated[0][1],perseed)[idx]*factor),
                      'combined_standard_error':float(integrated_stat_err[idx]),
                      'first_seed_to_combined_stderr_ratio':float(row_stderr(perseed_integrated[0][0],perseed_integrated[0][1],perseed)[idx]/integrated_stat_err[idx]) if integrated_stat_err[idx]>0 else None} 
                for idx,kind in enumerate(KINDS)
            },
        }
        for bin_idx,(lo,hi,width) in enumerate(zip(edges[:-1],edges[1:],widths)):
            row={
                'mchi_GeV':mchi,'sigma_p_cm2':sigma_p,
                'ER_low_keV':lo,'ER_high_keV':hi,
                'rate_unit':'events / kg / day / keV',
                'rate_collisionless':rates[0,bin_idx]/width,
                'rate_collisionless_stderr':errs[0,bin_idx]/width,
                'rate_ballistic':rates[1,bin_idx]/width,
                'rate_ballistic_stderr':errs[1,bin_idx]/width,
                'rate_transport':rates[2,bin_idx]/width,
                'rate_transport_stderr':errs[2,bin_idx]/width,
                'rate_scattered':rates[3,bin_idx]/width,
                'rate_scattered_stderr':errs[3,bin_idx]/width,
                'collisionless_bin_resolved':bool(bin_resolution[bin_idx]),
                'ratio_transport_to_collisionless':rates[2,bin_idx]/rates[0,bin_idx] if bin_resolution[bin_idx] and rates[0,bin_idx]>0 else '',
            }
            out_rows.append(row)
        per_sigma.append(sigma_info)
        print(f'mchi={mchi:g} sigma={sigma_p:.3g}: reached {status["reached"]}/{N}, scattered {n_scattered_reached}, rate full={integrated_rates[2]:.4g}/kg/day',flush=True)
    return out_rows,per_sigma,{'mass_GeV':mchi,'incoming_flux_cm2_s':flux,'incoming_flux_stderr_cm2_s':flux_err,'xe_mean_atom_mass_g':isotopic.mean_atomic_mass_g,'xenon_isotope_number_fractions':detector['xenon_isotope_number_fractions']}


def run():
    halo=load_json('halo_parameters.json');inter=load_json('interaction_parameters.json')
    det=load_json('detector_config.json');sim=load_json('simulation_config.json')
    atmosphere=load_csv('atmosphere_layers.csv');crust=load_csv('crust_composition.csv')
    validate_inputs(halo,inter,det,sim,atmosphere,crust)
    grid=halo.get('wimp_mass_grid_GeV',[halo['wimp_mass_GeV']]);hash_,files=config_fingerprint()
    rows=[];comparisons=[];fluxes=[]
    for mass in grid:
        mass_rows,mass_summaries,flux=compute_mass(float(mass),halo,inter,det,sim,atmosphere,crust)
        rows.extend(mass_rows);comparisons.extend(mass_summaries);fluxes.append(flux)
    dest=ROOT/'results.csv'
    with dest.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        writer.writeheader();writer.writerows(rows)
    summary={
        'benchmark_title':'Earth shielding of hypothetical WIMPs: transport and underground xenon recoils',
        'scientific_status':'Illustrative controlled reference benchmark, not measured laboratory constraints',
        'task_sections':'1-9 plus appendices A-B',
        'configuration_sha256':hash_, 'input_sha256':files,
        'reference_frame':halo['reference_frame'],
        'fixed_vertical_convention':halo['local_vertical_convention'],
        'seeds':sim['benchmark_seed_list'],
        'histories_per_seed':sim['benchmark_histories_per_seed'],
        'histories_per_mass_cross_section':len(sim['benchmark_seed_list'])*sim['benchmark_histories_per_seed'],
        'flux_estimation_samples':sim['benchmark_flux_estimation_samples'],
        'flux_sampling_seed':sim['benchmark_flux_sampling_seed'],
        'flux_normalization':'n_chi * E_halo[max((u-v_lab) dot inward,0)] * 1e5 cm/km, separately evaluated',
        'rate_estimator':'J_in * (86400 s/day) * (1000 g/kg)/m_Xe_mean * mean_history[(integral isotope-average dsigma/dER * efficiency dER) / cos_theta_at_first_downward_scoring_plane_crossing]',
        'rate_unit':'events / kg / day / keV',
        'integration_convention':'Rates integrated over each bin before dividing by bin width; one first downward scoring-plane crossing per initial history',
        'uncertainty_method':'Unbiased across-independent-initial-history sample variance, plus independent flux-estimation normalisation variance in quadrature per bin; bin uncertainties are correlated',
        'detector_exposure_convention':'Rate per kg xenon per day; no actual integrated exposure, finite geometry, instrumental backgrounds or resolution',
        'detector_isotopic_fractions':'Atomic (number) fractions; atom/kg uses fraction-weighted mean isotope mass',
        'flights':fluxes,
        'results':comparisons,
        'assumptions_and_limitations':[
            'Flat laterally infinite homogeneous-at-each-altitude layers; detector is an infinite scoring plane and only first downward crossing is recorded',
            'One fixed laboratory vertical artificially anti-parallel to galactic laboratory motion; sidereal rotation not modelled',
            'Stationary targets, elastic spin-independent isospin-conserving WIMP-nucleus interaction and Helm nuclear form factor',
            'Atmospheric layers, rock composition and effective Other rock component are provisional illustrative fixtures',
            'Natural xenon isotopic composition approximated by rounded atomic fractions; no detector backgrounds or resolution',
            'Finite Monte Carlo statistics may not resolve rare scattered arrivals at large optical depth',
            'The 1/cos(theta) planar crossing estimator has large variance for nearly grazing arrivals; report uncertainty and maximum inverse cosine',
            'Numerical interpolation of speed-dependent cross section and nuclear recoil CDF, with finite quadrature resolution',
            'The model includes only trajectories originating at top atmosphere; upward returns below detector depth and finite three-dimensional detector transport are not represented',
            'No claim of experimental discovery or exclusion of hypothetical parameter choices',
        ]
    }
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check',action='store_true',help='Run default frozen benchmark and invoke independent verifier')
    args=p.parse_args()
    run()
    if args.check:
        import subprocess,sys
        subprocess.run([sys.executable,str(ROOT/'verification.py')],cwd=ROOT,check=True)
