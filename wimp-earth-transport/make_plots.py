#!/usr/bin/env python3
"""Visualize only verified benchmark results, without modifying numerical output."""
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

P=Path(__file__).resolve().parent
rows=list(csv.DictReader((P/'results.csv').open()))
info=json.loads((P/'summary.json').read_text())
fig,axes=plt.subplots(2,2,figsize=(12.3,8.6),layout='constrained')
keys=['rate_collisionless','rate_ballistic','rate_transport']
labels=['Collisionless','Ballistic (zero scatter)','Full transport']
for ax,data in zip(axes.flat,info['results']):
    sigma=float(data['sigma_p_cm2']);mass=float(data['mchi_GeV'])
    group=sorted((r for r in rows if float(r['sigma_p_cm2'])==sigma and float(r['mchi_GeV'])==mass),key=lambda r:float(r['ER_low_keV']))
    x=np.asarray([(float(r['ER_low_keV'])+float(r['ER_high_keV']))/2 for r in group])
    for key,name in zip(keys,labels):
        y=np.asarray([float(r[key]) for r in group]);yerr=np.asarray([float(r[key+'_stderr']) for r in group])
        mask=y>0
        if mask.any():
            line,=ax.plot(x[mask],y[mask],marker='o',ms=3,lw=1.5,label=name)
            ax.fill_between(x[mask],np.maximum(1e-30,y[mask]-yerr[mask]),y[mask]+yerr[mask],alpha=.17,color=line.get_color())
    ax.set_title(rf'$m_\chi$ = {mass:g} GeV,  $\sigma_p$ = {sigma:.2g} cm$^2$')
    ax.set_xlabel('Nuclear recoil energy (keV)')
    ax.set_ylabel('dR/dE (events / kg / day / keV)')
    ax.set_yscale('log')
    ax.grid(alpha=.25)
axes[0,0].legend(fontsize=8)
fig.suptitle('Xenon nuclear recoil predictions: collisionless vs ballistic vs transport',fontsize=13)
fig.savefig(P/'recoil_spectra.png',dpi=170)
plt.close(fig)

summary=info['results']
sigma=np.array([r['sigma_p_cm2'] for r in summary],float)
frac_zero=np.array([r['zero_scatter_reached']/r['histories'] for r in summary])
frac_expected=np.array([r['zero_scatter_expected_count']/r['histories'] for r in summary])
frac_scattered=np.array([r['scattered_reached']/r['histories'] for r in summary])
mean_speeds=np.array([r['mean_arrival_speed_km_s'] for r in summary])
mean_in=np.array([r['mean_incoming_speed_km_s'] for r in summary])
integ=np.array([[r['integrated_rates_events_per_kg_day'][k] for k in ('collisionless','ballistic','transport','scattered')] for r in summary])
fig,axes=plt.subplots(2,2,figsize=(12.3,8.5),layout='constrained')
a=axes[0,0]
a.semilogx(sigma,frac_zero,'o-',label='Simulated zero-scatter arrivals')
a.semilogx(sigma,frac_expected,'s--',label='Ray-averaged exp(-tau)')
a.semilogx(sigma,frac_scattered,'d-',label='Scattered arrivals')
a.set_ylabel('Fraction of incident histories');a.set_xlabel('Proton cross section (cm²)');a.legend(fontsize=8);a.grid(alpha=.25)
a=axes[0,1]
a.plot(sigma,mean_in,'o--',label='Incident mean speed')
a.plot(sigma,mean_speeds,'s-',label='Mean speed of arrivals')
a.set_xscale('log');a.set_ylabel('Speed (km/s)');a.set_xlabel('Proton cross section (cm²)');a.legend(fontsize=8);a.grid(alpha=.25)
a=axes[1,0]
a.plot(sigma,integ[:,1]/integ[:,0],'s-',label='Ballistic / collisionless')
a.plot(sigma,integ[:,2]/integ[:,0],'o-',label='Full / collisionless')
a.plot(sigma,integ[:,3]/integ[:,0],'d-',label='Scattered / collisionless')
a.axhline(1,color='gray',lw=.8,ls='--')
a.set_xscale('log');a.set_ylabel('Integrated rate ratio');a.set_xlabel('Proton cross section (cm²)');a.legend(fontsize=8);a.grid(alpha=.25)
a=axes[1,1]
max_n=max((int(k) for r in summary for k in r['collision_multiplicity']),default=0)
max_plot=min(7,max_n)
for r in summary[-2:]:
    ys=[r['collision_multiplicity'].get(str(k),0) for k in range(max_plot+1)]
    a.plot(range(max_plot+1),ys,'o-',label=f"sigma={r['sigma_p_cm2']:.2g}")
a.set_xlabel('Number of collisions per history')
a.set_ylabel('Number of histories')
a.set_yscale('log');a.legend(fontsize=8);a.grid(alpha=.25)
fig.suptitle('Earth transport diagnostics and scoring-plane observations',fontsize=13)
fig.savefig(P/'transport_diagnostics.png',dpi=170)
plt.close(fig)
print('Wrote recoil_spectra.png and transport_diagnostics.png')
