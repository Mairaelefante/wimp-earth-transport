"""Reusable physics kernel for the WIMP Earth-shielding benchmark.

Geometry: infinite horizontal, piecewise homogeneous slabs, z increasing down.
The scoring convention records only the first downward crossing at detector depth.
The sampled incoming flux is number flux through the upper horizontal plane.
Units: GeV/c^2 for masses, km/s for speeds, cm for positions, keV for recoil.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / 'inputs'
AMU_GEV = 0.93149410242
AMU_G = 1.66053906660e-24
C_KM_S = 299792.458
SEC_PER_DAY = 86400.0
HBARC_GEV_FM = 0.1973269804


def load_json(name):
    with (INPUTS / name).open(encoding='utf-8') as fh:
        return json.load(fh)


def load_csv(name):
    with (INPUTS / name).open(newline='', encoding='utf-8') as fh:
        return list(csv.DictReader(fh))


def reduced_mass(mchi, mass):
    return mchi * mass / (mchi + mass)


def nuclear_sigma_zero(mchi, mp, A, sigma_p):
    mass = A * AMU_GEV
    return sigma_p * A ** 2 * (reduced_mass(mchi, mass) / reduced_mass(mchi, mp)) ** 2


def recoil_endpoint_keV(speed_kms, mchi, A):
    mass = A * AMU_GEV
    return 2 * reduced_mass(mchi, mass)**2 / mass * (np.asarray(speed_kms)/C_KM_S)**2 * 1e6


def helm_sq(energy_keV, A):
    """Squared Helm form factor, numerically regular at zero momentum transfer."""
    E = np.asarray(energy_keV, dtype=float)
    q = np.sqrt(np.maximum(0, 2 * A * AMU_GEV * E * 1e-6)) / HBARC_GEV_FM
    skin, a = 0.9, 0.52
    c = 1.23*A**(1/3) - 0.60
    rn = math.sqrt(max(0, c*c + 7*math.pi**2*a*a/3 - 5*skin*skin))
    x = q * rn
    with np.errstate(divide='ignore', invalid='ignore'):
        direct = 3*(np.sin(x)-x*np.cos(x))/x**3
    series = 1 - x*x/10 + x**4/280
    form = np.where(np.abs(x) < .025, series, direct)
    return np.maximum(0, form**2*np.exp(-(q*skin)**2))


def rotate_cm(direction, cos_theta, phi):
    """Return rotated *relative-velocity* unit vector in 3-D."""
    d = np.asarray(direction, float)
    d = d / np.linalg.norm(d)
    ref = np.array([0.,0.,1.]) if abs(d[2]) < 0.9 else np.array([1.,0.,0.])
    e1 = np.cross(ref, d)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(d, e1)
    s = math.sqrt(max(0., 1 - cos_theta*cos_theta))
    return cos_theta*d + s*(math.cos(phi)*e1 + math.sin(phi)*e2)


def scatter_elastic(velocity, mchi, A, recoil_fraction, azimuth):
    """Stationary nuclear target. Recoil fraction = E_R/E_Rmax.

    Returns (outgoing WIMP vector, E_R in keV, relative energy error,
    relative momentum error). Finite precision conservation is monitored.
    """
    v = np.asarray(velocity, float)
    mn = A * AMU_GEV
    speed = float(np.linalg.norm(v))
    cos_theta = 1 - 2 * recoil_fraction
    relative_final = speed * rotate_cm(v, cos_theta, azimuth)
    vcm = mchi/(mchi+mn) * v
    outgoing = vcm + mn/(mchi+mn)*relative_final
    nucleus_out = vcm - mchi/(mchi+mn)*relative_final
    recoil_keV = float(recoil_endpoint_keV(speed, mchi, A)*recoil_fraction)
    wimp_loss = .5*mchi*(np.dot(v,v)-np.dot(outgoing,outgoing))/(C_KM_S*C_KM_S)*1e6
    energy_error = abs(wimp_loss-recoil_keV)/max(1.,abs(recoil_keV))
    momentum_error = np.linalg.norm(mchi*v-mchi*outgoing-mn*nucleus_out)/max(1.,mchi*speed)
    return outgoing, recoil_keV, float(energy_error), float(momentum_error)


def validate_inputs(halo, interactions, detector, simulation, atmosphere, crust):
    """Raise descriptive ValueError for nonphysical values and inconsistent fixtures."""
    def pos(value, label, allow_zero=False):
        try:
            v = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f'{label}: expected a number') from exc
        if not math.isfinite(v) or (v < 0 if allow_zero else v <= 0):
            raise ValueError(f'{label}: expected finite {"nonnegative" if allow_zero else "positive"} value')
        return v

    for key in ('dark_matter_density_GeV_cm3','wimp_mass_GeV','halo_speed_scale_km_s','escape_speed_km_s'):
        pos(halo[key], key)
    lab=np.asarray(halo['lab_velocity_galactic_km_s'],float)
    if lab.shape!=(3,) or not np.all(np.isfinite(lab)) or np.linalg.norm(lab)<=0:
        raise ValueError('Laboratory velocity must be a finite nonzero three-vector')
    masses = list(halo.get('wimp_mass_grid_GeV', [halo['wimp_mass_GeV']]))
    if not masses or any(pos(m,'WIMP mass grid')<=0 for m in masses):
        raise ValueError('Empty WIMP mass grid')
    cross = interactions['sigma_p_cm2_grid']
    if not cross:
        raise ValueError('Empty cross-section grid')
    for sig in cross: pos(sig,'sigma_p_cm2')
    if interactions.get('isospin_conserving') is not True or interactions.get('spin_independent') is not True:
        raise ValueError('Benchmark physics requires isospin-conserving spin-independent scattering')
    for key in ('proton_mass_GeV','atomic_mass_unit_GeV','atomic_mass_unit_g','speed_of_light_km_s'):
        pos(interactions[key],key)
    depth=pos(detector['depth_km'],'detector depth')
    if not np.isclose(depth*1e5, pos(simulation['crust_depth_cm'],'crust depth'),rtol=1e-12):
        raise ValueError('Crust thickness must equal scoring depth')
    edges=np.asarray(detector['recoil_bin_edges_keV'],float)
    if len(edges)<2 or not np.all(np.isfinite(edges)) or edges[0]<0 or np.any(np.diff(edges)<=0):
        raise ValueError('Detector energy bin edges must be nonnegative, strictly increasing and finite')
    threshold=pos(detector['threshold_keV'],'threshold',allow_zero=True)
    if threshold>=edges[-1]:raise ValueError('Threshold must lie inside energy bin coverage')
    efficiency=pos(detector['efficiency_above_threshold'],'efficiency',allow_zero=True)
    if efficiency>1:raise ValueError('Efficiency cannot exceed one')
    isotopes=detector['xenon_isotope_number_fractions']
    if not isotopes or any(int(k)<100 or int(k)>150 for k in isotopes):
        raise ValueError('Invalid xenon isotope configuration')
    fractions=[pos(x,'xenon isotope fraction',allow_zero=True) for x in isotopes.values()]
    if abs(sum(fractions)-1)>1e-7:raise ValueError('Xenon isotope fractions must sum to one')
    if len(atmosphere)<1 or len(crust)<2:raise ValueError('Expected at least one atmospheric layer and multicomponent crust')
    alt=[]
    for row in atmosphere:
        lo=pos(row['altitude_bottom_km'],'atmosphere bottom',allow_zero=True)
        hi=pos(row['altitude_top_km'],'atmosphere top')
        if hi<=lo:raise ValueError('Atmosphere layer not positive thickness')
        alt.append((lo,hi))
        pos(row['density_g_cm3'],'atmospheric density')
        fr=[pos(row[k],k,allow_zero=True) for k in ('molecular_fraction_N2','molecular_fraction_O2','molecular_fraction_Ar')]
        if abs(sum(fr)-1)>1e-3:raise ValueError('Atmosphere molecular fractions exceed normalization tolerance')
    alt.sort()
    if alt[0][0]!=0 or any(abs(b[0]-a[1])>1e-10 for a,b in zip(alt[:-1],alt[1:])):
        raise ValueError('Atmospheric layers must be contiguous from ground to top')
    fractions=[]
    for row in crust:
        if int(row['A'])<1:raise ValueError('Bad nuclear mass number')
        pos(row['nuclear_mass_GeV'],'nuclear mass')
        fractions.append(pos(row['mass_fraction'],'crust mass fraction',allow_zero=True))
    if abs(sum(fractions)-1)>1e-7:raise ValueError('Crust mass fractions must sum to one')
    pos(simulation['crust_density_g_cm3'],'crust density')
    if len(simulation['benchmark_seed_list'])<2 or len(set(simulation['benchmark_seed_list']))!=len(simulation['benchmark_seed_list']):
        raise ValueError('Use at least two independent unique seeds')
    if int(simulation['benchmark_histories_per_seed'])<10:raise ValueError('Too few histories for benchmark')
    for key in ('benchmark_max_collisions_per_history','benchmark_flux_estimation_samples'):
        if int(simulation[key])<1:raise ValueError(f'{key} must be positive')


def incoming_halo_velocity(rng, halo):
    """Flux-weighted entrance distribution: proportional to positive laboratory v_z."""
    sigma = float(halo['halo_speed_scale_km_s'])/math.sqrt(2)
    vlab=np.asarray(halo['lab_velocity_galactic_km_s'],float)
    inward=-vlab/np.linalg.norm(vlab)
    esc=float(halo['escape_speed_km_s'])
    v_bound=esc+np.linalg.norm(vlab)
    ref=np.array([0.,0.,1.]) if abs(inward[2])<.9 else np.array([1.,0.,0.])
    ex=np.cross(ref,inward);ex/=np.linalg.norm(ex)
    ey=np.cross(inward,ex)
    while True:
        u=rng.normal(0,sigma,3)
        if np.linalg.norm(u)>esc:continue
        v=u-vlab
        vz=float(v@inward)
        if vz<=0:continue
        if rng.random()*v_bound>vz:continue
        return np.array([v@ex,v@ey,vz])


def estimate_incoming_flux(halo, seed, n):
    """Unbiased independent estimator of number flux and its standard error."""
    rng=np.random.default_rng(seed)
    std=float(halo['halo_speed_scale_km_s'])/math.sqrt(2)
    vlab=np.asarray(halo['lab_velocity_galactic_km_s'],float)
    inward=-vlab/np.linalg.norm(vlab)
    esc=float(halo['escape_speed_km_s'])
    values=[]
    while len(values)<n:
        u=rng.normal(0,std,(min(30000,max(1024,n-len(values))),3))
        u=u[np.linalg.norm(u,axis=1)<esc]
        vplus=np.maximum(0,(u-vlab)@inward)
        values.extend(vplus.tolist())
    vplus=np.asarray(values[:n])
    factor=float(halo['dark_matter_density_GeV_cm3'])/float(halo['wimp_mass_GeV'])*1e5
    return float(vplus.mean()*factor),float(vplus.std(ddof=1)/math.sqrt(n)*factor)


@dataclass
class Propagation:
    status: str
    velocity: np.ndarray
    collisions: int
    atmosphere_collisions: int
    crust_collisions: int
    max_energy_relative_error: float
    max_momentum_relative_error: float
    first_collision_path_cm: float | None = None


class Medium:
    def __init__(self, halo, interactions, simulation, sigma_p, atmosphere, crust):
        self.mchi=float(halo['wimp_mass_GeV'])
        self.mp=float(interactions['proton_mass_GeV'])
        self.depth=float(simulation['crust_depth_cm'])
        self.max_coll=int(simulation['benchmark_max_collisions_per_history'])
        self.layers=[]
        for row in atmosphere:
            f=np.array([float(row['molecular_fraction_N2']),float(row['molecular_fraction_O2']),float(row['molecular_fraction_Ar'])]);f/=f.sum()
            rho=float(row['density_g_cm3']); molecular_mass=np.dot(f,np.array([28,32,40]))
            molecules=rho/(molecular_mass*AMU_G)
            species=[(14,2*f[0]*molecules),(16,2*f[1]*molecules),(40,f[2]*molecules)]
            self.layers.append((-float(row['altitude_top_km'])*1e5,-float(row['altitude_bottom_km'])*1e5,tuple(species),'atmosphere'))
        crust_density=float(simulation['crust_density_g_cm3'])
        species=[(int(r['A']),crust_density*float(r['mass_fraction'])/(int(r['A'])*AMU_G)) for r in crust]
        self.layers.append((0.,self.depth,tuple(species),'crust'))
        self.layers.sort(key=lambda item:item[0]);self.top=self.layers[0][0]
        self.sigma_p=float(sigma_p)
        self.xgrid=np.linspace(0,1,81)
        self.vgrid=np.linspace(0,850,171)
        self.tables={}
        for A in {A for lo,hi,species,medium in self.layers for A,_ in species}:
            E=recoil_endpoint_keV(self.vgrid[:,None],self.mchi,A)*self.xgrid[None,:]
            F=helm_sq(E,A)
            cdf=np.zeros_like(F)
            cdf[:,1:]=np.cumsum(.5*(F[:,1:]+F[:,:-1])*np.diff(self.xgrid),axis=1)
            factors=cdf[:,-1].copy()
            cdf/=np.maximum(factors[:,None],1e-50)
            self.tables[A]=(factors,cdf)
        self.sigma0={A:nuclear_sigma_zero(self.mchi,self.mp,A,self.sigma_p) for A in self.tables}

    def sigma(self,A,speed):
        return self.sigma0[A]*float(np.interp(speed,self.vgrid,self.tables[A][0]))

    def rates(self,layer,speed):
        return [(A,n*self.sigma(A,speed)) for A,n in layer[2]]

    def recoil_fraction(self,A,speed,rng):
        speed=float(np.clip(speed,0,self.vgrid[-1]))
        hi=int(np.clip(np.searchsorted(self.vgrid,speed),1,len(self.vgrid)-1));lo=hi-1
        w=(speed-self.vgrid[lo])/(self.vgrid[hi]-self.vgrid[lo])
        cdf=(1-w)*self.tables[A][1][lo]+w*self.tables[A][1][hi]
        return float(np.interp(rng.random(),cdf,self.xgrid))

    def straight_tau(self,velocity):
        v=float(np.linalg.norm(velocity));dz=float(velocity[2]/v)
        return sum((hi-lo)/dz*sum(rate for _,rate in self.rates(layer,v)) for layer in self.layers for hi,lo in [(layer[1],layer[0])])

    def propagate(self,initial,rng):
        z=self.top
        velocity=np.asarray(initial,float).copy()
        ncoll=n_atm=n_crust=0
        max_e=max_p=0.
        accumulated_path=0.
        first_collision_path=None
        for _ in range(self.max_coll+1):
            speed=float(np.linalg.norm(velocity))
            if speed<=1e-12:
                return Propagation('truncated',velocity,ncoll,n_atm,n_crust,max_e,max_p,first_collision_path)
            d=velocity/speed
            if abs(d[2])<1e-12:
                return Propagation('truncated',velocity,ncoll,n_atm,n_crust,max_e,max_p,first_collision_path)
            target_tau=-math.log1p(-rng.random())
            while True:
                if z >= self.depth-1e-7 and d[2]>0:
                    return Propagation('reached',velocity,ncoll,n_atm,n_crust,max_e,max_p,first_collision_path)
                if z <= self.top+1e-7 and d[2]<0:
                    return Propagation('escaped',velocity,ncoll,n_atm,n_crust,max_e,max_p,first_collision_path)
                probe=z+(1e-6 if d[2]>0 else -1e-6)
                layer=next((L for L in self.layers if L[0]<=probe<L[1]),None)
                if layer is None:
                    return Propagation('reached' if d[2]>0 else 'escaped',velocity,ncoll,n_atm,n_crust,max_e,max_p,first_collision_path)
                boundary=layer[1] if d[2]>0 else layer[0]
                distance=max(0.,(boundary-z)/d[2])
                rates=self.rates(layer,speed)
                rate=sum(w for _,w in rates)
                optical=distance*rate
                if rate>0 and target_tau < optical:
                    distance_to_collision=target_tau/rate
                    z+=distance_to_collision*d[2]
                    if first_collision_path is None:first_collision_path=accumulated_path+distance_to_collision
                    r=rng.random()*rate
                    cumulative=0.;target=rates[-1][0]
                    for A,w in rates:
                        cumulative+=w
                        if r<cumulative:
                            target=A
                            break
                    x=self.recoil_fraction(target,speed,rng)
                    velocity,_,e,p=scatter_elastic(velocity,self.mchi,target,x,2*math.pi*rng.random())
                    max_e=max(max_e,e);max_p=max(max_p,p)
                    ncoll+=1
                    if layer[3]=='atmosphere':n_atm+=1
                    else:n_crust+=1
                    break
                target_tau-=optical
                z=boundary
                accumulated_path+=distance
            if ncoll>=self.max_coll:
                return Propagation('truncated',velocity,ncoll,n_atm,n_crust,max_e,max_p,first_collision_path)
        return Propagation('truncated',velocity,ncoll,n_atm,n_crust,max_e,max_p,first_collision_path)


class XenonResponse:
    """Bin-integrated detector cross section per xenon atom, as a speed lookup.

    Isotope number fractions (not mass fractions), averaged using xenon's mean
    atomic mass to obtain the number of atoms per kg. Efficiency is stepped.
    """
    def __init__(self, halo, interactions, detector):
        self.halo=halo
        self.interactions=interactions
        self.detector=detector
        self.mchi=float(halo['wimp_mass_GeV'])
        self.edges=np.array(detector['recoil_bin_edges_keV'],float)
        self.threshold=float(detector['threshold_keV'])
        self.eff=float(detector['efficiency_above_threshold'])
        self.isotopes={int(k):float(v) for k,v in detector['xenon_isotope_number_fractions'].items()}
        self.mean_atomic_mass_g=sum(A*weight for A,weight in self.isotopes.items())*AMU_G
        self.speed_grid=np.linspace(0.,850.,341)
        self.reference_sigma=1e-34
        self.table=np.zeros((len(self.speed_grid),len(self.edges)-1))
        nodes,weights=leggauss(16)
        for A,fraction in self.isotopes.items():
            Emax=recoil_endpoint_keV(self.speed_grid,self.mchi,A)
            sigma_zero=nuclear_sigma_zero(self.mchi,float(interactions['proton_mass_GeV']),A,self.reference_sigma)
            for j,(lo,hi) in enumerate(zip(self.edges[:-1],self.edges[1:])):
                low=max(lo,self.threshold)
                high=np.minimum(hi,Emax)
                valid=high>low
                if not np.any(valid):continue
                qmid=(high[valid]+low)/2
                qhalf=(high[valid]-low)/2
                E=qmid[:,None]+qhalf[:,None]*nodes[None,:]
                integral=(helm_sq(E,A) @ weights)*qhalf
                self.table[valid,j]+=fraction*sigma_zero*integral/Emax[valid]*self.eff

    def cross_section_bins(self, velocity, sigma_p):
        speed=float(np.linalg.norm(velocity))
        return np.array([np.interp(speed,self.speed_grid,self.table[:,j]) for j in range(self.table.shape[1])])*(sigma_p/self.reference_sigma)


def row_stderr(sum_values, sum_squares, N):
    """Unbiased sample standard error of mean for independent histories."""
    if N<2:return np.full_like(sum_values,float('nan'))
    sample_variance=np.maximum(0.,(sum_squares-sum_values*sum_values/N)/(N-1))
    return np.sqrt(sample_variance/N)
