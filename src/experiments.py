"""Reproducible experiments and raw evidence for a prescribed-potential RTD."""
from pathlib import Path
from dataclasses import replace, asdict
import json
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad, trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from negf import Device, transport, supply, integrate_current, current_prefactor
from reference import continuum, direct_banded, rectangular_exact

plt.rcParams.update({'font.size':10,'figure.dpi':120,'axes.grid':True,'grid.alpha':.22,
                     'savefig.dpi':170,'axes.spines.top':False,'axes.spines.right':False})


def csv(out,name,header,rows):
    np.savetxt(out/name,np.asarray(list(rows)),delimiter=',',header=header,comments='',fmt='%.12g')


def save(out,name,fig):
    fig.savefig(out/name,bbox_inches='tight'); plt.close(fig)


def first_resonance(device,bias=0.):
    # Ground resonance stays in this window for zero-bias geometries studied here.
    f=lambda e: float(transport(device,[e],bias)['T'][0])
    es=np.linspace(.005,.20,2001); ts=transport(device,es,bias)['T']
    peaks,_=find_peaks(ts)
    if not len(peaks): raise RuntimeError('No ground resonance found in search window')
    i=peaks[0]
    opt=minimize_scalar(lambda e:-f(e),bounds=(es[i-1],es[i+1]),method='bounded',options={'xatol':1e-14})
    ep=opt.x; tp=f(ep)
    # Half maximum above zero; this is a transmission width, not a fitted complex pole.
    le=brentq(lambda e:f(e)-tp/2,.000001,ep,xtol=1e-14)
    ri=brentq(lambda e:f(e)-tp/2,ep,min(.23,ep+.1),xtol=1e-14)
    return ep,tp,ri-le


def curve(device,biases,temperature=300.):
    return np.array([[v,*[r[k] for k in ('J_A_m2','estimated_error_A_m2','tail_bound_A_m2','evaluations')]]
                     for v in biases for r in [integrate_current(device,v,temperature)]])


def peak_valley(arr):
    # First local maximum, then first local minimum; never use global max at upper bias endpoint.
    p,_=find_peaks(arr[:,1]); m,_=find_peaks(-arr[:,1])
    if not len(p): return None
    later=m[m>p[0]]
    if not len(later): return None
    i=p[0]; j=later[0]
    return {'peak_V':float(arr[i,0]),'valley_V':float(arr[j,0]),
            'peak_A_m2':float(arr[i,1]),'valley_A_m2':float(arr[j,1]),'PVR':float(arr[i,1]/arr[j,1])}


def baseline(out,device):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    biases=np.linspace(0,.45,91); data=curve(device,biases)
    csv(out,'iv_300K.csv','bias_V,J_A_m2,quadrature_error_A_m2,tail_bound_A_m2,evaluations',data)
    pv=peak_valley(data)
    selected=[0.,pv['peak_V'],pv['valley_V'],.35]
    en=np.linspace(.00001,.40,4000)
    rows=[]; fig,ax=plt.subplots(1,2,figsize=(11,4))
    x,_=device.profile()
    for v in selected:
        _,u=device.profile(v); spec=transport(device,en,v)['T']
        ax[0].plot(x,u,label=f'{v:.3f} V'); ax[1].plot(en,spec,label=f'{v:.3f} V')
        rows.extend(zip(np.repeat(v,len(en)),en,spec))
    ax[0].set(xlabel='Position (nm)',ylabel='Conduction edge (eV)',title='Prescribed band-edge profiles')
    ax[1].set(xlabel='Longitudinal energy (eV)',ylabel='Transmission',title='Resonance alignment under bias')
    for a in ax:a.legend(fontsize=8)
    fig.tight_layout(); save(out,'bands_transmission.png',fig)
    csv(out,'transmission.csv','bias_V,energy_eV,T',rows)
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    ax[0].plot(data[:,0],data[:,1]/1e4); ax[0].plot([pv['peak_V'],pv['valley_V']],np.array([pv['peak_A_m2'],pv['valley_A_m2']])/1e4,'o')
    ax[0].set(xlabel='Bias (V)',ylabel='Current density (A/cm²)',title=f'First peak / valley = {pv["PVR"]:.2f} (5 mV grid)')
    slope=np.gradient(data[:,1],data[:,0]); ax[1].plot(data[1:-1,0],slope[1:-1]/1e4)
    ax[1].axhline(0,color='k',lw=.8); ax[1].set(xlabel='Bias (V)',ylabel='dJ/dV (A/cm²/V)',title='Negative differential resistance')
    fig.tight_layout(); save(out,'iv_ndr.png',fig)
    csv(out,'differential_conductance.csv','bias_V,dJdV_A_m2_V',zip(data[1:-1,0],slope[1:-1]))
    # Non-equilibrium local spectral evidence and energy-resolved current conservation.
    v=pv['peak_V']; enmap=np.linspace(.0001,.25,800)
    loc=transport(device,enmap,v,local=True)
    spec=loc['T']; ep=first_resonance(device,v)[0]
    one=transport(device,[ep],v,local=True)
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    im=ax[0].pcolormesh(x,enmap,np.log10(np.maximum(loc['LDOS'],1e-8)),shading='auto',cmap='magma')
    ax[0].set(xlabel='Position (nm)',ylabel='Energy (eV)',title=f'log10 LDOS per spin at {v:.3f} V')
    fig.colorbar(im,ax=ax[0],label='log10[1/(eV nm)]')
    al=one['A_L'][0]; ar=one['A_R'][0]
    ax[1].plot(x,al/(2*np.pi*device.spacing_nm),label='Left injectivity')
    ax[1].plot(x,ar/(2*np.pi*device.spacing_nm),label='Right injectivity')
    ax[1].set(xlabel='Position (nm)',ylabel='1/(eV nm), per spin',title=f'Contact-resolved LDOS at {ep:.5f} eV'); ax[1].legend()
    fig.tight_layout(); save(out,'local_spectra.png',fig)
    np.savez_compressed(out/'local_spectra.npz',x_nm=x,energy_eV=enmap,LDOS_per_eV_nm=loc['LDOS'],
                        bias_V=v,resonance_eV=ep,A_L_per_eV=one['A_L'],A_R_per_eV=one['A_R'])
    xb=(x[:-1]+x[1:])/2
    csv(out,'bond_spectrum.csv','position_nm,left_injected_flux,right_injected_flux,T',
        zip(xb,one['bond_left'][0],one['bond_right'][0],np.repeat(one['T'][0],len(xb))))
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    ax[0].plot(xb,one['bond_left'][0],label='Left injection'); ax[0].plot(xb,one['bond_right'][0],label='Right injection')
    ax[0].set(xlabel='Bond position (nm)',ylabel='Dimensionless spectral flux',title='Equal flux on every bond'); ax[0].legend()
    ax[1].plot(enmap,loc['R']+loc['T']-1)
    ax[1].set(xlabel='Energy (eV)',ylabel='R + T - 1',title='Scattering unitarity (no clipping)')
    fig.tight_layout(); save(out,'conservation.png',fig)
    errors={'unitarity_max':float(np.max(abs(loc['R']+loc['T']-1))),
            'bond_max_abs':float(np.max(abs(one['bond_left']-one['T'][:,None]))),
            'spectral_identity_max_per_eV':float(loc['spectral_identity_error'])}
    return data,pv,errors


def sensitivity(out,device):
    biases=np.linspace(0,.35,71)
    rows=[]; summary=[]; fig,ax=plt.subplots(1,2,figsize=(11,4))
    for w in [4.6,5.,5.4]:
        d=replace(device,well_nm=w); iv=curve(d,biases); pv=peak_valley(iv)
        ax[0].plot(iv[:,0],iv[:,1]/1e4,label=f'well {w:.1f} nm')
        rows.extend((w,*r) for r in iv); summary.append({'well_nm':w,**pv})
    csv(out,'well_sensitivity.csv','well_nm,bias_V,J_A_m2,quadrature_error_A_m2,tail_bound_A_m2,evaluations',rows)
    rows=[]
    for temp in [77.,150.,300.]:
        iv=curve(device,biases,temp); pv=peak_valley(iv)
        ax[1].plot(iv[:,0],iv[:,1]/1e4,label=f'{temp:.0f} K')
        rows.extend((temp,*r) for r in iv)
    csv(out,'temperature_sweep.csv','temperature_K,bias_V,J_A_m2,quadrature_error_A_m2,tail_bound_A_m2,evaluations',rows)
    for a in ax:a.set(xlabel='Bias (V)',ylabel='Current density (A/cm²)'); a.legend()
    ax[0].set_title('Well width: fixed prescribed-voltage rule'); ax[1].set_title('Fermi supply changes; no phonon scattering')
    fig.tight_layout(); save(out,'sensitivity.png',fig)
    return summary


def verification(out,device):
    # Continuum reference is exact for piecewise constant ZERO-BIAS layer potentials.
    energies=np.linspace(.015,.24,241)
    cont=np.array([continuum(e,device.layers,device.mass) for e in energies])
    rows=[]; mesh=[]; fig,ax=plt.subplots(1,2,figsize=(11,4))
    for spacing in [.2,.1,.05]:
        d=replace(device,spacing_nm=spacing); tr=transport(d,energies)['T']
        ep,tp,width=first_resonance(d); j=integrate_current(d,.1)['J_A_m2']
        mesh.append((spacing,ep,width,j,np.max(abs(tr-cont))))
        rows.extend(zip(np.repeat(spacing,len(energies)),energies,tr,cont))
        ax[0].plot(energies,tr,label=f'FD {spacing:g} nm')
    ax[0].plot(energies,cont,'k--',label='Continuum S cascade'); ax[0].legend(fontsize=8)
    ax[0].set(xlabel='Energy (eV)',ylabel='T',title='Independent continuum comparison (V=0)')
    continuum_peak=minimize_scalar(lambda e:-continuum(e,device.layers,device.mass),bounds=(.06,.10),method='bounded',options={'xatol':1e-14}).x
    ax[1].loglog(np.array(mesh)[:,0],np.abs(np.array(mesh)[:,1]-continuum_peak)*1e3,'o-',label='Ground resonance error')
    ax[1].set(xlabel='Mesh spacing (nm)',ylabel='Energy error (meV)',title='Resonance approaches continuum'); ax[1].legend()
    fig.tight_layout(); save(out,'independent_validation.png',fig)
    csv(out,'continuum_comparison.csv','spacing_nm,energy_eV,T_discrete,T_continuum',rows)
    csv(out,'mesh_convergence.csv','spacing_nm,resonance_eV,FWHM_eV,J_at_0p1V_A_m2,max_T_error_vs_continuum',mesh)
    # Exact one-barrier formula versus a separate continuum cascade.
    es=[.015,.05,.12,.25,.28,.35,.50]
    rectangular_errors=[abs(rectangular_exact(e,.28,3.)-continuum(e,[(3.,.28)])) for e in es if e!=.28]
    direct_errors=[abs(direct_banded(device,e,.12)-transport(device,[e],.12)['T'][0]) for e in es]
    quads=[]
    for panel in [.008,.004,.002]:
        for v in [.1,.2]:
            r=integrate_current(device,v,panel_eV=panel,rtol=1e-8)
            quads.append((panel,v,r['J_A_m2'],r['estimated_error_A_m2'],r['evaluations']))
    csv(out,'quadrature_convergence.csv','initial_panel_eV,bias_V,J_A_m2,estimated_error_A_m2,evaluations',quads)
    # A different integrator with an explicit dense set of initial subintervals.
    v=.1; hi=integrate_current(device,v)['Emax_eV']
    fun=lambda e:float(transport(device,[e],v)['T'][0]*supply(e,v,300.,device.mu_eV))
    qval,qerr=quad(fun,0.,hi,points=np.linspace(0.,hi,151)[1:-1],epsabs=1e-12,epsrel=2e-8,limit=1000)
    ref_j=qval*current_prefactor(device.mass)
    # Bias-grid error is independent of energy quadrature error.
    refined=curve(device,np.arange(.085,.131,.001))
    imax=np.argmax(refined[:,1]); peak=(refined[imax,0],refined[imax,1])
    refined_valley=curve(device,np.arange(.18,.241,.001)); imin=np.argmin(refined_valley[:,1]); valley=refined_valley[imin,:2]
    csv(out,'bias_refinement.csv','bias_V,J_A_m2,quadrature_error_A_m2,tail_bound_A_m2,evaluations',np.vstack([refined,refined_valley]))
    return {'continuum_resonance_eV':float(continuum_peak),'rectangular_cascade_max_error':float(max(rectangular_errors)),
            'recursive_direct_max_T_error':float(max(direct_errors)), 'mesh_records':mesh,
            'independent_quad_J_A_m2':float(ref_j),'independent_quad_error_A_m2':float(qerr*current_prefactor(device.mass)),
            'gauss_quad_relative_difference':float(abs(integrate_current(device,v)['J_A_m2']/ref_j-1)),
            'refined_peak_V':float(peak[0]),'refined_valley_V':float(valley[0]),'refined_PVR':float(peak[1]/valley[1])}


def sampling_failure(out,device):
    # Thicker barriers produce a narrower zero-bias resonance, even though Tpeak≈1.
    ds=[replace(device,barrier_nm=w) for w in [2.,3.,4.,5.]]
    res=np.array([(d.barrier_nm,*first_resonance(d)) for d in ds])
    d=ds[-1]; ep,tp,width=first_resonance(d)
    # Finite window centered on the actual resonance: quantity is spectral area,
    # not the full RTD current. Adaptive scalar quadrature is explicitly seeded.
    low=ep-.005; high=ep+.005
    f=lambda e:float(transport(d,[e])['T'][0])
    reference,error=quad(f,low,high,points=[ep-width,ep,ep+width],epsabs=1e-13,epsrel=1e-9,limit=500)
    rows=[]
    for step in [.002,.001,.0005,.0001,.00002]:
        for offset in [0.,.25,.5,.75]:
            # Always cover the identical window; include endpoints before trapezoid.
            es=np.unique(np.r_[low,np.arange(low+offset*step,high,step),high])
            approx=trapezoid(transport(d,es)['T'],es)
            rows.append((step,offset,approx,reference,approx/reference-1))
    csv(out,'sampling_failure.csv','step_eV,offset_fraction,area_eV,reference_area_eV,relative_error',rows)
    csv(out,'barrier_linewidth.csv','barrier_nm,resonance_eV,T_peak,FWHM_eV',res)
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    dense=np.linspace(low,high,10001); ax[0].plot((dense-ep)*1e3,transport(d,dense)['T'],label='Resolved resonance')
    coarse=np.arange(low,high,.001); ax[0].plot((coarse-ep)*1e3,transport(d,coarse)['T'],'o',label='1 meV samples')
    ax[0].set(xlabel='E - Eres (meV)',ylabel='T',title=f'5 nm barriers: FWHM {width*1e3:.4f} meV'); ax[0].legend()
    arr=np.array(rows)
    for offset in [0.,.25,.5,.75]:
        a=arr[arr[:,1]==offset]; ax[1].semilogx(a[:,0]*1e3,a[:,-1]*100,'o-',label=f'offset {offset:g}')
    ax[1].axhline(0,color='k',lw=.8); ax[1].set(xlabel='Energy step (meV)',ylabel='Spectral-area error (%)',title='Grid offset changes the apparent area'); ax[1].legend(fontsize=8)
    fig.tight_layout(); save(out,'sampling_failure.png',fig)
    return {'barrier_linewidth_records':res.tolist(),'narrow_window_reference_area_eV':reference,
            'reference_error_eV':error,'coarse_1meV_error_range':arr[arr[:,0]==.001,-1].tolist()}


def run_all(out,device=None):
    d=device or Device(); out=Path(out); out.mkdir(parents=True,exist_ok=True)
    data,pv,checks=baseline(out,d)
    sensitivity_records=sensitivity(out,d); verify=verification(out,d); failure=sampling_failure(out,d)
    # Eighth figure collects controlled numerical convergence, not duplicate science.
    mesh=np.array(verify['mesh_records']); q=np.loadtxt(out/'quadrature_convergence.csv',delimiter=',',skiprows=1)
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    ax[0].plot(mesh[:,0],(mesh[:,3]/mesh[-1,3]-1)*100,'o-')
    ax[0].set(xlabel='Mesh spacing (nm)',ylabel='J change vs finest (%)',title='Spatial discretization at 0.1 V')
    for v in [.1,.2]:
        a=q[np.isclose(q[:,1],v)]; ax[1].plot(a[:,0]*1e3,(a[:,2]/a[-1,2]-1)*1e6,'o-',label=f'{v:g} V')
    ax[1].set(xlabel='Initial Gauss panel (meV)',ylabel='J change vs smallest (ppm)',title='Independent energy refinement'); ax[1].legend()
    fig.tight_layout(); save(out,'convergence.png',fig)
    summary={'device':asdict(d),'model':'coherent constant-mass single-valley prescribed potential; no Poisson feedback',
             'baseline_5mV':pv,'well_sensitivity':sensitivity_records,'checks':checks,'validation':verify,'sampling':failure}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
