import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, awkward as ak, sys
import mplhep as hep
hep.style.use("CMS")

plt.rcParams['figure.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'
plt.rcParams['savefig.bbox']='tight'

plt.rcParams["figure.figsize"] = (7, 7)

config=sys.argv[1].split("/")[1]  #results/{config}/{benchmark_name}
outdir=sys.argv[1]+"/"
try:
    import os
    os.mkdir(outdir[:-1])
except:
    pass

import uproot as ur
arrays_sim={}
momenta=20, 30, 40, 50, 60,80,100
for p in momenta:
    arrays_sim[p] = ur.concatenate({
        f'sim_output/insert_tau/{config}_rec_tau-_{p}GeV_{index}.edm4eic.root': 'events'
        for index in range(5)
    })


for a in arrays_sim.values():
    #recon
    Etot=0
    px=0
    py=0
    pz=0
    for det in "HcalEndcapPInsert", "EcalEndcapP", "LFHCAL":
    
        E=a[f'{det}Clusters.energy']
        
        #todo apply corrections depending on whether this is an electromagnetic or hadronic shower.  
        
        x=a[f'{det}Clusters.position.x']
        y=a[f'{det}Clusters.position.y']
        z=a[f'{det}Clusters.position.z']
        r=np.sqrt(x**2+y**2+z**2)
        Etot=Etot+np.sum(E, axis=-1)
        px=px+np.sum(E*x/r,axis=-1)
        py=py+np.sum(E*y/r,axis=-1)
        pz=pz+np.sum(E*z/r,axis=-1)
    
    a['jet_p_recon']=np.sqrt(px**2+py**2+pz**2)
    a['jet_E_recon']=Etot
    
    a['jet_theta_recon']=np.arctan2(np.hypot(px*np.cos(-.025)-pz*np.sin(-.025),py), 
                                    pz*np.cos(-.025)+px*np.sin(-.025))
    
    #truth
    m=a['MCParticles.mass'][::,2:]
    px=a['MCParticles.momentum.x'][::,2:]
    py=a['MCParticles.momentum.y'][::,2:]
    pz=a['MCParticles.momentum.z'][::,2:]
    E=np.sqrt(m**2+px**2+py**2+pz**2)
    status=a['MCParticles.simulatorStatus'][::,2:]
    PDG=a['MCParticles.PDG'][::,2:]
    
    #find the hadronic part: initial-state tau - all leptons
    selection=1*(PDG==15)-1*(np.abs(PDG)==16)
    is_hadronic=1*(np.sum((PDG==-14)+(PDG==-12), axis=-1)==0)
    
    px_hfs, py_hfs, pz_hfs= np.sum(px*selection,axis=-1)*is_hadronic,np.sum(py*selection,axis=-1)*is_hadronic, np.sum(pz*selection,axis=-1)*is_hadronic
    
    a['hfs_p_truth']=np.sqrt(px_hfs**2+py_hfs**2+pz_hfs**2)
    a['hfs_E_truth']=np.sum(E*selection,axis=-1)*is_hadronic
    
    
    a['hfs_theta_truth']=np.arctan2(np.hypot(px_hfs*np.cos(-.025)-pz_hfs*np.sin(-.025),py_hfs), 
                                    pz_hfs*np.cos(-.025)+px_hfs*np.sin(-.025))
    a['hfs_eta_truth']=-np.log(np.tan(a['hfs_theta_truth']/2))
    a['n_mu']=np.sum(np.abs(PDG)==13, axis=-1)
    a['n_e']=np.sum(np.abs(PDG)==13, axis=-1)
    a['hfs_mass_truth']=np.sqrt(a['hfs_E_truth']**2-a['hfs_p_truth']**2)

for a in arrays_sim.values():
    selection=(a['hfs_eta_truth']>3.1) & (a['hfs_eta_truth']<3.8)\
            &(a['n_mu']==0)&(a['n_e']==0)&(a['hfs_mass_truth']>.140)&(a['jet_E_recon']>0)

Etruth=[]
Erecon=[]

theta_truth=[]
theta_recon=[]

eta_max=3.7
eta_min=3.3
for a in arrays_sim.values():
    selection=(a['hfs_eta_truth']>eta_min) & (a['hfs_eta_truth']<eta_max)\
            &(a['n_mu']==0)&(a['n_e']==0)&(a['hfs_mass_truth']>.140)&(a['jet_E_recon']>1)
    theta_truth=np.concatenate((theta_truth,a['hfs_theta_truth'][selection]))
    theta_recon=np.concatenate((theta_recon,a['jet_theta_recon'][selection]))
    Etruth=np.concatenate((Etruth,a['hfs_E_truth'][selection]))
    Erecon=np.concatenate((Erecon,a['jet_E_recon'][selection]))

plt.figure()
plt.scatter(theta_truth, theta_recon, 1)
plt.xlabel("$\\theta^{hfs}_{\\rm truth}$ [rad]")
plt.ylabel("$\\theta^{hfs}_{\\rm rec}$ [rad]")
plt.title(f"$E_{{\\tau}}$=20-100 GeV, ${eta_min}<\\eta_{{hfs}}<{eta_max}$")
plt.plot([0.04,0.1], [0.04, 0.1], color='tab:orange')
plt.ylim(0, 0.15)
plt.savefig(outdir+"/theta_scatter.pdf")

plt.figure()
plt.scatter(Etruth, Erecon, 1)
plt.xlabel("$E^{hfs}_{\\rm truth}$ [GeV]")
plt.ylabel("$E^{hfs}_{\\rm rec}$ [GeV]")
plt.title(f"$E_{{\\tau}}$=20-100 GeV, ${eta_min}<\\eta_{{hfs}}<{eta_max}$")
plt.plot((0,100), (0, 100), color='tab:orange')
plt.savefig(outdir+"/energy_scatter.pdf")

def gauss(x, A,mu, sigma):
    return A * np.exp(-(x-mu)**2/(2*sigma**2))
from scipy.optimize import curve_fit

res=[]
dres=[]
emid=[]
ew=[]
partitions=(20,30, 40, 60,80,100)
for emin, emax in zip(partitions[:-1], partitions[1:]):

    y,x = np.histogram((theta_recon-theta_truth)[(Etruth>emin)&(Etruth<emax)], bins=100, range=(-0.03,0.03))
    bc=(x[1:]+x[:-1])/2
    slc=abs(bc)<0.5
    # try:
    p0=(100, 0, 0.15)
    fnc=gauss
    sigma=np.sqrt(y[slc])+(y[slc]==0)

    coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0, sigma=list(sigma), maxfev=10000)
    res.append(abs(coeff[2]))
    dres.append(np.sqrt(var_matrix[2][2]))
    emid.append((emin+emax)/2)
    ew.append((emax-emin)/2)
plt.errorbar(emid, 1000*np.array(res),1000*np.array(dres), ew, ls='', label=f'{eta_min}<$\\eta_{{hfs}}$<{eta_max}')
plt.xlabel('$E_{hfs}$ [GeV]')
plt.ylabel('$\\theta$ resolution [mrad]')
plt.ylim(0)

fnc=lambda E,B:B/E
p0=[1,]
coeff, var_matrix = curve_fit(fnc, emid, res, p0=p0, sigma=list(dres), maxfev=10000)
xx=np.linspace(10, 100, 100)
plt.plot(xx, 1000*fnc(xx, *coeff), label=f"fit: ${coeff[0]:.2f}/E$ mrad")
plt.legend()
plt.savefig(outdir+"/theta_res.pdf")
