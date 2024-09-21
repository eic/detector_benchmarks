import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, awkward as ak, sys
import mplhep as hep
hep.style.use("CMS")

plt.rcParams['figure.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'
plt.rcParams['savefig.bbox']='tight'


outdir=sys.argv[1]+"/"
config=outdir.split("/")[1]
try:
    import os
    os.mkdir(outdir[:-1])
except:
    pass
    
def gauss(x, A,mu, sigma):
    return A * np.exp(-(x-mu)**2/(2*sigma**2))
    
#load the files
import uproot as ur
arrays_sim={}
momenta=20, 30, 50, 70, 100, 150, 200, 275
for p in momenta:
    filename=f'sim_output/photon_zdc/{config}_rec_photon_zdc_{p}GeV.edm4hep.root'
    print("opening file", filename)
    events = ur.open(filename+':events')
    arrays_sim[p] = events.arrays()#[:-1] #remove last event, which for some reason is blank
    import gc
    gc.collect()
    print("read", filename)
    
fig,axs=plt.subplots(1,3, figsize=(24, 8))
pvals=[]
resvals=[]
dresvals=[]
scalevals=[]
dscalevals=[]
for p in momenta:
    selection=[len(arrays_sim[p]["HcalFarForwardZDCClusters.energy"][i])==1 for i in range(len(arrays_sim[p]))]
    E=arrays_sim[p][selection]["HcalFarForwardZDCClusters.energy"]
    
    Etot=np.sum(E, axis=-1)
    #print(len(Etot))
    #print(p, res, mrecon)
    if p==100:
        plt.sca(axs[0])
        y, x, _=plt.hist(Etot, bins=100, range=(p*.75, p*1.25), histtype='step')
        plt.ylabel("events")
        plt.title(f"$p_{{\gamma}}$={p} GeV")
        plt.xlabel("$E^{\\gamma}_{recon}$ [GeV]")
    else:
        y, x = np.histogram(Etot, bins=100, range=(p*.75, p*1.25))
        
    bc=(x[1:]+x[:-1])/2
    from scipy.optimize import curve_fit
    slc=abs(bc-p)<10
    fnc=gauss
    p0=[100, p, 10]
    #print(list(y), list(x))
    coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                 sigma=list(np.sqrt(y[slc])+(y[slc]==0)))
    if p==100:
        xx=np.linspace(p*0.75,p*1.25, 100)
        plt.plot(xx, fnc(xx,*coeff))
    pvals.append(p)
    resvals.append(np.abs(coeff[2])/coeff[1])
    dresvals.append(np.sqrt(var_matrix[2][2])/coeff[1])
    scalevals.append(np.abs(coeff[1])/p)
    dscalevals.append(np.sqrt(var_matrix[2][2])/p)
    
plt.sca(axs[1])
plt.errorbar(pvals, resvals, dresvals, ls='', marker='o')
plt.ylim(0)
plt.ylabel("$\\sigma[E_{\\gamma}]/\\mu[E_{\\gamma}]$")
plt.xlabel("$p_{\\gamma}$ [GeV]")

xx=np.linspace(15, 275, 100)

fnc=lambda E,a: a/np.sqrt(E)
#pvals, resvals, dresvals
coeff, var_matrix = curve_fit(fnc, pvals, resvals, p0=(1,),
                                 sigma=dresvals)

xx=np.linspace(15, 275, 100)
plt.plot(xx, fnc(xx, *coeff), label=f'fit:  $\\frac{{{coeff[0]*100:.0f}\\%}}{{\\sqrt{{E}}}}$')
plt.legend()
plt.sca(axs[2])
plt.errorbar(pvals, scalevals, dscalevals, ls='', marker='o')
plt.ylim(0.8, 1.2)
plt.ylabel("$\\mu[E_{\\gamma}]/E_{\\gamma}$")
plt.xlabel("$p_{\\gamma}$ [GeV]")
plt.axhline(1, ls='--', alpha=0.7, color='0.5')
plt.tight_layout()
plt.savefig(outdir+"photon_energy_res.pdf")

#theta res
fig,axs=plt.subplots(1,2, figsize=(16, 8))
pvals=[]
resvals=[]
dresvals=[]
for p in momenta:
    selection=[len(arrays_sim[p]["HcalFarForwardZDCClusters.energy"][i])==1 for i in range(len(arrays_sim[p]))]
    x=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.x"][::,0]
    y=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.y"][::,0]
    z=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.z"][::,0]
    
    theta_recon=np.arctan2(np.hypot(x*np.cos(-.025)-z*np.sin(-.025), y), z*np.cos(-.025)+x*np.sin(-.025))
    
    px=arrays_sim[p][selection]["MCParticles.momentum.x"][::,2]
    py=arrays_sim[p][selection]["MCParticles.momentum.y"][::,2]
    pz=arrays_sim[p][selection]["MCParticles.momentum.z"][::,2]

    theta_truth=np.arctan2(np.hypot(px*np.cos(-.025)-pz*np.sin(-.025), py), pz*np.cos(-.025)+px*np.sin(-.025))
    
    Etot=np.sum(E, axis=-1)
    #print(p, res, mrecon)
    if p==100:
        plt.sca(axs[0])
        y, x, _=plt.hist(1000*(theta_recon-theta_truth), bins=100, range=(-0.5, 0.5), histtype='step')
        plt.ylabel("events")
        plt.title(f"$p_{{\gamma}}$={p} GeV")
        plt.xlabel("$\\theta^{\\gamma}_{recon}$ [mrad]")
    else:
        y, x = np.histogram(1000*(theta_recon-theta_truth), bins=100, range=(-0.5, 0.5))
        
    bc=(x[1:]+x[:-1])/2
    from scipy.optimize import curve_fit
    slc=abs(bc)<0.2#1.5*np.std(1000*(theta_recon-theta_truth))
    fnc=gauss
    p0=[100, 0, 0.1]
    #print(list(y), list(x))
    coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                 sigma=list(np.sqrt(y[slc])+(y[slc]==0)))
    if p==100:
        xx=np.linspace(-0.5,0.5, 100)
        plt.plot(xx, fnc(xx,*coeff))
    pvals.append(p)
    resvals.append(np.abs(coeff[2]))
    dresvals.append(np.sqrt(var_matrix[2][2]))
plt.sca(axs[1])
plt.errorbar(pvals, resvals, dresvals, ls='', marker='o')
#print(dresvals)

fnc=lambda E,a, b: np.hypot(a/np.sqrt(E), b)
#pvals, resvals, dresvals
coeff, var_matrix = curve_fit(fnc, pvals, resvals, p0=(1,.1),
                                 sigma=dresvals)

xx=np.linspace(15, 275, 100)

plt.plot(xx, fnc(xx, *coeff), label=f'fit:  $\\frac{{{coeff[0]:.2f}}}{{\\sqrt{{E}}}}\\oplus {coeff[1]:.3f}$ mrad')

plt.ylabel("$\\sigma[\\theta_{\\gamma}]$ [mrad]")
plt.xlabel("$p_{\\gamma}$ [GeV]")

plt.ylim(0, 0.1)
plt.legend()
plt.tight_layout()
plt.savefig(outdir+"photon_theta_res.pdf")
