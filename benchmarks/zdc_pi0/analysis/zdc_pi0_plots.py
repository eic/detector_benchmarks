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
    
import uproot as ur
arrays_sim={}
momenta=60, 80, 100, 130, 160,
for p in momenta:
    arrays_sim[p] = ur.concatenate({
        f'sim_output/zdc_pi0/{config}_rec_zdc_pi0_{p}GeV_{index}.edm4eic.root': 'events'
        for index in range(5)
    })

#energy res plot
fig,axs=plt.subplots(1,3, figsize=(24, 8))
pvals=[]
resvals=[]
dresvals=[]
scalevals=[]
dscalevals=[]
for p in momenta:
    selection=[len(arrays_sim[p]["HcalFarForwardZDCClusters.energy"][i])==2 for i in range(len(arrays_sim[p]))]
    E=arrays_sim[p][selection]["HcalFarForwardZDCClusters.energy"]
    
    Etot=np.sum(E, axis=-1)
    if len(Etot)<25:
        continue
    #print(p, res, mrecon)
    if p==100:
        plt.sca(axs[0])
        y, x, _=plt.hist(Etot, bins=100, range=(p*.5, p*1.5), histtype='step')
        plt.ylabel("events")
        plt.title(f"$p_{{\pi^0}}$={p} GeV")
        plt.xlabel("$E^{\\pi^{0}}_{recon}$ [GeV]")
    else:
        y, x = np.histogram(Etot, bins=100, range=(p*.5, p*1.5))
        
    bc=(x[1:]+x[:-1])/2
    from scipy.optimize import curve_fit
    slc=abs(bc-p)<10
    fnc=gauss
    p0=[100, p, 10]
    #print(list(y), list(x))
    coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                 sigma=list(np.sqrt(y[slc])+(y[slc]==0)), maxfev=10000)
    if p==100:
        xx=np.linspace(p*0.5,p*1.5, 100)
        plt.plot(xx, fnc(xx,*coeff))
    pvals.append(p)
    resvals.append(np.abs(coeff[2])/coeff[1])
    dresvals.append(np.sqrt(var_matrix[2][2])/coeff[1])
    scalevals.append(np.abs(coeff[1])/p)
    dscalevals.append(np.sqrt(var_matrix[2][2])/p)
    
plt.sca(axs[1])
plt.errorbar(pvals, resvals, dresvals, ls='', marker='o')

plt.ylabel("$\\sigma[E_{\\pi^0}]/\\mu[E_{\\pi^0}]$")
plt.xlabel("$p_{\\pi^0}$ [GeV]")

fnc=lambda E,a: a/np.sqrt(E)
#pvals, resvals, dresvals
coeff, var_matrix = curve_fit(fnc, pvals, resvals, p0=(1,),
                                 sigma=dresvals, maxfev=10000)
xx=np.linspace(55, 200, 100)
plt.plot(xx, fnc(xx, *coeff), label=f'fit:  $\\frac{{{coeff[0]:.2f}\\%}}{{\\sqrt{{E}}}}$')
plt.legend()
plt.ylim(0)
plt.sca(axs[2])
plt.errorbar(pvals, scalevals, dscalevals, ls='', marker='o')
plt.ylim(0.8, 1.2)
plt.ylabel("$\\mu[E_{\\pi^0}]/E_{\\pi^0}$")
plt.xlabel("$p_{\\pi^0}$ [GeV]")
plt.axhline(1, ls='--', alpha=0.7, color='0.5')
plt.tight_layout()
plt.savefig(outdir+"/pi0_energy_res.pdf")


fig,axs=plt.subplots(1,2, figsize=(16, 8))
pvals=[]
resvals=[]
dresvals=[]
for p in momenta:
    selection=[len(arrays_sim[p]["HcalFarForwardZDCClusters.energy"][i])==2 for i in range(len(arrays_sim[p]))]
    x=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.x"]
    y=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.y"]
    z=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.z"]
    E=arrays_sim[p][selection]["HcalFarForwardZDCClusters.energy"]
    r=np.sqrt(x**2+y**2+z**2)
    px=np.sum(E*x/r, axis=-1)
    py=np.sum(E*y/r, axis=-1)
    pz=np.sum(E*z/r, axis=-1)
    
    theta_recon=np.arctan2(np.hypot(px*np.cos(-.025)-pz*np.sin(-.025), py), pz*np.cos(-.025)+px*np.sin(-.025))
    if len(theta_recon)<25:
        continue
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
        plt.title(f"$p_{{\\pi^0}}$={p} GeV")
        plt.xlabel("$\\theta^{\\pi^0}_{recon}$ [mrad]")
    else:
        y, x = np.histogram(1000*(theta_recon-theta_truth), bins=100, range=(-0.5, 0.5))
        
    bc=(x[1:]+x[:-1])/2
    from scipy.optimize import curve_fit
    slc=abs(bc)<0.2#1.5*np.std(1000*(theta_recon-theta_truth))
    fnc=gauss
    p0=[100, 0, 0.1]
    #print(list(y), list(x))
    coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                 sigma=list(np.sqrt(y[slc])+(y[slc]==0)), maxfev=10000)
    if p==100:
        xx=np.linspace(-0.5,0.5, 100)
        plt.plot(xx, fnc(xx,*coeff))
    pvals.append(p)
    resvals.append(np.abs(coeff[2]))
    dresvals.append(np.sqrt(var_matrix[2][2]))
    
plt.sca(axs[1])
plt.errorbar(pvals, resvals, dresvals, ls='', marker='o')
#print(dresvals)

fnc=lambda E,a: a/np.sqrt(E)
#pvals, resvals, dresvals
coeff, var_matrix = curve_fit(fnc, pvals, resvals, p0=(1,),
                                 sigma=dresvals, maxfev=10000)

xx=np.linspace(55, 200, 100)

plt.plot(xx, fnc(xx, *coeff), label=f'fit:  $\\frac{{{coeff[0]:.2f}}}{{\\sqrt{{E}}}}$ mrad')

plt.ylabel("$\\sigma[\\theta_{\\pi^0}]$ [mrad]")
plt.xlabel("$p_{\\pi^0}$ [GeV]")

plt.ylim(0, 0.1)
plt.legend()
plt.tight_layout()
plt.savefig(outdir+"/pi0_theta_res.pdf")

fig,axs=plt.subplots(1,2, figsize=(16, 8))
pvals=[]
resvals=[]
dresvals=[]
for p in momenta:
    selection=[len(arrays_sim[p]["HcalFarForwardZDCClusters.energy"][i])==2 for i in range(len(arrays_sim[p]))]
    E=arrays_sim[p][selection]["HcalFarForwardZDCClusters.energy"]
    cx=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.x"]
    cy=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.y"]
    cz=arrays_sim[p][selection]["HcalFarForwardZDCClusters.position.z"]
    r=np.sqrt(cx**2+cy**2+cz**2)
    px=E*cx/r
    py=E*cy/r
    pz=E*cz/r
    
    cos_opening_angle=(cx/r)[::,0]*(cx/r)[::,1]+(cy/r)[::,0]*(cy/r)[::,1]+(cz/r)[::,0]*(cz/r)[::,1]
    mrecon=np.sqrt(2*E[::,0]*E[::,1]*(1-cos_opening_angle))
    
    if len(mrecon)<25:
        continue
    
    #print(p, res, mrecon)
    if p==100:
        plt.sca(axs[0])
        y, x, _=plt.hist(mrecon, bins=100, range=(0, 0.2), histtype='step')
        plt.ylabel("events")
        plt.title(f"$p_{{\pi^0}}$={p} GeV")
        plt.xlabel("$m^{\\pi^{0}}_{recon}$ [GeV]")
    else:
        #y, x, _=plt.hist(mrecon, bins=100, range=(0, 0.2), histtype='step')#y, x =np.histogram(mrecon, bins=100, range=(0, 0.2))
        y, x = np.histogram(mrecon, bins=100, range=(0, 0.2))
        
    bc=(x[1:]+x[:-1])/2
    from scipy.optimize import curve_fit
    slc=abs(bc-.135)<.1
    fnc=gauss
    p0=[100, .135, 0.2]
    #print(list(y), list(x))
    coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                 sigma=list(np.sqrt(y[slc])+(y[slc]==0)), maxfev=10000)
    if p==100:
        xx=np.linspace(0,0.2)
        plt.plot(xx, fnc(xx,*coeff))
    pvals.append(p)
    resvals.append(np.abs(coeff[2]))
    dresvals.append(np.sqrt(var_matrix[2][2]))
    
plt.sca(axs[1])
plt.errorbar(pvals, resvals, dresvals, ls='', marker='o')
plt.ylim(0)
plt.ylabel("$\\sigma[m_{\\pi^0}]$ [GeV]")
plt.xlabel("$p_{\\pi^0}$ [GeV]")

fnc=lambda E,a,b: a+b*E
#pvals, resvals, dresvals
coeff, var_matrix = curve_fit(fnc, pvals, resvals, p0=(1,1),
                                 sigma=dresvals, maxfev=10000)
xx=np.linspace(55, 200, 100)
#plt.plot(xx, fnc(xx, *coeff), label=f'fit:  ${coeff[0]*1000:.1f}+{coeff[1]*1000:.4f}\\times E$ MeV')
plt.plot(xx, fnc(xx, *coeff), label=f'fit:  $({coeff[0]*1000:.1f}+{coeff[1]*1000:.4f}\\times [E\,in\,GeV])$ MeV')
plt.legend()


plt.tight_layout()
plt.savefig(outdir+"/pi0_mass_res.pdf")
