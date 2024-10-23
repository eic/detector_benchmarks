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

def Landau(x, normalization,location,stdev):
    #print(type(x))
    u=(x-location)*3.591/stdev/2.355
    renormalization = 1.64872*normalization
    return renormalization * np.exp(-u/2 - np.exp(-u)/2)

import uproot as ur
arrays_sim={}
momenta=50,
for p in momenta:
    arrays_sim[p] = ur.concatenate({
        f'sim_output/insert_muon/{config}_sim_mu-_{p}GeV_{index}.edm4hep.root': 'events'
        for index in range(5)
    })

for array in arrays_sim.values():
    tilt=-0.025
    px=array['MCParticles.momentum.x'][:,2]
    py=array['MCParticles.momentum.y'][:,2]
    pz=array['MCParticles.momentum.z'][:,2]
    p=np.sqrt(px**2+py**2+pz**2)
    
    pxp=px*np.cos(tilt)-pz*np.sin(tilt)
    pyp=py
    pzp=pz*np.cos(tilt)+px*np.sin(tilt)
    
    array['eta_truth']=1/2*np.log((p+pzp)/(p-pzp))
    array['phi_truth']=np.arctan2(pyp,pxp)
    
for p in 50,:
    E=arrays_sim[p]["HcalEndcapPInsertHits.energy"]
    y, x,_=plt.hist(1e3*ak.flatten(E),bins=100, range=(0, 1.2), histtype='step')
    bc=(x[1:]+x[:-1])/2
    from scipy.optimize import curve_fit
    slc=abs(bc-.48)<.15
    fnc=Landau
    p0=[100, .5, .05]
    #print(list(y), list(x))
    coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                 sigma=list(np.sqrt(y[slc])+(y[slc]==0)), maxfev=10000)
    print(coeff)
    xx=np.linspace(0,.7, 100)
    MIP=coeff[1]/1000
    plt.plot(xx, fnc(xx,*coeff), label=f'Landau fit:\nMIP={coeff[1]*1e3:.0f}$\\pm${1e3*np.sqrt(var_matrix[1][1]):.0f} keV')
    plt.xlabel("hit energy [MeV]")
    plt.ylabel("hits")
    plt.title(f"$E_{{\\mu^-}}=${p} GeV")
    plt.legend(fontsize=20)
    plt.savefig(outdir+"/MIP.pdf")

    plt.figure(figsize=(10,7))
    array=arrays_sim[p]
    bins=30; r=((-np.pi, np.pi),(2.8, 4.2))
    selection=np.sum(array["HcalEndcapPInsertHits.energy"]>0.5*MIP,axis=-1)>0
    h1, xedges, yedges = np.histogram2d(list(array[selection]['phi_truth']),list(array[selection]['eta_truth']), bins=bins, range=r)
    h2, xedges, yedges = np.histogram2d(list(array['phi_truth']),list(array['eta_truth']), bins=bins, range=r)

    h = h1 / h2
    pc=plt.pcolor(xedges, yedges, h.T,linewidth=0)
    plt.xlabel("$\\phi^*$ [rad]")
    plt.ylabel("$\\eta^*$")
    cb = plt.colorbar(pc)
    cb.set_label("acceptance")
    plt.title(f"$E_{{\\mu^-}}=${p} GeV")
    plt.savefig(outdir+"/acceptance.pdf")
