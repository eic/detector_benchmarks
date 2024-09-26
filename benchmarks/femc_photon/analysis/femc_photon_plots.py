import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, awkward as ak, sys, uproot as ur
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
arrays_sim={p:ur.open(f'sim_output/femc_photon/{config}_rec_e-_{p}GeV.edm4hep.root:events').arrays() for p in (20, 30, 40, 50, 60,70,80)}

for p in arrays_sim:
    array=arrays_sim[p]
    tilt=-.025
    px=array['MCParticles.momentum.x'][:,2]
    py=array['MCParticles.momentum.y'][:,2]
    pz=array['MCParticles.momentum.z'][:,2]
    p=np.sqrt(px**2+py**2+pz**2)
    
    pxp=px*np.cos(tilt)-pz*np.sin(tilt)
    pyp=py
    pzp=pz*np.cos(tilt)+px*np.sin(tilt)
    
    array['eta_truth']=1/2*np.log((p+pzp)/(p-pzp))
    array['nclust_endcap']=[len(array['EcalEndcapPClusters.energy'][i]) for i in range(len(array))]
    
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

#number of clusters
plt.figure()
for eta_min, eta_max, field in (1.5, 2.8, 'nclust_endcap'),:
    for p in arrays_sim:
        array=arrays_sim[p]
        plt.hist(array[field][(array['eta_truth']>eta_min)&(array['eta_truth']<eta_max)],
                 bins=np.linspace(0,10,11), histtype='step', label=f'{p} GeV', density=True)
    plt.ylabel("events")
    plt.xlabel("# of Ecal clusters")
    plt.legend()
    plt.savefig(outdir+f"/{field}.pdf")

fig, axs=plt.subplots(1,2, figsize=(16,8))
avgs=[]
stds=[]
pvals=[]

#number of hits per cluster
fig, axs=plt.subplots(1,2, figsize=(16,8))
avgs=[]
stds=[]
pvals=[]

for p in arrays_sim:

    a=arrays_sim[p]
    n=[]
    nn=-a['EcalEndcapPClusters.hits_begin']+a['EcalEndcapPClusters.hits_end']
    E=a['EcalEndcapPClusters.energy']
    for evt in range(len(array)):
        maxE=np.max(E[evt])
        found=False
        for i in range(len(E[evt])):
            if E[evt][i]==maxE:
                n.append(nn[evt][i])
                found=True
                break
        #if not found:
        #    n.append(0)
    
    if p ==50:
        plt.sca(axs[0])
        y,x,_=plt.hist(n, range=(0,100), bins=100, histtype='step', label=f"E={p} GeV")
        plt.ylabel("events")
        plt.xlabel("# hits in cluster")
        plt.title(f"e-, E={p} GeV")
    pvals.append(p)
    avgs.append(np.mean(n))
    stds.append(np.std(n))

plt.sca(axs[1])
plt.errorbar(pvals, avgs, stds, marker='o',ls='')
plt.xlabel("E [GeV]")
plt.ylabel("# hits in cluster [mean$\\pm$std]")
plt.ylim(0)
plt.savefig(outdir+"/nhits_per_cluster.pdf")


#energy resolution
def gauss(x, A,mu, sigma):
    return A * np.exp(-(x-mu)**2/(2*sigma**2))
from scipy.optimize import curve_fit

fig, axs=plt.subplots(1,3, figsize=(24,8))
pvals=[]
res=[]
dres=[]
scale=[]
dscale=[]
for p in arrays_sim:
    bins=np.linspace(15*p/20,22*p/20, 50)
    if p==50:
        plt.sca(axs[0])
        plt.title(f"E={p} GeV")
        y,x,_=plt.hist(ak.flatten(arrays_sim[p]['EcalEndcapPClusters.energy']), bins=bins, histtype='step')
        plt.ylabel("events")
        plt.xlabel("$E^{rec}_e$ [GeV]")
    else:
        y,x=np.histogram(ak.flatten(arrays_sim[p]['EcalEndcapPClusters.energy']), bins=bins)
    bcs=(x[1:]+x[:-1])/2

    fnc=gauss
    slc=abs(bcs-p)<3
    sigma=np.sqrt(y[slc])+0.5*(y[slc]==0)
    p0=(100, p, 3)

    coeff, var_matrix = curve_fit(fnc, list(bcs[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
    #res=np.abs(coeff[2]/coeff[1])
    if p==50:
        xx=np.linspace(15*p/20,22*p/20, 100)

        plt.plot(xx, fnc(xx,*coeff), label=f"$\\sigma_E/E={abs(coeff[2])/coeff[1]*100:.1f}\%$")
        plt.axvline(p, color='g', ls='--', alpha=0.7)
        plt.legend()
        #plt.xlim(0,60)
    #plt.show()
    pvals.append(p)
    res.append(abs(coeff[2])/coeff[1])
    dres.append(np.sqrt(var_matrix[2][2])/coeff[1])
    scale.append(abs(coeff[1])/p)
    dscale.append(np.sqrt(var_matrix[1][1])/p)
plt.sca(axs[1])
plt.errorbar(pvals, 100*np.array(res), 100*np.array(dres), ls='', marker='o')
fnc = lambda E, a, b: np.hypot(a,b/np.sqrt(E))
p0=(.05, .12)
coeff, var_matrix = curve_fit(fnc, pvals, res, p0=p0,sigma=dres)
xx=np.linspace(15, 85, 100)
plt.plot(xx, 100*fnc(xx,*coeff), label=f'fit:{100*coeff[0]:.0f}%$\\oplus\\frac{{{100*coeff[1]:.1f}\\%}}{{\\sqrt{{E}}}}$')
plt.legend()
plt.ylim(0)
plt.ylabel("E resolution [%]")
plt.xlabel("E truth [GeV]")
plt.sca(axs[2])

plt.errorbar(pvals, 100*np.array(scale), 100*np.array(dscale), ls='', marker='o')
plt.ylabel("energy scale [%]")
plt.xlabel("E truth [GeV]")
plt.axhline(100, color='0.5', alpha=0.5, ls='--')
plt.ylim(0, 110)
plt.tight_layout()
plt.savefig(outdir+"/energy_res.pdf")
