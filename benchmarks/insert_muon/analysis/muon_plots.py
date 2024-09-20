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
momenta=20, 30, 40, 50, 60,70, 80
for p in momenta:
    filename=f'sim_output/insert_muon/{config}_sim_mu-_{p}GeV.edm4hep.root'
    print("opening file", filename)
    events = ur.open(filename+':events')
    arrays_sim[p] = events.arrays()
    import gc
    gc.collect()
    print("read", filename)

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
                                 sigma=list(np.sqrt(y[slc])+(y[slc]==0)))
    print(coeff)
    xx=np.linspace(0,.7, 100)
    plt.plot(xx, fnc(xx,*coeff), label=f'Landau fit:\nMIP={coeff[1]*1e3:.0f}$\\pm${1e3*np.sqrt(var_matrix[1][1]):.0f} keV')
    plt.xlabel("hit energy [MeV]")
    plt.ylabel("hits")
    plt.title(f"$E_{{\\mu^-}}=${p} GeV")
    plt.legend(fontsize=20)
    plt.savefig(outdir+"/MIP.pdf")

for p in 50,:
    array=arrays_sim[p]
    bins=30
    selection=np.sum(array["HcalEndcapPInsertHits.energy"],axis=-1)>0
    h1, xedges, yedges = np.histogram2d(list(array[selection]['phi_truth']),list(array[selection]['eta_truth']), bins=bins)
    h2, xedges, yedges = np.histogram2d(list(array['phi_truth']),list(array['eta_truth']), bins=bins)

    h = h1 / h2
    pc=plt.pcolor(xedges, yedges, h.T)
    plt.xlabel("$\\phi^*$ [rad]")
    plt.ylabel("$\\eta^*$ [rad]")
    cb = plt.colorbar(pc)
    cb.set_label("acceptance")
    plt.title(f"$E_{{\\mu^-}}=${p} GeV")
    plt.savefig(outdir+"/acceptance.pdf")
