import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, awkward as ak, sys
import mplhep as hep
hep.style.use("CMS")

plt.rcParams['figure.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'
plt.rcParams['savefig.bbox']='tight'

plt.rcParams["figure.figsize"] = (7, 7)

outdir=sys.argv[1]+"/"
config=outdir.split("/")[1]
try:
    import os
    os.mkdir(outdir[:-1])
except:
    pass
import uproot as ur
arrays_sim={}
momenta=100, 125, 150, 175,200,225,250,275
for p in momenta:
    arrays_sim[p] = ur.concatenate({
        f'sim_output/zdc_lambda/{config}_rec_lambda_dec_{p}GeV_{index}.edm4eic.root': 'events'
        for index in range(5)
    })

def gauss(x, A,mu, sigma):
    return A * np.exp(-(x-mu)**2/(2*sigma**2))

#keep track of the number of clusters per event
nclusters={}

for p in momenta:
    plt.figure()
    nclusters[p]=[]
    for i in range(len(arrays_sim[p])):
        nclusters[p].append(len(arrays_sim[p]["HcalFarForwardZDCClusters.position.x"][i]))
    nclusters[p]=np.array(nclusters[p])
    plt.hist(nclusters[p],bins=20, range=(0,20))
    plt.xlabel("number of clusters")
    plt.yscale('log')
    plt.title(f"$p_\\Lambda={p}$ GeV")
    plt.ylim(1)
    plt.savefig(outdir+f"nclust_{p}GeV_recon.pdf")
    print("saved file ", outdir+f"nclust_{p}GeV_recon.pdf")



pt_truth={}
theta_truth={}

for p in momenta:
    #get the truth value of theta* and pt*
    px=arrays_sim[p]["MCParticles.momentum.x"][:,2]
    py=arrays_sim[p]["MCParticles.momentum.y"][:,2]
    pz=arrays_sim[p]["MCParticles.momentum.z"][:,2]
    tilt=-0.025
    pt_truth[p]=np.hypot(px*np.cos(tilt)-pz*np.sin(tilt), py)
    theta_truth[p]=np.arctan2(pt_truth[p],pz*np.cos(tilt)+px*np.sin(tilt))

if "ReconstructedFarForwardZDCLambdas.momentum.x" not in arrays_sim[momenta[0]].fields:
    print("ReconstructedFarForwardZDCLambdas collection is not available (needs EICrecon 1.23)")
else:
    theta_recon={}
    E_recon={}
    zvtx_recon={}
    mass_recon={}
    print(arrays_sim[p].fields)
    for p in momenta:

        px,py,pz,m=(arrays_sim[p][f"ReconstructedFarForwardZDCLambdas.{a}"] for a in "momentum.x momentum.y momentum.z mass".split())
        theta_recon[p]=np.arctan2(np.hypot(px*np.cos(tilt)-pz*np.sin(tilt), py),pz*np.cos(tilt)+px*np.sin(tilt))
        E_recon[p]=np.sqrt(px**2+py**2+pz**2+m**2)
        zvtx_recon[p]=arrays_sim[p][f"ReconstructedFarForwardZDCLambdas.referencePoint.z"]*np.cos(tilt)+arrays_sim[p][f"ReconstructedFarForwardZDCLambdas.referencePoint.x"]*np.sin(tilt)
        mass_recon[p]=m

    #theta plots
    fig,axs=plt.subplots(1,3, figsize=(24, 8))
    plt.sca(axs[0])
    plt.title(f"$E_{{\\Lambda}}=100-275$ GeV")
    x=[]
    y=[]
    import awkward as ak
    for p in momenta:
        x+=list(ak.flatten(theta_truth[p]+0*theta_recon[p])*1000)
        y+=list(ak.flatten(theta_recon[p]*1000))
    plt.scatter(x,y)
    plt.xlabel("$\\theta^{*\\rm truth}_{\\Lambda}$ [mrad]")
    plt.ylabel("$\\theta^{*\\rm recon}_{\\Lambda}$ [mrad]")
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)

    plt.sca(axs[1])
    plt.title(f"$E_{{\\Lambda}}=100-275$ GeV")
    y,x,_=plt.hist(y-np.array(x), bins=50, range=(-1,1))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc)<0.3
    fnc=gauss
    p0=[100, 0, 0.05]
    coeff, var_matrix = curve_fit(fnc, bc[slc], y[slc], p0=p0,
                                     sigma=np.sqrt(y[slc])+(y[slc]==0), maxfev=10000)
    x=np.linspace(-1, 1)
    plt.plot(x, gauss(x, *coeff), color='tab:orange')
    plt.xlabel("$\\theta^{*\\rm recon}_{\\Lambda}-\\theta^{*\\rm truth}_{\\Lambda}$ [mrad]")
    plt.ylabel("events")

    plt.sca(axs[2])
    sigmas=[]
    dsigmas=[]
    xvals=[]
    for p in momenta:
        y,x=np.histogram(ak.flatten(theta_recon[p]-theta_truth[p])*1000, bins=100, range=(-1,1))
        bc=(x[1:]+x[:-1])/2

        from scipy.optimize import curve_fit
        slc=abs(bc)<0.3
        fnc=gauss
        p0=(100, 0, 0.06)
        sigma=np.sqrt(y[slc])+(y[slc]==0)
        try:
            coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0, sigma=list(sigma), maxfev=10000)
            sigmas.append(coeff[2])
            dsigmas.append(np.sqrt(var_matrix[2][2]))
            xvals.append(p)
        except:
            print("fit failed")
    plt.ylim(0, 0.3)

    plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', color='k')
    x=np.linspace(100, 275, 100)
    plt.plot(x, 3/np.sqrt(x), color='tab:orange')
    plt.text(170, .23, "YR requirement:\n   3 mrad/$\\sqrt{E}$")
    plt.xlabel("$E_{\\Lambda}$ [GeV]")
    plt.ylabel("$\\sigma[\\theta^*_{\\Lambda}]$ [mrad]")
    plt.tight_layout()
    plt.savefig(outdir+"thetastar_recon.pdf")
    #plt.show()

#vtx z
fig,axs=plt.subplots(1,3, figsize=(24, 8))
plt.sca(axs[0])
plt.title(f"$E_{{\\Lambda}}=100-275$ GeV")
x=[]
y=[]
for p in momenta:
    x+=list(ak.flatten(arrays_sim[p]['MCParticles.vertex.z'][:,3]+0*zvtx_recon[p])/1000)
    y+=list(ak.flatten(zvtx_recon[p])/1000)
plt.scatter(x,y)
#print(x,y)
plt.xlabel("$z^{\\rm truth}_{\\rm vtx}$ [m]")
plt.ylabel("$z^{\\rm recon}_{\\rm vtx}$  [m]")
plt.xlim(0,40)
plt.ylim(0,40)

plt.sca(axs[1])
plt.title(f"$E_{{\\Lambda}}=100-275$ GeV")
y,x,_=plt.hist(y-np.array(x), bins=50, range=(-10,10))
bc=(x[1:]+x[:-1])/2

from scipy.optimize import curve_fit
slc=abs(bc)<5
fnc=gauss
p0=[100, 0, 1]
coeff, var_matrix = curve_fit(fnc, bc[slc], y[slc], p0=p0,
                                 sigma=np.sqrt(y[slc])+(y[slc]==0), maxfev=10000)
x=np.linspace(-5, 5)
plt.plot(x, gauss(x, *coeff), color='tab:orange')
print(coeff[2], np.sqrt(var_matrix[2][2]))
plt.xlabel("$z^{*\\rm recon}_{\\rm vtx}-z^{*\\rm truth}_{\\rm vtx}$ [m]")
plt.ylabel("events")

plt.sca(axs[2])
sigmas=[]
dsigmas=[]
xvals=[]
for p in momenta:
    

    a=ak.flatten((zvtx_recon[p]-arrays_sim[p]['MCParticles.vertex.z'][:,3])/1000)
    y,x=np.histogram(a, bins=100, range=(-10,10))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc)<5
    fnc=gauss
    p0=(100, 0, 1)
    #print(bc[slc],y[slc])
    sigma=np.sqrt(y[slc])+(y[slc]==0)
    try:
        coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0, sigma=list(sigma), maxfev=10000)
        sigmas.append(coeff[2])
        dsigmas.append(np.sqrt(var_matrix[2][2]))
        xvals.append(p)
    except:
        print("fit failed")
plt.ylim(0, 2)

plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', color='k')
x=np.linspace(100, 275, 100)

avg=np.sum(sigmas/np.array(dsigmas)**2)/np.sum(1/np.array(dsigmas)**2)
plt.axhline(avg, color='tab:orange')
plt.text(150, 1.25,f"$\\sigma\\approx${avg:.1f} m")

plt.xlabel("$E_{\\Lambda}$ [GeV]")
plt.ylabel("$\\sigma[z_{\\rm vtx}]$ [m]")
plt.tight_layout()
plt.savefig(outdir+"zvtx_recon.pdf")
#plt.show()

p=100
fig,axs=plt.subplots(1,2, figsize=(16, 8))
plt.sca(axs[0])
lambda_mass=1.115683
vals=[]
for p in momenta:
    vals+=list(ak.flatten(mass_recon[p]))

y,x,_= plt.hist(vals, bins=100, range=(1.0, 1.25))
bc=(x[1:]+x[:-1])/2
plt.axvline(lambda_mass, ls='--', color='tab:green', lw=3)
plt.text(lambda_mass+.01, np.max(y)*1.05, "PDG mass", color='tab:green')
plt.xlabel("$m_{\\Lambda}^{\\rm recon}$ [GeV]")
plt.ylim(0, np.max(y)*1.2)
plt.xlim(1.0, 1.25)

from scipy.optimize import curve_fit
slc=abs(bc-lambda_mass)<0.05
fnc=gauss
p0=[100, lambda_mass, 0.04]
coeff, var_matrix = curve_fit(fnc, bc[slc], y[slc], p0=p0,
                                 sigma=np.sqrt(y[slc])+(y[slc]==0), maxfev=10000)
x=np.linspace(0.8, 1.3, 200)
plt.plot(x, gauss(x, *coeff), color='tab:orange')
print(coeff[2], np.sqrt(var_matrix[2][2]))
plt.xlabel("$m^{\\rm recon}_{\\Lambda}$ [GeV]")
plt.ylabel("events")
plt.title(f"$E_{{\\Lambda}}=100-275$ GeV")

plt.sca(axs[1])
xvals=[]
sigmas=[]
dsigmas=[]
for p in momenta:
    y,x= np.histogram(ak.flatten(mass_recon[p]), bins=100, range=(0.6,1.4))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc-lambda_mass)<0.05
    fnc=gauss
    p0=[100, lambda_mass, 0.05]
    try:
        coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                       sigma=list(np.sqrt(y[slc])+(y[slc]==0)), maxfev=10000)
        x=np.linspace(0.8, 1.3, 200)
        sigmas.append(coeff[2])
        dsigmas.append(np.sqrt(var_matrix[2][2]))
        xvals.append(p)
    except:
        print("fit failed")
    
plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', color='k')
avg=np.sum(sigmas/np.array(dsigmas)**2)/np.sum(1/np.array(dsigmas)**2)
plt.axhline(avg, color='tab:orange')
plt.text(150, 0.01,f"$\\sigma\\approx${avg*1000:.0f} MeV")
plt.xlabel("$E_{\\Lambda}$ [GeV]")
plt.ylabel("$\\sigma[m_{\\Lambda}]$ [GeV]")
plt.ylim(0, 0.02)
plt.tight_layout()
plt.savefig(outdir+"lambda_mass_rec.pdf")


#now for the CM stuff:
phi_residuals=[]
theta_residuals=[]
for p in momenta:
    isNeutron=arrays_sim[p]['ReconstructedFarForwardZDCLambdaDecayProductsCM.PDG']==2112
    pxcm=arrays_sim[p]['ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.x']
    pycm=arrays_sim[p]['ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.y']
    pzcm=arrays_sim[p]['ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.z']


    import ROOT
    px,py,pz,E=arrays_sim[p]['MCParticles.momentum.x'], arrays_sim[p]['MCParticles.momentum.y'], arrays_sim[p]['MCParticles.momentum.z'],np.sqrt(arrays_sim[p]['MCParticles.momentum.x']**2+arrays_sim[p]['MCParticles.momentum.y']**2+arrays_sim[p]['MCParticles.momentum.z']**2\
                +arrays_sim[p]['MCParticles.mass']**2)
    phicmtruth=list(np.repeat(-9999, len(arrays_sim[p])))
    thetacmtruth=list(np.repeat(-9999, len(arrays_sim[p])))
    for i in range(len(arrays_sim[p])):
        l=ROOT.TLorentzVector(px[i,2], py[i,2],  pz[i,2], E[i,2])
        n=ROOT.TLorentzVector(px[i,3], py[i,3],  pz[i,3], E[i,3])
        ncm=n.Clone();
        ncm.Boost(-l.BoostVector())
        phicmtruth[i]=ncm.Phi()
        thetacmtruth[i]=ncm.Theta()
        
    arrays_sim[p]["phicmtruth"]=phicmtruth
    arrays_sim[p]["thetacmtruth"]=thetacmtruth

    phicmtruth=arrays_sim[p]["phicmtruth"]
    thetacmtruth=arrays_sim[p]["thetacmtruth"]
    phi_residuals=np.concatenate((phi_residuals, ak.flatten((np.arctan2(pycm,pxcm)[isNeutron]-phicmtruth)*np.sin(thetacmtruth))))
    theta_residuals=np.concatenate((theta_residuals, ak.flatten(np.arctan2(np.hypot(pycm,pxcm),pzcm)[isNeutron]-thetacmtruth)))
plt.figure()
plt.hist(phi_residuals*1000, bins=100, range=(-300, 300))
plt.xlabel("$(\\phi^n_{cm,rec}-\\phi^n_{cm,truth})\\times\\sin\\theta^n_{cm,truth} [mrad]$")
plt.savefig(outdir+"neutron_phi_cm_res.pdf")

plt.figure()
plt.hist(1000*theta_residuals, bins=100, range=(-1000, 1000))
plt.xlabel("$\\theta^n_{cm,rec}-\\theta^n_{cm,truth}$ [mrad]")
plt.savefig(outdir+"neutron_theta_cm_res.pdf")
