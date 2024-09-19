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
    filename=f'sim_output/zdc_sigma/{config}_rec_sigma_dec_{p}GeV.edm4hep.root'
    print("opening file", filename)
    events = ur.open(filename+':events')
    arrays_sim[p] = events.arrays()[:-1] #remove last event, which for some reason is blank
    import gc
    gc.collect()
    print("read", filename)

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
    plt.title(f"$p_\Sigma={p}$ GeV")
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

#create an array with the same shape as the cluster-level arrays
is_neutron_cand={}
for p in momenta:
    is_neutron_cand[p]=(0*arrays_sim[p][f"HcalFarForwardZDCClusters.energy"]).to_list()
    
    #largest_eigenvalues
    for i in range(len(arrays_sim[p])):
        pars=arrays_sim[p]["_HcalFarForwardZDCClusters_shapeParameters"][i]
        index_of_max=-1
        max_val=0
        eigs=[]
        #Must make sure this doesn't get messed up if someone changes the number of shape parameters in EICrecon.
        nClust=nclusters[p][i]
        nShapePars=len(pars)//nClust
        for j in range(nClust):
            largest_eigenvalue=max(pars[nShapePars*j+4:nShapePars*j+7])
            eigs.append(largest_eigenvalue)
            if(largest_eigenvalue>max_val):
                max_val=largest_eigenvalue
                index_of_max=j
        if index_of_max >=0:
            is_neutron_cand[p][i][index_of_max]=1
        eigs.sort()
        
    is_neutron_cand[p]=ak.Array(is_neutron_cand[p])

import ROOT

lambda_mass=1.115683
pi0_mass=0.1349768
pt_recon_corr={}
theta_recon_corr={}
mass_recon_corr={}
mass_lambda_recon_corr={}
mass_pi0_recon_corr={}
pi0_converged={}
zvtx_recon={}

#run this event-by-event:
maxZ=35800
for p in momenta:
    pt_recon_corr[p]=[]
    theta_recon_corr[p]=[]
    mass_recon_corr[p]=[]
    mass_lambda_recon_corr[p]=[]
    mass_pi0_recon_corr[p]=[]
    zvtx_recon[p]=[]
    for evt in range(len(arrays_sim[p])):
        if nclusters[p][evt]!=4:
            nan=-1
            pt_recon_corr[p].append(nan)
            theta_recon_corr[p].append(nan)
            mass_recon_corr[p].append(nan)
            mass_lambda_recon_corr[p].append(nan)
            mass_pi0_recon_corr[p].append(nan)
            zvtx_recon[p].append(nan)
            continue
        xc=arrays_sim[p][f"HcalFarForwardZDCClusters.position.x"][evt]
        yc=arrays_sim[p][f"HcalFarForwardZDCClusters.position.y"][evt]
        zc=arrays_sim[p][f"HcalFarForwardZDCClusters.position.z"][evt]
        E=arrays_sim[p][f"HcalFarForwardZDCClusters.energy"][evt]
        
        #apply correction to the neutron candidates only
        A,B,C=-0.0756, -1.91,  2.30
        neutron_corr=(1-is_neutron_cand[p][evt])+is_neutron_cand[p][evt]/(1+A+B/np.sqrt(E)+C/E)
        E=E*neutron_corr

        pabs=np.sqrt(E**2-is_neutron_cand[p][evt]*.9406**2)
        tilt=-0.025
        xcp=xc*np.cos(tilt)-zc*np.sin(tilt)
        ycp=yc
        zcp=zc*np.cos(tilt)+xc*np.sin(tilt)
        
        #search for the combination of photons that would give the best lambda mass
        pt_best=-999
        theta_best=-999
        mass_lambda_best=-999
        mass_sigma_best=-999
        mass_pi0_best=-999
        zvtx_best=-999
        for hypothesis in range(4):
            if is_neutron_cand[p][evt][hypothesis]:
                continue
            
            xvtx=0
            yvtx=0
            zvtx=0
            #find the vertex position that reconstructs the pi0 mass
            for iteration in range(20):
                tot=ROOT.TLorentzVector(0,0,0,0)
                Lambda=ROOT.TLorentzVector(0,0,0,0)
                pi0=ROOT.TLorentzVector(0,0,0,0)

                for i in range(4):

                    if i!=hypothesis:
                        ux=xcp[i]-xvtx
                        uy=ycp[i]-yvtx
                        uz=zcp[i]-zvtx
                    else:
                        ux=xcp[i]
                        uy=ycp[i]
                        uz=zcp[i]
                    u=np.sqrt(ux**2+uy**2+uz**2)
                    ux/=u
                    uy/=u
                    uz/=u

                    P=ROOT.TLorentzVector(pabs[i]*ux, pabs[i]*uy, pabs[i]*uz, E[i])
                    tot+=P
                    if not is_neutron_cand[p][evt][i] and i!=hypothesis:
                        pi0+=P
                    if i!=hypothesis:
                        Lambda+=P
                alpha=1
                if iteration==0:
                    zeta=1/2
                    zvtx=maxZ*np.power(zeta,alpha)
                    xvtx=Lambda.X()/Lambda.Z()*zvtx
                    yvtx=Lambda.Y()/Lambda.Z()*zvtx
                else :
                    s=2*(pi0.M()<pi0_mass)-1
                    zeta=np.power(zvtx/maxZ, 1/alpha)
                    zeta=zeta+s*1/2**(1+iteration)
                    zvtx=maxZ*np.power(zeta,alpha)
                    xvtx=Lambda.X()/Lambda.Z()*zvtx
                    yvtx=Lambda.Y()/Lambda.Z()*zvtx

            if abs(Lambda.M()-lambda_mass)< abs(mass_lambda_best-lambda_mass):
                pt_best=tot.Pt()
                theta_best=tot.Theta()
                mass_lambda_best=Lambda.M()
                mass_sigma_best=tot.M()
                mass_pi0_best=pi0.M()
                zvtx_best=zvtx
                
        pt_recon_corr[p].append(pt_best)
        theta_recon_corr[p].append(theta_best)
        mass_recon_corr[p].append(mass_sigma_best)
        mass_lambda_recon_corr[p].append(mass_lambda_best)
        mass_pi0_recon_corr[p].append(mass_pi0_best)
        zvtx_recon[p].append(zvtx_best)
    pt_recon_corr[p]=ak.Array(pt_recon_corr[p])
    theta_recon_corr[p]=ak.Array(theta_recon_corr[p])
    mass_recon_corr[p]=ak.Array(mass_recon_corr[p])
    mass_lambda_recon_corr[p]=ak.Array(mass_lambda_recon_corr[p])
    mass_pi0_recon_corr[p]=ak.Array(mass_pi0_recon_corr[p])
    zvtx_recon[p]=ak.Array(zvtx_recon[p])
        
#now make plots

#reconstructed vertex position plot
fig,axs=plt.subplots(1,3, figsize=(24, 8))
plt.sca(axs[0])
plt.title(f"$E_{{\\Sigma}}=100-275$ GeV")
x=[]
y=[]
for p in momenta:
    accept=(nclusters[p]==4)# &(pi0_converged[p])
    x+=list(theta_truth[p][accept]*1000)
    y+=list(theta_recon_corr[p][accept]*1000)
#print(x)
plt.scatter(x,y)
plt.xlabel("$\\theta^{*\\rm truth}_{\\Sigma}$ [mrad]")
plt.ylabel("$\\theta^{*\\rm recon}_{\\Sigma}$ [mrad]")
plt.xlim(0,3.2)
plt.ylim(0,3.2)

plt.sca(axs[1])
plt.title(f"$E_{{\\Sigma}}=100-275$ GeV")
y,x,_=plt.hist(y-np.array(x), bins=50, range=(-1,1))
bc=(x[1:]+x[:-1])/2

from scipy.optimize import curve_fit
slc=abs(bc)<0.6
fnc=gauss
p0=[100, 0, 0.5]
coeff, var_matrix = curve_fit(fnc, bc[slc], y[slc], p0=p0,
                                 sigma=np.sqrt(y[slc])+(y[slc]==0))
x=np.linspace(-1, 1)
plt.plot(x, gauss(x, *coeff), color='tab:orange')
plt.xlabel("$\\theta^{*\\rm recon}_{\\Sigma}-\\theta^{*\\rm truth}_{\\Sigma}$ [mrad]")
plt.ylabel("events")

plt.sca(axs[2])
sigmas=[]
dsigmas=[]
xvals=[]
for p in momenta:
    
    accept=(nclusters[p]==4)
    y,x=np.histogram((theta_recon_corr[p]-theta_truth[p])[accept]*1000, bins=100, range=(-0.5,0.5))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc)<0.3
    fnc=gauss
    p0=(100, 0, 0.15)
    #print(bc[slc],y[slc])
    sigma=np.sqrt(y[slc])+(y[slc]==0)
    try:
        coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
        sigmas.append(np.abs(coeff[2]))
        dsigmas.append(np.sqrt(var_matrix[2][2]))
        xvals.append(p)
    except:
        print(f"fit failed for p={p}")
print(xvals)
print(sigmas)
plt.ylim(0, 0.3)

plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', color='k')
x=np.linspace(100, 275, 100)
plt.plot(x, 3/np.sqrt(x), color='tab:orange')
plt.text(170, .23, "YR requirement:\n   3 mrad/$\\sqrt{E}$")
plt.xlabel("$E_{\\Sigma}$ [GeV]")
plt.ylabel("$\\sigma[\\theta^*_{\\Sigma}]$ [mrad]")
plt.tight_layout()
plt.savefig(outdir+"thetastar_recon.pdf")

#reconstructed vertex position plot
fig,axs=plt.subplots(1,3, figsize=(24, 8))
plt.sca(axs[0])
plt.title(f"$E_{{\\Sigma}}=100-275$ GeV")
x=[]
y=[]
for p in momenta:
    accept=(nclusters[p]==4)&(abs(mass_pi0_recon_corr[p]-pi0_mass)<.01)
    tilt=-0.025
    x+=list(arrays_sim[p]['MCParticles.vertex.z'][:,5][accept]*np.cos(tilt)/1000
            +np.sin(tilt)*arrays_sim[p]['MCParticles.vertex.z'][:,5][accept]/1000)
    y+=list(zvtx_recon[p][accept]/1000)
plt.scatter(x,y)
#print(x,y)
plt.xlabel("$z^{\\rm truth}_{\\rm vtx}$ [m]")
plt.ylabel("$z^{\\rm recon}_{\\rm vtx}$  [m]")
plt.xlim(0,40)
plt.ylim(0,40)

plt.sca(axs[1])
plt.title(f"$E_{{\\Sigma}}=100-275$ GeV")
y,x,_=plt.hist(y-np.array(x), bins=50, range=(-10,10))
bc=(x[1:]+x[:-1])/2

from scipy.optimize import curve_fit
slc=abs(bc)<5
fnc=gauss
p0=[100, 0, 1]
coeff, var_matrix = curve_fit(fnc, bc[slc], y[slc], p0=p0,
                                 sigma=np.sqrt(y[slc])+(y[slc]==0))
x=np.linspace(-5, 5)
plt.plot(x, gauss(x, *coeff), color='tab:orange')
plt.xlabel("$z^{*\\rm recon}_{\\rm vtx}-z^{*\\rm truth}_{\\rm vtx}$ [m]")
plt.ylabel("events")

plt.sca(axs[2])
sigmas=[]
dsigmas=[]
xvals=[]
for p in momenta:
    
    accept=(nclusters[p]==4)&(abs(mass_pi0_recon_corr[p]-pi0_mass)<.01)
    a=list((zvtx_recon[p]-arrays_sim[p]['MCParticles.vertex.z'][:,5])[accept]/1000)
    y,x=np.histogram(a, bins=100, range=(-10,10))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc)<5
    fnc=gauss
    p0=(100, 0, 1)
    #print(bc[slc],y[slc])
    sigma=np.sqrt(y[slc])+(y[slc]==0)
    try:
        coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
        sigmas.append(abs(coeff[2]))
        dsigmas.append(np.sqrt(var_matrix[2][2]))
        xvals.append(p)
    except:
        print(f"fit failed for p={p}")
plt.ylim(0, 3)

plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', color='k')
x=np.linspace(100, 275, 100)

avg=np.sum(sigmas/np.array(dsigmas)**2)/np.sum(1/np.array(dsigmas)**2)
plt.axhline(avg, color='tab:orange')
plt.text(150, 1.25,f"$\\sigma\\approx${avg:.1f} m")

plt.xlabel("$E_{\\Sigma}$ [GeV]")
plt.ylabel("$\\sigma[z_{\\rm vtx}]$ [m]")
plt.tight_layout()
plt.savefig(outdir+"zvtx_recon.pdf")
        
#lambda mass reconstruction
fig,axs=plt.subplots(1,2, figsize=(16, 8))
plt.sca(axs[0])
lambda_mass=1.115683
vals=[]
for p in momenta:
    accept=(nclusters[p]==4)&(abs(mass_pi0_recon_corr[p]-pi0_mass)<.01)
    vals+=list(mass_lambda_recon_corr[p][accept])

y,x,_= plt.hist(vals, bins=100, range=(0.9, 1.3))
bc=(x[1:]+x[:-1])/2
plt.axvline(lambda_mass, ls='--', color='tab:green', lw=3)
plt.text(lambda_mass+.01, np.max(y)*1.05, "PDG mass", color='tab:green')
plt.xlabel("$m_{\\Lambda}^{\\rm recon}$ [GeV]")
plt.ylim(0, np.max(y)*1.2)
plt.xlim(0.9, 1.3)

from scipy.optimize import curve_fit
slc=abs(bc-lambda_mass)<0.05
fnc=gauss
p0=[100, lambda_mass, 0.03]
coeff, var_matrix = curve_fit(fnc, bc[slc], y[slc], p0=p0,
                                 sigma=np.sqrt(y[slc])+(y[slc]==0))
x=np.linspace(0.8, 1.3, 200)
plt.plot(x, gauss(x, *coeff), color='tab:orange')
print(coeff[2], np.sqrt(var_matrix[2][2]))
plt.xlabel("$m^{\\rm recon}_{\\Lambda}$ [GeV]")
plt.ylabel("events")
plt.title(f"$E_{{\\Sigma}}=100-275$ GeV")

plt.sca(axs[1])
xvals=[]
sigmas=[]
dsigmas=[]
for p in momenta:
    accept=(nclusters[p]==4)&(abs(mass_pi0_recon_corr[p]-pi0_mass)<.01)
    y,x= np.histogram(mass_lambda_recon_corr[p][accept], bins=100, range=(0.6,1.4))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc-lambda_mass)<0.05
    fnc=gauss
    p0=[100, lambda_mass, 0.03]
    try:
        coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                       sigma=list(np.sqrt(y[slc])+(y[slc]==0)))
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
plt.xlabel("$E_{\\Sigma}$ [GeV]")
plt.ylabel("$\\sigma[m_{\\Lambda}]$ [GeV]")
plt.ylim(0, 0.02)
plt.tight_layout()
plt.savefig(outdir+"lambda_mass_rec_from_sigma_decay.pdf")

#sigma mass reconstruction
p=100
fig,axs=plt.subplots(1,2, figsize=(16, 8))
plt.sca(axs[0])
sigma_mass=1.192
vals=[]
for p in momenta:
    accept=(nclusters[p]==4)&(abs(mass_pi0_recon_corr[p]-pi0_mass)<.01)
    vals+=list(mass_recon_corr[p][accept])

y,x,_= plt.hist(vals, bins=100, range=(1.0, 1.4))
bc=(x[1:]+x[:-1])/2
plt.axvline(sigma_mass, ls='--', color='tab:green', lw=3)
plt.text(sigma_mass+.01, np.max(y)*1.05, "PDG mass", color='tab:green')
plt.xlabel("$m_{\\Sigma}^{\\rm recon}$ [GeV]")
plt.ylim(0, np.max(y)*1.2)
plt.xlim(1.0, 1.45)

from scipy.optimize import curve_fit
slc=abs(bc-sigma_mass)<0.02
fnc=gauss
p0=[100, sigma_mass, 0.03]
coeff, var_matrix = curve_fit(fnc, bc[slc], y[slc], p0=p0,
                                 sigma=np.sqrt(y[slc])+(y[slc]==0))
x=np.linspace(0.8, 1.3, 200)
plt.plot(x, gauss(x, *coeff), color='tab:orange')
print(coeff[2], np.sqrt(var_matrix[2][2]))
plt.xlabel("$m^{\\rm recon}_{\\Sigma}$ [GeV]")
plt.ylabel("events")
plt.title(f"$E_{{\\Sigma}}=100-275$ GeV")

plt.sca(axs[1])
xvals=[]
sigmas=[]
dsigmas=[]
for p in momenta:
    accept=(nclusters[p]==4)&(abs(mass_pi0_recon_corr[p]-pi0_mass)<.01)
    y,x= np.histogram(mass_recon_corr[p][accept], bins=100, range=(1.0,1.4))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc-sigma_mass)<0.02
    fnc=gauss
    p0=[100, sigma_mass, 0.03]
    try:
        coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,
                                       sigma=list(np.sqrt(y[slc])+(y[slc]==0)))
        sigmas.append(abs(coeff[2]))
        dsigmas.append(np.sqrt(var_matrix[2][2]))
        xvals.append(p)
    except:
        print("fit failed")
    
plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', color='k')
avg=np.sum(sigmas/np.array(dsigmas)**2)/np.sum(1/np.array(dsigmas)**2)
plt.axhline(avg, color='tab:orange')
plt.text(150, 0.01,f"$\\sigma\\approx${avg*1000:.0f} MeV")
plt.xlabel("$E_{\\Sigma}$ [GeV]")
plt.ylabel("$\\sigma[m_{\\Sigma}]$ [GeV]")
plt.ylim(0, 0.1)
plt.tight_layout()
plt.savefig(outdir+"sigma_mass_rec.pdf")
