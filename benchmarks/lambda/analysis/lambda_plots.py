import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, awkward as ak, sys
import mplhep as hep
hep.style.use("CMS")

plt.rcParams['figure.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'
plt.rcParams['savefig.bbox']='tight'

plt.rcParams["figure.figsize"] = (7, 7)

outdir=sys.argv[1]+"/"
try:
    import os
    os.mkdir(outdir[:-1])
except:
    pass
import uproot as ur
arrays_sim={}
momenta=100, 125, 150, 175,200,225,250,275
for p in momenta:
    filename=f'sim_output/lambda/epic_zdc_sipm_on_tile_only_rec_lambda_dec_{p}GeV.edm4hep.root'
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
    plt.title(f"$p_\Lambda={p}$ GeV")
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
        shape_params_per_cluster=7
        for j in range(len(pars)//shape_params_per_cluster):
            largest_eigenvalue=max(pars[shape_params_per_cluster*j+4:shape_params_per_cluster*j+7])
            eigs.append(largest_eigenvalue)
            if(largest_eigenvalue>max_val):
                max_val=largest_eigenvalue
                index_of_max=j
        if index_of_max >=0:
            is_neutron_cand[p][i][index_of_max]=1
        eigs.sort()
        
    is_neutron_cand[p]=ak.Array(is_neutron_cand[p])


#with the position of the vertex determined by assuming the mass of the pi0
#corrected pt* and theta* recon
pt_recon_corr={}
theta_recon_corr={}
mass_recon_corr={}
mass_pi0_recon_corr={}
pi0_converged={}
zvtx_recon={}

maxZ=35800
for p in momenta:
    xvtx=0
    yvtx=0
    zvtx=0
    
    for iteration in range(20):
    
        #compute the value of theta* using the clusters in the ZDC
        xc=arrays_sim[p][f"HcalFarForwardZDCClusters.position.x"]
        yc=arrays_sim[p][f"HcalFarForwardZDCClusters.position.y"]
        zc=arrays_sim[p][f"HcalFarForwardZDCClusters.position.z"]
        E=arrays_sim[p][f"HcalFarForwardZDCClusters.energy"]
        #apply correction to the neutron candidates only
        A,B,C=-0.0756, -1.91,  2.30
        neutron_corr=(1-is_neutron_cand[p])+is_neutron_cand[p]/(1+A+B/np.sqrt(E)+C/E)
        E=E*neutron_corr

        E_recon=np.sum(E, axis=-1)
        pabs=np.sqrt(E**2-is_neutron_cand[p]*.9406**2)
        tilt=-0.025
        xcp=xc*np.cos(tilt)-zc*np.sin(tilt)
        ycp=yc
        zcp=zc*np.cos(tilt)+xc*np.sin(tilt)
        rcp=np.sqrt(xcp**2+ycp**2+zcp**2)
        
        ux=(xcp-xvtx)
        uy=(ycp-yvtx)
        uz=(zcp-zvtx)
        
        norm=np.sqrt(ux**2+uy**2+uz**2)
        ux=ux/norm
        uy=uy/norm
        uz=uz/norm
        
        px_recon,py_recon,pz_recon=np.sum(pabs*ux, axis=-1),np.sum(pabs*uy, axis=-1),np.sum(pabs*uz, axis=-1)

        pt_recon_corr[p]=np.hypot(px_recon,py_recon)
        theta_recon_corr[p]=np.arctan2(pt_recon_corr[p], pz_recon)
        
        mass_recon_corr[p]=np.sqrt((E_recon)**2\
                                -(px_recon)**2\
                                -(py_recon)**2\
                                -(pz_recon)**2)
        mass_pi0_recon_corr[p]=np.sqrt(np.sum(pabs*(1-is_neutron_cand[p]), axis=-1)**2\
                                    -np.sum(pabs*ux*(1-is_neutron_cand[p]), axis=-1)**2\
                                    -np.sum(pabs*uy*(1-is_neutron_cand[p]), axis=-1)**2\
                                    -np.sum(pabs*uz*(1-is_neutron_cand[p]), axis=-1)**2)
        alpha=1
        if iteration==0:
            u=np.sqrt(px_recon**2+py_recon**2+pz_recon**2)
            ux=px_recon/u
            uy=py_recon/u
            uz=pz_recon/u
            zeta=1/2
            zvtx=maxZ*np.power(zeta,alpha)
            xvtx=ux/uz*zvtx
            yvtx=uy/uz*zvtx
        else :
            u=np.sqrt(px_recon**2+py_recon**2+pz_recon**2)
            ux=px_recon/u
            uy=py_recon/u
            uz=pz_recon/u
            s=2*(mass_pi0_recon_corr[p]<0.135)-1
            zeta=np.power(zvtx/maxZ, 1/alpha)
            zeta=zeta+s*1/2**(1+iteration)
            zvtx=maxZ*np.power(zeta,alpha)
            xvtx=ux/uz*zvtx
            yvtx=uy/uz*zvtx
        #print(zvtx)
    pi0_converged[p]=np.abs(mass_pi0_recon_corr[p]-0.135)<0.01
    zvtx_recon[p]=zvtx
        
fig,axs=plt.subplots(1,3, figsize=(24, 8))
plt.sca(axs[0])
plt.title(f"$E_{{\\Lambda}}=100-275$ GeV")
x=[]
y=[]
for p in momenta:
    accept=(nclusters[p]==3) &(pi0_converged[p])
    x+=list(theta_truth[p][accept]*1000)
    y+=list(theta_recon_corr[p][accept]*1000)
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
                                 sigma=np.sqrt(y[slc])+(y[slc]==0))
x=np.linspace(-1, 1)
plt.plot(x, gauss(x, *coeff), color='tab:orange')
plt.xlabel("$\\theta^{*\\rm recon}_{\\Lambda}-\\theta^{*\\rm truth}_{\\Lambda}$ [mrad]")
plt.ylabel("events")

plt.sca(axs[2])
sigmas=[]
dsigmas=[]
xvals=[]
for p in momenta:
    
    accept=(nclusters[p]==3) &(pi0_converged[p])
    y,x=np.histogram((theta_recon_corr[p]-theta_truth[p])[accept]*1000, bins=100, range=(-1,1))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc)<0.3
    fnc=gauss
    p0=(100, 0, 0.06)
    #print(bc[slc],y[slc])
    sigma=np.sqrt(y[slc])+(y[slc]==0)
    try:
        coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
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


fig,axs=plt.subplots(1,3, figsize=(24, 8))
plt.sca(axs[0])
plt.title(f"$E_{{\\Lambda}}=100-275$ GeV")
x=[]
y=[]
for p in momenta:
    accept=(nclusters[p]==3) &(pi0_converged[p])
    x+=list(arrays_sim[p]['MCParticles.vertex.z'][:,3][accept]/1000)
    y+=list(zvtx_recon[p][accept]/1000)
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
                                 sigma=np.sqrt(y[slc])+(y[slc]==0))
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
    
    accept=(nclusters[p]==3) &(pi0_converged[p])
    a=list((zvtx_recon[p]-arrays_sim[p]['MCParticles.vertex.z'][:,3])[accept]/1000)
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
    accept=(nclusters[p]==3) &(pi0_converged[p])
    vals+=list(mass_recon_corr[p][accept])

y,x,_= plt.hist(vals, bins=100, range=(1.0, 1.25))
bc=(x[1:]+x[:-1])/2
plt.axvline(lambda_mass, ls='--', color='tab:green', lw=3)
plt.text(lambda_mass+.01, np.max(y)*1.05, "PDG mass", color='tab:green')
plt.xlabel("$m_{\\Lambda}^{\\rm recon}$ [GeV]")
plt.ylim(0, np.max(y)*1.2)
plt.xlim(1.0, 1.25)

from scipy.optimize import curve_fit
slc=abs(bc-lambda_mass)<0.07
fnc=gauss
p0=[100, lambda_mass, 0.04]
coeff, var_matrix = curve_fit(fnc, bc[slc], y[slc], p0=p0,
                                 sigma=np.sqrt(y[slc])+(y[slc]==0))
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
    accept=(nclusters[p]==3) &(pi0_converged[p])
    y,x= np.histogram(mass_recon_corr[p][accept], bins=100, range=(0.6,1.4))
    bc=(x[1:]+x[:-1])/2

    from scipy.optimize import curve_fit
    slc=abs(bc-lambda_mass)<0.07
    fnc=gauss
    p0=[100, lambda_mass, 0.05]
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
plt.xlabel("$E_{\\Lambda}$ [GeV]")
plt.ylabel("$\\sigma[m_{\\Lambda}]$ [GeV]")
plt.ylim(0, 0.02)
plt.tight_layout()
plt.savefig(outdir+"lambda_mass_rec.pdf")
