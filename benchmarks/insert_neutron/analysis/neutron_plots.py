import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl, awkward as ak, sys, uproot as ur
import mplhep as hep
hep.style.use("CMS")

plt.rcParams['figure.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'
plt.rcParams['savefig.bbox']='tight'

plt.rcParams["figure.figsize"] = (7, 7)
config=sys.argv[1].split("/")[1]  #results/{config}/neutron
outdir=sys.argv[1]+"/"
try:
    import os
    os.mkdir(outdir[:-1])
except:
    pass

#read files
arrays_sim={}
for p in 20,30,40,50,60,70,80:
    arrays_sim[p] = ur.open(f'sim_output/insert_neutron/{config}_rec_neutron_{p}GeV.edm4eic.root:events')\
                    .arrays()

def gauss(x, A,mu, sigma):
    return A * np.exp(-(x-mu)**2/(2*sigma**2))

#get the truth pseudorapidity and energy
for array in arrays_sim.values():
    tilt=-0.025
    px=array['MCParticles.momentum.x'][:,2]
    py=array['MCParticles.momentum.y'][:,2]
    pz=array['MCParticles.momentum.z'][:,2]
    p=np.sqrt(px**2+py**2+pz**2)
    
    pxp=px*np.cos(tilt)-pz*np.sin(tilt)
    pyp=py
    pzp=pz*np.cos(tilt)+px*np.sin(tilt)
    
    array['E_truth']=np.hypot(p, 0.9406)
    array['eta_truth']=1/2*np.log((p+pzp)/(p-pzp))
    array['theta_truth']=np.arccos(pzp/p)

#
# get reconstructed theta as avg of theta of cluster centers, weighted by energy
for array in arrays_sim.values():
    tilt=-0.025
    x=array['HcalEndcapPInsertClusters.position.x']
    y=array['HcalEndcapPInsertClusters.position.y']
    z=array['HcalEndcapPInsertClusters.position.z']
    E=array['HcalEndcapPInsertClusters.energy']
    r=np.sqrt(x**2+y**2+z**2)
    
    xp=x*np.cos(tilt)-z*np.sin(tilt)
    yp=y
    zp=z*np.cos(tilt)+x*np.sin(tilt)
    
    w=E
    
    array['theta_recon']=np.sum(np.arccos(zp/r)*w, axis=-1)/np.sum(w, axis=-1)
    array['eta_recon']=-np.log(np.tan(array['theta_recon']/2))
    

#plot theta residuals:
print("making theta recon plot")
from scipy.optimize import curve_fit

fig, axs=plt.subplots(1,2, figsize=(16,8))
plt.sca(axs[0])
p=40
eta_min=3.4; eta_max=3.6
y,x,_=plt.hist(1000*(arrays_sim[p]['theta_recon']-arrays_sim[p]['theta_truth'])\
               [(arrays_sim[p]['eta_truth']>eta_min)&(arrays_sim[p]['eta_truth']<eta_max)], bins=50,
                    range=(-10,10), histtype='step')
bc=(x[1:]+x[:-1])/2
slc=abs(bc)<3
# try:
fnc=gauss
sigma=np.sqrt(y[slc])+(y[slc]==0)
p0=(100, 0, 5)
coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
xx=np.linspace(-5,5,100)
plt.plot(xx,fnc(xx,*coeff))
# except:
#     pass
plt.xlabel("$\\theta_{rec}-\\theta_{truth}$ [mrad]")
plt.ylabel("events")
plt.title(f"$p={p}$ GeV, ${eta_min}<\\eta<{eta_max}$")

r=[3.2,3.4,3.6,3.8,4.0]
for eta_min, eta_max in zip(r[:-1],r[1:]):
    xvals=[]
    sigmas=[]
    dsigmas=[]
    for p in 20,30,40, 50, 60, 70, 80:
        y,x=np.histogram(1000*(arrays_sim[p]['theta_recon']-arrays_sim[p]['theta_truth'])\
                         [(arrays_sim[p]['eta_truth']>eta_min)&(arrays_sim[p]['eta_truth']<eta_max)],
                         bins=50, range=(-10,10))
        bc=(x[1:]+x[:-1])/2
        slc=abs(bc)<3
        fnc=gauss
        p0=(100, 0, 5)
        #print(bc[slc],y[slc])
        sigma=np.sqrt(y[slc])+(y[slc]==0)
        try:
            coeff, var_matrix = curve_fit(fnc, list(bc[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
            sigmas.append(np.abs(coeff[2]))
            dsigmas.append(np.sqrt(var_matrix[2][2]))
            xvals.append(p)
        except:
            pass
    plt.sca(axs[1])
    plt.errorbar(xvals, sigmas, dsigmas, ls='', marker='o', label=f"${eta_min}<\\eta<{eta_max}$")
plt.xlabel("$p_{n}$ [GeV]")
plt.ylabel("$\\sigma[\\theta]$ [mrad]")
plt.ylim(0)
plt.legend()
plt.tight_layout()
plt.savefig(outdir+"neutron_theta_recon.pdf")

#now determine the energy recon parameters
pvals=[]
resvals=[]
reserrs=[]
wvals=[]
svals=[]
Euncorr=[]

print("determining the energy recon correction parameters")
fig,axs=plt.subplots(1,2, figsize=(20,10))
eta_min=3.4;eta_max=3.6
for p in 20, 30,40,50,60,70, 80:
    best_res=1000
    res_err=1000
    best_s=1000
    wrange=np.linspace(30, 70, 41)*0.0257
    coeff_best=None
    
    wbest=0
    a=arrays_sim[p]
    h=np.sum(a[f'HcalEndcapPInsertClusters.energy'], axis=-1)
    e=np.sum(a[f'EcalEndcapPInsertClusters.energy'], axis=-1)
    for w in wrange:
        
        r=(e/w+h)[(h>0)&(a['eta_truth']>eta_min)&(a['eta_truth']<eta_max)]
        y,x=np.histogram(r,bins=50)
        bcs=(x[1:]+x[:-1])/2
        fnc=gauss
        slc=abs(bcs-np.mean(r)*1.25)<2*np.std(r)
        sigma=np.sqrt(y[slc])+0.5*(y[slc]==0)
        p0=(100, np.mean(r), np.std(r))
        try:
            coeff, var_matrix = curve_fit(fnc, list(bcs[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
            res=np.abs(coeff[2]/coeff[1])
            
            if res<best_res:
                best_res=res
                coeff_best=coeff
                res_err=res*np.hypot(np.sqrt(var_matrix[2][2])/coeff[2], np.sqrt(var_matrix[1][1])/coeff[1])
                wbest=w
                best_s=np.hypot(p,0.9406)/coeff[1]
                Euncorr_best=coeff[1]
        except :
            print("fit failed")
    
    if p==50:
        r=(e/wbest+h)[(h>0)&(a['eta_truth']>3.4)&(a['eta_truth']<3.6)]
        plt.sca(axs[0])
        y, x, _= plt.hist(r, histtype='step', bins=50)
        xx=np.linspace(20, 55, 100)
        plt.plot(xx,fnc(xx, *coeff_best), ls='-')
        plt.xlabel("$E_{uncorr}=E_{Hcal}+E_{Ecal}/w$ [GeV]")
        plt.title(f"p=50 GeV, ${eta_min}<\\eta<{eta_max}$, w={wbest:.2f}")
        plt.axvline(np.sqrt(50**2+.9406**2), color='g', ls=':')
        plt.text(40, max(y)*0.9, "generated\nenergy", color='g', fontsize=20)
        
    Euncorr.append(Euncorr_best)
    resvals.append(best_res)
    reserrs.append(res_err)
    pvals.append(p)
    svals.append(best_s)
    wvals.append(wbest)

pvals=np.array(pvals)
svals=np.array(svals)
Euncorr=np.array(Euncorr)
plt.sca(axs[1])
plt.plot(Euncorr,wvals, label="w")
w_avg=np.mean(wvals)
plt.axhline(w_avg, label=f'w avg: {w_avg:.2f}', ls=':')
plt.plot(Euncorr,svals, label="s")
m=(np.sum(svals*Euncorr)*len(Euncorr)-np.sum(Euncorr)*np.sum(svals))/(np.sum(Euncorr**2)*len(Euncorr)-np.sum(Euncorr)**2)
b=np.mean(svals)-np.mean(Euncorr)*m
plt.plot(Euncorr,Euncorr*m+b, label=f"s fit: ${m:.4f}E_{{uncorr}}+{b:.2f}$", ls=':')

plt.xlabel("$E_{uncorr}=E_{Hcal}+E_{Ecal}/w$ [GeV]")
plt.title("$E_{n,recon}=s\\times(E_{Hcal}+E_{Ecal}/w)$")
plt.ylabel('parameter values')
plt.legend()
plt.ylim(0)
plt.savefig(outdir+"neutron_energy_params.pdf")

#now make the energy plot
print("making energy recon plot")
fig, axs=plt.subplots(1,3, figsize=(24,8))
partitions=[3.2,3.4, 3.6, 3.8, 4.0]
for eta_min, eta_max in zip(partitions[:-1],partitions[1:]):
    pvals=[]
    resvals=[]
    reserrs=[]
    scalevals=[]
    scaleerrs=[]
    for p in 20, 30,40,50,60,70, 80:
        best_res=1000
        res_err=1000


        wrange=np.linspace(30, 70, 30)*0.0257
        
        w=w_avg
        a=arrays_sim[p]
        h=np.sum(a[f'HcalEndcapPInsertClusters.energy'], axis=-1)
        e=np.sum(a[f'EcalEndcapPInsertClusters.energy'], axis=-1)
        #phi=a['phi_truth']
        uncorr=(e/w+h)
        s=-0.0064*uncorr+1.80
        r=uncorr*s #reconstructed energy with correction
        r=r[(h>0)&(a['eta_truth']>eta_min)&(a['eta_truth']<eta_max)]#&(abs(phi)>np.pi/2)]
        y,x=np.histogram(r,bins=50)
        bcs=(x[1:]+x[:-1])/2
        fnc=gauss
        slc=abs(bcs-np.mean(r)*1.25)<2*np.std(r)
        sigma=np.sqrt(y[slc])+0.5*(y[slc]==0)
        p0=(100, np.mean(r), np.std(r))
        try:
            coeff, var_matrix = curve_fit(fnc, list(bcs[slc]), list(y[slc]), p0=p0,sigma=list(sigma))
            res=np.abs(coeff[2]/coeff[1])

            if res<best_res:
                best_res=res
                res_err=res*np.hypot(np.sqrt(var_matrix[2][2])/coeff[2], np.sqrt(var_matrix[1][1])/coeff[1])
                wbest=w
                scale=coeff[1]/np.sqrt(p**2+0.9406**2)
                dscale=np.sqrt(var_matrix[1][1]/np.sqrt(p**2+0.9406**2))
        except :
            print("fit failed")
        if p==50 and eta_min==3.4:
            plt.sca(axs[0])
            plt.errorbar(bcs, y, np.sqrt(y)+(y==0),marker='o', ls='', )
            plt.title(f'p={p} GeV, ${eta_min}<\\eta<{eta_max}$')
            plt.xlabel("$E_{recon}$ [GeV]")
            plt.ylabel("events")
            #plt.ylim(0)
            xx=np.linspace(40, 70, 50)
            plt.plot(xx, fnc(xx, *coeff))
        resvals.append(best_res)
        reserrs.append(res_err)
        scalevals.append(scale)
        scaleerrs.append(dscale)
        pvals.append(p)
    plt.sca(axs[1])
    plt.errorbar(pvals, resvals, reserrs, marker='o', ls='', label=f"${eta_min}<\\eta<{eta_max}$")
    #plt.ylim(0)
    plt.ylabel("$\\sigma[E]/\\mu[E]$")
    plt.xlabel("$p_{n}$ [GeV]")

    plt.sca(axs[2])
    plt.errorbar(pvals, scalevals, scaleerrs, marker='o', ls='', label=f"${eta_min}<\\eta<{eta_max}$")
    
    
    plt.ylabel("$\\mu[E]/E$")


axs[2].set_xlabel("$p_n$ [GeV]")
axs[2].axhline(1, ls='--', color='0.5', alpha=0.7)
axs[0].set_ylim(0)
axs[1].set_ylim(0, 0.35)
axs[2].set_ylim(0)
axs[1].legend()
axs[2].legend()
plt.tight_layout()
plt.savefig(outdir+"neutron_energy_recon.pdf")
