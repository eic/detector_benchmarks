import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot
import pandas as pd
import os
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

plt.figure()
hep.set_style(hep.style.CMS)
hep.set_style("CMS")

def gaussian(x, amp, mean, sigma):
    return amp * np.exp( -(x - mean)**2 / (2*sigma**2) ) 

def rotateY(xdata, zdata, angle):
    s = np.sin(angle)
    c = np.cos(angle)
    rotatedz = c*zdata - s*xdata
    rotatedx = s*zdata + c*xdata
    return rotatedx, rotatedz
    
Energy = [0.005, 0.01, 0.05, 0.1, 0.5, 1.0]
q0 = [5, 10, 40, 90, 400, 700]
q1 = [0.5, 0.5, 0.9, 5, 10, 20]

df = pd.DataFrame({})
for eng in Energy:
    tree = uproot.open(f'sim_output/zdc_lyso/{os.environ["DETECTOR_CONFIG"]}_gamma_{eng}GeV_theta_0deg_thru_0.3deg.eicrecon.tree.edm4eic.root')['events']
    ecal_reco_energy = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.energy'].array()
    #hcal_reco_energy = tree['HcalFarForwardZDCClusters/HcalFarForwardZDCClusters.energy'].array()
    
    tree = uproot.open(f'sim_output/zdc_lyso/{os.environ["DETECTOR_CONFIG"]}_gamma_{eng}GeV_theta_0deg_thru_0.3deg.edm4hep.root')['events']
    ecal_sim_energy = tree['EcalFarForwardZDCHits/EcalFarForwardZDCHits.energy'].array()
    hcal_sim_energy = tree['HcalFarForwardZDCHits/HcalFarForwardZDCHits.energy'].array()

    par_x = tree['MCParticles/MCParticles.momentum.x'].array()[:,2]
    par_y = tree['MCParticles/MCParticles.momentum.y'].array()[:,2]
    par_z = tree['MCParticles/MCParticles.momentum.z'].array()[:,2]
    
    eng = int(eng*1000)

    ecal_reco_energy = pd.DataFrame({f'ecal_reco_energy_{eng}': np.array(ecal_reco_energy.tolist(), dtype=object)})
    #hcal_reco_energy = pd.DataFrame({f'hcal_reco_energy_{eng}': np.array(hcal_reco_energy.tolist(), dtype=object)})
    ecal_sim_energy = pd.DataFrame({f'ecal_sim_energy_{eng}': np.array(ecal_sim_energy.tolist(), dtype=object)})
    hcal_sim_energy = pd.DataFrame({f'hcal_sim_energy_{eng}': np.array(hcal_sim_energy.tolist(), dtype=object)})

    par_x = pd.DataFrame({f'par_x_{eng}': np.array(par_x.tolist(), dtype=object)})
    par_y = pd.DataFrame({f'par_y_{eng}': np.array(par_y.tolist(), dtype=object)})
    par_z = pd.DataFrame({f'par_z_{eng}': np.array(par_z.tolist(), dtype=object)})

    df = pd.concat([df,ecal_reco_energy,ecal_sim_energy,hcal_sim_energy,par_x,par_y,par_z],axis=1)

mu = []
sigma = []
resolution = []

fig1, ax = plt.subplots(3,2,figsize=(20,10))
#fig1.suptitle('ZDC ECal Cluster Energy Reconstruction')

plt.tight_layout()
for i in range(6):
    plt.sca(ax[i%3,i//3])
    eng = int(Energy[i]*1000)
    plt.title(f'Gamma Energy: {eng} MeV')
    temp = np.array([sum(item) if (item != 0) else 0 for item in df[f'ecal_reco_energy_{eng}']])
    hist, x = np.histogram(np.array(temp)*1000,bins=30)
    x = x[1:]/2 + x[:-1]/2
    plt.errorbar(x,hist,yerr=np.sqrt(hist),fmt='-o')
    temp = np.array([item[0] for item in df[f'ecal_reco_energy_{eng}'] if item])
    hist, x = np.histogram(np.array(temp)*1000,bins=30)
    x = x[1:]/2 + x[:-1]/2
    coeff, covar = curve_fit(gaussian,x,hist,p0=(200,q0[i],q1[i]),maxfev = 80000)
    plt.plot(np.linspace(coeff[1]-5*coeff[2],coeff[1]+5*coeff[2],50),gaussian(np.linspace(coeff[1]-5*coeff[2],coeff[1]+5*coeff[2],50),*coeff)
            ,label = f'$\mu$ = {coeff[1]:.3f} $\pm$ {covar[1][1]:.3f}\n$\sigma$ = {np.abs(coeff[2]):.3f} $\pm$ {covar[2][2]:.3f}\nResolution = {np.abs(coeff[2])*100/coeff[1]:.2f}%')
    
    plt.xlabel('Energy (MeV)')
    plt.legend()
    mu.append(coeff[1])
    sigma.append(coeff[2])
    resolution.append(np.abs(coeff[2])*100/coeff[1])
#plt.savefig('results/Energy_reconstruction_cluster.pdf')
#plt.show()

fig2, (ax1,ax2) = plt.subplots(2,1,figsize=(15,10),sharex=True)

plt.tight_layout()
# Plot data on primary axis
ax1.scatter(np.array(Energy)*1000, mu, c='b')
ax1.plot([4.5,1000],[4.5,1000],c='black',label='x=y')
ax1.set_ylabel('Reconstructed Energy (MeV)')
ax1.set_yscale('log')
ax1.legend()
ax1.set_title('ECal Craterlake Cluster Energy Reconstruction')

ax2.plot(np.array(Energy)*1000, resolution, c='r')
ax2.scatter(np.array(Energy)*1000, resolution, c='r')
ax2.set_ylabel('Resolution (%)')
ax2.set_xlabel('Gamma Energy (MeV)')
ax2.set_xscale('log')
#plt.savefig('results/Energy_resolution.pdf')
#plt.show()

htower = []
herr = []
hmean = []
hhits = []
hhits_cut = []
emean = []
ehits = []
etower = []
eerr = []
ehits_cut = []

fig3, ax = plt.subplots(2,3,figsize=(20,10))
fig3.suptitle('ZDC Simulation Energy Reconstruction')
for i in range(6):
    plt.sca(ax[i//3,i%3])
    eng = int(Energy[i]*1000)

    x = df[f'par_x_{eng}'].astype(float).to_numpy()
    y = df[f'par_y_{eng}'].astype(float).to_numpy()
    z = df[f'par_z_{eng}'].astype(float).to_numpy()
    x, z = rotateY(x,z, 0.025)
    theta = np.arccos(z/np.sqrt((x**2+y**2+z**2)))*1000
    condition = theta <= 3.5

    plt.title(f'Gamma Energy: {eng} MeV')
    energy1 = np.array([sum(item) if (item != 0) else 0 for item in df[f'hcal_sim_energy_{eng}']])#df.eval(f'hcal_sim_energy_{eng}').apply(lambda row: sum(row))
    hist, x = np.histogram(energy1*1000,bins=np.logspace(0,3,200))
    x = x[1:]/2 + x[:-1]/2
    plt.plot(x,hist,marker='o',label="HCal")
    hhits.append(len(energy1[energy1!=0]))
    condition1 = energy1!=0
    hhits_cut.append(len(energy1[condition & condition1])/len(condition[condition==True]))
    energy = np.array([sum(item) if (item != 0) else 0 for item in df[f'ecal_sim_energy_{eng}']])#df.eval(f'ecal_sim_energy_{eng}').apply(lambda row: sum(row))
    hist, x = np.histogram(energy*1000,bins=np.logspace(0,3,200))
    x = x[1:]/2 + x[:-1]/2
    plt.plot(x,hist,marker='o',label="ECal")
    emean.append(sum(energy[energy!=0])*1000/len(energy[energy!=0]))
    hmean.append(sum(energy1[energy!=0])*1000/len(energy[energy!=0]))
    condition1 = energy!=0
    ehits_cut.append(len(energy[condition & condition1])/len(condition[condition==True]))
    ehits.append(len(energy[energy!=0]))
    plt.legend()
    plt.xscale('log')
    plt.xlim(1e0,1e3)

plt.xlabel('Energy (MeV)')
#plt.savefig('results/Energy_deposition.pdf')
#plt.show()

fig4, ax = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios': [2,1]})
plt.sca(ax[0])
plt.errorbar(np.array(Energy)*1000,np.array(hmean)*47.619+np.array(emean),label='HCal+ECal',fmt='-o')
plt.errorbar(np.array(Energy)*1000,emean,label='ECal',fmt='-o')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Simulation Energy (MeV)')
plt.sca(ax[1])
plt.errorbar(np.array(Energy)*1000,(1 - np.array(emean)/(np.array(hmean)*47.619+np.array(emean)))*100,label='Total/ECal',fmt='-o')
plt.legend()
plt.ylabel('Fraction of energy\n deposited in Hcal (%)')
plt.xlabel('Truth Energy (MeV)')
#plt.savefig('results/Energy_ratio_and_Leakage.pdf')
plt.tight_layout()
#plt.show()

fig5 = plt.figure()
plt.errorbar(np.array(Energy)*1000,np.array(hhits)/1000*100,label='HCal Hits',fmt='-o')
plt.errorbar(np.array(Energy)*1000,np.array(ehits)/1000*100,label='ECal Hits',fmt='-o')
#plt.errorbar(np.array(Energy)*1000,np.array(hhits)/np.array(ehits)*100,label='HCal / ECal',fmt='-o',c='b')

plt.errorbar(np.array(Energy)*1000,np.array(hhits_cut)*100,label='HCal Hits with 3.5 mRad cut',fmt='-^')
plt.errorbar(np.array(Energy)*1000,np.array(ehits_cut)*100,label='ECal Hits with 3.5 mRad cut',fmt='-^')
#plt.errorbar(np.array(Energy)*1000,np.array(hhits_cut)/np.array(ehits_cut)*100,label='HCal / ECal with 3.5 mRad cut',fmt='-^',c='b')
### 3mrad cuts

plt.legend()
plt.xlabel('Simulation Truth Gamma Energy (MeV)')
plt.ylabel('Fraction of Events with non-zero energy (%)')
#plt.savefig('results/Hits.pdf')
plt.xscale('log')
#plt.show()

#pdfs = ['results/Energy_reconstruction_cluster.pdf','results/Energy_resolution.pdf','results/Energy_deposition.pdf','results/Energy_ratio_and_Leakage.pdf','results/Hits.pdf']
with PdfPages(f'results/{os.environ["DETECTOR_CONFIG"]}/zdc_lyso/plots.pdf') as pdf:
    pdf.savefig(fig1)
    pdf.savefig(fig2)
    pdf.savefig(fig3)
    pdf.savefig(fig4)
    pdf.savefig(fig5)

