import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import uproot
import pandas as pd
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

plt.figure()
hep.set_style(hep.style.CMS)
hep.set_style("CMS")

def gaussian(x, amp, mean, sigma):
    return amp * np.exp( -(x - mean)**2 / (2*sigma**2) ) 

Energy = [0.005, 0.01, 0.05, 0.1, 0.5, 1]

df = pd.DataFrame({})
for eng in Energy:
    tree = uproot.open(f'data/epic_craterlake_reco_gamma_{eng}GeV_theta_0deg_thru_0.3deg.edm4hep.root')['events']
    ecal_reco_energy = tree['EcalFarForwardZDCClusters/EcalFarForwardZDCClusters.energy'].array()
    #hcal_reco_energy = tree['HcalFarForwardZDCClusters/HcalFarForwardZDCClusters.energy'].array()

    tree = uproot.open(f'data/epic_craterlake_sim_gamma_{eng}GeV_theta_0deg_thru_0.3deg.edm4hep.root')['events']
    ecal_sim_energy = tree['EcalFarForwardZDCHits/EcalFarForwardZDCHits.energy'].array()
    hcal_sim_energy = tree['HcalFarForwardZDCHits/HcalFarForwardZDCHits.energy'].array()

    eng = int(eng*1000)

    ecal_reco_energy = pd.DataFrame({f'ecal_reco_energy_{eng}': np.array(ecal_reco_energy.tolist(), dtype=object)})
    #hcal_reco_energy = pd.DataFrame({f'hcal_reco_energy_{eng}': np.array(hcal_reco_energy.tolist(), dtype=object)})
    ecal_sim_energy = pd.DataFrame({f'ecal_sim_energy_{eng}': np.array(ecal_sim_energy.tolist(), dtype=object)})
    hcal_sim_energy = pd.DataFrame({f'hcal_sim_energy_{eng}': np.array(hcal_sim_energy.tolist(), dtype=object)})

    df = pd.concat([df,ecal_reco_energy,ecal_sim_energy,hcal_sim_energy],axis=1)


mu = []
sigma = []
resolution = []
fig1, ax = plt.subplots(3,2,figsize=(20,10))
plt.tight_layout()
for i in range(6):
    plt.sca(ax[i%3,i//3])
    eng = int(Energy[i]*1000)
    plt.title(f'Gamma Energy: {eng} MeV')
    temp = [item[0] for item in df[f'ecal_reco_energy_{eng}'] if item]
    hist, x = np.histogram(temp,bins=15)
    x = x[1:]/2 + x[:-1]/2
    plt.scatter(x,hist,marker='o')
    coeff, covar = curve_fit(gaussian,x,hist,p0=(200,eng,eng/2))
    plt.plot(np.linspace(coeff[1]-5*coeff[2],coeff[1]+5*coeff[2],50),gaussian(np.linspace(coeff[1]-5*coeff[2],coeff[1]+5*coeff[2],50),*coeff)
            ,label = f'$\mu$ = {coeff[1]:.3f} $\pm$ {covar[1][1]:.3f}\n$\sigma$ = {np.abs(coeff[2]):.3f} $\pm$ {covar[2][2]:.3f}\nResolution = {np.abs(coeff[2])*100/coeff[1]:.2f}%')
    plt.xlabel('Energy (MeV)')
    plt.legend()
    mu.append(coeff[1])
    sigma.append(coeff[2])
    resolution.append(np.abs(coeff[2])*100/coeff[1])
plt.savefig('results/Energy_reconstruction_cluster.pdf')

plt.show()

fig2, (ax1,ax2) = plt.subplots(2,1,figsize=(15,10),sharex=True)
plt.tight_layout()
# Plot data on primary axis
ax1.scatter(Energy, mu, c='b')
ax1.plot([0.0045,1],[0.0045,1],c='black',label='x=y')
ax1.set_ylabel('Reconstructed Energy (GeV)')
ax1.legend()
ax1.set_title('ECal Craterlake Energy Reconstruction')

ax2.plot(Energy, resolution, c='r')
ax2.scatter(Energy, resolution, c='r')
ax2.set_ylabel('Resolution (%)')
ax2.set_xlabel('Gamma Energy (GeV)')
ax2.set_xscale('log')
plt.savefig('results/Energy_resolution.pdf')

plt.show()


htower = []
herr = []
hmean = []
hhits = []
emean = []
ehits = []
etower = []
eerr = []
fig3, ax = plt.subplots(2,3,figsize=(20,10))
for i in range(6):
    plt.sca(ax[i//3,i%3])
    
    eng = int(Energy[i]*1000)

    energy = np.array([sum(item) for item in df[f'hcal_sim_energy_{eng}'] if item])#df.eval(f'hcal_sim_energy_{eng}').apply(lambda row: sum(row))
    hist, x = np.histogram(energy*1000,bins=np.logspace(0,3,200))
    x = x[1:]/2 + x[:-1]/2
    plt.plot(x,hist,marker='o',label="HCal")
    hmean.append(sum(energy)*1000/len(energy))
    hhits.append(len(energy[energy!=0]))
    energy = np.array([sum(item) for item in df[f'ecal_sim_energy_{eng}'] if item])#df.eval(f'ecal_sim_energy_{eng}').apply(lambda row: sum(row))
    hist, x = np.histogram(energy*1000,bins=np.logspace(0,3,200))
    x = x[1:]/2 + x[:-1]/2
    plt.plot(x,hist,marker='o',label="ECal")
    emean.append(sum(energy)*1000/len(energy))
    ehits.append(len(energy[energy!=0]))
    plt.legend()
    plt.xscale('log')
    plt.xlim(1e0,1e3)
plt.savefig('results/Energy_deposition.pdf')
plt.show()

fig4, ax = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios': [2,1]})
plt.sca(ax[0])
plt.errorbar(np.array(Energy)*1000,hmean,label='HCal Active Layers',fmt='-o')
plt.errorbar(np.array(Energy)*1000,np.array(hmean)*47.619,label='HCal',fmt='-o')
plt.errorbar(np.array(Energy)*1000,np.array(hmean)*47.619+emean,label='HCal+ECal',fmt='-o')
plt.errorbar(np.array(Energy)*1000,emean,label='ECal',fmt='-o')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Simulation Energy (MeV)')
plt.sca(ax[1])
plt.errorbar(np.array(Energy)*1000,np.array(hmean)/np.array(emean),label='HCal Active Layers/ECal',fmt='-o')
plt.errorbar(np.array(Energy)*1000,np.array(hmean)/np.array(emean)*47.619,label='HCal/ECal',fmt='-o')
plt.legend()
plt.ylabel('Leakage Ratio')
plt.xlabel('Truth Energy (MeV)')
plt.savefig('results/Energy_ratio_and_Leakage.pdf')
plt.tight_layout()
plt.show()

fig5 = plt.figure()
plt.scatter(Energy,np.array(hhits)/1000*100,label='HCal Hits')
plt.scatter(Energy,np.array(ehits)/1000*100,label='ECal Hits')
plt.plot(Energy,np.array(hhits)/1000*100)
plt.plot(Energy,np.array(ehits)/1000*100)
plt.plot(Energy,np.array(hhits)/np.array(ehits)*100,label='HCal / ECal')
plt.scatter(Energy,np.array(hhits)/np.array(ehits)*100)

plt.legend()
plt.xlabel('Simulation Truth Gamma Energy (GeV)')
plt.ylabel('Simulation Hits at ZDC (%)')
plt.savefig('results/Hits.pdf')
plt.show()

pdfs = ['results/Energy_reconstruction_cluster.pdf','results/Energy_resolution.pdf','results/Energy_deposition.pdf','results/Energy_ratio_and_Leakage.pdf','results/Hits.pdf']
with PdfPages("results/plots.pdf") as pdf:
    pdf.savefig(fig1)
    pdf.savefig(fig2)
    pdf.savefig(fig3)
    pdf.savefig(fig4)
    pdf.savefig(fig5)
