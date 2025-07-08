#!/usr/bin/env python
# Copyright 2023, Christopher Dilks
# Subject to the terms in the LICENSE file found in the top-level directory.

import ROOT as r
import sys, getopt, pathlib, math

# suppress graphics
r.gROOT.SetBatch(True)

# GLOBAL SETTINGS
################################################################
RESID_MAX = 10 # cherenkov angle |residual| maximum
RADIATORS = {
        'Aerogel': { 'p_max': 30, 'p_rebin': 12, 'rindex_ref': 1.0190},
        'Gas':     { 'p_max': 70, 'p_rebin': 30, 'rindex_ref': 1.00076},
        'Merged':  { 'p_max': 70, 'p_rebin': 30, },
        }
PARTICLE_MASS = {
        'e':  0.00051,
        'pi': 0.13957,
        'K':  0.49368,
        'p':  0.93827,
        }

# ARGUMENTS
################################################################

ana_file_name = 'out/ana.edm4hep.root'
output_dir    = 'out/ana.plots'

helpStr = f'''
{sys.argv[0]} [OPTIONS]

    -i <input file>: specify an input file, e.g., hepmc
       default: {ana_file_name}

    -o <output dir>: specify an output directory
       default: {output_dir}
    
    -h: show this usage guide

    '''

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:o:h')
except getopt.GetoptError:
    print('\n\nERROR: invalid argument\n', helpStr, file=sys.stderr)
    sys.exit(2)
for opt, arg in opts:
    if(opt == '-i'): ana_file_name = arg.lstrip()
    if(opt == '-o'): output_dir    = arg.lstrip()
    if(opt == '-h'):
        print(helpStr)
        sys.exit(2)
print(f'''
ana_file_name = {ana_file_name}
output_dir    = {output_dir}
''')


# PLOTTING
################################################################

# make canvases
ana_file = r.TFile.Open(ana_file_name, "READ")
def make_canv(name, nx, ny):
    canv = r.TCanvas(name, name, nx*600, ny*400)
    canv.Divide(nx,ny)
    return canv
canv_dict = {
    "photon_spectra": make_canv("photon_spectra", 2, 2),
    "digitization":   make_canv("digitization",   2, 2),
    "pidAerogel":     make_canv("pidAerogel",     5, 3),
    "pidGas":         make_canv("pidGas",         5, 3),
    "pidMerged":      make_canv("pidMerged",      2, 3),
}

# helper functions
### minimum momentum for Cherenkov radiation
def calculate_mom_min(mass,n):
    return mass / math.sqrt(n**2-1)
### Cherenkov angle at saturation
def calculate_theta_max(n):
    return 1000 * math.acos(1/n)
### execute `op` for every pad of `canv`
def loop_pads(canv, op):
    for pad in canv.GetListOfPrimitives():
        if pad.InheritsFrom(r.TVirtualPad.Class()):
            op(pad)
### draw a profile on 2D histogram `hist`
def draw_profile(hist, rad, style="BOX"):
    prof = hist.ProfileX("_pfx", 1, -1, "i")
    prof.SetMarkerStyle(r.kFullCircle)
    prof.SetMarkerColor(r.kRed)
    prof.SetMarkerSize(2)
    hist_rebinned = hist.Clone(f'{hist.GetName()}_rebin')
    if 'vs_p' in hist.GetName():
        hist_rebinned.RebinX(RADIATORS[rad]['p_rebin'])
    else:
        hist_rebinned.RebinX(4)
    hist_rebinned.SetLineColor(r.kBlue)
    hist_rebinned.SetFillColor(r.kBlue)
    for p in [hist_rebinned, prof]:
        if 'vs_p' in hist.GetName():
            p.GetXaxis().SetRangeUser(0, RADIATORS[rad]['p_max'])
        if 'rindex_ref' in RADIATORS[rad] and 'theta_vs_p' in hist.GetName():
            p.GetYaxis().SetRangeUser(0, 1.5*calculate_theta_max(RADIATORS[rad]['rindex_ref']))
        if 'thetaResid_vs' in hist.GetName():
            p.GetYaxis().SetRangeUser(-RESID_MAX, RESID_MAX)
    hist_rebinned.Draw(style)
    prof.Draw("SAME")

# Cherenkov angle vs. p functions
for rad, radHash in RADIATORS.items():
    if rad == "Merged":
        continue
    radHash['cherenkov_curves'] = []
    n = radHash['rindex_ref']
    for part, mass in PARTICLE_MASS.items():
        ftn = r.TF1(
                f'ftn_theta_{part}_{rad}',
                f'1000*TMath::ACos(TMath::Sqrt(x^2+{mass}^2)/({n}*x))',
                calculate_mom_min(mass, n),
                radHash['p_max']
                )
        ftn.SetLineColor(r.kBlack)
        ftn.SetLineWidth(1)
        radHash['cherenkov_curves'].append(ftn)

# draw photon spectra
canv = canv_dict["photon_spectra"]
def canv_op(pad):
    pad.SetGrid(1,1)
    pad.SetLogy()
loop_pads(canv, canv_op)
canv.cd(1)
ana_file.Get("phot/phot_spectrum_sim").Draw()
canv.cd(2)
ana_file.Get("phot/nphot_dist").Draw()
canv.cd(3)
ana_file.Get("digi/phot_spectrum_rec").Draw()
canv.cd(4)
ana_file.Get("digi/nhits_dist").Draw()

# draw digitization
canv = canv_dict["digitization"]
def canv_op(pad):
    pad.SetGrid(1,1)
    if pad.GetNumber() < 3:
        pad.SetLogy()
    else:
        pad.SetLogz()
loop_pads(canv, canv_op)
canv.cd(1)
ana_file.Get("digi/adc_dist").Draw()
canv.cd(2)
ana_file.Get("digi/tdc_dist").Draw()
canv.cd(3)
ana_file.Get("digi/tdc_vs_adc").Draw("COLZ")

# draw CherenkovPID for each radiator
for rad in ["Aerogel", "Gas"]:
    pid_name = f'pid{rad}'
    canv = canv_dict[pid_name]
    def canv_op(pad):
        pad.SetGrid(1,1)
    loop_pads(canv, canv_op)
    canv.cd(1)
    ana_file.Get(f'{pid_name}/npe_dist_{rad}').Draw()
    canv.cd(2)
    ana_file.Get(f'{pid_name}/theta_dist_{rad}').Draw()
    canv.cd(3)
    hist = ana_file.Get(f'{pid_name}/thetaResid_dist_{rad}')
    hist.DrawCopy()
    canv.cd(4)
    hist.SetTitle(hist.GetTitle() + " - ZOOM")
    hist.GetXaxis().SetRangeUser(-RESID_MAX, RESID_MAX)
    hist.Draw()
    resid_mode = hist.GetBinCenter(hist.GetMaximumBin())
    resid_dev  = hist.GetStdDev()
    nsigma = 3
    hist.Fit("gaus", "", "", resid_mode - nsigma*resid_dev, resid_mode + nsigma*resid_dev)
    if hist.GetFunction("gaus"):
        hist.GetFunction("gaus").SetLineWidth(5)
    canv.cd(5)
    ana_file.Get(f'{pid_name}/highestWeight_dist_{rad}').Draw()
    canv.cd(6)
    hist = ana_file.Get(f'{pid_name}/npe_vs_p_{rad}')
    draw_profile(hist, rad)
    canv.cd(7)
    hist = ana_file.Get(f'{pid_name}/theta_vs_p_{rad}')
    draw_profile(hist, rad)
    if 'cherenkov_curves' in RADIATORS[rad]:
        for ftn in RADIATORS[rad]['cherenkov_curves']:
            ftn.Draw("SAME")
    canv.cd(8)
    hist = ana_file.Get(f'{pid_name}/thetaResid_vs_p_{rad}')
    hist.GetYaxis().SetRangeUser(-RESID_MAX, RESID_MAX)
    draw_profile(hist, rad)
    canv.cd(9)
    ana_file.Get(f'{pid_name}/photonTheta_vs_photonPhi_{rad}').Draw("COLZ")
    canv.cd(10)
    ana_file.Get(f'{pid_name}/mcRindex_{rad}').Draw()
    canv.cd(11)
    hist = ana_file.Get(f'{pid_name}/npe_vs_eta_{rad}')
    draw_profile(hist, rad)
    canv.cd(12)
    hist = ana_file.Get(f'{pid_name}/theta_vs_eta_{rad}')
    draw_profile(hist, rad)
    canv.cd(13)
    hist = ana_file.Get(f'{pid_name}/thetaResid_vs_eta_{rad}')
    hist.GetYaxis().SetRangeUser(-RESID_MAX, RESID_MAX)
    draw_profile(hist, rad)

# draw CherenkovPID for merged radiators
pid_name = f'pidMerged'
canv = canv_dict[pid_name]
def canv_op(pad):
    pad.SetGrid(1,1)
loop_pads(canv, canv_op)
canv.cd(1)
ana_file.Get(f'{pid_name}/npe_dist_Merged').Draw()
canv.cd(2)
ana_file.Get(f'{pid_name}/highestWeight_dist_Merged').Draw()
canv.cd(3)
hist = ana_file.Get(f'{pid_name}/npe_vs_p_Merged')
draw_profile(hist, 'Merged')
canv.cd(4)
ana_file.Get(f'{pid_name}/photonTheta_vs_photonPhi_Merged').Draw("COLZ")
canv.cd(5)
hist = ana_file.Get(f'{pid_name}/npe_vs_eta_Merged')
draw_profile(hist, 'Merged')

# FINISH
################################################################

pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
for name, canvas in canv_dict.items():
    canvas.SaveAs(f'{output_dir}/{name}.png')
ana_file.Close()
