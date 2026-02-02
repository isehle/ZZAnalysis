import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT
import uproot as up
import numpy as np

import json

from NanoAnalysis.scripts.ZpX_estimation import ZpX

procs       = ["ggZZ", "ZpX", "VVV", "H", "ZLZL", "ZLZT", "ZTZT", "ZZ_LO", "ZZ_NLO"]
pol_samples = ['ZLZL', 'ZLZT', 'ZTZT', 'ZZ_LO']
fstates     = ['fs_4l', 'fs_2x2e', 'fs_2x2mu']
#fstates     = ['fs_2x2e', 'fs_2x2mu']
zpx_procs   = ["DY", "TT", "WZ"]

properties = [
    "cosTheta1",
    "cosTheta3",
    "cosThetaStar",
    #"delPhiStar",
    #"delPhi",
    "delRapidity"
]

variations = ["", "_LHEScaleWeightUp", "_LHEScaleWeightDown", "_LHEPdfWeightUp", "_LHEPdfWeightDown", "_puWeightUp", "_puWeightDown", "_lepIDRecUp", "_lepIDRecDown"]

def get_vars(proc):
    if proc in pol_samples + ["H"]:
        return ["", "_LHEScaleWeightUp", "_LHEScaleWeightDown", "_LHEPdfWeightUp", "_LHEPdfWeightDown", "_puWeightUp", "_puWeightDown", "_lepIDRecUp", "_lepIDRecDown"]
    elif proc in ["ggZZ", "VVV"]:
        return ["", "_puWeightUp", "_puWeightDown", "_lepIDRecUp", "_lepIDRecDown"]
    else:
        return [""]

norm_hist = lambda hist, val=1.: hist.Scale(val/hist.Integral())

get_key = lambda key: key.replace(";1","")

def fill_empty_hist(hist, content=1e-5, error=1e-3):
    for bin_idx in range(hist.GetNbinsX()):
        hist.SetBinContent(bin_idx+1, content)
        hist.SetBinError(bin_idx+1, error)
    return hist

def sum_hist_list(hists, subprocs, fs_key, var, proc):
    hist_list = [hists[subproc][fs_key][var] for subproc in subprocs]
    new_hist = hist_list[0].Clone(proc)
    for hist in hist_list[1:]:
        new_hist.Add(hist)
    
    # Protect against absolute zero yields
    if new_hist.Integral() == 0:
        new_hist = fill_empty_hist(new_hist)

    return new_hist

# def zpx_shape(zpx_hist, year):

#     # best fit values from (Data - notZpX) in HighMassSSRelaxed CR (4l, Full2022 or Full2023)
#     if year == 2022:
#         params = (349.358, 63.3026, -48.312, 10.9672)
#     elif year == 2023:
#         params = (309.166, -70.8494, -29.6992, 7.13474)

#     def fit_func(x, *popt):
#         return popt[0] + popt[1]*x + popt[2]*(x**2) + popt[3]*(x**3)
#         #return 349.358 - 63.3026*x - 48.312*(x**2) + 10.9672*(x**3)

#     for bin_idx in range(zpx_hist.GetNbinsX()):
#         bin_center = zpx_hist.GetBinCenter(bin_idx+1)
#         bin_val    = fit_func(bin_center, *params)

#         zpx_hist.SetBinContent(bin_idx+1, bin_val)
#         zpx_hist.SetBinError(bin_idx+1, 0.)
    
#     zpx_hist.Scale(1./zpx_hist.Integral())
    
#     return zpx_hist

def get_zpx_hist(full_hist):
    hists = {}
    with up.open(full_hist) as HistFile:
        for proc in HistFile["delRapidity"]["HighMassSSRelaxed"]["fs_4l"].keys():
            key = proc.replace(";1", "")
            hist = HistFile["delRapidity"]["HighMassSSRelaxed"]["fs_4l"][proc].to_pyroot()
            #hists[key] = hist.Rebin(8, "dely_rebin", np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]))
            hists[key] = hist.Rebin(2, "dely_rebin")
    
    new_hist = hists["Data"].Clone("ZpX")
    for proc in hists.keys():
        if proc not in zpx_procs + ["Data", "ZLZL", "ZLZT", "ZTZT", "ZZ_LO"]:
            new_hist.Add(hists[proc], -1)

    return new_hist

def get_fit_pars(new_hist, fit_func):
    #f = ROOT.TF1("f", fit_func, 0.0, 4.0, npar=4, ndim=1)
    f = ROOT.TF1("f", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0.0, 4.0)

    fit_result = new_hist.Fit("f", "LES")

    pars = []
    errs = []
    for i in range(4):
        pars.append(fit_result.Parameter(i))
        errs.append(fit_result.ParError(i))

    nom_pars = np.array(pars)
    up_pars  = nom_pars + np.array(errs)
    dn_pars  = nom_pars - np.array(errs)

    return nom_pars, up_pars, dn_pars

def zpx_shape(full_hist, zpx_norm):
    def fit_func(x, *popt):
        return popt[0] + popt[1]*x + popt[2]*(x**2) + popt[3]*(x**3)

    new_hist = get_zpx_hist(full_hist)

    nom_pars, up_pars, dn_pars = get_fit_pars(new_hist, fit_func)

    # nom_hist = ROOT.TH1D("delRapidity", "delRapidity", 8, 0., 4.)
    # up_hist  = ROOT.TH1D("delRapidity", "delRapidityUp", 8, 0., 4.)
    # dn_hist  = ROOT.TH1D("delRapidity", "delRapidityDown", 8, 0., 4.)

    nom_hist = ROOT.TH1D("delRapidity", "delRapidity", 20, 0., 4.)
    up_hist  = ROOT.TH1D("delRapidity", "delRapidityUp", 20, 0., 4.)
    dn_hist  = ROOT.TH1D("delRapidity", "delRapidityDown", 20, 0., 4.)

    for bin_idx in range(20):
        bin_center = nom_hist.GetBinCenter(bin_idx+1)

        nom_val = fit_func(bin_center, *nom_pars)
        up_val  = fit_func(bin_center, *up_pars)
        dn_val  = fit_func(bin_center, *dn_pars)

        nom_hist.SetBinContent(bin_idx+1, max(nom_val, 0))
        up_hist.SetBinContent(bin_idx+1, max(up_val, 0))
        dn_hist.SetBinContent(bin_idx+1, max(dn_val, 0))

        nom_hist.SetBinError(bin_idx+1, 0.)
        up_hist.SetBinError(bin_idx+1, 0.)
        dn_hist.SetBinError(bin_idx+1, 0.)

    nom_count = nom_hist.Integral()
    # up_count  = up_hist.Integral()
    # dn_count  = dn_hist.Integral()

    nom_hist.Scale(zpx_norm/nom_count)
    # nom_hist.Scale(1./nom_count)

    up_hist.Scale(zpx_norm/up_hist.Integral())
    dn_hist.Scale(zpx_norm/dn_hist.Integral())

    # up_hist.Scale(1./nom_count)
    # dn_hist.Scale(1./nom_count)

    # nom_hist.Scale(1./nom_hist.Integral())
    # up_hist.Scale(1./up_hist.Integral())
    # dn_hist.Scale(1./dn_hist.Integral())
    
    return nom_hist, up_hist, dn_hist
    

if __name__ == "__main__":
    from tqdm import tqdm

    # For Polarization variation normalizations
    # zz_nlo_var_norms = "/eos/user/i/iehle/ZZ_NLO_normVars.json"
    zz_nlo_var_norms = "/eos/user/i/iehle/ZZ_NLO_normVars_reweightNonDegEvents22_23_09_25.json"

    with open(zz_nlo_var_norms, "r") as zz_nlo_norms:
        the_dict = json.load(zz_nlo_norms)
        zz_2022, zz_2023 = the_dict["2022"], the_dict["2023"]

    # Used to get ZpX shape
    #full_hist = "/eos/user/i/iehle/Analysis/histograms/Full/hists_dataInCRs.root"
    full_hist = "/eos/user/i/iehle/Analysis/histograms/Full/hists_delRapidity_HMSSRelaxed_reweightNonDegEvents22_23_09_25.root"

    # For ZpX normalization
    #zpx_2022  = "/afs/cern.ch/user/i/iehle/cmssw/CMSSW_14_1_6/src/ZZAnalysis/ZpX_info_2022_Full_2x2e_2x2mu_4l.json"
    zpx_2022  = "/eos/user/i/iehle/ZpX_info_2022_Full_reweightNonDegEvents_22_09_25.json"
    zpx_2023  = "/afs/cern.ch/user/i/iehle/cmssw/CMSSW_14_1_6/src/ZZAnalysis/ZpX_info_2023_Full_2x2e_2x2mu_4l.json"

    with open(zpx_2022, "r") as ZpX22, open(zpx_2023, "r") as ZpX23:
        dict_2022, dict_2023 = json.load(ZpX22), json.load(ZpX23)
        zpx_info_2022 = dict_2022["N_ZpX_MidMass"]
        zpx_info_2023 = dict_2023["N_ZpX_MidMass"]

    # hist_file = "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_24_06_25.root"
    # outfile   = "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_zpxForPlots_11_07_25.root"

    #hist_file = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_sigOnTop_allFStates.root"
    #outfile   = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_delRapidity_GoodNorms_03_07_25.root"

    #hist_file = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_sigOnTop_allFStates.root"

    # hist_file = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_reweightNonDegEvents_22_09_25.root"
    # outfile   = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_goodZpX_reweightNonDegEvents_22_09_25.root"

    hist_file = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_reweightNonDegEvents_22_09_25.root"
    outfile   = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_ZpXForPlots_reweightNonDegEvents_22_09_25.root"

    with up.open(hist_file) as HistFile:
        with up.recreate(outfile) as OutFile:
            for prop in tqdm(properties, desc = "Properties", position = 0):
                for fs in tqdm(fstates, desc = "Final States", position = 1, leave = False):
                    zpx_norm = zpx_info_2022[fs][0]
                    #zz_norms = zz_2022[fs]
                    for proc in tqdm(procs, desc = "Processes", position = 2, leave = False):

                        if proc != "ZpX":
                            proc_vars = get_vars(proc)
                            for var in tqdm(proc_vars, desc = "Variations", position = 3 , leave = False):
                                
                                key = f"{prop}/SR/{fs}/{proc}{var}"
                                
                                hist = HistFile[key].to_pyroot()

                                if prop == "delRapidity":
                                    hist = hist.Rebin(2, "delRapidity")
                                
                                # if ("ZLZL" in proc) or ("ZLZT" in proc) or ("ZTZT" in proc):
                                #     if var == "":
                                #         hist.Scale(1./hist.Integral())
                                #     else:
                                #         norm_val = zz_norms[var.replace("_", "")]
                                #         hist.Scale(norm_val/hist.Integral())
                                    #norm_hist(hist)
                                
                                OutFile[key] = hist
                        
                        else:
                            if prop == "delRapidity":
                                hist = hist.Rebin(2, "delRapidity")
                            
                            # use previous hist as template to find nbins, low and upper edge
                            nbins = hist.GetNbinsX()
                            low, high = hist.GetBinLowEdge(1), hist.GetBinLowEdge(nbins+1)

                            # Want integral to be 1 so we can uniformly scale it with a rateParam in combine
                            #bin_yield = 1./nbins
                            # For plots we directly put the norm
                            bin_yield = zpx_norm/nbins
                            
                            zpx_hist = ROOT.TH1D(prop, prop, nbins, low, high)

                            if prop == "delRapidity":
                                zpx_hist, zpx_histUp, zpx_histDn = zpx_shape(full_hist, zpx_norm)

                                # zpx_norm_up, zpx_norm_dn = zpx_histUp.Integral()/zpx_norm, zpx_histDn.Integral()/zpx_norm
                                # zpx_hist.Scale(1./zpx_hist.Integral())
                                # zpx_histUp.Scale(zpx_norm_up/zpx_histUp.Integral())
                                # zpx_histDn.Scale(zpx_norm_dn/zpx_histDn.Integral())

                                OutFile[f"{prop}/SR/{fs}/{proc}"]         = zpx_hist
                                OutFile[f"{prop}/SR/{fs}/{proc}_zpxUp"]   = zpx_histUp
                                OutFile[f"{prop}/SR/{fs}/{proc}_zpxDown"] = zpx_histDn
                            
                            else:
                                for bin_idx in range(nbins):
                                    zpx_hist.SetBinContent(bin_idx+1, bin_yield)
                                    zpx_hist.SetBinError(bin_idx+1, 0.) # error will be handled through the rateParam                            

                                OutFile[f"{prop}/SR/{fs}/{proc}"]         = zpx_hist