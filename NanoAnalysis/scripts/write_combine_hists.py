import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT
import uproot as up
import numpy as np

import json

from NanoAnalysis.scripts.ZpX_estimation import ZpX

procs       = ["ggZZ", "ZpX", "VVV", "H", "ZLZL", "ZLZT", "ZTZT"]
pol_samples = ['ZLZL', 'ZLZT', 'ZTZT']
fstates     = ['fs_2x2e', 'fs_2x2mu']

properties = [
    "cosTheta1",
    "cosTheta3",
    "cosThetaStar",
    "delPhiStar",
    "delPhi",
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

if __name__ == "__main__":
    from tqdm import tqdm

    hist_file = "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_24_06_25.root"
    outfile   = "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_normPolsFlatZpX_24_06_25.root"

    with up.open(hist_file) as HistFile:
        with up.recreate(outfile) as OutFile:
            for prop in tqdm(properties, desc = "Properties", position = 0):
                for fs in tqdm(fstates, desc = "Final States", position = 1, leave = False):
                    
                    for proc in tqdm(procs, desc = "Processes", position = 2, leave = False):

                        if proc != "ZpX":
                            proc_vars = get_vars(proc)
                            for var in tqdm(proc_vars, desc = "Variations", position = 3 , leave = False):
                                
                                key = f"{prop}/SR/{fs}/{proc}{var}"
                                
                                hist = HistFile[key].to_pyroot()
                                
                                if ("ZLZL" in proc) or ("ZLZT" in proc) or ("ZTZT" in proc):
                                    norm_hist(hist)
                                
                                OutFile[key] = hist
                        
                        else:
                            
                            # use previous hist as template to find nbins, low and upper edge
                            nbins = hist.GetNbinsX()
                            low, high = hist.GetBinLowEdge(1), hist.GetBinLowEdge(nbins+1)

                            # Want integral to be 1 so we can uniformly scale it with a rateParam in combine
                            bin_yield = 1./nbins
                            
                            zpx_hist = ROOT.TH1D(prop, prop, nbins, low, high)
                            for bin_idx in range(nbins):
                                zpx_hist.SetBinContent(bin_idx+1, bin_yield)
                                zpx_hist.SetBinError(bin_idx+1, 0.) # error will be handled through the rateParam                            

                            OutFile[f"{prop}/SR/{fs}/{proc}"] = zpx_hist