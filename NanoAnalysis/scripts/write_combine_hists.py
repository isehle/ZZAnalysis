import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT
import uproot as up
import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep

from NanoAnalysis.scripts.histReader import HistReader
from NanoAnalysis.scripts.ZpX_estimation import ZpX

procs = dict(
    ggZZ = ['ggTo2e2mu', 'ggTo2e2tau', 'ggTo2mu2tau', 'ggTo4e', 'ggTo4mu', 'ggTo4tau'],
    ZpX  = ['DY', 'TT', 'WZ'],
    VVV  = ['WWZ', 'WZZ', 'ZZZ'],
    H125 = 'H',
    ZLZL = 'ZLZL',
    ZLZT = 'ZLZT',
    ZTZT = 'ZTZT'
)

#fstates = ['fs_4e', 'fs_4mu', ['fs_2e2mu', 'fs_2mu2e']]
fstates = ['fs_4e', 'fs_4mu', 'fs_2e2mu', 'fs_2mu2e']

pol_samples = ['ZLZL', 'ZLZT', 'ZTZT']
norm_samples = pol_samples + ["ZpX"]

get_vars = lambda proc: ["", "_LHEScaleWeightUp", "_LHEScaleWeightDown", "_LHEPdfWeightUp", "_LHEPdfWeightDown", "_puWeightUp", "_puWeightDown", "_lepIDRecUp", "_lepIDRecDown"] if proc in pol_samples else [""]

norm_hist = lambda hist: hist.Scale(1./hist.Integral())

def th1_to_np(hist):
    #hist = hist.GetValue()

    bin_contents = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
    bin_edges    = np.array([hist.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)])

    return bin_contents, bin_edges

def get_hists(pol_hist_file):
    with up.open(pol_hist_file) as HistFile:
        nominal = HistFile["ZLZL"]["SR"]["delPhi"]["fs_4mu"].to_pyroot()

        QCDUp   = HistFile["ZLZL"]["SR"]["delPhi"]["fs_4mu_LHEScaleWeightUp"].to_pyroot()
        QCDDn   = HistFile["ZLZL"]["SR"]["delPhi"]["fs_4mu_LHEScaleWeightDown"].to_pyroot()

        PDFUp   = HistFile["ZLZL"]["SR"]["delPhi"]["fs_4mu_LHEPdfWeightUp"].to_pyroot()
        PDFDn   = HistFile["ZLZL"]["SR"]["delPhi"]["fs_4mu_LHEPdfWeightDown"].to_pyroot()

        norm_hist(nominal)

        norm_hist(QCDUp)
        norm_hist(QCDDn)

        norm_hist(PDFUp)
        norm_hist(PDFDn)

    return th1_to_np(nominal), th1_to_np(PDFUp), th1_to_np(PDFDn)

def writeFigs(hist_file):
    labels = [r"$q \bar{q} \rightarrow Z_L Z_L$", r"$(q \bar{q} \rightarrow Z_L Z_L)_{QCDScale\_Up}$", r"$(q \bar{q} \rightarrow Z_L Z_L)_{QCDScale\_Dn}$"]

    nom_hist, pdf_up, pdf_dn = get_hists(hist_file)

    hep.style.use("CMS")
    ratio_fig_style = {
        "figsize": [10, 10],
        "gridspec_kw": {
            "height_ratios": [3, 1]
        }
    }

    fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **ratio_fig_style)
    fig.subplots_adjust(hspace=0.07)

    hep.histplot(
        [nom_hist, pdf_up, pdf_dn],
        histtype="step",
        label = labels,
        ax = ax,
        color = ["black", "blue", "orange"]
    )

    ax.legend()
    ax.set_ylabel("Events")

    bin_centers = (nom_hist[1][:-1] + nom_hist[1][1:])/2

    rax.scatter(bin_centers,
                pdf_up[0]/nom_hist[0],
                color = "blue",
                label = "QCDScale_Up/Nom")

    rax.scatter(bin_centers,
            pdf_dn[0]/nom_hist[0],
            color = "orange",
            label = "QCDScale_Dn/Nom")

    #rax.legend()

    rax.axhline(y=1,linestyle="--")
    rax.set_ylim(0.97, 1.03)

    #rax.set_ylabel("gg/qq")
    rax.set_xlabel(r"$cos(\theta_1)$")

    fig.suptitle('QCD Variations, 4mu, 2022EE')

    outfile = f"delPhi_ZLZL_qcdVarRatio_4mu_2022EE.png"
    fig.savefig(outfile)

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

def write_zpx_hist(hist_info, norm, prop):
    hist = ROOT.TH1D(prop, prop, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"]))

    counts_per_bin = norm/float(hist_info["nbinsx"])

    for bin_idx in range(int(hist_info["nbinsx"])):
        hist.SetBinContent(bin_idx+1, counts_per_bin)
    
    return hist

if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--reg", choices=("SR", "OS_NoSIP_HighMass", "OS_NoSIP_MidMass", "OS_NoSIP_LowMass", "SS_NoSIP_HighMass", "SS_NoSIP_MidMass", "SS_NoSIP_LowMass"), default="SR")
    parser.add_argument("--prop", default="mass")
    parser.add_argument("--year", choices=(2022, 2023), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "Full"), default="Full")
    parser.add_argument("--tag", default="")
    parser.add_argument("--infile", default="")
    parser.add_argument("--lumi_tag", default=0, type=int)
    args = vars(parser.parse_args())

    with open("/afs/cern.ch/user/i/iehle/cmssw/CMSSW_14_1_6/src/ZZAnalysis/NanoAnalysis/scripts/hist_config.yaml") as config:
        cfg = yaml.safe_load(config)

    histReader = HistReader(cfg, args)

    hist_info = cfg["hist_info"]["delPhi"]

    #hist_file = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/hists_Z1pt.root"
    hist_file = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/hists_Z1pt_v2.root"
    outpath   = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/delPhi_hists_NormZpX.root"

    # hist_file = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/hists_newLepPtReqs_v3_cp.root"
    # outpath   = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/delPhi_hists_newLepPtReqs_v3_fixPDFVars.root"

    #writeFigs(pol_hist_file)

    # zpx = ZpX()
    # hists, counts, errors = histReader.read_hists_and_counts(hist_file)
    # zpx_info = zpx.get_zpx(hists, errors, fstates)

    with up.open(hist_file) as HistFile:
        with up.recreate(outpath) as OutFile:
            for proc, subprocs in procs.items():
                if isinstance(subprocs, list):
                    variations = get_vars(proc[0])
                    if proc == "ZpX":
                        #hists = zpx.write_hists(zpx_info, hist_info, "delPhi")
                        nbins = int(hist_info["nbinsx"])
                        count_per_bin = 1./nbins
                        for fs in fstates:
                            hist = ROOT.TH1D("delPhi", "delPhi", nbins, float(hist_info["xlow"]), float(hist_info["xhigh"]))
                            for bin_idx in range(nbins):
                                hist.SetBinContent(bin_idx+1, count_per_bin)
                                hist.SetBinError(bin_idx+1, 0.)
                            OutFile[f"{fs}/{proc}"] = hist

                            # OutFile[f"{fs}/{proc}"]       = hists[fs]["Nominal"]
                            # OutFile[f"{fs}/{proc}_sipUp"] = hists[fs]["Up"]
                            # OutFile[f"{fs}/{proc}_sipDown"] = hists[fs]["Down"]
                    else:
                        hists = {}
                        for subproc in subprocs:
                            hists[subproc] = {}
                            for fs in fstates:
                                if type(fs) == list:
                                    hists[subproc][fs[0]] = {}
                                    for var in variations:
                                        hist_1 = HistFile[subproc]["SR"]["delPhi"][f"{fs[0]}{var}"].to_pyroot()
                                        hist_2 = HistFile[subproc]["SR"]["delPhi"][f"{fs[1]}{var}"].to_pyroot()

                                        hist = hist_1 + hist_2

                                        if not hist.Integral() ==  0 and proc in pol_samples: norm_hist(hist)
                                        hists[subproc][fs[0]][var] = hist

                                else:
                                    hists[subproc][fs] = {}
                                    for var in variations:
                                        hist = HistFile[subproc]["SR"]["delPhi"][f"{fs}{var}"].to_pyroot()

                                        if not hist.Integral() ==  0 and proc in pol_samples: norm_hist(hist)
                                        hists[subproc][fs][var] = hist
                    
                        for fs in fstates:
                            fs_key = fs if type(fs) != list else fs[0]
                            for var in variations:
                                new_hist = sum_hist_list(hists, subprocs, fs_key, var, proc)
                                # new_hist.Rebin(4)
                                if new_hist.Integral() == 0:
                                    breakpoint()
                                OutFile[f"{fs_key}/{proc}{var}"] = new_hist
                else:
                    variations = get_vars(proc)
                    for fs in fstates:
                        if type(fs) == list:
                            for var in variations:
                                hist_1 = HistFile[subprocs]["SR"]["delPhi"][f"{fs[0]}{var}"].to_pyroot()
                                hist_2 = HistFile[subprocs]["SR"]["delPhi"][f"{fs[1]}{var}"].to_pyroot()

                                hist = hist_1 + hist_2
                                # if not hist.Integral() ==  0 and proc in pol_samples:
                                #     norm_hist(hist)
                                if hist.Integral() == 0:
                                    hist = fill_empty_hist(hist)
                                elif proc in pol_samples:
                                    norm_hist(hist)
                                
                                # hist.Rebin(4)
                                OutFile[f"{fs[0]}/{proc}{var}"] = hist
                        else:
                            for var in variations:
                                hist = HistFile[subprocs]["SR"]["delPhi"][f"{fs}{var}"].to_pyroot()
                                    
                                # if not hist.Integral() ==  0 and proc in pol_samples:
                                #     norm_hist(hist)
                                if hist.Integral() == 0:
                                    hist = fill_empty_hist(hist)
                                elif proc in pol_samples:
                                    norm_hist(hist)
                                
                                # hist.Rebin(4)
                                OutFile[f"{fs}/{proc}{var}"] = hist
