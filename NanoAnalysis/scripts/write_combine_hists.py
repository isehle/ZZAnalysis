import ROOT
import uproot as up
import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep

procs = dict(
    ggZZ = ['ggTo2e2mu', 'ggTo2e2tau', 'ggTo2mu2tau', 'ggTo4e', 'ggTo4mu', 'ggTo4tau'],
    ZpX  = ['DY', 'TT', 'WZ'],
    VVV  = ['WWZ', 'WZZ', 'ZZZ'],
    H125 = 'H',
    ZLZL = 'ZLZL',
    ZLZT = 'ZLZT',
    ZTZT = 'ZTZT'
)

fstates = ['fs_4e', 'fs_4mu', ['fs_2e2mu', 'fs_2mu2e']]

pol_samples = ['ZLZL', 'ZLZT', 'ZTZT']

get_vars = lambda proc: ["", "_LHEScaleWeightUp", "_LHEScaleWeightDown", "_LHEPdfWeightUp", "_LHEPdfWeightDown"] if proc in pol_samples else [""]

norm_hist = lambda hist: hist.Scale(1./hist.Integral())

def th1_to_np(hist):
    #hist = hist.GetValue()

    bin_contents = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
    bin_edges    = np.array([hist.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)])

    return bin_contents, bin_edges

def get_hists(pol_hist_file):
    with up.open(pol_hist_file) as InFile:
        nominal = InFile["ZLZL"]["SR"]["cosTheta1"]["fs_4mu"].to_pyroot()

        QCDUp   = InFile["ZLZL"]["SR"]["cosTheta1"]["fs_4mu_LHEScaleWeightUp"].to_pyroot()
        QCDDn   = InFile["ZLZL"]["SR"]["cosTheta1"]["fs_4mu_LHEScaleWeightDown"].to_pyroot()

        PDFUp   = InFile["ZLZL"]["SR"]["cosTheta1"]["fs_4mu_LHEPdfWeightUp"].to_pyroot()
        PDFDn   = InFile["ZLZL"]["SR"]["cosTheta1"]["fs_4mu_LHEPdfWeightDown"].to_pyroot()

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

    outfile = f"cosTheta1_ZLZL_qcdVarRatio_4mu_2022EE.png"
    fig.savefig(outfile)


def sum_hist_list(hists, subprocs, fs_key, var, proc):
    hist_list = [hists[subproc][fs_key][var] for subproc in subprocs]
    new_hist = hist_list[0].Clone(proc)
    for hist in hist_list[1:]:
        new_hist.Add(hist)

    return new_hist

if __name__ == "__main__":
    pol_hist_file = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/hists_newLepPtReqs_v2.root"
    outpath       = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta3_hists_newLepPtReqs.root"

    #writeFigs(pol_hist_file)

    with up.open(pol_hist_file) as InFile:
        with up.recreate(outpath) as OutFile:
            for proc, subprocs in procs.items():
                if isinstance(subprocs, list):
                    variations = get_vars(proc[0])
                    hists = {}
                    for subproc in subprocs:
                        hists[subproc] = {}
                        for fs in fstates:
                            if type(fs) == list:
                                hists[subproc][fs[0]] = {}
                                for var in variations:
                                    hist_1 = InFile[subproc]["SR"]["cosTheta3"][f"{fs[0]}{var}"].to_pyroot()
                                    hist_2 = InFile[subproc]["SR"]["cosTheta3"][f"{fs[1]}{var}"].to_pyroot()

                                    hist = hist_1 + hist_2

                                    if not hist.Integral() ==  0 and proc in pol_samples: norm_hist(hist)
                                    hists[subproc][fs[0]][var] = hist

                            else:
                                hists[subproc][fs] = {}
                                for var in variations:
                                    hist = InFile[subproc]["SR"]["cosTheta3"][f"{fs}{var}"].to_pyroot()

                                    if not hist.Integral() ==  0 and proc in pol_samples: norm_hist(hist)
                                    hists[subproc][fs][var] = hist
                    
                    for fs in fstates:
                        fs_key = fs if type(fs) != list else fs[0]
                        for var in variations:
                            new_hist = sum_hist_list(hists, subprocs, fs_key, var, proc)
                            OutFile[f"{fs_key}/{proc}{var}"] = new_hist
                else:
                    variations = get_vars(proc)
                    for fs in fstates:
                        if type(fs) == list:
                            for var in variations:
                                hist_1 = InFile[subprocs]["SR"]["cosTheta3"][f"{fs[0]}{var}"].to_pyroot()
                                hist_2 = InFile[subprocs]["SR"]["cosTheta3"][f"{fs[1]}{var}"].to_pyroot()

                                hist = hist_1 + hist_2
                                if not hist.Integral() ==  0 and proc in pol_samples: norm_hist(hist)

                                OutFile[f"{fs[0]}/{proc}{var}"] = hist
                        else:
                            for var in variations:
                                hist = InFile[subprocs]["SR"]["cosTheta3"][f"{fs}{var}"].to_pyroot()
                                if not hist.Integral() ==  0 and proc in pol_samples: norm_hist(hist)

                                OutFile[f"{fs}/{proc}{var}"] = hist