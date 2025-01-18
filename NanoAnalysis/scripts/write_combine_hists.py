import ROOT
import uproot as up
import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep

pol_hist_file = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/hists_shapeVars.root"

# cosTheta1_path    = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnlyQCDScaleTestv3.root"
# cosTheta3_path    = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnly.root"
# cosThetaStar_path = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnly.root"
# delRapidity_path  = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnly.root"

#out_cosTheta1    = ROOT.TFile.Open(cosTheta1_path, "RECREATE")
# out_cosTheta3    = ROOT.TFile.Open(cosTheta3_path, "RECREATE")
# out_cosThetaStar = ROOT.TFile.Open(cosThetaStar_path, "RECREATE")
# out_delRapidity  = ROOT.TFile.Open(delRapidity_path, "RECREATE")

norm_hist = lambda hist: hist.Scale(1./hist.Integral())

def th1_to_np(hist):
    #hist = hist.GetValue()

    bin_contents = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
    bin_edges    = np.array([hist.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)])

    return bin_contents, bin_edges

def get_hists(pol_hist_file):
    with up.open(pol_hist_file) as InFile:
        nominal = InFile["ZLZL"]["SR"]["cosTheta1"]["fs_4mu"].to_pyroot()
        PDFUp   = InFile["ZLZL"]["SR"]["cosTheta1"]["fs_4mu_LHEScaleWeightUp"].to_pyroot()
        PDFDn   = InFile["ZLZL"]["SR"]["cosTheta1"]["fs_4mu_LHEScaleWeightDown"].to_pyroot()

        norm_hist(nominal)
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

if __name__ == "__main__":
    pol_hist_file = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/hists_shapeVars.root"
    cosTheta1_path    = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnlyTheoryVars.root"

    #writeFigs(pol_hist_file)

    with up.open(pol_hist_file) as InFile:
        with up.recreate(cosTheta1_path) as OutFile:
            for proc in ["ZLZL", "ZLZT", "ZTZT"]:
                for fs in ["fs_4e", "fs_4mu", ["fs_2e2mu", "fs_2mu2e"]]:
                    if type(fs) == list:
                        for var in ["", "_LHEScaleWeightUp", "_LHEScaleWeightDown", "_LHEPdfWeightUp", "_LHEPdfWeightDown"]:
                            hist_1 = InFile[proc]["SR"]["cosTheta1"][f"{fs[0]}{var}"].to_pyroot()
                            hist_2 = InFile[proc]["SR"]["cosTheta1"][f"{fs[1]}{var}"].to_pyroot()

                            hist = hist_1 + hist_2
                            norm_hist(hist)

                            OutFile[f"{fs[0]}/{proc}{var}"] = hist
                    else:
                        for var in ["", "_LHEScaleWeightUp", "_LHEScaleWeightDown", "_LHEPdfWeightUp", "_LHEPdfWeightDown"]:
                            hist = InFile[proc]["SR"]["cosTheta1"][f"{fs}{var}"].to_pyroot()
                            norm_hist(hist)

                            OutFile[f"{fs}/{proc}{var}"] = hist