import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT
import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep

base_path = "/eos/user/i/iehle/Analysis/rootFiles/PrivateProduction/2022/EraEFG"

all_paths = dict(
    qq = dict(
        ll = os.path.join(base_path, "qqZLZL_goodSeeds_hadd_Skim.root"),
        tt = os.path.join(base_path, "qqZTZT_goodSeeds_hadd_Skim.root")
    ),
    gg = dict(
        ll = os.path.join(base_path, "ggZLZL_4l_5f_Skim.root"),
        tt = os.path.join(base_path, "ggZTZT_4l_5f_Skim.root")
    )
)

x_labels = dict(
    cosTheta1 = r"$cos(\theta_1)$",
    cosTheta3 = r"$cos(\theta_3)$",
    cosThetaStar = r"$cos(\theta^*)$",
    delRapidity = r"$|\Delta (y_{Z Z'})|$"
)

def get_df(path):

    df = ROOT.RDataFrame("Events", path)
    df = df.Filter("HLT_passZZ4l")

    Runs = ROOT.RDataFrame("Runs", path)

    return df.Define("genEventSumw",str(Runs.Sum("genEventSumw").GetValue()))

def write_weight(df, reg):
    cand     = lambda reg: "ZZCand" if reg=="SR" else "ZLLCand"
    reg_filt = lambda reg: "{}.at(0) == 1".format(reg)
    
    # dataMCWeight stored as RVec but only ever has one entry
    cand_weight = "{}_dataMCWeight.at(0)*".format(cand(reg))

    df = df.Define("weight", "{}overallEventWeight/genEventSumw".format(cand_weight))

    return df.Filter(reg_filt(reg))

def get_hists():
    hists = {}
    for state in ["ll", "tt"]:
        hists[state] = {}
        
        qq_path, gg_path = all_paths["qq"][state], all_paths["gg"][state]

        qq_df = get_df(qq_path)
        gg_df = get_df(gg_path)

        qq_df = write_weight(qq_df, "SR")
        gg_df = write_weight(gg_df, "SR")

        #for var in ["cosTheta1", "cosTheta3", "cosThetaStar"]:
        for var in ["delRapidity"]:
            #qq_hist = qq_df.Histo1D(("","",21,-1.,1.1), "ZZCand_"+var, "weight")
            #gg_hist = gg_df.Histo1D(("","",21,-1.,1.1), "ZZCand_"+var, "weight")

            qq_hist = qq_df.Histo1D(("","",41,0.,4.1), "ZZCand_"+var, "weight")
            gg_hist = gg_df.Histo1D(("","",41,0.,4.1), "ZZCand_"+var, "weight")

            qq_hist.Scale(1./qq_hist.Integral())
            gg_hist.Scale(1./gg_hist.Integral())

            hists[state][var] = dict(
                qq = th1_to_np(qq_hist),
                gg = th1_to_np(gg_hist)
            )
    
    return hists

def th1_to_np(hist):
    hist = hist.GetValue()

    bin_contents = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])
    bin_edges    = np.array([hist.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)])

    return bin_contents, bin_edges

def writeFigs(hists):
    ll_label = [r"$gg \rightarrow Z_L Z_L$", r"$qq \rightarrow Z_L Z_L$"]
    tt_label = [r"$gg \rightarrow Z_T Z_T$", r"$qq \rightarrow Z_T Z_T$"]
    for state in hists.keys():
        labels = ll_label if state == "ll" else tt_label
        for var in hists[state].keys():
            qq_hist = hists[state][var]['qq']
            gg_hist = hists[state][var]['gg']

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
                [gg_hist, qq_hist],
                histtype="step",
                label = labels,
                ax = ax,
            )

            ax.legend()
            ax.set_ylabel("Events")

            bin_centers = (qq_hist[1][:-1] + qq_hist[1][1:])/2

            rax.scatter(bin_centers,
                        gg_hist[0]/qq_hist[0],
                        color = "black")

            rax.axhline(y=1,linestyle="--")
            rax.set_ylim(0, 2)

            rax.set_ylabel("gg/qq")
            rax.set_xlabel(x_labels[var])

            outfile = f"{var}_gg_qq_{state}_ratio_2022EE.png"
            fig.savefig(outfile)


            



if __name__ == "__main__":
    hists = get_hists()
    writeFigs(hists)
        