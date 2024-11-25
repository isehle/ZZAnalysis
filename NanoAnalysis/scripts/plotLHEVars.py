import matplotlib.pyplot as plt
import mplhep as hep

import uproot as up
import numpy as np

import yaml

import os
import datetime
from pathlib import Path

import helper_functions

def adjustBinEdges(hist_list, norm=False):
    new_list = []
    for hist, edges in hist_list:
        edges[0] = edges[1] - (edges[2] - edges[1])
        edges[-1] = edges[-2] + (edges[-2] - edges[-3])
        if norm: hist = hist/hist.sum()
        new_list.append((hist, edges))
    return new_list

def plotter(var, pol_hists, outfile, **kwargs):    

    pol_labels = [
        r"$q\bar{q} \rightarrow Z_L Z_L$",
        r"$q\bar{q} \rightarrow Z_L Z_T$",
        r"$q\bar{q} \rightarrow Z_T Z_L$",
        r"$q\bar{q} \rightarrow Z_{L(T)} Z_{T(L)}$",
        r"$q\bar{q} \rightarrow Z_T Z_T$",
    ]

    hep.style.use("CMS")

    fig, ax = plt.subplots()

    hep.histplot(
        pol_hists,
        histtype="step",
        label = pol_labels,
        ax = ax,
        color = ["red", "orange", "green", "black", "blue"]
    )

    ax.legend(ncol=2, fontsize="x-small")

    ax.set_xlabel(var)

    hep.rescale_to_axessize(ax, 10, 10/1.62)

    fig.savefig(outfile)

if __name__ == "__main__":
    infile = "qqZZ_2e2mu_LHE_vars.root"

    states = ["ZLZL", "ZLZT", "ZTZL", "ZTZT"]

    variables = ["cosTheta1", "cosTheta3", "cosThetaStar", "delRapidity", "delPhiStar"]

    with up.open(infile) as MyFile:
        for var in variables:
            zlzl, zlzt, ztzl, ztzt = MyFile[var+"_ZLZL"].to_numpy(), MyFile[var+"_ZLZT"].to_numpy(), MyFile[var+"_ZTZL"].to_numpy(), MyFile[var+"_ZTZT"].to_numpy()

            mixed = (zlzt[0]+ztzl[0])/2 # For normalized only!!

            mixed_pol = (mixed, zlzt[1])

            outfile = var + "_LHE_Normalized_wMixed_v2.png"
            plotter(var, [zlzl, zlzt, ztzl, mixed_pol, ztzt], outfile)
