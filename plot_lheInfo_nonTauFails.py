import json
import re

import uproot as up
import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep

with open('lheInfo_cloneBranchFailures_noTau_v4.json', 'r') as infile:
    fail_data = json.load(infile)

with open('lheInfo_cloneBranchPasses_noTau.json', 'r') as inPasses:
    pass_data = json.load(inPasses)

# with up.open("AnalysisStep/data/LeptonEffScaleFactors/SF2D_postEE_RMS.root") as egamma_id:
#     egamma_id_hist = egamma_id["EGamma_SF2D"].to_numpy()

x_labels = dict(
    el_pt = r"${p_T}^e$",
    el_eta = r"${\eta}^e$",
    mu_pt = r"${p_T}^{\mu}$",
    mu_eta = r"${\eta}^{\mu}$"
)

hep.style.use("CMS")
# for part in ["el", "mu"]:
#     for var in ["pt", "eta"]:
#         fails  = fail_data[f"{part}_{var}"]
#         passes = pass_data[f"{part}_{var}"]

#         if "pt" in var:
#             fail_hist = np.histogram(fails, bins=100)
#             pass_hist = np.histogram(passes, bins=100)
#         else:
#             fail_hist = np.histogram(fails, bins=20)
#             pass_hist = np.histogram(passes, bins=20)

#         ratios = [p/(f+p) for f, p in zip(fail_hist[0], pass_hist[0])]

#         fig, (ax, rax) = plt.subplots(2, 1, sharex=True, figsize=(10,10), gridspec_kw = {"height_ratios": (3, 1)})
#         hep.cms.label(year = 2022, com = 13.6, data = False, ax = ax)

#         hep.histplot(
#             [fail_hist[0], pass_hist[0]],
#             bins = fail_hist[1],
#             histtype = "fill",
#             label = ["Failures", "Passes"],
#             ax = ax,
#             alpha = 0.6
#         )

#         ax.legend()

#         bin_centers = (fail_hist[1][:-1] + fail_hist[1][1:])/2

#         rax.scatter(bin_centers, ratios)
#         rax.set_ylim(0, 1.1)
#         rax.axhline(1, linestyle="--")
#         rax.set_xlabel(x_labels[f"{part}_{var}"])
#         rax.set_ylabel("Pass/Total")

#         fig.savefig(f"{part}_{var}_zzFiller_ratioHists.png")

pass_hist_pt_eta = np.histogram2d(pass_data["mu_eta"], pass_data["mu_pt"], bins=(20, 100))
fail_hist_pt_eta = np.histogram2d(fail_data["mu_eta"], fail_data["mu_pt"], bins=(20, 100))

ratios = pass_hist_pt_eta[0]/(pass_hist_pt_eta[0] + fail_hist_pt_eta[0])
ratio_hist_pt_eta = (ratios, pass_hist_pt_eta[1], pass_hist_pt_eta[2])

fig, ax = plt.subplots()
hep.cms.label(year = 2022, com = 13.6, data = False, ax = ax)
outfile = "mu_pT_eta_2DHist_lheInfo_zzFillerEff_noTau.png"

hep.hist2dplot(
    ratio_hist_pt_eta,
    ax = ax,
    cmin = 0.,
    cmax = 1.,
)

ax.set_xlabel(r"$\eta_{\mu}$")
ax.set_ylabel(r"${p_T}^{\mu}$")

fig.savefig(outfile)


# with up.open("AnalysisStep/data/LeptonEffScaleFactors/SF2D_postEE_RMS.root") as egamma_id:
#     egamma_id_hist = egamma_id["EGamma_SF2D"].to_numpy()

# hep.style.use("CMS")

# hist_pt_eta = np.histogram2d(data["el_eta"], data["el_pt"], (egamma_id_hist[1], egamma_id_hist[2]))

# fig, ax = plt.subplots()
# hep.cms.label(year = 2022, com=13.6, data=False, ax=ax)

# outfile = "el_pT_eta_2hist_lheInfo_cloneBranchFailures_noTau_goodBins.png"

# hep.hist2dplot(
#     hist_pt_eta,
#     ax = ax
# )

# ax.set_xlabel(r"$\eta_{e}$")
# ax.set_ylabel(r"${p_T}^{e}$")

# fig.savefig(outfile)



