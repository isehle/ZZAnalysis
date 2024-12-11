import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import datetime
from pathlib import Path

import ROOT
import uproot as up
from tqdm import tqdm

import matplotlib.pyplot as plt
import mplhep as hep

from NanoAnalysis.scripts.helper_functions import *

class HistPlotter:
    def __init__(self, cfg, args):
        self.cfg  = cfg
        self.args = args

        if args["infile"] == "":
            self.infile = os.path.join(cfg["output"]["base_dir"], "rootFiles", str(args["year"]), args["era"], "hists{}.root".format(args["tag"]))
        else:
            self.infile = args["infile"]

        self.fill_colors = [cfg["plot_styling"]["mc_colors"][proc]["fill"] for proc in cfg["datasets"]["MC_Procs"].keys()]
        self.line_colors = [cfg["plot_styling"]["mc_colors"][proc]["line"] for proc in cfg["datasets"]["MC_Procs"].keys()]

        self.ratio_fig_style = cfg["plot_styling"]["ratio_fig_style"]
        self.hatch_style     = cfg["plot_styling"]["hatch_style"]
        self.errorbar_style  = cfg["plot_styling"]["errorbar_style"]

        self.bin_edges   = []
        self.bin_centers = []

        self.is_empty = lambda a_dict, key: a_dict[key] == {}

        self.labels = dict(
            MC = dict(
                ggZZ   = r"$gg \rightarrow ZZ$",
                ZZ_NLO = r"$(q\bar{q} \rightarrow ZZ)_{NLO}$",
                DY     = r"$DY$",
                TT     = r"$t\bar{t}$",
                WZ     = r"$WZ$",
                H      = r"$H$",
                VVV    = r"$VVV$"
        ),
            Pol = dict(
                ZLZL  = r"$q\bar{q} \rightarrow Z_L Z_L$",
                ZLZT  = r"$q\bar{q} \rightarrow Z_L Z_T$",
                ZTZT  = r"$q\bar{q} \rightarrow Z_T Z_T$",
                ZZ_LO = r"$(q\bar{q} \rightarrow ZZ)_{LO}$"                
        ),
            Data = "Data"
        )

        self.pol_colors = ["magenta", "lime", "red", "black"]

    def set_outfile(self, reg, fstate, prop):
        if self.args["year"] != -1:
            outdir = os.path.join(self.cfg["output"]["plot_dir"], str(self.args["year"]), self.args["era"], reg, fstate)
        else:
            outdir = os.path.join(self.cfg["output"]["plot_dir"], "Full", reg, fstate)
        Path(outdir).mkdir(parents=True, exist_ok=True)
        self.outfile = os.path.join(outdir, prop+self.args["tag"]+"_countTest.png")

    def set_lumi_tag(self):
        if self.args["lumi_tag"] == 0:
            lumi_tag = self.cfg["datasets"]["year_"+str(self.args["year"])][self.args["era"]]["Lumi"]
            self.lumi_tag = round(float(lumi_tag)*1e-3)
        else:
            self.lumi_tag = self.args["lumi_tag"]

    def set_prop_info(self, prop, reg):
        self.prop_info = self.cfg["hist_info"][prop]
        
        if "mass" in prop:
            mass_reg = reg.split("_")[-1] if reg != "SR" else "SR"
            if mass_reg == "HighMass": mass_reg = "SR"
            self.ylabel = to_raw_string(self.prop_info[mass_reg]["ylabel"])
            self.prop_info = self.cfg["hist_info"][prop][mass_reg]
        else:
            self.ylabel = to_raw_string(self.prop_info["ylabel"])
        
        self.xlabel = to_raw_string(self.prop_info["xlabel"])

    def setBins(self, hists):
        if self.bin_edges == []:
            key = list(hists.keys())[0]
            self.bin_edges = hists[key][1]
            self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:])/2

    def adjustBinEdges(self, hist_list):
        new_list = []
        for hist, edges in hist_list:
            edges[0] = edges[1] - (edges[2] - edges[1])
            edges[-1] = edges[-2] + (edges[-2] - edges[-3])
            new_list.append((hist, edges))
        return new_list

    def addCounts(self, labels, counts, errors=[]):

        if len(errors) == 0:
            err = np.sqrt(counts)
            new_lab = labels + ": " + str(counts) + r" $\pm$ " + str(round(err, 2))
            return new_lab
        else:
            new_labels = []
            for proc in labels.keys():
                init_lab = labels[proc]
                count    = counts[proc]
                bin_errs = errors[proc]
                err      = np.sqrt(np.square(bin_errs).sum())
                
                new_lab = init_lab + ": " + str(count) + r" $\pm$ " + str(round(err, 2))
                new_labels.append(new_lab)
            return new_labels

    def writeFigs(self):
        hep.style.use("CMS")

        if not self.prop_info["Blind"]:
            self.fig, (self.ax, self.rax) = plt.subplots(2, 1, sharex=True, **self.ratio_fig_style)
            self.fig.subplots_adjust(hspace=0.07)
        else:
            self.fig, self.ax = plt.subplots()

    def cms_label(self):
        if self.args["year"] != -1:
            hep.cms.label(label="Work in Progress", year=self.args["year"], lumi = self.lumi_tag, com=13.6, data=True, ax=self.ax)
        else:
            hep.cms.label(label="Work in Progress", lumi = self.lumi_tag, com=13.6, data=True, ax=self.ax)

    def draw_mc_hists(self, hists, labels):
        ordered_hists = [hists[proc] for proc in self.labels["MC"].keys()]
        adjusted_hists = self.adjustBinEdges(ordered_hists)
        hep.histplot(
            adjusted_hists,
            stack=True,
            histtype='fill',
            label = labels,
            ax = self.ax,
            color = self.fill_colors,
            edgecolor = self.line_colors
        )

    def draw_pol_hists(self, hists, labels):
        ordered_hists = [hists[proc] for proc in self.labels["Pol"].keys()]
        adjusted_hists = self.adjustBinEdges(ordered_hists)
        hep.histplot(
            adjusted_hists,
            histtype="step",
            label = labels,
            ax = self.ax,
            color = self.pol_colors
        )

    def draw_data_hist(self, hist, label):
        self.ax.errorbar(self.bin_centers,
                    hist[0],
                    yerr = np.sqrt(hist[0]),
                    color = "black",
                    fmt = "o",
                    label = label,
                    markersize = 3)

    def draw_ratio(self, data_hist):
        self.rax.fill_between(x=self.bin_centers, y1= 1 - self.total_mc_errors/self.total_mc_counts, y2 = 1 + self.total_mc_errors/self.total_mc_counts, step='mid', **self.hatch_style)
        self.rax.errorbar(x=self.bin_centers, y=data_hist[0]/self.total_mc_counts, yerr=np.sqrt(data_hist[0])/self.total_mc_counts, **self.errorbar_style)

        self.rax.set_ylim(0, 2)
        self.rax.set_ylabel('Data / MC')
        self.rax.set_xlabel(self.xlabel)
        self.rax.autoscale(axis='x', tight=True)

    def set_stackCountsErrs(self, mc_hists, mc_errors):
        mc_count_arr = np.array([mc_hists[key][0] for key in self.cfg["plot_styling"]["mc_colors"].keys()])
        self.total_mc_counts = np.sum(mc_count_arr, axis=0)
        self.total_mc_counts[ self.total_mc_counts < 0] = 0
        self.total_mc_errors = np.sqrt(np.sum([err**2 for err in mc_errors.values()], axis=0))

    def plotter(self, prop, reg, fs, hists, counts, errors):
        # Set up figures, styling, names etc
        self.set_prop_info(prop, reg)
        self.set_lumi_tag()    
        self.writeFigs()
        self.cms_label()
        self.set_outfile(reg, fs, prop)

        self.draw_mc   = not self.is_empty(hists, "MC")
        self.draw_pol  = not self.is_empty(hists, "Pol")
        self.draw_data = not self.is_empty(hists, "Data") and not self.prop_info["Blind"]

        if self.draw_mc:
            mc_hists  = hists["MC"]
            mc_counts = counts["MC"]
            mc_errors = errors["MC"]

            mc_labels = self.addCounts(self.labels["MC"], mc_counts, mc_errors)

            # Bin centers and edges
            self.setBins(mc_hists)

            # Draw MC hists
            self.draw_mc_hists(mc_hists, mc_labels)

            # Draw MC Err
            self.set_stackCountsErrs(mc_hists, mc_errors)
            self.ax.fill_between(x=self.bin_centers, y1 = self.total_mc_counts - self.total_mc_errors, y2 = self.total_mc_counts + self.total_mc_errors, label = "Stat. Unc.", step='mid', **self.hatch_style)
        
        if self.draw_pol:
            pol_hists  = hists["Pol"]
            pol_counts = counts["Pol"]
            pol_errors = errors["Pol"]

            pol_labels = self.addCounts(self.labels["Pol"], pol_counts, pol_errors)

            self.setBins(pol_hists)

            self.draw_pol_hists(pol_hists, pol_labels)

        if self.draw_data:
            data_hist   = hists["Data"]
            data_counts = counts["Data"]

            data_label = self.addCounts(self.labels["Data"], data_counts)

            self.setBins(data_hist)

            self.draw_data_hist(data_hist, data_label)

            # Adjust plot size to fit legend
            max_bin_count = max([max(hist[0]) for hist in mc_hists.values()]+[max(data_hist[0])])
            self.ax.set_ylim(0, max_bin_count*2.5) 

        self.ax.legend(ncol=2, fontsize="x-small")

        # Labels and ratio (xlabel drawn on rax if unblinded)
        self.ax.set_ylabel(self.ylabel)
        if self.draw_data:
            self.draw_ratio(data_hist)
        else:
            self.ax.set_xlabel(self.xlabel)

        hep.rescale_to_axessize(self.ax, 10, 10/1.62)

        self.fig.savefig(self.outfile)

        self.bin_edges   = []
        self.bin_centers = []


if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--reg", choices=("SR", "OS_NoSIP_HighMass", "OS_NoSIP_MidMass", "OS_NoSIP_LowMass", "SS_NoSIP_HighMass", "SS_NoSIP_MidMass", "SS_NoSIP_LowMass"), default="SR")
    parser.add_argument("--prop", default="mass")
    parser.add_argument("--year", choices=(2022, 2023), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "B"), default="EFG")
    parser.add_argument("--tag", default="")
    args = vars(parser.parse_args())

    with open("/afs/cern.ch/user/i/iehle/cmssw/CMSSW_13_3_3/src/ZZAnalysis/NanoAnalysis/scripts/hist_config.yaml") as config:
        cfg = yaml.safe_load(config)

    plots = HistPlotter(cfg, args)