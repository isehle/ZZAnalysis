import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import datetime
from pathlib import Path

import yaml
import json

import ROOT
import uproot as up
from tqdm import tqdm

import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep

from copy import deepcopy

plot_cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/plot_cfg.yaml")

zpx_paths = {
    "2022": "ZpX_info_2022_Full_reweightNonDegEvents_22_09_25.json",
    "2023": "ZpX_info_2023_Full_2x2e_2x2mu_4l.json"
}

class HistPlotter:
    def __init__(self, year, era, tag, lumi, all_hists, all_counts, all_errors, zpx_from_mc, cfg_path=plot_cfg_path):
        self.year = year
        self.era  = era
        self.tag  = tag

        self.lumi = lumi

        self.all_hists  = all_hists
        self.all_counts = all_counts
        self.all_errors = all_errors

        self.zpx_from_mc = zpx_from_mc
        self._read_zpx()

        self._set_cfg(cfg_path)

        self.pol_colors      = self.cfg["pol_colors"]

        self.ratio_fig_style = self.cfg["ratio_fig_style"]
        self.hatch_style     = self.cfg["hatch_style"]
        self.errorbar_style  = self.cfg["errorbar_style"]

        self.pol_ratio_style = deepcopy(self.errorbar_style)
        del self.pol_ratio_style["color"]

        self.ang_vars = ["cosTheta1", "cosTheta3", "cosThetaStar", "delRapidity", "delPhi", "delPhiStar"]

        self.bin_edges   = []
        self.bin_centers = []

        self.is_empty = lambda a_dict, key: a_dict[key] == {}

        self.sr_labels = dict(
            MC = dict(
                ZpX    = r"$Z+X$",
                VVV    = r"$VVV$",
                ggZZ   = r"$gg \rightarrow ZZ$",
                ZZ_NLO = r"$(q\bar{q} \rightarrow ZZ)_{NNLO}$",
        ),
            Pol = dict(
                ZLZL  = r"$q\bar{q} \rightarrow Z_L Z_L$",
                ZLZT  = r"$q\bar{q} \rightarrow Z_L Z_T$",
                ZTZT  = r"$q\bar{q} \rightarrow Z_T Z_T$",
                ZZ_LO = r"$(q\bar{q} \rightarrow ZZ)_{LO}$"                
        ),
            Data = "Data"
        )

        self.lm_labels = dict(
            MC = dict(
                ZpX    = r"$Z+X$",
                H      = r"$H$",
                VVV    = r"$VVV$",
                ggZZ   = r"$gg \rightarrow ZZ$",
                ZZ_NLO = r"$(q\bar{q} \rightarrow ZZ)_{NNLO}$",
        ),
            Pol = dict(
                ZLZL  = r"$q\bar{q} \rightarrow Z_L Z_L$",
                ZLZT  = r"$q\bar{q} \rightarrow Z_L Z_T$",
                ZTZT  = r"$q\bar{q} \rightarrow Z_T Z_T$",
                ZZ_LO = r"$(q\bar{q} \rightarrow ZZ)_{LO}$"                
        ),
            Data = "Data"
        )

        self.cr_labels = dict(
            MC = dict(
                ZpX    = r"$Z+X$",
                VVV    = r"$VVV$",
                ggZZ   = r"$gg \rightarrow ZZ$",
                ZZ_NLO = r"$(q\bar{q} \rightarrow ZZ)_{NNLO}$",
        ),
            Pol = dict(
                ZLZL  = r"$q\bar{q} \rightarrow Z_L Z_L$",
                ZLZT  = r"$q\bar{q} \rightarrow Z_L Z_T$",
                ZTZT  = r"$q\bar{q} \rightarrow Z_T Z_T$",
                ZZ_LO = r"$(q\bar{q} \rightarrow ZZ)_{LO}$"                
        ),
            Data = "Data"
        )

    def _set_procs(self, reg):
        if reg == "SR":
            self.labels = self.sr_labels
        elif "LowMass" in reg:
            self.labels = self.lm_labels
        else:
            self.labels = self.cr_labels

        self.fill_colors = [self.cfg["proc_info"][proc]["fill"] for proc in self.labels["MC"].keys()]
        self.line_colors = [self.cfg["proc_info"][proc]["line"] for proc in self.labels["MC"].keys()]

    def _set_cfg(self, cfg_path):
        with open(cfg_path) as config:
            self.cfg = yaml.safe_load(config)

    def _to_raw_string(self, s):
        return s.encode('unicode_escape').decode('utf-8')

    def _set_prop_info(self, prop, reg):
        hist_cfg = os.path.join(parent_dir, "NanoAnalysis/scripts/hist_cfg.yaml")
        with open(hist_cfg) as config:
            hist_cfg = yaml.safe_load(config)

        self.prop_info = hist_cfg[prop]
        if "mass" in prop:
            mass_reg = "SR" if reg == "SR" or "HighMass" in reg else reg.replace("OSSIP","").replace("SSSIP", "") 
            self.prop_info = self.prop_info[mass_reg]

        self.xlabel = self._to_raw_string(self.prop_info["xlabel"])
        self.ylabel = self._to_raw_string(self.prop_info["ylabel"])
        
        self.blind  = self.prop_info["Blind"] and reg not in ["SS", "HighMassSSRelaxed"]

    def _setBins(self, hists):
        if isinstance(hists, dict):
            key = list(hists.keys())[0]
            self.bin_edges = hists[key][1]
        else:
            self.bin_edges = hists[1]

        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:])/2

        self.adjusted_bin_centers = deepcopy(self.bin_centers)
        if self.bin_centers[-1] == np.inf:
            # Only works when bins are equal length!
            self.adjusted_bin_centers[-1] = self.bin_centers[-2] + (self.bin_centers[-2] - self.bin_centers[-3])

    def _read_zpx(self):
        zpx_info = {}
        for year, path in zpx_paths.items():
            with open(zpx_paths[year]) as json_file:
                zpx_info[year] = json.load(json_file)

        self.zpx_info = zpx_info

    def adjustBinEdges(self, hist_list):
        new_list = []
        for hist, edges in hist_list:
            edges[0] = edges[1] - (edges[2] - edges[1])
            edges[-1] = edges[-2] + (edges[-2] - edges[-3])
            new_list.append((hist, edges))
        return new_list

    def addCounts(self, labels, counts, errors):
        new_labels = []
        for proc in labels.keys():
            init_lab = labels[proc]
            count    = counts[proc]

            if proc == "ZpX" and not self.zpx_from_mc:
                if self.year != "Full":
                    count, err = self.zpx_info[self.year]["N_ZpX_MidMass"][self.fs]
                else:
                    count_22, err_22 = self.zpx_info["2022"]["N_ZpX_MidMass"][self.fs]
                    count_23, err_23 = self.zpx_info["2023"]["N_ZpX_MidMass"][self.fs]
                    
                    count, err       = count_22+count_23, err_22+err_23
            else:
                bin_errs = errors[proc]
                err      = np.sqrt(np.square(bin_errs).sum())

            new_lab = init_lab + ": " + str(round(count,2)) + r" $\pm$ " + str(round(err, 2))
            new_labels.append(new_lab)
        return new_labels

    def writeFigs(self):
        hep.style.use("CMS")

        self.fig, (self.ax, self.rax) = plt.subplots(2, 1, sharex=True, **self.ratio_fig_style)
        self.fig.subplots_adjust(hspace=0.07)

    def cms_label(self):
        if self.year != "Full":
            hep.cms.label(label="Private Work", year=self.year, lumi = round(self.lumi*1e-3), com=13.6, data= not self.blind, ax=self.ax)
        else:
            hep.cms.label(label="Private Work", lumi = round(self.lumi*1e-3), com=13.6, data= not self.blind, ax=self.ax)
 
    def draw_mc_hists(self, hists, labels):

        ordered_hists = [hists[proc] for proc in self.labels["MC"].keys()]

        new_hists = []
        for hist in ordered_hists:
            hist[1][hist[1] == -np.inf] = 0.
            hist[1][hist[1] == np.inf]  = hist[1][-2] + (hist[1][-2] - hist[1][-3])
            new_hists.append(hist)

        hep.histplot(
            new_hists,
            stack=True,
            histtype='fill',
            label = labels,
            ax = self.ax,
            color = self.fill_colors,
            edgecolor = self.line_colors,
            flow = "sum"
        )

    def add_pols(self, hists, counts, errors):
        hists_sum  = (hists["ZLZL"][0]+hists["ZLZT"][0]+hists["ZTZT"][0], hists["ZLZL"][1])
        counts_sum = counts["ZLZL"] + counts["ZLZT"] + counts["ZTZT"]
        errs_sum = np.sqrt(errors["ZLZL"]**2 + errors["ZLZT"]**2 + errors["ZTZT"]**2)

        return hists_sum, counts_sum, errs_sum

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

    def draw_pol_sum(self, hists, counts, errors):
        hist_sum, count_sum, err_sum = self.add_pols(hists["Pol"], counts["Pol"], errors["Pol"])
        count_err = np.sqrt((err_sum**2).sum())
        
        label_sum = r"$q\bar{q} \rightarrow \sum_{\lambda \lambda '}{Z_{\lambda}Z_{\lambda '}}$"
        label_sum += ": " + str(round(count_sum, 2)) + r" $\pm$ " + str(round(count_err, 2))
        
        adjusted_hist_sum = self.adjustBinEdges([hist_sum])
        hep.histplot(
            adjusted_hist_sum,
            histtype="step",
            label = label_sum,
            ax = self.ax,
            color = "black",
            linestyle = "--"
        )

        self.pol_sum_hist, self.pol_sum_err = hist_sum, err_sum

    def draw_data_hist(self, hist, label):
        self.ax.errorbar(self.adjusted_bin_centers,
                    hist[0],
                    yerr = np.sqrt(hist[0]),
                    color = "black",
                    fmt = "o",
                    label = label,
                    markersize = 3)

    def set_stackCountsErrs(self, mc_hists, mc_errors):
        mc_count_arr = np.array([mc_hists[key][0] for key in self.labels["MC"].keys()])
        mc_errors = np.array([mc_errors[key] for key in self.labels["MC"].keys()])

        self.total_mc_counts = np.sum(mc_count_arr, axis=0)
        self.total_mc_counts[ self.total_mc_counts < 0] = 0

        self.total_mc_errors = np.sqrt(np.sum([err**2 for err in mc_errors], axis=0))

    def get_pol_ratios(self, hists, errors):
        zz_lo, zz_nnlo         = hists["Pol"]["ZZ_LO"][0], hists["MC"]["ZZ_NLO"][0]
        zz_lo_err, zz_nnlo_err = errors["Pol"]["ZZ_LO"], errors["MC"]["ZZ_NLO"]

        ratio_lo_nnlo     = zz_lo/zz_nnlo
        ratio_lo_nnlo_err = ratio_lo_nnlo*np.sqrt((zz_lo_err/zz_lo)**2 + (zz_nnlo_err/zz_nnlo)**2)

        ratio_polsum     = self.pol_sum_hist[0]/zz_lo
        ratio_polsum_err = ratio_polsum*np.sqrt((self.pol_sum_err/self.pol_sum_hist[0])**2 + (zz_lo_err/zz_lo)**2)

        return ratio_lo_nnlo, ratio_lo_nnlo_err, ratio_polsum, ratio_polsum_err

    def draw_ratio(self, hists, errors):

        if self.draw_data:
            data_hist = hists["Data"]

            self.rax.fill_between(x=self.bin_edges[:-1], y1= 1 - self.total_mc_errors/self.total_mc_counts, y2 = 1 + self.total_mc_errors/self.total_mc_counts, step='post', **self.hatch_style)
            self.rax.errorbar(x=self.adjusted_bin_centers, y=data_hist[0]/self.total_mc_counts, yerr=np.sqrt(data_hist[0])/self.total_mc_counts, **self.errorbar_style)
            self.rax.set_ylabel('Data / MC')
        else:
            ratio_lo_nnlo, ratio_lo_nnlo_err, ratio_polsum, ratio_polsum_err = self.get_pol_ratios(hists, errors)

            self.rax.errorbar(x=self.adjusted_bin_centers, y = ratio_lo_nnlo, yerr = ratio_lo_nnlo_err, label = r"LO/NNLO", **self.pol_ratio_style)
            self.rax.errorbar(x=self.adjusted_bin_centers, y = ratio_polsum, yerr = ratio_polsum_err, label = r"$\sum Z_{\lambda} Z_{\lambda '}/ZZ_{U}$", **self.pol_ratio_style)

            self.rax.legend(
                bbox_to_anchor = (1.05, 1),
                loc            = "upper left",
                borderaxespad  = 0.,
                fontsize       = "x-small"
            )

        self.rax.axhline(y=1, linestyle="--", color="black")

        self.rax.set_ylim(0, 2)
        self.rax.set_xlabel(self.xlabel)
        self.rax.autoscale(axis='x', tight=True)

    def plotter(self, prop, reg, fs, norm=False):
        # Set up figures, styling, names etc
        self._set_procs(reg)
        self._set_prop_info(prop, reg)   
        self.writeFigs()
        self.cms_label()

        hists  = self.all_hists[prop][reg][fs]
        counts = self.all_counts[prop][reg][fs]
        errors = self.all_errors[prop][reg][fs]

        self.draw_mc   = not self.is_empty(hists, "MC")
        self.draw_pol  = not self.is_empty(hists, "Pol") and (prop in self.ang_vars)
        self.draw_data = not self.is_empty(hists, "Data") and not self.blind

        max_bin_counts = []

        if self.draw_mc:
            mc_hists  = hists["MC"]
            mc_counts = counts["MC"]
            mc_errors = errors["MC"]

            mc_labels = self.addCounts(self.labels["MC"], mc_counts, mc_errors)

            # Bin centers and edges
            self._setBins(mc_hists)

            # Draw MC hists
            self.draw_mc_hists(mc_hists, mc_labels)

            # Draw MC Err
            self.set_stackCountsErrs(mc_hists, mc_errors)

            self.ax.fill_between(x=self.bin_edges[:-1], y1 = self.total_mc_counts - self.total_mc_errors, y2 = self.total_mc_counts + self.total_mc_errors, label = "Stat. Unc.", step="post", **self.hatch_style)

            max_bin_counts.append(max([max(hist[0]) for hist in mc_hists.values()]))

        if self.draw_pol:
            pol_hists  = hists["Pol"]
            pol_counts = counts["Pol"]
            pol_errors = errors["Pol"]

            if norm:
                bins = pol_hists["ZLZL"][1]
                pol_hists = dict(
                    ZLZL = (pol_hists["ZLZL"][0]/pol_hists["ZLZL"][0].sum(), bins),
                    ZLZT = (pol_hists["ZLZT"][0]/pol_hists["ZLZT"][0].sum(), bins),
                    ZTZT = (pol_hists["ZTZT"][0]/pol_hists["ZTZT"][0].sum(), bins),
                )
                pol_labels = self.labels["Pol"]
            else:
                pol_labels = self.addCounts(self.labels["Pol"], pol_counts, pol_errors)

            self._setBins(pol_hists)

            self.draw_pol_hists(pol_hists, pol_labels)
            self.draw_pol_sum(hists, counts, errors)
            
            max_bin_counts.append(max([max(hist[0]) for hist in pol_hists.values()]))

        if self.draw_data:
            data_hist   = hists["Data"]
            data_counts = counts["Data"]

            data_label = f"{self.labels['Data']}: {data_counts} " + r"$\pm$" + str(round(np.sqrt(data_counts), 2))

            self._setBins(data_hist)

            self.draw_data_hist(data_hist, data_label)

            max_bin_counts.append(max(data_hist[0]))


        if not self.draw_pol:
            # Adjust plot size to fit legend
            max_bin_count = max(max_bin_counts)
            self.ax.set_ylim(0, max_bin_count*1.75) 
            self.ax.legend(ncol=2, fontsize="x-small")
        else:
            # Place legend outside of main axis
            self.ax.legend(
                bbox_to_anchor = (1.05, 1),
                loc            = "upper left",
                borderaxespad  = 0.,
                fontsize       = "x-small"
            )

        # Labels and ratio (xlabel drawn on rax if unblinded)
        self.ax.set_ylabel(self.ylabel)

        self.draw_ratio(hists, errors)

        hep.rescale_to_axessize(self.ax, 10, 10/1.62)

        return self.fig

    def main(self):
        figs = {}
        for prop in self.all_hists.keys():
            figs[prop] = {}
            for reg in self.all_hists[prop].keys():
                figs[prop][reg] = {}
                for fs in self.all_hists[prop][reg].keys():
                    self.fs = fs
                    figs[prop][reg][fs] = self.plotter(prop, reg, fs)
        return figs