import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import datetime
from pathlib import Path

import yaml

import ROOT
import uproot as up
from tqdm import tqdm

import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep

plot_cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/plot_cfg.yaml")

class HistPlotter:
    def __init__(self, year, era, tag, lumi, all_hists, all_counts, all_errors, cfg_path=plot_cfg_path):
        self.year = year
        self.era  = era
        self.tag  = tag

        self.lumi = lumi

        self.all_hists  = all_hists
        self.all_counts = all_counts
        self.all_errors = all_errors

        self._set_cfg(cfg_path)

        self.pol_colors      = self.cfg["pol_colors"]

        self.ratio_fig_style = self.cfg["ratio_fig_style"]
        self.hatch_style     = self.cfg["hatch_style"]
        self.errorbar_style  = self.cfg["errorbar_style"]

        self.fill_colors = [self.cfg["proc_info"][proc]["fill"] for proc in self.cfg["proc_info"].keys()]
        self.line_colors = [self.cfg["proc_info"][proc]["line"] for proc in self.cfg["proc_info"].keys()]

        self.bin_edges   = []
        self.bin_centers = []

        self.is_empty = lambda a_dict, key: a_dict[key] == {}

        self.labels = dict(
            MC = dict(
                ggZZ   = r"$gg \rightarrow ZZ$",
                ZZ_NLO = r"$(q\bar{q} \rightarrow ZZ)_{NNLO}$",
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
            #mass_reg = reg.split("_")[-1] if reg != "SR" or "HighMass" in reg else reg.replace("OSSIP","").replace("SSSIP", "")
            self.prop_info = self.prop_info[mass_reg]

        self.xlabel = self._to_raw_string(self.prop_info["xlabel"])
        self.ylabel = self._to_raw_string(self.prop_info["ylabel"])
        
        self.blind  = self.prop_info["Blind"]

    def _setBins(self, hists):
        if isinstance(hists, dict):
            key = list(hists.keys())[0]
            self.bin_edges = hists[key][1]
        else:
            self.bin_edges = hists[1]
        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:])/2
        # if self.bin_edges == []:
        #     key = list(hists.keys())[0]
        #     self.bin_edges = hists[key][1]
        #     self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:])/2

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
            bin_errs = errors[proc]
            err      = np.sqrt(np.square(bin_errs).sum())
            
            new_lab = init_lab + ": " + str(round(count,2)) + r" $\pm$ " + str(round(err, 2))
            new_labels.append(new_lab)
        return new_labels

    def writeFigs(self):
        hep.style.use("CMS")

        if not self.blind:
            self.fig, (self.ax, self.rax) = plt.subplots(2, 1, sharex=True, **self.ratio_fig_style)
            self.fig.subplots_adjust(hspace=0.07)
        else:
            self.fig, self.ax = plt.subplots()

    def cms_label(self):
        hep.cms.label(label="Work in Progress", year=self.year, lumi = round(self.lumi*1e-3), com=13.6, data=True, ax=self.ax)

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

    def draw_data_hist(self, hist, label):
        self.ax.errorbar(self.bin_centers,
                    hist[0],
                    yerr = np.sqrt(hist[0]),
                    color = "black",
                    fmt = "o",
                    label = label,
                    markersize = 3)

    def set_stackCountsErrs(self, mc_hists, mc_errors):
        #mc_count_arr = np.array([mc_hists[key][0] for key in self.cfg["plot_styling"]["mc_colors"].keys()])
        mc_count_arr = np.array([mc_hists[key][0] for key in mc_hists.keys()])
        self.total_mc_counts = np.sum(mc_count_arr, axis=0)
        self.total_mc_counts[ self.total_mc_counts < 0] = 0
        self.total_mc_errors = np.sqrt(np.sum([err**2 for err in mc_errors.values()], axis=0))

    def draw_ratio(self, data_hist):
        self.rax.fill_between(x=self.bin_centers, y1= 1 - self.total_mc_errors/self.total_mc_counts, y2 = 1 + self.total_mc_errors/self.total_mc_counts, step='mid', **self.hatch_style)
        self.rax.errorbar(x=self.bin_centers, y=data_hist[0]/self.total_mc_counts, yerr=np.sqrt(data_hist[0])/self.total_mc_counts, **self.errorbar_style)

        self.rax.set_ylim(0, 2)
        self.rax.set_ylabel('Data / MC')
        self.rax.set_xlabel(self.xlabel)
        self.rax.autoscale(axis='x', tight=True)

    def plotter(self, prop, reg, fs, norm=False):
        # Set up figures, styling, names etc
        self._set_prop_info(prop, reg)
        #if not norm: self.set_lumi_tag()    
        self.writeFigs()
        self.cms_label()

        hists  = self.all_hists[prop][reg][fs]
        counts = self.all_counts[prop][reg][fs]
        errors = self.all_errors[prop][reg][fs]

        self.draw_mc   = not self.is_empty(hists, "MC")
        self.draw_pol  = not self.is_empty(hists, "Pol")
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
            try:
                self.ax.fill_between(x=self.bin_centers[1:-1], y1 = self.total_mc_counts[1:-1] - self.total_mc_errors[1:-1], y2 = self.total_mc_counts[1:-1] + self.total_mc_errors[1:-1], label = "Stat. Unc.", step='mid', **self.hatch_style)
            except:
                breakpoint()

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

            #data_label = self.addCounts(self.labels["Data"], data_counts, data_errors)
            #data_label = self.addCounts(self.labels, data_counts, data_errors)

            self._setBins(data_hist)

            self.draw_data_hist(data_hist, data_label)

            max_bin_counts.append(max(data_hist[0]))

        # Adjust plot size to fit legend
        max_bin_count = max(max_bin_counts)
        self.ax.set_ylim(0, max_bin_count*2.5) 

        self.ax.legend(ncol=2, fontsize="x-small")

        # Labels and ratio (xlabel drawn on rax if unblinded)
        self.ax.set_ylabel(self.ylabel)
        if self.draw_data:
            self.draw_ratio(data_hist)
        else:
            self.ax.set_xlabel(self.xlabel)

        hep.rescale_to_axessize(self.ax, 10, 10/1.62)

        return self.fig

    def main(self):
        figs = {}
        for prop in self.all_hists.keys():
            figs[prop] = {}
            for reg in self.all_hists[prop].keys():
                figs[prop][reg] = {}
                for fs in self.all_hists[prop][reg].keys():
                    #if fs != "fs_4l": continue
                    figs[prop][reg][fs] = self.plotter(prop, reg, fs)
        return figs