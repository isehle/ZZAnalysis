import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

from pathlib import Path

import ROOT

from tqdm import tqdm
import uproot as up
import numpy as np

from NanoAnalysis.scripts.fileHandler import FileHandler
from NanoAnalysis.scripts.histWriter import HistWriter
from NanoAnalysis.scripts.histPlotter import HistPlotter
from NanoAnalysis.scripts.ZpX_estimation import ZpX

import matplotlib.pyplot as plt

from copy import deepcopy

import json

class HistManager:
    def __init__(self, cfg, args):
        self.cfg          = cfg

        self.year         = args["year"] if args["year"] != -1 else "Full"
        self.era          = args["era"]

        self.tag          = args["tag"]

        self.file_handler = FileHandler(self.year, self.era, self.tag)
        self.lumi         = self.file_handler.lumi

        self.get_cfg_path = lambda step: os.path.join(parent_dir, "NanoAnalysis/scripts", self.cfg[step]["cfg_file"])

        self.zpx = ZpX()

        self.hists = {}

    def _get_cfg(self, step):
        cfg_path = self.get_cfg_path(step)
        with open(cfg_path) as config:
            step_cfg = yaml.safe_load(config)
        return step_cfg, self.cfg[step]

    def _get_hist_structure(self, Hists):
        """Helper function to determine directory
        structure of Hists (from uproot). Returns
        list of distinct properties, regions, and final
        states."""
        
        props   = []
        regions = []
        fstates = []

        for key in Hists.keys():
            full_key = key.split("/")
            # Write dict structure
            if len(full_key) == 3:
                prop, reg, fs = full_key
                fs = fs.replace(";1", '')
                if prop not in props: props.append(prop)
                if reg not in regions:regions.append(reg)
                if fs not in fstates: fstates.append(fs)

        return props, regions, fstates

    def write_hists(self):
        """Wrapper for histogram writing."""

        write_cfg, step_cfg = self._get_cfg("writing")

        regions = step_cfg["regions"]
        fstates = step_cfg["fstates"]
        col_tag = step_cfg["col_tag"]

        props   = write_cfg.keys()

        hist_writer = HistWriter(write_cfg, self.lumi, col_tag)

        hists = hist_writer.main(self.file_handler.file_paths, regions, props, fstates)
        self.file_handler.write_hists(hists)

        return hists

    def rebin_hist(self, hist, new_bins, overflow, name=""):
        """Rebins histogram, saves new histogram to numpy hist.
        Returns numpy hist, total count, and rebinned bin errors."""

        pyhist = hist.to_pyroot()

        nbins = len(new_bins) - 1 
        
        rebinned = pyhist.Rebin(nbins, name, np.array(new_bins))

        max_bin = nbins + 2 if overflow else nbins+1

        bin_contents = np.array([rebinned.GetBinContent(i) for i in range(1, max_bin)])
        bin_errors   = np.array([rebinned.GetBinError(i) for i in range(1, max_bin)])
        
        new_bins = np.array(new_bins)
        if overflow:
            new_bins = np.append(new_bins, np.inf)

        count = bin_contents.sum()

        return (bin_contents, new_bins), count, bin_errors

    def combine_procs(self, hists, counts, errors, **kwargs):
        """Combines processes into groups, for example
        ZpX = DY+TT+WZ."""
        
        for prop in hists.keys():
            for reg in hists[prop].keys():
                for fs in hists[prop][reg].keys():
                    mc_hists  = hists[prop][reg][fs]["MC"]
                    mc_counts = counts[prop][reg][fs]["MC"]
                    mc_errors = errors[prop][reg][fs]["MC"]

                    for group, procs in kwargs.items():
                        group_count = np.array([mc_counts[proc] for proc in procs]).sum()

                        group_hist, edges = mc_hists[procs[0]]

                        group_error = np.square(mc_errors[procs[0]])

                        for proc in procs[1:]:
                            group_hist += mc_hists[proc][0]
                            
                            group_error += np.square(mc_errors[proc])

                        group_hist = (group_hist, edges)
                        group_error = np.sqrt(group_error)

                        hists[prop][reg][fs]["MC"][group] = group_hist
                        counts[prop][reg][fs]["MC"][group] = group_count
                        errors[prop][reg][fs]["MC"][group] = group_error

        return hists, counts, errors
                    

    def read_hists(self, path, **kwargs):
        """Stores histograms, counts, and bin errors
        in dictionaries with the structure:
        <property>/<region>/<final_state>/<MC/Data/Pol>/<process>.
        Also handles optional rebinning and histogram grouping."""
        hists, counts, errors = {}, {}, {}
        with up.open(path) as Hists:
            props, regions, fstates = self._get_hist_structure(Hists)
            for prop in props:
                hists[prop] = {}
                counts[prop] = {}
                errors[prop] = {}
                for reg in regions:
                    hists[prop][reg] = {}
                    counts[prop][reg] = {}
                    errors[prop][reg] = {}
                    for fs in fstates:
                        hists[prop][reg][fs] = dict(
                            MC   = {},
                            Data = {},
                            Pol  = {},
                        )
                        counts[prop][reg][fs] = dict(
                            MC   = {},
                            Data = {},
                            Pol  = {},
                        )
                        errors[prop][reg][fs] = dict(
                            MC   = {},
                            Data = {},
                            Pol  = {},
                        )
                        for proc, hist in Hists[prop][reg][fs].items():
                            proc = proc.replace(";1","")

                            if len(kwargs) > 0:
                                if kwargs["rebin"] and prop in kwargs:
                                    if prop == "mass":
                                        if reg == "SR" or "HighMass" in reg:
                                            bin_info = kwargs[prop]["HighMass"]
                                        elif "MidMass" in reg:
                                            bin_info = kwargs[prop]["MidMass"]
                                        elif "LowMass" in reg:
                                            bin_info = kwargs[prop]["LowMass"]
                                    else:
                                        bin_info = kwargs[prop]

                                    this_hist, count, err = self.rebin_hist(hist, bin_info["bins"], bin_info["overflow"])

                                else:
                                    this_hist = hist.to_numpy(flow=True)
                                    this_hist = (this_hist[0][1:], this_hist[1][1:])
                                    count     = this_hist[0].sum()
                                    err       = hist.errors(flow=True)[1:]
                            else:
                                this_hist = hist.to_numpy(flow=True)
                                this_hist = (this_hist[0][1:], this_hist[1][1:])
                                count     = this_hist[0].sum()
                                err       = hist.errors(flow=True)[1:]

                            if proc in list(self.file_handler.mc_procs.keys()) + ["ZpX"]:

                                hists[prop][reg][fs]["MC"][proc]  = this_hist
                                counts[prop][reg][fs]["MC"][proc] = count
                                errors[prop][reg][fs]["MC"][proc] = err

                            elif proc == "Data":
                                hists[prop][reg][fs]["Data"]  = this_hist
                                counts[prop][reg][fs]["Data"] = count
                                errors[prop][reg][fs]["Data"] = err

                            elif proc in self.file_handler.pol_procs:
                                hists[prop][reg][fs]["Pol"][proc]  = this_hist
                                counts[prop][reg][fs]["Pol"][proc] = count
                                errors[prop][reg][fs]["Pol"][proc] = err

        if "group" in kwargs and kwargs["group"]:
            hists, counts, errors = self.combine_procs(hists, counts, errors, **kwargs["groups"])
            
        return hists, counts, errors
    
    def plot_hists(self, path="", **kwargs):
        """Plots histograms stored in <path>. Histograms will first
        be written to <path> if it does not yet exist."""
        if path=="": path = self.file_handler.hist_path
        if not os.path.exists(path):
            self.write_hists()

        hists, counts, errors = self.read_hists(path, **kwargs)

        zpx_from_mc = "ZpX" in kwargs["groups"] and kwargs["group"]

        hist_plotter = HistPlotter(self.year, self.era, self.tag, self.lumi, hists, counts, errors, zpx_from_mc)
        figs = hist_plotter.main()

        self.file_handler.write_plots(figs)

    def combine_eras(self):
        """Combines histograms from two different ROOT files, outputting new
        new ROOT file in appropriate directory. Input file paths are in gen_cfg.yaml."""
        
        if self.year != -1:
            year_tag = f"year_{str(self.year)}"        
        else:
            year_tag = f"year_22_23"

        infiles = list(self.cfg["combine_eras"][year_tag].values())
        eras    = list(self.cfg["combine_eras"][year_tag].keys())

        base_dir, filename = os.path.split(infiles[0])
        
        if self.year != -1:
            # .../2022/CD/... --> .../2022/Full/...
            outdir = base_dir.replace(eras[0], "Full")
        else:
            # .../2022/Full/... --> .../Full/...
            base_dir.replace("/2022", "")
        
        Path(outdir).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(outdir, filename)

        with up.open(infiles[0]) as Hists_1, up.open(infiles[1]) as Hists_2, up.recreate(outfile) as NewHists:
            for key in tqdm(Hists_1.keys()):
                if key.count("/")==3 and key in Hists_2.keys():
                    hist_1, hist_2 = Hists_1[key], Hists_2[key]
                    NewHists[key.replace(";1", "")] = hist_1.to_pyroot() + hist_2.to_pyroot()

    def plot_zpx(self, recalc):
        """Plot Z+X info for 2022, 2023 or 2022+2023.
        By default will use saved info from json files (paths
        set in gen_cfg.yaml), but with recalc==True will recalculate
        the estimations from the histograms."""

        if not recalc:
            if self.year != "Full":

                year_tag = f"year_{str(self.year)}"
                with open(self.cfg["zpx"][year_tag], "r") as InFile:
                    zpx_info = json.load(InFile)
                
                for step in zpx_info.keys():
                    self.zpx.plot_zpx(zpx_info, self.year, step, tag=self.tag)
            
            else:

                with open(self.cfg["zpx"]["year_2022"], "r") as InFile:
                    zpx_info_22 = json.load(InFile)
                with open(self.cfg["zpx"]["year_2023"], "r") as InFile:
                    zpx_info_23 = json.load(InFile)

                for step in zpx_info_22.keys():
                    self.zpx.plot_zpx_years(zpx_info_22, zpx_info_23, step, tag=self.tag)

        else:
            assert self.year != -1, "Can only calculate ZpX for one year at a time."

            hist_path = self.file_handler.hist_path
            
            if not os.path.exists(hist_path):
                self.write_hists()
            
            hists, counts, errors = self.read_hists(hist_path)

            zpx_info = self.zpx.get_zpx(hists, errors)
            
            outpath = f"ZpX_info_{self.year}_{self.era}_{self.tag}.json"
            with open(outpath, "w") as outfile:
                json.dump(zpx_info, outfile, indent=4)
            
            for step in zpx_info.keys():
                self.zpx.plot_zpx(zpx_info, self.year, step, tag=self.tag)

if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--year", choices=(2022, 2023, -1), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "Full"), default="EFG")
    parser.add_argument("--tag", default="")
    parser.add_argument("--mode", choices=("plot_hists", "plot_zpx", "combine_eras"), default="plot_hists")
    parser.add_argument("--recalculate_zpx", action='store_true')
    args = vars(parser.parse_args())    

    cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/gen_cfg.yaml")
    with open(cfg_path) as config:
        cfg = yaml.safe_load(config)
    
    hist_manager = HistManager(cfg, args)

    if args["mode"] == "plot_hists":
        plot_cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/plot_cfg.yaml")
        with open(plot_cfg_path) as plot_config:
            plot_cfg = yaml.safe_load(plot_config)

        hist_manager.plot_hists(**plot_cfg["extra"])
    
    elif args["mode"] == "plot_zpx":
        hist_manager.plot_zpx(recalc=args["recalculate_zpx"])

    elif args["mode"] == "combine_eras":
        hist_manager.combine_eras()