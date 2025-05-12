import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

from pathlib import Path

import yaml
import uproot as up

default_cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/files_cfg.yaml")

class FileHandler:
    def __init__(self, year, era, tag, cfg_path=default_cfg_path):
        self.year = year
        self.era  = era

        self.tag  = tag
        
        self._set_cfg(cfg_path)

        self.central_base = self.cfg["eos_base"]
        self.store        = self.cfg["store"]
        self.mc_file_name = self.cfg["mc_file_name"]
        self.mc_procs     = self.cfg["MC_Procs"]
        self.era_info     = self.cfg["year_"+str(year)][era]

        self.lumi = self.era_info["lumi"]

        self.file_paths = {}

        self.mc_samples   = {}
        self.data_samples = {}
        self.pol_samples  = {}

        self._set_file_paths()

        self.output    = self.cfg["output"]
        self.hist_path = os.path.join(self.output["histograms"], str(self.year), self.era, f"hists{self.tag}.root")


    def _set_cfg(self, cfg_path):
        with open(cfg_path) as config:
            self.cfg = yaml.safe_load(config)

    def _set_file_paths(self):
        
        if "MC" in self.era_info:
            mc_path = self.store + os.path.join(self.central_base, self.era_info["MC"])
            for cat, procs in self.mc_procs.items():
                if isinstance(procs, dict):
                    self.mc_samples[cat] = {key: os.path.join(mc_path, val, self.mc_file_name) for key, val in procs.items()}
                else:
                    path = os.path.join(mc_path, procs, self.mc_file_name)
                    self.mc_samples[cat] = path
        
        if "Data" in self.era_info:
            data_path = self.era_info["Data"]
            path = self.store + os.path.join(self.central_base, data_path)
            self.data_samples = dict(
                Data = path
            )

        if "Pol" in self.era_info:
            self.pol_samples = {cat: self.store + os.path.join(self.central_base, pol_file) for cat, pol_file in self.era_info["Pol"].items()}

        self.file_paths = self.mc_samples | self.data_samples | self.pol_samples

    def write_hists(self, hists):
        print(f"\nWriting to {self.hist_path}... ")
        with up.recreate(self.hist_path) as OutFile:
            for prop in hists.keys():
                for reg in hists[prop].keys():
                    for fs in hists[prop][reg].keys():
                        for proc_type in hists[prop][reg][fs].keys():
                            #hist = hists[prop][reg][fs][proc_type].GetValue() if "Data" not in proc_type else hists[prop][reg][fs][proc_type]
                            OutFile[f"{prop}/{reg}/{fs}/{proc_type}"] = hists[prop][reg][fs][proc_type].GetValue()

    def write_plots(self, figs):
        base_dir = os.path.join(self.output["plots"], str(self.year), self.era)
        for prop in figs.keys():
            for reg in figs[prop].keys():
                for fs, fig in figs[prop][reg].items():
                    outdir = os.path.join(base_dir, reg, fs)
                    #outdir = os.path.join(self.output["plots"], self.year, self.era, prop, reg, fs)
                    Path(outdir).mkdir(parents=True, exist_ok=True)
                    outfile = os.path.join(outdir, f"{prop}{self.tag}.png")
                    print("Writing: ", outfile)
                    fig.savefig(outfile)
                    # try:
                    #     fig.savefig(outfile)
                    # except ValueError:
                    #     breakpoint()
                    #     print("Ugh")


