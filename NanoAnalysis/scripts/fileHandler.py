import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import yaml

default_cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/files_cfg.yaml")

class FileHandler:
    def __init__(self, year, era, cfg_path=default_cfg_path):
        self.year = year
        self.era  = era
        
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


