import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import multiprocessing as mp
# mp.set_start_method('spawn')

from multiprocessing.dummy import Pool

from copy import deepcopy

from NanoAnalysis.scripts.fileHandler import FileHandler
from NanoAnalysis.scripts.histWriter import HistWriter

class HistManager:
    def __init__(self, cfg, args):
        self.cfg          = cfg

        self.year         = args["year"]
        self.era          = args["era"]

        self.file_handler = FileHandler(self.year, self.era)
        self.lumi         = self.file_handler.lumi

        self.get_cfg_path = lambda step: os.path.join(parent_dir, "NanoAnalysis/scripts", self.cfg[step]["cfg_file"])

        self.hists = {}

    def _get_cfg(self, step):
        cfg_path = self.get_cfg_path(step)
        with open(cfg_path) as config:
            step_cfg = yaml.safe_load(config)
        return step_cfg

    def write_hists(self):
        step_cfg = self._get_cfg("writing")
        histWriter = HistWriter(step_cfg)

        regions = self.cfg["writing"]["regions"]
        fstates = self.cfg["writing"]["fstates"]
        props   = step_cfg.keys()

        hists = {}

        df = histWriter.get_df("ZZ4l_NLO_HZZSelection_Skim.root")

        for reg in regions:
            hists[reg] = {}
            reg_df = histWriter.init_defs(df, reg)
            for prop in props:
                hists[reg][prop] = {}
                reg_df = histWriter.define_cols(reg_df, prop)
                for fs in fstates:
                    fs_df = histWriter.fs_filt(reg_df, fs)
                    hists[reg][prop][fs] = histWriter.write_hist(fs_df, prop, reg, self.lumi)

        # Attempt to parallelize —— since the histWriter's class methods are all trying to access the same object at once,
        # we're currently getting a segfault. Attempt to use deepcopy(df) did not solve the issue.
        # with mp.Pool(processes=len(regions)) as pool:
        #     reg_dfs = pool.starmap(histWriter.init_defs, [(df, reg) for reg in regions])

        # Loop over all the processed ntuples (CJLST output) (commented because only have local file atm)
        # for proc, subprocs in self.file_handler.file_paths.items():
        #     hists[proc] = {}
        #     if isinstance(subprocs, dict):
        #         for path in subprocs.values():
        #             df = histWriter.get_df(path)
        #             breakpoint()

        return hists

    

if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--year", choices=(2022, 2023), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "Full"), default="EFG")
    parser.add_argument("--tag", default="")
    args = vars(parser.parse_args())    

    cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/gen_cfg.yaml")
    with open(cfg_path) as config:
        cfg = yaml.safe_load(config)
    
    hist_manager = HistManager(cfg, args)
    hist_manager.write_hists()