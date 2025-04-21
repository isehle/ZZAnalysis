import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT

from tqdm import tqdm
import uproot as up

from NanoAnalysis.scripts.fileHandler import FileHandler
from NanoAnalysis.scripts.histWriter import HistWriter

def get_df(path):
    df = ROOT.RDataFrame("Events", path)
    df = df.Filter("HLT_passZZ4l")

    if "Data" in path:
        return df

    Runs = ROOT.RDataFrame("Runs", path)
    return df.Define("genEventSumw", str(Runs.Sum("genEventSumw").GetValue())) 

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
        return step_cfg, self.cfg[step]

    def write_hists(self):
        write_cfg, step_cfg = self._get_cfg("writing")

        regions = step_cfg["regions"]
        fstates = step_cfg["fstates"]
        col_tag = step_cfg["col_tag"]

        props   = write_cfg.keys()

        histWriter = HistWriter(write_cfg, self.lumi, col_tag)

        file_paths = self.file_handler.file_paths

        hists = histWriter.main(file_paths, regions, props, fstates)
        #hists = histWriter.main(self.file_handler, regions, props, fstates)

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
    hists = hist_manager.write_hists()

    print("\nWriting...")
    with up.recreate("hist_output_test_v5.root") as OutFile:
        for prop in tqdm(hists.keys(), desc = "Properties", position = 0):
            for reg in tqdm(hists[prop].keys(), desc = "Regions", position = 1, leave = False):
                for fs in tqdm(hists[prop][reg].keys(), desc = "Final States", position = 2, leave = False):
                    for proc_type in tqdm(hists[prop][reg][fs].keys(), desc = "Processes", position = 3, leave = False):
                        hist = hists[prop][reg][fs][proc_type].GetValue() if "Data" not in proc_type else hists[prop][reg][fs][proc_type]
                        OutFile[f"{prop}/{reg}/{fs}/{proc_type}"] = hist