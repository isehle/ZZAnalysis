import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

from NanoAnalysis.scripts.fileHandler import FileHandler

import ROOT
import uproot as up
import numpy as np
from tqdm import tqdm

class HistWriter:
    def __init__(self, cfg):
        self.cfg          = cfg
    
        self.col_tag = "good"

        self.lep_idx = lambda z, l: f"{self.col_tag}_{z}{l}Idx"
        self.lep_col = lambda lep_prop, z, l: f"{self.col_tag}_{lep_prop}_{z}{l}"

        self.pdgs = dict(
            fs_4e    = (-121, -121),
            fs_4mu   = (-169, -169),
            fs_2e2mu = (-121, -169),
            fs_2mu2e = (-169, -121)
        )

        self.isData = False

    def get_hist_info(self, reg, prop):
        hist_info = self.cfg[prop]
        if "mass" in prop:
            if reg == "SR" or "HighMass" in reg:
                hist_info = hist_info["SR"]
            else:
                reg = reg.split("_")
                hist_info = hist_info[reg[-1]]

        return hist_info

    def get_df(self, path):
        df = ROOT.RDataFrame("Events", path)
        df = df.Filter("HLT_passZZ4l")

        if "Data" in path:
            self.isData = True
            return df

        Runs = ROOT.RDataFrame("Runs", path)
        return df.Define("genEventSumw", str(Runs.Sum("genEventSumw").GetValue()))    

    def init_defs(self, df, reg):
        """Defines columns that aren't necessarily plotted
        but are useful later, i.e. which candidate idx to use
        and lepton indices."""
        if reg == "SR":
            self.cand = "ZZCand"
            self.reg_idx  = "bestCandIdx"

        else:
            self.cand = "ZLLCand"
            self.reg_idx  = f"ZLLbest{reg}Idx"

        if not self.isData:
            df = df.Define("weight", f"{self.cand}_dataMCWeight[{self.reg_idx}]*overallEventWeight/genEventSumw")

        for Z in ["Z1", "Z2"]:
            for L in ["l1", "l2"]:
                df = df.Define(self.lep_idx(Z, L), f"{self.cand}_{Z}{L}Idx[{self.reg_idx}]")

        return df

    def define_cols(self, df, prop):
        if "Lepton" not in prop:
            col = f"{self.col_tag}_{prop}"
            return df.Define(col, f"{self.cand}_{prop}[{self.reg_idx}]")
        else:
            for Z in ["Z1", "Z2"]:
                for L in ["l1", "l2"]:
                    lep_idx = self.lep_idx(Z, L)
                    lep_col = self.lep_col(prop, Z, L)
                    df = df.Define(lep_col, f"{prop}[{lep_idx}]")

        return df

    def fs_filt(self, df, fs):
        if "4l" in fs:
            return df
        
        z1flav, z2flav = self.pdgs[fs]

        return df.Filter(f"{self.cand}_Z1flav[{self.reg_idx}] == {z1flav}").Filter(f"{self.cand}_Z2flav[{self.reg_idx}] == {z2flav}")

    def write_hist(self, df, col, reg, lumi, weight_col="weight"):
        hist_info = self.get_hist_info(reg, col)
        if not self.isData:
            hist = df.Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), self.col_tag + "_" + col, weight_col)
            hist.Scale(lumi)
            hist = hist.GetValue()
        else:
            hist = df.Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), self.col_tag + "_" + col, weight_col)
            hist.SetBinErrorOption(ROOT.TH1.kPoisson)
        return hist
    
if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--year", choices=(2022, 2023), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "Full"), default="EFG")
    parser.add_argument("--tag", default="")
    args = vars(parser.parse_args())

    cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/hist_cfg.yaml")
    with open(cfg_path) as config:
        cfg = yaml.safe_load(config)

    histWriter = HistWriter(cfg)
    df = histWriter.main("ZZ4l_NLO_HZZSelection_Skim.root")