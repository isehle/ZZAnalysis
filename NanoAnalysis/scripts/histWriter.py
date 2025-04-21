import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

from NanoAnalysis.scripts.fileHandler import FileHandler

import ROOT
ROOT.ROOT.EnableImplicitMT()

import uproot as up
import numpy as np
from tqdm import tqdm

ROOT.gInterpreter.Declare(
    """
using ROOT::RVecF;
using ROOT::RVecI;
float propByRegFS(const RVecF &prop, Short_t reg_idx, const RVecI &z1_flav, const RVecI &z2_flav, int pdg_1, int pdg_2){
    if (reg_idx >= 0) {
        if ((z1_flav[reg_idx] == pdg_1) && (z2_flav[reg_idx] == pdg_2)) {
            return prop[reg_idx];
        }
    }
    return -1.;
}
"""
)

class HistWriter:
    def __init__(self, cfg, lumi, col_tag):
        self.cfg  = cfg
        self.lumi = lumi

        self.col_tag = col_tag

        self.lep_idx = lambda z, l, reg: f"{self.col_tag}_{z}{l}Idx_{reg}"
        self.lep_col = lambda lep_prop, z, l, reg: f"{self.col_tag}_{lep_prop}_{z}{l}_{reg}"

        self.pdgs = dict(
            fs_4e    = (-121, -121),
            fs_4mu   = (-169, -169),
            fs_2e2mu = (-121, -169),
            fs_2mu2e = (-169, -121)
        )

        self.isData = False

        self.mc_hists   = {}
        self.pol_hists  = {}
        self.data_hists = {}

    def get_hist_info(self, reg, prop):
        hist_info = self.cfg[prop]
        if "mass" in prop:
            if reg == "SR" or "HighMass" in reg:
                hist_info = hist_info["SR"]
            elif "MidMass" in reg:
                hist_info = hist_info["MidMass"]
            elif "LowMass" in reg:
                hist_info = hist_info["LowMass"]

        return hist_info
    
    def _get_full_df(self, file_paths):
        spec = ROOT.RDF.Experimental.RDatasetSpec()

        for proc, subprocs in file_paths.items():
            if isinstance(subprocs, dict):
                for subproc, subproc_path in subprocs.items():
                    Runs = ROOT.RDataFrame("Runs", subproc_path)
                    genEventSumw = Runs.Sum("genEventSumw").GetValue()

                    meta = ROOT.RDF.Experimental.RMetaData()
                    meta.Add("name", subproc)
                    meta.Add("sample_type", proc)
                    meta.Add("genEventSumw", float(genEventSumw))
                    
                    sample = ROOT.RDF.Experimental.RSample(subproc, "Events", subproc_path, meta)
                    spec.AddSample(sample)
            
            elif "Data" not in proc:
                Runs = ROOT.RDataFrame("Runs", subprocs)
                genEventSumw = Runs.Sum("genEventSumw").GetValue()

                meta = ROOT.RDF.Experimental.RMetaData()
                meta.Add("name", proc)
                meta.Add("sample_type", proc)
                meta.Add("genEventSumw", float(genEventSumw))

                sample = ROOT.RDF.Experimental.RSample(proc, "Events", subprocs, meta)
                spec.AddSample(sample)
            else:
                meta = ROOT.RDF.Experimental.RMetaData()
                meta.Add("name", proc)
                meta.Add("sample_type", "Data")
                meta.Add("genEventSumw", float(1.0))

                sample = ROOT.RDF.Experimental.RSample(proc, "Events", subprocs, meta)
                spec.AddSample(sample)

        df = ROOT.RDataFrame(spec)
        ROOT.RDF.Experimental.AddProgressBar(df)
        
        df = df.DefinePerSample("sample_type",'rdfsampleinfo_.GetS("sample_type")')
        df = df.DefinePerSample("genEventSumw", 'rdfsampleinfo_.GetD("genEventSumw")')

        return df

    def _get_reg(self, reg):
        if reg == "SR":
            return "ZZCand", "bestCandIdx"
        else:
            return "ZLLCand", f"ZLLbest{reg}Idx"

    def define_cols(self, df):
        for prop in ["dataMCWeight"] + list(self.properties):
            for reg in self.regions:
                cand, reg_idx = self._get_reg(reg)
                for fs in self.final_states:
                    if fs == "fs_4l": continue
                    new_col = f"{prop}_{reg}_{fs}"
                    old_col = f"{cand}_{prop}"
                    z1flav = f"{cand}_Z1flav"
                    z2flav = f"{cand}_Z2flav"
                    pdg1, pdg2 = self.pdgs[fs]
                    if "SSSIP" in reg: pdg2 *= -1
                    df = df.Define(new_col, f"propByRegFS({old_col}, {reg_idx}, {z1flav}, {z2flav}, {pdg1}, {pdg2})")
                    if prop != "dataMCWeight":
                        weight = f"weight_{reg}_{fs}"
                        if weight not in df.GetDefinedColumnNames():
                            df = df.Define(f"weight_{reg}_{fs}", f"dataMCWeight_{reg}_{fs}*overallEventWeight/genEventSumw")

        return df
    
    def write_hists(self, df):
        mc_hist_list   = []
        data_hist_list = []
        
        # Book hists
        hists = {}
        for prop in tqdm(self.properties, desc = "Properties", position= 0):
            hists[prop] = {}
            for reg in tqdm(self.regions, desc = "Regions", position= 1, leave = False):
                hist_info = self.get_hist_info(reg, prop)
                hists[prop][reg] = {}
                for fs in tqdm(self.final_states, desc = "Final States", position= 2, leave = False):
                    hists[prop][reg][fs] = {}
                    if fs == "fs_4l": continue
                    col = f"{prop}_{reg}_{fs}"
                    weight = f"weight_{reg}_{fs}"

                    for proc_type in tqdm(self.file_paths, desc = "Process Types", position= 3, leave = False):
                        if proc_type != "Data":
                            hist = df.Filter(f'sample_type=="{proc_type}"').Filter(f"{col} >= 0").Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), col, weight)
                            mc_hist_list.append(hist)
                        else:
                            hist = df.Filter(f'sample_type=="{proc_type}"').Filter(f"{col} >= 0").Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), col)
                            data_hist_list.append(hist)

                        hists[prop][reg][fs][proc_type] = hist
        
        # Trigger event loop
        ROOT.RDF.RunGraphs(mc_hist_list)
        ROOT.RDF.RunGraphs(data_hist_list)

        for hist in mc_hist_list:
            hist.Scale(self.lumi)
        for hist in data_hist_list:
            hist.SetBinErrorOption(ROOT.TH1.kPoisson)
        
        return hists

    def main(self, file_paths, regions, properties, final_states):
        self.file_paths   = file_paths

        self.regions      = regions
        self.properties   = properties
        self.final_states = final_states

        df = self._get_full_df(self.file_paths)

        df = self.define_cols(df)

        hists = self.write_hists(df)

        return hists
    
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