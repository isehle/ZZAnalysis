import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

from NanoAnalysis.scripts.fileHandler import FileHandler

import ROOT
#ROOT.ROOT.EnableImplicitMT()

import uproot as up
import numpy as np
from tqdm import tqdm

ROOT.gInterpreter.Declare(
    """
using ROOT::RVecF;
using ROOT::RVecI;
float propByRegFS(const RVecF &prop, Short_t reg_idx, const RVecI &z1_flav, const RVecI &z2_flav, int pdg_1, int pdg_2){
    if (reg_idx >= 0) {
        if ((pdg_1 == -1) && (pdg_2 == -1)){
            return prop[reg_idx];
        }
        else if ((z1_flav[reg_idx] == pdg_1) && (z2_flav[reg_idx] == pdg_2)) {
            return prop[reg_idx];
        }
    }
    return -1.;
}
"""
)

ROOT.gInterpreter.Declare(
"""
using ROOT::RVecF;
using ROOT::RVecI;
float lepPropByRegFS(const RVecF &lep_prop, const RVecI &cand_lep_idx, const RVecI &z_flav, Short_t reg_idx, int pdg){
    if (reg_idx >= 0){
        if (pdg==-1){
            return lep_prop[cand_lep_idx[reg_idx]];
        }
        else if (z_flav[reg_idx] == pdg){
            return lep_prop[cand_lep_idx[reg_idx]];
        }
    }
    return -1.;
}
"""
)

ROOT.gInterpreter.Declare(
"""
using ROOT::RVecF;
ROOT::VecOps::RVec<double> weight(RVecF dataMCWeight, Float_t overallEventWeight, double genEventSumw){
    return dataMCWeight*overallEventWeight/genEventSumw;
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
            fs_4l    = (-1, -1),
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
        dfs = {}
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
                data_df = ROOT.RDataFrame("Events", subprocs)
                dfs["Data"] = data_df.Filter("HLT_passZZ4l")

        df = ROOT.RDataFrame(spec)
        ROOT.RDF.Experimental.AddProgressBar(df)
        
        df = df.DefinePerSample("sample_type",'rdfsampleinfo_.GetS("sample_type")')
        df = df.DefinePerSample("genEventSumw", 'rdfsampleinfo_.GetD("genEventSumw")')

        dfs["MC"] = df.Filter("HLT_passZZ4l")

        return dfs

    def _get_reg(self, reg):
        if reg == "SR":
            return "ZZCand", "bestCandIdx"
        else:
            return "ZLLCand", f"ZLLbest{reg}Idx"

    def define_lep_cols(self, df, prop, reg, fs):
        base_lep_prop = "_".join(prop.split("_")[:2]) # We usually want Lepton_sip3d_Z2. The base corresponds to Lepton_sip3d.

        if not base_lep_prop in df.GetColumnNames():
            prop = base_lep_prop.split("_")[-1] # Just sip3d
            df = df.Define(base_lep_prop, f"ROOT::VecOps::Concatenate(Electron_{prop},Muon_{prop})")

        cand, reg_idx = self._get_reg(reg)
        for Z in ["Z1", "Z2"]:
            zflav = f"{cand}_{Z}flav"
            pdg = self.pdgs[fs][0] if Z=="Z1" else self.pdgs[fs][1]
            if "SSSIP" in reg and Z=="Z2": pdg *= -1
            for L in ["l1", "l2"]:
                this_lep     = f"{base_lep_prop}_{Z}{L}_{reg}_{fs}"
                cand_lep_idx = f"{cand}_{Z}{L}Idx"
                
                df = df.Define(this_lep, f"lepPropByRegFS({base_lep_prop}, {cand_lep_idx}, {zflav}, {reg_idx}, {pdg})")
            
            # Combine the two leptons
            lep_by_z     = f"{base_lep_prop}_{Z}_{reg}_{fs}"
            def_lep_by_z = "ROOT::VecOps::RVec<float> {" + f"{base_lep_prop}_{Z}l1_{reg}_{fs}, {base_lep_prop}_{Z}l2_{reg}_{fs}" + "}"
            
            df = df.Define(lep_by_z, def_lep_by_z)
        
        # Combine the four leptons
        df = df.Define(f"{base_lep_prop}_{reg}_{fs}", f"ROOT::VecOps::Concatenate({base_lep_prop}_Z1_{reg}_{fs}, {base_lep_prop}_Z2_{reg}_{fs})")

        return df

    def define_cols(self, dfs):
        new_dfs = {}
        for key, df in dfs.items():
            props = ["dataMCWeight"] + list(self.properties) if key != "Data" else list(self.properties)
            
            for prop in props:
                for reg in self.regions:
                    cand, reg_idx = self._get_reg(reg)
                    for fs in self.final_states:
                        
                        if "Lepton" in prop:
                            df = self.define_lep_cols(df, prop, reg, fs)
                        else:    
                            new_col = f"{prop}_{reg}_{fs}"
                            old_col = f"{cand}_{prop}"
                            z1flav = f"{cand}_Z1flav"
                            z2flav = f"{cand}_Z2flav"
                            pdg1, pdg2 = self.pdgs[fs]
                            
                            if "SSSIP" in reg and "4l" not in fs: pdg2 *= -1
                            
                            df = df.Define(new_col, f"propByRegFS({old_col}, {reg_idx}, {z1flav}, {z2flav}, {pdg1}, {pdg2})")
                        
                        if prop != "dataMCWeight" and key != "Data":
                            weight = f"weight_{reg}_{fs}"
                            if weight not in df.GetDefinedColumnNames():
                                df = df.Define(f"weight_{reg}_{fs}", f"dataMCWeight_{reg}_{fs}*overallEventWeight/genEventSumw")
            
            new_dfs[key] = df

        return new_dfs
    
    def write_hists(self, dfs):
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

                    col = f"{prop}_{reg}_{fs}"
                    weight = f"weight_{reg}_{fs}"

                    for proc_type in tqdm(self.file_paths, desc = "Process Types", position= 3, leave = False):
                        # Lepton Columns are RVecs, need to handle the filtering differently
                        if ("Lepton" in col) and not (prop.endswith("l1") or prop.endswith("l2")):
                            if "Z" in col:
                                l1_col, l2_col = f"{prop}l1_{reg}_{fs}", f"{prop}l2_{reg}_{fs}"
                                if proc_type != "Data":
                                    df = dfs["MC"]
                                    hist = df.Filter(f'sample_type=="{proc_type}"').Filter(f"{l1_col} >= 0.").Filter(f"{l2_col} >= 0.").Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), col, weight)
                                    mc_hist_list.append(hist)
                                else:
                                    df = dfs["Data"]
                                    hist = df.Filter(f"{l1_col} >= 0.").Filter(f"{l2_col} >= 0.").Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), col)
                                    data_hist_list.append(hist)
                            else:
                                z1l1_col, z1l2_col, z2l1_col, z2l2_col = f"{prop}_Z1l1_{reg}_{fs}", f"{prop}_Z1l2_{reg}_{fs}", f"{prop}_Z2l1_{reg}_{fs}", f"{prop}_Z2l2_{reg}_{fs}"
                                if proc_type != "Data":
                                    df = dfs["MC"]
                                    hist = df.Filter(f'sample_type=="{proc_type}"').Filter(f"{z1l1_col} >= 0.").Filter(f"{z1l2_col} >= 0.").Filter(f"{z2l1_col} >= 0.").Filter(f"{z2l2_col} >= 0.").Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), col, weight)
                                    mc_hist_list.append(hist)
                                else:
                                    df = dfs["Data"]
                                    hist = df.Filter(f"{z1l1_col} >= 0.").Filter(f"{z1l2_col} >= 0.").Filter(f"{z2l1_col} >= 0.").Filter(f"{z2l2_col} >= 0.").Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), col)
                                    data_hist_list.append(hist)
                        else:
                            if proc_type != "Data":
                                df = dfs["MC"]
                                hist = df.Filter(f'sample_type=="{proc_type}"').Filter(f"{col} >= 0").Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), col, weight)
                                mc_hist_list.append(hist)
                            else:
                                df = dfs["Data"]
                                hist = df.Filter(f"{col} >= 0").Histo1D((col, col, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), col)
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

        dfs = self._get_full_df(self.file_paths)

        dfs = self.define_cols(dfs)

        hists = self.write_hists(dfs)

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