import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT
import uproot as up
from tqdm import tqdm

from NanoAnalysis.scripts.helper_functions import *

class HistWriter:
    def __init__(self, cfg, args):
        self.cfg     = cfg
        self.args    = args
        self.samples = self._get_samples(args["year"], args["era"])
        self.lumi    = float(cfg["datasets"]["year_"+str(args["year"])][args["era"]]["Lumi"])

        self.outfile = os.path.join("/eos/user/i/iehle/Analysis/rootFiles", str(args["year"]), args["era"], "hists{}.root".format(args["tag"]))

        self.regions = cfg["regions"]
        self.fstates = dict(
            fs_4e    = (-121, -121),
            fs_4mu   = (-169, -169),
            fs_2e2mu = (-121, -169),
            fs_2mu2e = (-169, -121)        
        )

        self.props    = list(cfg["hist_info"].keys())

        self.cand     = lambda reg: "ZZCand" if reg=="SR" else "ZLLCand"
        self.reg_filt = lambda reg: "{}.at(0) == 1".format(reg)
        self.z1_flav  = lambda reg, fs: "{}_Z1flav.at(0)=={}".format(self.cand(reg), self.fstates[fs][0])
        self.z2_flav  = lambda reg, fs: "{}_Z2flav.at(0)=={}".format(self.cand(reg), self.fstates[fs][1])
        self.lep_idx  = lambda reg, z, l: "{}_Z{}l{}Idx.at(0)".format(self.cand(reg), z, l)

    def z_flav(self, reg, fs, z):
        pid = self.fstates[fs][z-1]
        if "SS" in reg and z==2:
            pid *= -1
        return "{}_Z{}flav.at(0)=={}".format(self.cand(reg), z, pid)

    def _get_samples(self, year, era):
        central_base = self.cfg["datasets"]["eos_base"]

        mc_sub_path  = self.cfg["datasets"]["year_"+str(year)][era]["MC"]

        central_mc_procs = self.cfg["datasets"]["MC_Procs"]
        
        central_mc_path = os.path.join(central_base, mc_sub_path)

        central_mc_samples = {}
        for cat, procs in central_mc_procs.items():
            if isinstance(procs, dict):
                for key, val in procs.items():
                    central_mc_samples[key] = os.path.join(central_mc_path, val, "ZZ4lAnalysis.root")
            else:
                central_mc_samples[cat] = os.path.join(central_mc_path, procs, "ZZ4lAnalysis.root")

        if "Data" in self.cfg["datasets"]["year_"+str(year)][era]:
            data_sub_path = self.cfg["datasets"]["year_"+str(year)][era]["Data"]
            data_sample = dict(Data = os.path.join(central_base, data_sub_path))
        else:
            data_sample = {}

        if "Pol" in self.cfg["datasets"]["year_"+str(year)][era]:
            pol_samples = {cat: os.path.join(central_base, pol_file) for cat, pol_file in self.cfg["datasets"]["year_"+str(year)][era]["Pol"].items()}
        else:
            pol_samples = {}

        return central_mc_samples | data_sample | pol_samples

    def get_df(self, path):

        df = ROOT.RDataFrame("Events", path)
        df = df.Filter("HLT_passZZ4l")

        if self.isData:
            return df

        Runs = ROOT.RDataFrame("Runs", path)

        return df.Define("genEventSumw",str(Runs.Sum("genEventSumw").GetValue()))

    def write_weight(self, df, reg):
        # dataMCWeight stored as RVec but only ever has one entry
        cand_weight = "{}_dataMCWeight.at(0)*".format(self.cand(reg))

        df = df.Define("weight", "{}overallEventWeight/genEventSumw".format(cand_weight))

        return df.Filter(self.reg_filt(reg))

    def fs_df(self, df, fs, reg):
        return df.Filter(self.z_flav(reg, fs, 1)).Filter(self.z_flav(reg, fs, 2))
        #return df.Filter(self.z1_flav(reg, fs)).Filter(self.z2_flav(reg, fs))

    def lep_df(self, df, reg, lep_prop):
        prop = lep_prop.split("_")[1]

        lep_col = "Lepton_{}".format(prop)
        df = df.Define(lep_col, "ROOT::VecOps::Concatenate(Electron_{},Muon_{})".format(prop, prop))

        for Z in [1, 2]:
            for L in [1, 2]:
                col = lep_col + "_Z{}l{}".format(Z, L)
                df = df.Define(col, "{}.at({})".format(lep_col, self.lep_idx(reg, Z, L)))
                
        def_by_z = lambda z: "ROOT::VecOps::RVec<float> {" + lep_col +"_Z{}l1, ".format(z) + lep_col + "_Z{}l2".format(z) + "}"

        df = df.Define("{}_Z1".format(lep_col), def_by_z(1))
        df = df.Define("{}_Z2".format(lep_col), def_by_z(2))

        df = df.Redefine(lep_col, "ROOT::VecOps::Concatenate({}_Z1, {}_Z2)".format(lep_col, lep_col))

        return df

    def get_hist_info(self, reg, prop):
        hist_info = self.cfg["hist_info"][prop]
        if "mass" in prop:
            if reg == "SR" or "HighMass" in reg:
                hist_info = hist_info["SR"]
            else:
                reg = reg.split("_")
                hist_info = hist_info[reg[-1]]
        
        return hist_info
    
    def get_column_name(self, reg, prop):
        if "Lepton" in prop:
            return prop
        else:
            return "{}_{}".format(self.cand(reg), prop)

    def write_hist(self, df, reg, prop, hist_info):
        column = self.get_column_name(reg, prop)

        if not self.isData:
            hist = df.Histo1D((prop, column, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), column, "weight")
            hist.Scale(self.lumi)
            hist = hist.GetValue()
        else:
            hist = df.Histo1D((prop, column, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), column).GetValue()
            hist.SetBinErrorOption(ROOT.TH1.kPoisson)

        return hist

    def write_hists(self):
        with up.recreate(self.outfile) as OutFile:
            hists = {}
            for sample, sample_path in tqdm(self.samples.items(), desc = "Processes", position = 0):
                if "Data" in sample_path: self.isData = True
                else: self.isData = False

                df = self.get_df(sample_path)

                hists[sample] = {}
                for reg in tqdm(self.regions, desc = "Regions", position = 1, leave = False):
                    if not self.isData:
                        df_reg = self.write_weight(df, reg)
                    else: 
                        df_reg = df.Filter(self.reg_filt(reg))
                    
                    hists[sample][reg] = {}               
                    for prop in tqdm(self.props, desc = "Properties", position = 2, leave = False):
                        hist_info = self.get_hist_info(reg, prop)

                        # Define lepton columns
                        if "Lepton" in prop:
                            df_4l = self.lep_df(df_reg, reg, prop)
                        else:
                            df_4l = df_reg

                        hists[sample][reg][prop] = {}

                        fs_dfs = dict(
                            fs_4l    = df_4l,
                            fs_4e    = self.fs_df(df_4l, "fs_4e", reg),
                            fs_4mu   = self.fs_df(df_4l, "fs_4mu", reg),
                            fs_2e2mu = self.fs_df(df_4l, "fs_2e2mu", reg),
                            fs_2mu2e = self.fs_df(df_4l, "fs_2mu2e", reg)
                        )

                        for fs, df_fs in tqdm(fs_dfs.items(), desc = "Final States", position = 3, leave = False):
                            OutFile["{}/{}/{}/{}".format(sample, reg, prop, fs)] = self.write_hist(df_fs, reg, prop, hist_info)
        
        return hists

if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--reg", choices=("SR", "OS_NoSIP_HighMass", "OS_NoSIP_MidMass", "OS_NoSIP_LowMass", "SS_NoSIP_HighMass", "SS_NoSIP_MidMass", "SS_NoSIP_LowMass"), default="SR")
    parser.add_argument("--prop", default="mass")
    parser.add_argument("--year", choices=(2022, 2023), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "B"), default="EFG")
    parser.add_argument("--tag", default="")
    args = vars(parser.parse_args())

    with open("/afs/cern.ch/user/i/iehle/cmssw/CMSSW_13_3_3/src/ZZAnalysis/NanoAnalysis/scripts/hist_config.yaml") as config:
        cfg = yaml.safe_load(config)

    zzNLO = HistWriter(cfg, args)
    hists = zzNLO.write_hists()