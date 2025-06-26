import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

from NanoAnalysis.scripts.fileHandler import FileHandler
from NanoAnalysis.scripts.histSystematics import HistSystematics

import ROOT

import uproot as up
import numpy as np
from tqdm import tqdm

ROOT.gInterpreter.Declare("""

#include <string>

string get_event_idx(int run, int lumi, long event) {

    string run_string = to_string(run);
    string lumi_string = to_string(lumi);
    string event_string = to_string(event);

    string event_idx = run_string + "_" + lumi_string + "_" + event_string;

    return event_idx;
}

""")

# See: https://gist.github.com/eguiraud/77a0ca3566e66bc6b8cd0f9e156c983b
# and https://root-forum.cern.ch/t/select-unique-candidates-based-on-their-id-and-variable/59668/3
ROOT.gInterpreter.Declare("""
// A thread-safe stateful filter that lets only one event pass for each value of
// "category" (where "category" is a random character).
// It is using gCoreMutex, which is a read-write lock, to have a bit less contention between threads.

class FilterOnePerKind {
  std::unordered_set<string> _seenCategories;
  
public:
  bool operator()(string category) {
    {
      R__READ_LOCKGUARD(ROOT::gCoreMutex); // many threads can take a read lock concurrently
      if (_seenCategories.count(category) == 1)
        return false;
    }
    // if we are here, `category` was not already in _seenCategories
    R__WRITE_LOCKGUARD(ROOT::gCoreMutex); // only one thread at a time can take the write lock
    _seenCategories.insert(category);
    return true;
  }
};
""")

class HistWriter:
    def __init__(self, cfg, lumi, col_tag):
        self.cfg  = cfg
        self.lumi = lumi

        self.cand    = lambda reg: "ZZCand" if reg == "SR" else "ZLLCand"
        self.reg_idx = lambda reg: "bestCandIdx" if reg == "SR" else f"ZLLbest{reg}Idx"

        self.reg_prop = lambda prop, reg: f"{self.cand(reg)}_{prop}"

        self.lep_idx = lambda z, l, reg: f"{self.col_tag}_{z}{l}Idx_{reg}"
        self.lep_col = lambda lep_prop, z, l, reg: f"{self.col_tag}_{lep_prop}_{z}{l}_{reg}"

        self.zpx     = ["DY", "TT", "WZ"]

        self.syst_variables = ["cosTheta1", "cosTheta3", "cosThetaStar", "delRapidity", "delPhi", "delPhiStar"]

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

    def run_systs(self, hist_info, reg):
        return ("systematics" in hist_info) and (reg=="SR") and (self.category not in self.zpx) and (not self.data)

    def _get_df(self, process, filepath):

        df = ROOT.RDataFrame("Events", filepath)

        ROOT.RDF.Experimental.AddProgressBar(df)

        df = df.Filter("HLT_passZZ4l")

        # Only needed temporarily, filter events that belong to at least one reg
        reg_bools = [f"({self.reg_idx(reg)} >= 0)" for reg in self.regions]
        reg_filt  = " || ".join(reg_bools)
        df = df.Filter(reg_filt)

        if process == "Data":
            return df

        Runs = ROOT.RDataFrame("Runs", filepath)

        df = df.Define("genEventSumw",str(Runs.Sum("genEventSumw").GetValue()))

        # Temporarily needed to handle duplicate events!
        df = df.Define("event_idx", "get_event_idx(run, luminosityBlock, event)")
        
        event_idx = ROOT.std.vector['string'](["event_idx"])
        
        return df.Filter(ROOT.FilterOnePerKind(), event_idx)

    def _combine_procs(self, category, sub_hists):
        processes = list(sub_hists.keys())

        combined_hists = {}
        for prop in self.properties:
            combined_hists[prop] = {}
            for reg in self.regions:
                combined_hists[prop][reg] = {}

                for fs in sub_hists[processes[0]][prop][reg].keys():
                # for fs in self.final_states:

                    new_hist = sub_hists[processes[0]][prop][reg][fs].Clone(category)
                    for proc in processes[1:]:
                        if isinstance(sub_hists[proc][prop][reg][fs], ROOT.TH1D):
                            new_hist.Add(sub_hists[proc][prop][reg][fs])
                        else:
                            new_hist.Add(sub_hists[proc][prop][reg][fs].GetValue())

                    combined_hists[prop][reg][fs] = new_hist
        
        return combined_hists

    def all_hists(self):
        self.hists = {}

        for key in self.file_paths: #MC, Pol, or Data
            if len(self.file_paths[key]) == 0: continue
            if key == "Data":
                self.data = True
                
                df = self._get_df(key, self.file_paths[key])
                
                self.hists[key] = self.write_hists(df)
            
            else:
                self.data = False

                key_paths = self.file_paths[key]

                if len(key_paths) != 0:
                    
                    for category in key_paths:

                        cat_paths = key_paths[category]

                        # Save to use for systematics
                        self.category = category
                        
                        if isinstance(cat_paths, dict):
                            sub_hists = {}    
                            for process, filepath in cat_paths.items():
                                
                                df = self._get_df(process, filepath)

                                sub_hists[process] = self.write_hists(df)
                            
                            self.hists[category] = self._combine_procs(category, sub_hists)

                        else:
                            
                            df = self._get_df(category, key_paths[category])

                            self.hists[category] = self.write_hists(df)
                            print("Wow!")

    def def_lep_id_cols(self, df, reg):
        for Z in ["Z1", "Z2"]:
            for l in ["l1", "l2"]:
                lep_idx          = self.reg_prop(f"{Z}{l}Idx", reg)
                good_lep_idx     = f"{lep_idx}_{reg}"
                good_lep_idx_def = f"{lep_idx}[{self.good_reg_idx}]"

                df = df.Define(good_lep_idx, good_lep_idx_def)
        
        return df

    def define_reg_cols(self, df, reg):
        z1_flav = self.reg_prop("Z1flav", reg)
        z2_flav = self.reg_prop("Z2flav", reg)

        self.good_reg_idx = self.reg_idx(reg)
        
        if not self.data:
            dMC_wgt           = self.reg_prop("dataMCWeight", reg)
            self.good_dMC_wgt = f"{dMC_wgt}_{reg}"
            df                = df.Define(self.good_dMC_wgt, f"{dMC_wgt}[{self.good_reg_idx}]")

        self.good_z1_flav = f"{z1_flav}_{reg}"
        self.good_z2_flav = f"{z2_flav}_{reg}"

        df = df.Define(self.good_z1_flav, f"{z1_flav}[{self.good_reg_idx}]").Define(self.good_z2_flav, f"{z2_flav}[{self.good_reg_idx}]")

        df = self.def_lep_id_cols(df, reg)
        
        return df

    def lep_cols(self, df, reg, lep_prop):
        # Ex: lep_prop is usally Lepton_sip3d_Z2, need to extract just sip3d
        prop = lep_prop.split("_")[1]

        lep_col = "Lepton_{}".format(prop)
        if not lep_col in df.GetColumnNames():
            df = df.Define(lep_col, f"ROOT::VecOps::Concatenate(Electron_{prop}, Muon_{prop})")

        for Z in ["Z1", "Z2"]:
            for l in ["l1", "l2"]:
                lep_idx          = self.reg_prop(f"{Z}{l}Idx", reg)
                this_lep_idx     = f"{lep_idx}_{reg}"

                this_lep_col     = f"{lep_col}_{Z}{l}_{reg}"
                this_lep_col_def = f"{lep_col}[{this_lep_idx}]"

                df = df.Define(this_lep_col, this_lep_col_def)

            good_lep_col_Z     = f"{lep_col}_{Z}_{reg}"
            good_lep_col_Z_def = "ROOT::VecOps::RVec<float> {" + f"{lep_col}_{Z}l1_{reg}, " + f"{lep_col}_{Z}l2_{reg}" + "}"
            
            df = df.Define(good_lep_col_Z, good_lep_col_Z_def)
        
        return df.Define(f"{lep_col}_{reg}", f"ROOT::VecOps::Concatenate({lep_col}_Z1_{reg}, {lep_col}_Z2_{reg})")

    def fs_filt(self, df, reg, fs):

        if "4l" in fs:
            return df.Filter(f"({self.good_z1_flav} == -121) || ({self.good_z1_flav} == -169)")

        elif "2x2e" in fs:
            pdg2 = -121
        elif "2x2mu" in fs:
            pdg2 = -169
        
        if ("SSSIP" in reg) or ("SSRelaxed" in reg):
            pdg2 *= -1
        
        return df.Filter(f"{self.good_z2_flav}=={pdg2}")

    def write_hists(self, df):
        hists = {}

        hist_list = []
        
        for reg in tqdm(self.regions, desc = "Regions", position = 0):
            
            df = self.define_reg_cols(df, reg) # Define Z1flav, Z2flav, and dataMCWeight columns by region
            
            if not self.data:
                weight_col     = f"weight_{reg}"
                weight_col_def = f"{self.good_dMC_wgt}*overallEventWeight/genEventSumw"
                df             = df.Define(weight_col, weight_col_def)

            hists[reg] = {}
            for prop in tqdm(self.properties, desc = "Properties", position = 1, leave = False):
                hists[reg][prop] = {}

                hist_info = self.get_hist_info(reg, prop)


                if "Lepton" not in prop:
                    branch_name = self.reg_prop(prop, reg)
                    good_branch = f"{branch_name}_{reg}"

                    good_branch_def = f"{branch_name}[{self.good_reg_idx}]"
                    df              = df.Define(good_branch, good_branch_def)

                else:
                    good_branch = f"{prop}_{reg}"
                    df = self.lep_cols(df, reg, prop)

                for fs in tqdm(self.final_states, desc = "Final States", position = 2, leave = False):
                    df_fs = self.fs_filt(df, reg, fs)

                    if not self.data:
                        hist = df_fs.Histo1D((good_branch+"_nom", good_branch, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), good_branch, weight_col) 
                    else:
                        hist = df_fs.Histo1D((good_branch+"_nom", good_branch, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), good_branch)

                    hist_list.append(hist)

                    hists[reg][prop][fs] = hist

                    # if not self.data and (self.category not in self.zpx) and reg == "SR":
                    if self.run_systs(hist_info, reg):

                        systWriter = HistSystematics(hist_info, self.lumi, good_branch, weight_col, self.good_reg_idx)

                        # if self.category == "ggZZ":
                        #     variations = ["puWeight", "lepIDRec"]
                        # elif self.category == "VVV":
                        #     variations = ["puWeight", "lepIDRec", "LHEScaleWeight"]

                        if self.category == "ggZZ" or self.category == "VVV":
                            variations = ["puWeight", "lepIDRec"]
                        else:
                            variations = ["puWeight", "lepIDRec", "LHEScaleWeight", "LHEPdfWeight"]

                        for var in variations:
                            syst_vars = systWriter.get_vars(df_fs, var)
                            if "LHE" not in var:
                                up_hist, down_hist = syst_vars
                                
                                hists[reg][prop][f"{fs}_{var}Up"]   = up_hist
                                hists[reg][prop][f"{fs}_{var}Down"] = down_hist
                            else:
                                hists[reg][prop][f"{fs}_{var}"] = syst_vars # contains all (booked) variations, need to calculate Up/Down after event loop

                                # qcd_up, qcd_down = lhe_vars["QCD"]
                                # pdf_up, pdf_down = lhe_vars["PDF"]

                                # hists[reg][prop][f"{fs}_LHEScaleWeightUp"]   = qcd_up
                                # hists[reg][prop][f"{fs}_LHEScaleWeightDown"] = qcd_down

                                # hists[reg][prop][f"{fs}_LHEPdfWeightUp"]   = pdf_up
                                # hists[reg][prop][f"{fs}_LHEPdfWeightDown"] = pdf_down

        try:
            ROOT.RDF.RunGraphs(hist_list)
        except:
            breakpoint()
        
        final_hists = {}
        for prop in self.properties:
            final_hists[prop] = {}
            for reg in self.regions:
                final_hists[prop][reg] = {}
                
                hist_info = self.get_hist_info(reg, prop)

                # for fs in self.final_states:
                #     hist = hists[reg][prop][fs]
                    
                #     if not self.data:
                #         hist.Scale(self.lumi)
                #     else:
                #         hist.SetBinErrorOption(ROOT.TH1.kPoisson)       
                    
                #     final_hists[prop][reg][fs] = hist.GetValue()

                #     if self.run_systs(hist_info, reg):
                #         for var in hist_info["systematics"]["vars"]:
                #             final_hists[prop][reg][f"{fs}_{var}Up"] = hists[reg][prop][f"{fs}_{var}Up"]
                #             final_hists[prop][reg][f"{fs}_{var}Down"] = hists[reg][prop][f"{fs}_{var}Down"]

                for key, val in hists[reg][prop].items():
                    
                    if "LHEScaleWeight" in key:
                        up_hist, down_hist = systWriter.qcd_up_down(val, hist_info)

                        up_hist.Scale(self.lumi)
                        down_hist.Scale(self.lumi)

                        final_hists[prop][reg][key+"Up"]   = up_hist
                        final_hists[prop][reg][key+"Down"] = down_hist
                    
                    elif "LHEPdfWeight" in key:
                        up_hist, down_hist = systWriter.pdf_up_down(val, hist_info)

                        up_hist.Scale(self.lumi)
                        down_hist.Scale(self.lumi)

                        final_hists[prop][reg][key+"Up"]   = up_hist
                        final_hists[prop][reg][key+"Down"] = down_hist
                    
                    else:
                        hist = val
                    
                        if not self.data:
                            hist.Scale(self.lumi)
                        else:
                            hist.SetBinErrorOption(ROOT.TH1.kPoisson)  

                        final_hists[prop][reg][key] = hist     
                    
                    # if key in self.final_states:
                    #     final_hists[prop][reg][key] = hist.GetValue()
                    # else:
                    #     final_hists[prop][reg][key] = hist
        
        return final_hists

    def main(self, file_paths, regions, properties, final_states):
        self.file_paths   = file_paths

        self.regions      = regions
        self.properties   = properties
        self.final_states = final_states

        self.all_hists()

        return self.hists
    
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