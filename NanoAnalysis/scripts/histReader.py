import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT
import uproot as up
from tqdm import tqdm

from NanoAnalysis.scripts.helper_functions import *

class HistReader:
    def __init__(self, cfg, args):
        self.cfg  = cfg
        self.args = args

        self.regions = cfg["regions"]
        self.props   = cfg["hist_info"]
        self.procs   = get_samples(cfg, args["year"], args["era"])

        self.fstates = dict(
            fs_4l    = (),
            fs_4e    = (-121, -121),
            fs_4mu   = (-169, -169),
            fs_2e2mu = (-121, -169),
            fs_2mu2e = (-169, -121)        
        )

    def get_fileDepth(self, hist_keys):
        """Iterate through list of histogram keys,
        where directories are separated by slashes.
        Returns the total depth, i.e.:
        <Process>/<Region>/<Final State>/<Variable> --> 4"""
        for i, key in enumerate(hist_keys):
            current_depth = key.count("/") + 1
            next_depth = hist_keys[i+1].count("/") + 1
            if current_depth == next_depth:
                return current_depth

    def get_hist_order(self, structure):
        template = []
        for el in structure:
            if el in self.procs:
                template.append("Procs")
            elif el in self.regions:
                template.append("Regions")
            elif el in self.fstates:
                template.append("FStates")
            elif el in self.props:
                template.append("Props")
            else:
                raise Exception("Unrecognized key: ", el)
        return template

    def get_final_key(self, order, reg, proc, prop, fs):
        """Returns key corresponding to reg, proc, prop, fs
        in the order specified by the passed order template."""
        key = []
        for el in order:
            key.append(el.replace("Procs",proc).replace("Regions",reg).replace("Props",prop).replace("FStates",fs))
        return "/".join(key)
            
    def read_hists_and_counts(self, infile="", procs={}):
        """Reads hists from file and arranges them in
        dictionary(ies) corresponding to
        <Regions>/<Properties>/<Final States>/<Processes>."""
        if procs != {}: self.procs = procs["MC"] | procs["Pol"] | procs["Data"]
        with up.open(infile) as Hists:
            all_keys     = list(Hists.keys())
            depth        = self.get_fileDepth(all_keys)
            ex_structure = all_keys[depth].replace(";1", "").split("/")
            order = self.get_hist_order(ex_structure)

            hists  = {}
            counts = {}
            errors = {}
            for reg in self.regions:
                hists[reg]  = {}
                counts[reg] = {}
                errors[reg] = {}
                for prop in self.props:
                    hists[reg][prop]  = {}
                    counts[reg][prop] = {}
                    errors[reg][prop] = {}
                    for fs in self.fstates:
                        hists[reg][prop][fs]  = dict(MC = {}, Pol = {})
                        counts[reg][prop][fs] = dict(MC = {}, Pol = {})
                        errors[reg][prop][fs] = dict(MC = {}, Pol = {})
                        for proc in self.procs:
                            key = self.get_final_key(order, reg, proc, prop, fs)
                            hist = Hists[key]

                            np_hist = hist.to_numpy(flow = True)
                            cnts  = round(np.sum(np_hist[0]), 3)
                            err     = np.sqrt(hist.variances(flow=True))

                            if proc == "Data":
                                hists[reg][prop][fs]["Data"]  = np_hist
                                counts[reg][prop][fs]["Data"] = cnts
                                errors[reg][prop][fs]["Data"] = err

                            elif proc in ["ZLZL", "ZLZT", "ZTZT", "ZZ_LO"]:
                                hists[reg][prop][fs]["Pol"][proc]  = np_hist
                                counts[reg][prop][fs]["Pol"][proc] = cnts
                                errors[reg][prop][fs]["Pol"][proc] = err

                            else:
                                hists[reg][prop][fs]["MC"][proc]  = np_hist
                                counts[reg][prop][fs]["MC"][proc] = cnts
                                errors[reg][prop][fs]["MC"][proc] = err

        return hists, counts, errors

if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--reg", choices=("SR", "OS_NoSIP_HighMass", "OS_NoSIP_MidMass", "OS_NoSIP_LowMass", "SS_NoSIP_HighMass", "SS_NoSIP_MidMass", "SS_NoSIP_LowMass"), default="SR")
    parser.add_argument("--prop", default="mass")
    parser.add_argument("--year", choices=(2022, 2023), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "B"), default="EFG")
    parser.add_argument("--tag", default="")
    parser.add_argument("--infile", default="")
    parser.add_argument("--lumi_tag", default=0, type=int)
    args = vars(parser.parse_args())

    with open("/afs/cern.ch/user/i/iehle/cmssw/CMSSW_13_3_3/src/ZZAnalysis/NanoAnalysis/scripts/hist_config.yaml") as config:
        cfg = yaml.safe_load(config)
  
    histManager = HistManager(cfg, args)
    histManager.plot_hists()