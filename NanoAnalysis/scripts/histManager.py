import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT

from tqdm import tqdm
import uproot as up
import numpy as np

from NanoAnalysis.scripts.fileHandler import FileHandler
from NanoAnalysis.scripts.histWriter import HistWriter
from NanoAnalysis.scripts.histPlotter import HistPlotter
from NanoAnalysis.scripts.ZpX_estimation import ZpX

class HistManager:
    def __init__(self, cfg, args):
        self.cfg          = cfg

        self.year         = args["year"]
        self.era          = args["era"]

        self.tag          = args["tag"]

        self.file_handler = FileHandler(self.year, self.era, self.tag)
        self.lumi         = self.file_handler.lumi

        self.get_cfg_path = lambda step: os.path.join(parent_dir, "NanoAnalysis/scripts", self.cfg[step]["cfg_file"])

        self.zpx = ZpX()

        self.hists = {}

    def _get_cfg(self, step):
        cfg_path = self.get_cfg_path(step)
        with open(cfg_path) as config:
            step_cfg = yaml.safe_load(config)
        return step_cfg, self.cfg[step]

    def _get_hist_structure(self, Hists):
        props   = []
        regions = []
        fstates = []

        for key in Hists.keys():
            full_key = key.split("/")
            # Write dict structure
            if len(full_key) == 3:
                prop, reg, fs = full_key
                fs = fs.replace(";1", '')
                if prop not in props: props.append(prop)
                if reg not in regions:regions.append(reg)
                if fs not in fstates: fstates.append(fs)

        return props, regions, fstates

    def write_hists(self):
        write_cfg, step_cfg = self._get_cfg("writing")

        regions = step_cfg["regions"]
        fstates = step_cfg["fstates"]
        col_tag = step_cfg["col_tag"]

        props   = write_cfg.keys()

        hist_writer = HistWriter(write_cfg, self.lumi, col_tag)

        hists = hist_writer.main(self.file_handler.file_paths, regions, props, fstates)
        self.file_handler.write_hists(hists)

        return hists

    def read_hists(self, path):
        props   = []
        regions = []
        fstates = []

        hists, counts, errors = {}, {}, {}
        with up.open(path) as Hists:
            props, regions, fstates = self._get_hist_structure(Hists)
            for prop in props:
                hists[prop] = {}
                counts[prop] = {}
                errors[prop] = {}
                for reg in regions:
                    hists[prop][reg] = {}
                    counts[prop][reg] = {}
                    errors[prop][reg] = {}
                    for fs in fstates:
                        hists[prop][reg][fs] = dict(
                            MC   = {},
                            Data = {},
                            Pol  = {},
                        )
                        counts[prop][reg][fs] = dict(
                            MC   = {},
                            Data = {},
                            Pol  = {},
                        )
                        errors[prop][reg][fs] = dict(
                            MC   = {},
                            Data = {},
                            Pol  = {},
                        )
                        for proc, hist in Hists[prop][reg][fs].items():
                            proc = proc.replace(";1","")

                            if proc in self.file_handler.mc_samples.keys():
                                hists[prop][reg][fs]["MC"][proc]  = hist.to_numpy(flow=True)
                                counts[prop][reg][fs]["MC"][proc] = np.sum(hist.values(flow=True))
                                errors[prop][reg][fs]["MC"][proc] = hist.errors(flow=True)

                            elif proc in self.file_handler.data_samples.keys():
                                hists[prop][reg][fs]["Data"][proc]  = hist.to_numpy(flow=True)
                                counts[prop][reg][fs]["Data"][proc] = np.sum(hist.values(flow=True))
                                errors[prop][reg][fs]["Data"][proc] = hist.errors(flow=True)

                            elif proc in self.file_handler.pol_samples.keys():
                                hists[prop][reg][fs]["Pol"][proc]  = hist.to_numpy(flow=True)
                                counts[prop][reg][fs]["Pol"][proc] = np.sum(hist.values(flow=True))
                                errors[prop][reg][fs]["Pol"][proc] = hist.errors(flow=True)
            
        return hists, counts, errors
    
    def plot_hists(self, path=""):
        if path=="": path = self.file_handler.hist_path
        if not os.path.exists(path):
            self.write_hists()

        hists, counts, errors = self.read_hists(path)

        hist_plotter = HistPlotter(self.year, self.era, self.tag, self.lumi, hists, counts, errors)
        figs = hist_plotter.main()
        breakpoint()

        self.file_handler.write_plots(figs)
        #fig = hist_plotter.plotter("mass", "SR", "fs_4mu")
        print("Yay!")

    def write_zpx(self, path=""):
        if path=="": path = self.file_handler.hist_path
        if not os.path.exists(path):
            self.write_hists()

        hists, counts, errors = self.read_hists(path)

        yields = self.zpx.get_yields(hists, errors, ["fs_4l", "fs_4e", "fs_4mu", "fs_2e2mu", "fs_2mu2e"])
        breakpoint()
        print("Wow!")

   # def plot_zpx(self, infile_1, years, eras, infile_2=""):
    def plot_zpx(self, infile_1, infile_2=""):
        self.zpx = ZpX()
        
        hists_1, counts_1, errors_1 = self.read_hists(infile_1)

        zpx_info_1 = self.zpx.get_zpx(hists_1, errors_1, ["fs_4e", "fs_4mu", "fs_2e2mu", "fs_2mu2e"])
        if infile_2 != "":
            hists_2, counts_2, errors_2 = self.histReader.read_hists_and_counts(infile_2)
            hists_2, counts_2, errors_2 = self.combine_processes(hists_2, counts_2, errors_2)
            zpx_info_2 = self.zpx.get_zpx(hists_2, errors_2, self.fstates)

            categories = self.fstates
            counts = {
                "2022": [zpx_info_1["N_ZpX_MidMass"][fs][0] for fs in self.fstates],
                "2023": [zpx_info_2["N_ZpX_MidMass"][fs][0] for fs in self.fstates],
            }
            errs = {
                "2022": [zpx_info_1["N_ZpX_MidMass"][fs][1] for fs in self.fstates],
                "2023": [zpx_info_2["N_ZpX_MidMass"][fs][1] for fs in self.fstates],
            }
            
            lumi_2022 = self.cfg["datasets"]["year_2022"]["Full"]["Lumi"]
            lumi_2023 = self.cfg["datasets"]["year_2023"]["Full"]["Lumi"]

            norm_counts = {
                "2022": [cnt/(lumi_2022*1e-3) for cnt in counts["2022"]],
                "2023": [cnt/(lumi_2023*1e-3) for cnt in counts["2023"]]
            }
            # Assuming a 1.5% lumi err for full 2022 and 2023
            norm_errs = {
                "2022": [abs(nm_cnt)*np.sqrt((err/cnt)**2 + (0.015)**2) for nm_cnt, err, cnt in zip(norm_counts["2022"], errs["2022"], counts["2022"])],
                "2023": [abs(nm_cnt)*np.sqrt((err/cnt)**2 + (0.015)**2) for nm_cnt, err, cnt in zip(norm_counts["2023"], errs["2023"], counts["2023"])]
            }

            x = np.arange(len(self.fstates))
            group_width = 0.5
            offset_step = group_width/len(counts)

            fig, ax = plt.subplots(layout='constrained')
            for i, (key, val) in enumerate(counts.items()):
                offset = x - (group_width - offset_step)/2 + i*offset_step
                ax.errorbar(
                    offset,
                    norm_counts[key],
                    yerr=norm_errs[key],
                    fmt="o",
                    label=key
                )
            
            ax.legend()

            ax.set_xticks(x)
            ax.set_xticklabels(self.fstates)

            ax.set_ylabel(r"$N_{Z+X}/{fb^{-1}}$")
            title = "N_ZpX/fb^-1 Full 2022, 2023"
            outfile = "N_ZpX_Full_2022_2023_perInvFb_Z1pt"
            ax.set_title(title)
            fig.savefig(outfile+".png")
        
        else:
            for step in zpx_info_1.keys():
                #era, year = eras[0], years[0]
                era, year = "EFG", 2022
                fig = self.zpx.plot_zpx(zpx_info_1, step, year, era)
                fig.savefig(f"zpx_test_{step}_{year}_{era}.png")

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
    #hist_manager.write_hists()
    #hist_manager.plot_hists()
    #hist_manager.plot_hists("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists.root")
    #hist_manager.write_zpx("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_wLepCols.root")
    hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_wLepCols.root")