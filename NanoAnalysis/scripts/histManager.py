import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

from pathlib import Path

import ROOT

from tqdm import tqdm
import uproot as up
import numpy as np

from NanoAnalysis.scripts.fileHandler import FileHandler
from NanoAnalysis.scripts.histWriter import HistWriter
from NanoAnalysis.scripts.histPlotter import HistPlotter
from NanoAnalysis.scripts.ZpX_estimation import ZpX

import matplotlib.pyplot as plt

from copy import deepcopy

class HistManager:
    def __init__(self, cfg, args):
        self.cfg          = cfg

        self.year         = args["year"] if args["year"] != -1 else "Full"
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

    def rebin_hist(self, hist, new_bins, overflow, name=""):
        pyhist = hist.to_pyroot()

        nbins = len(new_bins) - 1 
        
        rebinned = pyhist.Rebin(nbins, name, np.array(new_bins))

        max_bin = nbins + 2 if overflow else nbins+1
        #count   = rebinned.Integral(0, max_bin)

        bin_contents = np.array([rebinned.GetBinContent(i) for i in range(1, max_bin)])
        bin_errors   = np.array([rebinned.GetBinError(i) for i in range(1, max_bin)])
        
        new_bins = np.array(new_bins)
        if overflow:
            new_bins = np.append(new_bins, np.inf)

        count = bin_contents.sum()

        return (bin_contents, new_bins), count, bin_errors

    def combine_procs(self, hists, counts, errors, **kwargs):
        for prop in hists.keys():
            for reg in hists[prop].keys():
                for fs in hists[prop][reg].keys():
                    mc_hists  = hists[prop][reg][fs]["MC"]
                    mc_counts = counts[prop][reg][fs]["MC"]
                    mc_errors = errors[prop][reg][fs]["MC"]

                    for group, procs in kwargs.items():
                        group_count = np.array([mc_counts[proc] for proc in procs]).sum()

                        group_hist, edges = mc_hists[procs[0]]

                        group_error = np.square(mc_errors[procs[0]])

                        # del mc_hists[procs[0]]
                        # del mc_counts[procs[0]]
                        # del mc_errors[procs[0]]

                        for proc in procs[1:]:
                            group_hist += mc_hists[proc][0]
                            
                            group_error += np.square(mc_errors[proc])

                            # del mc_hists[proc]
                            # del mc_counts[proc]
                            # del mc_errors[proc]

                        group_hist = (group_hist, edges)
                        group_error = np.sqrt(group_error)

                        hists[prop][reg][fs]["MC"][group] = group_hist
                        counts[prop][reg][fs]["MC"][group] = group_count
                        errors[prop][reg][fs]["MC"][group] = group_error

        return hists, counts, errors
                    

    def read_hists(self, path, **kwargs):
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
                #for reg in self.cfg["writing"]["regions"]:
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

                            # We want to keep overflow but not underflow bin
                            if len(kwargs) > 0:
                                if kwargs["rebin"] and prop in kwargs:
                                    if prop == "mass":
                                        if reg == "SR" or "HighMass" in reg:
                                            bin_info = kwargs[prop]["HighMass"]
                                        elif "MidMass" in reg:
                                            bin_info = kwargs[prop]["MidMass"]
                                        elif "LowMass" in reg:
                                            bin_info = kwargs[prop]["LowMass"]
                                    else:
                                        bin_info = kwargs[prop]

                                    this_hist, count, err = self.rebin_hist(hist, bin_info["bins"], bin_info["overflow"])

                                else:
                                    this_hist = hist.to_numpy(flow=True)
                                    this_hist = (this_hist[0][1:], this_hist[1][1:])
                                    count     = this_hist[0].sum()
                                    err       = hist.errors(flow=True)[1:]
                            else:
                                this_hist = hist.to_numpy(flow=True)
                                this_hist = (this_hist[0][1:], this_hist[1][1:])
                                count     = this_hist[0].sum()
                                err       = hist.errors(flow=True)[1:]

                            # this_hist = hist.to_numpy(flow=True)
                            # count     = np.sum(hist.values(flow=True))
                            # err       = hist.errors(flow=True)

                            # this_hist = hist.to_numpy(flow=False)
                            # count     = np.sum(hist.values(flow=False))
                            # err       = hist.errors(flow=False)

                            if proc in list(self.file_handler.mc_procs.keys()) + ["ZpX"]:

                                hists[prop][reg][fs]["MC"][proc]  = this_hist
                                counts[prop][reg][fs]["MC"][proc] = count
                                errors[prop][reg][fs]["MC"][proc] = err

                            elif proc == "Data":
                                hists[prop][reg][fs]["Data"]  = this_hist
                                counts[prop][reg][fs]["Data"] = count
                                errors[prop][reg][fs]["Data"] = err

                            elif proc in self.file_handler.pol_procs:
                                hists[prop][reg][fs]["Pol"][proc]  = this_hist
                                counts[prop][reg][fs]["Pol"][proc] = count
                                errors[prop][reg][fs]["Pol"][proc] = err

        if "group" in kwargs and kwargs["group"]:
            hists, counts, errors = self.combine_procs(hists, counts, errors, **kwargs["groups"])
            
        return hists, counts, errors
    
    def plot_hists(self, path="", **kwargs):
        if path=="": path = self.file_handler.hist_path
        if not os.path.exists(path):
            self.write_hists()

        hists, counts, errors = self.read_hists(path, **kwargs)

        hist_plotter = HistPlotter(self.year, self.era, self.tag, self.lumi, hists, counts, errors)
        figs = hist_plotter.main()

        self.file_handler.write_plots(figs)
        #fig = hist_plotter.plotter("mass", "SR", "fs_4mu")
        print("Yay!")

    def combine_eras(self, infile_1, infile_2, years, eras):

        if years[0]==years[1]:
            base_dir, filename = os.path.split(infile_1)
            outdir = base_dir.replace(eras[0], "Full")
            Path(outdir).mkdir(parents=True, exist_ok=True)
            outfile = os.path.join(outdir, filename)
        else:
            # NEED TO CHANGE THIS
            outfile = "/eos/user/i/iehle/Analysis/histograms/Full/hists_zpxForPlots_11_07_25.root"

        with up.open(infile_1) as Hists_1, up.open(infile_2) as Hists_2, up.recreate(outfile) as NewHists:
            for key in tqdm(Hists_1.keys()):
                if key.count("/")==3 and key in Hists_2.keys():
                    hist_1, hist_2 = Hists_1[key], Hists_2[key]
                    NewHists[key.replace(";1", "")] = hist_1.to_pyroot() + hist_2.to_pyroot()
            print("yay!")

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

        #fstates = ["fs_4e", "fs_4mu", "fs_2e2mu", "fs_2mu2e"]
        fstates = ["fs_4l", "fs_2x2e", "fs_2x2mu"]

        zpx_info_1 = self.zpx.get_zpx(hists_1, errors_1, fstates)

        import json
        with open("ZpX_info_2022_Full_24_06_25.json", "w") as myfile:
            json.dump(zpx_info_1, myfile, indent=4)

        if infile_2 != "":
            hists_2, counts_2, errors_2 = self.read_hists(infile_2)

            zpx_info_2 = self.zpx.get_zpx(hists_2, errors_2, fstates)

            categories = fstates
            counts = {
                "2022": [zpx_info_1["N_ZpX_MidMass"][fs][0] for fs in fstates],
                "2023": [zpx_info_2["N_ZpX_MidMass"][fs][0] for fs in fstates],
            }
            errs = {
                "2022": [zpx_info_1["N_ZpX_MidMass"][fs][1] for fs in fstates],
                "2023": [zpx_info_2["N_ZpX_MidMass"][fs][1] for fs in fstates],
            }
            
            #lumi_2022 = self.cfg["datasets"]["year_2022"]["Full"]["Lumi"]
            #lumi_2023 = self.cfg["datasets"]["year_2023"]["Full"]["Lumi"]

            lumi_2022 = 34.6532e3
            lumi_2023 = 27.245e3

            norm_counts = {
                "2022": [cnt/(lumi_2022*1e-3) for cnt in counts["2022"]],
                "2023": [cnt/(lumi_2023*1e-3) for cnt in counts["2023"]]
            }

            norm_errs = {
                "2022": [abs(nm_cnt)*np.sqrt((err/cnt)**2 + (0.014)**2) for nm_cnt, err, cnt in zip(norm_counts["2022"], errs["2022"], counts["2022"])],
                "2023": [abs(nm_cnt)*np.sqrt((err/cnt)**2 + (0.013)**2) for nm_cnt, err, cnt in zip(norm_counts["2023"], errs["2023"], counts["2023"])]
            }

            x = np.arange(len(fstates))
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
            ax.set_xticklabels(fstates)

            ax.set_ylabel(r"$N_{Z+X}/{fb^{-1}}$")
            title = "N_ZpX/fb^-1 Full 2022, 2023"
            outfile = "N_ZpX_Full_2022_2023_perInvFb_noSmartCut_2x2l"
            ax.set_title(title)
            fig.savefig(outfile+".png")
        
        else:
            for step in zpx_info_1.keys():
                #era, year = eras[0], years[0]
                #era, year = "CD", 2022
                #era, year = "EFG", 2022
                era, year = "Full", 2022
                #era, year = "C", 2023
                #era, year = "D", 2023
                #era, year = "Full", 2023
                fig = self.zpx.plot_zpx(zpx_info_1, step, year, era)
                fig.savefig(f"zpx_test_{step}_{year}_{era}_24_06_25.png")

if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--year", choices=(2022, 2023, -1), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "Full"), default="EFG")
    parser.add_argument("--tag", default="")
    args = vars(parser.parse_args())    

    cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/gen_cfg.yaml")
    with open(cfg_path) as config:
        cfg = yaml.safe_load(config)
    
    hist_manager = HistManager(cfg, args)

    # hist_manager.combine_eras(
    #     #infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_allShapeVar_24_06_25.root",
    #     #infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_sigOnTop_allFStates.root",
    #     #infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_24_06_25.root",
    #     # infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_dataInCRs.root",
    #     # infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_dataInCRs.root",
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_delRapidityForPlots_03_07_25_v2.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_zpxForPlots_11_07_25.root",
    #     years    = (2022, 2023),
    #     eras     = ("Full", "Full")
    # )

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/CD/hists_sigOnTop_allFStates.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_sigOnTop_allFStates.root",
    #     years    = (2022, 2022),
    #     eras     = ("CD", "EFG")
    # )

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2023/C/hists_allFstats_SROnly_massOnly.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/D/hists_allFstats_SROnly_massOnly.root",
    #     years    = (2023, 2023),
    #     eras     = ("C", "D")
    # )

    hist_manager.plot_hists(
        "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_allFstats_SROnly_massOnly.root",
        rebin=True,
        mass = dict(
            HighMass = dict(
                bins     = [180., 200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400., 420.],
                overflow = True,          
            )
        ),
        group = True,
        groups = dict(
            ZpX = ["DY", "TT", "WZ"]
        )
    )

    #hist_manager.write_hists()
    # #hist_manager.plot_hists()
    # hist_manager.plot_hists(
    #     # "/eos/user/i/iehle/Analysis/histograms/Full/hists_4l_02_07_25.root",
    #     #"/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_2x2e_2x2mu_4l.root",
    #     #"/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_24_06_25.root",
    #     #"/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_dataInCRs.root",
    #     #"/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_SSTest_28_06_25.root",
    #     #"/eos/user/i/iehle/Analysis/histograms/Full/hists_zpxForPlots_03_07_25_v2.root",
    #     #"/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_sigOnTop_allFStates.root",
    #     "/eos/user/i/iehle/Analysis/histograms/Full/hists_zpxForPlots_11_07_25.root",
    #     rebin=True,
    #     mass = dict(
    #         HighMass = dict(
    #             bins     = [180., 200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400., 420.],
    #             overflow = True,          
    #         ),
    #         MidMass = dict(
    #             bins     = [140., 150., 160., 170., 180.],
    #             overflow = False,          
    #         ), 
    #         LowMass = dict(
    #             bins     = [105., 115., 125., 135., 140.],
    #             overflow = False,          
    #         ), 
    #     ),
    #     Z1mass = dict(
    #         bins = [81., 83., 85., 87., 89., 91., 93., 95., 97., 99., 101., 103.],
    #         overflow = False
    #     ),
    #     Z2mass = dict(
    #         bins = [81., 83., 85., 87., 89., 91., 93., 95., 97., 99., 101., 103.],
    #         overflow = False
    #     ),
    #     Lepton_sip3d_Z2 = dict(
    #         bins     = [0.0, 4.0, 20.0],
    #         overflow = True,
    #     ),
    #     # delRapidity = dict(
    #     #     bins = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0],
    #     #     overflow = True,
    #     # ),
    #     # group = True,
    #     # groups = dict(
    #     #     ZpX = ["DY", "TT", "WZ"]
    #     # )
    # )
    # #hist_manager.plot_hists("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists.root")
    #hist_manager.write_zpx("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_wLepCols.root")
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_wLepCols.root")

    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_fixCROverlap_wData.root")
    
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_NewAlg_noSmartCut.root")

    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/CD/hists_sipShape.root")
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_sipShape.root")
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_sipShape.root")

    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/C/hists_sipShape.root")
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/D/hists_sipShape.root")
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_sipShape.root")

    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_sipShape.root", "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_sipShape.root")

    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_sipShape.root", "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_sipShape.root")

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/CD/hists_newAlg_dropDuplicates.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_newAlg_dropDuplicates.root",
    #     years    = (2022, 2022),
    #     eras     = ("CD", "EFG")
    # )

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/CD/hists_2x2e_2x2mu.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_2x2e_2x2mu.root",
    #     years    = (2022, 2022),
    #     eras     = ("CD", "EFG")
    # )

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/CD/hists_allShapeVar_24_06_25.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_allShapeVar_24_06_25.root",
    #     years    = (2022, 2022),
    #     eras     = ("CD", "EFG")
    # )
    # hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_allShapeVar_24_06_25.root")

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/CD/hists_2x2e_2x2mu_4l.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_2x2e_2x2mu_4l.root",
    #     years    = (2022, 2022),
    #     eras     = ("CD", "EFG")
    # )

    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_2x2e_2x2mu_4l.root"
    # )
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_2x2e_2x2mu_4l.root")

    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_2x2e_2x2mu.root"
    # )
    # hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_2x2e_2x2mu.root")

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2022/CD/hists_newAlg_lepPts_20_10x3.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2022/EFG/hists_newAlg_lepPts_20_10x3.root",
    #     years    = (2022, 2022),
    #     eras     = ("CD", "EFG")
    # )

    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_newAlg_lepPts_20_10x3.root"
    # )
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_newAlg_lepPts_20_10x3.root")

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2023/C/hists_2x2e_2x2mu.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/D/hists_2x2e_2x2mu.root",
    #     years    = (2023, 2023),
    #     eras     = ("C", "D")
    # )
    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_2x2e_2x2mu.root"
    # )
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_2x2e_2x2mu.root")

    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_2x2e_2x2mu.root", "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_2x2e_2x2mu.root")

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2023/C/hists_newNanos.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/D/hists_newNanos.root",
    #     years    = (2023, 2023),
    #     eras     = ("C", "D")
    # )
    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_2x2e_2x2mu.root"
    # )
    #hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_2x2e_2x2mu.root")

    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_newNanos.root"
    # )
    # hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_newNanos.root")

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2023/C/hists_newAlg_lepPts_20_10x3.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/D/hists_newAlg_lepPts_20_10x3.root",
    #     years    = (2023, 2023),
    #     eras     = ("C", "D")
    # )
    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_newAlg_lepPts_20_10x3.root"
    # )
    # hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_newAlg_lepPts_20_10x3.root")

    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_newAlg_dropDuplicates.root"
    # )

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2023/C/hists_newAlg_dropDuplicates_v2.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/D/hists_newAlg_dropDuplicates_v2.root",
    #     years    = (2023, 2023),
    #     eras     = ("C", "D")
    # )

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2023/C/hists_dataInCRs.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/D/hists_dataInCRs.root",
    #     years    = (2023, 2023),
    #     eras     = ("C", "D")
    # )

    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_newAlg_dropDuplicates_v2.root"
    # )

    # hist_manager.plot_hists(
    #     path = "/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_2x2e_2x2mu.root"
    # )

    # hist_manager.plot_hists(
    #     path = "/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_2x2e_2x2mu.root"
    # )

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2023/C/hists_MuonSF_ScaleSmear_noSysts_noPols.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/D/hists_MuonSF_ScaleSmear_noSysts_noPols.root",
    #     years    = (2023, 2023),
    #     eras     = ("C", "D")
    # )
    # hist_manager.plot_hists(
    #     path="/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_ScaleSmear_noSysts_noPols.root"
    # )
    # hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_ScaleSmear_noSysts_noPols.root")

    # hist_manager.combine_eras(
    #     infile_1 = "/eos/user/i/iehle/Analysis/histograms/2023/C/hists_MuonSF_24_06_25.root",
    #     infile_2 = "/eos/user/i/iehle/Analysis/histograms/2023/D/hists_MuonSF_24_06_25.root",
    #     years    = (2023, 2023),
    #     eras     = ("C", "D")
    # )
    # hist_manager.plot_hists(
    #     #path="/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_24_06_25.root",
    #     #path="/eos/user/i/iehle/Analysis/histograms/2022/Full/hists_2x2e_2x2mu.root",
    #     rebin=True,
    #     mass = dict(
    #         HighMass = dict(
    #             bins     = [180., 200., 220., 240., 260., 280., 300., 320., 340., 360., 380., 400., 420.],
    #             overflow = True,          
    #         ),
    #         MidMass = dict(
    #             bins     = [140., 150., 160., 170., 180.],
    #             overflow = False,          
    #         ), 
    #         LowMass = dict(
    #             bins     = [105., 115., 125., 135., 140.],
    #             overflow = False,          
    #         ), 
    #     ),
    #     Lepton_sip3d_Z2 = dict(
    #         bins = [0., 4., 20.],
    #         overflow = True,
    #     ),
    #     Z1mass = dict(
    #         bins = [81., 83., 85., 87., 89., 91., 93., 95., 97., 99., 101., 103.],
    #         overflow = False
    #     ),
    #     Z2mass = dict(
    #         bins = [81., 83., 85., 87., 89., 91., 93., 95., 97., 99., 101., 103.],
    #         overflow = False
    #     ),
    #     group = True,
    #     groups = dict(
    #         ZpX = ["DY", "TT", "WZ"]
    #     )
    # )
    # hist_manager.plot_zpx("/eos/user/i/iehle/Analysis/histograms/2023/Full/hists_MuonSF_24_06_25.root")