import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

from pathlib import Path

import collections

import ROOT
import uproot as up
import numpy as np
from tqdm import tqdm
import pprint

import matplotlib.pyplot as plt

from NanoAnalysis.scripts.helper_functions import *
from NanoAnalysis.scripts.histWriter import HistWriter
from NanoAnalysis.scripts.histReader import HistReader
from NanoAnalysis.scripts.histPlotter import HistPlotter
from NanoAnalysis.scripts.ZpX_estimation import ZpX

class HistManager:
    def __init__(self, cfg, args):
        self.cfg  = cfg
        self.args = args

        self.year = args["year"]
        self.era  = args["era"]

        #self.samples = self._get_samples(args["year"], args["era"])

        self.regions = cfg["regions"]
        self.fstates = cfg["fstates"]
        self.props   = cfg["hist_info"].keys()
        
        self.mc_procs  = cfg["datasets"]["MC_Procs"]
        #self.procs = self.mc_procs | self.samples["Data"] | {"ZLZL": "ZLZL", "ZLZT": "ZLZT", "ZTZT": "ZTZT", "ZZ_LO": "ZZ_LO"}
        #self.procs = self.mc_procs | {"Data": "Data"} | {"ZLZL": "ZLZL", "ZLZT": "ZLZT", "ZTZT": "ZTZT", "ZZ_LO": "ZZ_LO"}
        self.procs = self.mc_procs | {"Data": "Data"}

        if args["infile"] == "":
            self.infile = os.path.join(cfg["output"]["base_dir"], "rootFiles", str(args["year"]), args["era"], "hists{}.root".format(args["tag"]))
        else:
            self.infile = args["infile"]

        self.lep_reqs = cfg["lep_reqs"]

        self.histReader = HistReader(self.cfg, self.args)

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

        return dict(
            MC = central_mc_samples,
            Data = data_sample,
            Pol = pol_samples
        )

    '''def _get_new(self, info_1, info_2, reg, prop, fs, i, cat):
        add_hists = lambda hist_1, hist_2: (hist_1[0]+hist_2[0], hist_1[1])
        add_errs  = lambda err_1, err_2: np.sqrt(err_1**2 + err_2**2)

        if cat == "Data":
            val_1 = info_1[i][reg][prop][fs][cat]
            val_2 = info_2[i][reg][prop][fs][cat]

            if i == 0:
                return add_hists(val_1, val_2)
            elif i == 1:
                return val_1 + val_2
            elif i == 2:
                return add_errs(val_1, val_2)

        new_vals = {}

        procs   = info_1[i][reg][prop][fs][cat].keys()
        procs_2 = info_2[i][reg][prop][fs][cat].keys()

        if procs == procs_2:
            vals_1 = list(info_1[i][reg][prop][fs][cat].values())
            vals_2 = list(info_2[i][reg][prop][fs][cat].values())
            
            if i==0:
                generator = (add_hists(val_1, val_2) for val_1, val_2 in zip(vals_1, vals_2))
            elif i==1:
                generator = (val_1 + val_2 for val_1, val_2 in zip(vals_1, vals_2))
            elif i==2:
                generator = (add_errs(val_1, val_2) for val_1, val_2 in zip(vals_1, vals_2))

            for proc, new_val in zip(procs, generator):
                new_vals[proc] = new_val
        
        return new_vals'''

    def _get_new(self, info_1, info_2, reg, prop, fs, cat):
        add_hists = lambda hist_1, hist_2: (hist_1[0]+hist_2[0], hist_1[1])

        if cat == "Data":
            val_1 = info_1[reg][prop][fs][cat]
            val_2 = info_2[reg][prop][fs][cat]

            return add_hists(val_1, val_2)

        new_vals = {}

        procs   = info_1[reg][prop][fs][cat].keys()
        procs_2 = info_2[reg][prop][fs][cat].keys()

        if procs == procs_2:
            vals_1 = list(info_1[reg][prop][fs][cat].values())
            vals_2 = list(info_2[reg][prop][fs][cat].values())
            
            generator = (add_hists(val_1, val_2) for val_1, val_2 in zip(vals_1, vals_2))

            for proc, new_val in zip(procs, generator):
                new_vals[proc] = new_val
        
        return new_vals

    def write_hists(self):
        self.histWriter = HistWriter(self.cfg, self.args)
        self.histWriter.write_hists(self.lep_reqs)

    def combine_processes(self, hists, counts, errors):
        for category in self.procs.keys():
            for reg in hists.keys():
                for prop in hists[reg].keys():
                    for fs in hists[reg][prop].keys():
                        
                        all_procs = list(hists[reg][prop][fs]["MC"]) + list(hists[reg][prop][fs]["Pol"]) + ["Data"]
                        if category not in all_procs:
                            procs = list(self.cfg["datasets"]["MC_Procs"][category].keys())
                            
                            hist_counts, hist_edges = hists[reg][prop][fs]["MC"][procs[0]]
                            hist_errors = [errors[reg][prop][fs]["MC"][procs[0]]]

                            del hists[reg][prop][fs]["MC"][procs[0]]
                            del counts[reg][prop][fs]["MC"][procs[0]]
                            del errors[reg][prop][fs]["MC"][procs[0]]
                            
                            for proc in procs[1:]:
                                hist = hists[reg][prop][fs]["MC"][proc]
                                err  = errors[reg][prop][fs]["MC"][proc]

                                del hists[reg][prop][fs]["MC"][proc]
                                del counts[reg][prop][fs]["MC"][proc]
                                del errors[reg][prop][fs]["MC"][proc]

                                hist_counts += hist[0]
                                hist_errors.append(err)

                            hists[reg][prop][fs]["MC"][category] = (hist_counts, hist_edges)
                            counts[reg][prop][fs]["MC"][category] = round(np.sum(hist_counts), 3)
                            errors[reg][prop][fs]["MC"][category] = np.sqrt(np.sum(np.array(hist_errors)**2, axis=0))
        
        return hists, counts, errors

    def get_hists(self, infile="", procs=[]):
        in_file = self.infile if infile=="" else infile
        if procs==[]:
            info = self.histReader.read_hists_and_counts(in_file)
        else:
            info = self.histReader.read_hists_and_counts(in_file, procs)
        return self.combine_processes(*info)

    def plot_zpx(self, infile_1, years, eras, infile_2=""):
        self.zpx = ZpX()

        hists_1, counts_1, errors_1 = self.histReader.read_hists_and_counts(infile_1)
        hists_1, counts_1, errors_1 = self.combine_processes(hists_1, counts_1, errors_1)
        zpx_info_1 = self.zpx.get_zpx(hists_1, errors_1, self.fstates)
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
                era, year = eras[0], years[0]
                self.zpx.plot_zpx(zpx_info_1, step, year, era)

    def plot_hists(self):
        self.histReader = HistReader(self.cfg, self.args)
        self.histPlotter = HistPlotter(self.cfg, self.args)


        all_hists, all_counts, all_errors = self.histReader.read_hists_and_counts(self.infile)
        all_hists, all_counts, all_errors = self.combine_processes(all_hists, all_counts, all_errors)


        for reg in self.regions:
            for prop in self.props:
                for fs in self.fstates:
                    hists  = all_hists[reg][prop][fs]
                    counts = all_counts[reg][prop][fs]
                    errors = all_errors[reg][prop][fs]

                    self.histPlotter.plotter(prop, reg, fs, hists, counts, errors)

    def plot_counts(self):
        all_hists, all_counts, all_errors = self.histReader.read_hists_and_counts(self.infile)
        all_hists, all_counts, all_errors = self.combine_processes(all_hists, all_counts, all_errors)

        categories = ("Data", "MC")
        counts = {}
        for reg in self.regions:
            counts[reg] = {}
            for fs in self.fstates:
                count_dict = all_counts[reg]["mass"][fs]
            
                n_data    = count_dict["Data"]
                n_mc   = np.array([count for count in count_dict["MC"].values()]).sum()

                counts[reg][fs] = (n_data, n_mc)
            
            x = np.arange(len(categories))  # the label locations
            width = 0.15  # the width of the bars
            multiplier = 0

            fig, ax = plt.subplots(layout='constrained')

            for fs, count in counts[reg].items():
                offset = width * multiplier
                label = fs.replace("fs_", "")
                ax.bar(x + offset, count, width, label=label, bottom = bottom)
                multiplier += 1
            
            ax.set_ylabel("Counts")
            ax.set_title("{} Counts by Final State ()".format(reg))
            ax.set_xticks(x+width, categories)
            ax.legend()
            
            fig.savefig("NotZpXCountsByFS_{}.png".format(reg))

    def combine_eras(self, infile_1, infile_2, years, eras):

        if years[0]==years[1]:
            base_dir, filename = os.path.split(infile_1)
            outdir = base_dir.replace(eras[0], "Full")
            Path(outdir).mkdir(parents=True, exist_ok=True)
            outfile = os.path.join(outdir, filename)
        else:
            # NEED TO CHANGE THIS
            outfile = "/eos/user/i/iehle/Analysis/rootFiles/hists_full.root"

        with up.open(infile_1) as Hists_1, up.open(infile_2) as Hists_2, up.recreate(outfile) as NewHists:
            for key in tqdm(Hists_1.keys()):
                if key.count("/")==3 and key in Hists_2.keys():
                    hist_1, hist_2 = Hists_1[key], Hists_2[key]
                    NewHists[key.replace(";1", "")] = hist_1.to_pyroot() + hist_2.to_pyroot()

    def insert_hists(self, infile_1, infile_2, new_file=False):
        def get_histograms(directory, path=""):
            histograms = {}
    
            # Get list of keys in the current directory
            keys = directory.GetListOfKeys()
            if not keys:
                return histograms

            for key in keys:
                obj = key.ReadObj()
                
                # If it's a TDirectory, recurse into it
                if isinstance(obj, ROOT.TDirectory):
                    subdir_name = obj.GetName()
                    subdir_path = f"{path}/{subdir_name}" if path else subdir_name
                    histograms.update(get_histograms(obj, subdir_path))
                
                # If it's a TH1 (or other histogram type), save it
                elif isinstance(obj, ROOT.TH1):
                    hist_name = obj.GetName()
                    full_hist_path = f"{path}/{hist_name}" if path else hist_name
                    histograms[full_hist_path] = obj
            
            return histograms

        if new_file:
            tag = self.args["tag"]
            outfile = infile_1.replace(".root", f"_{tag}.root")
            with up.open(infile_1) as Hists_1, up.open(infile_2) as Hists_2, up.recreate(outfile) as NewHists:
                for key in tqdm(Hists_1.keys()):
                    pass
        else:
            with ROOT.TFile(infile_1, "UPDATE") as File1, ROOT.TFile(infile_2, "READ") as File2:
                file_2_hists = get_histograms(File2)
                for key, hist in tqdm(file_2_hists.items()):
                    dir_list = key.split("/")
                    dir_name, hist_name = "/".join(dir_list[:-1]), dir_list[-1]

                    File1.cd(dir_name)
                    hist.SetDirectory(File1)
                    hist.Write(hist_name, ROOT.TObject.kOverwrite)
                    #File1.WriteObject(hist, hist_name)

    def write_zpx_hists(self, infile, year, era):
        self.zpx = ZpX()

        hists, counts, errors = self.histReader.read_hists_and_counts(infile)
        hists, counts, errors = self.combine_processes(hists, counts, errors)
        
        #zpx_info = self.zpx.get_zpx(hists, errors, self.fstates)
        yields = self.zpx.get_yields(hists, errors, self.fstates)

        with up.recreate("SS_NoSIP_HighMass_OneBin.root") as OutFile:
            for fs in yields.keys():
                for proc in yields[fs].keys():
                    count, err = yields[fs][proc]

                    hist = ROOT.TH1D(proc, proc, 1, 0., 1.)
                    hist.SetBinContent(1, count)
                    hist.SetBinError(1, err)

                    OutFile[f"{fs}/{proc}"] = hist



        #pprint.pprint(yields)

        # my_dict = {}
        # for fs in ["fs_4e", "fs_4mu", "fs_2e2mu", "fs_2mu2e"]:
        #     key = fs
        #     n_zpp_ss = max(0, zpx_info["N_ZPP_SS"][fs][0])
        #     r_os_ss  = zpx_info["r_OS_SS_MidMass"][fs][0]
            
        #     my_dict[key] = dict(
        #         N_ZPP_SS = n_zpp_ss,
        #         r_OS_SS  = r_os_ss
        #     )

        # import json
        # with open("ZpX_info_2022_EFG_20_15_15_15.json", "w") as myfile:
        #     json.dump(my_dict, myfile, indent=4)
        
        #for prop in self.props:
        # for prop in ["cosTheta1", "cosTheta3", "cosThetaStar", "delRapidity", "delPhi"]:
        #     outfile = f"{prop}_{year}_{era}_ZpX_up_down_good.root"

        #     hist_info = self.cfg["hist_info"][prop]
            
        #     hists = self.zpx.write_hists(zpx_info, hist_info, prop)

        #     with up.recreate(outfile) as OutFile:
        #         for fs, hist_dict in hists.items():
        #             OutFile[f"{fs}/ZpX"]         = hist_dict["Nominal"]
        #             OutFile[f"{fs}/ZpX_sipUp"]   = hist_dict["Up"]
        #             OutFile[f"{fs}/ZpX_sipDown"] = hist_dict["Down"]

if __name__ == "__main__":
    import yaml
    from argparse import ArgumentParser
    parser = ArgumentParser(description="")
    parser.add_argument("--reg", choices=("SR", "OS_NoSIP_HighMass", "OS_NoSIP_MidMass", "OS_NoSIP_LowMass", "SS_NoSIP_HighMass", "SS_NoSIP_MidMass", "SS_NoSIP_LowMass"), default="SR")
    parser.add_argument("--prop", default="mass")
    parser.add_argument("--year", choices=(2022, 2023), default=2022, type=int)
    parser.add_argument("--era", choices=("C", "D", "CD", "EFG", "Full"), default="Full")
    parser.add_argument("--tag", default="")
    parser.add_argument("--infile", default="")
    parser.add_argument("--lumi_tag", default=0, type=int)
    args = vars(parser.parse_args())

    with open("/afs/cern.ch/user/i/iehle/cmssw/CMSSW_14_1_6/src/ZZAnalysis/NanoAnalysis/scripts/hist_config.yaml") as config:
        cfg = yaml.safe_load(config)
  
    histManager = HistManager(cfg, args)

    # base_dir = "/eos/user/i/iehle/Analysis"
    # infile_1 = os.path.join(base_dir, "rootFiles/2022/Full/hists_Z1pt.root") # Both with no Tau pols
    # infile_2 = os.path.join(base_dir, "rootFiles/2023/Full/hists_Z1pt.root") # Both with no Tau pols
    #histManager.plot_zpx(infile_1, years=(2022, 2023), eras=("Full", "Full"), infile_2=infile_2)
    #histManager.plot_zpx(infile_1, years=(2022, 2022), eras=("Full", "Full"))
    #histManager.plot_zpx(infile_2, years=(2023, 2023), eras=("Full", "Full"))

    # base_dir = "/eos/user/i/iehle/Analysis"
    # infile_1 = os.path.join(base_dir, "rootFiles/2022/EFG/hists_Z1pt_v2.root")
    # histManager.plot_zpx(infile_1, years=(2022,2022), eras=("EFG", "EFG"))

    # base_dir = "/eos/user/i/iehle/Analysis"
    # infile = os.path.join(base_dir, "rootFiles/2022/EFG/hists_Z1pt_v2.root")
    # histManager.write_zpx_hists(infile, year=2022, era="EFG")

    # base_dir = "/eos/user/i/iehle/Analysis"
    # infile_1 = os.path.join(base_dir, "rootFiles/2022/EFG/hists_Z1pt_v2.root") # Both with no Tau pols
    # infile_2 = os.path.join(base_dir, "rootFiles/2022/EFG/hists_Z1pt_polOnly_puAndLepIDRecoVars_v2.root") # Both with no Tau pols
    # histManager.insert_hists(infile_1, infile_2)

    base_dir = "/eos/user/i/iehle/Analysis"
    infile_1 = os.path.join(base_dir, "rootFiles/2022/EFG/delPhi_hists_NormZpX_wAsimov.root") # Both with no Tau pols
    infile_2 = "/afs/cern.ch/user/i/iehle/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit/asimov_data_delPhi_test.root" # Both with no Tau pols
    histManager.insert_hists(infile_1, infile_2)

    #histManager.write_hists()

    # histManager.plot_hists()

    # base_dir = "/eos/user/i/iehle/Analysis"
    # infile_1 = os.path.join(base_dir, "rootFiles/2023/C/hists_Z1pt.root")
    # infile_2 = os.path.join(base_dir, "rootFiles/2023/D/hists_Z1pt.root")
    # histManager.combine_eras(infile_1, infile_2, years=[2023, 2023], eras=["C", "D"])
   
    # histManager.combine_eras(infile_1, infile_2, years=[2022, 2022], eras=["CD", "EFG"])

    # base_dir = "/eos/user/i/iehle/Analysis"
    # infile_1 = os.path.join(base_dir, "rootFiles/2022/Full/hists_goodSeeds.root")
    # infile_2 = os.path.join(base_dir, "rootFiles/2023/Full/hists_17Jan.root")
    # histManager.plot_zpx(infile_2, years=(2023, 2023), eras=("Full", "Full"))

    # histManager.combine_eras(infile_1, infile_2, years=[2023, 2023], eras=["C", "D"])
    # histManager.plot_hists()
    
    #histManager.write_hists()
    #histManager.plot_hists(combine_eras=True, infile_1=infile_1, infile_2=infile_2, year=2022, eras=["CD", "EFG"])
    #histManager.plot_zpx(combine_eras=True, infile_1=infile_1, infile_2=infile_2, year=2023, eras=["C", "D"])
    #histManager.plot_counts()
