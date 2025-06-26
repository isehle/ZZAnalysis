import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT
import uproot as up
import numpy as np
from tqdm import tqdm

import yaml

cfg_path = os.path.join(parent_dir, "NanoAnalysis/scripts/systs_cfg.yaml")
with open(cfg_path) as config:
    cfg = yaml.safe_load(config)

class HistSystematics:
    def __init__(self, hist_info, lumi, column, weight, reg_idx, pdf_prescription="hessian"):
        self.cfg              = cfg

        self.hist_info        = hist_info

        self.lumi             = lumi

        self.column           = column
        self.weight           = weight
        self.reg_idx          = reg_idx

        self.pdf_prescription = pdf_prescription

        self.th1_model        = lambda var: ROOT.RDF.TH1DModel(f"{self.column}_{var}", self.column, int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))

        self.variations       = {}

    def get_vars(self, df, var):
        if var == "LHEScaleWeight":
            #return self.lhe_scale_weight(df)
            return self.qcd_vars(df)
        elif var == "LHEPdfWeight":
            #return self.lhe_pdf_weight(df)
            return self.pdf_vars(df)
        elif var == "puWeight":
            return self.pu_weight(df)
        elif var == "lepIDRec":
            return self.lepIDReco(df)
        elif var == "LHE":
            return self.lhe_vars(df)
        else:
            raise ValueError("{} is not a supported varation. Options are: [LHEScaleWeight, LHEPdfWeight, puWeight, lepIDRec].".format(var))

    def lhe_scale_weight(self, df):
        # Indices where ratio between mu_f, mu_r variations is less than 4
        phys_variations = [0, 1, 3, 4, 5, 7, 8]

        nominal      = df.Vary(f"{self.weight}", f"{self.weight}*LHEScaleWeight", [f"qcd_{i}" for i in range(9)], "QCDScale").Histo1D(self.th1_model("LHEScaleWeight"), self.column, f"{self.weight}")
        hists_varied = ROOT.RDF.Experimental.VariationsFor(nominal)

        up_hist = ROOT.TH1D(self.column+"_LHEScaleWeightUp", self.column+"_Up", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))
        dn_hist = ROOT.TH1D(self.column+"_LHEScaleWeightDown", self.column+"_Down", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))

        # Triggers the event loop!
        for bin_idx in range(int(self.hist_info["nbinsx"])):
            bin_counts    = [hists_varied[f"QCDScale:qcd_{i}"].GetBinContent(bin_idx+1) for i in phys_variations]
            
            max_bin_count = np.max(bin_counts)
            min_bin_count = np.min(bin_counts)

            up_hist.SetBinContent(bin_idx+1, max_bin_count)
            dn_hist.SetBinContent(bin_idx+1, min_bin_count)

        return up_hist, dn_hist

    def lhe_pdf_weight(self, df):
        nominal      = df.Vary(f"{self.weight}",f"{self.weight}*LHEPdfWeight",[f"pdf_{i}" for i in range(103)], "PDF").Histo1D(self.th1_model("LHEPdfWeight"), self.column, f"{self.weight}")
        hists_varied = ROOT.RDF.Experimental.VariationsFor(nominal) 

        up_hist = ROOT.TH1D(self.column+"_LHEPdfWeightUp", self.column+"_Up", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))
        dn_hist = ROOT.TH1D(self.column+"_LHEPdfWeightDown", self.column+"_Down", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))

        # Triggers the event loop!
        for bin_idx in range(int(self.hist_info["nbinsx"])):
            nom_bin_count = hists_varied["nominal"].GetBinContent(bin_idx+1)

            # Hessian representation
            if self.pdf_prescription == "hessian":
                pdf_nom  = np.full(102, nom_bin_count)
                #pdf_vars = np.array([hists_varied[f"{self.weight}:pdf_{i}"].GetBinContent(bin_idx+1)*self.lumi for i in range(1,103)])
                pdf_vars = np.array([hists_varied[f"PDF:pdf_{i}"].GetBinContent(bin_idx+1) for i in range(1,103)])
                delta = np.sqrt(np.sum(np.square(pdf_nom-pdf_vars)))

                up_hist.SetBinContent(bin_idx+1, nom_bin_count + delta)
                dn_hist.SetBinContent(bin_idx+1, nom_bin_count - delta)

            # This is the MC uncertainty prescription (sec. 6.3.1)
            elif self.pdf_prescription == "mc":
                #bin_std       = np.std([hists_varied[f"{self.weight}:pdf_{i}"].GetBinContent(bin_idx+1)*self.lumi for i in range(103)])
                bin_std       = np.std([hists_varied[f"{self.weight}:pdf_{i}"].GetBinContent(bin_idx+1) for i in range(103)])

                up_hist.SetBinContent(bin_idx+1, nom_bin_count + bin_std)
                dn_hist.SetBinContent(bin_idx+1, nom_bin_count - bin_std)

        return up_hist, dn_hist

    # def qcd_up_down(self, hists_varied):
    #     # Indices where ratio between mu_f, mu_r variations is less than 4
    #     phys_variations = [0, 1, 3, 4, 5, 7, 8]

    #     up_hist = ROOT.TH1D(self.column+"_LHEScaleWeightUp", self.column+"_Up", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))
    #     dn_hist = ROOT.TH1D(self.column+"_LHEScaleWeightDown", self.column+"_Down", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))

    #     # Triggers the event loop!
    #     for bin_idx in range(int(self.hist_info["nbinsx"])):
    #         bin_counts    = [hists_varied[f"QCDScale:qcd_{i}"].GetBinContent(bin_idx+1) for i in phys_variations]
            
    #         max_bin_count = np.max(bin_counts)
    #         min_bin_count = np.min(bin_counts)

    #         up_hist.SetBinContent(bin_idx+1, max_bin_count)
    #         dn_hist.SetBinContent(bin_idx+1, min_bin_count)

    #     return up_hist, dn_hist

    def qcd_up_down(self, hists_varied, hist_info):
        # Indices where ratio between mu_f, mu_r variations is less than 4
        phys_variations = [0, 1, 3, 4, 5, 7, 8]

        up_hist = ROOT.TH1D("", "", int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"]))
        dn_hist = ROOT.TH1D("", "", int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"]))

        # Triggers the event loop!
        for bin_idx in range(int(hist_info["nbinsx"])):
            bin_counts    = [hists_varied[f"QCDScale:qcd_{i}"].GetBinContent(bin_idx+1) for i in phys_variations]
            
            max_bin_count = np.max(bin_counts)
            min_bin_count = np.min(bin_counts)

            up_hist.SetBinContent(bin_idx+1, max_bin_count)
            dn_hist.SetBinContent(bin_idx+1, min_bin_count)

        return up_hist, dn_hist

    def pdf_up_down(self, hists_varied, hist_info):
        up_hist = ROOT.TH1D("", "", int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"]))
        dn_hist = ROOT.TH1D("", "", int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"]))

        # Triggers the event loop!
        for bin_idx in range(int(hist_info["nbinsx"])):
            nom_bin_count = hists_varied["nominal"].GetBinContent(bin_idx+1)

            # Hessian representation
            if self.pdf_prescription == "hessian":
                pdf_nom  = np.full(102, nom_bin_count)
                #pdf_vars = np.array([hists_varied[f"{self.weight}:pdf_{i}"].GetBinContent(bin_idx+1)*self.lumi for i in range(1,103)])
                pdf_vars = np.array([hists_varied[f"PDF:pdf_{i}"].GetBinContent(bin_idx+1) for i in range(1,103)])
                delta = np.sqrt(np.sum(np.square(pdf_nom-pdf_vars)))

                up_hist.SetBinContent(bin_idx+1, nom_bin_count + delta)
                dn_hist.SetBinContent(bin_idx+1, nom_bin_count - delta)

            # This is the MC uncertainty prescription (sec. 6.3.1)
            elif self.pdf_prescription == "mc":
                #bin_std       = np.std([hists_varied[f"{self.weight}:pdf_{i}"].GetBinContent(bin_idx+1)*self.lumi for i in range(103)])
                bin_std       = np.std([hists_varied[f"PDF:pdf_{i}"].GetBinContent(bin_idx+1) for i in range(103)])

                up_hist.SetBinContent(bin_idx+1, nom_bin_count + bin_std)
                dn_hist.SetBinContent(bin_idx+1, nom_bin_count - bin_std)

        return up_hist, dn_hist

    def qcd_vars(self, df):
        qcd_nom = df.Vary(f"{self.weight}", f"{self.weight}*LHEScaleWeight", [f"qcd_{i}" for i in range(9)], "QCDScale").Histo1D(self.th1_model("LHEScaleWeight"), self.column, f"{self.weight}")
        return  ROOT.RDF.Experimental.VariationsFor(qcd_nom)

    def pdf_vars(self, df):
        pdf_nom = df.Vary(f"{self.weight}",f"{self.weight}*LHEPdfWeight",[f"pdf_{i}" for i in range(103)], "PDF").Histo1D(self.th1_model("LHEPdfWeight"), self.column, f"{self.weight}")
        return ROOT.RDF.Experimental.VariationsFor(pdf_nom)

    def lhe_vars(self, df):
        # Must book all variations before triggering the event loop, unbooked variations will always return zero otherwise

        # QCD Scale
        qcd_nom    = df.Vary(f"{self.weight}", f"{self.weight}*LHEScaleWeight", [f"qcd_{i}" for i in range(9)], "QCDScale").Histo1D(self.th1_model("LHEScaleWeight"), self.column, f"{self.weight}")
        qcd_varies = ROOT.RDF.Experimental.VariationsFor(qcd_nom)

        # PDF
        pdf_nom    = df.Vary(f"{self.weight}",f"{self.weight}*LHEPdfWeight",[f"pdf_{i}" for i in range(103)], "PDF").Histo1D(self.th1_model("LHEPdfWeight"), self.column, f"{self.weight}")
        pdf_varies = ROOT.RDF.Experimental.VariationsFor(pdf_nom)

        return qcd_varies, pdf_varies

        # qcd_up, qcd_down = self.qcd_up_down(qcd_varies)
        # pdf_up, pdf_down = self.pdf_up_down(pdf_varies)

        # return dict(
        #     QCD = (qcd_up, qcd_down),
        #     PDF = (pdf_up, pdf_down)
        # )

    def pu_weight(self, df):
        df = df.Define(f"{self.weight}_pu_up", f"{self.weight}*(puWeightUp/puWeight)")
        df = df.Define(f"{self.weight}_pu_dn", f"{self.weight}*(puWeightDn/puWeight)")

        up_hist = df.Histo1D(self.th1_model("puWeight"), self.column, f"{self.weight}_pu_up")
        dn_hist = df.Histo1D(self.th1_model("puWeight"), self.column, f"{self.weight}_pu_dn")

        return up_hist, dn_hist

    def lepIDReco(self, df):
        df = df.Define(f"{self.weight}_lepIDRecUp", f"(ZZCand_dataMCWeight[{self.reg_idx}] + ZZCand_dataMCWeight_err[{self.reg_idx}])*(overallEventWeight/genEventSumw)")
        df = df.Define(f"{self.weight}_lepIDRecDn", f"(ZZCand_dataMCWeight[{self.reg_idx}] - ZZCand_dataMCWeight_err[{self.reg_idx}])*(overallEventWeight/genEventSumw)")

        up_hist = df.Histo1D(self.th1_model("lepIDRec"), self.column, f"{self.weight}_lepIDRecUp")
        dn_hist = df.Histo1D(self.th1_model("lepIDRec"), self.column, f"{self.weight}_lepIDRecDn")

        return up_hist, dn_hist