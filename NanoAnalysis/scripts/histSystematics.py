import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT
import uproot as up
import numpy as np
from tqdm import tqdm

class HistSystematics:
    #def __init__(self, hist_info, lumi, prop, column, pdf_prescription="hessian"):
    def __init__(self, hist_info, lumi, column, weight, pdf_prescription="hessian"):
        self.hist_info        = hist_info

        self.lumi             = lumi

        # self.prop             = prop
        # self.column           = column

        self.column           = column
        self.weight           = weight

        self.pdf_prescription = pdf_prescription
        
        self.th1_model        = ROOT.RDF.TH1DModel("", "", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))

    def get_vars(self, df, var):
        if var == "LHEScaleWeight":
            return self.lhe_scale_weight(df)
        elif var == "LHEPdfWeight":
            return self.lhe_pdf_weight(df)
        elif var == "puWeight":
            return self.pu_weight(df)
        elif var == "lepIDRec":
            return self.lepIDReco(df)
        else:
            raise ValueError("{} is not a supported varation. Options are: [LHEScaleWeight, LHEPdfWeight, puWeight, lepIDRec].".format(var))

    def lhe_scale_weight(self, df):
        # Indices where ratio between mu_f, mu_r variations is less than 4
        phys_variations = [0, 1, 3, 4, 5, 7, 8]

        nominal      = df.Vary(f"{self.weight}", f"{self.weight}*LHEScaleWeight", [f"qcd_{i}" for i in range(9)]).Histo1D(self.th1_model, self.column, f"{self.weight}")
        hists_varied = ROOT.RDF.Experimental.VariationsFor(nominal)

        up_hist = ROOT.TH1D(self.column+"_Up", self.column+"_Up", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))
        dn_hist = ROOT.TH1D(self.column+"_Down", self.column+"_Down", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))

        for bin_idx in range(int(self.hist_info["nbinsx"])):
            max_bin_count = np.max([hists_varied[f"{self.weight}:qcd_{i}"].GetBinContent(bin_idx+1) for i in phys_variations])
            min_bin_count = np.min([hists_varied[f"{self.weight}:qcd_{i}"].GetBinContent(bin_idx+1) for i in phys_variations])

            up_hist.SetBinContent(bin_idx+1, max_bin_count)
            dn_hist.SetBinContent(bin_idx+1, min_bin_count)

        up_hist.Scale(self.lumi)
        dn_hist.Scale(self.lumi)

        return up_hist, dn_hist

    def lhe_pdf_weight(self, df):
        nominal      = df.Vary(f"{self.weight}",f"{self.weight}*LHEPdfWeight",[f"pdf_{i}" for i in range(103)]).Histo1D(self.th1_model, self.column, f"{self.weight}")
        hists_varied = ROOT.RDF.Experimental.VariationsFor(nominal) 

        up_hist = ROOT.TH1D(self.column+"_Up", self.column+"_Up", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))
        dn_hist = ROOT.TH1D(self.column+"_Down", self.column+"_Down", int(self.hist_info["nbinsx"]), float(self.hist_info["xlow"]), float(self.hist_info["xhigh"]))

        for bin_idx in range(int(self.hist_info["nbinsx"])):
            nom_bin_count = hists_varied["nominal"].GetBinContent(bin_idx+1) * self.lumi

            # Hessian representation
            if self.pdf_prescription == "hessian":
                pdf_nom  = np.full(102, nom_bin_count)
                pdf_vars = np.array([hists_varied[f"{self.weight}:pdf_{i}"].GetBinContent(bin_idx+1)*self.lumi for i in range(1,103)])
                delta = np.sqrt(np.sum(np.square(pdf_nom-pdf_vars)))

                up_hist.SetBinContent(bin_idx+1, nom_bin_count + delta)
                dn_hist.SetBinContent(bin_idx+1, nom_bin_count - delta)

            # This is the MC uncertainty prescription (sec. 6.3.1)
            elif self.pdf_prescription == "mc":
                bin_std       = np.std([hists_varied[f"{self.weight}:pdf_{i}"].GetBinContent(bin_idx+1)*self.lumi for i in range(103)])

                up_hist.SetBinContent(bin_idx+1, nom_bin_count + bin_std)
                dn_hist.SetBinContent(bin_idx+1, nom_bin_count - bin_std)

        return up_hist, dn_hist

    def pu_weight(self, df):
        df = df.Define(f"{self.weight}_pu_up", f"{self.weight}*(puWeightUp/puWeight)")
        df = df.Define(f"{self.weight}_pu_dn", f"{self.weight}*(puWeightDn/puWeight)")

        up_hist = df.Histo1D(self.th1_model, self.column, f"{self.weight}_pu_up").GetValue()
        dn_hist = df.Histo1D(self.th1_model, self.column, f"{self.weight}_pu_dn").GetValue()

        up_hist.Scale(self.lumi)
        dn_hist.Scale(self.lumi)

        return up_hist, dn_hist

    def lepIDReco(self, df):
        df = df.Define("weight_lepIDRecUp", "(ZZCand_dataMCWeight.at(0) + ZZCand_lepSF_err.at(0))*(overallEventWeight/genEventSumw)")
        df = df.Define("weight_lepIDRecDn", "(ZZCand_dataMCWeight.at(0) - ZZCand_lepSF_err.at(0))*(overallEventWeight/genEventSumw)")

        up_hist = df.Histo1D(self.th1_model, self.column, "weight_lepIDRecUp").GetValue()
        dn_hist = df.Histo1D(self.th1_model, self.column, "weight_lepIDRecDn").GetValue()

        up_hist.Scale(self.lumi)
        dn_hist.Scale(self.lumi)

        return up_hist, dn_hist