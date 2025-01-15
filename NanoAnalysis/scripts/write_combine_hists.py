import ROOT
import uproot as up

pol_hist_file = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/hists_shapeVars.root"

cosTheta1_path    = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnlyQCDScaleTestv3.root"
# cosTheta3_path    = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnly.root"
# cosThetaStar_path = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnly.root"
# delRapidity_path  = "/eos/user/i/iehle/Analysis/rootFiles/2022/EFG/cosTheta1_hists_polOnly.root"

#out_cosTheta1    = ROOT.TFile.Open(cosTheta1_path, "RECREATE")
# out_cosTheta3    = ROOT.TFile.Open(cosTheta3_path, "RECREATE")
# out_cosThetaStar = ROOT.TFile.Open(cosThetaStar_path, "RECREATE")
# out_delRapidity  = ROOT.TFile.Open(delRapidity_path, "RECREATE")

norm_hist = lambda hist: hist.Scale(1./hist.Integral())

with up.open(pol_hist_file) as InFile:
    with up.recreate(cosTheta1_path) as OutFile:
        for proc in ["ZLZL", "ZLZT", "ZTZT"]:
            for fs in ["fs_4e", "fs_4mu", ["fs_2e2mu", "fs_2mu2e"]]:
                if type(fs) == list:
                    for var in ["", "_LHEScaleWeightUp", "_LHEScaleWeightDown"]:
                        hist_1 = InFile[proc]["SR"]["cosTheta1"][f"{fs[0]}{var}"].to_pyroot()
                        hist_2 = InFile[proc]["SR"]["cosTheta1"][f"{fs[1]}{var}"].to_pyroot()

                        hist = hist_1 + hist_2
                        norm_hist(hist)

                        OutFile[f"{fs[0]}/{proc}{var}"] = hist
                else:
                    for var in ["", "_LHEScaleWeightUp", "_LHEScaleWeightDown"]:
                        hist = InFile[proc]["SR"]["cosTheta1"][f"{fs}{var}"].to_pyroot()
                        norm_hist(hist)

                        OutFile[f"{fs}/{proc}{var}"] = hist