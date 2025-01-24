import os
import sys

parent_dir = os.path.abspath(__file__ + 3 * "/..")
sys.path.insert(0, parent_dir)

import ROOT

from NanoAnalysis.scripts.RDFFunctions import *

import awkward as ak
import numpy as np
import uproot as up

from tqdm import tqdm

def to_raw_string(s):
    return s.encode('unicode_escape').decode('utf-8')

'''def get_samples(cfg, year, era):
    central_base = cfg["datasets"]["eos_base"]

    if era == "Full":
        era = "CD" if year==2022 else "C"

    mc_sub_path  = cfg["datasets"]["year_"+str(year)][era]["MC"]

    central_mc_procs = cfg["datasets"]["MC_Procs"]
    
    central_mc_path = os.path.join(central_base, mc_sub_path)

    central_mc_samples = {}
    for cat, procs in central_mc_procs.items():
        if isinstance(procs, dict):
            for key, val in procs.items():
                central_mc_samples[key] = os.path.join(central_mc_path, val, "ZZ4lAnalysis.root")
        else:
            central_mc_samples[cat] = os.path.join(central_mc_path, procs, "ZZ4lAnalysis.root")

    if "Data" in cfg["datasets"]["year_"+str(year)][era]:
        data_sub_path = cfg["datasets"]["year_"+str(year)][era]["Data"]
        data_sample = dict(Data = os.path.join(central_base, data_sub_path))
    else:
        data_sample = {}

    if "Pol" in cfg["datasets"]["year_"+str(year)][era]:
        pol_samples = {cat: os.path.join(central_base, pol_file) for cat, pol_file in cfg["datasets"]["year_"+str(year)][era]["Pol"].items()}
    else:
        pol_samples = {}

    return central_mc_samples | data_sample | pol_samples'''

def get_samples(cfg, year, era):
    central_base = cfg["datasets"]["eos_base"]

    era_dict = cfg["datasets"]["year_"+str(year)][era]

    central_mc_samples = {}
    if "MC" in era_dict:
        central_mc_path = os.path.join(central_base, era_dict["MC"])
        central_mc_procs = cfg["datasets"]["MC_Procs"]
        for cat, procs in central_mc_procs.items():
            if isinstance(procs, dict):
                for key, val in procs.items():
                    central_mc_samples[key] = os.path.join(central_mc_path, val, "ZZ4lAnalysis.root")
            else:
                central_mc_samples[cat] = os.path.join(central_mc_path, procs, "ZZ4lAnalysis.root")

    if "Data" in era_dict:
        data_sub_path = era_dict["Data"]
        data_sample = dict(Data = os.path.join(central_base, data_sub_path))
    else:
        data_sample = {}

    if "Pol" in era_dict:
        pol_samples = {cat: os.path.join(central_base, pol_file) for cat, pol_file in era_dict["Pol"].items()}
    else:
        pol_samples = {}

    return central_mc_samples | data_sample | pol_samples

def get_hist_info(cfg, reg, prop):
    hist_info = cfg["hist_info"][prop]
    if "mass" in prop:
        if reg == "SR" or "HighMass" in reg:
            hist_info = hist_info["SR"]
        else:
            reg = reg.split("_")
            hist_info = hist_info[reg[-1]]
    
    return hist_info

def hist_and_errors(hist, reg, prop, fstate, proc, flow=True):
    hist_arrs = hist[reg][prop][fstate][proc].to_numpy(flow=flow)
    hist_errs = np.sqrt(hist[reg][prop][fstate][proc].variances(flow=flow))

    return hist_arrs, hist_errs

def read_hists_and_counts(hist_path, **kwargs):
    hists  = {}
    counts = {}
    errors = {}

    with up.open(hist_path) as theHist:
        for reg in kwargs["regions"]:
            hists[reg]  = {}
            counts[reg] = {}
            errors[reg] = {}
            for prop in kwargs["props"]:
                hists[reg][prop] = {}
                counts[reg][prop] = {}
                errors[reg][prop] = {}
                for fstate in kwargs["fstates"]:
                    hists[reg][prop][fstate]  = dict(MC = {}, Pol = {})
                    counts[reg][prop][fstate] = dict(MC = {}, Pol = {})
                    errors[reg][prop][fstate] = dict(MC = {}, Pol = {})
                    for proc in kwargs["procs"]:
                        if proc == "Data":
                            hists[reg][prop][fstate]["Data"]  = theHist[reg+"/"+prop+"/"+fstate+"/"+proc].to_numpy(flow=True)
                            counts[reg][prop][fstate]["Data"] = round(np.sum(hists[reg][prop][fstate]["Data"][0]), 3)

                            errors[reg][prop][fstate]["Data"] = np.sqrt(theHist[reg+"/"+prop+"/"+fstate+"/"+proc].variances(flow=True))
                        elif proc in ["ZLZL", "ZLZT", "ZTZT", "ZZ_LO"]:
                                hists[reg][prop][fstate]["Pol"][proc]  = theHist[reg+"/"+prop+"/"+fstate+"/"+proc].to_numpy(flow=True)
                                counts[reg][prop][fstate]["Pol"][proc] = round(np.sum(hists[reg][prop][fstate]["Pol"][proc][0]), 3)

                                errors[reg][prop][fstate]["Pol"][proc] = np.sqrt(theHist[reg+"/"+prop+"/"+fstate+"/"+proc].variances(flow=True))
                        else:
                            if proc + ";1" in theHist[reg][prop][fstate].keys():
                                hists[reg][prop][fstate]["MC"][proc]  = theHist[reg+"/"+prop+"/"+fstate+"/"+proc].to_numpy(flow=True)
                                counts[reg][prop][fstate]["MC"][proc] = round(np.sum(hists[reg][prop][fstate]["MC"][proc][0]), 3)

                                errors[reg][prop][fstate]["MC"][proc] = np.sqrt(theHist[reg+"/"+prop+"/"+fstate+"/"+proc].variances(flow=True))
                            else:
                                sub_procs = list(kwargs['procs'][proc].keys())

                                hist_counts, hist_edges = theHist[reg][prop][fstate][sub_procs[0]].to_numpy(flow=True)

                                hist_errors = [np.sqrt(theHist[reg][prop][fstate][sub_procs[0]].variances(flow=True))]

                                for sub_proc in sub_procs[1:]:
                                    next_hist = theHist[reg][prop][fstate][sub_proc].to_numpy(flow=True)
                                    hist_counts += next_hist[0]

                                    hist_errors.append(np.sqrt(theHist[reg][prop][fstate][sub_proc].variances(flow=True)))

                                hists[reg][prop][fstate]["MC"][proc] = (hist_counts, hist_edges)
                                counts[reg][prop][fstate]["MC"][proc] = round(np.sum(hist_counts), 3)
                                errors[reg][prop][fstate]["MC"][proc] = np.sqrt(np.sum(np.array(hist_errors)**2, axis=0))
        
    return hists, counts, errors

def fill_empty_bins(hi):
    for bin in range(hi.GetNbinsX()):
        if hi.GetBinContent(bin) == 0:
            hi.SetBinContent(bin, 1e-8)
            hi.SetBinError(bin, 1e-4)
    return hi

def get_df(path, reg_idx="SR"):
    cand = "ZZCand" if reg_idx=="SR" else "ZLLCand"

    df = ROOT.RDataFrame("Events", path)
    df = df.Filter("HLT_passZZ4l")

    if "Data" in path:
        return df

    Runs = ROOT.RDataFrame("Runs", path)
    
    df = df.Define("genEventSumw", str(Runs.Sum("genEventSumw").GetValue()))
    
    # dataMCWeight stored as RVec but only ever has one entry
    cand_weight = "{}_dataMCWeight.at(0)*".format(cand)
    df = df.Define("weight", "{}overallEventWeight/genEventSumw".format(cand_weight))

    return df

def reg_df(df, reg):
    cand = "ZZCand" if reg=="SR" else "ZLLCand"
    return df.Filter("{}.at(0)==1".format(reg)).Filter("{}_Z1flav.size()==1".format(cand)).Filter("{}_Z2flav.size()==1".format(cand))

def fs_df(df, fs, reg):
    cand = "ZZCand" if reg=="SR" else "ZLLCand"
    fs_dict = dict(
        fs_4e    = (-121, -121),
        fs_4mu   = (-169, -169),
        fs_2e2mu = (-121, -169),
        fs_2mu2e = (-169, -121)
    )
    return df.Filter("{}_Z1flav.at(0)=={}".format(cand, fs_dict[fs][0])).Filter("{}_Z2flav.at(0)=={}".format(cand, fs_dict[fs][1]))

def fillLepCols(df, reg, lep_prop):
    cand = "ZZCand" if reg=="SR" else "ZLLCand"
    idx = lambda z, l: "{}_Z{}l{}Idx.at(0)".format(cand, z, l)

    lep_col = "Lepton_{}".format(lep_prop)

    def my_define(df, new_col, comm):
        if new_col not in df.GetDefinedColumnNames():
            return df.Define(new_col, comm)
        else:
            return df.Redefine(new_col, comm) 

    #df = df.Define(lep_col, "ROOT::VecOps::Concatenate(Electron_{},Muon_{})".format(lep_prop, lep_prop))

    df = my_define(df, lep_col, "ROOT::VecOps::Concatenate(Electron_{},Muon_{})".format(lep_prop, lep_prop))

    for Z in [1, 2]:
        for L in [1, 2]:
            col = lep_col + "_Z{}l{}".format(Z, L)
            #df = df.Define(col, "{}.at({})".format(lep_col, idx(Z, L)))
            df = my_define(df, col, "{}.at({})".format(lep_col, idx(Z, L)))

    def_by_z = lambda z: "ROOT::VecOps::RVec<float> {" + lep_col +"_Z{}l1, ".format(z) + lep_col + "_Z{}l2".format(z) + "}"

    #df = df.Define("{}_Z1".format(lep_col), def_by_z(1))
    #df = df.Define("{}_Z2".format(lep_col), def_by_z(2))

    df = my_define(df, "{}_Z1".format(lep_col), def_by_z(1))
    df = my_define(df, "{}_Z2".format(lep_col), def_by_z(2))

    df = df.Redefine(lep_col, "ROOT::VecOps::Concatenate({}_Z1, {}_Z2)".format(lep_col, lep_col))

    return df

def scale_hists(zlzl, zlzt, ztzt):
    zlzl_norm = zlzl.Integral()
    zlzl.Scale(1/zlzl_norm)

    zxzt       = zlzt.Clone()
    ztzt_clone = ztzt.Clone()

    zxzt.Add(ztzt_clone)

    zxzt_norm = zxzt.Integral()
    zxzt.Scale(1/zxzt_norm)

    return zlzl, zxzt

def envelope(hists, proc, prop, vary):
    up_counts  = [np.max([hist.GetBinContent(i+1) for hist in hists]) for i in range(10)]
    down_counts = [np.min([hist.GetBinContent(i+1) for hist in hists]) for i in range(10)]

    return up_counts, down_counts
        
def rms_uncertainty(hists, nbins=10):
    if nbins > -1:
        up_counts, down_counts = [], []
        sig_bins = []
        for i in range(nbins):
            all_bin_yields = np.array([hist.GetBinContent(i+1) for hist in hists])

            nom_bin_yield  = all_bin_yields[0]

            nom_bin_yields = np.full(len(hists), nom_bin_yield)
            sig_bin        = np.sqrt(np.sum((all_bin_yields - nom_bin_yields)**2))/nom_bin_yield
            sig_bins.append(sig_bin)

            '''upYield   = nom_bin_yield + sig_bin
            downYield = nom_bin_yield - sig_bin

            up_counts.append(upYield)
            down_counts.append(downYield)'''
        
        return sig_bins
        #return up_counts, down_counts
    else:
        all_bin_yields = np.array([hist.Integral() for hist in hists])

        nom_bin_yield  = all_bin_yields[0]

        nom_bin_yields = np.full(len(hists), nom_bin_yield)
        sig_bin        = np.sqrt(np.sum((all_bin_yields - nom_bin_yields)**2))/nom_bin_yield

        return sig_bin

def get_UpDownHists(proc, prop, vary, hists):

    if vary == "LHEScaleWeight":
        up_counts, down_counts = envelope(hists, proc, prop, vary)
    else:
        up_counts, down_counts = rms_uncertainty(hists)

    name = "{}_{}".format(proc, vary)
    
    upHist   = ROOT.TH1D(name+"Up", prop, 10, -1, 1)
    downHist = ROOT.TH1D(name+"Down", prop, 10, -1, 1)

    for bin in range(upHist.GetNbinsX()):
        upHist.SetBinContent(bin, up_counts[bin])
        downHist.SetBinContent(bin, down_counts[bin])

    return upHist, downHist

# def get_variation(df, proc, prop, prop_fs, vary, lumi=27.007e3):
#     hists = []

#     print("Varying ", proc)
    
#     if vary == "LHEPdfWeight":
#         for i in tqdm(range(101)):
#             df = df.Define("weight_{}".format(i), "weight*LHEPdfWeight.at({})".format(i))
#             hist = df.Histo1D((proc, prop, int(2/0.2), -1, 1), prop_fs, "weight_{}".format(i))
#             hist.Scale(lumi)
#             hists.append(hist)
    
#     elif vary == "LHEScaleWeight":
#         # Nominal value not one of the 8 options in LHEScaleWeight (for LHEPdfWeight it's i=0)
#         hist = df.Histo1D((proc, prop, int(2/0.2), -1, 1), prop_fs, "weight")
#         hist.Scale(lumi)
#         hists.append(hist)
#         for i in tqdm(range(8)):
#             if i in [2, 5]: continue # (muF, muR) = (0.5, 2) or (2, 0.5) which are unphysical
#             df = df.Define("weight_{}".format(i), "weight*LHEScaleWeight.at({})".format(i))
#             hist = df.Histo1D((proc, prop, int(2/0.2), -1, 1), prop_fs, "weight_{}".format(i))
#             hist.Scale(lumi)
#             hists.append(hist)
#     print("")

#     return hists

def get_varHists(df, prop, proc, hist_info, fs, lumi, syst_proc_info):
    prop_fs = prop + "_fs"

    df = df.Define(prop_fs, "propByFS({}, Z1flav, Z2flav, {})".format(prop, fs))

    var_hists ={}
    for var_type, var_info in syst_proc_info.items():

        var_hists[var_type] = {}
        if var_type == "Envelope":
            for var_name in var_info:
                env_hist_vars = []
                
                for var_index in var_info[var_name]:
                    new_weight      = "weight*{}.at({})".format(var_name, var_index)
                    new_weight_name = "weight_{}_{}".format(var_name, var_index)
                    df = df.Define(new_weight_name, new_weight)
                    hist = df.Histo1D((proc, prop, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), prop_fs, new_weight_name)
                    hist.Scale(lumi)
                    env_hist_vars.append(hist)
                
                hist_0 = env_hist_vars[0]
                hist_1 = env_hist_vars[1]
                hist_2 = env_hist_vars[2]

                name = "{}_{}".format(proc, var_name)
                upHist = ROOT.TH1D(name+"Up", prop, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])) 
                dnHist = ROOT.TH1D(name+"Down", prop, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"]))                  
                
                for bin_idx in range(hist_info["nbinsx"]):
                    count_0, count_1, count_2 = hist_0.GetBinContent(bin_idx), hist_1.GetBinContent(bin_idx), hist_2.GetBinContent(bin_idx)
                    up_count = max(count_0, count_1, count_2)
                    dn_count = min(count_0, count_1, count_2)
                    
                    upHist.SetBinContent(bin_idx, up_count)
                    dnHist.SetBinContent(bin_idx, dn_count)
                
                var_hists[var_type][var_name] = {"Down": dnHist, "Up": upHist}

        elif var_type == "RMS":
            pass

            # hist = df.Histo1D((proc, prop, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), prop_fs, "weight")
            # hist.Scale(lumi)

    return var_hists

def get_hist(df, prop, proc, hist_info, fs, lumi):
    prop_fs = prop + "_fs"

    df = df.Define(prop_fs, "propByFS({}, Z1flav, Z2flav, {})".format(prop, fs))

    if proc != "Data":
        hist = df.Histo1D((proc, prop, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), prop_fs, "weight")
        hist.Scale(lumi)
    else:
        hist = df.Histo1D((proc, prop, int(hist_info["nbinsx"]), float(hist_info["xlow"]), float(hist_info["xhigh"])), prop_fs)

    return hist

def get_hists(df, prop, proc, fs, vary, lumi=27.007e3):
    fstates = {
        "fs_4e": 14641,
        "fs_4mu": 28561,
        "fs_2e2mu": 20449
    }

    prop_fs = "{}_{}".format(prop, fs)
    df = df.Define(prop_fs, "propByFS({}, Z1flav, Z2flav, {})".format(prop, fstates[fs]))

    hists = get_variation(df, proc, prop, prop_fs, vary)

    return hists
    '''hist_nom = hists[0]
    hist_up, hist_down = get_UpDownHists(proc, prop, vary, hists)

    return dict(
        nominal = hist_nom,
        up      = hist_up,
        down    = hist_down,
    )'''

def def_cols(df, reg_idx, props, proc):

    cand = "ZZCand" if reg_idx == "SR" else "ZLLCand"

    zz_props  = [prop for prop in props if "Lepton" not in prop]
    lep_props = [prop.split("_")[-1] for prop in props if "Lepton" in prop]

    def reg_prop(prop):
        if prop not in df.GetColumnNames():
            return "applyRegFilt({}_{}, {})".format(cand, prop, reg_idx)
        else:
            return "applyRegFilt({}, {})".format(prop, reg_idx)

    lep_from_cand = lambda prop, z: "lepFromCand(Electron_{}, Muon_{}, {}_{}l1Idx, {}_{}l2Idx, {})".format(prop, prop, cand, z, cand, z, reg_idx)

    for prop in zz_props:
        tag = "" if prop not in df.GetColumnNames() else "_SR"
        df = df.Define(prop+tag, reg_prop(prop))

    for prop in lep_props:
        df = df.Define("Lepton_Z1_{}".format(prop), lep_from_cand(prop, "Z1"))
        df = df.Define("Lepton_Z2_{}".format(prop), lep_from_cand(prop, "Z2"))
        df = df.Define("Lepton_{}".format(prop), "ROOT::VecOps::Concatenate(Lepton_Z1_{}, Lepton_Z2_{})".format(prop, prop))        

    return df

def get_arrs(path, props, reg="SR"):
    df = get_df(path)
    df = def_cols(df, reg, props)

    arrs = ak.from_rdataframe(df, columns=props)

    good_idx = ak.flatten(np.argwhere(ak.num(arrs[props[0]]) == 1))

    good_arrs = {}
    for prop in arrs.fields:
        if prop=="weight":
            good_arrs[prop] = arrs[prop][good_idx]
        else:
            good_arrs[prop] = ak.flatten(arrs[prop][good_idx])

    return good_arrs

def get_bin_err(prop_arr, weight_arr, bins):
        
    digits = np.digitize(prop_arr, bins)
    
    nbins = len(bins)
    sum_w2 = np.zeros([nbins], dtype=np.float32)
    for bin in range(nbins):
        bin_weight = weight_arr[np.where(digits==bin)[0]]
        sum_w2[bin] = np.sum(bin_weight**2)
    
    return np.sqrt(sum_w2)

def count_and_err(prop_arr, weight_arr, bins):

    counts, _ = np.histogram(prop_arr, bins=bins, weights=weight_arr)
    bin_errs  = get_bin_err(prop_arr, weight_arr, bins)

    tot = np.sum(counts)
    err = np.sqrt(np.sum(bin_errs**2))

    return tot, err

def ratio_and_err(x, x_err, y, y_err):
    ratio = x/y
    err   = ratio*np.sqrt((x_err/x)**2 + (y_err/y)**2)

    return ratio, err


if __name__=="__main__":
    cms_store = "root://cms-xrd-global.cern.ch/"
    zz_nano   = "/store/mc/Run3Summer22EENanoAODv12/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/50000/dea56a0a-462c-4690-953c-c7db96dd3ab5.root"

    df = get_df("pol_skims/wCosThetaStar/ZZ_hadd_Skim.root")