#!/usr/bin/env python3
###
# Example for running the analysis locally, after customizing variables.
# Run with: 
# python runLocal.py
###
from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertAfter

# Check that the checkout recipe has been properly updated 
from ZZAnalysis.AnalysisStep.validateCheckout import validateCheckout 
if not validateCheckout() :
    exit(1)

#SampleToRun = "MCsync_Rereco"
#SampleToRun = "MCsync_UL"
#SampleToRun = "Data2022"
#SampleToRun = "MC2022"

#SampleToRun = "ZLZLTo4L_2022CD"
#SampleToRun = "ZLZTTo4L_2022CD"
#SampleToRun = "ZTZTTo4L_2022CD"
#SampleToRun = "ZZTo4L_2022CD_LO"

#SampleToRun = "ZZTo4L_2022EE" # Centrally produced Powheg NLO
#SampleToRun = "ZZTo4L_2022EE_MG" # Privately produced MadGraph LO
#SampleToRun = "ZLZLTo4L_2022EE"
#SampleToRun = "ZLZTTo4L_2022EE"
#SampleToRun = "ZTZTTo4L_2022EE"

#SampleToRun = "ZLZLTo4L_2023C"
#SampleToRun = "ZLZTTo4L_2023C"
#SampleToRun = "ZTZTTo4L_2023C"
SampleToRun = "ZZTo4L_2023C_LO"

#SampleToRun = "ZLZLTo4L_2023D"
#SampleToRun = "ZLZTTo4L_2023D"
#SampleToRun = "ZTZTTo4L_2023D"
#SampleToRun = "ZZTo4L_2023D_LO"

#SampleToRun = "ggZLZLTo4L_2022EE" # Privately produced gg-->ZZ (box diagram) MadGraph
#SampleToRun = "ggZLZTTo4L_2022EE"
#SampleToRun = "ggZTZTTo4L_2022EE"

#SampleToRun = "ggZZ_2022EE"

### Customize processing variables
#setConf("runMELA", False)
#setConf("bestCandByMELA", False)
#setConf("APPLYMUCORR", False)
#setConf("APPLYELECORR", False)


## Force filling K factors and weights (default: all off)
#setConf("APPLY_K_NNLOQCD_ZZGG", 1) # 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
#setConf("APPLY_K_NNLOQCD_ZZQQB", True)
#setConf("APPLY_K_NNLOEW_ZZQQB", True)
#setConf("APPLY_QCD_GGF_UNCERT", True)

setConf("PROCESS_CR", True)
setConf("PROCESS_ZL", True)
setConf("DEBUG", False)
setConf("SYNCMODE", True) # Force muon resolution correction with fixed +1 sigma smearing
#setConf("ADD_ALLEVENTS", True) # Add extra tree of gen info for all events

json = None #replace this if needed

################################################################################
if SampleToRun == "Data2022" :
    # 2022 data sample from /MuonEG/Run2022D-PromptNanoAODv10_v1-v1/NANOAOD
    setConf("IsMC", False)
    setConf("LEPTON_SETUP", 2022)
    setConf("PD", "any")
    setConf("SAMPLENAME", "test")
    setConf("TRIGPASSTHROUGH", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/data/Run2022D/MuonEG/NANOAOD/PromptNanoAODv10_v2-v1/50000/68f42f42-3274-46ec-b23d-bfadc13012c2.root",
        ])


################################################################################
elif SampleToRun == "ggh125_UL" : ### 2018 UL test sample
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2018)
    setConf("DATA_TAG", "UL")
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("fileNames",[
        "/store/mc/RunIISummer20UL18NanoAODv2/WplusH_HToZZTo4L_M125_TuneCP5_13TeV_powheg2-minlo-HWJ_JHUGenV7011_pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/3B6A5CB5-2B7C-924D-85B4-FC3B0C1F4909.root",
        ])

################################################################################
elif SampleToRun == "MCsync_UL" :
    # Custom-reprocessed Rereco nanoAOD file with updated FSR and electron MVA,
    # no packing for genparticle p3; 26000 events
    # corresponding to:/store/mc/RunIISummer20UL17MiniAODv2/GluGluHToZZTo4L_M125_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/130000/3E4E8D55-3993-2B43-AF3B-7AB45BBE0BDA.root
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("LEPTON_SETUP", 2017)
    setConf("NANOVERSION", 10) # variable defined as per nanoAOD v10 (notably electron_mvaHZZIso)
    setConf("DATA_TAG", "UL")
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR.root"])
#    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_2017UL_fixedFSR_nopacking.root"]) # with no packing of muon eta, phi, mass


################################################################################
elif SampleToRun == "MCsync_Rereco" :
     # Custom-reprocessed Rereco nanoAOD file with updated FSR,
     # corresponding to:/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root
    setConf("APPLYMUCORR", True)
    setConf("SAMPLENAME", "ggH125")
    setConf("XSEC", 48.58*0.0002745)
    setConf("NANOVERSION", 9)
    setConf("store","")
    setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_fixedFSR.root"])


################################################################################
elif SampleToRun == "MC2022" :
    # 2022 MC sample
    setConf("SAMPLENAME", "ggH125")
    setConf("DATA_TAG", "post_EE")
    setConf("XSEC", 52.23*0.0002745)
    setConf("LEPTON_SETUP", 2022)
    setConf("IsMC", True)
    setConf("store","root://cms-xrd-global.cern.ch/")
    setConf("APPLY_QCD_GGF_UNCERT", True) # for ggH
    setConf("fileNames",[
        "/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2540000/25c8f5ff-9de0-4a0c-9e2f-757332ad392f.root",
#        "/store/mc/Run3Summer22EENanoAODv12/GluGluHtoZZto4L_M-125_TuneCP5_13p6TeV_powheg2-JHUGenV752-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/2530000/8f306f2b-1284-41b8-a98f-744267f64b9c.root",
        ])
#    json = {"1": [[1245, 1245],[1306, 1306],[1410, 1410],[1692, 1692],[1903, 1903],[1910, 1910],[1915, 1915],[1927, 1927],[1939, 1939],[1940, 1940],[1944, 1944],[1945, 1945],[1956, 1956],[1960, 1960],[1965, 1965],[1967, 1967],[1968, 1968],[1969, 1969],[2104, 2104]]}

elif SampleToRun == "ZZTo4L_2022EE":
    setConf("SAMPLENAME", "ZZTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", True)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 1.39)
    setConf("LEPTON_SETUP", 2022)
    setConf("DATA_TAG", "post_EE")
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("runMELA", False)
    setConf("bestCandByMELA", False)
    setConf("store", "root://cms-xrd-global.cern.ch/")
    setConf("fileNames", [
        "/store/mc/Run3Summer22EENanoAODv12/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/50000/dea56a0a-462c-4690-953c-c7db96dd3ab5.root",
    ])

elif SampleToRun == "ZLZLTo4L_2022CD":
    setConf("SAMPLENAME", "ZLZLTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 0.003963)
    setConf("LEPTON_SETUP", 2022)
    setConf("DATA_TAG", "pre_EE")
    setConf("NANOVERSION", 11)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZLZLTo4l_5f_Summer22EraCD_onlyCERN/crab_qqzlzl-4l-v3/241216_101620/0000/qqZLZL_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZLZTTo4L_2022CD":
    setConf("SAMPLENAME", "ZLZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 0.01551)
    setConf("LEPTON_SETUP", 2022)
    setConf("DATA_TAG", "pre_EE")
    setConf("NANOVERSION", 11)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZLZTTo4l_5f_Summer22EraCD_onlyCERN/crab_qqzlzt-4l-v3/241216_161621/0000/qqZLZT_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZTZTTo4L_2022CD":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 0.045113)
    setConf("LEPTON_SETUP", 2022)
    setConf("DATA_TAG", "pre_EE")
    setConf("NANOVERSION", 11)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZTZTTo4l_5f_Summer22EraCD_onlyCERN/crab_qqztzt-4l-v3/241216_162933/0000/qqZTZT_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZZTo4L_2022CD_LO":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 0.06644)
    setConf("LEPTON_SETUP", 2022)
    setConf("DATA_TAG", "pre_EE")
    setConf("NANOVERSION", 11)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZZTo4l_5f_Summer22EraCD_onlyCERN/crab_qqzz-4l-v3/241216_163742/0000/qqZZ_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZLZLTo4L_2022EE":
    setConf("SAMPLENAME", "ZLZLTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 0.003963)
    setConf("LEPTON_SETUP", 2022)
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZLZLTo4l_5f_Summer22EraEFG_onlyCERN/crab_qqzlzl-4l-cern/241214_150850/0000/qqZLZL_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZLZTTo4L_2022EE":
    setConf("SAMPLENAME", "ZLZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 0.01551)
    setConf("LEPTON_SETUP", 2022)
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZLZTTo4l_5f_Summer22EraEFG_onlyCERN/crab_qqzlzt-4l-v3/241215_112122/0000/qqZLZT_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZTZTTo4L_2022EE":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 0.045113)
    setConf("LEPTON_SETUP", 2022)
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZTZTTo4l_5f_Summer22EraEFG_onlyCERN/crab_qqztzt-4l-v3/241215_112239/0000/qqZTZT_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZZTo4L_2022EE_MG":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("XSEC", 0.06644)
    setConf("LEPTON_SETUP", 2022)
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZZTo4l_5f_Summer22EraEFG_onlyCERN/crab_qqzz-4l-v3/241215_112331/0000/qqZZ_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ggZLZLTo4L_2022EE":
    setConf("SAMPLENAME", "ggZLZLTo4L")
    setConf("APPLY_K_NNLOQCD_ZZGG", 2)
    setConf("XSEC", 0.001274)
    setConf("LEPTON_SETUP", 2022)
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/ggZLZL_4l_5f.root"
    ])

# elif SampleToRun == "ggZLZTTo4L_2022EE":
#     setConf("SAMPLENAME", "ZLZTTo4L")
#     setConf("APPLY_K_NNLOQCD_ZZGG", 2)
#     #setConf("XSEC", 0.0281)
#     setConf("XSEC", 0.01551)
#     setConf("LEPTON_SETUP", 2022)
#     setConf("NANOVERSION", 12)
#     setConf("IsMC", True)
#     setConf("store", "root://eos.grif.fr/")
#     setConf("fileNames", [
#         #"/eos/grif/cms/llr/store/user/iehle/MyZLZTTo4l/crab_zlzt-4l-tarBallTest_v1/240430_100126/0000/ZLZT_hadd.root"
#         "/eos/grif/cms/llr/store/user/iehle/ggZLZT_4l_5f.root"
#     ])

elif SampleToRun == "ggZTZTTo4L_2022EE":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZGG", 2)
    setConf("XSEC", 0.006941)
    setConf("LEPTON_SETUP", 2022)
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        #"/eos/grif/cms/llr/store/user/iehle/MyZTZTTo4l/crab_ztzt-4l-tarBallTest_v1/240430_154317/0000/ZTZT_hadd.root"
        "/eos/grif/cms/llr/store/user/iehle/ggZTZT_4l_5f.root"
    ])

elif SampleToRun == "ggZZ_2022EE":
    setConf("SAMPLENAME", "ggTo4e")
    setConf("APPLY_K_NNLOQCD_ZZGG", 1)
    setConf("XSEC", 0.00305851)
    setConf("LEPTON_SETUP", 2022)
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://cms-xrd-global.cern.ch/")
    setConf("fileNames", [
        "/store/mc/Run3Summer22EENanoAODv12/GluGlutoContinto2Zto4E_TuneCP5_13p6TeV_mcfm-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v2/40000/3606c6cd-0b9c-4413-8338-4ac70c369c3a.root"
    ])

elif SampleToRun == "ZLZLTo4L_2023C":
    setConf("SAMPLENAME", "ZLZLTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLYELECORR", False) # TEMPORARY UNTIL THEY ARE READY
    setConf("APPLYJETCORR", False)
    setConf("XSEC", 0.003939)
    setConf("LEPTON_SETUP", 2023)
    setConf("NANOVERSION", 12)
    setConf("DATA_TAG", "pre_BPix")
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZLZLTo4l_5f_Summer23EraC_onlyCERN_good/crab_qqzlzl-4l-v5/241220_112409/0000/qqZLZL_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZLZTTo4L_2023C":
    setConf("SAMPLENAME", "ZLZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLYELECORR", False) # TEMPORARY UNTIL THEY ARE READY
    setConf("APPLYJETCORR", False)
    setConf("XSEC", 0.01592) # Mean value of LHEWeight_originalXWGTUP
    setConf("LEPTON_SETUP", 2023)
    setConf("DATA_TAG", "pre_BPix")
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZLZTTo4l_5f_Summer23EraC_onlyCERN_good/crab_qqzlzt-4l-v2/241220_112440/0000/qqZLZT_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZTZTTo4L_2023C":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLYELECORR", False) # TEMPORARY UNTIL THEY ARE READY
    setConf("APPLYJETCORR", False)
    #setConf("XSEC", 0.045113)
    setConf("XSEC", 0.04569)
    setConf("LEPTON_SETUP", 2023)
    setConf("DATA_TAG", "pre_BPix")
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZTZTTo4l_5f_Summer23EraC_onlyCERN_good/crab_qqztzt-4l-v2/241220_112509/0000/qqZTZT_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZZTo4L_2023C_LO":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLYELECORR", False) # TEMPORARY UNTIL THEY ARE READY
    setConf("APPLYJETCORR", False)
    setConf("XSEC", 0.0664)
    #setConf("XSEC", 0.06761)
    setConf("LEPTON_SETUP", 2023)
    setConf("DATA_TAG", "pre_BPix")
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZZTo4l_5f_Summer23EraC_onlyCERN_good/crab_qqzz-4l-v3/241220_112540/0000/qqZZ_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZLZLTo4L_2023D":
    setConf("SAMPLENAME", "ZLZLTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLYELECORR", False) # TEMPORARY UNTIL THEY ARE READY
    setConf("APPLYJETCORR", False)
    setConf("XSEC", 0.003939)
    setConf("LEPTON_SETUP", 2023)
    setConf("NANOVERSION", 12)
    setConf("DATA_TAG", "post_BPix")
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZLZLTo4l_5f_Summer23EraD_goodSeeds/crab_qqzlzl-4l-v3/241220_103638/0000/qqZLZL_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZLZTTo4L_2023D":
    setConf("SAMPLENAME", "ZLZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLYELECORR", False) # TEMPORARY UNTIL THEY ARE READY
    setConf("APPLYJETCORR", False)
    setConf("XSEC", 0.01592) # Mean value of LHEWeight_originalXWGTUP
    setConf("LEPTON_SETUP", 2023)
    setConf("DATA_TAG", "post_BPix")
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZLZTTo4l_5f_Summer23EraD_goodSeeds/crab_qqzlzt-4l-v3/241220_103717/0000/qqZLZT_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZTZTTo4L_2023D":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLYELECORR", False) # TEMPORARY UNTIL THEY ARE READY
    setConf("APPLYJETCORR", False)
    #setConf("XSEC", 0.045113)
    setConf("XSEC", 0.04569)
    setConf("LEPTON_SETUP", 2023)
    setConf("DATA_TAG", "post_BPix")
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZTZTTo4l_5f_Summer23EraD_goodSeeds/crab_qqztzt-4l-v3/241220_103744/0000/qqZTZT_goodSeeds_hadd.root"
    ])

elif SampleToRun == "ZZTo4L_2023D_LO":
    setConf("SAMPLENAME", "ZTZTTo4L")
    setConf("APPLY_K_NNLOQCD_ZZQQB", False)
    setConf("APPLY_K_NNLOEW_ZZQQB", True)
    setConf("APPLYELECORR", False) # TEMPORARY UNTIL THEY ARE READY
    setConf("APPLYJETCORR", False)
    setConf("XSEC", 0.0664)
    #setConf("XSEC", 0.06761)
    setConf("LEPTON_SETUP", 2023)
    setConf("DATA_TAG", "post_BPix")
    setConf("NANOVERSION", 12)
    setConf("IsMC", True)
    setConf("store", "root://eos.grif.fr/")
    setConf("fileNames", [
        "/eos/grif/cms/llr/store/user/iehle/qqZZTo4l_5f_Summer23EraD_goodSeeds/crab_qqzz-4l-v3/241220_103811/0000/qqZZ_goodSeeds_hadd.root"
    ])

#####################################################################
### This import should be done AFTER all customizations (setConf calls)
from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *
######################################################################

### Tweak postprocessor parameters as necessary
p.prefetch=True # Prefetch remote files
p.longTermCache=True # keep prefetched files (useful for rerunning tests several times)
if len(p.inputFiles) == 1 :
    p.haddFileName = None # Skip final hadd
#p.maxEntries = 10000

### Select specific events to debug
#p.cut = "run==316239  && luminosityBlock==226 && event==284613817"

### Print out detailed candidate information for debug purposes
#from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
#p.cut = None # Remove preselction
#insertAfter(p.modules,"lepFiller",dumpEvents(level=-1),getConf("NANOVERSION", 11)) 

#p.branchsel=None #Read all branches
#p.outputbranchsel=None #Output all branches

#replace JSON
p.json = json

### Run the postprocessor
p.run()