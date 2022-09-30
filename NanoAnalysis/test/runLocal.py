#!/usr/bin/env python
###
# Example for running the analysis locally, after customizing variables.
# Run with: 
# python runLocal.py
###
from __future__ import print_function
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf, insertAfter


### Customize processing variables

### K factors and weights (default: all off)
#setConf("APPLY_K_NNLOQCD_ZZGG", 1) # 0:None; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
#setConf("APPLY_K_NNLOQCD_ZZQQB", True)
#setConf("APPLY_K_NNLOEW_ZZQQB", True)
#setConf("APPLY_QCD_GGF_UNCERT", True)

### ggH125
setConf("SAMPLENAME", "ggH125")
setConf("XSEC", 48.58*0.0002745)
#setConf("runMELA", False)
#setConf("bestCandByMELA", False)
setConf("syncMode", True)

# setConf("store","/eos/cms")
# setConf("fileNames",[
#          "/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root",
#          "/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/445F4788-F322-C443-AB54-699C7716976A.root",
#          "/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/1DA1D554-E919-0146-9873-77E0B1B08DB5.root",
#          ])

# Custom-reprocessed nanoAOD file with updated FSR, corresponding to:/store/mc/RunIIAutumn18NanoAODv7/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/BA6D7F40-ED5E-7D4E-AB14-CE8A9C5DE7EC.root
setConf("store","")
setConf("fileNames",["/eos/user/n/namapane/H4lnano/ggH125_fixedFSR.root"])


### ZH125
# setConf("SAMPLENAME", "ZH125")
# setConf("store","root://cms-xrd-global.cern.ch/")
# setConf("fileNames",["/store/mc/RunIIAutumn18NanoAODv7/ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/C670B342-C02E-1848-AE6A-0B5550E3DFE3.root"])

### qqZZ: /ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/NANOAODSIM
# setConf("SAMPLENAME", "qqZZ")
# setConf("XSEC", 1.256)
# setConf("fileNames",["/store/mc/RunIIAutumn18NanoAODv7/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/3BD1674D-804E-0E45-B106-EAB3C0D59E32.root"])

### High mass sample
# setConf("SAMPLENAME", "VBFH3000")
# setConf("fileNames",["/store/mc/RunIIAutumn18NanoAODv7/VBF_HToZZTo4L_M3000_13TeV_powheg2_JHUGenV7011_pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A82C31DF-0A6F-B44F-857B-89BFD72CEFEA.root"])


#Options to run on DATA
#setConf("IsMC", False)
#setConf("PD", "any")
#setConf("SAMPLENAME", "test")
#setConf("store","root://cms-xrd-global.cern.ch/")
#setConf("fileNames",["/eos/cms/store/data/Run2018A/MuonEG/NANOAOD/02Apr2020-v1/50000/420527E8-C6CA-9745-AD5D-6ADF63B808B5.root"])
#setConf("fileNames",["/store/data/Run2018A/DoubleMuon/NANOAOD/02Apr2020-v1/30000/CD25F914-CA18-2847-989C-91E6BF71E22F.root"])


# Select specific events to debug
#setConf("preselection","run==316239  && luminosityBlock==226 && event==284613817")


### This import should be done AFTER all customizations (setConf calls)
from ZZAnalysis.NanoAnalysis.nanoZZ4lAnalysis import *

### Tweak postprocessor parameters as necessary
p.prefetch=False # Skip local prefetching
p.haddFileName=None # Skip final hadd
#p.maxEntries = 1000

### Print out detailed candidate information for debug purposes
#p.cut = None # Remove preselction
#from ZZAnalysis.NanoAnalysis.dumpEvents import dumpEvents
#insertAfter(p.modules,"lepFiller",dumpEvents(level=-1)) 

#p.branchsel=None #Read all branches
#p.outputbranchsel=None #Output all branches


### Run the postprocessor
p.run()
