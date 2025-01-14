from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import sys

class triggerAndSkim(Module):
    def __init__(self, isMC=True, PD="", era=2018, passThru=False):
        '''Apply good PV filter + Check trigger requirements and PD precedence rules for data.
        Parameters:
        PD: for data, should be "Muon", "EGamma", "MuonEG" to apply the PD precedence rules, or 
            "any" to select events passing any of the selected triggers regardless of PD precedence rules.
        passThru: do not filter events failing triggers (just store the module's variables). 
             Note that the good PV filter is still applied even if passThru = True.
        '''
        self.writeHistFile = False
        self.isMC = isMC
        self.PD = PD
        self.era = era
        self.passThru = passThru
        self.debug = False
        print("***triggerAndSkim: IsMC:", self.isMC, "PD:", self.PD, "era:", self.era, "passThru:", passThru, "- This module filters events.", flush=True)
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("HLT_passZZ4lEle", "O")   # pass Ele triggers
        self.out.branch("HLT_passZZ4lMu", "O")    # pass Muon triggers
        self.out.branch("HLT_passZZ4lMuEle", "O") # pass MuEle triggers
        self.out.branch("HLT_passZZ4l", "O")      # pass trigger requirements for the given PD (including PD precedence vetos) 


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        PD = self.PD

        ### Good PV filter
        if not event.Flag_goodVertices : return False # equivalent to: event.PV_npvsGood == 0


        ### Trigger requirements
        passTrigger = False
        if self.era == 2017 :
            passSingleEle = event.HLT_Ele35_WPTight_Gsf or event.HLT_Ele38_WPTight_Gsf or event.HLT_Ele40_WPTight_Gsf
            passSingleMu = event.HLT_IsoMu27
            passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_DoubleEle33_CaloIdL_MW
            passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 or event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8
            passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ or event.HLT_Mu8_DiEle12_CaloIdL_TrackIdL or event.HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ
            passTriEle = event.HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL
            passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        elif self.era == 2018 :
            passSingleEle = event.HLT_Ele32_WPTight_Gsf
            passSingleMu = event.HLT_IsoMu24
            passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_DoubleEle25_CaloIdL_MW
            passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ or event.HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ
            passTriEle = False
            passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        elif self.era == 2022 : # Checked that these are unprescaled in run 359751
            passSingleEle = event.HLT_Ele30_WPTight_Gsf #Note: we used Ele32 in 2018! 
            passSingleMu = event.HLT_IsoMu24
            passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_DoubleEle25_CaloIdL_MW
            passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ or event.HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ
            passTriEle = False
            passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        elif self.era == 2023 : # Checked that these are unprescaled, reference twikis for 2023 Eg & Muon Triggers https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIIISummary & https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2023
            passSingleEle = event.HLT_Ele30_WPTight_Gsf
            passSingleMu = event.HLT_IsoMu24
            passDiEle = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL
            passDiMu = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
            passMuEle = event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
            passTriEle = False
            passTriMu = event.HLT_TripleMu_10_5_5_DZ or event.HLT_TripleMu_12_10_5
        else:
            sys.exit("ERROR: era not supported: ", self.era)

        
        if self.isMC or PD == "any" :
            passTrigger = passDiEle or passDiMu or passMuEle or passTriEle or passTriMu or passSingleEle or passSingleMu

        else: # Data: ensure each event is taken only from a single PD
            if PD == "" : sys.exit("ERROR: PD must be set in data") # we may want to merge triggers for test runs 
            if ((PD=="DoubleEle" or PD=="DoubleEG"  or PD=="EGamma" ) and (passDiEle or passTriEle)) or \
               ((PD=="Muon" or PD=="DoubleMu"  or PD=="DoubleMuon") and (passDiMu or passTriMu) and not passDiEle and not passTriEle) or \
               ((PD=="MuEG" or PD=="MuonEG") and passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle) or \
               ((PD=="SingleElectron" or PD=="EGamma") and passSingleEle and not passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle) or \
               ((PD=="SingleMuon" or PD=="Muon") and passSingleMu and not passSingleEle and not passMuEle and not passDiMu and not passTriMu and not passDiEle and not passTriEle) :
                   passTrigger = True
            
        if self.debug :
            print(PD)
            print((passDiEle or passTriEle))
            print ((passDiMu or passTriMu) and not passDiEle and not passTriEle)
            print(passDiEle, passDiMu , passMuEle , passTriEle , passTriMu , passSingleEle , passSingleMu)
            print(passTrigger)

        self.out.fillBranch("HLT_passZZ4lEle", passSingleEle or passDiEle or passTriEle)
        self.out.fillBranch("HLT_passZZ4lMu", passSingleMu or passDiMu or passTriMu)
        self.out.fillBranch("HLT_passZZ4lMuEle", passMuEle)
        self.out.fillBranch("HLT_passZZ4l", passTrigger)

        if self.passThru :
            return True
        else :
            return passTrigger



