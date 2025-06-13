from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR

from functools import cmp_to_key
# from ROOT import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar, TLorentzVector
from ROOT import TLorentzVector

from ROOT.Math import LorentzVector, PxPyPzE4D, PtEtaPhiM4D, Boost, LorentzRotation
import ROOT
import numpy as np

from ctypes import c_float
import Mela

class AngularVars:
    def __init__(self, ZZ):
        self.ZZ = ZZ
        self.Z  = {"1": ZZ.Z1, "2": ZZ.Z2}

        # Cache useful p4s
        self.ZZ_p4 = self.ZZ.p4
        self.Z1_p4 = self.ZZ.Z1.p4
        self.Z2_p4 = self.ZZ.Z2.p4
        
        # Define boosts
        self.boost_to_4l = self._getBoost(self._getLorentzVec(self.ZZ_p4))
        self.boost_to_z1 = self._getBoost(self._getLorentzVec(self.Z1_p4))
        self.boost_to_z2 = self._getBoost(self._getLorentzVec(self.Z2_p4))


    def _getLorentzVec(self, tLorVec):
        return LorentzVector(PxPyPzE4D('double'))(
            tLorVec.Px(),
            tLorVec.Py(),
            tLorVec.Pz(),
            tLorVec.E()
        )

    def _getPosLep(self, z_cand):
        lep_p4 = z_cand.l1DressedP4 if z_cand.l1.charge == 1 else z_cand.l2DressedP4
        return self._getLorentzVec(lep_p4)

    def _getNegLep(self, z_cand):
        lep_p4 = z_cand.l1DressedP4 if z_cand.l1.charge == -1 else z_cand.l2DressedP4
        return self._getLorentzVec(lep_p4)

    def _getBoost(self, p4):
        return Boost(p4.BoostToCM())

    def cosTheta(self, z):
        lab_z_cand = self.Z[z]
        
        # Four vectors in the lab frame for 4l system, Z, and its associated l+
        lab_pos_lep = self._getPosLep(lab_z_cand)

        zcand_4l_com  = self.boost_to_4l*self.Z1_p4 if z=="1" else self.boost_to_4l*self.Z2_p4
        pos_lep_z_com = self.boost_to_z1*lab_pos_lep if z=="1" else self.boost_to_z2*lab_pos_lep

        return ROOT.Math.VectorUtil.CosTheta(pos_lep_z_com, zcand_4l_com)

    def cosThetaStar(self, z="1"):
        # THIS DEFINITION IS STRONGLY DEPENDENT ON DEFINTION OF Z1 : SORTING BY pT GIVES VERY ASSYMETRIC DISTRIBUTION
        
        zcand_4l_com = self.boost_to_4l*self.Z1_p4 if z=="1" else self.boost_to_4l*self.Z2_p4
        return ROOT.Math.VectorUtil.CosTheta(zcand_4l_com, self.ZZ_p4)

    def delPhiStar(self):
        lab_z1_cand, lab_z2_cand = self.Z["1"], self.Z["2"]

        # l+ p4's in lab frame
        lab_pos_lep_z1 = self._getPosLep(lab_z1_cand)
        lab_pos_lep_z2 = self._getPosLep(lab_z2_cand)

        pos_lep_z1_cm = self.boost_to_z1*lab_pos_lep_z1
        pos_lep_z2_cm = self.boost_to_z2*lab_pos_lep_z2

        z1_4l_p4 = self.boost_to_4l*self.Z1_p4
        z2_4l_p4 = self.boost_to_4l*self.Z2_p4

        # Calculate phi*'s (angle between l+ in its Z's rest frame and the Z's direction of flight in the event CM (4l frame))
        phiStar_z1 = ROOT.Math.VectorUtil.Angle(z1_4l_p4, pos_lep_z1_cm)
        phiStar_z2 = ROOT.Math.VectorUtil.Angle(z2_4l_p4, pos_lep_z2_cm)

        delPhiStar = np.abs(phiStar_z1 - phiStar_z2)

        return min(delPhiStar, 2*np.pi - delPhiStar)

    def delPhi(self):
        lab_z1_cand, lab_z2_cand = self.Z["1"], self.Z["2"]

        # l+(-) p4's in lab frame
        pos_lep_z1 = self._getPosLep(lab_z1_cand)
        neg_lep_z2 = self._getNegLep(lab_z2_cand)

        return np.abs(ROOT.Math.VectorUtil.DeltaPhi(pos_lep_z1, neg_lep_z2))

    def delRapidity(self):
        z1_p4, z2_p4 = self._getLorentzVec(self.Z["1"].p4), self._getLorentzVec(self.Z["2"].p4)

        return np.abs(z1_p4.Rapidity() - z2_p4.Rapidity())

class candProps:
    def __init__(self, final_cands, region_filters):
        self.final_cands    = final_cands

        self.region_bools = {reg: [] for reg in region_filters}
        
        self.prop_names = ("mass", "pt", "eta", "phi", "massPreFSR", "Z1mass", "Z1pt", "Z1eta", "Z1phi", "Z1flav",
                            "Z1pt", "Z1eta", "Z1phi", "Z2mass", "Z2flav", "Z1l1Idx", "Z1l2Idx", "Z2l1Idx", "Z2l2Idx")

        self.props = dict(
            mass         = lambda cand: cand.M,
            pt           = lambda cand: cand.p4.Pt(),
            eta          = lambda cand: cand.p4.Eta(),
            phi          = lambda cand: cand.p4.Phi(),
            cosTheta1    = lambda cand: cand.cosTheta1,
            cosTheta3    = lambda cand: cand.cosTheta3,
            cosThetaStar = lambda cand: cand.cosThetaStar,
            delPhiStar   = lambda cand: cand.delPhiStar,
            delPhi       = lambda cand: cand.delPhi,
            delRapidity  = lambda cand: cand.delRapidity,
            massPreFSR   = lambda cand: cand.massPreFSR(),
            Z1mass       = lambda cand: cand.Z1.M,
            Z1pt         = lambda cand: cand.Z1.p4.Pt(),
            Z1eta        = lambda cand: cand.Z1.p4.Eta(),
            Z1phi        = lambda cand: cand.Z1.p4.Phi(),
            Z1flav       = lambda cand: cand.Z1.finalState(),
            Z2mass       = lambda cand: cand.Z2.M,
            Z2pt         = lambda cand: cand.Z2.p4.Pt(),
            Z2eta        = lambda cand: cand.Z2.p4.Eta(),
            Z2phi        = lambda cand: cand.Z2.p4.Phi(),
            Z2flav       = lambda cand: cand.Z2.finalState(),
            Z1l1Idx      = lambda cand: cand.Z1.l1Idx,
            Z1l2Idx      = lambda cand: cand.Z1.l2Idx,
            Z2l1Idx      = lambda cand: cand.Z2.l1Idx,
            Z2l2Idx      = lambda cand: cand.Z2.l2Idx
        )

        self.branches = {prop: [] for prop in self.props.keys()}

        self._fill_props()

        self.branches.update(self.region_bools)
    
    def _fill_props(self):
        for passing_region, cand in self.final_cands.items():
            self.region_bools[passing_region].append(True)
            for reg in self.region_bools:
                if reg==passing_region: continue
                else: self.region_bools[reg].append(False)

            for prop, prop_list in self.branches.items():
                prop_list.append(self.props[prop](cand))

class ZZFiller(Module):

    def __init__(self, bestCandByMELA, MELA, isMC, year, data_tag, processCR=False, addZL=False, filter='Cands', candsToStore='BestCandOnly',  debug=False):
        """Build candidates:
        -ZZCand: SR candidates. The index of the best ZZ candidate in each event is stored as bestCandIdx
        -ZLLCand: CR candidates (SS, 3P1F, 2P2F, SIP CRs, with indices: ZZLLbestSSIdx, ZLLbest3P1FIdx, ZLLbest2P2FIdx, ZLLbestSIPCRIdx)
        -ZLCand: Z+L CR, for fake rate computation (only the index of the additional lepton is stored)
        -ZCand: Zs referenced in the above collections (note that when filter='Cands', are applied, only events with at least one candidate in the above collections are retained, so the tree will not contain all events with one Z)
        Parameters:
          bestCandByMELA: True = select best candidate by KD; False = with best Z1/highest-pTsum Z2 
          year: data taking year        
          MELA: The MELA object passed from nanoZZ4lAnalysis.py 
          processCR: add ZLLCand CR collections
          data_tag: subperiod or processing, e.g. "pre_EE" (currently unused)
          addZL: add ZL CR
          filter: criteria to keep or skip events:
                  'Cands' = keep only events with at least one ZZ, ZLL, or ZL candidate are kept
                  'Z' = keep any event that has at least one good Z candidate (passing the analysis Z selection criteria, and 12<mll<120)
                  '3L_20_10' = keep all events with 3 good leptons, pt1>20, pt2>10; useful for trigger studies
                  'NoFilter' = don't filter events
          candsToStore: which candidates should be stored in the ZZCand collection:
                  'BestCandOnly' = only the best SR candidate in the event is saved (default)
                  'AllCands' = keep all SR candidates passing the full selection and analysis cuts (including permutations of leptons).
                  'AllWithRelaxedMuId' = keep any SR candidate that can be made, even if leptons don't pass ID cuts (useful for ID cut optimization studies).
                   Note that this option does not affect the ZLLCand collection: for each CR that is activated, only the best candidate is stored.
        """
        print("***ZZFiller: isMC:", isMC, "year:", year, "data_tag:", data_tag, "bestCandByMELA:", bestCandByMELA, "filter:", filter, "candsToStore:", candsToStore, ("- This module filters events." if filter!='NoFilter' else ""),  flush=True)
        self.writeHistFile = False
        self.mela = MELA
        self.isMC = isMC
        self.year = year
        if bestCandByMELA :
            self.bestCandCmp = self.bestCandByDbkgKin
        else:
            self.bestCandCmp = self.bestCandByZ1Z2

        self.addSSCR = processCR
        self.addOSCR = processCR
        self.addSIPCR = processCR
        self.addZLCR = addZL

        self.DATA_TAG = data_tag

        # Translate option strings into enums, for efficiency
        self.noFilter, self.filterOnCands, self.filterOnZ, self.filter_3L_20_10 = range(0,4)
        filters = {'NoFilter':self.noFilter,
                   'Cands':self.filterOnCands,
                   'Z':self.filterOnZ,
                   '3L_20_10':self.filter_3L_20_10}
        try:
            self.filterType = filters[filter]
        except :
            raise ValueError("ZZFiller: filter =", filter, "not supported")

        self.BestCandOnly, self.AllCands, self.AllWithRelaxedMuId =  range(0,3)
        storeOptions = {'BestCandOnly':self.BestCandOnly,
                        'AllCands':self.AllCands,
                        'AllWithRelaxedMuId':self.AllWithRelaxedMuId}
        try:
            self.candsToStore = storeOptions[candsToStore]
        except :
            raise ValueError("ZZFiller: candsToStore =", candsToStore, "not supported")

        self.DEBUG = debug
        self.ZmassValue = 91.1876;
        self.shell_cond = 10 # on-shell condition for polarization analysis (requires |mll - zmass| < shell_cond)

        # Pre-selection of leptons used to reduce combinatorial when building Z and LL candidates.
        # Note that the actual lepton selection cuts for SR and CR are applied later; this preselection only affects what
        # leptons are considered in making the combinatorial (ie processing speed)
        
        if self.candsToStore != self.AllWithRelaxedMuId :
            # Normal case: is the full ID + iso if only the SR is considered, or the relaxed ID if CRs are also filled.
            if self.addSIPCR or self.addOSCR or self.addSSCR or self.addZLCR:
                self.leptonPresel = (lambda l : l.ZZRelaxedIdNoSIP) # minimal selection good for all CRs: no SIP, no ID, no iso
            else : # SR only
                self.leptonPresel = (lambda l : l.ZZFullSel)
        else :
            # Use relaxed muon preselection for muon ID studies: fully relax muon ID.
            # Electron ID is unchanged; FullSel electrons are preselected in this case, so the effect of cut variations can be studied for muons only.
            # for this reason, CRs cannot be properly built.
            if self.addSIPCR or self.addOSCR or self.addSSCR or self.addZLCR:
                raise Exception("WARNING: CRs are not supported when candsToStore==AllWithRelaxedMuId")
            self.leptonPresel = (lambda l : (abs(l.pdgId)==13 and l.pt>5 and abs(l.eta) < 2.4) or (abs(l.pdgId)==11 and l.ZZFullSel))

        
        

        # Example of adding control histograms (requires self.writeHistFile = True)
        # def beginJob(self,histFile=None, histDirName=None):
        #    Module.beginJob(self, histFile, histDirName+"_ZZFiller")
        #    self.histFile=None # Hack to prevent histFile to be closed before other modules write their histograms
        #    self.h_ZZMass = ROOT.TH1F('ZZMass','ZZMass',130,70,200)
        #    self.addObject(self.h_ZZMass)


    def endJob(self):
         print("", flush=True)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("nZCand", "I", title="Z candidates passing the full H4l selection")
        self.out.branch("ZCand_mass", "F", lenVar="nZCand", title="mass")
        self.out.branch("ZCand_pt", "F", lenVar="nZCand")
        self.out.branch("ZCand_eta", "F", lenVar="nZCand", limitedPrecision=16)
        self.out.branch("ZCand_rapidity", "F", lenVar="nZCand", limitedPrecision=12)
        self.out.branch("ZCand_phi", "F", lenVar="nZCand", limitedPrecision=16)
        self.out.branch("ZCand_flav", "I", lenVar="nZCand", title="Product of the pdgIds of the 2 daughters")
        self.out.branch("ZCand_l1Idx", "S", lenVar="nZCand", title="index of 1st daughter in Electron+Muon merged collection")
        self.out.branch("ZCand_l2Idx", "S", lenVar="nZCand", title="index of 2nd daughter in Electron+Muon merged collection")
        self.out.branch("ZCand_fsr1Idx", "S", lenVar="nZCand", title="index of FSR associated to l1 (-1 if none)")
        self.out.branch("ZCand_fsr2Idx", "S", lenVar="nZCand", title="index of FSR associated to l2 (-1 if none)") 
        self.out.branch("bestZIdx", "S", title="Best Z in the event (mass closest to mZ)")

        self.out.branch("nZZCand", "I", title="ZZ candidates passing the full H4l selection")
        self.out.branch("ZZCand_mass", "F", lenVar="nZZCand", title="mass")
        self.out.branch("ZZCand_pt", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_eta", "F", lenVar="nZZCand", limitedPrecision=16)
        self.out.branch("ZZCand_rapidity", "F", lenVar="nZZCand", limitedPrecision=12)
        self.out.branch("ZZCand_phi", "F", lenVar="nZZCand", limitedPrecision=16)
        self.out.branch("ZZCand_massPreFSR", "F", lenVar="nZZCand", title="mass without FSR photons")
        self.out.branch("ZZCand_Z1mass", "F", lenVar="nZZCand", title="Z1 mass")
        self.out.branch("ZZCand_Z1flav", "I", lenVar="nZZCand", title="Product of the pdgIds of the 2 Z1 daughters")
        self.out.branch("ZZCand_Z1pt", "F", lenVar="nZZCand", title="Z1 pt")
        self.out.branch("ZZCand_Z1eta", "F", lenVar="nZZCand", title="Z1 eta")
        self.out.branch("ZZCand_Z1phi", "F", lenVar="nZZCand", title="Z1 phi")
        self.out.branch("ZZCand_Z1rapidity", "F", lenVar="nZZCand", title="Z1 rapidity")
        self.out.branch("ZZCand_Z2mass", "F", lenVar="nZZCand", title="Z2 mass")
        self.out.branch("ZZCand_Z2flav", "I", lenVar="nZZCand", title="Product of the pdgIds of the 2 Z2 daughters")
        self.out.branch("ZZCand_Z2pt", "F", lenVar="nZZCand", title="Z2 pt")
        self.out.branch("ZZCand_Z2eta", "F", lenVar="nZZCand", title="Z2 eta")
        self.out.branch("ZZCand_Z2phi", "F", lenVar="nZZCand", title="Z2 phi")
        self.out.branch("ZZCand_Z2rapidity", "F", lenVar="nZZCand", title="Z2 rapidity")
        self.out.branch("ZZCand_KD", "F", lenVar="nZZCand", title="Kinematic discriminant for the choice of best candidate")
        self.out.branch("ZZCand_Z2sumpt", "F", lenVar="nZZCand", title="sum of Z2 daughter pts (used in the choice of best candidate)")
        # Note: lepton indices are numbered for leps=list(electrons)+list(muons) and run up to nlep=len(leps);
        # no special ordering of l1, l2 is applied
        self.out.branch("ZZCand_Z1l1Idx", "S", lenVar="nZZCand", title="Index of 1st Z1 daughter in the Electron+Muon merged collection")
        self.out.branch("ZZCand_Z1l2Idx", "S", lenVar="nZZCand", title="Index of 2nd Z1 daughter in the Electron+Muon merged collection")
        self.out.branch("ZZCand_Z2l1Idx", "S", lenVar="nZZCand", title="Index of 1st Z2 daughter in the Electron+Muon merged collection")
        self.out.branch("ZZCand_Z2l2Idx", "S", lenVar="nZZCand", title="Index of 2nd Z2 daughter in the Electron+Muon merged collection")

        self.out.branch("ZZCand_cosTheta1", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_cosTheta3", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_cosThetaStar", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_delPhiStar","F", lenVar="nZZCand")
        self.out.branch("ZZCand_delPhi", "F", lenVar="nZZCand")
        self.out.branch("ZZCand_delRapidity", "F", lenVar="nZZCand")

        self.out.branch("bestCandIdx", "S", title="Index of seleced ZZCand candidate in the event")

        if self.addSSCR or self. addOSCR or self.addSIPCR :
            self.out.branch("nZLLCand", "I", title="Z+LL control region candidates")
            self.out.branch("ZLLCand_mass", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_massPreFSR", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_pt", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_eta", "F", lenVar="nZLLCand", limitedPrecision=16)
            self.out.branch("ZLLCand_rapidity", "F", lenVar="nZLLCand", limitedPrecision=12)
            self.out.branch("ZLLCand_phi", "F", lenVar="nZLLCand", limitedPrecision=16)
            self.out.branch("ZLLCand_Z1mass", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z1flav", "I", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z1pt", "F", lenVar="nZLLCand", title="Z1 pt")
            self.out.branch("ZLLCand_Z1eta", "F", lenVar="nZLLCand", title="Z1 eta")
            self.out.branch("ZLLCand_Z1phi", "F", lenVar="nZLLCand", title="Z1 phi")
            self.out.branch("ZLLCand_Z1rapidity", "F", lenVar="nZLLCand", title="Z1 rapidity")
            self.out.branch("ZLLCand_Z2mass", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z2flav", "S", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z2pt", "F", lenVar="nZLLCand", title="Z2 pt")
            self.out.branch("ZLLCand_Z2eta", "F", lenVar="nZLLCand", title="Z2 eta")
            self.out.branch("ZLLCand_Z2phi", "F", lenVar="nZLLCand", title="Z2 phi")
            self.out.branch("ZLLCand_Z2rapidity", "F", lenVar="nZLLCand", title="Z2 rapidity")
            self.out.branch("ZLLCand_Z1l1Idx", "S", lenVar="nZLLCand") 
            self.out.branch("ZLLCand_Z1l2Idx", "S", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z2l1Idx", "S", lenVar="nZLLCand")
            self.out.branch("ZLLCand_Z2l2Idx", "S", lenVar="nZLLCand")
            self.out.branch("ZLLCand_KD", "F", lenVar="nZLLCand", limitedPrecision=12)
            self.out.branch("ZLLbestSSIdx", "S", title="best candidate for the SS CR")
            self.out.branch("ZLLbest2P2FIdx", "S", title="best candidate for the 2P2F CR")
            self.out.branch("ZLLbest3P1FIdx", "S", title="best candidate for the 3P1F CR")
            self.out.branch("ZLLbestSIPCRIdx", "S", title="best candidate for the SIP CR")
            self.out.branch("ZLLbestHighMassOSSIPIdx", "S", title="best candidate for the High Mass OS SIP CR")
            self.out.branch("ZLLbestMidMassOSSIPIdx", "S", title="best candidate for the Mid Mass OS SIP CR")
            self.out.branch("ZLLbestLowMassOSSIPIdx", "S", title="best candidate for the Low Mass OS SIP CR")
            self.out.branch("ZLLbestHighMassSSSIPIdx", "S", title="best candidate for the High Mass OS SIP CR")
            self.out.branch("ZLLbestMidMassSSSIPIdx", "S", title="best candidate for the Mid Mass OS SIP CR")
            self.out.branch("ZLLbestLowMassSSSIPIdx", "S", title="best candidate for the Low Mass OS SIP CR")
            self.out.branch("ZLLbestHighMassSSRelaxedIdx", "S", title="best candidate for the High Mass SS Relaxed Sel CR (For Z+X shape)")

            self.out.branch("ZLLCand_cosTheta1", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_cosTheta3", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_cosThetaStar", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_delPhiStar","F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_delPhi", "F", lenVar="nZLLCand")
            self.out.branch("ZLLCand_delRapidity", "F", lenVar="nZLLCand")

        if self.addZLCR :            
            self.out.branch("ZLCand_lepIdx", "S", title="Index of extra lep for the ZL CR")


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if self.DEBUG : print ('Event {}:{}:{}'.format(event.run,event.luminosityBlock,event.event))
 
        # Collections
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        leps = list(electrons)+list(muons)
        nlep=len(leps)

        
        ### Apply initial event filter.
        if self.filterType == self.filterOnCands : # Filter on ZZ, ZLL, ZL candidates will happen later. Here We just apply a minimal pre-filter (>=3L if ZL is included, 4 otherwise)
            if nlep < 3 or (not self.addZLCR and nlep < 4) :
                return False
        elif self.filterType == self.filterOnZ : # Filter on Z candidates will happen later. Here We just apply a minimal pre-filter (>=2L)
            if event.nMuon < 2 and event.nElectron < 2 :
                return False
        elif self.filterType == self.filter_3L_20_10 : # Select events with 3 good leptons, pt>20/10 - no further filter will be applied later
            nGoodLeps = 0
            nLeps10 = 0
            nLeps20 = 0
            for lep in leps :
                if lep.ZZFullSel :
                    nGoodLeps += 1
                    if lep.pt>10. : nLeps10 +=1
                    if lep.pt>20. : nLeps20 +=1
            if nGoodLeps < 3 or nLeps20 < 1 or nLeps10 < 2 :
                return False
        
        Zs = [] # all Z cands, used to build SR and CRs
        selZs = [] # Selected Zs, to be written out
        bestZIdx = -1 # index of the best Z in the event, among selZs
        ZZs = []
        bestCandIdx = -1
        ZLLs = [] 
        ZLLsTemp = []
        best3P1FCRIdx = -1
        best2P2FCRIdx = -1
        bestSSCRIdx = -1
        bestSIPCRIdx = -1
        bestHighMassOSSIPIdx = -1
        bestMidMassOSSIPIdx = -1
        bestLowMassOSSIPIdx = -1
        bestHighMassSSSIPIdx = -1
        bestMidMassSSSIPIdx = -1
        bestLowMassSSSIPIdx = -1
        bestHighMassSSRelaxedIdx = -1
        ZLCand_lIdx = -1 #index of the additional lepton for the Z+L CR

        ### Z combinatorial over selected leps (after FSR-corrected ISO cut for muons)
        for i,l1 in enumerate(leps):
            if self.leptonPresel(l1) :
                for j in range(i+1,nlep):
                    l2 = leps[j]
                    nPassLep = 0
                    isOSSF = False
                    isSR   = False
                    is1FCR = False
                    is2FCR = False
                    isSSCR = False
                    isSIPCR = False
                    isOSSIPCR = False
                    isSSSIPCR = False

                    isSSRelaxedCR = False
                    if self.leptonPresel(l2):
                        if l1.pdgId == -l2.pdgId : #OS,SF: candidate for signal region 
                            isOSSF = True
                                
                            if l1.ZZRelaxedId and l2.ZZRelaxedId:
                                nPassLep = int(l1.ZZFullSel) + int(l2.ZZFullSel) # 0,1,2 leptons passing full sel
                                if nPassLep == 2 : isSR = True # SR (to be used to set the best candidate)
                                if self.addOSCR :
                                    if nPassLep == 0 : is2FCR = True
                                    elif nPassLep == 1 : is1FCR = True

                            # For OS/SS transfer function using SIP method
                            if l1.ZZFullSelNoSIP and l2.ZZFullSelNoSIP:
                                isOSSIPCR = True

                        elif l1.pdgId == l2.pdgId : # SS control regions
                            if self.addSSCR and l1.ZZRelaxedId and l2.ZZRelaxedId : isSSCR = True
                            if self.addSIPCR and l1.ZZFullSelNoSIP and l2.ZZFullSelNoSIP : isSIPCR = True
                            if l1.ZZFullSelNoSIP and l2.ZZFullSelNoSIP : isSSSIPCR = True
                            if (l1.ZZRelaxedId and l1.passIso) and (l2.ZZRelaxedId and l2.passIso): isSSRelaxedCR = True
                            if not (isSSCR or isSIPCR) : continue
                        else:
                            continue
                        
                        # Set a default order for OS leptons in the Z candidate: 1=l+, 2=l-
                        idx1, idx2 = i, j
                        if (l1.pdgId*l2.pdgId<0 and l1.pdgId>0) :
                            idx1, idx2 = idx2, idx1
                        
                        aZ = self.ZCand(idx1, idx2, leps, fsrPhotons)
                        aZ.on_shell = abs(aZ.M - self.ZmassValue) < self.shell_cond
                        aZ.isOSSF = isOSSF
                        aZ.isSR   = isSR and aZ.on_shell # Polarization requirement
                        aZ.is1FCR = is1FCR
                        aZ.is2FCR = is2FCR
                        aZ.isSSCR = isSSCR
                        aZ.isSIPCR = isSIPCR
                        aZ.isOSSIPCR = isOSSIPCR
                        aZ.isSSSIPCR = isSSSIPCR
                        aZ.isSSRelaxedCR = isSSRelaxedCR
                        
                        zmass = aZ.M
                        if self.DEBUG: print('Z={:.4g} pt1={:.3g} pt2={:.3g} fsr1={} fsr2={} SR={} 1F={} 2F={} SS={} OSSIP={}'.format(zmass, l1.pt, l2.pt, l1.fsrPhotonIdx,  l2.fsrPhotonIdx, isSR, is1FCR, is2FCR, isSSCR, isSIPCR))
                        if (zmass>12. and zmass<120.):
                            Zs.append(aZ)
                            if self.candsToStore == self.AllWithRelaxedMuId : # For ID studies, store any OS Z made with preselected leptons
                                if aZ.isOSSF :
                                    selZs.append(aZ)
                            else : # Default: store only Zs passing the SR lepton selection
                                if aZ.isSR :
                                    selZs.append(aZ)

                            # Choose best Z among SR Zs and store index wihin the selZs list
                            if aZ.isSR :
                                if (bestZIdx<0 or abs(zmass-self.ZmassValue)<abs(selZs[bestZIdx].M-self.ZmassValue)) :
                                    bestZIdx = len(selZs)-1
                                

        ### Build ZZ and ZLL combinations passing the ZZ selection
        if len(Zs) >= 2:
            ### Signal region
            for iZ,aZ in enumerate(Zs):
                for jZ in range(iZ+1, len(Zs)):                    
                    if Zs[iZ].isOSSF and Zs[jZ].isOSSF : # 2 OSSSF Zs
                        ZZ = self.makeCand(Zs[iZ],Zs[jZ])
                        if ZZ == None: continue
                        ZZs.append(ZZ)

                        if self.DEBUG : print("ZZ:", len(ZZs)-1, ZZ.p4.M(), ZZ.Z1.M, ZZ.Z2.M, ZZ.Z2.sumpt(), ZZ.finalState(), ZZ.p_GG_SIG_ghg2_1_ghz1_1_JHUGen, ZZ.p_QQB_BKG_MCFM, ZZ.KD, (ZZ.Z1.isSR and ZZ.Z2.isSR))

                        #Search for the the best cand in the SR (ie among those passing full ID cuts)
                        if ZZ.Z1.isSR and ZZ.Z2.isSR and ZZ.HighMass:
                            if bestCandIdx<0 or self.bestCandCmp(ZZ,ZZs[bestCandIdx]) < 0: bestCandIdx = len(ZZs)-1

            if self.DEBUG : print("bestCand:", bestCandIdx)


            ### ZLL combinations for control regions, made of 1 good Z + 1 ll pair;
            ### these are considered only if no SR candidate is present in the event
            if bestCandIdx < 0 and (self.addSSCR or self. addOSCR or self.addSIPCR) :
                for iZ1,Z1 in enumerate(Zs):
                    if Z1.isSR : 
                        for iZ2,Z2 in enumerate(Zs):
                            if Z2.is1FCR or Z2.is2FCR or Z2.isSSCR or Z2.isOSSIPCR or Z2.isSSSIPCR:
                                ZLL = self.makeCand(Z1, Z2, sortZsByMass=False, fillIDVars=False)
                                if ZLL == None: continue
                                if ZLL.Z2.is2FCR  and (best2P2FCRIdx<0 or self.bestCandCmp(ZLL,ZLLsTemp[best2P2FCRIdx]) < 0) : best2P2FCRIdx = len(ZLLsTemp)
                                if ZLL.Z2.is1FCR  and (best3P1FCRIdx<0 or self.bestCandCmp(ZLL,ZLLsTemp[best3P1FCRIdx]) < 0) : best3P1FCRIdx = len(ZLLsTemp)
                                if ZLL.Z2.isSSCR  and (bestSSCRIdx<0 or self.bestCandCmp(ZLL,ZLLsTemp[bestSSCRIdx]) < 0) : bestSSCRIdx = len(ZLLsTemp)
                                if ZLL.Z2.isOSSIPCR:
                                    # if ZLL.HighMass and ZLL.Z2.on_shell and (bestHighMassOSSIPIdx < 0 or self.bestCandCmp(ZLL, ZLLsTemp[bestHighMassOSSIPIdx]) < 0) : bestHighMassOSSIPIdx = len(ZLLsTemp)
                                    if ZLL.MidMass and (bestMidMassOSSIPIdx < 0 or self.bestCandCmp(ZLL, ZLLsTemp[bestMidMassOSSIPIdx]) < 0) : bestMidMassOSSIPIdx = len(ZLLsTemp)
                                    if ZLL.LowMass and (bestLowMassOSSIPIdx < 0 or self.bestCandCmp(ZLL, ZLLsTemp[bestLowMassOSSIPIdx]) < 0) : bestLowMassOSSIPIdx = len(ZLLsTemp)
                                if ZLL.Z2.isSSSIPCR:
                                    if ZLL.HighMass and ZLL.Z2.on_shell and (bestHighMassSSSIPIdx < 0 or self.bestCandCmp(ZLL, ZLLsTemp[bestHighMassSSSIPIdx]) < 0) : bestHighMassSSSIPIdx = len(ZLLsTemp)
                                    if ZLL.MidMass and (bestMidMassSSSIPIdx < 0 or self.bestCandCmp(ZLL, ZLLsTemp[bestMidMassSSSIPIdx]) < 0) : bestMidMassSSSIPIdx = len(ZLLsTemp)
                                    if ZLL.LowMass and (bestLowMassSSSIPIdx < 0 or self.bestCandCmp(ZLL, ZLLsTemp[bestLowMassSSSIPIdx]) < 0) : bestLowMassSSSIPIdx = len(ZLLsTemp)
                                if ZLL.Z2.isSSRelaxedCR:
                                    if ZLL.HighMass and ZLL.Z2.on_shell and (bestHighMassSSRelaxedIdx < 0 or self.bestCandCmp(ZLL, ZLLsTemp[bestHighMassSSRelaxedIdx]) < 0): bestHighMassSSRelaxedIdx = len(ZLLsTemp)
                                        
                                if ZLL.Z2.isSIPCR and (bestSIPCRIdx<0 or self.bestCandCmp(ZLL,ZLLsTemp[bestSIPCRIdx]) < 0) : bestSIPCRIdx = len(ZLLsTemp)
                                ZLLsTemp.append(ZLL)
                                
                # Check that only one candidate is selected in each event for the 2 CRs of the OS method?
                # Actually not needed, at the overlap is accounted for in the method
#                if best2P2FCRIdx >= 0 and best3P1FCRIdx >= 0 :
#                    print ('WARNING: event {}:{}:{} has CR candidates in both 2P2F and 3P1F regions'.format(event.run,event.luminosityBlock,event.event))

                # Different SIP CRs can contain overlapping candidates, choose the best among them
                sip_crs = [bestMidMassOSSIPIdx, bestLowMassOSSIPIdx, bestHighMassSSSIPIdx, bestMidMassSSSIPIdx, bestLowMassSSSIPIdx]
                n_crs = len(sip_crs)
                if sum([idx >= 0 for idx in sip_crs]) > 1:
                    for i, za_idx in enumerate(sip_crs):
                        for j in range(i+1, n_crs):
                            zb_idx = sip_crs[j]
                            if zb_idx == -1: continue
                            za, zb = ZLLsTemp[za_idx], ZLLsTemp[zb_idx]

                            # Check if leps overlap
                            if bool(set(za.leps()) & set(zb.leps())):
                                if self.bestCandCmp(za, zb) < 0:
                                    if   j == 1: bestLowMassOSSIPIdx  = -1
                                    elif j == 2: bestHighMassSSSIPIdx = -1
                                    elif j == 3: bestMidMassSSSIPIdx  = -1
                                    elif j == 4: bestLowMassSSIPIdx   = -1
                                else:
                                    if   i == 0: bestMidMassOSSIPIdx  = -1
                                    elif i == 1: bestLowMassOSSIPIdx  = -1
                                    elif i == 2: bestHighMassSSSIPIdx = -1
                                    elif i == 3: bestMidMassSSSIPIdx  = -1

                # Store only ZLL candidates that belong to at least 1 CR
                for iZLL, ZLL in enumerate(ZLLsTemp) :
                    select = False
                    if iZLL == best2P2FCRIdx :
                        best2P2FCRIdx = len(ZLLs)
                        select = True
                    if iZLL == best3P1FCRIdx :
                        best3P1FCRIdx = len(ZLLs)
                        select = True
                    if iZLL == bestSSCRIdx :
                        bestSSCRIdx = len(ZLLs)
                        select = True
                    if iZLL == bestSIPCRIdx :
                        bestSIPCRIdx = len(ZLLs)
                        select = True
                    if iZLL == bestHighMassOSSIPIdx:
                        bestHighMassOSSIPIdx = len(ZLLs)
                        select = True
                    if iZLL == bestMidMassOSSIPIdx:
                        bestMidMassOSSIPIdx = len(ZLLs)
                        select = True
                    if iZLL == bestLowMassOSSIPIdx:
                        bestLowMassOSSIPIdx = len(ZLLs)
                        select = True
                    if iZLL == bestHighMassSSSIPIdx:
                        bestHighMassSSSIPIdx = len(ZLLs)
                        select = True
                    if iZLL == bestMidMassSSSIPIdx:
                        bestMidMassSSSIPIdx = len(ZLLs)
                        select = True
                    if iZLL == bestLowMassSSSIPIdx:
                        bestLowMassSSSIPIdx = len(ZLLs)
                        select = True
                    if iZLL == bestHighMassSSRelaxedIdx:
                        bestHighMassSSRelaxedIdx = len(ZLLs)
                        select = True

                    if select : ZLLs.append(ZLL)
                    if self.DEBUG: print("ZLL:", iZLL, ZLL.p4.M(), ZLL.Z1.M, ZLL.Z2.M, ZLL.Z2.sumpt(), ZLL.finalState(), ZLL.p_GG_SIG_ghg2_1_ghz1_1_JHUGen, ZLL.p_QQB_BKG_MCFM, ZLL.KD,
                                         "2P2F:", int(iZLL==best2P2FCRIdx), "3P1F:", int(iZLL==best3P1FCRIdx), "SS:", int(iZLL == bestSSCRIdx), "SIP:", int(iZLL == bestSIPCRIdx))
                    
            if self.candsToStore == self.BestCandOnly : # keep only the best cand as single element of the ZZ collection
                if bestCandIdx >= 0 :
                    ZZs = [ZZs[bestCandIdx]]
                    bestCandIdx = 0
                else :
                    ZZs = []

        ### Z+L CR, for fake rate. This is considered only for events with a Z + exactly 1 additional lepton passing the relaxed selection.
        ### This ensures that there is no overlap with the SR and OS, 3P1F and 2P2F CRs (there can be an overlap with the SIP CR
        ### as that keeps leptons failing SIP)
        if self.addZLCR and bestZIdx >= 0 and selZs[bestZIdx].M > 40 and selZs[bestZIdx].M < 120:
            aZ = selZs[bestZIdx]
            for i,aL in enumerate(leps):
                # Search for additional lepton, with ghost suppression DR cut
                if i != aZ.l1Idx and i!= aZ.l2Idx and aL.ZZRelaxedId and \
                   deltaR(aL.eta, aL.phi, aZ.l1.eta, aZ.l1.phi) > 0.02 and \
                   deltaR(aL.eta, aL.phi, aZ.l2.eta, aZ.l2.phi) > 0.02 :
                    if ZLCand_lIdx < 0 :
                        ZLCand_lIdx = i
                    else : # more than 1 additional lepton, drop the CR
                        ZLCand_lIdx = -1
                        break
            if ZLCand_lIdx >= 0 :
                aL = leps[ZLCand_lIdx]
                #Apply QCD suppression cut (mLL>4 cut on all OS pairs) 
                if (aL.charge != aZ.l1.charge and (aL.p4()+aZ.l1.p4()).M() <= 4) or \
                   (aL.charge != aZ.l2.charge and (aL.p4()+aZ.l2.p4()).M() <= 4) :
                    ZLCand_lIdx = -1
                        
        ### Filter events with no candidates
        if self.filterType == self.filterOnCands and len(ZZs) == 0 and len(ZLLs) == 0 and ZLCand_lIdx < 0: return False
        if self.filterType == self.filterOnZ and len(selZs) == 0 : return False

        ### Now fill the variables to be stored as output
        # Fill selected Zs
        ZCand_mass = [0.]*len(selZs)
        ZCand_pt = [0.]*len(selZs)
        ZCand_eta = [0.]*len(selZs)
        ZCand_rapidity = [0.]*len(selZs)
        ZCand_phi = [0.]*len(selZs)
        ZCand_flav = [0.]*len(selZs)
        ZCand_l1Idx = [-1]*len(selZs)
        ZCand_l2Idx = [-1]*len(selZs)
        ZCand_fsr1Idx = [-1]*len(selZs)
        ZCand_fsr2Idx = [-1]*len(selZs)

        for iZ, aZ in enumerate(selZs) :
            ZCand_mass[iZ] = aZ.p4.M()
            ZCand_pt[iZ] = aZ.p4.Pt()
            ZCand_eta[iZ] = aZ.p4.Eta()
            ZCand_rapidity[iZ] = aZ.p4.Rapidity()
            ZCand_phi[iZ] = aZ.p4.Phi()
            ZCand_flav[iZ] = aZ.finalState()
            ZCand_l1Idx[iZ] = aZ.l1Idx
            ZCand_l2Idx[iZ] = aZ.l2Idx
            ZCand_fsr1Idx[iZ] = aZ.fsr1Idx
            ZCand_fsr2Idx[iZ] = aZ.fsr2Idx

        self.out.fillBranch("ZCand_mass", ZCand_mass)
        self.out.fillBranch("ZCand_pt", ZCand_pt)
        self.out.fillBranch("ZCand_eta", ZCand_eta)
        self.out.fillBranch("ZCand_rapidity", ZCand_rapidity)
        self.out.fillBranch("ZCand_phi", ZCand_phi)
        self.out.fillBranch("ZCand_flav", ZCand_flav)
        self.out.fillBranch("ZCand_l1Idx", ZCand_l1Idx)
        self.out.fillBranch("ZCand_l2Idx", ZCand_l2Idx)
        self.out.fillBranch("ZCand_fsr1Idx", ZCand_fsr1Idx)
        self.out.fillBranch("ZCand_fsr2Idx", ZCand_fsr2Idx)
        self.out.fillBranch("bestZIdx", bestZIdx)

        ### Fill ZZ candidates
        ZZCand_mass = [0.]*len(ZZs)
        ZZCand_massPreFSR = [0.]*len(ZZs)
        ZZCand_pt = [0.]*len(ZZs)
        ZZCand_eta = [0.]*len(ZZs)
        ZZCand_rapidity = [0.]*len(ZZs)
        ZZCand_phi = [0.]*len(ZZs)
        ZZCand_Z1mass = [0.]*len(ZZs)
        ZZCand_Z1flav = [0.]*len(ZZs)
        ZZCand_Z1pt   = [0.]*len(ZZs)
        ZZCand_Z1eta  = [0.]*len(ZZs)
        ZZCand_Z1phi  = [0.]*len(ZZs)
        ZZCand_Z1rapidity = [0.]*len(ZZs)
        ZZCand_Z2mass = [0.]*len(ZZs)
        ZZCand_Z2flav = [0.]*len(ZZs)
        ZZCand_Z2pt   = [0.]*len(ZZs)
        ZZCand_Z2eta  = [0.]*len(ZZs)
        ZZCand_Z2phi  = [0.]*len(ZZs)
        ZZCand_Z2rapidity = [0.]*len(ZZs)
        ZZCand_Z1l1Idx = [-1]*len(ZZs)
        ZZCand_Z1l2Idx = [-1]*len(ZZs)
        ZZCand_Z2l1Idx = [-1]*len(ZZs)
        ZZCand_Z2l2Idx = [-1]*len(ZZs)
        ZZCand_KD = [0.]*len(ZZs)
        ZZCand_Z2sumpt = [0.]*len(ZZs)

        ZZCand_cosTheta1    = [0.]*len(ZZs)
        ZZCand_cosTheta3    = [0.]*len(ZZs)
        ZZCand_cosThetaStar = [0.]*len(ZZs)
        ZZCand_delPhi       = [0.]*len(ZZs)
        ZZCand_delPhiStar   = [0.]*len(ZZs)
        ZZCand_delRapidity  = [0.]*len(ZZs)

        for iZZ, ZZ in enumerate(ZZs) :
            ang_vars = AngularVars(ZZ)

            ZZCand_cosTheta1[iZZ]    = ang_vars.cosTheta("1")
            ZZCand_cosTheta3[iZZ]    = ang_vars.cosTheta("2")
            ZZCand_cosThetaStar[iZZ] = ang_vars.cosThetaStar()
            ZZCand_delPhiStar[iZZ]   = ang_vars.delPhiStar()
            ZZCand_delPhi[iZZ]       = ang_vars.delPhi()
            ZZCand_delRapidity[iZZ]  = ang_vars.delRapidity()

            ZZCand_mass[iZZ] = ZZ.p4.M()
            ZZCand_massPreFSR[iZZ] = ZZ.massPreFSR()
            ZZCand_pt[iZZ] = ZZ.p4.Pt()
            ZZCand_eta[iZZ] = ZZ.p4.Eta()
            ZZCand_rapidity[iZZ] = ZZ.p4.Rapidity()
            ZZCand_phi[iZZ] = ZZ.p4.Phi()
            ZZCand_Z1mass[iZZ] = ZZ.Z1.M
            ZZCand_Z1flav[iZZ] = ZZ.Z1.finalState()
            ZZCand_Z1pt[iZZ] = ZZ.Z1.p4.Pt()
            ZZCand_Z1eta[iZZ] = ZZ.Z1.p4.Eta()
            ZZCand_Z1phi[iZZ] = ZZ.Z1.p4.Phi()
            ZZCand_Z1rapidity[iZZ] = ZZ.Z1.p4.Rapidity()
            ZZCand_Z2mass[iZZ] = ZZ.Z2.M
            ZZCand_Z2flav[iZZ] = ZZ.Z2.finalState()
            ZZCand_Z2pt[iZZ] = ZZ.Z2.p4.Pt()
            ZZCand_Z2eta[iZZ] = ZZ.Z2.p4.Eta()
            ZZCand_Z2phi[iZZ] = ZZ.Z2.p4.Phi()
            ZZCand_Z2rapidity[iZZ] = ZZ.Z2.p4.Rapidity()
            ZZCand_Z1l1Idx[iZZ] = ZZ.Z1.l1Idx
            ZZCand_Z1l2Idx[iZZ] = ZZ.Z1.l2Idx
            ZZCand_Z2l1Idx[iZZ] = ZZ.Z2.l1Idx
            ZZCand_Z2l2Idx[iZZ] = ZZ.Z2.l2Idx
            ZZCand_KD[iZZ] = ZZ.KD
            ZZCand_Z2sumpt[iZZ] = ZZ.Z2.sumpt()

        self.out.fillBranch("ZZCand_mass", ZZCand_mass)
        self.out.fillBranch("ZZCand_massPreFSR", ZZCand_massPreFSR)
        self.out.fillBranch("ZZCand_pt", ZZCand_pt)
        self.out.fillBranch("ZZCand_eta", ZZCand_eta)
        self.out.fillBranch("ZZCand_rapidity", ZZCand_rapidity)
        self.out.fillBranch("ZZCand_phi", ZZCand_phi)
        self.out.fillBranch("ZZCand_Z1mass", ZZCand_Z1mass)
        self.out.fillBranch("ZZCand_Z1flav", ZZCand_Z1flav)
        self.out.fillBranch("ZZCand_Z1pt", ZZCand_Z1pt)
        self.out.fillBranch("ZZCand_Z1eta", ZZCand_Z1eta)
        self.out.fillBranch("ZZCand_Z1phi", ZZCand_Z1phi)
        self.out.fillBranch("ZZCand_Z1rapidity", ZZCand_Z1rapidity)
        self.out.fillBranch("ZZCand_Z2mass", ZZCand_Z2mass)
        self.out.fillBranch("ZZCand_Z2flav", ZZCand_Z2flav)
        self.out.fillBranch("ZZCand_Z2pt", ZZCand_Z2pt)
        self.out.fillBranch("ZZCand_Z2eta", ZZCand_Z2eta)
        self.out.fillBranch("ZZCand_Z2phi", ZZCand_Z2phi)
        self.out.fillBranch("ZZCand_Z2rapidity", ZZCand_Z2rapidity)
        self.out.fillBranch("ZZCand_KD", ZZCand_KD)
        self.out.fillBranch("ZZCand_Z2sumpt", ZZCand_Z2sumpt)
        self.out.fillBranch("ZZCand_Z1l1Idx", ZZCand_Z1l1Idx)
        self.out.fillBranch("ZZCand_Z1l2Idx", ZZCand_Z1l2Idx)
        self.out.fillBranch("ZZCand_Z2l1Idx", ZZCand_Z2l1Idx)
        self.out.fillBranch("ZZCand_Z2l2Idx", ZZCand_Z2l2Idx)

        self.out.fillBranch("ZZCand_cosTheta1", ZZCand_cosTheta1)
        self.out.fillBranch("ZZCand_cosTheta3", ZZCand_cosTheta3)
        self.out.fillBranch("ZZCand_cosThetaStar", ZZCand_cosThetaStar)
        self.out.fillBranch("ZZCand_delPhiStar", ZZCand_delPhiStar)
        self.out.fillBranch("ZZCand_delPhi", ZZCand_delPhi)
        self.out.fillBranch("ZZCand_delRapidity", ZZCand_delRapidity)

        self.out.fillBranch("bestCandIdx", bestCandIdx)

        if self.addSSCR or self. addOSCR or self.addSIPCR :
            ZLLCand_mass   = [0.]*len(ZLLs)
            ZLLCand_massPreFSR = [0.]*len(ZLLs)
            ZLLCand_pt     = [0.]*len(ZLLs)
            ZLLCand_eta    = [0.]*len(ZLLs)
            ZLLCand_rapidity = [0.]*len(ZLLs)
            ZLLCand_phi    = [0.]*len(ZLLs)
            ZLLCand_Z1mass = [0.]*len(ZLLs)
            ZLLCand_Z1flav = [0.]*len(ZLLs)
            ZLLCand_Z1pt   = [0.]*len(ZLLs)
            ZLLCand_Z1eta  = [0.]*len(ZLLs)
            ZLLCand_Z1phi  = [0.]*len(ZLLs)
            ZLLCand_Z1rapidity = [0.]*len(ZLLs)
            ZLLCand_Z2mass = [0.]*len(ZLLs)
            ZLLCand_Z2flav = [0.]*len(ZLLs)
            ZLLCand_Z2pt   = [0.]*len(ZLLs)
            ZLLCand_Z2eta  = [0.]*len(ZLLs)
            ZLLCand_Z2phi  = [0.]*len(ZLLs)
            ZLLCand_Z2rapidity = [0.]*len(ZLLs)
            ZLLCand_Z1l1Idx = [-1]*len(ZLLs)
            ZLLCand_Z1l2Idx = [-1]*len(ZLLs)
            ZLLCand_Z2l1Idx = [-1]*len(ZLLs)
            ZLLCand_Z2l2Idx = [-1]*len(ZLLs)
            ZLLCand_KD     = [0.]*len(ZLLs)

            ZLLCand_cosTheta1    = [0.]*len(ZLLs)
            ZLLCand_cosTheta3    = [0.]*len(ZLLs)
            ZLLCand_cosThetaStar = [0.]*len(ZLLs)
            ZLLCand_delPhi       = [0.]*len(ZLLs)
            ZLLCand_delPhiStar   = [0.]*len(ZLLs)
            ZLLCand_delRapidity  = [0.]*len(ZLLs)

            for iZLL, ZLL in enumerate(ZLLs) :
                ang_vars = AngularVars(ZLL)

                ZLLCand_cosTheta1[iZLL]    = ang_vars.cosTheta("1")
                ZLLCand_cosTheta3[iZLL]    = ang_vars.cosTheta("2")
                ZLLCand_cosThetaStar[iZLL] = ang_vars.cosThetaStar()
                ZLLCand_delPhiStar[iZLL]   = ang_vars.delPhiStar()
                ZLLCand_delPhi[iZLL]       = ang_vars.delPhi()
                ZLLCand_delRapidity[iZLL]  = ang_vars.delRapidity()

                ZLLCand_mass[iZLL] = ZLL.p4.M()
                ZLLCand_massPreFSR[iZLL] = ZLL.massPreFSR()
                ZLLCand_pt[iZLL] = ZLL.p4.Pt()
                ZLLCand_eta[iZLL] = ZLL.p4.Eta()
                ZLLCand_rapidity[iZLL] = ZLL.p4.Rapidity()
                ZLLCand_phi[iZLL] = ZLL.p4.Phi()
                ZLLCand_Z1mass[iZLL] = ZLL.Z1.M
                ZLLCand_Z1flav[iZLL] = ZLL.Z1.finalState()
                ZLLCand_Z1pt[iZLL] = ZLL.Z1.p4.Pt()
                ZLLCand_Z1eta[iZLL] = ZLL.Z1.p4.Eta()
                ZLLCand_Z1phi[iZLL] = ZLL.Z1.p4.Phi()
                ZLLCand_Z1rapidity[iZLL] = ZLL.Z1.p4.Rapidity()
                ZLLCand_Z2mass[iZLL] = ZLL.Z2.M
                ZLLCand_Z2flav[iZLL] = ZLL.Z2.finalState()
                ZLLCand_Z2pt[iZLL] = ZLL.Z2.p4.Pt()
                ZLLCand_Z2eta[iZLL] = ZLL.Z2.p4.Eta()
                ZLLCand_Z2phi[iZLL] = ZLL.Z2.p4.Phi()
                ZLLCand_Z2rapidity[iZLL] = ZLL.Z2.p4.Rapidity()
                ZLLCand_Z1l1Idx[iZLL] = ZLL.Z1.l1Idx
                ZLLCand_Z1l2Idx[iZLL] = ZLL.Z1.l2Idx
                ZLLCand_Z2l1Idx[iZLL] = ZLL.Z2.l1Idx
                ZLLCand_Z2l2Idx[iZLL] = ZLL.Z2.l2Idx
                ZLLCand_KD[iZLL] = ZLL.KD

            self.out.fillBranch("ZLLCand_mass",   ZLLCand_mass)
            self.out.fillBranch("ZLLCand_massPreFSR",   ZLLCand_massPreFSR)
            self.out.fillBranch("ZLLCand_pt",     ZLLCand_pt)
            self.out.fillBranch("ZLLCand_eta",    ZLLCand_eta)
            self.out.fillBranch("ZLLCand_rapidity",    ZLLCand_rapidity)
            self.out.fillBranch("ZLLCand_phi",    ZLLCand_phi)
            self.out.fillBranch("ZLLCand_Z1mass", ZLLCand_Z1mass)
            self.out.fillBranch("ZLLCand_Z1flav", ZLLCand_Z1flav)
            self.out.fillBranch("ZLLCand_Z1pt", ZLLCand_Z1pt)
            self.out.fillBranch("ZLLCand_Z1eta", ZLLCand_Z1eta)
            self.out.fillBranch("ZLLCand_Z1phi", ZLLCand_Z1phi)
            self.out.fillBranch("ZLLCand_Z1rapidity", ZLLCand_Z1rapidity)
            self.out.fillBranch("ZLLCand_Z2mass", ZLLCand_Z2mass)
            self.out.fillBranch("ZLLCand_Z2flav", ZLLCand_Z2flav)
            self.out.fillBranch("ZLLCand_Z2pt", ZLLCand_Z2pt)
            self.out.fillBranch("ZLLCand_Z2eta", ZLLCand_Z2eta)
            self.out.fillBranch("ZLLCand_Z2phi", ZLLCand_Z2phi)
            self.out.fillBranch("ZLLCand_Z2rapidity", ZLLCand_Z2rapidity)
            self.out.fillBranch("ZLLCand_Z1l1Idx", ZLLCand_Z1l1Idx)
            self.out.fillBranch("ZLLCand_Z1l2Idx", ZLLCand_Z1l2Idx)
            self.out.fillBranch("ZLLCand_Z2l1Idx", ZLLCand_Z2l1Idx)
            self.out.fillBranch("ZLLCand_Z2l2Idx", ZLLCand_Z2l2Idx)
            self.out.fillBranch("ZLLCand_KD",     ZLLCand_KD)

            self.out.fillBranch("ZLLCand_cosTheta1", ZLLCand_cosTheta1)
            self.out.fillBranch("ZLLCand_cosTheta3", ZLLCand_cosTheta3)
            self.out.fillBranch("ZLLCand_cosThetaStar", ZLLCand_cosThetaStar)
            self.out.fillBranch("ZLLCand_delPhiStar", ZLLCand_delPhiStar)
            self.out.fillBranch("ZLLCand_delPhi", ZLLCand_delPhi)
            self.out.fillBranch("ZLLCand_delRapidity", ZLLCand_delRapidity)

            if self.addSSCR :
                self.out.fillBranch("ZLLbestSSIdx",  bestSSCRIdx)
            if self.addOSCR :
                self.out.fillBranch("ZLLbest2P2FIdx", best2P2FCRIdx)
                self.out.fillBranch("ZLLbest3P1FIdx", best3P1FCRIdx)
            if self.addSIPCR :
                self.out.fillBranch("ZLLbestSIPCRIdx", bestSIPCRIdx)
                self.out.fillBranch("ZLLbestHighMassOSSIPIdx", bestHighMassOSSIPIdx)
                self.out.fillBranch("ZLLbestMidMassOSSIPIdx", bestMidMassOSSIPIdx)
                self.out.fillBranch("ZLLbestLowMassOSSIPIdx", bestLowMassOSSIPIdx)
                self.out.fillBranch("ZLLbestHighMassSSSIPIdx", bestHighMassSSSIPIdx)
                self.out.fillBranch("ZLLbestMidMassSSSIPIdx", bestMidMassSSSIPIdx)
                self.out.fillBranch("ZLLbestLowMassSSSIPIdx", bestLowMassSSSIPIdx)
                self.out.fillBranch("ZLLbestHighMassSSRelaxedIdx", bestHighMassSSRelaxedIdx)
        if self.addZLCR :
            self.out.fillBranch("ZLCand_lepIdx", ZLCand_lIdx)

        ### Fill control plot (example)
        # self.h_ZZMass.Fill(ZZCand_mass[bestCandIdx])

        return True


    # Temporary class to store information on a ZZ candidate.
    # NOTE: We may need to move Zs into the Event as persistent objects, and to be re-used for CRs. They would be built by a separate module in that case.
    class ZCand: 
        def __init__(self, l1Idx, l2Idx, leps, fsrPhotons):
            self.l1Idx = l1Idx
            self.l2Idx = l2Idx
            self.l1 = leps[l1Idx]
            self.l2 = leps[l2Idx]
            self.fsr1Idx = self.l1.fsrPhotonIdx
            self.fsr2Idx = self.l2.fsrPhotonIdx
    
            self.l1DressedP4 = self.l1.p4()
            self.l2DressedP4 = self.l2.p4()
            if self.fsr1Idx>=0 : self.l1DressedP4 += fsrPhotons[self.fsr1Idx].p4()
            if self.fsr2Idx>=0 : self.l2DressedP4 += fsrPhotons[self.fsr2Idx].p4()
    
            self.p4 = self.l1DressedP4 + self.l2DressedP4
    
            self.M = self.p4.M() # cache the mass as it is used often
    
        def sumpt(self) : # sum of lepton pTs, used to sort candidates
            return self.l1.pt+self.l2.pt
    
        def finalState(self) :
            return self.l1.pdgId*self.l2.pdgId
    
    
    # Temporary class to store information on a ZZ candidate.
    class ZZCand:
        def __init__(self, Z1, Z2, p_GG_SIG_ghg2_1_ghz1_1_JHUGen=0., p_QQB_BKG_MCFM=1.):
            self.Z1 = Z1
            self.Z2 = Z2
            self.p4 = Z1.p4+Z2.p4
            self.M  = self.p4.M()
            self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen = p_GG_SIG_ghg2_1_ghz1_1_JHUGen
            self.p_QQB_BKG_MCFM = p_QQB_BKG_MCFM
            self.KD = p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(p_GG_SIG_ghg2_1_ghz1_1_JHUGen+p_QQB_BKG_MCFM) # without c-constants, for candidate sorting

        def finalState(self) :
            return self.Z1.finalState()*self.Z2.finalState()

        def massPreFSR(self) :
            return (self.Z1.l1.p4()+self.Z1.l2.p4()+self.Z2.l1.p4()+self.Z2.l2.p4()).M()

        def leps(self) :
            return([self.Z1.l1, self.Z1.l2, self.Z2.l1, self.Z2.l2])


    ### Comparators to select the best candidate in the event. Return -1 if a is better than b, +1 otherwise
    # Choose by abs(MZ1-MZ), or sum(PT) if same Z1
    def bestCandByZ1Z2(self,a,b): 
        if a.Z1.l1Idx == b.Z1.l1Idx and a.Z1.l2Idx == b.Z1.l2Idx :
            # Same Z1, choose by sum of Z2 pTs.
            # Note that leptons are ordered (1=+, 2=-), there is no need to check the alternative pairing
            if a.Z2.sumpt() > b.Z2.sumpt() :
                return -1
            else :
                return 1
        else : # choose based on Z1 masses
            if abs(a.Z1.M-self.ZmassValue) < abs(b.Z1.M-self.ZmassValue) :
                return -1 
            else :
                return 1

    # Choose by DbkgKin
    def bestCandByDbkgKin(self,a,b): 
        if set([a.Z1.l1Idx, a.Z1.l2Idx, a.Z2.l1Idx, a.Z2.l2Idx]) == \
           set([b.Z1.l1Idx, b.Z1.l2Idx, b.Z2.l1Idx, b.Z2.l2Idx]) :
            # Equivalent: same masss (tolerance 100 keV) and same FS -> different permutation of the same leptons.
            # Note that this can only happen in SR, not in CRs where the Z1 is always the best Z in the event.
            return self.bestCandByZ1Z2(a,b)
        if a.KD > b.KD : return -1 # choose by best dbkgkin
        else: return 1
        

    def makeCand(self, Za, Zb, sortZsByMass=True, fillIDVars=True) :
        '''Build a ZZ object from given Za, Zb pair, if it passes selection cuts; None is returned otherwise.
        All relevant candidate variables are computed for the candidate. Options:
        sortZsByMass : set Z1 and Z2 according to closest-Mz criteria (for SR); otherwise, specified order is
                       kept (useful for CRs)
        '''
        
        # check that these Zs are mutually exclusive (not sharing the same lepton) 
        if Za.l1Idx==Zb.l1Idx or Za.l2Idx==Zb.l2Idx or Za.l2Idx==Zb.l1Idx or Za.l2Idx==Zb.l2Idx: return None
        
        # set Z1 and Z2
        Z1, Z2 = Za, Zb
        if sortZsByMass and abs(Zb.M-self.ZmassValue) < abs(Za.M-self.ZmassValue):
            Z1, Z2 = Z2, Z1

        # Z1 mass cut (regardless of region!)
        if Z1.M <= 40. : return None


        zzleps = [Z1.l1, Z1.l2, Z2.l1, Z2.l2]
        lepPts = []
        # QCD suppression on all OS pairs, regardelss of flavour
        passQCD = True    # QCD suppression on all OS pairs, regardelss of flavour
        passDeltaR = True # DeltaR>0.02 cut among all leptons to protect against split tracks
        for k in range(4):
            lepPts.append(zzleps[k].pt)
            for l in range (k+1,4):
                if zzleps[k].charge!=zzleps[l].charge and (zzleps[k].p4()+zzleps[l].p4()).M()<=4.:
                    passQCD = False
                    break
                if deltaR(zzleps[k].eta, zzleps[k].phi, zzleps[l].eta, zzleps[l].phi) <= 0.02 :
                    passDeltaR = False
                    break

        if self.DEBUG : print(f"ZZ: Z1: {Z1.M}, Z2: {Z2.M}, pTs: {lepPts}, passDR: {passDeltaR} passQCD: {passQCD}")
        if not (passQCD and passDeltaR) : return None

        # trigger acceptance cuts (20,10 GeV)
        lepPts.sort()
        if not (lepPts[3]>20. and lepPts[2]>10.) : return None

        #Compute D_bkg^kin
        p_GG_SIG_ghg2_1_ghz1_1_JHUGen = 0.
        p_QQB_BKG_MCFM = 1.
        if self.mela != None:
            daughters = Mela.SimpleParticleCollection_t()
            daughters.add_particle(Mela.SimpleParticle_t(Z1.l1.pdgId, Z1.l1DressedP4.Px(), Z1.l1DressedP4.Py(), Z1.l1DressedP4.Pz(), Z1.l1DressedP4.E()))
            daughters.add_particle(Mela.SimpleParticle_t(Z1.l2.pdgId, Z1.l2DressedP4.Px(), Z1.l2DressedP4.Py(), Z1.l2DressedP4.Pz(), Z1.l2DressedP4.E()))
            daughters.add_particle(Mela.SimpleParticle_t(Z2.l1.pdgId, Z2.l1DressedP4.Px(), Z2.l1DressedP4.Py(), Z2.l1DressedP4.Pz(), Z2.l1DressedP4.E()))
            daughters.add_particle(Mela.SimpleParticle_t(Z2.l2.pdgId, Z2.l2DressedP4.Px(), Z2.l2DressedP4.Py(), Z2.l2DressedP4.Pz(), Z2.l2DressedP4.E()))
            self.mela.setInputEvent(daughters, None, None, 0)
            self.mela.setProcess(Mela.Process.HSMHiggs, Mela.MatrixElement.JHUGen, Mela.Production.ZZGG)

            p_GG_SIG_ghg2_1_ghz1_1_JHUGen = self.mela.computeP(True)

            self.mela.setProcess(Mela.Process.bkgZZ, Mela.MatrixElement.MCFM, Mela.Production.ZZQQB)

            p_QQB_BKG_MCFM = self.mela.computeP(True)

            self.mela.resetInputEvent()

        if (p_GG_SIG_ghg2_1_ghz1_1_JHUGen+p_QQB_BKG_MCFM == 0.) :
            print ("ERROR", p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM)
            p_QQB_BKG_MCFM = 1. # FIXME: fix for error with message: "TUtil::CheckPartonMomFraction: At least one of the parton momentum fractions is greater than 1."
        ZZ = self.ZZCand(Z1, Z2, p_GG_SIG_ghg2_1_ghz1_1_JHUGen, p_QQB_BKG_MCFM)

        ZZ.HighMass = ZZ.M > 180
        ZZ.MidMass  = (ZZ.M > 140) & (ZZ.M < 180)
        ZZ.LowMass  = (ZZ. M > 105) & (ZZ.M < 140)

        if not (ZZ.HighMass or ZZ.MidMass or ZZ.LowMass):
            return None
        
        return ZZ
