from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR

# For debugging
from ZZAnalysis.NanoAnalysis.tools import Mother, getParentID

#from ROOT import Mela, SimpleParticle_t, SimpleParticleCollection_t, TVar, TLorentzVector, TMath, LeptonSFHelper
from ROOT import LeptonSFHelper
from ROOT.Math import LorentzVector, PxPyPzE4D, PtEtaPhiM4D, Boost, LorentzRotation
import ROOT

from ctypes import c_float

from itertools import combinations

import csv
import numpy as np
from json import dump

class StoreOption:
    # Enum used to define which candidates should be stored in the ZZCand collection.
    # BestCandOnly = only the best SR candidate in the event is retained. This is normally the default for production jobs.
    # AllCands   = keep all SR candidates passing the full selection and analysis cuts (including permutations of leptons).
    # AllWithRelaxedMuId = keep any SR candidate that can be made, even if leptons don't pass ID cuts (useful for ID cut optimization studies).
    # Note: This does not affect the ZLLCand collection; for each CR that is activated, only the best candidate is stored.
    BestCandOnly, AllCands, AllWithRelaxedMuId = range(0,3)

'''class AngularVars:
    def __init__(self, ZZ):
        self.ZZ = ZZ
        self.Z  = {"1": ZZ.Z1, "2": ZZ.Z2}

    def _boost_to_cm(self, v1, v2):
        """Returns v1 4-vector in the v2 rest frame."""
        boosted_v1 = v1.Clone()
        boosted_v1.Boost(-v2.BoostVector())
        return boosted_v1

    def cosTheta(self, z):
        """Return cos(θ) where θ is defined
        as the angle between the negatively
        charged lepton in its parent Z's rest
        frame and the Z in the four lepton rest
        frame."""

        Z = self.Z[z]

        lep_p4_zRest = self._boost_to_cm(Z.l2DressedP4, Z.p4)
        z_p4_m4lRest = self._boost_to_cm(Z.p4, self.ZZ.p4)

        return TMath.Cos(lep_p4_zRest.Angle(z_p4_m4lRest.Vect()))

    def cosThetaStar(self):
        """Return cos(θ*) where θ* is defined
        as the production angle of the Z1 in the four
        lepton rest frame."""

        z_p4_m4lRest = self._boost_to_cm(self.ZZ.Z1.p4, self.ZZ.p4)

        return z_p4_m4lRest.CosTheta()

    def delPhi(self, z):
        # Verify this!!
        """Return Δ(φ_12) where φ_12 is defined
        as the azimuthal separation of the two
        leptons from the Z defined in the four-lepton
        rest frame."""

        Z = self.Z[z]        

        lep1_p4_m4lRest = self._boost_to_cm(Z.l1DressedP4, self.ZZ.p4)
        lep2_p4_m4lRest = self._boost_to_cm(Z.l2DressedP4, self.ZZ.p4)

        return lep1_p4_m4lRest.DeltaPhi(lep2_p4_m4lRest)'''

class lheParts:
    def __init__(self, lhepart):
        self.lhepart = lhepart
        self.z1_flav = 0
        self.z2_flav = 0
        self.zz_flav = 0

    def _getLorentzVec(self, idx):
        return LorentzVector(PtEtaPhiM4D('double'))(
            self.lhepart[idx].pt,
            self.lhepart[idx].eta,
            self.lhepart[idx].phi,
            self.lhepart[idx].mass
        )

    def _setInfo(self):
        z1 = self._getLorentzVec(2)
        z2 = self._getLorentzVec(3)

        zz = z1+z2

        ll_z1 = self._getLorentzVec(4) + self._getLorentzVec(5)
        ll_z2 = self._getLorentzVec(6) + self._getLorentzVec(7)

        if abs(z1.M() - ll_z1.M()) > abs(z1.M() - ll_z2.M()):
            ll_z1, ll_z2 = ll_z2, ll_z1
            self.z1_flav = self.lhepart[6].pdgId*self.lhepart[7].pdgId
            self.z2_flav = self.lhepart[4].pdgId*self.lhepart[5].pdgId

            self.lep_1_pt  = self.lhepart[6].pt
            self.lep_1_eta = self.lhepart[6].eta

            self.lep_2_pt  = self.lhepart[7].pt
            self.lep_2_eta = self.lhepart[7].eta

            self.lep_3_pt  = self.lhepart[4].pt
            self.lep_3_eta = self.lhepart[4].eta

            self.lep_4_pt  = self.lhepart[5].pt
            self.lep_4_eta = self.lhepart[5].eta
        else:
            self.z1_flav = self.lhepart[4].pdgId*self.lhepart[5].pdgId
            self.z2_flav = self.lhepart[6].pdgId*self.lhepart[7].pdgId

            self.lep_1_pt  = self.lhepart[4].pt
            self.lep_1_eta = self.lhepart[4].eta

            self.lep_2_pt  = self.lhepart[5].pt
            self.lep_2_eta = self.lhepart[5].eta

            self.lep_3_pt  = self.lhepart[6].pt
            self.lep_3_eta = self.lhepart[6].eta

            self.lep_4_pt  = self.lhepart[7].pt
            self.lep_4_eta = self.lhepart[7].eta

        self.zz_flav = self.z1_flav * self.z2_flav
        self.has_tau = self.zz_flav % 15**2 == 0

        self.z1_mass = z1.M()
        self.z2_mass = z2.M()
        self.zz_mass = zz.M()

    def printInfo(self):
        self._setInfo()
        print("\nLHE Info:")
        print("===========================")
        print("\nZZ Flavor: ", self.zz_flav)
        print("\nZ1 Flavor: ", self.z1_flav)
        print("\nZ2 Flavor: ", self.z2_flav)
        print("===========================")
        print("\nZZ Mass: ", self.zz_mass)
        print("\nZ1 Mass: ", self.z1_mass)
        print("\nZ2 Mass: ", self.z2_mass)
        print("===========================")
        print("\nLep 1 pT: ", self.lep_1_pt)
        print("\nLep 2 pT: ", self.lep_2_pt)
        print("\nLep 3 pT: ", self.lep_3_pt)
        print("\nLep 4 pT: ", self.lep_4_pt)
        print("===========================")
        print("\nLep 1 eta: ", self.lep_1_eta)
        print("\nLep 2 eta: ", self.lep_2_eta)
        print("\nLep 3 eta: ", self.lep_3_eta)
        print("\nLep 4 eta: ", self.lep_4_eta)

class AngularVars:
    def __init__(self, ZZ):
        self.ZZ = ZZ
        self.Z  = {"1": ZZ.Z1, "2": ZZ.Z2}

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

    def _getBoost(self, p4):
        return Boost(p4.BoostToCM())

    def cosTheta(self, z):
        lab_z_cand = self.Z[z]
        
        # Four vectors in the lab frame for 4l system, Z, and its associated l+
        lab_4l_p4   = self._getLorentzVec(self.ZZ.p4)
        lab_z_p4    = self._getLorentzVec(lab_z_cand.p4)
        
        lab_pos_lep = self._getPosLep(lab_z_cand)

        # Boost vectors to 4l rest frame and Z rest frame
        boost_to_4l  = self._getBoost(lab_4l_p4)
        boost_to_ZCM = self._getBoost(lab_z_p4)

        zcand_4l_com  = boost_to_4l*lab_z_p4
        pos_lep_z_com = boost_to_ZCM*lab_pos_lep

        return ROOT.Math.VectorUtil.CosTheta(pos_lep_z_com, zcand_4l_com)

    def cosThetaStar(self, z="1"):
        lab_z_cand = self.Z[z]

        # Four vectors in the lab frame for 4l system and Z
        lab_4l_p4   = self._getLorentzVec(self.ZZ.p4)
        lab_z_p4    = self._getLorentzVec(lab_z_cand.p4)

        # Boost vectors to 4l rest frame and Z rest frame
        boost_to_4l  = self._getBoost(lab_4l_p4)

        zcand_4l_com = boost_to_4l*lab_z_p4

        return ROOT.Math.VectorUtil.CosTheta(zcand_4l_com, lab_4l_p4)

    def delPhiStar(self):
        lab_z1_cand, lab_z2_cand = self.Z["1"], self.Z["2"]

        lab_pos_lep_z1 = self._getPosLep(lab_z1_cand)
        lab_pos_lep_z2 = self._getPosLep(lab_z2_cand)

        lab_z1_p4 = self._getLorentzVec(lab_z1_cand.p4)
        lab_z2_p4 = self._getLorentzVec(lab_z2_cand.p4)

        boost_to_z1 = self._getBoost(lab_z1_p4)
        boost_to_z2 = self._getBoost(lab_z2_p4)

        pos_lep_z1_cm = boost_to_z1*lab_pos_lep_z1
        pos_lep_z2_cm = boost_to_z2*lab_pos_lep_z2

        az_diff = np.abs(pos_lep_z1_cm.Phi() - pos_lep_z2_cm.Phi())

        return min(az_diff, 2*np.pi - az_diff)

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

    def __init__(self, runMELA, bestCandByMELA, isMC, year, data_tag, processCR=False, addZL=False, filter='Cands', debug=False):
        """Build candidates:
        -ZZCand: SR candidates. The index of the best ZZ candidate in each event is stored as bestCandIdx
        -ZLLCand: CR candidates (SS, 3P1F, 2P2F, SIP CRs, with indices: ZZLLbestSSIdx, ZLLbest3P1FIdx, ZLLbest2P2FIdx, ZLLbestSIPCRIdx)
        -ZLCand: Z+L CR, for fake rate computation (only the index of the additional lepton is stored)
        -ZCand: Zs referenced in the above collections (note that when filter='Cands', are applied, only events with at least one candidate in the above collections are retained, so the tree will not contain all events with one Z)
        Parameters:
          runMELA: compute MELA KD
          bestCandByMELA: True = select best candidate by KD; False = with best Z1/highest-pTsum Z2 
          year: data taking year        
          processCR: add ZLLCand CR collections
          data_tag: subperiod or processing, e.g. "pre_EE" (currently unused)
          addZL: add ZL CR
          filter: criteria to keep or skip events:
                  'Cands' = store only events with at least one ZZ, ZLL, or ZL candidate are kept
                  'Z' = store any event that has at least one good Z candidate (passing the analysis Z selection criteria, and 12<mll<120)
                  '3L_20_10' = filter on events with 3 good leptons, pt1>20, pt2>10; useful for trigger studies
                  'NoFilter' = don't skip events
        """
        print("***ZZFiller: isMC:", isMC, "year:", year, "data_tag:", data_tag, "bestCandByMELA:", bestCandByMELA, "filter:", filter, "- This module filters events.",  flush=True)
        self.writeHistFile = False
        self.isMC = isMC
        self.year = year

        if self.isMC:
            self.lepSFHelper = LeptonSFHelper(data_tag)

        self.addSSCR = processCR
        self.addOSCR = processCR
        self.addSIPCR = processCR
        self.addZLCR = addZL

        self.DATA_TAG = data_tag

        # Translate filter string into an enum, for efficiency
        self.noFilter, self.filterOnCands, self.filterOnZ, self.filter_3L_20_10 = range(0,4)
        filters = {'NoFilter':self.noFilter,
                   'Cands':self.filterOnCands,
                   'Z':self.filterOnZ,
                   '3L_20_10':self.filter_3L_20_10}
        try:
            self.filterType = filters[filter]
        except :
            raise ValueError("ZZFiller: filter =", filter, "not supported")
        self.DEBUG = debug
        self.ZmassValue = 91.1876

        self.candsToStore = StoreOption.BestCandOnly # Only store the best candidate for the SR

        # DeltaR>0.02 cut among all leptons to protect against split tracks
        self.passDeltaR = lambda l1, l2: deltaR(l1.eta, l1.phi, l2.eta, l2.phi) > 0.02

        self.lepsExclusive = lambda z1, z2: z1.l1 != z2.l1 and z1.l1 != z2.l2 and z1.l2 != z2.l1 and z1.l2 != z2.l2

        # Pre-selection of leptons used to reduce combinatorial when building Z and LL candidates.
        # Normally it is the full ID + iso if only the SR is considered, or the relaxed ID if CRs are also filled,
        # but can be relaxed further to also build loose candidates for ID studies (see below).
        # Note that the actual lepton selection cuts for SR and CR are applied later; this preselection only affects what
        # leptons are considered in making the combinatorial (ie processing speed)
        if self.addSIPCR or self.addOSCR or self.addSSCR :
            self.leptonPresel = (lambda l : l.ZZRelaxedIdNoSIP) # minimal selection good for all CRs: no SIP, no ID, no iso
        else : # SR only
            self.leptonPresel = (lambda l : l.ZZFullSel)

        self.muonIDs=[]
        self.muonIDVars=[]
        
        # Use relaxed muon preselection for muon ID studies: fully relax muon ID. Electron ID is unchanged.
        if self.candsToStore == StoreOption.AllWithRelaxedMuId :
            if self.addSIPCR or self.addOSCR or self.addSSCR :
                self.leptonPresel = (lambda l : (abs(l.pdgId)==13 and l.pt>5 and abs(l.eta) < 2.4) or (abs(l.pdgId)==11 and l.ZZRelaxedIdNoSIP))
            else : 
                self.leptonPresel = (lambda l : (abs(l.pdgId)==13 and l.pt>5 and abs(l.eta) < 2.4) or (abs(l.pdgId)==11 and l.ZZFullSel))
            
            # Add flags for muon ID studies. Each Flag will be set to true for a candidate if all of its muons pass the specified ID.
            self.muonIDs=[dict(name="ZZFullSel", sel=lambda l : l.ZZFullId and l.passIso), # Standard ZZ selection; this is used for setting default bestCandIdx
                          dict(name="ZZRelaxedIDNoSIP",sel=lambda l : l.pt>5 and abs(l.eta)<2.4 and (l.isGlobal or (l.isTracker and l.nStations>0))),# ZZ relaxed mu ID without dxy, dz, SIP cuts (for optimization). Note: this is looser than nanoAOD presel.
                          dict(name="ZZFullIDNoSIP",   sel=lambda l : l.pt>5 and abs(l.eta)<2.4 and (l.isGlobal or (l.isTracker and l.nStations>0)) and (l.isPFcand or (l.highPtId>0 and l.pt>200.))),# ZZ full ID without dxy, dz, SIP, and isolation cuts (for optimization)                      
                          dict(name="looseId", sel=lambda l : l.looseId),   # POG CutBasedIdLoose
                          dict(name="mediumId", sel=lambda l : l.mediumId), # POG CutBasedIdMedium
                          dict(name="mediumPromptId", sel=lambda l : l.mediumPromptId), # POG CutBasedIdMediumPrompt (=mediumId + tighter dxy, dz cuts)
                          dict(name="tightId", sel=lambda l : l.tightId), # POG CutBasedIdTight
                          dict(name="highPtId", sel=lambda l : l.highPtId>0), # >0 = POG tracker high pT; 2 = global high pT, which includes the former
                          dict(name="isPFcand", sel=lambda l : l.isPFcand),
                          dict(name="isGlobal", sel=lambda l : l.isGlobal),  # Note: this is looser than nanoAOD presel.
                          dict(name="isTracker", sel=lambda l : l.isTracker),# Note: this is looser than nanoAOD presel.
                          dict(name="isTrackerArb", sel=lambda l : l.isTracker and l.nStations>0), # Arbitrated tracker muon. Note: this is looser than nanoAOD presel.
                          dict(name="inTimeMuon", sel=lambda l : l.inTimeMuon), # 
                          ]

            # Add variable to store the worst value of a given quantity among the 4 leptons of a candidate, for optimization studies.
            # Worst is intended as lowest value (as for an MVA), unless the variable's name starts with "max".
            self.muonIDVars=[dict(name="maxsip3d", sel=lambda l : l.sip3d if (abs(l.dxy)<0.5 and abs(l.dz) < 1) else 999.), # dxy, dz cuts included with SIP
                             dict(name="maxpfRelIso03FsrCorr", sel=lambda l : l.pfRelIso03FsrCorr), # FSR-corrected iso, DR=0.3
                             dict(name="maxpfRelIso03_all", sel=lambda l : l.pfRelIso03_all),
                             dict(name="maxpfRelIso04_all", sel=lambda l : l.pfRelIso04_all),
                             dict(name="maxminiPFRelIso_all", sel=lambda l : l.miniPFRelIso_all), # miniIso
                             dict(name="mvaLowPt", sel=lambda l : l.mvaLowPt if (l.looseId and l.sip3d<4. and l.dxy<0.5 and l.dz < 1) else -2.), # additional presel is required, cf: https://github.com/cms-sw/cmssw/blob/90f498af750cf4271c0a988fef352b0698012a40/PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc#L762-L764
#                             dict(name="promptMVA", sel=lambda l : l.promptMVA if (l.looseId and l.sip3d<4. and l.dxy<0.5 and l.dz < 1) else -2.), # former mvaTTH. adding H4l preselection for consistencty with mvaLowPt; this is looser than the original recommendation (https://twiki.cern.ch/twiki/bin/viewauth/CMS/LeptonMVA). Variable retrained in v14 (and renamed)
                             dict(name="mvaMuID", sel=lambda l : l.mvaMuID if (l.looseId and l.sip3d<4. and l.dxy<0.5 and l.dz < 1) else -2.), # muon MVA from 22-001. FIXME: Was retrained in v14; using H4l preselection for consistency, see above
                             ]


        # Data-MC SFs. 
        # NanoAODTools provides a module based on LeptonEfficiencyCorrector.cc, but that does not seem to be flexible enough for us:
        # https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/lepSFProducer.py
        #self.lepSFHelper = LeptonSFHelper(self.DATA_TAG) # FIXME for 2016 UL samples: requires passing bool preVFP

        # Example of adding control histograms (requires self.writeHistFile = True)
        # def beginJob(self,histFile=None, histDirName=None):
        #    Module.beginJob(self, histFile, histDirName+"_ZZFiller")
        #    self.histFile=None # Hack to prevent histFile to be closed before other modules write their histograms
        #    self.h_ZZMass = ROOT.TH1F('ZZMass','ZZMass',130,70,200)
        #    self.addObject(self.h_ZZMass)

    def endJob(self):
        with open('lheInfo_zzFailures_noTau_v2.json', 'w') as outfile:
            dump(self.failed_lhe_info, outfile, indent=4)
        with open('lheInfo_zzPasses_noTau.json', 'w') as outfile:
            dump(self.passed_lhe_info, outfile, indent=4)
        print("", flush=True)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.sign_reg  = ("OS", "SS")
        self.mass_reg  = ("HighMass", "MidMass", "LowMass")
        #self.sip_reg   = ("PassSIP", "FailSIP", "NoSIP")
        self.sip_reg   = ("PassSIP", "FailSIP")
        #self.sip_reg   = ("NoSIP")
        
        self.props = (("mass", "F"), ("pt", "F"), ("eta", "F"), ("phi", "F"), ("cosTheta1", "F"), ("cosTheta3", "F"), ("cosThetaStar", "F"), ("delPhiStar", "F"), ("delRapidity", "F"), ("massPreFSR", "F"),
                          ("Z1mass", "F"), ("Z1pt", "F"), ("Z1eta", "F"), ("Z1phi", "F"), ("Z1flav", "I"),
                          ("Z2mass", "F"), ("Z2pt", "F"), ("Z2eta", "F"), ("Z2phi", "F"), ("Z2flav", "I"),
                          ("KD", "F"), ("Z1l1Idx", "S"), ("Z1l2Idx", "S"), ("Z2l1Idx", "S"), ("Z2l2Idx", "S"))

        self.out.branch("nZZCand", "I",  title="ZZ candidates passing the full selection")
        self.out.branch("nZLLCand", "I", title="Z+LL control region candidates")
        
        self.out.branch("SR", "O", lenVar="nZZCand", title="Event candidate is in SR") # SR Filter Branch
        for prop, typ in self.props:
            sr_branch = "ZZCand_{}".format(prop)
            cr_branch = "ZLLCand_{}".format(prop)
            
            self.out.branch(sr_branch, typ, lenVar="nZZCand")
            self.out.branch(cr_branch, typ, lenVar="nZLLCand")
            
        if self.isMC:
            self.out.branch("ZZCand_dataMCWeight", "F", lenVar="nZZCand", title="data/MC efficiency correction weight")
            self.out.branch("ZLLCand_dataMCWeight", "F", lenVar="nZLLCand", title="data/MC efficiency correction weight")

        # CR Filter branches
        #self.cr_regions = ["{}_{}_{}".format(sign, sip, mass) for sign in self.sign_reg for sip in self.sip_reg for mass in self.mass_reg]
        self.cr_regions = ["{}_NoSIP_{}".format(sign, mass) for sign in self.sign_reg for mass in self.mass_reg]
        for branch in self.cr_regions:
            self.out.branch(branch, "O", lenVar="nZLLCand")
        
        self.filt_regions = ["SR"] + self.cr_regions

        self.failed_lhe_info = dict(
            el_pt    = [],
            el_eta   = [],
            mu_pt    = [],
            mu_eta   = []
        )
        self.passed_lhe_info = dict(
            el_pt    = [],
            el_eta   = [],
            mu_pt    = [],
            mu_eta   = []
        )

        # self.non_tau_fails = dict(
        #     LHEInfo = dict(
        #         ZZ = {"mass": []},
        #         Z1 = {
        #             "Leps": {
        #                 "pT": [],
        #                 "eta": []
        #             },
        #             "mass": [],
        #             "flavor": [],
        #         },
        #         Z2 = {
        #             "Leps": {
        #                 "pT": [],
        #                 "eta": []
        #             },
        #             "mass": [],
        #             "flavor": []
        #         }
        #     )
        # )
    def fill_dict(self, the_dict, lhepart):
        for i in range(4, 8):
            lep = lhepart[i]
            if abs(lep.pdgId) == 11:                    
                the_dict["el_pt"].append(lep.pt)
                the_dict["el_eta"].append(lep.eta)
            elif abs(lep.pdgId) == 13:
                the_dict["mu_pt"].append(lep.pt)
                the_dict["mu_eta"].append(lep.eta)

        return the_dict

    def buildZs(self, lep_idx_pairs, lep_pairs, fsrPhotons, looseMass=(12.,120.), tightMass=10.):
        """Builds Z candidates from pairs of leptons and their associated indices. Both leptons
        must pass the preselection, be SF and pass a deltaR cut to be considered. Lepton pairs
        which pass all of the above are used to instiate a ZCand object which constructs the Z
        from the leptons and their associated FSR photons. For the object to be saved to the output
        cands dictionary, it must pass the looseMass criteria and both of its leptons must pass
        FullSelNoSIP. Z1 candidates are defined as those which are OS+PassSIP+OnShell."""

        cands = dict(
            Z1_Cands   = [],
            #OS_PassSIP = [], # Only difference from above is that the Z2 is off shell
            #OS_FailSIP = [],
            OS_NoSIP   = [],
            #SS_PassSIP = [],
            #SS_FailSIP = [],
            SS_NoSIP   = [],
        )

        for lep_idx, lep_pair in zip(lep_idx_pairs, lep_pairs):
            l1, l2 = lep_pair
            if self.leptonPresel(l1) and self.leptonPresel(l2):

                if abs(l1.pdgId) == abs(l2.pdgId) and self.passDeltaR(l1, l2):
                    # Set a default order for OS leptons in the Z candidate: 1=l+, 2=l-
                    idx1, idx2 = lep_idx
                    if (l1.pdgId*l2.pdgId<0 and l1.pdgId>0):
                        idx1, idx2 = idx2, idx1
                        l1, l2     = l2, l1
                    
                    aZ = self.ZCand(idx1, idx2, l1, l2, fsrPhotons, tightMass=tightMass)

                    if aZ.PassSelNoSIP and (aZ.M>looseMass[0] and aZ.M<looseMass[1]):
                        if aZ.Z1Cand:
                            cands["Z1_Cands"].append(aZ)
                        else:
                            key = "OS" if aZ.OS else "SS"
                            
                            cands[key+"_NoSIP"].append(aZ)
                            '''if aZ.PassSIP:
                                cands[key+"_PassSIP"].append(aZ)
                            else:
                                cands[key+"_FailSIP"].append(aZ)'''
                            
                            '''key += "_PassSIP" if aZ.PassSIP else "_FailSIP"
                            cands[key].append(aZ)'''
        
        return cands

    def getSortedZPairs(self, z_cands, reg, z1=None):
        """Accepts a dictionary Z_Cands which contains zed candidates for each
        region, the region in question, and a Z1 candidate if sorting a CR. Returns
        a list of tuples (Z1, Z2) sorted in two ways with Z1 defined as:
        1) closest to Z1 mass if both candidates are possible Z1s (SR)
        2) the passed value of z1. This should only be used in two cases:
            a) there is only one Z1 candidate in the event.
            b) all combinations of (Z1_a, Z1_b) failed to pass the SR criteria.
        
        The final returned list is sorted by closest Z1 (Z2) mass for the SR (CRs)."""
        
        cands = z_cands[reg]
        which_z_idx = 0 if reg == "Z1_Cands" else 1
        if reg == "Z1_Cands":
            z_pairs = combinations(cands, 2)
            # Define Z1 as closest mass to mZ, ex [(88, 90), (84, 91)] --> [(90, 88), (91, 84)]
            z_pairs = [(zp[0], zp[1]) if abs(zp[0].M - self.ZmassValue) < abs(zp[1].M - self.ZmassValue) else (zp[1], zp[0]) for zp in z_pairs]
        else:
            z_pairs = [(z1, ll) for ll in cands]
        
        # Sort total list by Z1 mass if both Zs are Z1Cands, Z2 mass if only Z1 is a Z1Cand
        return sorted(z_pairs, key = lambda x: abs(x[which_z_idx].M - self.ZmassValue))

    def bestInSR(self, Z_Cands, sr_cands):
        """Accepts a dictionary Z_Cands which contains zed candidates for each
        region, and an empty list sr_cands. Loops through a sorted list of
        (Za, Zb) pairs (both of which must pass Z1 criteria), appending the
        first one that passes SR criteria."""
        
        sorted_pairs = self.getSortedZPairs(Z_Cands, "Z1_Cands")
        for zz_pair in sorted_pairs:
            zz = self.ZZCand_Base(*zz_pair)
            if zz.SR:
                sr_cands.append(zz)
                break
        
        return sr_cands        

    def bestInCR(self, Z_Cands, z1_cands, temp_zlls):
        """Accepts a dictionary Z_Cands which contains zed candidates for each
        region, a list of possible z1_cands (candidates that pass all
        selections), and an empty dictionary.For each region and each z1 candidate,
        constructs an ordered list of best Z_LL pairs, and saves the first pair that
        passes the baseline selection to the dictionary."""

        ll_regions = list(Z_Cands.keys())
        ll_regions.remove("Z1_Cands")
        for region in ll_regions:
            
            if len(Z_Cands[region]) < 1: continue
            
            for z1 in z1_cands:
                sorted_pairs = self.getSortedZPairs(Z_Cands, region, z1)
                for zz_pair in sorted_pairs:
                    zz = self.ZZCand_Base(*zz_pair)
                    if zz.passBaseline:
                        if zz.HighMass:
                            if zz.Z2OnShell : temp_zlls[region+"_HighMass"] = zz
                            else: continue
                        elif zz.MidMass: temp_zlls[region+"_MidMass"] = zz
                        elif zz.LowMass: temp_zlls[region+"_LowMass"] = zz
                        break
                        
        return temp_zlls

    def bestZLLs(self, temp_zlls, final_cands):
        """Takes a dictionary of temporary ZLL candidates (temp_zlls)
        containing the best candidate for each CR. CRs for which no candidate
        passes the baseline selection should not be in temp_zlls.

        For each pair of regions, check if candidate leptons overlap:
        --> Overlap:    choose best cand by Z2 mass.
        --> No overlap: keep both candidates.
        """

        # Candidate leptons for each region
        lep_lists = {reg: temp_zlls[reg].leps() for reg in temp_zlls}
        # Pairs of regions
        reg_pairs = combinations(lep_lists.keys(), 2)
        
        for reg_pair in reg_pairs:
            reg_1, reg_2 = reg_pair
            
            # NoSIP region _does_ overlap with Pass/Fail Regions
            if ("NoSIP" in reg_1) or ("NoSIP" in reg_2):
                # Only overlap allowed is in SIP Flags
                if np.array_equal(np.char.equal(reg_1.split("_"), reg_2.split("_")), [True, False, True]):
                    final_cands[reg_1] = temp_zlls[reg_1]
                    final_cands[reg_2] = temp_zlls[reg_2]
                    continue

            cand_1_leps, cand_2_leps = lep_lists[reg_1], lep_lists[reg_2]
            
            leps_overlap = bool(set(cand_1_leps) & set(cand_2_leps))

            if leps_overlap:
                cand_1_z2, cand_2_z2 = temp_zlls[reg_1].Z2, temp_zlls[reg_2].Z2
                
                if abs(cand_1_z2.M - self.ZmassValue) < abs(cand_2_z2.M - self.ZmassValue):
                    final_cands[reg_1] = temp_zlls[reg_1]
                else:
                    final_cands[reg_2] = temp_zlls[reg_2]
            
            else:
                final_cands[reg_1] = temp_zlls[reg_1]
                final_cands[reg_2] = temp_zlls[reg_2]
        
        return final_cands

    def write_branches(self, final_cands):
        
        cand_props = candProps(final_cands, self.filt_regions)
        
        reg_idx = [key for key in final_cands][0]

        if self.isMC:
            mcWeight = [self.getDataMCWeight(final_cands[reg_idx].leps())]
            if cand_props.branches["SR"][0]:
                self.out.fillBranch("ZZCand_dataMCWeight", mcWeight)
            else:
                self.out.fillBranch("ZLLCand_dataMCWeight", mcWeight)
        
        prep = "ZZCand_" if cand_props.branches["SR"][0] else "ZLLCand_"

        for branch, vals in cand_props.branches.items():
            if not (branch == "SR" or "OS" in branch or "SS" in branch):
                branch = prep + branch
            self.out.fillBranch(branch, vals)
        


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if self.DEBUG : print ('Event {}:{}:{}'.format(event.run,event.luminosityBlock,event.event))

        # Collections
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        leps = list(electrons)+list(muons)
        nlep=len(leps)

        # For Debugging
        genpart = Collection(event, "GenPart")
        lhepart = Collection(event, "LHEPart")
        #lhePart = lheParts(lhepart)

        ### Skip events with too few leptons (min 3 if ZL is included, 4 otherwise) 
        if nlep < 3 or (not self.addZLCR and nlep < 4) :
            return False

        lep_idxs      = [*range(0,nlep,1)]

        lep_idx_pairs = combinations(lep_idxs, 2)
        lep_pairs     = combinations(leps, 2)

        # {Z1_cands: [...], OS_PassSIP = [...], OS_FailSIP = [...], SS_PassSIP = [...], SS_FailSIP = [...]}
        Z_Cands = self.buildZs(lep_idx_pairs, lep_pairs, fsrPhotons)
        z1_cands = Z_Cands["Z1_Cands"]
        
        temp_zlls = {} # Save best ZLL candidate in each region separately
        sr_cands = []
        

        # More than one Z1_Cands --> Possible SR
        if len(z1_cands) > 1:

            sr_cands = self.bestInSR(Z_Cands, sr_cands)

            # >= 2 Z1 Cands but no ZZ Cand passing SR criteria --> check for CRs
            if len(sr_cands) == 0:
                temp_zlls = self.bestInCR(Z_Cands, z1_cands, temp_zlls)


        elif len(z1_cands) == 1:
            
            if sum([len(Z_Cands[reg]) for reg in Z_Cands.keys()]) < 2:
                self.fill_dict(self.failed_lhe_info, lhepart)
                #print("\nFail: <2 Tot Cands")
                # lhePart._setInfo()
                # if not lhePart.has_tau:
                #     self.non_tau_fails["LHEInfo"]["ZZ"]["mass"].append(lhePart.zz_mass)
                #     self.non_tau_fails["LHEInfo"]["Z1"]["mass"].append(lhePart.z1_mass)
                #     self.non_tau_fails["LHEInfo"]["Z1"]["flavor"].append(lhePart.z1_flav)
                #     self.non_tau_fails["LHEInfo"]["Z2"]["mass"].append(lhePart.z2_mass)
                #     self.non_tau_fails["LHEInfo"]["Z2"]["flavor"].append(lhePart.z2_flav)
                #     self.non_tau_fails["LHEInfo"]["Z1"]["Leps"]["pT"].append((lhePart.lep_1_pt, lhePart.lep_2_pt))
                #     self.non_tau_fails["LHEInfo"]["Z1"]["Leps"]["eta"].append((lhePart.lep_1_eta, lhePart.lep_2_eta))
                #     self.non_tau_fails["LHEInfo"]["Z2"]["Leps"]["pT"].append((lhePart.lep_3_pt, lhePart.lep_4_pt))
                #     self.non_tau_fails["LHEInfo"]["Z2"]["Leps"]["eta"].append((lhePart.lep_3_eta, lhePart.lep_4_eta))
                return False
            else:
                temp_zlls = self.bestInCR(Z_Cands, z1_cands, temp_zlls)
        
        else:
            #print("\nFail: No Z Cands")
            # test = [(gp.pdgId, getParentID(gp, genpart)) for gp in genpart]
            # test_2 = [(gp.pdgId, Mother(gp, genpart)) for gp in genpart]
            # lhePart._setInfo()
            # if not lhePart.has_tau:
            #     self.non_tau_fails["LHEInfo"]["ZZ"]["mass"].append(lhePart.zz_mass)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["mass"].append(lhePart.z1_mass)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["flavor"].append(lhePart.z1_flav)
            #     self.non_tau_fails["LHEInfo"]["Z2"]["mass"].append(lhePart.z2_mass)
            #     self.non_tau_fails["LHEInfo"]["Z2"]["flavor"].append(lhePart.z2_flav)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["Leps"]["pT"].append((lhePart.lep_1_pt, lhePart.lep_2_pt))
            #     self.non_tau_fails["LHEInfo"]["Z1"]["Leps"]["eta"].append((lhePart.lep_1_eta, lhePart.lep_2_eta))
            #     self.non_tau_fails["LHEInfo"]["Z2"]["Leps"]["pT"].append((lhePart.lep_3_pt, lhePart.lep_4_pt))
            #     self.non_tau_fails["LHEInfo"]["Z2"]["Leps"]["eta"].append((lhePart.lep_3_eta, lhePart.lep_4_eta))
            self.fill_dict(self.failed_lhe_info, lhepart)
            return False
        


        final_cands = {}
        if (len(sr_cands) == 0) and len(temp_zlls) == 0:
            #print("\nFail: len(sr)==0 and len(temp_zlls)==0")
            #test = [(gp.pdgId, getParentID(gp, genpart)) for gp in genpart]
            #test_2 = [(gp.pdgId, Mother(gp, genpart)) for gp in genpart]
            # lhePart._setInfo()
            # if not lhePart.has_tau:
            #     self.non_tau_fails["LHEInfo"]["ZZ"]["mass"].append(lhePart.zz_mass)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["mass"].append(lhePart.z1_mass)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["flavor"].append(lhePart.z1_flav)
            #     self.non_tau_fails["LHEInfo"]["Z2"]["mass"].append(lhePart.z2_mass)
            #     self.non_tau_fails["LHEInfo"]["Z2"]["flavor"].append(lhePart.z2_flav)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["Leps"]["pT"].append((lhePart.lep_1_pt, lhePart.lep_2_pt))
            #     self.non_tau_fails["LHEInfo"]["Z1"]["Leps"]["eta"].append((lhePart.lep_1_eta, lhePart.lep_2_eta))
            #     self.non_tau_fails["LHEInfo"]["Z2"]["Leps"]["pT"].append((lhePart.lep_3_pt, lhePart.lep_4_pt))
            #     self.non_tau_fails["LHEInfo"]["Z2"]["Leps"]["eta"].append((lhePart.lep_3_eta, lhePart.lep_4_eta))
            self.fill_dict(self.failed_lhe_info, lhepart)
            return False
        
        else:
            # sr_cands can have at most one candidate inside
            if len(sr_cands) > 0:
                final_cands["SR"] = sr_cands[0]
            
            elif len(temp_zlls) == 1:
                region = list(temp_zlls.keys())[0]
                final_cands[region] = temp_zlls[region]
            else:
                # fill final_cands dict with best of non-overlapping cands
                final_cands = self.bestZLLs(temp_zlls, final_cands)

        if len(final_cands) == 0:
            #print("\nFail: len(final_cands)==0")
            #test = [(gp.pdgId, getParentID(gp, genpart)) for gp in genpart]
            #test_2 = [(gp.pdgId, Mother(gp, genpart)) for gp in genpart]
            # lhePart._setInfo()
            # if not lhePart.has_tau:
            #     self.non_tau_fails["LHEInfo"]["ZZ"]["mass"].append(lhePart.zz_mass)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["mass"].append(lhePart.z1_mass)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["flavor"].append(lhePart.z1_flav)
            #     self.non_tau_fails["LHEInfo"]["Z2"]["mass"].append(lhePart.z2_mass)
            #     self.non_tau_fails["LHEInfo"]["Z2"]["flavor"].append(lhePart.z2_flav)
            #     self.non_tau_fails["LHEInfo"]["Z1"]["Leps"]["pT"].append((lhePart.lep_1_pt, lhePart.lep_2_pt))
            #     self.non_tau_fails["LHEInfo"]["Z1"]["Leps"]["eta"].append((lhePart.lep_1_eta, lhePart.lep_2_eta))
            #     self.non_tau_fails["LHEInfo"]["Z2"]["Leps"]["pT"].append((lhePart.lep_3_pt, lhePart.lep_4_pt))
            #     self.non_tau_fails["LHEInfo"]["Z2"]["Leps"]["eta"].append((lhePart.lep_3_eta, lhePart.lep_4_eta))
            self.fill_dict(self.failed_lhe_info, lhepart)
            return False

        self.write_branches(final_cands)
        self.fill_dict(self.passed_lhe_info, lhepart)
            
        return True
    
    class ZCand: 
        def __init__(self, l1Idx, l2Idx, l1, l2, fsrPhotons, tightMass=10):
            self.l1Idx = l1Idx
            self.l2Idx = l2Idx
            self.l1    = l1
            self.l2    = l2
            self.fsr1Idx = self.l1.fsrPhotonIdx
            self.fsr2Idx = self.l2.fsrPhotonIdx
    
            self.l1DressedP4 = self.l1.p4()
            self.l2DressedP4 = self.l2.p4()
            if self.fsr1Idx>=0 : self.l1DressedP4 += fsrPhotons[self.fsr1Idx].p4()
            if self.fsr2Idx>=0 : self.l2DressedP4 += fsrPhotons[self.fsr2Idx].p4()
    
            self.p4 = self.l1DressedP4 + self.l2DressedP4
    
            self.M = self.p4.M() # cache the mass as it is used often

            self.ZmassValue = 91.1876

            self.OS = self.l1.pdgId == -self.l2.pdgId
            self.SS = self.l1.pdgId == self.l2.pdgId

            self.PassSelNoSIP = False
            self.PassSIP      = False
            if self.l1.ZZFullSelNoSIP and self.l2.ZZFullSelNoSIP:
                self.PassSelNoSIP = True
            if self.l1.ZZFullSel and self.l2.ZZFullSel:
                self.PassSIP = True

            self.OnShell = abs(self.M - self.ZmassValue) < tightMass

            self.Z1Cand  = self.OS and self.PassSIP and self.OnShell
    
        def sumpt(self) : # sum of lepton pTs, used to sort candidates
            return self.l1.pt+self.l2.pt
    
        def finalState(self) :
            return self.l1.pdgId*self.l2.pdgId
    
    # Temporary class to store information on a ZZ candidate.
    class ZZCand_Base:
        def __init__(self, Z1, Z2):
            self.Z1 = Z1
            self.Z2 = Z2

            self.lepsExclusive = self.Z1.l1 != self.Z2.l1 and self.Z1.l1 != self.Z2.l2 and self.Z1.l2 != self.Z2.l1 and self.Z1.l2 != self.Z2.l2

            self.p4 = Z1.p4+Z2.p4

            self.M = self.p4.M()

            angVars = AngularVars(self)
            self.cosTheta1 = angVars.cosTheta("1")
            self.cosTheta3 = angVars.cosTheta("2")
            self.cosThetaStar = angVars.cosThetaStar()
            self.delPhiStar = angVars.delPhiStar()
            self.delRapidity = angVars.delRapidity()

            self.HighMass = self.M > 180.
            self.MidMass  = 140. < self.M < 180.
            self.LowMass  = 105. < self.M < 140.

            # Minimum requirements to belong to any region
            self.passBaseline = self.Z1.Z1Cand and self.lepsExclusive and self.passLepPts() and self.passQCDandDeltaR() and self.M > 105
            
            self.SR = self.passBaseline and self.Z2.Z1Cand and self.HighMass

            self.OS = self.Z2.OS
            self.SS = self.Z2.SS

            self.PassSIP = self.Z2.PassSIP
            self.Z2OnShell = self.Z2.OnShell

        def finalState(self) :
            return self.Z1.finalState()*self.Z2.finalState()

        def massPreFSR(self) :
            return (self.Z1.l1.p4()+self.Z1.l2.p4()+self.Z2.l1.p4()+self.Z2.l2.p4()).M()

        def leps(self) :
            return([self.Z1.l1, self.Z1.l2, self.Z2.l1, self.Z2.l2])

        def passLepPts(self, lead=20., sublead=10.):
            # Original / HZZ Requirements
            # lep_pts = [l.pt for l in self.leps()]
            # lep_pts.sort()

            # return lep_pts[3] > lead and lep_pts[2] > sublead

            # To kill Z+X background
            lep_reqs = [20., 15., 15., 15.]
            return sum([l_pt > req for lep, req in zip(lep_pts, lep_reqs)]) == 4

        def passQCDandDeltaR(self):
            zzleps = self.leps()
            # QCD suppression on all OS pairs, regardelss of flavour
            passQCD = True    # QCD suppression on all OS pairs, regardelss of flavour
            passDeltaR = True # DeltaR>0.02 cut among all leptons to protect against split tracks
            for k in range(4):
                for l in range (k+1,4):
                    if zzleps[k].charge!=zzleps[l].charge and (zzleps[k].p4()+zzleps[l].p4()).M()<=4.:
                        passQCD = False
                        break
                    if deltaR(zzleps[k].eta, zzleps[k].phi, zzleps[l].eta, zzleps[l].phi) <= 0.02 :
                        passDeltaR = False
                        break

            return passQCD and passDeltaR

    ### Compute lepton efficiency scale factor
    def getDataMCWeight(self, leps) :
       if self.year > 2023:
           return 1.
       dataMCWeight = 1.
       for lep in leps:           
           myLepID = abs(lep.pdgId)
           mySCeta = lep.eta
           isCrack = False
           if myLepID==11 :
               mySCeta = lep.eta + lep.deltaEtaSC
               # FIXME: isGap() is not available in nanoAODs, and cannot be recomputed easily based on eta, phi. We thus use the non-gap SFs for all electrons.

           # Deal with very rare cases when SCeta is out of 2.5 bounds
           mySCeta = min(mySCeta,2.49)
           mySCeta = max(mySCeta,-2.49)

           SF = self.lepSFHelper.getSF(self.year, myLepID, lep.pt, lep.eta, mySCeta, isCrack)
#           SF_Unc = self.lepSFHelper.getSFError(year, myLepID, lep.pt, lep.eta, mySCeta, isCrack)
           dataMCWeight *= SF

       return dataMCWeight