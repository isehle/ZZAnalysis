import ROOT
from ROOT.Math import LorentzVector,PtEtaPhiM4D

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR

from itertools import combinations

from tqdm import tqdm

inpath = "root://eos.grif.fr//eos/grif/cms/llr/store/user/iehle/qqZLZLTo4l_noTau_5f_Summer22EraEFG_onlyCERN/crab_qqzlzl-4l-noTau-cern/250120_103551/0000/qqZLZL_noTau_hadd.root"

class lheParts:
    def __init__(self, lhepart):
        self.ZmassValue = 91.1876

        self.lhepart = lhepart

        self.leps = [lhepart[i] for i in range(4,8)]
        self.lep_pairs = combinations(self.leps, 2)

        self.z1_flav = self.leps[0].pdgId*self.leps[1].pdgId
        self.z2_flav = self.leps[2].pdgId*self.leps[3].pdgId
        self.zz_flav = self.z1_flav*self.z2_flav

        self._setP4s()

        self.pass_eta = lambda l: abs(l.eta) < 2.4 if abs(l.pdgId) == 13 else abs(l.eta) < 2.5
        self.presel   = True
        
        for lep in self.leps:
            self.presel *= self.pass_eta(lep)

        self.passDeltaR = lambda l1, l2: deltaR(l1.eta, l1.phi, l2.eta, l2.phi) > 0.02

        self.os_pair    = lambda l1, l2: l1.pdgId/abs(l1.pdgId) == -1*l2.pdgId/abs(l2.pdgId)

    def _getLorentzVec(self, idx):
        if type(idx) == int:
            return LorentzVector(PtEtaPhiM4D('double'))(
            self.lhepart[idx].pt,
            self.lhepart[idx].eta,
            self.lhepart[idx].phi,
            self.lhepart[idx].mass
        )
        else:
            return LorentzVector(PtEtaPhiM4D('double'))(
                idx.pt,
                idx.eta,
                idx.phi,
                idx.mass
            )

    def _setP4s(self):
        self.z1_p4 = self._getLorentzVec(2)
        self.z2_p4 = self._getLorentzVec(3)

        self.zz_p4 = self.z1_p4+self.z2_p4

        self.lep_1_p4 = self._getLorentzVec(4)
        self.lep_2_p4 = self._getLorentzVec(5)
        self.lep_3_p4 = self._getLorentzVec(6)
        self.lep_4_p4 = self._getLorentzVec(7)

        self.z1_ll = self.lep_1_p4 + self.lep_2_p4
        self.z2_ll = self.lep_3_p4 + self.lep_4_p4

        self.zz_4l = self.z1_ll + self.z2_ll

    def _checkLL(self):
        self.all_pass_dR = True
        self.n_pass_mZ = 0

        for lep_pair in self.lep_pairs:
            lep_1, lep_2 = lep_pair

            if lep_1.pdgId != -1*lep_2.pdgId:
                continue

            self.all_pass_dR *= self.passDeltaR(lep_1, lep_2)

            ll_p4 = self._getLorentzVec(lep_1) + self._getLorentzVec(lep_2)

            if abs(ll_p4.M() - self.ZmassValue) < 10:
                self.n_pass_mZ += 1

    def _check4l(self):
        el_pt_reqs = [20, 10, 7, 7]
        mu_pt_reqs = [20, 10, 5, 5]

        self.pass_lep_pt = True

        for i in range(4):
            lep = self.leps[i]
            if abs(lep.pdgId) == 11:
                self.pass_lep_pt *= lep.pt >= el_pt_reqs[i]
            if abs(lep.pdgId) == 13:
                self.pass_lep_pt *= lep.pt >= mu_pt_reqs[i]

        self.pass_os_mass = True

        for lep_pair in self.lep_pairs:
            lep_1, lep_2 = lep_pair
            
            self.all_pass_dR *= self.passDeltaR(lep_1, lep_2)
            
            if self.os_pair(lep_1, lep_2):
                ll_p4 = self._getLorentzVec(lep_1) + self._getLorentzVec(lep_2)

                self.pass_os_mass *= ll_4.M() >= 4

        self.pass_m4l = self.zz_4l.M() >= 180

        

infile = ROOT.TFile.Open(inpath)

event = infile.Events
nEntries = event.GetEntries()

# Preselection (Detetector acceptance (eta))
fail_presel   = 0

# Building Z Cands
fail_first_dR = 0
fail_mll_mz   = 0

# Building ZZCands
fail_lep_pts  = 0
fail_os_mass  = 0
fail_scnd_dR  = 0
fail_m4l      = 0

iEntry        = 0
pbar = tqdm(total = nEntries+1)
while iEntry < nEntries and event.GetEntry(iEntry):
    iEntry += 1
    pbar.update(1)
    lhepart = Collection(event, "LHEPart")
    
    lhePart = lheParts(lhepart)

    # Check Detector acceptance (eta)
    if not lhePart.presel:
        fail_presel += 1
        continue

    # Check building Z Cands
    lhePart._checkLL()
    if not lhePart.all_pass_dR:
        fail_first_dR += 1
        continue
    if lhePart.n_pass_mZ != 2:
        fail_mll_mz += 1
        continue

    # Check building ZZ Cands
    lhePart._check4l()
    if not lhePart.pass_lep_pt:
        fail_lep_pts += 1
        continue
    if not lhePart.all_pass_dR:
        fail_scnd_dR += 1
        continue
    if not lhePart.pass_os_mass:
        fail_os_mass += 1
        continue
    if not lhePart.pass_m4l:
        fail_m4l += 1
        continue
pbar.close()

infile.Close()

print("Total Failures: ", fail_presel+fail_first_dR+fail_mll_mz+fail_lep_pts+fail_os_mass+fail_scnd_dR+fail_m4l)
print("\nDetector Acceptance: \n")
print("Presel:   ", fail_presel/nEntries)
print("\nZ Cand Building: \n")
print("First dR: ", fail_first_dR/nEntries)
print("mLL:      ", fail_mll_mz/nEntries)
print("\nZZ Cand Building: \n")
print("Lep pTs:  ", fail_lep_pts/nEntries)
print("OS Mass:  ", fail_os_mass/nEntries)
print("Scnd dR:  ", fail_scnd_dR/nEntries)
print("m4l:      ", fail_m4l/nEntries)