import os

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

import ROOT
from ROOT.Math import LorentzVector, PtEtaPhiM4D, Boost, LorentzRotation

import numpy as np

xsecs = dict(
    ZLZL = 0.0008994,
    ZLZT = 0.001773,
    ZTZL = 0.001758,
    ZTZT = 0.01025
)

def getIdx(lhe_parts, pdgId):
    return [lhe.pdgId for lhe in lhe_parts].index(pdgId)

def getLorentzVec(lhe_parts, pdgId, idx=-1):
    if idx == -1:
        idx = getIdx(lhe_parts, pdgId)
    return LorentzVector(PtEtaPhiM4D('double'))(
        lhe_parts[idx].pt,
        lhe_parts[idx].eta,
        lhe_parts[idx].phi,
        lhe_parts[idx].mass
    )

def getBoost(p4):
    return Boost(p4.BoostToCM())

def getIdx(lhe_parts, pdgId):
    return [lhe.pdgId for lhe in lhe_parts].index(pdgId)

def getZ(lhe_parts, pdgId):
    z1 = getLorentzVec(lhe_parts, 23, 2)
    z2 = getLorentzVec(lhe_parts, 23, 3)

    ll = getLorentzVec(lhe_parts, -1*pdgId) + getLorentzVec(lhe_parts, pdgId)

    z1_mass, z2_mass, ll_mass = z1.M(), z2.M(), ll.M()

    if np.abs(z1_mass - ll_mass) < np.abs(z2_mass - ll_mass):
        return z1
    else:
        return z2

def cosTheta(lhe_parts, z):
    z_ee_lab   = getZ(lhe_parts, 11)
    z_mumu_lab = getZ(lhe_parts, 13)

    lab_4l = z_ee_lab + z_mumu_lab

    lab_pos_lep = getLorentzVec(lhe_parts, -11) if z==1 else getLorentzVec(lhe_parts, -13)

    boost_to_4l = getBoost(lab_4l)
    boost_to_ZCM = getBoost(z_ee_lab) if z==1 else getBoost(z_mumu_lab)

    z_4l_cm = boost_to_4l*z_ee_lab if z==1 else boost_to_4l*z_mumu_lab
    pos_lep_z_cm = boost_to_ZCM*lab_pos_lep

    return ROOT.Math.VectorUtil.CosTheta(pos_lep_z_cm, z_4l_cm)

def getCosThetaStar(lhe_parts):
    z1_lab = getLorentzVec(lhe_parts, 23, 2)
    z2_lab = getLorentzVec(lhe_parts, 23, 3)

    lab_4l = z1_lab + z2_lab

    boost_to_4l = getBoost(lab_4l)

    z1_4l_cm = boost_to_4l*z1_lab

    return ROOT.Math.VectorUtil.CosTheta(z1_4l_cm, lab_4l)

def delRapidity(lhe_parts):
    z1_lab = getLorentzVec(lhe_parts, 23, 2)
    z2_lab = getLorentzVec(lhe_parts, 23, 3)
    
    z1_rapidity = 0.5*np.log((z1_lab.E() + z1_lab.Pz())/(z1_lab.E() - z1_lab.Pz()))
    z2_rapidity = 0.5*np.log((z2_lab.E() + z2_lab.Pz())/(z2_lab.E() - z2_lab.Pz()))

    return np.abs(z1_rapidity - z2_rapidity)

def delPhiStar(lhe_parts):
    pos_el_p4 = getLorentzVec(lhe_parts, -11)
    pos_mu_p4 = getLorentzVec(lhe_parts, -13)

    z_ee_lab   = getZ(lhe_parts, 11)
    z_mumu_lab = getZ(lhe_parts, 13)

    z_ee_boost   = getBoost(z_ee_lab)
    z_mumu_boost = getBoost(z_mumu_lab)

    pos_el_zCM = z_ee_boost*pos_el_p4
    pos_mu_zCM = z_mumu_boost*pos_mu_p4

    # Above is good, need Phi w.r.t. its parent Z spatial direction in the CM frame

    lab_4l = z_ee_lab + z_mumu_lab
    
    boost_to_4l = getBoost(lab_4l)

    z_ee_4l   = boost_to_4l*z_ee_lab
    z_mumu_4l = boost_to_4l*z_mumu_lab

    phiStar_e  = ROOT.Math.VectorUtil.Angle(pos_el_zCM, z_ee_4l)
    phiStar_mu = ROOT.Math.VectorUtil.Angle(pos_mu_zCM, z_mumu_4l)

    # Test with z_ll_cm (its own rest frame)

    # z_ee_cm   = z_ee_boost*z_ee_lab
    # z_mumu_cm = z_mumu_boost*z_mumu_lab

    # phiStar_e = ROOT.Math.VectorUtil.Angle(pos_el_zCM, z_ee_cm)
    # phiStar_mu = ROOT.Math.VectorUtil.Angle(pos_mu_zCM, z_mumu_cm)

    # Test with z_ll_lab

    # phiStar_e = ROOT.Math.VectorUtil.Angle(pos_el_zCM, z_ee_lab)
    # phiStar_mu = ROOT.Math.VectorUtil.Angle(pos_mu_zCM, z_mumu_lab)

    diff = np.abs(phiStar_e - phiStar_mu)
    return np.rad2deg(min(diff, 2*np.pi - diff))

    # az_diff = np.abs(pos_el_zCM.Phi() - pos_mu_zCM.Phi())

    # return min(az_diff, 2*np.pi - az_diff)


def make_vars(filename, state):

    cosTheta1Hist = ROOT.TH1F("cosTheta1_"+state, "cosTheta1_"+state, 20, -1., 1.)
    cosTheta1Hist.GetXaxis().SetTitle("cos(#theta_1)")
    cosTheta1Hist.GetYaxis().SetTitle("Events")

    cosTheta3Hist = ROOT.TH1F("cosTheta3_"+state, "cosTheta3_"+state, 20, -1., 1.)
    cosTheta3Hist.GetXaxis().SetTitle("cos(#theta_3)")
    cosTheta3Hist.GetYaxis().SetTitle("Events")

    cosThetaStarHist = ROOT.TH1F("cosThetaStar_"+state, "cosThetaStar_"+state, 20, -1., 1.)
    cosThetaStarHist.GetXaxis().SetTitle("cos(#theta^*)")
    cosThetaStarHist.GetYaxis().SetTitle("Events")

    delRapHist = ROOT.TH1F("delRapidity_"+state, "delRapidity_"+state, 40, 0., 4.0)
    delRapHist.GetXaxis().SetTitle("Delta(y_(z z')")
    delRapHist.GetYaxis().SetTitle("Events")

    delPhiHist = ROOT.TH1F("delPhiStar_"+state, "delPhiStar_"+state, 36, 0., 180.)
    delPhiHist.GetXaxis().SetTitle("Delta(Phi_(e+ mu+)")
    delPhiHist.GetYaxis().SetTitle("Events")

    f = ROOT.TFile.Open(filename)

    event = f.Events

    nEntries = event.GetEntries()

    iEntry=0
    while iEntry<nEntries and event.GetEntry(iEntry):
        iEntry+=1
        
        LHEParts = Collection(event, "LHEPart")

        cosTheta1 = cosTheta(LHEParts, 1)
        cosTheta3 = cosTheta(LHEParts, 2)
        cosThetaStar = getCosThetaStar(LHEParts)
        delRap = delRapidity(LHEParts)
        delPhi = delPhiStar(LHEParts)

        cosTheta1Hist.Fill(cosTheta1)
        cosTheta3Hist.Fill(cosTheta3)
        cosThetaStarHist.Fill(cosThetaStar)
        delRapHist.Fill(delRap)
        delPhiHist.Fill(delPhi)

    f.Close()

    xsec = xsecs[state]

    cosTheta1Hist.Scale(xsec)
    cosTheta3Hist.Scale(xsec)
    cosThetaStarHist.Scale(xsec)
    delRapHist.Scale(xsec)
    delPhiHist.Scale(xsec)

    # cosTheta1Hist.Scale(1.)
    # cosTheta3Hist.Scale(1.)
    # cosThetaStarHist.Scale(1.)
    # delRapHist.Scale(1.)
    # delPhiHist.Scale(1.)
    
    return [cosTheta1Hist, cosTheta3Hist, cosThetaStarHist, delRapHist, delPhiHist]     

if __name__=="__main__":

    base = "/afs/cern.ch/user/i/iehle/cmssw/CMSSW_13_0_14/src/LHEprod/LHEDumper/"

    files = [os.path.join(base, "qq{}_2e2mu_LHE.root".format(pol)) for pol in ["ZLZL", "ZLZT", "ZTZT"]]

    variables = ["LHEPart_{}".format(var) for var in ["pdgID", "status", "pt", "spin", "eta", "phi", "mass"]]

    outfile = "qqZZ_2e2mu_LHE_vars_byXSec.root"

    OutFile = ROOT.TFile.Open(outfile, "recreate")

    for state in ["ZLZL", "ZLZT", "ZTZL", "ZTZT"]:
        file = os.path.join(base, "qq{}_2e2mu_LHE.root".format(state))
        hists = make_vars(file, state)
        for hist in hists:
            OutFile.WriteObject(hist, hist.GetName())

    OutFile.Close()
