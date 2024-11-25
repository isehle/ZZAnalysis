import ROOT

ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
using ROOT::RVecB;
RVecF applyRegFilt(RVecF &prop, RVecB &regMask){
    RVecF reg_prop;
    for (int i=0; i<prop.size(); i++){
        if (regMask.at(i) == true){
            reg_prop.push_back(prop.at(i));
        }
    }
    return reg_prop;
}
""")

ROOT.gInterpreter.Declare("""
using ROOT::RVecI;
using ROOT::RVecB;
RVecI applyRegFilt(RVecI &prop, RVecB &regMask){
    RVecI reg_prop;
    for (int i=0; i<prop.size(); i++){
        if (regMask.at(i) == true){
            reg_prop.push_back(prop.at(i));
        }
    }
    return reg_prop;
}
""")
                          
ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
using ROOT::RVecB;
float applyRegFilt(RVecF &pdfWeights, RVecB &regMask, int pdfBin){
    float w_pdf;
    for (int i=0; i<regMask.size(); i++){
        if (regMask.at(i) == true){
            w_pdf = pdfWeights.at(pdfBin);
            break;
        }
    }
    return w_pdf;
}
""")

ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
using ROOT::RVecB;
RVecF lepFromCand(RVecF &elProp, RVecF &muProp, ROOT::VecOps::RVec<short> &l1Idx, ROOT::VecOps::RVec<short> &l2Idx, RVecB &regMask){
    RVecF lepProp = ROOT::VecOps::Concatenate(elProp, muProp);

    RVecF goodLepProp;

    for (int i=0; i<regMask.size(); i++){
        if (regMask.at(i) == true){
            int l1_idx = l1Idx.at(i);
            int l2_idx = l2Idx.at(i);

            goodLepProp.push_back(lepProp.at(l1_idx));
            goodLepProp.push_back(lepProp.at(l2_idx));
        }
    }

    return goodLepProp;
}
""")

ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
RVecF deltaR(RVecF &lepEta, RVecF &lepPhi){
    RVecF deltaR;

    for (int i=0; i<lepEta.size(); i++){
        float l1_eta = lepEta.at(i);
        float l1_phi = lepPhi.at(i);
        for (int j=1; j<lepEta.size(); j++){
            float l2_eta = lepEta.at(j);
            float l2_phi = lepPhi.at(j);
            float dR = TMath::Power(TMath::Power(l2_eta - l1_eta, 2.) + TMath::Power(l2_phi - l1_phi, 2.), 0.5);
            deltaR.push_back(dR);
        }
    }

    return deltaR;
}
""")

ROOT.gInterpreter.Declare("""
using ROOT::RVecB;
ROOT::VecOps::RVec<short> applyRegFilt(ROOT::VecOps::RVec<short> &prop, RVecB &regMask){
    ROOT::VecOps::RVec<short> reg_prop;
    for (int i=0; i<prop.size(); i++){
        if (regMask.at(i) == true){
            reg_prop.push_back(prop.at(i));
        }
    }
    return reg_prop;
}
""")

ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
RVecF lepFromCand(RVecF &elProp, RVecF &muProp, short Z1l1Idx, short Z1l2Idx, short Z2l1Idx, short Z2l2Idx){
    RVecF lepProp = ROOT::VecOps::Concatenate(elProp, muProp);

    RVecF goodLepProp = ROOT::RVecF {lepProp.at(Z1l1Idx), lepProp.at(Z1l2Idx), lepProp.at(Z2l1Idx), lepProp.at(Z2l2Idx)};

    return goodLepProp;
}
""")
                          
ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
RVecF lepFromCand(RVecF &elProp, RVecF &muProp, short l1Idx, short l2Idx){
    RVecF lepProp = ROOT::VecOps::Concatenate(elProp, muProp);

    RVecF goodLepProp = ROOT::RVecF {lepProp.at(l1Idx), lepProp.at(l2Idx)};

    return goodLepProp;
}
""")

ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
bool passSIP(RVecF &lepSIP){
    bool pass = true;

    for (int i=0; i<lepSIP.size(); i++){
        if (lepSIP.at(i) > 4.) {
            pass = false;
            break;
        }
    }

    return pass;
}
""")
                          
ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
using ROOT::RVecB;
RVecF defWeight(RVecF &ZZCand_dataMCWeight, RVecB &regMask, float overallEventWeight, float genEventSumw){
    RVecF weight;
    float wgt;

    RVecF reg_dataMCWeight = applyRegFilt(ZZCand_dataMCWeight, regMask);
    
    for (int i=0; i<reg_dataMCWeight.size(); i++){
        wgt = reg_dataMCWeight.at(i)*(overallEventWeight/genEventSumw);
        weight.push_back(wgt);
    }

    return weight;
}
""")

ROOT.gInterpreter.Declare("""
using ROOT::RVecF;
using ROOT::RVecI;
RVecF propByFS(RVecF &prop, RVecI Z1flav, RVecI Z2flav, int fs){
    RVecF prop_fs;
                          
    if (fs == -1){
        prop_fs = prop;                      
    }
    else{                 
        for (int i=0; i<prop.size(); i++){
            int cand_fs = Z1flav.at(i)*Z2flav.at(i);
            if (cand_fs == fs){
                prop_fs.push_back(prop.at(i));
            }
        }
    }

    return prop_fs;
}
""")