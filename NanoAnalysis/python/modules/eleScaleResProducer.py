import os

# Set up the NATModules eleScaleRes module
def getEleScaleRes(era, tag, is_mc, overwritePt=True) :
    from PhysicsTools.NATModules.modules.eleScaleRes import eleScaleRes

    # Check for supported eras
    if era not in [2022, 2023]:
        raise ValueError(f"getEleScaleRes: Era {era} not supported")

    if era == 2022:
        scaleKey = "Scale"
        smearKey = "Smearing" if is_mc else None

        if "pre_EE" in tag :
            fname = "electronSS_preEE.json.gz"
        else:
            fname = "electronSS_postEE.json.gz"

    elif era == 2023:
        if "preBPix" in tag:
            scaleKey = "2023PromptC_ScaleJSON"
            smearKey = "2023PromptC_SmearingJSON" if is_mc else None
            fname = "electronSS_preBPix.json.gz"
        else:
            scaleKey = "2023PromptD_ScaleJSON"
            smearKey = "2023PromptD_SmearingJSON" if is_mc else None
            fname = "electronSS_postBPix.json.gz"
 
    json = "%s/src/ZZAnalysis/NanoAnalysis/data/%s" % (os.environ['CMSSW_BASE'], fname)

    print("***eleScaleRes: era:", era, "tag:", tag, "is MC:", is_mc, "overwritePt:", overwritePt, "json:", json)
    return eleScaleRes(json, scaleKey, smearKey, overwritePt)
