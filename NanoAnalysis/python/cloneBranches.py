from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
# Debugging
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from json import dump

class cloneBranches(Module):

    def __init__(self, treeName, varlist, continueFor = lambda evt : (True)) :
        '''Clone the branches listed in "varlist" to a new tree named "treeName" in the output file,
        and stop processing events that do not satisfy the condition "continueFor".
        This is useul in particular to store gen variables for all events in a separate tree.
        '''

        self.treeName = treeName
        self.varlist = varlist
        self.continueFor = continueFor
        print("***cloneBranches: cloning into tree:", treeName, "- This module filters events.", flush=True )

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.inputTree = inputTree
        # Clone the requested branches into the new tree.
        eventsTree = wrappedOutputTree.tree()
        eventsTree.SetBranchStatus("*", 0)
        for var in self.varlist :
            eventsTree.SetBranchStatus(var, 1)
        self.newTree = eventsTree.CloneTree(0)
        self.newTree.SetName(self.treeName)
        eventsTree.SetBranchStatus("*", 1)
        self.failed_lhe_info = dict(
            el_pt    = [],
            el_eta   = [],
            mu_pt    = [],
            mu_eta   = []
        )
        self.passed_lhe_info = dict(
            el_pt  = [],
            el_eta = [],
            mu_pt  = [],
            mu_eta = []
        )

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
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree) :
        # print("\n========================\n")
        # print("Total cloneBranches Failures: ", self.all_failures)
        # print("Tau Fail Percentage: ", self.tau_failures/self.all_failures)
        # print("\n========================\n")
        # with open('lheInfo_cloneBranchFailures_noTau_v4.json', 'w') as outfile:
        #     dump(self.failed_lhe_info, outfile, indent=4)
        # with open('lheInfo_cloneBranchPasses_noTau.json', 'w') as outfile:
        #     dump(self.passed_lhe_info, outfile, indent=4)
        self.newTree.Write()
        
    def analyze(self, event) :
        self.inputTree.readAllBranches()        
        self.newTree.Fill()

        # lhepart = Collection(event, "LHEPart")
        # if not self.continueFor(event):
        #     self.fill_dict(self.failed_lhe_info, lhepart)
        #     # lhepart = Collection(event, "LHEPart")
        #     # for i in range(4, 8):
        #     #     lep = lhepart[i]
        #     #     if abs(lep.pdgId) == 11:                    
        #     #         self.failed_lhe_info["el_pt"].append(lep.pt)
        #     #         self.failed_lhe_info["el_eta"].append(lep.eta)
        #     #     elif abs(lep.pdgId) == 13:
        #     #         self.failed_lhe_info["mu_pt"].append(lep.pt)
        #     #         self.failed_lhe_info["mu_eta"].append(lep.eta)
        #         # lep_idx = i+4
        #         # pt_key = "lep" + str(i+1) + "_pt"
        #         # eta_key = "lep" + str(i+1) + "_eta"
        #         # self.failed_lhe_info[pt_key].append(lhepart[lep_idx].pt)
        #         # self.failed_lhe_info[eta_key].append(lhepart[lep_idx].eta)
        # else:
        #     self.fill_dict(self.passed_lhe_info, lhepart)
  
        return self.continueFor(event) 
