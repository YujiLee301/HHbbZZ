from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import yaml
import os
from Helper import *
ROOT.PyConfig.IgnoreCommandLineOptions = True


class HZZAnalysisCppProducer(Module):
    def __init__(self,year,cfgFile,isMC,isFSR):
        base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/el9_amd64_gcc12/libJHUGenMELAMELA.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/el9_amd64_gcc12/libjhugenmela.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/el9_amd64_gcc12/libmcfm_710.so" % base)
        ROOT.gSystem.Load("%s/JHUGenMELA/MELA/data/el9_amd64_gcc12/libcollier.so" % base)
        if "/GenAnalysis_cc.so" not in ROOT.gSystem.GetLibraries():
            print("Load GenAnalysis C++ module")
            base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
            if base:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/GenAnalysis.cc+O" % base)
            else:
                base = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/interface/GenAnalysis.h" % base)
        if "/H4LTools_cc.so" not in ROOT.gSystem.GetLibraries():
            print("Load H4LTools C++ module")
            base = "$CMSSW_BASE/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim"
            if base:
                ROOT.gROOT.ProcessLine(
                    ".L %s/src/H4LTools.cc+O" % base)
            else:
                base = "$CMSSW_BASE//src/PhysicsTools/NanoAODTools"
                ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")
                ROOT.gROOT.ProcessLine(
                    ".L %s/interface/H4LTools.h" % base)
        self.year = year
        self.isMC = isMC
        self.genworker = ROOT.GenAnalysis()
        with open(cfgFile, 'r') as ymlfile:
          cfg = yaml.full_load(ymlfile)
          self.worker = ROOT.H4LTools(self.year,self.isMC)
          self.worker.InitializeElecut(cfg['Electron']['pTcut'],cfg['Electron']['Etacut'],cfg['Electron']['Sip3dcut'],cfg['Electron']['Loosedxycut'],cfg['Electron']['Loosedzcut'],
                                       cfg['Electron']['Isocut'],cfg['Electron']['BDTWP']['LowEta']['LowPT'],cfg['Electron']['BDTWP']['MedEta']['LowPT'],cfg['Electron']['BDTWP']['HighEta']['LowPT'],
                                       cfg['Electron']['BDTWP']['LowEta']['HighPT'],cfg['Electron']['BDTWP']['MedEta']['HighPT'],cfg['Electron']['BDTWP']['HighEta']['HighPT'])
          self.worker.InitializeMucut(cfg['Muon']['pTcut'],cfg['Muon']['Etacut'],cfg['Muon']['Sip3dcut'],cfg['Muon']['Loosedxycut'],cfg['Muon']['Loosedzcut'],cfg['Muon']['Isocut'],
                                       cfg['Muon']['Tightdxycut'],cfg['Muon']['Tightdzcut'],cfg['Muon']['TightTrackerLayercut'],cfg['Muon']['TightpTErrorcut'],cfg['Muon']['HighPtBound'])
          self.worker.InitializeFsrPhotonCut(cfg['FsrPhoton']['pTcut'],cfg['FsrPhoton']['Etacut'],cfg['FsrPhoton']['Isocut'],cfg['FsrPhoton']['dRlcut'],cfg['FsrPhoton']['dRlOverPtcut'])
          self.worker.InitializeJetcut(cfg['Jet']['pTcut'],cfg['Jet']['Etacut'],cfg['Jet']['Ncut'])
          self.worker.InitializeEvtCut(cfg['MZ1cut'],cfg['MZZcut'],cfg['Higgscut']['down'],cfg['Higgscut']['up'],cfg['Zmass'],cfg['MZcut']['down'],cfg['MZcut']['up'])
          #self.PUweightfile = cfg["outputdataNPV"]
          #self.PUweighthisto = cfg["PUweightHistoName"]
        #PUinput_file = ROOT.TFile.Open(self.PUweightfile)
        #PUinput_hist = PUinput_file.Get(self.PUweighthisto)
        #self.PUweight_list = []
        #for i in range(1, PUinput_hist.GetNbinsX() + 1):
        #    self.PUweight_list.append(PUinput_hist.GetBinContent(i))
        #PUinput_file.Close()
        self.passtrigEvts = 0
        self.passZZEvts = 0
        self.cfgFile = cfgFile
        self.worker.isFSR = isFSR
        self.print_count = 0
        pass
    def beginJob(self):
        pass

    def endJob(self):
        print(("PassTrig: "+str(self.passtrigEvts)+" Events"))
        print(("Pass4eCut: "+str(self.worker.cut4e)+" Events"))
        print(("Pass4eGhostRemoval: "+str(self.worker.cutghost4e)+" Events"))
        print(("Pass4eLepPtCut: "+str(self.worker.cutLepPt4e)+" Events"))
        print(("Pass4eQCDSupress: "+str(self.worker.cutQCD4e)+" Events"))
        print(("PassmZ1mZ2Cut_4e: "+str(self.worker.cutZZ4e)+" Events"))
        print(("Passm4l_105_160_Cut_4e: "+str(self.worker.cutm4l4e)+" Events"))
        print(("Pass4muCut: "+str(self.worker.cut4mu)+" Events"))
        print(("Pass4muGhostRemoval: "+str(self.worker.cutghost4mu)+" Events"))
        print(("Pass4muLepPtCut: "+str(self.worker.cutLepPt4mu)+" Events"))
        print(("Pass4muQCDSupress: "+str(self.worker.cutQCD4mu)+" Events"))
        print(("PassmZ1mZ2Cut_4mu: "+str(self.worker.cutZZ4mu)+" Events"))
        print(("Passm4l_105_160_Cut_4mu: "+str(self.worker.cutm4l4mu)+" Events"))
        print(("Pass2e2muCut: "+str(self.worker.cut2e2mu)+" Events"))
        print(("Pass2e2muGhostRemoval: "+str(self.worker.cutghost2e2mu)+" Events"))
        print(("Pass2e2muLepPtCut: "+str(self.worker.cutLepPt2e2mu)+" Events"))
        print(("Pass2e2muQCDSupress: "+str(self.worker.cutQCD2e2mu)+" Events"))
        print(("PassmZ1mZ2Cut_2e2mu: "+str(self.worker.cutZZ2e2mu)+" Events"))
        print(("Passm4l_105_160_Cut_2e2mu: "+str(self.worker.cutm4l2e2mu)+" Events"))
        print(("PassZZSelection: "+str(self.passZZEvts)+" Events"))
        if self.isMC:
            print(("PassGEN4eCut: "+str(self.genworker.nGEN4e)+" Events"))
            print(("PassGEN4eZ1Cut: "+str(self.genworker.nGEN4epassZ1)+" Events"))
            print(("PassGEN4efidCut: "+str(self.genworker.nGEN4epassFid)+" Events"))
            print(("PassGEN2e2muCut: "+str(self.genworker.nGEN2e2mu)+" Events"))
            print(("PassGEN2e2muZ1Cut: "+str(self.genworker.nGEN2e2mupassZ1)+" Events"))
            print(("PassGEN2e2mufidCut: "+str(self.genworker.nGEN2e2mupassFid)+" Events"))
            print(("PassGEN4muCut: "+str(self.genworker.nGEN4mu)+" Events"))
            print(("PassGEN4muZ1Cut: "+str(self.genworker.nGEN4mupassZ1)+" Events"))
            print(("PassGEN4mufidCut: "+str(self.genworker.nGEN4mupassFid)+" Events"))
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree)  # initReaders must be called in beginFile
        self.out = wrappedOutputTree
        self.out.branch("mass4l",  "F")
        self.out.branch("GENmass4l",  "F")
        self.out.branch("mass4e",  "F")
        self.out.branch("mass4mu",  "F")
        self.out.branch("mass2e2mu",  "F")
        self.out.branch("pT4l",  "F")
        self.out.branch("GENpT4l",  "F")
        self.out.branch("rapidity4l",  "F")
        self.out.branch("njets_pt30_eta4p7", "I")
        self.out.branch("finalState", "I")
        self.out.branch("GENnjets_pt30_eta4p7", "I")
        self.out.branch("GENrapidity4l",  "F")
        self.out.branch("eta4l",  "F")
        self.out.branch("phi4l",  "F")
        self.out.branch("massZ1",  "F")
        self.out.branch("pTZ1",  "F")
        self.out.branch("etaZ1",  "F")
        self.out.branch("phiZ1",  "F")
        self.out.branch("massZ2",  "F")
        self.out.branch("pTZ2",  "F")
        self.out.branch("etaZ2",  "F")
        self.out.branch("phiZ2",  "F")
        self.out.branch("D_CP",  "F")
        self.out.branch("D_0m",  "F")
        self.out.branch("D_0hp",  "F")
        self.out.branch("D_int",  "F")
        self.out.branch("D_L1",  "F")
        self.out.branch("D_L1Zg",  "F")

        self.out.branch("massL1",  "F")
        self.out.branch("pTL1",  "F")
        self.out.branch("etaL1",  "F")
        self.out.branch("phiL1",  "F")
        self.out.branch("massL2",  "F")
        self.out.branch("pTL2",  "F")
        self.out.branch("etaL2",  "F")
        self.out.branch("phiL2",  "F")
        self.out.branch("massL3",  "F")
        self.out.branch("pTL3",  "F")
        self.out.branch("etaL3",  "F")
        self.out.branch("phiL3",  "F")
        self.out.branch("massL4",  "F")
        self.out.branch("pTL4",  "F")
        self.out.branch("etaL4",  "F")
        self.out.branch("phiL4",  "F")
        self.out.branch("mj1",  "F")
        self.out.branch("pTj1",  "F")
        self.out.branch("etaj1",  "F")
        self.out.branch("phij1",  "F")
        self.out.branch("pTj2",  "F")
        self.out.branch("etaj2",  "F")
        self.out.branch("phij2",  "F")
        self.out.branch("mj2",  "F")
        self.out.branch("btagger1_DJ", "F")
        self.out.branch("btagger1_PN", "F")
        self.out.branch("btagger1_RPT", "F")
        self.out.branch("btagger2_DJ", "F")
        self.out.branch("btagger2_PN", "F")
        self.out.branch("btagger2_RPT", "F")
        self.out.branch("invjj", "F")
        self.out.branch("GENmass2j", "F")
        self.out.branch("GENpTj1", "F")
        self.out.branch("GENetaj1", "F")
        self.out.branch("GENphij1", "F")
        self.out.branch("GENmj1", "F")
        self.out.branch("GENpTj2", "F")
        self.out.branch("GENetaj2", "F")
        self.out.branch("GENphij2", "F")
        self.out.branch("GENmj2", "F")
        self.out.branch("EvtNum",  "I")
        self.out.branch("Weight",  "F")
        self.out.branch("pileupWeight",  "F")
        self.out.branch("dataMCWeight_new",  "F")
        self.out.branch("prefiringWeight",  "F")
        self.out.branch("passedTrig",  "O")
        self.out.branch("passedFullSelection",  "O")
        self.out.branch("passedZ4lSelection",  "O")
        self.out.branch("passedZ4lZ1LSelection",  "O")
        self.out.branch("passedZ4lZXCRSelection",  "O")
        self.out.branch("passedZXCRSelection",  "O")
        self.out.branch("passedFiducialSelection",  "O")
        GENHlepNum = 4
        GENZNum = 2
        GENHjetNum = 2
        self.out.branch("GENlep_MomId",  "I", lenVar = "nGenPart")
        self.out.branch("GENlep_MomMomId",  "I", lenVar = "nGenPart")
        self.out.branch("GENZ_MomId",  "I", lenVar = "nVECZ")
        self.out.branch("GENZ_DaughtersId",  "I", lenVar = "GENZNum")
        self.out.branch("GENlep_Hindex",  "I", lenVar = "GENHlepNum")
        self.out.branch("lep_Hindex",  "I", lenVar = "GENHlepNum")
        self.out.branch("GENlep_id",  "I", lenVar = "nGENLeptons")
        self.out.branch("lep_genindex",  "I", lenVar = "Lepointer")
        self.out.branch("Electron_Fsr_pt",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_eta",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Electron_Fsr_phi",  "F", lenVar = "nElectron_Fsr")
        self.out.branch("Muon_Fsr_pt",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_eta",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("Muon_Fsr_phi",  "F", lenVar = "nMuon_Fsr")
        self.out.branch("goodJet_pt", "F", lenVar="ngoodJets")
        self.out.branch("goodJet_eta", "F", lenVar="ngoodJets")
        self.out.branch("goodJet_phi", "F", lenVar="ngoodJets")
        self.out.branch("goodJet_mass", "F", lenVar="ngoodJets")
        self.out.branch("goodJet_btagRPT", "F", lenVar="ngoodJets")
        self.out.branch("GENjet_hadronFlavour",  "I", lenVar = "nGenJet")
        self.out.branch("GENjet_Hindex",  "I", lenVar = "GENHjetNum")
        self.out.branch("nTightEle", "I")
        self.out.branch("nTightMu", "I")

        with open("SyncLepton2018GGH.txt", 'w') as f:
            f.write("Sync data list:"+"\n")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    # this function gets the pointers to Value and ArrayReaders and sets
    # them in the C++ worker class
    def initReaders(self, tree):
        # self._ttreereaderversion must be set AFTER all calls to
        # tree.valueReader or tree.arrayReader
        self._ttreereaderversion = tree._ttreereaderversion

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail,
        go to next event)"""
        # do this check at every event, as other modules might have read
        # further branches
        #if event._tree._ttreereaderversion > self._ttreereaderversion:
        #    self.initReaders(event._tree)
        # do NOT access other branches in python between the check/call to
        # initReaders and the call to C++ worker code
        self.worker.Initialize()
        isMC = self.isMC
        self.worker.SetObjectNum(event.nElectron,event.nMuon,event.nJet,event.nFsrPhoton)
        if isMC:
            self.worker.SetObjectNumGen(event.nGenPart)
            self.genworker.Initialize()
            self.genworker.SetObjectNumGen(event.nGenPart, event.nGenJet)
        keepIt = False
        Lepointer = 0
        EvtNum = 0
        Weight = 1
        passedTrig=False
        passedFullSelection=False
        passedZ4lSelection=False
        passedQCDcut=False
        passedZ1LSelection=False
        passedZ4lZ1LSelection=False
        passedZ4lZXCRSelection=False
        passedZXCRSelection=False
        passedFiducialSelection=False
        nZXCRFailedLeptons=0
        prefiringWeight = 1
        dataMCWeight_new = 1
        pileupWeight = 1
        mass4e=0
        mass2e2mu=0
        mass4mu=0
        finalState=-1
        GENmass4l = -99
        GENpT4l = -99
        nVECZ = 2
        GENrapidity4l = -99
        GENnjets_pt30_eta4p7 = -1
        nGENLeptons = 0
        nGenPart = 0
        nGenJet = 0
        pTZ1 = -99
        etaZ1 = -99
        phiZ1 = -99
        massZ1 = 0
        pTZ2 = -99
        etaZ2 = -99
        phiZ2 = -99
        massZ2 = -99
        pT4l = -99
        eta4l = -99
        phi4l = -99
        mass4l = 0
        rapidity4l = -99
        btagger1_DJ = -99
        btagger1_PN = -99
        btagger1_RPT = -99
        btagger2_DJ = -99
        btagger2_PN = -99
        btagger2_RPT = -99
        invjj = -99
        GENmass2j = -99
        GENpTj1 = -99
        GENetaj1 = -99
        GENphij1 = -99
        GENmj1 = -99
        GENpTj2 = -99
        GENetaj2 = -99
        GENphij2 = -99
        GENmj2 = -99
        passedTrig = PassTrig(event, self.cfgFile)
        nTightEle = 0
        nTightMu = 0
        if (passedTrig==True):
            self.passtrigEvts += 1
        else:
            return keepIt
        if(isMC):
            pileupWeight = 1.0
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        fsrPhotons = Collection(event, "FsrPhoton")
        jets = Collection(event, "Jet")
        if isMC:
            nGenPart = event.nGenPart
            genparts = Collection(event, "GenPart")
            genjets = Collection(event, "GenJet")
            for xj in genjets:
                self.genworker.SetGenJets(xj.pt,xj.eta,xj.phi,xj.mass,xj.hadronFlavour)
            for xg in genparts:
                self.worker.SetGenParts(xg.pt)
                self.genworker.SetGenParts(xg.pt,xg.eta,xg.phi,xg.mass,xg.pdgId,xg.status,xg.statusFlags,xg.genPartIdxMother)
            for xm in muons:
                self.worker.SetMuonsGen(xm.genPartIdx)
            for xe in electrons:
                self.worker.SetElectronsGen(xe.genPartIdx)
                
        branches = [b.GetName() for b in event._tree.GetListOfBranches()]        
        for xe in electrons:
            hasFall17 = "Electron_mvaFall17V2Iso_WP80" in branches
            hasRun3   = "Electron_mvaIso_WP80" in branches
            
            wp80_Fall17 = xe.mvaFall17V2Iso_WP80 if hasFall17 else False
            wp90_Fall17 = xe.mvaFall17V2Iso_WP90 if hasFall17 else False
            wpl_Fall17  = xe.mvaFall17V2Iso_WPL  if hasFall17 else False
            wp80_Run3   = xe.mvaIso_WP80 if hasRun3 else False
            wp90_Run3   = xe.mvaIso_WP90 if hasRun3 else False
            self.worker.SetElectrons(xe.pt, xe.eta, xe.phi, xe.mass, xe.dxy,xe.dz, xe.sip3d, wp80_Fall17, wp90_Fall17, wpl_Fall17,
                                     wp80_Run3, wp90_Run3, xe.pdgId, xe.pfRelIso03_all)
        
        branches = [b.GetName() for b in event._tree.GetListOfBranches()]        
        for xm in muons:
            hasCutBased = "Muon_looseId" in branches
            hasLowPtMVA = "Muon_mvaLowPtId" in branches #run2
            hasMVA_Run2 = "Muon_mvaId" in branches #run2
            hasMVA_Run3 = "Muon_mvaLowPt" in branches #run3
            
            looseId  = xm.looseId  if hasCutBased else False
            mediumId = xm.mediumId if hasCutBased else False
            tightId  = xm.tightId  if hasCutBased else False
            
            mvaLowPtId = xm.mvaLowPtId if hasLowPtMVA else 0
            mvaId      = xm.mvaId      if hasMVA_Run2 else 0
            mvaWP      = xm.mvaMuID_WP if hasMVA_Run3 else 0  # 0=fail,1=MVAIDwpMedium,2=MVAIDwpTight
            self.worker.SetMuons(xm.pt, xm.eta, xm.phi, xm.mass, xm.isGlobal, xm.isTracker,
                                xm.dxy, xm.dz, xm.sip3d, xm.ptErr, looseId, mediumId, tightId,
                                mvaLowPtId, mvaId, mvaWP,
                                xm.nTrackerLayers, xm.isPFcand, xm.pdgId, xm.charge, xm.pfRelIso03_all)
            
        for xf in fsrPhotons:
            self.worker.SetFsrPhotons(xf.dROverEt2,xf.eta,xf.phi,xf.pt,xf.relIso03,xf.electronIdx,xf.muonIdx)
        for xj in jets:
            self.worker.SetJets(xj.pt,xj.eta,xj.phi,xj.mass,xj.jetId,xj.btagDeepFlavB,xj.btagPNetB,xj.btagRobustParTAK4B, 0.8, 7)
        self.worker.BatchFsrRecovery_Run3()
        
        self.worker.LeptonSelection()
        if ((self.worker.nTightEle<2)&(self.worker.nTightMu<2)):
            pass
        if isMC:
            self.genworker.SetGenVariables()
            GENmass4l = self.genworker.GENmass4l
            GENpT4l = self.genworker.GENpT4l
            GENrapidity4l = self.genworker.GENrapidity4l
            GENnjets_pt30_eta4p7 = self.genworker.GENnjets_pt30_eta4p7
            nGENLeptons = self.genworker.nGENLeptons
            
            GENmass2j = self.genworker.GENmass2j
            GENpTj1 = self.genworker.GENpTj1
            GENetaj1 = self.genworker.GENetaj1
            GENphij1 = self.genworker.GENphij1
            GENmj1 = self.genworker.GENmj1
            GENpTj2 = self.genworker.GENpTj2
            GENetaj2 = self.genworker.GENetaj2
            GENphij2 = self.genworker.GENphij2
            GENmj2 = self.genworker.GENmj2
        
        passedFiducialSelection = self.genworker.passedFiducialSelection

        Electron_Fsr_pt_vec = self.worker.ElectronFsrPt()
        Electron_Fsr_eta_vec = self.worker.ElectronFsrEta()
        Electron_Fsr_phi_vec = self.worker.ElectronFsrPhi()
        Muon_Fsr_pt_vec = self.worker.MuonFsrPt()
        Muon_Fsr_eta_vec = self.worker.MuonFsrEta()
        Muon_Fsr_phi_vec = self.worker.MuonFsrPhi()
        Electron_Fsr_pt = []
        Electron_Fsr_eta = []
        Electron_Fsr_phi = []
        Muon_Fsr_pt = []
        Muon_Fsr_eta = []
        Muon_Fsr_phi = []
        
        if len(Electron_Fsr_pt_vec)>0:
            for i in range(len(Electron_Fsr_pt_vec)):
                Electron_Fsr_pt.append(Electron_Fsr_pt_vec[i])
                Electron_Fsr_eta.append(Electron_Fsr_eta_vec[i])
                Electron_Fsr_phi.append(Electron_Fsr_phi_vec[i])
        if len(Muon_Fsr_pt_vec)>0:
            for i in range(len(Muon_Fsr_pt_vec)):
                Muon_Fsr_pt.append(Muon_Fsr_pt_vec[i])
                Muon_Fsr_eta.append(Muon_Fsr_eta_vec[i])
                Muon_Fsr_phi.append(Muon_Fsr_phi_vec[i])
                
        goodJet_idxs = self.worker.jetidx
        goodJet_pt = []
        goodJet_eta = []
        goodJet_phi = []
        goodJet_mass = []
        goodJet_btagRPT = []

        for idx in goodJet_idxs:
            goodJet_pt.append(jets[idx].pt)
            goodJet_eta.append(jets[idx].eta)
            goodJet_phi.append(jets[idx].phi)
            goodJet_mass.append(jets[idx].mass)
            goodJet_btagRPT.append(jets[idx].btagRobustParTAK4B)

        #if self.print_count < 20:
        #    print(f"[Python] Event {self.print_count} - Selected Jets (index | pt | eta | phi | btag):")
        #    for idx in goodJet_idxs:
        #        jet = jets[idx]
        #        print(f"  idx = {idx}, pt = {jet.pt:.2f}, eta = {jet.eta:.2f}, phi = {jet.phi:.2f}, btag = {jet.btagRobustParTAK4B:.3f}")
        #    self.print_count += 1
    
        GENlep_id = []
        GENlep_Hindex = []
        GENZ_DaughtersId = []
        GENZ_MomId = []
        GENlep_MomId = []
        GENlep_MomMomId = []
        GENjet_hadronFlavour = []
        GENjet_Hindex = []
        if isMC:
            GENlep_id_vec = self.genworker.GENlep_id
            if len(GENlep_id_vec)>0:
                for i in range(len(GENlep_id_vec)):
                    GENlep_id.append(GENlep_id_vec[i])
            GENlep_Hindex_vec = self.genworker.GENlep_Hindex
            if len(GENlep_Hindex_vec)>0:
                for i in range(len(GENlep_Hindex_vec)):
                    GENlep_Hindex.append(GENlep_Hindex_vec[i])   
            GENZ_DaughtersId_vec = self.genworker.GENZ_DaughtersId
            if len(GENZ_DaughtersId_vec)>0:
                for i in range(len(GENZ_DaughtersId_vec)):
                    GENZ_DaughtersId.append(GENZ_DaughtersId_vec[i])
            nVECZ = self.genworker.nVECZ
            GENZ_MomId_vec = self.genworker.GENZ_MomId
            if len(GENZ_MomId_vec)>0:
                for i in range(len(GENZ_MomId_vec)):
                    GENZ_MomId.append(GENZ_MomId_vec[i])
            GENlep_MomId_vec = self.genworker.GENlep_MomId
            if len(GENlep_MomId_vec)>0:
                for i in range(len(GENlep_MomId_vec)):
                    GENlep_MomId.append(GENlep_MomId_vec[i])
            GENlep_MomMomId_vec = self.genworker.GENlep_MomMomId
            if len(GENlep_MomMomId_vec)>0:
                for i in range(len(GENlep_MomMomId_vec)):
                    GENlep_MomMomId.append(GENlep_MomMomId_vec[i])
            GENjet_hadronFlavour_vec = self.genworker.GENjet_hadronFlavour
            if len(GENjet_hadronFlavour_vec)>0:
                for i in range(len(GENjet_hadronFlavour_vec)):
                    GENjet_hadronFlavour.append(GENjet_hadronFlavour_vec[i])
            GENjet_Hindex_vec = self.genworker.GENjet_Hindex
            if len(GENjet_Hindex_vec)>0:
                for i in range(len(GENjet_Hindex_vec)):
                    GENjet_Hindex.append(GENjet_Hindex_vec[i])

        foundZZCandidate = self.worker.ZZSelection()
        passedFullSelection=foundZZCandidate
        Lepointer = self.worker.Lepointer
        lep_Hindex = []
        lep_Hindex_vec = self.worker.lep_Hindex
        if len(lep_Hindex_vec)>0:
            for i in range(len(lep_Hindex_vec)):
                lep_Hindex.append(lep_Hindex_vec[i])
        lep_genindex = []
        if isMC:
            lep_genindex_vec = self.worker.lep_genindex
            if len(lep_genindex_vec)>0:
                for i in range(len(lep_genindex_vec)):
                    lep_genindex.append(lep_genindex_vec[i])
        if (foundZZCandidate):
            self.passZZEvts += 1
        #if (foundZZCandidate | passedFiducialSelection):
        if (foundZZCandidate):
        #if (passedFiducialSelection):
            EvtNum += 1
            keepIt = True
        if self.worker.RecoFourMuEvent: finalState = 1
        if self.worker.RecoFourEEvent: finalState = 2
        if self.worker.RecoTwoETwoMuEvent: finalState = 3
        if self.worker.RecoTwoMuTwoEEvent: finalState = 4
        pTZ1 = self.worker.Z1.Pt()
        etaZ1 = self.worker.Z1.Eta()
        phiZ1 = self.worker.Z1.Phi()
        massZ1 = self.worker.Z1.M()
        pTZ2 = self.worker.Z2.Pt()
        etaZ2 = self.worker.Z2.Eta()
        phiZ2 = self.worker.Z2.Phi()
        massZ2 = self.worker.Z2.M()
        D_CP = self.worker.D_CP
        D_0m = self.worker.D_0m
        D_0hp = self.worker.D_0hp
        D_int = self.worker.D_int
        D_L1 = self.worker.D_L1
        D_L1Zg = self.worker.D_L1Zg

        pTL1 = self.worker.pTL1
        etaL1 = self.worker.etaL1
        phiL1 = self.worker.phiL1
        massL1 = self.worker.massL1
        pTL2 = self.worker.pTL2
        etaL2 = self.worker.etaL2
        phiL2 = self.worker.phiL2
        massL2 = self.worker.massL2
        pTL3 = self.worker.pTL3
        etaL3 = self.worker.etaL3
        phiL3 = self.worker.phiL3
        massL3 = self.worker.massL3
        pTL4 = self.worker.pTL4
        etaL4 = self.worker.etaL4
        phiL4 = self.worker.phiL4
        massL4 = self.worker.massL4
        pTj1 = self.worker.pTj1
        etaj1 = self.worker.etaj1
        phij1 = self.worker.phij1
        mj1 = self.worker.mj1
        pTj2 = self.worker.pTj2
        etaj2 = self.worker.etaj2
        phij2 = self.worker.phij2
        mj2 = self.worker.mj2
        btagger1_DJ = self.worker.btagger1_DJ
        btagger1_PN = self.worker.btagger1_PN
        btagger1_RPT = self.worker.btagger1_RPT
        btagger2_DJ = self.worker.btagger2_DJ
        btagger2_PN = self.worker.btagger2_PN
        btagger2_RPT = self.worker.btagger2_RPT
        invjj = self.worker.invjj

        if pTL2>pTL1:
            pTL1, pTl2 = pTL2, pTL1
            etaL1, etaL2 = etaL2, etaL1
            phiL1, phiL2 = phiL2, phiL1
            massL1,massL2 = massL2, massL1
        if pTL4>pTL3:
            pTL3, pTL4 = pTL4, pTL3
            etaL3, etaL4 = etaL4, etaL3
            phiL3, phiL4 = phiL4, phiL3
            massL3, massL4 = massL4, massL3
        if passedFullSelection: 
            pT4l = self.worker.ZZsystem.Pt()
            eta4l = self.worker.ZZsystem.Eta()
            phi4l = self.worker.ZZsystem.Phi()
            mass4l = self.worker.ZZsystem.M()
            rapidity4l = self.worker.ZZsystem.Rapidity()
        njets_pt30_eta4p7 = self.worker.njets_pt30_eta4p7
        if self.worker.flag4e:
            mass4e = mass4l
        if self.worker.flag2e2mu:
            mass2e2mu = mass4l
        if self.worker.flag4mu:
            mass4mu = mass4l
        if (self.worker.isFSR==False & passedFullSelection):
            pT4l = self.worker.ZZsystemnofsr.Pt()
            eta4l = self.worker.ZZsystemnofsr.Eta()
            phi4l = self.worker.ZZsystemnofsr.Phi()
            mass4l = self.worker.ZZsystemnofsr.M()
            rapidity4l = self.worker.ZZsystemnofsr.Rapidity()
        Weight = event.genWeight * pileupWeight * dataMCWeight_new * prefiringWeight
        self.out.fillBranch("mass4l",mass4l)
        self.out.fillBranch("GENmass4l",GENmass4l)
        self.out.fillBranch("mass4e",mass4e)
        self.out.fillBranch("mass2e2mu",mass2e2mu)
        self.out.fillBranch("mass4mu",mass4mu)
        self.out.fillBranch("pT4l",pT4l)
        self.out.fillBranch("GENpT4l",GENpT4l)
        self.out.fillBranch("rapidity4l",rapidity4l)
        self.out.fillBranch("GENrapidity4l",GENrapidity4l)
        self.out.fillBranch("njets_pt30_eta4p7",njets_pt30_eta4p7)
        self.out.fillBranch("finalState",finalState)
        self.out.fillBranch("GENnjets_pt30_eta4p7",GENnjets_pt30_eta4p7)
        self.out.fillBranch("eta4l",eta4l)
        self.out.fillBranch("phi4l",phi4l)
        self.out.fillBranch("massZ1",massZ1)
        self.out.fillBranch("pTZ1",pTZ1)
        self.out.fillBranch("etaZ1",etaZ1)
        self.out.fillBranch("phiZ1",phiZ1)
        self.out.fillBranch("massZ2",massZ2)
        self.out.fillBranch("pTZ2",pTZ2)
        self.out.fillBranch("etaZ2",etaZ2)
        self.out.fillBranch("phiZ2",phiZ2)
        self.out.fillBranch("D_CP",D_CP)
        self.out.fillBranch("D_0m",D_0m)
        self.out.fillBranch("D_0hp",D_0hp)
        self.out.fillBranch("D_int",D_int)
        self.out.fillBranch("D_L1",D_L1)
        self.out.fillBranch("D_L1Zg",D_L1Zg)
        self.out.fillBranch("passedTrig",  passedTrig)
        self.out.fillBranch("passedFullSelection",  passedFullSelection)
        self.out.fillBranch("passedZ4lSelection", passedZ4lSelection)
        self.out.fillBranch("passedZ4lZ1LSelection",  passedZ4lZ1LSelection)
        self.out.fillBranch("passedZ4lZXCRSelection",  passedZ4lZXCRSelection)
        self.out.fillBranch("passedZXCRSelection",  passedZXCRSelection)
        self.out.fillBranch("passedFiducialSelection",  passedFiducialSelection)
        self.out.fillBranch("EvtNum",EvtNum)
        self.out.fillBranch("massL1",massL1)
        self.out.fillBranch("pTL1",pTL1)
        self.out.fillBranch("etaL1",etaL1)
        self.out.fillBranch("phiL1",phiL1)
        self.out.fillBranch("massL2",massL2)
        self.out.fillBranch("pTL2",pTL2)
        self.out.fillBranch("etaL2",etaL2)
        self.out.fillBranch("phiL2",phiL2)
        self.out.fillBranch("massL3",massL3)
        self.out.fillBranch("pTL3",pTL3)
        self.out.fillBranch("etaL3",etaL3)
        self.out.fillBranch("phiL3",phiL3)
        self.out.fillBranch("massL4",massL4)
        self.out.fillBranch("pTL4",pTL4)
        self.out.fillBranch("etaL4",etaL4)
        self.out.fillBranch("phiL4",phiL4)

        self.out.fillBranch("mj1",mj1)
        self.out.fillBranch("pTj1",pTj1)
        self.out.fillBranch("etaj1",etaj1)
        self.out.fillBranch("phij1",phij1)
        self.out.fillBranch("mj2",mj2)
        self.out.fillBranch("pTj2",pTj2)
        self.out.fillBranch("etaj2",etaj2)
        self.out.fillBranch("phij2",phij2)
        self.out.fillBranch("btagger1_DJ",btagger1_DJ)
        self.out.fillBranch("btagger1_PN",btagger1_PN)
        self.out.fillBranch("btagger1_RPT",btagger1_RPT)
        self.out.fillBranch("btagger2_DJ",btagger2_DJ)
        self.out.fillBranch("btagger2_PN",btagger2_PN)
        self.out.fillBranch("btagger2_RPT",btagger2_RPT)
        self.out.fillBranch("invjj",invjj)
        self.out.fillBranch("GENmass2j",GENmass2j)
        self.out.fillBranch("GENpTj1",GENpTj1)
        self.out.fillBranch("GENetaj1",GENetaj1)
        self.out.fillBranch("GENphij1",GENphij1)
        self.out.fillBranch("GENmj1",GENmj1)
        self.out.fillBranch("GENpTj2",GENpTj2)
        self.out.fillBranch("GENetaj2",GENetaj2)
        self.out.fillBranch("GENphij2",GENphij2)
        self.out.fillBranch("GENmj2",GENmj2)
        
        self.out.fillBranch("pileupWeight",pileupWeight)
        self.out.fillBranch("dataMCWeight_new",dataMCWeight_new)
        self.out.fillBranch("prefiringWeight",prefiringWeight)
        self.out.fillBranch("Weight",Weight)
        # self.out.fillBranch("nElectron_Fsr", len(electrons))
        # self.out.fillBranch("nMuon_Fsr", len(muons))
        
        self.out.fillBranch("GENlep_id",GENlep_id)
        self.out.fillBranch("GENlep_Hindex",GENlep_Hindex)
        self.out.fillBranch("GENZ_DaughtersId",GENZ_DaughtersId)
        self.out.fillBranch("GENZ_MomId",GENZ_MomId)
        self.out.fillBranch("GENlep_MomId",GENlep_MomId)
        self.out.fillBranch("GENlep_MomMomId",GENlep_MomMomId)
        self.out.fillBranch("Electron_Fsr_pt",Electron_Fsr_pt)
        self.out.fillBranch("Electron_Fsr_eta",Electron_Fsr_eta)
        self.out.fillBranch("Electron_Fsr_phi",Electron_Fsr_phi)

        self.out.fillBranch("lep_Hindex",lep_Hindex)
        self.out.fillBranch("lep_genindex",lep_genindex)
        self.out.fillBranch("Muon_Fsr_pt",Muon_Fsr_pt)
        self.out.fillBranch("Muon_Fsr_eta",Muon_Fsr_eta)
        self.out.fillBranch("Muon_Fsr_phi",Muon_Fsr_phi)
        self.out.fillBranch("goodJet_pt", goodJet_pt)
        self.out.fillBranch("goodJet_eta", goodJet_eta)
        self.out.fillBranch("goodJet_phi", goodJet_phi)
        self.out.fillBranch("goodJet_mass", goodJet_mass)
        self.out.fillBranch("goodJet_btagRPT", goodJet_btagRPT)
        self.out.fillBranch("GENjet_hadronFlavour",GENjet_hadronFlavour)
        self.out.fillBranch("GENjet_Hindex",GENjet_Hindex)
        self.out.fillBranch("nTightEle", self.worker.nTightEle)
        self.out.fillBranch("nTightMu", self.worker.nTightMu)

        """with open("SyncLepton2018GGH.txt", 'a') as f:
            if(foundZZCandidate):
                f.write(str('%.4f' % event.run)+":"+str('%.4f' % event.luminosityBlock)+":"+str('%.4f' % event.event)+":" \
                        +str('%.4f' % self.worker.pTL1)+":"+str('%.4f' % self.worker.etaL1)+":"+str('%.4f' % self.worker.phiL1)+":"+str('%.4f' % self.worker.massL1)+":" \
                        +str('%.4f' % self.worker.pTL2)+":"+str('%.4f' % self.worker.etaL2)+":"+str('%.4f' % self.worker.phiL2)+":"+str('%.4f' % self.worker.massL2)+":" \
                        +str('%.4f' % self.worker.pTL3)+":"+str('%.4f' % self.worker.etaL3)+":"+str('%.4f' % self.worker.phiL3)+":"+str('%.4f' % self.worker.massL3)+":" \
                        +str('%.4f' % self.worker.pTL4)+":"+str('%.4f' % self.worker.etaL4)+":"+str('%.4f' % self.worker.phiL4)+":"+str('%.4f' % self.worker.massL4)+"\n")
            else:
                f.write(str('%.4f' % event.run)+":"+str('%.4f' % event.luminosityBlock)+":"+str('%.4f' % event.event)+":" \
                        +str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":" \
                        +str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":" \
                        +str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":" \
                        +str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+":"+str('%.4f'%-1.0000)+"\n")"""




        return keepIt


# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed

#H4LCppModule() = lambda: HZZAnalysisCppProducer(year)
