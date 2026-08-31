#ifndef H4LTools_h
#define H4LTools_h

#include <utility>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TSpline.h>
#include <vector>
//#include "yaml-cpp/yaml.h"
#include "../JHUGenMELA/MELA/interface/Mela.h"

class H4LTools {
    public:
      H4LTools(int year, bool isMC_);
      float elePtcut, MuPtcut, eleEtacut, MuEtacut, elesip3dCut, Musip3dCut,Zmass,MZ1cut,MZcutup,MZcutdown,MZZcut,HiggscutUp,HiggscutDown;
      float eleLoosedxycut,eleLoosedzcut,MuLoosedxycut,MuLoosedzcut,MuTightdxycut,MuTightdzcut,MuTightTrackerLayercut,MuTightpTErrorcut,MuHighPtBound,eleIsocut,MuIsocut;
      float fsrphotonPtcut,fsrphotonEtacut,fsrphotonIsocut,fsrphotondRlcut,fsrphotondRlOverPtcut,JetPtcut,JetEtacut,JetNcut;
      float btagger1_DJ,btagger1_PN,btagger1_RPT,btagger1_UPT,btagger2_DJ,btagger2_PN,btagger2_RPT,btagger2_UPT,invjj;
      float eleBDTWPLELP,eleBDTWPMELP,eleBDTWPHELP,eleBDTWPLEHP,eleBDTWPMEHP,eleBDTWPHEHP;
      bool RecoFourMuEvent, RecoFourEEvent, RecoTwoETwoMuEvent, RecoTwoMuTwoEEvent;

      int cut4e, cut4mu, cut2e2mu;
      int cutZZ4e, cutZZ4mu, cutZZ2e2mu, cutm4l4e, cutm4l4mu, cutm4l2e2mu, cutghost2e2mu, cutQCD2e2mu, cutLepPt2e2mu, cutghost4e, cutQCD4e, cutLepPt4e, cutghost4mu, cutQCD4mu, cutLepPt4mu, Z1flav; //2l2j
    
      // 2l2j cutflow counters
      int passTwoTightLeps;
      int passZCand;
      int passAtLeastTwoRawJets;
      int passAtLeastTwoPtEtaJets;
      int passAtLeastTwoJetIdJets;
      int passAtLeastTwoPuIdJets;
      int passTwoGoodJets;
      int passDijet;
      int passFinal;

      // event-level jet multiplicities for cumulative jet cutflow
      int nRawJetsThisEvent;
      int nPtEtaJetsThisEvent;
      int nJetIdJetsThisEvent;
      int nPuIdJetsThisEvent;

      // event-level 2l2j cutflow flags
      bool eventPassTwoTightLeps;
      bool eventPassZCand;
      bool eventPassAtLeastTwoRawJets;
      bool eventPassAtLeastTwoPtEtaJets;
      bool eventPassAtLeastTwoJetIdJets;
      bool eventPassAtLeastTwoPuIdJets;
      bool eventPassTwoGoodJets;
      bool eventPassDijet;
      bool eventPassFinal;

      void InitializeElecut(float elePtcut_,float eleEtacut_,float elesip3dCut_,float eleLoosedxycut_,float eleLoosedzcut_,float eleIsocut_,float eleBDTWPLELP_,float eleBDTWPMELP_, float eleBDTWPHELP_,float eleBDTWPLEHP_,float eleBDTWPMEHP_,float eleBDTWPHEHP_){
        elePtcut = elePtcut_;
        eleEtacut = eleEtacut_;
        elesip3dCut = elesip3dCut_;
        eleLoosedxycut = eleLoosedxycut_;
        eleLoosedzcut = eleLoosedzcut_;
        eleIsocut = eleIsocut_;
        eleBDTWPLELP = eleBDTWPLELP_;
        eleBDTWPMELP = eleBDTWPMELP_;
        eleBDTWPHELP = eleBDTWPHELP_;
        eleBDTWPLEHP = eleBDTWPLEHP_;
        eleBDTWPMEHP = eleBDTWPMEHP_;
        eleBDTWPHEHP = eleBDTWPHEHP_;
      }
      void InitializeMucut(float MuPtcut_,float MuEtacut_,float Musip3dCut_,float MuLoosedxycut_,float MuLoosedzcut_,float MuIsocut_,float MuTightdxycut_,float MuTightdzcut_,float MuTightTrackerLayercut_,float MuTightpTErrorcut_,float MuHighPtBound_){
        MuPtcut = MuPtcut_;
        MuEtacut = MuEtacut_;
        Musip3dCut = Musip3dCut_;
        MuLoosedxycut = MuLoosedxycut_;
        MuLoosedzcut = MuLoosedzcut_;
        MuIsocut = MuIsocut_;
        MuTightdxycut = MuTightdxycut_;
        MuTightdzcut = MuTightdzcut_;
        MuTightTrackerLayercut = MuTightTrackerLayercut_;
        MuTightpTErrorcut = MuTightpTErrorcut_;
        MuHighPtBound = MuHighPtBound_;
      }
      void InitializeFsrPhotonCut(float fsrphotonPtcut_, float fsrphotonEtacut_, float fsrphotonIsocut_, float fsrphotondRlcut_, float fsrphotondRlOverPtcut_){
        fsrphotonPtcut = fsrphotonPtcut_;
        fsrphotonEtacut = fsrphotonEtacut_;
        fsrphotonIsocut = fsrphotonIsocut_;
        fsrphotondRlcut = fsrphotondRlcut_;
        fsrphotondRlOverPtcut = fsrphotondRlOverPtcut_;
      }
      void InitializeJetcut(float JetPtcut_, float JetEtacut_, float JetNcut_){
        JetPtcut = JetPtcut_;
        JetEtacut = JetEtacut_;
        JetNcut = JetNcut_;
      }
      void InitializeEvtCut(float MZ1cut_,float MZZcut_,float HiggscutDown_,float HiggscutUp_,float Zmass_,float MZcutdown_, float MZcutup_){
        MZ1cut = MZ1cut_;
        MZZcut = MZZcut_;
        HiggscutDown = HiggscutDown_;
        HiggscutUp = HiggscutUp_;
        Zmass = Zmass_;
        MZcutdown = MZcutdown_;
        MZcutup = MZcutup_;
      }
      //void SetElectrons(float Electron_pt_, float Electron_eta_, float Electron_phi_, float Electron_mass_, float Electron_dxy_, float Electron_dz_, float Electron_sip3d_, 
      //                  float Electron_deltaEtaSC_, float Electron_mvaHZZIso_, bool Electron_mvaIso_WPHZZ_, int Electron_pdgId_, float Electron_pfRelIso03_all_){
      void SetElectrons(float Electron_pt_, float Electron_eta_, float Electron_phi_, float Electron_mass_, float Electron_dxy_, float Electron_dz_, float Electron_sip3d_, 
                        float Electron_deltaEtaSC_, float Electron_mvaNoIso_, bool Electron_mvaIso_WPHZZ_, int Electron_pdgId_, float Electron_pfRelIso03_all_){
        Electron_pt.push_back(Electron_pt_); 
        Electron_phi.push_back(Electron_phi_);
        Electron_eta.push_back(Electron_eta_);
        Electron_mass.push_back(Electron_mass_);
        Electron_dxy.push_back(Electron_dxy_);
        Electron_dz.push_back(Electron_dz_);
        Electron_sip3d.push_back(Electron_sip3d_);
        Electron_deltaEtaSC.push_back(Electron_deltaEtaSC_);
        //Electron_mvaHZZIso.push_back(Electron_mvaHZZIso_);
        Electron_mvaNoIso.push_back(Electron_mvaNoIso_);
        Electron_mvaIso_WPHZZ.push_back(Electron_mvaIso_WPHZZ_);
        Electron_pdgId.push_back(Electron_pdgId_);
        Electron_pfRelIso03_all.push_back(Electron_pfRelIso03_all_);
      }

      void SetJets(float Jet_pt_, float Jet_eta_, float Jet_phi_, float Jet_mass_, int Jet_jetId_, float Jet_btagDeepFlavB_, float Jet_btagPNetB_, float Jet_btagRobustParTAK4B_, float Jet_btagUParTAK4B_, float Jet_btagDeepC_,
                         int Jet_puId_){
        Jet_pt.push_back(Jet_pt_); 
        Jet_phi.push_back(Jet_phi_);
        Jet_eta.push_back(Jet_eta_);
        Jet_mass.push_back(Jet_mass_);
        Jet_btagDeepC.push_back(Jet_btagDeepC_);
        Jet_jetId.push_back(Jet_jetId_);
        Jet_puId.push_back(Jet_puId_); //1 or 0?
        Jet_btagDeepFlavB.push_back(Jet_btagDeepFlavB_);
        Jet_btagPNetB.push_back(Jet_btagPNetB_);
        // NanoAODv12 / 2022 / 2023
        Jet_btagRobustParTAK4B.push_back(Jet_btagRobustParTAK4B_);

        // NanoAODv15 / 2024
        Jet_btagUParTAK4B.push_back(Jet_btagUParTAK4B_);
      }
    
      
      void SetMuons(float Muon_pt_, float Muon_eta_, float Muon_phi_, float Muon_mass_, bool Muon_isGlobal_, bool Muon_isTracker_,
                        float Muon_dxy_, float Muon_dz_,float Muon_sip3d_, float Muon_ptErr_,
                        bool Muon_looseId_, bool Muon_mediumId_, bool Muon_tightId_, 
                        unsigned char Muon_mvaLowPtId_, unsigned char Muon_mvaId_, float Muon_mvaLowPt_,
                        int Muon_nTrackerLayers_, bool Muon_isPFcand_, int Muon_pdgId_,int Muon_charge_, float Muon_pfRelIso03_all_
                        ){
        Muon_pt.push_back(Muon_pt_); 
        Muon_phi.push_back(Muon_phi_);
        Muon_eta.push_back(Muon_eta_);
        Muon_mass.push_back(Muon_mass_);
        Muon_isGlobal.push_back(Muon_isGlobal_);
        Muon_isTracker.push_back(Muon_isTracker_);
        Muon_dxy.push_back(Muon_dxy_);
        Muon_dz.push_back(Muon_dz_);
        Muon_sip3d.push_back(Muon_sip3d_);
        Muon_ptErr.push_back(Muon_ptErr_);
        Muon_nTrackerLayers.push_back(Muon_nTrackerLayers_);
        Muon_isPFcand.push_back(Muon_isPFcand_);
        Muon_pdgId.push_back(Muon_pdgId_);
        Muon_charge.push_back(Muon_charge_);
        Muon_pfRelIso03_all.push_back(Muon_pfRelIso03_all_);
        Muon_looseId.push_back(Muon_looseId_);
        Muon_mediumId.push_back(Muon_mediumId_);
        Muon_tightId.push_back(Muon_tightId_);
        Muon_mvaLowPtId.push_back(Muon_mvaLowPtId_);
        Muon_mvaId.push_back(Muon_mvaId_);
        Muon_mvaLowPt.push_back(Muon_mvaLowPt_);
      }
      void SetMuonsGen(int Muon_genPartIdx_){
        Muon_genPartIdx.push_back(Muon_genPartIdx_);
      }
      void SetElectronsGen(int Electron_genPartIdx_){
        Electron_genPartIdx.push_back(Electron_genPartIdx_);
      }
      /*void SetMuons(TTreeReaderArray<float> *Muon_pt_, TTreeReaderArray<float> *Muon_eta_,
                        TTreeReaderArray<float> *Muon_phi_, TTreeReaderArray<float> *Muon_mass_, TTreeReaderArray<bool> *Muon_isGlobal_, TTreeReaderArray<bool> *Muon_isTracker_,
                        TTreeReaderArray<float> *Muon_dxy_, TTreeReaderArray<float> *Muon_dz_,TTreeReaderArray<float> *Muon_sip3d_, TTreeReaderArray<float> *Muon_ptErr_,
                        TTreeReaderArray<int> *Muon_nTrackerLayers_, TTreeReaderArray<bool> *Muon_isPFcand_, TTreeReaderArray<int> *Muon_pdgId_,TTreeReaderArray<int> *Muon_charge_, TTreeReaderArray<float> *Muon_pfRelIso03_all_,
                        TTreeReaderArray<int> *Muon_genPartIdx_){
        Muon_pt = Muon_pt_; 
        Muon_phi = Muon_phi_;
        Muon_eta = Muon_eta_;
        Muon_mass = Muon_mass_;
        Muon_isGlobal = Muon_isGlobal_;
        Muon_isTracker = Muon_isTracker_;
        Muon_dxy = Muon_dxy_;
        Muon_dz = Muon_dz_;
        Muon_sip3d = Muon_sip3d_;
        Muon_ptErr = Muon_ptErr_;
        Muon_nTrackerLayers = Muon_nTrackerLayers_;
        Muon_isPFcand = Muon_isPFcand_;
        Muon_pdgId = Muon_pdgId_;
        Muon_charge = Muon_charge_;
        Muon_pfRelIso03_all = Muon_pfRelIso03_all_;
        Muon_genPartIdx = Muon_genPartIdx_;
      }*/
      void SetFsrPhotons(float FsrPhoton_dROverEt2_, float FsrPhoton_eta_,
                        float FsrPhoton_phi_, float FsrPhoton_pt_, float FsrPhoton_relIso03_, int FsrPhoton_electronIdx_, int FsrPhoton_muonIdx_){
        FsrPhoton_dROverEt2.push_back(FsrPhoton_dROverEt2_); 
        FsrPhoton_phi.push_back(FsrPhoton_phi_);
        FsrPhoton_eta.push_back(FsrPhoton_eta_);
        FsrPhoton_pt.push_back(FsrPhoton_pt_);
        FsrPhoton_relIso03.push_back(FsrPhoton_relIso03_);
        FsrPhoton_electronIdx.push_back(FsrPhoton_electronIdx_);
        FsrPhoton_muonIdx.push_back(FsrPhoton_muonIdx_);
      }
      /*void SetFsrPhotons(TTreeReaderArray<float> *FsrPhoton_dROverEt2_, TTreeReaderArray<float> *FsrPhoton_eta_,
                        TTreeReaderArray<float> *FsrPhoton_phi_, TTreeReaderArray<float> *FsrPhoton_pt_, 
                        TTreeReaderArray<float> *FsrPhoton_relIso03_){
        FsrPhoton_dROverEt2 = FsrPhoton_dROverEt2_; 
        FsrPhoton_phi = FsrPhoton_phi_;
        FsrPhoton_eta = FsrPhoton_eta_;
        FsrPhoton_pt = FsrPhoton_pt_;
        FsrPhoton_relIso03 = FsrPhoton_relIso03_; 
      }*/
      void SetGenParts(float GenPart_pt_){
        GenPart_pt.push_back(GenPart_pt_);
      }
      /*void SetGenParts(TTreeReaderArray<float> *GenPart_pt_){
        GenPart_pt = GenPart_pt_;
      }*/
      void SetObjectNum(unsigned nElectron_,unsigned nMuon_,unsigned nJet_,unsigned nFsrPhoton_){
        nElectron = nElectron_; 
        nMuon = nMuon_;
        nJet = nJet_;
        nFsrPhoton = nFsrPhoton_;
      }
      void SetObjectNumGen(unsigned nGenPart_){
        nGenPart = nGenPart_;
      }
      bool isMC;
      std::string analysisMode;
      void SetAnalysisMode(const std::string& mode);
      void SetNanoVersion(int v){ nanoVersion = v; }
      bool BuildZZCandidate();
      bool BuildBestDijet();
      std::vector<unsigned int> goodLooseElectrons2012();
      std::vector<unsigned int> goodLooseMuons2012();
      std::vector<unsigned int> goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex);
      std::vector<unsigned int> goodElectrons2015_noIso_noBdt(std::vector<unsigned int> Electronindex);
      std::vector<bool> pass_Ele_Id(int nanoVersion);
      std::vector<bool> pass_Mu_Id(const std::string& era, const std::string& method, const std::string& wp);
      std::vector<unsigned int> goodFsrPhotons();
      unsigned doFsrRecovery(TLorentzVector Lep);
      unsigned doFsrRecovery_Run3(std::vector<unsigned int> goodfsridx, unsigned lepidx, int lepflavor);//lepflavor 11 or 13
      std::vector<TLorentzVector> BatchFsrRecovery(std::vector<TLorentzVector> LepList);
      void BatchFsrRecovery_Run3();
      std::vector<TLorentzVector> ElectronFsr();
      std::vector<TLorentzVector> MuonFsr();
      std::vector<float> ElectronFsrPt();
      std::vector<float> ElectronFsrEta();
      std::vector<float> ElectronFsrPhi();
      std::vector<float> MuonFsrPt();
      std::vector<float> MuonFsrEta();
      std::vector<float> MuonFsrPhi();
      std::vector<unsigned int> SelectedJets(std::vector<unsigned int> ele, std::vector<unsigned int> mu);
      std::vector<TLorentzVector> Electrondressed_Run3;
      std::vector<TLorentzVector> Muondressed_Run3;
      std::vector<TLorentzVector> Zlist;
      std::vector<TLorentzVector> Zlistnofsr;
      std::vector<int> Zflavor; //mu->13, e->11
      std::vector<int> Zlep1index;
      std::vector<int> Zlep2index;
      std::vector<int> Zlep1lepindex;
      std::vector<int> Zlep2lepindex;
      std::vector<float> Zlep1pt;
      std::vector<float> Zlep1eta;
      std::vector<float> Zlep1phi;
      std::vector<float> Zlep1mass;
      std::vector<float> Zlep1chg;
      std::vector<float> Zlep2pt;
      std::vector<float> Zlep2eta;
      std::vector<float> Zlep2phi;
      std::vector<float> Zlep2mass;
      std::vector<float> Zlep2chg;
      std::vector<float> Zlep1ptNoFsr;
      std::vector<float> Zlep1etaNoFsr;
      std::vector<float> Zlep1phiNoFsr;
      std::vector<float> Zlep1massNoFsr;
      std::vector<float> Zlep2ptNoFsr;
      std::vector<float> Zlep2etaNoFsr;
      std::vector<float> Zlep2phiNoFsr;
      std::vector<float> Zlep2massNoFsr;
      std::vector<unsigned int> jetidx;
      int nTightEle;
      int nTightMu;
      int nTightEleChgSum;
      int nTightMuChgSum;
      int njets_pt30_eta4p7;
      int Lepointer;
    
      bool flag4e;
      bool flag4mu;
      bool flag2e2mu;

      void LeptonSelection();
      std::vector<unsigned int> step1Ele,step1Mu,bestEle,bestMu, tighteleforjetidx, tightmuforjetidx;
      std::vector<unsigned int> Electronindex;
      std::vector<unsigned int> Muonindex;
      std::vector<bool> AllEid;
      std::vector<bool> AllMuid;
      std::vector<TLorentzVector> Elelist;
      std::vector<TLorentzVector> Mulist;
      std::vector<TLorentzVector> ElelistFsr;
      std::vector<TLorentzVector> MulistFsr;
      std::vector<int> Elechg;
      std::vector<int> Muchg;
      std::vector<float> Muiso,Eiso;
      std::vector<bool> Eid;
      std::vector<bool> muid;
      std::vector<int> lep_genindex;
      std::vector<int> TightElelep_index;
      std::vector<int> TightMulep_index;
      int lep_Hindex[4];
      TLorentzVector Z1;
      TLorentzVector Z1nofsr;
      TLorentzVector Z2;
      TLorentzVector Z2nofsr;
      TLorentzVector ZZsystem;
      TLorentzVector ZZsystemnofsr;

      std::vector<int> TightEleindex;
      std::vector<int> TightMuindex;
      void Initialize(){
        step1Ele.clear();
        step1Mu.clear();
        bestEle.clear();
        bestMu.clear();
        tighteleforjetidx.clear();
        tightmuforjetidx.clear();
        Electronindex.clear();
        Muonindex.clear();
        AllEid.clear();
        AllMuid.clear();
        Elelist.clear();
        Mulist.clear();
        ElelistFsr.clear();
        MulistFsr.clear();
        Electron_pt.clear();Electron_phi.clear();Electron_eta.clear();Electron_mass.clear();Electron_dxy.clear();Electron_dz.clear();Electron_sip3d.clear();
        Electron_deltaEtaSC.clear();Electron_mvaIso_WPHZZ.clear();
        //Electron_mvaHZZIso.clear();
        Electron_mvaNoIso.clear();
        Electron_pdgId.clear();Electron_genPartIdx.clear();Electron_pfRelIso03_all.clear();
        Muon_pt.clear();Muon_phi.clear();Muon_eta.clear();Muon_mass.clear();Muon_dxy.clear();Muon_dz.clear();Muon_sip3d.clear();Muon_ptErr.clear();Muon_pfRelIso03_all.clear();
        Muon_nTrackerLayers.clear();Muon_genPartIdx.clear();Muon_pdgId.clear();Muon_charge.clear();
        Muon_isTracker.clear();Muon_isGlobal.clear();Muon_isPFcand.clear();
        Muon_looseId.clear();Muon_mediumId.clear();Muon_tightId.clear();Muon_mvaLowPtId.clear();Muon_mvaId.clear();Muon_mvaLowPt.clear();
        Jet_pt.clear();Jet_phi.clear();Jet_eta.clear();Jet_mass.clear();Jet_btagDeepC.clear();Jet_btagDeepFlavB.clear();Jet_btagPNetB.clear();Jet_btagRobustParTAK4B.clear();Jet_btagUParTAK4B.clear();
        Jet_jetId.clear();Jet_puId.clear(); Zlep1lepindex.clear();Zlep2lepindex.clear(); 
        FsrPhoton_dROverEt2.clear();FsrPhoton_phi.clear();FsrPhoton_eta.clear();FsrPhoton_pt.clear();FsrPhoton_relIso03.clear(); FsrPhoton_electronIdx.clear(); FsrPhoton_muonIdx.clear();
        GenPart_pt.clear();
        Zlist.clear();
        Zlistnofsr.clear();
        Zflavor.clear();
        Electrondressed_Run3.clear();
        Muondressed_Run3.clear();
        Zlep1index.clear();
        Zlep2index.clear();
        Zlep1pt.clear(); Zlep1eta.clear(); Zlep1phi.clear(); Zlep1mass.clear();
        Zlep2pt.clear(); Zlep2eta.clear(); Zlep2phi.clear(); Zlep2mass.clear();
        Zlep1chg.clear(); Zlep2chg.clear();
        Zlep1ptNoFsr.clear(); Zlep1etaNoFsr.clear(); Zlep1phiNoFsr.clear(); Zlep1massNoFsr.clear();
        Zlep2ptNoFsr.clear(); Zlep2etaNoFsr.clear(); Zlep2phiNoFsr.clear(); Zlep2massNoFsr.clear();
        jetidx.clear(); lep_genindex.clear(); TightElelep_index.clear();TightMulep_index.clear();
        step1Ele.clear(); step1Mu.clear(); bestEle.clear(); bestMu.clear();  tighteleforjetidx.clear();  tightmuforjetidx.clear(); 
        Electronindex.clear();  Muonindex.clear(); AllEid.clear(); AllMuid.clear(); Elelist.clear(); Mulist.clear(); ElelistFsr.clear(); Mulist.clear(); 
        Elechg.clear(); Muchg.clear(); Muiso.clear();Eiso.clear(); Eid.clear(); muid.clear(); TightEleindex.clear(); TightMuindex.clear();
        for (int i=0; i<4; i++) {lep_Hindex[i]=-1;}
        Z1.SetPtEtaPhiM(0,0,0,0);
        Z1nofsr.SetPtEtaPhiM(0,0,0,0);
        Z2.SetPtEtaPhiM(0,0,0,0);
        Z2nofsr.SetPtEtaPhiM(0,0,0,0);
        ZZsystem.SetPtEtaPhiM(0,0,0,0);
        ZZsystemnofsr.SetPtEtaPhiM(0,0,0,0);
        nElectron = 0; nMuon = 0; nJet = 0; nFsrPhoton = 0; nGenPart = 0; nGenJet = 0;
        nTightEle = 0; nTightMu = 0; nTightEleChgSum = 0; nTightMuChgSum = 0; Z1flav = 0; Zsize = 0;
        Lepointer = 0; 
        
        nRawJetsThisEvent = 0;
        nPtEtaJetsThisEvent = 0;
        nJetIdJetsThisEvent = 0;
        nPuIdJetsThisEvent = 0;

        pTL1 = -999; etaL1 = -999; phiL1 = -999; massL1 = -999;
        pTL2 = -999; etaL2 = -999; phiL2 = -999; massL2 = -999;
        pTL3 = -999; etaL3 = -999; phiL3 = -999; massL3 = -999;
        pTL4 = -999; etaL4 = -999; phiL4 = -999; massL4 = -999;

        pTj1 = -99;  etaj1 = -99;  phij1 = -99;  mj1 = -99; btagger1_DJ = -99; btagger1_PN = -99; btagger1_RPT = -99; btagger1_UPT = -99;
        pTj2 = -99;  etaj2 = -99;  phij2 = -99;  mj2 = -99; btagger2_DJ = -99; btagger2_PN = -99; btagger2_RPT = -99; btagger2_UPT = -99;
        invjj = -99;
        njets_pt30_eta4p7 = 0;
        RecoFourMuEvent=false; RecoFourEEvent=false; RecoTwoETwoMuEvent=false; RecoTwoMuTwoEEvent=false;
        flag4e=false; flag4mu=false; flag2e2mu=false;

        eventPassTwoTightLeps = false;
        eventPassZCand = false;
        eventPassAtLeastTwoRawJets = false;
        eventPassAtLeastTwoPtEtaJets = false;
        eventPassAtLeastTwoJetIdJets = false;
        eventPassAtLeastTwoPuIdJets = false;
        eventPassTwoGoodJets = false;
        eventPassDijet = false;
        eventPassFinal = false;

      }
      bool isFSR=true;
      unsigned int Zsize=0;
      TSpline *spline_g4;
      TSpline *spline_g2;
      TSpline *spline_L1;
      TSpline *spline_L1Zgs;
      bool findZCandidate();
      bool ZZSelection();

      Mela* mela;
      float me_0plus_JHU, me_qqZZ_MCFM, p0plus_m4l, bkg_m4l;
      float D_bkg_kin, D_bkg, D_g4, D_g1g4, D_CP, D_int, D_L1_int, D_L1Zgint;
      float D_0m, D_0hp, D_L1, D_L1Zg;
      float D_bkg_kin_vtx_BS;
      float p0minus_VAJHU, Dgg10_VAMCFM, pg1g4_VAJHU;
      float p0plus_VAJHU, p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, pDL1_VAJHU, pD_L1Zgint; //, p0plus_VAJHU;
      float getDg4Constant(float ZZMass);
      float getDg2Constant(float ZZMass);
      float getDL1Constant(float ZZMass);
      float getDL1ZgsConstant(float ZZMass); 

      float pTj1, etaj1, phij1, mj1, pTj2, etaj2, phij2, mj2;
      float pTL1, etaL1, phiL1, massL1;
      float pTL2, etaL2, phiL2, massL2;
      float pTL3, etaL3, phiL3, massL3;
      float pTL4, etaL4, phiL4, massL4;


      
    private:
      std::vector<float> Electron_pt,Electron_phi,Electron_eta,Electron_mass,Electron_dxy,Electron_dz,Electron_sip3d;
      std::vector<float> Electron_pfRelIso03_all;
      std::vector<float> Electron_deltaEtaSC;
      //std::vector<float> Electron_mvaHZZIso;
      std::vector<float> Electron_mvaNoIso;
      std::vector<bool> Electron_mvaIso_WPHZZ;
      std::vector<int> Electron_pdgId,Electron_genPartIdx;

      std::vector<float> Jet_pt,Jet_phi,Jet_eta,Jet_mass,Jet_btagDeepC,Jet_btagDeepFlavB,Jet_btagPNetB,Jet_btagRobustParTAK4B,Jet_btagUParTAK4B;
      std::vector<int> Jet_jetId,Jet_puId;

      std::vector<float> Muon_pt,Muon_phi,Muon_eta,Muon_mass,Muon_dxy,Muon_dz,Muon_sip3d,Muon_ptErr,Muon_pfRelIso03_all;
      std::vector<int> Muon_nTrackerLayers,Muon_genPartIdx,Muon_pdgId,Muon_charge;
      std::vector<bool> Muon_isTracker,Muon_isGlobal,Muon_isPFcand;
      std::vector<bool> Muon_looseId,Muon_mediumId,Muon_tightId;
      std::vector<unsigned char> Muon_mvaLowPtId,Muon_mvaId;
      std::vector<float> Muon_mvaLowPt;

      std::vector<float> FsrPhoton_dROverEt2, FsrPhoton_phi, FsrPhoton_pt, FsrPhoton_relIso03, FsrPhoton_eta, FsrPhoton_muonIdx, FsrPhoton_electronIdx;

      std::vector<float> GenPart_pt;
      
      unsigned nElectron,nMuon,nJet,nGenPart,nFsrPhoton,nGenJet;
      int nanoVersion = 12;


};

H4LTools::H4LTools(int year, bool isMC_){
  isMC = isMC_;
  std::cout<<"year"<<" "<<year<<std::endl;
  mela = new Mela(13.0, 125.0, TVar::SILENT);
  mela->setCandidateDecayMode(TVar::CandidateDecay_ZZ);  
  TFile *gConstant_g4 = TFile::Open("CoupleConstantsForMELA/gConstant_HZZ2e2mu_g4.root");
  spline_g4 = (TSpline*) gConstant_g4->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g4");
  gConstant_g4->Close();
  delete gConstant_g4;
  TFile *gConstant_g2 = TFile::Open("CoupleConstantsForMELA/gConstant_HZZ2e2mu_g2.root");
  spline_g2 = (TSpline*) gConstant_g2->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_g2");
  gConstant_g2->Close();
  delete gConstant_g2;
  TFile *gConstant_L1 = TFile::Open("CoupleConstantsForMELA/gConstant_HZZ2e2mu_L1.root");
  spline_L1 = (TSpline*) gConstant_L1->Get("sp_tgfinal_HZZ2e2mu_SM_over_tgfinal_HZZ2e2mu_L1");
  gConstant_L1->Close();
  delete gConstant_L1;
  TFile *gConstant_L1Zgs = TFile::Open("CoupleConstantsForMELA/gConstant_HZZ2e2mu_L1Zgs.root");
  spline_L1Zgs = (TSpline*) gConstant_L1Zgs->Get("sp_tgfinal_HZZ2e2mu_SM_photoncut_over_tgfinal_HZZ2e2mu_L1Zgs");
  gConstant_L1Zgs->Close();
  delete gConstant_L1Zgs;

  cut2e2mu = 0;
  cut4e = 0;
  cut4mu = 0;
  cutghost2e2mu = 0;
  cutghost4e = 0;
  cutghost4mu = 0;
  cutLepPt2e2mu = 0;
  cutLepPt4e = 0;
  cutLepPt4mu = 0;
  cutQCD2e2mu = 0;
  cutQCD4e = 0;
  cutQCD4mu = 0;
  cutZZ2e2mu = 0;
  cutZZ4e = 0;
  cutZZ4mu = 0;
  cutm4l2e2mu = 0;
  cutm4l4e = 0;
  cutm4l4mu = 0;

  passTwoTightLeps = 0;
  passZCand = 0;
  passTwoGoodJets = 0;
  passDijet = 0;
  passFinal = 0;

  passAtLeastTwoRawJets = 0;
  passAtLeastTwoPtEtaJets = 0;
  passAtLeastTwoJetIdJets = 0;
  passAtLeastTwoPuIdJets = 0;
  
  nRawJetsThisEvent = 0;
  nPtEtaJetsThisEvent = 0;
  nJetIdJetsThisEvent = 0;
  nPuIdJetsThisEvent = 0;
}
#endif

