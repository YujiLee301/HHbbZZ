#include "../interface/H4LTools.h"
#include "../interface/H4LSelection.h"
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <vector>
#include <iostream>

void H4LTools::SetAnalysisMode(const std::string& mode){
    analysisMode = mode;
}

std::vector<unsigned int> H4LTools::goodLooseElectrons2012(){
    std::vector<unsigned int> LooseElectronindex;
    for (unsigned int i=0; i<Electron_pt.size(); i++){
        if ((Electron_pt[i]>elePtcut)&&(fabs(Electron_eta[i])<eleEtacut)){
            LooseElectronindex.push_back(i);
        }
    }

    return LooseElectronindex;
}

std::vector<unsigned int> H4LTools::goodLooseMuons2012(){
    std::vector<unsigned int> LooseMuonindex;
    for (unsigned int i=0; i<Muon_eta.size(); i++){
        if ((Muon_pt[i]>MuPtcut)&&(fabs(Muon_eta[i])<MuEtacut)&&((Muon_isGlobal[i]||Muon_isTracker[i]||Muon_isPFcand[i]))){
            LooseMuonindex.push_back(i);
        }
    }

    return LooseMuonindex;
}

std::vector<unsigned int> H4LTools::goodMuons2015_noIso_noPf(std::vector<unsigned int> Muonindex){
    std::vector<unsigned int> bestMuonindex;
    for (unsigned int i=0; i<Muonindex.size(); i++){
        if ((Muon_pt[Muonindex[i]]>MuPtcut)&&(fabs(Muon_eta[Muonindex[i]])<MuEtacut)&&(Muon_isGlobal[Muonindex[i]]||Muon_isTracker[Muonindex[i]])){
            if (Muon_sip3d[Muonindex[i]]<Musip3dCut){
                if((fabs(Muon_dxy[Muonindex[i]])<MuLoosedxycut)&&(fabs(Muon_dz[Muonindex[i]])<MuLoosedzcut)){
                    bestMuonindex.push_back(Muonindex[i]);
                }
            }
        }
    }

    return bestMuonindex;
}

std::vector<unsigned int> H4LTools::goodElectrons2015_noIso_noBdt(std::vector<unsigned int> Electronindex){
    std::vector<unsigned int> bestElectronindex;
    for (unsigned int i=0; i<Electronindex.size(); i++){
        if ((Electron_pt[Electronindex[i]])>elePtcut){
            if(Electron_sip3d[Electronindex[i]]<elesip3dCut){
                if((fabs(Electron_dxy[Electronindex[i]])<eleLoosedxycut)&&(fabs(Electron_dz[Electronindex[i]])<eleLoosedzcut)){
                    bestElectronindex.push_back(Electronindex[i]);
                }
            }
        }
    }

    return bestElectronindex;
}

std::vector<bool> H4LTools::pass_Ele_Id(int nanoVersion){
    std::vector<bool> passid;
    float cutVal = 1000.;
    float mvaVal = -999.;

    for (unsigned int i = 0; i < Electron_pt.size(); i++){

        // Use SC eta as in ZZAnalysis-style electron BDT definition
        float fSCeta = fabs(Electron_eta[i] + Electron_deltaEtaSC[i]);

        //if (nanoVersion < 14){
            // Run3 NanoAODv12/v13 style:
            // Electron_mvaHZZIso + hand cuts
            
            if(Electron_pt[i]<10){
                if(fSCeta<0.8) cutVal = eleBDTWPLELP;
                if((fSCeta>=0.8)&&(fSCeta<1.479)) cutVal = eleBDTWPMELP;
                if(fSCeta>=1.479) cutVal = eleBDTWPHELP;
            }
            else{
                if(fSCeta<0.8) cutVal = eleBDTWPLEHP;
                if((fSCeta>=0.8)&&(fSCeta<1.479)) cutVal = eleBDTWPMEHP;
                if(fSCeta>=1.479) cutVal = eleBDTWPHEHP;
            }
            
            //mvaVal = Electron_mvaHZZIso[i];
            mvaVal = Electron_mvaNoIso[i];
            passid.push_back(mvaVal > cutVal);
            //passid.push_back(Electron_mvaIso_WP90[i]);
            //passid.push_back(Electron_mvaNoIso_WP90[i]);
        //}
        //else{
            // Run3 NanoAODv14+ style:
            // use predefined HZZ working point
            //passid.push_back(Electron_mvaIso_WPHZZ[i]);
            //passid.push_back(Electron_mvaIso_WP90[i]);
            //passid.push_back(Electron_mvaNoIso_WP90[i]);
        //}
    }

    return passid;
}


std::vector<bool> H4LTools::pass_Mu_Id(const std::string& era, const std::string& method, const std::string& wp) {
    std::vector<bool> passId;

    for (unsigned int i = 0; i < Muon_pt.size(); i++) {

        if (method == "PF Cut-Based") {
            // PF Cut-Based selection (same for Run2 and Run3)
            if (wp == "loose")      passId.push_back(Muon_looseId[i]);
            else if (wp == "medium") passId.push_back(Muon_mediumId[i]);
            else if (wp == "tight")  passId.push_back(Muon_tightId[i]);
            else {
                std::cerr << "[MuonID] Unknown WP for PF Cut-Based: " << wp << std::endl;
                passId.push_back(false);
            }
        }
        else if (method == "MVA-Based") {
            if (era == "Run2") {
                // Run2 NanoAOD MVA ID
                if (Muon_pt[i] < 3.0) {
                    // soft muon, use Muon_mvaLowPtId
                    if      (wp == "loose")  passId.push_back(Muon_mvaLowPtId[i] >= 1);
                    else if (wp == "medium") passId.push_back(Muon_mvaLowPtId[i] >= 2);
                    else if (wp == "tight")  passId.push_back(Muon_mvaLowPtId[i] >= 3);
                    else {
                        std::cerr << "[MuonID] Unknown WP for Run2 low-pT MVA: " << wp << std::endl;
                        passId.push_back(false);
                    }
                } else {
                    // normal muon, use Muon_mvaId
                    if      (wp == "loose")  passId.push_back(Muon_mvaId[i] >= 1);
                    else if (wp == "medium") passId.push_back(Muon_mvaId[i] >= 2);
                    else if (wp == "tight")  passId.push_back(Muon_mvaId[i] >= 3);
                    else if (wp == "vtight") passId.push_back(Muon_mvaId[i] >= 4);
                    else if (wp == "vvtight") passId.push_back(Muon_mvaId[i] >= 5);
                    else {
                        std::cerr << "[MuonID] Unknown WP for Run2 MVA-Based: " << wp << std::endl;
                        passId.push_back(false);
                    }
                }
            }

            else if (era == "Run3") {
            // Cut-based + low-pT MVA + loose SIP cut

                bool passSIP = fabs(Muon_sip3d[i]) < 8; // 3D impact parameter significance wrt first PV

                bool passCutBased = false;
                if      (wp == "loose")  passCutBased = Muon_looseId[i];
                else if (wp == "medium") passCutBased = Muon_mediumId[i];
                else if (wp == "tight")  passCutBased = Muon_tightId[i];
                else {
                    std::cerr << "[MuonID] Unknown WP for Run3 Cut-Based: " << wp << std::endl;
                    passId.push_back(false);
                    continue;
                }

                bool passLowPtMVA = (Muon_mvaLowPt[i] > -0.6);

                passId.push_back(passCutBased && passLowPtMVA && passSIP);
            }
        }
    }

    return passId;
}

std::vector<unsigned int> H4LTools::goodFsrPhotons(){
    std::vector<unsigned int> goodFsrPhoton;
    //unsigned nFsr = (*nFsrPhoton).Get()[0];
    for (unsigned int i=0; i<FsrPhoton_pt.size(); i++){
        if((FsrPhoton_pt[i]>fsrphotonPtcut)&&(fabs(FsrPhoton_eta[i])<fsrphotonEtacut)&&(FsrPhoton_relIso03[i]<fsrphotonIsocut)){
            goodFsrPhoton.push_back(i);
        }
    }
    return goodFsrPhoton;
}

std::vector<unsigned int> H4LTools::SelectedJets(std::vector<unsigned int> ele, std::vector<unsigned int> mu){
    std::vector<unsigned int> goodJets;

    nRawJetsThisEvent = Jet_pt.size();
    nPtEtaJetsThisEvent = 0;
    nJetIdJetsThisEvent = 0;
    nPuIdJetsThisEvent = 0;

    for(unsigned int i=0; i<Jet_pt.size(); i++){

        // step 1: pt/eta
        if((Jet_pt[i] > JetPtcut) && (fabs(Jet_eta[i]) < JetEtacut)){
            nPtEtaJetsThisEvent++;

            // step 2: jetId
            if(Jet_jetId[i] > 0){
                nJetIdJetsThisEvent++;

                // step 3: lepton cleaning
                int overlaptag = 0;
                TLorentzVector jettest;
                jettest.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

                for(unsigned int ie=0; ie<ele.size(); ie++){
                    TLorentzVector eletest;
                    eletest.SetPtEtaPhiM(Electron_pt[ele[ie]], Electron_eta[ele[ie]], Electron_phi[ele[ie]], Electron_mass[ele[ie]]);
                    if(eletest.DeltaR(jettest) < 0.4) overlaptag++;
                }

                for(unsigned int im=0; im<mu.size(); im++){
                    TLorentzVector mutest;
                    mutest.SetPtEtaPhiM(Muon_pt[mu[im]], Muon_eta[mu[im]], Muon_phi[mu[im]], Muon_mass[mu[im]]);
                    if(mutest.DeltaR(jettest) < 0.4) overlaptag++;
                }

                if(overlaptag == 0) goodJets.push_back(i);
                
            }
        }
    }

    njets_pt30_eta4p7 = goodJets.size();
    return goodJets;
}

unsigned H4LTools::doFsrRecovery(TLorentzVector Lep){
    // This Function returns the index for the possible FsrPhoton
    unsigned int FsrIdx = 999; //only Idx>0 works, pay attention!
    std::vector<unsigned int> BestFsrPhotons;
    BestFsrPhotons = goodFsrPhotons();
    float dRl,dRlOverPt;
    dRl = 999;
    dRlOverPt = 999;
    for(unsigned int i=0;i<BestFsrPhotons.size();i++){
        TLorentzVector fsrcand;
        fsrcand.SetPtEtaPhiM(FsrPhoton_pt[BestFsrPhotons[i]],FsrPhoton_eta[BestFsrPhotons[i]],FsrPhoton_phi[BestFsrPhotons[i]],0);
        float dRlC,dRlOverPtC;
        dRlC = fsrcand.DeltaR(Lep);
        if ((dRlC<fsrphotondRlcut)&&(FsrPhoton_dROverEt2[BestFsrPhotons[i]]<fsrphotondRlOverPtcut)){
            if(FsrPhoton_dROverEt2[BestFsrPhotons[i]]<dRlOverPt){
                dRl = dRlC;
                dRlOverPt = FsrPhoton_dROverEt2[BestFsrPhotons[i]];
                FsrIdx = BestFsrPhotons[i];
            }
        }
    }

    return FsrIdx;

}

unsigned H4LTools::doFsrRecovery_Run3(std::vector<unsigned int> goodfsridx, unsigned lepidx, int lepflavor){//lepflavor 11 or 13

    unsigned matchedfsridx = 999;
    if(lepflavor == 11){
        for(unsigned fsridx=0; fsridx<goodfsridx.size(); fsridx++){
            if(FsrPhoton_electronIdx[goodfsridx[fsridx]] == lepidx){
                matchedfsridx = goodfsridx[fsridx];
                break;
            }
        }
    }
    if(lepflavor == 13){
        for(unsigned fsridx=0; fsridx<goodfsridx.size(); fsridx++){
            if(FsrPhoton_muonIdx[goodfsridx[fsridx]] == lepidx){
                matchedfsridx = goodfsridx[fsridx];
                break;
            }
        }
    }
    return matchedfsridx;
}

void H4LTools::BatchFsrRecovery_Run3(){
    unsigned fsridx;
    std::vector<unsigned> fsrlist;
    fsrlist = goodFsrPhotons();
    for(unsigned int i=0; i<Electron_pt.size(); i++){
        TLorentzVector fsr,lep;
        lep.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
        fsridx = doFsrRecovery_Run3(fsrlist,i,11);
        if(fsridx<900){
            fsr.SetPtEtaPhiM(
                FsrPhoton_pt[fsridx], 
                FsrPhoton_eta[fsridx], 
                FsrPhoton_phi[fsridx], 
                0
            );
            lep = lep + fsr;
            Electrondressed_Run3.push_back(lep);
        }
        else{
            Electrondressed_Run3.push_back(lep);
        }
    }
    for(unsigned int j=0; j<Muon_pt.size(); j++){
        TLorentzVector fsr,lep;
        lep.SetPtEtaPhiM(Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_mass[j]);
        fsridx = doFsrRecovery_Run3(fsrlist,j,13);
        if(fsridx<900){
            fsr.SetPtEtaPhiM(
                FsrPhoton_pt[fsridx], 
                FsrPhoton_eta[fsridx], 
                FsrPhoton_phi[fsridx], 
                0
            );
            lep = lep + fsr;
            Muondressed_Run3.push_back(lep);
        }
        else{
            Muondressed_Run3.push_back(lep);
        }
    }
}

std::vector<TLorentzVector> H4LTools::BatchFsrRecovery(std::vector<TLorentzVector> LepList){
    std::vector<TLorentzVector> LepFsrList;

    for(unsigned int i=0;i<LepList.size();i++){
        int Fsrtag;
        Fsrtag = doFsrRecovery(LepList[i]);
        if (Fsrtag>900){
            LepFsrList.push_back(LepList[i]);
            continue;
        }
        TLorentzVector FsrPhoton;
        FsrPhoton.SetPtEtaPhiM(FsrPhoton_pt[Fsrtag],FsrPhoton_eta[Fsrtag],FsrPhoton_phi[Fsrtag],0);
        TLorentzVector LepFsrRecovery;
        LepFsrRecovery = FsrPhoton + LepList[i];
        LepFsrList.push_back(LepFsrRecovery);
    }
    return LepFsrList;
}

std::vector<TLorentzVector> H4LTools::ElectronFsr(){
    std::vector<TLorentzVector> leplist;
    std::vector<TLorentzVector> leplistfsr;
    //unsigned nlep = (*nElectron).Get()[0];
    for(unsigned int i=0;i<Electron_pt.size();i++){
        TLorentzVector Lep;
        Lep.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
        leplist.push_back(Lep);
    }
    leplistfsr = BatchFsrRecovery(leplist);
    return leplistfsr;
}

std::vector<TLorentzVector> H4LTools::MuonFsr(){
    std::vector<TLorentzVector> leplist;
    std::vector<TLorentzVector> leplistfsr;
    //unsigned nlep = (*nMuon).Get()[0];
    for(unsigned int i=0;i<Muon_pt.size();i++){
        TLorentzVector Lep;
        Lep.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],Muon_mass[i]);
        leplist.push_back(Lep);
    }
    leplistfsr = BatchFsrRecovery(leplist);
    return leplistfsr;
}

std::vector<float> H4LTools::ElectronFsrPt(){
    std::vector<float> lepPt;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepPt.push_back(Electrondressed_Run3[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::ElectronFsrEta(){
    std::vector<float> lepEta;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepEta.push_back(Electrondressed_Run3[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::ElectronFsrPhi(){
    std::vector<float> lepPhi;
    for (unsigned int i=0;i<Electrondressed_Run3.size();i++){
        lepPhi.push_back(Electrondressed_Run3[i].Phi());
    }
    return lepPhi;
}

std::vector<float> H4LTools::MuonFsrPt(){
    std::vector<float> lepPt;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepPt.push_back(Muondressed_Run3[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::MuonFsrEta(){
    std::vector<float> lepEta;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepEta.push_back(Muondressed_Run3[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::MuonFsrPhi(){
    std::vector<float> lepPhi;
    for (unsigned int i=0;i<Muondressed_Run3.size();i++){
        lepPhi.push_back(Muondressed_Run3[i].Phi());
    }
    return lepPhi;
}

/*std::vector<float> H4LTools::ElectronFsrPt(){
    std::vector<float> lepPt;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPt.push_back(leplistfsr[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::ElectronFsrEta(){
    std::vector<float> lepEta;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepEta.push_back(leplistfsr[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::ElectronFsrPhi(){
    std::vector<float> lepPhi;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = ElectronFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPhi.push_back(leplistfsr[i].Phi());
    }
    return lepPhi;
}

std::vector<float> H4LTools::MuonFsrPt(){
    std::vector<float> lepPt;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPt.push_back(leplistfsr[i].Pt());
    }
    return lepPt;
}

std::vector<float> H4LTools::MuonFsrEta(){
    std::vector<float> lepEta;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepEta.push_back(leplistfsr[i].Eta());
    }
    return lepEta;
}

std::vector<float> H4LTools::MuonFsrPhi(){
    std::vector<float> lepPhi;
    std::vector<TLorentzVector> leplistfsr;
    leplistfsr = MuonFsr();
    for (unsigned int i=0;i<leplistfsr.size();i++){
        lepPhi.push_back(leplistfsr[i].Phi());
    }
    return lepPhi;
}*/

void H4LTools::LeptonSelection(){

    step1Ele = goodLooseElectrons2012();
    step1Mu = goodLooseMuons2012();
    bestEle = goodElectrons2015_noIso_noBdt(step1Ele);
    bestMu = goodMuons2015_noIso_noPf(step1Mu);
    Electronindex = bestEle;
    Muonindex = bestMu;
    AllEid = pass_Ele_Id(nanoVersion);
    AllMuid = pass_Mu_Id("Run3","MVA-Based","loose");
    for (unsigned int iuj=0;iuj<step1Ele.size();iuj++){
        if(AllEid[step1Ele[iuj]]) tighteleforjetidx.push_back(step1Ele[iuj]);
    }
    for (unsigned int juj=0;juj<step1Mu.size();juj++){
        if(AllMuid[step1Mu[juj]]) tightmuforjetidx.push_back(step1Mu[juj]);
    }
    jetidx = SelectedJets(tighteleforjetidx,tightmuforjetidx);

    for(unsigned int ie=0; ie<Electronindex.size();ie++){
        if(Electron_pdgId[Electronindex[ie]]>0){
            Elechg.push_back(-1);
        }
        else{
            Elechg.push_back(1);
        }
        TLorentzVector Ele;
        Ele.SetPtEtaPhiM(Electron_pt[Electronindex[ie]],Electron_eta[Electronindex[ie]],Electron_phi[Electronindex[ie]],Electron_mass[Electronindex[ie]]);
        Elelist.push_back(Ele);
        ElelistFsr.push_back(Electrondressed_Run3[Electronindex[ie]]);
        Eiso.push_back(Electron_pfRelIso03_all[Electronindex[ie]]);
        Eid.push_back(AllEid[Electronindex[ie]]);
    }

    for(unsigned int imu=0; imu<Muonindex.size();imu++){
        if(Muon_pdgId[Muonindex[imu]]>0){
            Muchg.push_back(-1);
        }
        else{
            Muchg.push_back(1);
        }
        TLorentzVector Mu;
        Mu.SetPtEtaPhiM(Muon_pt[Muonindex[imu]],Muon_eta[Muonindex[imu]],Muon_phi[Muonindex[imu]],Muon_mass[Muonindex[imu]]);
        Mulist.push_back(Mu);
        MulistFsr.push_back(Muondressed_Run3[Muonindex[imu]]);
        muid.push_back(AllMuid[Muonindex[imu]]);
        Muiso.push_back(Muon_pfRelIso03_all[Muonindex[imu]]);
    }
    //ElelistFsr = BatchFsrRecovery(Elelist);
    //MulistFsr = BatchFsrRecovery(Mulist);
    for(unsigned int ae=0; ae<Eid.size();ae++){
        float RelEleIsoNoFsr;
        RelEleIsoNoFsr = Eiso[ae];
        unsigned FsrEleidx;
        FsrEleidx = doFsrRecovery_Run3(goodFsrPhotons(), Electronindex[ae], 11);
        if (isFSR && (FsrEleidx < 900)){
            TLorentzVector fsrele;
            fsrele.SetPtEtaPhiM(FsrPhoton_pt[FsrEleidx],FsrPhoton_eta[FsrEleidx],FsrPhoton_phi[FsrEleidx],0);
            if(Elelist[ae].DeltaR(fsrele)>0.01){
                RelEleIsoNoFsr = RelEleIsoNoFsr - FsrPhoton_pt[FsrEleidx]/Elelist[ae].Pt();
            }
          //FsrEleidx = doFsrRecovery(Elelist[ae]);
          /*if(FsrEleidx<900){
              TLorentzVector fsrele;
              fsrele.SetPtEtaPhiM(FsrPhoton_pt[FsrEleidx],FsrPhoton_eta[FsrEleidx],FsrPhoton_phi[FsrEleidx],0);
              std::cout<<"Ele correction: "<< std::endl;
              if(Elelist[ae].DeltaR(fsrele)>0.01){
                RelEleIsoNoFsr = RelEleIsoNoFsr - FsrPhoton_pt[FsrEleidx]/Elelist[ae].Pt();
              }
          }*/
        }
        //if((Eid[ae]==true)&&(RelEleIsoNoFsr<0.35)){
        if(Eid[ae]==true){
            nTightEle++;
            TightEleindex.push_back(ae);
            nTightEleChgSum += Elechg[ae];
            TightElelep_index.push_back(Lepointer);
            Lepointer++;
            if (isMC) lep_genindex.push_back(Electron_genPartIdx[Electronindex[ae]]);
            else lep_genindex.push_back(-1);
        }
    }
    for(unsigned int amu=0; amu<muid.size();amu++){
        float RelIsoNoFsr;
        RelIsoNoFsr = Muiso[amu];
        unsigned int FsrMuonidx;
        FsrMuonidx = doFsrRecovery_Run3(goodFsrPhotons(), Muonindex[amu], 13);
        if (isFSR && (FsrMuonidx < 900)){
            TLorentzVector fsrmuon;
            fsrmuon.SetPtEtaPhiM(FsrPhoton_pt[FsrMuonidx],FsrPhoton_eta[FsrMuonidx],FsrPhoton_phi[FsrMuonidx],0);
            if(Mulist[amu].DeltaR(fsrmuon)>0.01){
                RelIsoNoFsr = RelIsoNoFsr - FsrPhoton_pt[FsrMuonidx]/Mulist[amu].Pt();
            }
          /*if(FsrMuonidx<900){
              TLorentzVector fsrmuon;
              fsrmuon.SetPtEtaPhiM(FsrPhoton_pt[FsrMuonidx],FsrPhoton_eta[FsrMuonidx],FsrPhoton_phi[FsrMuonidx],0);
              std::cout<<"muon FSR recovered"<<endl;
              if(Mulist[amu].DeltaR(fsrmuon)>0.01){
                RelIsoNoFsr = RelIsoNoFsr - FsrPhoton_pt[FsrMuonidx]/Mulist[amu].Pt();
              }
          }*/
        }
        if((muid[amu]==true)&&(RelIsoNoFsr<0.35)){
            nTightMu++;
            TightMuindex.push_back(amu);
            nTightMuChgSum += Muchg[amu];
            TightMulep_index.push_back(Lepointer);
            Lepointer++;
            if (isMC) lep_genindex.push_back(Muon_genPartIdx[Muonindex[amu]]);
            else lep_genindex.push_back(-1);
        }
    }
    if (analysisMode == "2l2j") {
        if (nTightEle >= 2 || nTightMu >= 2) {
            eventPassTwoTightLeps = true;
            passTwoTightLeps++;
        }
    }
}

bool H4LTools::findZCandidate(){
    TLorentzVector z1,z2;
    if (nTightEle>=4) {
        cut4e++;
        flag4e = true;
    }
    else if (nTightMu>=4){
        cut4mu++;
        flag4mu = true;
    }
    else if ((nTightMu>=2)&&(nTightEle>=2)){
        cut2e2mu++;
        flag2e2mu = true;
    }
    //Just record 3 cases, but do not cut anything here.

    if(TightEleindex.size()>1){
        for(unsigned int ke=0; ke<(TightEleindex.size()-1);ke++){
            for(unsigned int je=ke+1;je<TightEleindex.size();je++){
                if ((Elechg[TightEleindex[ke]]+Elechg[TightEleindex[je]])==0){
                    TLorentzVector Zcan;
                    Zcan = ElelistFsr[TightEleindex[ke]] + ElelistFsr[TightEleindex[je]];
                    if((Zcan.M()>MZcutdown)&&(Zcan.M()<MZcutup)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(TightEleindex[ke]);
                        Zlep2index.push_back(TightEleindex[je]);
                        Zlep1lepindex.push_back(TightElelep_index[ke]);
                        Zlep2lepindex.push_back(TightElelep_index[je]);
                        Zflavor.push_back(11);
                        Zlep1pt.push_back(ElelistFsr[TightEleindex[ke]].Pt());
                        Zlep2pt.push_back(ElelistFsr[TightEleindex[je]].Pt());
                        Zlep1eta.push_back(ElelistFsr[TightEleindex[ke]].Eta());
                        Zlep2eta.push_back(ElelistFsr[TightEleindex[je]].Eta());
                        Zlep1phi.push_back(ElelistFsr[TightEleindex[ke]].Phi());
                        Zlep2phi.push_back(ElelistFsr[TightEleindex[je]].Phi());
                        Zlep1mass.push_back(ElelistFsr[TightEleindex[ke]].M());
                        Zlep2mass.push_back(ElelistFsr[TightEleindex[je]].M());
                        Zlep1ptNoFsr.push_back(Elelist[TightEleindex[ke]].Pt());
                        Zlep2ptNoFsr.push_back(Elelist[TightEleindex[je]].Pt());
                        Zlep1etaNoFsr.push_back(Elelist[TightEleindex[ke]].Eta());
                        Zlep2etaNoFsr.push_back(Elelist[TightEleindex[je]].Eta());
                        Zlep1phiNoFsr.push_back(Elelist[TightEleindex[ke]].Phi());
                        Zlep2phiNoFsr.push_back(Elelist[TightEleindex[je]].Phi());
                        Zlep1massNoFsr.push_back(Elelist[TightEleindex[ke]].M());
                        Zlep2massNoFsr.push_back(Elelist[TightEleindex[je]].M());
                        Zlep1chg.push_back(Elechg[TightEleindex[ke]]);
                        Zlep2chg.push_back(Elechg[TightEleindex[je]]);
                    }
                }
            }
        }
    }

    if(TightMuindex.size()>1){
        for(unsigned int kmu=0; kmu<(TightMuindex.size()-1);kmu++){
            for(unsigned int jmu=kmu+1;jmu<TightMuindex.size();jmu++){
                if ((Muchg[TightMuindex[kmu]]+Muchg[TightMuindex[jmu]])==0){
                    TLorentzVector Zcan;
                    Zcan = MulistFsr[TightMuindex[kmu]] + MulistFsr[TightMuindex[jmu]];
                    if((Zcan.M()>MZcutdown)&&(Zcan.M()<MZcutup)){
                        Zlist.push_back(Zcan);
                        Zlep1index.push_back(TightMuindex[kmu]);
                        Zlep2index.push_back(TightMuindex[jmu]);
                        Zlep1lepindex.push_back(TightMulep_index[kmu]);
                        Zlep2lepindex.push_back(TightMulep_index[jmu]);
                        Zflavor.push_back(13);
                        Zlep1pt.push_back(MulistFsr[TightMuindex[kmu]].Pt());
                        Zlep2pt.push_back(MulistFsr[TightMuindex[jmu]].Pt());
                        Zlep1eta.push_back(MulistFsr[TightMuindex[kmu]].Eta());
                        Zlep2eta.push_back(MulistFsr[TightMuindex[jmu]].Eta());
                        Zlep1phi.push_back(MulistFsr[TightMuindex[kmu]].Phi());
                        Zlep2phi.push_back(MulistFsr[TightMuindex[jmu]].Phi());
                        Zlep1mass.push_back(MulistFsr[TightMuindex[kmu]].M());
                        Zlep2mass.push_back(MulistFsr[TightMuindex[jmu]].M());
                        Zlep1ptNoFsr.push_back(Mulist[TightMuindex[kmu]].Pt());
                        Zlep2ptNoFsr.push_back(Mulist[TightMuindex[jmu]].Pt());
                        Zlep1etaNoFsr.push_back(Mulist[TightMuindex[kmu]].Eta());
                        Zlep2etaNoFsr.push_back(Mulist[TightMuindex[jmu]].Eta());
                        Zlep1phiNoFsr.push_back(Mulist[TightMuindex[kmu]].Phi());
                        Zlep2phiNoFsr.push_back(Mulist[TightMuindex[jmu]].Phi());
                        Zlep1massNoFsr.push_back(Mulist[TightMuindex[kmu]].M());
                        Zlep2massNoFsr.push_back(Mulist[TightMuindex[jmu]].M());
                        Zlep1chg.push_back(Muchg[TightMuindex[kmu]]);
                        Zlep2chg.push_back(Muchg[TightMuindex[jmu]]);
                    }
                }
            }
        }
    }
    for (unsigned int znofsr = 0; znofsr<Zlist.size(); znofsr++){
        TLorentzVector Zlep1nofsr,Zlep2nofsr,Zcannofsr;
        Zlep1nofsr.SetPtEtaPhiM(Zlep1ptNoFsr[znofsr],Zlep1etaNoFsr[znofsr],Zlep1phiNoFsr[znofsr],Zlep1massNoFsr[znofsr]);
        Zlep2nofsr.SetPtEtaPhiM(Zlep2ptNoFsr[znofsr],Zlep2etaNoFsr[znofsr],Zlep2phiNoFsr[znofsr],Zlep2massNoFsr[znofsr]);
        Zcannofsr = Zlep1nofsr + Zlep2nofsr;
        Zlistnofsr.push_back(Zcannofsr);
    }
    Zsize = Zlist.size();

    if (analysisMode == "2l2j" && Zsize > 0) {
        eventPassZCand = true;
        passZCand++;
    }

    if (Zsize>0){
        return true;
    }
    else{
        return false;
    }
}

bool H4LTools::BuildBestDijet(){

    if (analysisMode == "2l2j") {
        if ((int)jetidx.size() >= JetNcut) {
            eventPassTwoGoodJets = true;
            passTwoGoodJets++;
        }
    }

    if(jetidx.size() < JetNcut){
        return false;
    }

    // NanoAODv12 / 2022 / 2023: use RobustParT
    // NanoAODv15 / 2024:        use UParT
    bool useUPT = (nanoVersion >= 15);

    const std::vector<float>& Jet_btagParTAK4B =
        useUPT ? Jet_btagUParTAK4B : Jet_btagRobustParTAK4B;

    unsigned int jet1index = 99;
    unsigned int jet2index = 99;

    if(jetidx.size() == 2)
    {
        jet1index = jetidx[0];
        jet2index = jetidx[1];

        if(Jet_btagParTAK4B[jetidx[1]] > Jet_btagParTAK4B[jetidx[0]])
        {
            jet1index = jetidx[1];
            jet2index = jetidx[0];
        }
    }

    if(jetidx.size() > 2)
    {
        jet1index = jetidx[0];
        jet2index = jetidx[1];

        if(Jet_btagParTAK4B[jetidx[1]] > Jet_btagParTAK4B[jetidx[0]])
        {
            jet1index = jetidx[1];
            jet2index = jetidx[0];
        }

        for (unsigned int pj = 2; pj < jetidx.size(); pj++){
            if(
                (Jet_btagParTAK4B[jetidx[pj]] > Jet_btagParTAK4B[jet1index]) &&
                (Jet_btagParTAK4B[jetidx[pj]] > Jet_btagParTAK4B[jet2index])
            ){
                jet2index = jet1index;
                jet1index = jetidx[pj];
            }
            else if(
                (Jet_btagParTAK4B[jetidx[pj]] > Jet_btagParTAK4B[jet2index]) &&
                (Jet_btagParTAK4B[jetidx[pj]] < Jet_btagParTAK4B[jet1index])
            ){
                jet2index = jetidx[pj];
            }
        }
    }

    TLorentzVector Jet1, Jet2;

    if(jetidx.size() > 1){
        Jet1.SetPtEtaPhiM(
            Jet_pt[jet1index],
            Jet_eta[jet1index],
            Jet_phi[jet1index],
            Jet_mass[jet1index]
        );

        pTj1  = Jet1.Pt();
        etaj1 = Jet1.Eta();
        phij1 = Jet1.Phi();
        mj1   = Jet1.M();

        btagger1_DJ = Jet_btagDeepFlavB[jet1index];
        btagger1_PN = Jet_btagPNetB[jet1index];

        Jet2.SetPtEtaPhiM(
            Jet_pt[jet2index],
            Jet_eta[jet2index],
            Jet_phi[jet2index],
            Jet_mass[jet2index]
        );

        pTj2  = Jet2.Pt();
        etaj2 = Jet2.Eta();
        phij2 = Jet2.Phi();
        mj2   = Jet2.M();

        btagger2_DJ = Jet_btagDeepFlavB[jet2index];
        btagger2_PN = Jet_btagPNetB[jet2index];

        if(useUPT){
            // NanoAODv15 / 2024
            btagger1_RPT = -999.;
            btagger2_RPT = -999.;

            btagger1_UPT = Jet_btagUParTAK4B[jet1index];
            btagger2_UPT = Jet_btagUParTAK4B[jet2index];
        }
        else{
            // NanoAODv12 / 2022 / 2023
            btagger1_RPT = Jet_btagRobustParTAK4B[jet1index];
            btagger2_RPT = Jet_btagRobustParTAK4B[jet2index];

            btagger1_UPT = -999.;
            btagger2_UPT = -999.;
        }

        invjj = (Jet1 + Jet2).M();

        if (analysisMode == "2l2j") {
            eventPassDijet = true;
            passDijet++;
        }

        return true;
    }

    return false;
}

bool H4LTools::BuildZZCandidate(){

    bool foundZZCandidate = false;

    if((nTightMu+nTightEle)<4){
        return foundZZCandidate;
    }
    if((abs(nTightEleChgSum)+abs(nTightMuChgSum))>(nTightMu+nTightEle-4)){
        return foundZZCandidate;
    }
    if(Zsize<2){
        return foundZZCandidate;
    }

    //Find ZZ candidate
    std::vector<int> Z1CanIndex;
    std::vector<int> Z2CanIndex;
    int ghosttag = 0, QCDtag=0, lepPtTag = 0;
    for (unsigned int m=0; m<(Zlist.size()-1); m++){
        for (unsigned int n=m+1; n<Zlist.size(); n++){
            if (Zflavor[m]==Zflavor[n]){
               if ((Zlep1index[m] == Zlep1index[n])||(Zlep2index[m] == Zlep1index[n])) continue;  //non-overlapping
               if ((Zlep1index[m] == Zlep2index[n])||(Zlep2index[m] == Zlep2index[n])) continue;
            }

            // Apply ghost removal and the low-mass pair veto to the four bare leptons.
            std::array<h4l::CandidateLepton, 4> fourLeptons;
            fourLeptons[0].bareP4.SetPtEtaPhiM(
                Zlep1ptNoFsr[m], Zlep1etaNoFsr[m], Zlep1phiNoFsr[m], Zlep1massNoFsr[m]
            );
            fourLeptons[0].charge = Zlep1chg[m];
            fourLeptons[1].bareP4.SetPtEtaPhiM(
                Zlep2ptNoFsr[m], Zlep2etaNoFsr[m], Zlep2phiNoFsr[m], Zlep2massNoFsr[m]
            );
            fourLeptons[1].charge = Zlep2chg[m];
            fourLeptons[2].bareP4.SetPtEtaPhiM(
                Zlep1ptNoFsr[n], Zlep1etaNoFsr[n], Zlep1phiNoFsr[n], Zlep1massNoFsr[n]
            );
            fourLeptons[2].charge = Zlep1chg[n];
            fourLeptons[3].bareP4.SetPtEtaPhiM(
                Zlep2ptNoFsr[n], Zlep2etaNoFsr[n], Zlep2phiNoFsr[n], Zlep2massNoFsr[n]
            );
            fourLeptons[3].charge = Zlep2chg[n];

            if (!h4l::passesGhostRemoval(fourLeptons)) continue;
            ghosttag++;
            bool nPassPt20;
            int nPassPt10;
            nPassPt20 = (Zlep1pt[m]>20) || (Zlep2pt[m]>20) || (Zlep1pt[n]>20) || (Zlep2pt[n]>20);
            nPassPt10 = 0;
            if (Zlep1pt[m]>10) nPassPt10 += 1;
            if (Zlep2pt[m]>10) nPassPt10 += 1;
            if (Zlep1pt[n]>10) nPassPt10 += 1;
            if (Zlep2pt[n]>10) nPassPt10 += 1;
            if (nPassPt10 < 2) continue;
            if (nPassPt20 == false) continue; //lep Pt requirements
            lepPtTag++;
            if (!h4l::passesOppositeSignPairMass(fourLeptons)) continue;
            QCDtag++;

            const std::pair<unsigned int, unsigned int> orderedZs =
                h4l::orderZCandidates(m, n, Zlist[m], Zlist[n], Zmass);
            const unsigned int z1Index = orderedZs.first;
            const unsigned int z2Index = orderedZs.second;
            const TLorentzVector& zZ1 = Zlist[z1Index];
            const TLorentzVector& zZ2 = Zlist[z2Index];
            if (!h4l::passesZ1Mass(zZ1, MZ1cut)) continue;

            bool passSmartCut = true;
            if (Zflavor[m]==Zflavor[n]){
                TLorentzVector Za,Zb,lepM1,lepM2,lepN1,lepN2;
                int lepM1chg,lepN1chg;
                lepM1.SetPtEtaPhiM(Zlep1pt[m],Zlep1eta[m],Zlep1phi[m],Zlep1mass[m]);
                lepM2.SetPtEtaPhiM(Zlep2pt[m],Zlep2eta[m],Zlep2phi[m],Zlep2mass[m]);
                lepN1.SetPtEtaPhiM(Zlep1pt[n],Zlep1eta[n],Zlep1phi[n],Zlep1mass[n]);
                lepN2.SetPtEtaPhiM(Zlep2pt[n],Zlep2eta[n],Zlep2phi[n],Zlep2mass[n]);
                lepM1chg = Zlep1chg[m];
                lepN1chg = Zlep1chg[n];
                if(lepM1chg == lepN1chg){
                    Za = lepM1 + lepN2;
                    Zb = lepN1 + lepM2;
                }
                else{
                    Za = lepM1 + lepN1;
                    Zb = lepN2 + lepM2;
                }
                if (fabs(Za.M()-Zmass)<fabs(Zb.M()-Zmass)){
                    if ( (fabs(Za.M()-Zmass)<fabs(zZ1.M()-Zmass)) && (Zb.M()<MZcutdown) ) passSmartCut=false;
                }
                else{
                    if ( (fabs(Zb.M()-Zmass)<fabs(zZ1.M()-Zmass)) && (Za.M()<MZcutdown) ) passSmartCut=false;
                }
            }
            if (passSmartCut==false) continue ;
            // MZZcut is the lower threshold on the four-lepton invariant mass.
            if (!h4l::passesFourLeptonMass(zZ1, zZ2, MZZcut)) continue;
            foundZZCandidate = true;
            Z1CanIndex.push_back(z1Index);
            Z2CanIndex.push_back(z2Index);
        }
    }
    if(ghosttag){
        if (flag2e2mu) {
            cutghost2e2mu++;
        }
        if (flag4e){
            cutghost4e++;
        }
        if (flag4mu) {
            cutghost4mu++;
        }
    }
    if(lepPtTag){
        if (flag2e2mu) {
            cutLepPt2e2mu++;
        }
        if (flag4e) {
            cutLepPt4e++;
        }
        if (flag4mu) {
            cutLepPt4mu++;
        }
    }
    if(QCDtag){
        if (flag2e2mu) {
            cutQCD2e2mu++;
        }
        if (flag4e) {
            cutQCD4e++;
        }
        if (flag4mu) {
            cutQCD4mu++;
        }
    }
    if(foundZZCandidate == false){
        return foundZZCandidate;
    }
    if (flag2e2mu) {
        cutZZ2e2mu++;
    }
    if (flag4e) {
        cutZZ4e++;
    }
    if (flag4mu) {
        cutZZ4mu++;
    }
    int Z1index = Z1CanIndex[0];
    int Z2index = Z2CanIndex[0];
    double Z1Distance = fabs(Zlist[Z1index].M()-Zmass);
    double Z2Ptsum = Zlep1pt[Z2index] + Zlep2pt[Z2index];
    if(Z1CanIndex.size()>1){
        for(unsigned int iz=1;iz<Z1CanIndex.size();iz++){
            const double candidateZ1Distance = fabs(Zlist[Z1CanIndex[iz]].M()-Zmass);
            const double candidateZ2Ptsum =
                Zlep1pt[Z2CanIndex[iz]] + Zlep2pt[Z2CanIndex[iz]];
            if(h4l::isBetterZZCandidate(
                candidateZ1Distance, candidateZ2Ptsum, Z1Distance, Z2Ptsum
            )){
                Z1index = Z1CanIndex[iz];
                Z2index = Z2CanIndex[iz];
                Z1Distance = candidateZ1Distance;
                Z2Ptsum = candidateZ2Ptsum;
            }
        }
    }
    Z1 = Zlist[Z1index];
    Z2 = Zlist[Z2index];

    Z1nofsr = Zlistnofsr[Z1index];
    Z2nofsr = Zlistnofsr[Z2index];
    ZZsystem = Z1+Z2;
    ZZsystemnofsr = Z1nofsr+Z2nofsr;

    /*if(abs(ZZsystemnofsr.M()-ZZsystem.M())>0.000001){
        std::cout<<"FSR works "<<abs(ZZsystemnofsr.M()-ZZsystem.M())<<std::endl;
        std::cout<<"FSR: "<<ZZsystem.M()<<" noFSR:"<<ZZsystemnofsr.M()<<std::endl;
    }*/
    float massZZ;
    if (isFSR) massZZ = ZZsystem.M();
    else massZZ = ZZsystemnofsr.M();
    if ((massZZ>HiggscutDown)&&(massZZ<HiggscutUp)){
        if (flag2e2mu) cutm4l2e2mu++;
        if (flag4e) cutm4l4e++;
        if (flag4mu) cutm4l4mu++;
    }//It doesn’t define PassZZSelection; it just selects events within it that fall in the Higgs mass window.

    SimpleParticleCollection_t associated;

    if(analysisMode == "4l2j"){
        unsigned int jet1index = 99;
        unsigned int jet2index = 99;

        // Use the same b-tag score as BuildBestDijet().
        bool useUPT = (nanoVersion >= 15);

        const std::vector<float>& Jet_btagParTAK4B =
            useUPT ? Jet_btagUParTAK4B : Jet_btagRobustParTAK4B;

        if(jetidx.size() == 2)
        {
            jet1index = jetidx[0];
            jet2index = jetidx[1];

            if(Jet_btagParTAK4B[jetidx[1]] > Jet_btagParTAK4B[jetidx[0]])
            {
                jet1index = jetidx[1];
                jet2index = jetidx[0];
            }
        }

        if(jetidx.size() > 2)
        {
            jet1index = jetidx[0];
            jet2index = jetidx[1];

            if(Jet_btagParTAK4B[jetidx[1]] > Jet_btagParTAK4B[jetidx[0]])
            {
                jet1index = jetidx[1];
                jet2index = jetidx[0];
            }

            for (unsigned int pj = 2; pj < jetidx.size(); pj++){
                if(
                    (Jet_btagParTAK4B[jetidx[pj]] > Jet_btagParTAK4B[jet1index]) &&
                    (Jet_btagParTAK4B[jetidx[pj]] > Jet_btagParTAK4B[jet2index])
                ){
                    jet2index = jet1index;
                    jet1index = jetidx[pj];
                }
                else if(
                    (Jet_btagParTAK4B[jetidx[pj]] > Jet_btagParTAK4B[jet2index]) &&
                    (Jet_btagParTAK4B[jetidx[pj]] < Jet_btagParTAK4B[jet1index])
                ){
                    jet2index = jetidx[pj];
                }
            }
        }

        if(jetidx.size() > 1){
            TLorentzVector Jet1, Jet2;

            Jet1.SetPtEtaPhiM(
                Jet_pt[jet1index],
                Jet_eta[jet1index],
                Jet_phi[jet1index],
                Jet_mass[jet1index]
            );
            associated.push_back(SimpleParticle_t(0, Jet1));

            Jet2.SetPtEtaPhiM(
                Jet_pt[jet2index],
                Jet_eta[jet2index],
                Jet_phi[jet2index],
                Jet_mass[jet2index]
            );
            associated.push_back(SimpleParticle_t(0, Jet2));
        }
    }

    SimpleParticleCollection_t daughters;
    TLorentzVector Lep1,Lep2,Lep3,Lep4;

    Lep1.SetPtEtaPhiM(Zlep1pt[Z1index],Zlep1eta[Z1index],Zlep1phi[Z1index],Zlep1mass[Z1index]);
    Lep2.SetPtEtaPhiM(Zlep2pt[Z1index],Zlep2eta[Z1index],Zlep2phi[Z1index],Zlep2mass[Z1index]);
    Lep3.SetPtEtaPhiM(Zlep1pt[Z2index],Zlep1eta[Z2index],Zlep1phi[Z2index],Zlep1mass[Z2index]);
    Lep4.SetPtEtaPhiM(Zlep2pt[Z2index],Zlep2eta[Z2index],Zlep2phi[Z2index],Zlep2mass[Z2index]);

    if ((Zflavor[Z1index]==11)&&(Zflavor[Z2index]==11)) RecoFourEEvent=true;
    if ((Zflavor[Z1index]==13)&&(Zflavor[Z2index]==13)) RecoFourMuEvent=true;
    if ((Zflavor[Z1index]==11)&&(Zflavor[Z2index]==13)) RecoTwoETwoMuEvent=true;
    if ((Zflavor[Z1index]==13)&&(Zflavor[Z2index]==11)) RecoTwoMuTwoEEvent=true;
    lep_Hindex[0] = Zlep1lepindex[Z1index];
    lep_Hindex[1] = Zlep2lepindex[Z1index];
    lep_Hindex[2] = Zlep1lepindex[Z2index];
    lep_Hindex[3] = Zlep2lepindex[Z2index];
    pTL1 = Lep1.Pt();
    etaL1 = Lep1.Eta();
    phiL1 = Lep1.Phi();
    massL1 = Lep1.M();
    pTL2 = Lep2.Pt();
    etaL2 = Lep2.Eta();
    phiL2 = Lep2.Phi();
    massL2 = Lep2.M();
    pTL3 = Lep3.Pt();
    etaL3 = Lep3.Eta();
    phiL3 = Lep3.Phi();
    massL3 = Lep3.M();
    pTL4 = Lep4.Pt();
    etaL4 = Lep4.Eta();
    phiL4 = Lep4.Phi();
    massL4 = Lep4.M();

    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z1index]*Zlep1chg[Z1index], Lep1));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z1index]*Zlep2chg[Z1index], Lep2));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z2index]*Zlep1chg[Z2index], Lep3));
    daughters.push_back(SimpleParticle_t((-1)*Zflavor[Z2index]*Zlep2chg[Z2index], Lep4));
    me_0plus_JHU=999.0; me_qqZZ_MCFM=999.0; p0plus_m4l=999.0; bkg_m4l=999.0; D_bkg_kin=999.0; D_bkg=999.0;
    D_bkg_kin_vtx_BS=999.0;

    p0minus_VAJHU=999.0; pg1g4_VAJHU=999.0; Dgg10_VAMCFM=999.0; D_g4=999.0; D_g1g4=999.0; D_CP=999.0; D_int=999.0; D_L1_int=999.0; D_L1Zgint=999.0;
    D_0m=999.0; D_0hp=999.0; D_L1=999.0; D_L1Zg=999.0;
    p0plus_VAJHU=9999.0; p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen=999.0; pDL1_VAJHU=999.0; pD_L1Zgint=999.0; p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen=999.0; p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen=999.0, p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen=999.0, p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen=999.0, p0plus_VAJHU=999.0;

    mela->setInputEvent(&daughters, &associated, 0, 0);
    mela->setCurrentCandidateFromIndex(0);
    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(me_0plus_JHU, true);

    mela->setProcess(TVar::H0minus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0minus_VAJHU, true);
    // additional probabilities   GG_SIG_ghg2_1_ghz2_1_JHUGen
    mela->setProcess(TVar::H0hplus, TVar::JHUGen, TVar::ZZGG);
    mela->computeP(p0plus_VAJHU, true);

    // p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, Couplings:ghg2=1,0;ghz1=1,0;ghz2=1,0 Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz2_1_JHUGen
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_2][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen, true);    //FIXME

    p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen -= p0plus_VAJHU+me_0plus_JHU;


    // p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen
    // Couplings:ghg2=1,0;ghz1_prime2=10000,0
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1_PRIME2][0]=10000.;
    // (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen, true);    //FIXME

    // p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen/1e8, ghg2=1,0;ghz1=1,0;ghz1_prime2=10000,0, Options:SubtractP=GG_SIG_ghg2_1_ghz1_1_JHUGen,GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1_PRIME2][0]=10000.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen, true);    //FIXME
    p_GG_SIG_ghg2_1_ghz1_1_ghz1prime2_1E4_JHUGen -= p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen+me_0plus_JHU;

    // p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, ghg2=1,0;ghzgs1_prime2=10000,0
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_ZA_1_PRIME2][0]=10000.;
    // (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen, true);    //FIXME

    // p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, ghg2=1,0;ghzgs1_prime2=10000,0
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][gHIGGS_GG_2][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_VV_1][0]=1.;
    (mela->selfDHzzcoupl)[0][gHIGGS_ZA_1_PRIME2][0]=10000.;
    mela->computeP(p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen, true);    //FIXME

    p_GG_SIG_ghg2_1_ghz1_1_ghza1prime2_1E4_JHUGen -= me_0plus_JHU+p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen;

    // pg1g4_VAJHU=0.0;
    mela->setProcess(TVar::SelfDefine_spin0, TVar::JHUGen, TVar::ZZGG);
    (mela->selfDHggcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][0][0]=1.;
    (mela->selfDHzzcoupl)[0][3][0]=1.;
    mela->computeP(pg1g4_VAJHU, true);

    pg1g4_VAJHU -= me_0plus_JHU+p0minus_VAJHU;

    mela->setProcess(TVar::bkgZZ, TVar::MCFM, TVar::ZZQQB);
    mela->computeP(me_qqZZ_MCFM, true);

    mela->computeD_gg(TVar::MCFM, TVar::D_gg10, Dgg10_VAMCFM);

    mela->setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::ZZGG);
    mela->computePM4l(TVar::SMSyst_None, p0plus_m4l);

    mela->setProcess(TVar::bkgZZ, TVar::JHUGen, TVar::ZZGG);
    mela->computePM4l(TVar::SMSyst_None, bkg_m4l);
    D_0m = me_0plus_JHU / (me_0plus_JHU + (p0minus_VAJHU * pow(getDg4Constant(massZZ),2)));
    D_CP = pg1g4_VAJHU / (2 * sqrt(me_0plus_JHU * p0minus_VAJHU ));
    D_0hp = me_0plus_JHU / (me_0plus_JHU + (p0plus_VAJHU * pow(getDg2Constant(massZZ),2)));
    D_int = p_GG_SIG_ghg2_1_ghz1_1_ghz2_1_JHUGen / (2 * sqrt(me_0plus_JHU * p0plus_VAJHU));
    D_L1 = me_0plus_JHU / (me_0plus_JHU + ((p_GG_SIG_ghg2_1_ghz1prime2_1E4_JHUGen/1e8) * pow(getDL1Constant(massZZ),2)));
    D_L1Zg = me_0plus_JHU / (me_0plus_JHU + ((p_GG_SIG_ghg2_1_ghza1prime2_1E4_JHUGen/1e8) * pow(getDL1ZgsConstant(massZZ),2)));
    mela->resetInputEvent();

    return foundZZCandidate;
}

bool H4LTools::ZZSelection(){

    bool foundZZCandidate = false;
    if(!findZCandidate()){
        return foundZZCandidate;
    }

    if(analysisMode == "4l"){
        foundZZCandidate = BuildZZCandidate();
        return foundZZCandidate;
    }

    if(analysisMode == "2l2j"){

        if(Zsize<1){
            return false;
        }
        // --------------------------------------------------
        // Define a 2l2j Z candidate as soon as findZCandidate()
        // Choose the candidate whose mass is closest to nominal Z mass.
        // --------------------------------------------------
        unsigned int bestZIdx = 0;
        float bestDm = fabs(Zlist[0].M() - Zmass);

        for(unsigned int iz=1; iz<Zlist.size(); iz++){
            float dm = fabs(Zlist[iz].M() - Zmass);
            if(dm < bestDm){
                bestDm = dm;
                bestZIdx = iz;
            }
        }

        Z1 = Zlist[bestZIdx];
        Z1nofsr = Zlistnofsr[bestZIdx];
        Z1flav= Zflavor[bestZIdx];

        pTL1 = Zlep1pt[bestZIdx];
        etaL1 = Zlep1eta[bestZIdx];
        phiL1 = Zlep1phi[bestZIdx];
        massL1 = Zlep1mass[bestZIdx];

        pTL2 = Zlep2pt[bestZIdx];
        etaL2 = Zlep2eta[bestZIdx];
        phiL2 = Zlep2phi[bestZIdx];
        massL2 = Zlep2mass[bestZIdx];

        // -------Check data/MC agreement in 2l2j mode------
        if(Z1.Pt() <= 40.0){
            return false;
        }
        // --------------------------------------------------

        if(nRawJetsThisEvent >= 2){
            eventPassAtLeastTwoRawJets = true;
            passAtLeastTwoRawJets++;
        }
        if(nPtEtaJetsThisEvent >= 2){
            eventPassAtLeastTwoPtEtaJets = true;
            passAtLeastTwoPtEtaJets++;
        }
        if(nJetIdJetsThisEvent >= 2){
            eventPassAtLeastTwoJetIdJets = true;
            passAtLeastTwoJetIdJets++;
        }
        if(nPuIdJetsThisEvent >= 2){
            eventPassAtLeastTwoPuIdJets = true;
            passAtLeastTwoPuIdJets++;
        }

        foundZZCandidate = BuildBestDijet();

        if (foundZZCandidate) {
            eventPassFinal = true;
            passFinal++;
        }

        return foundZZCandidate;
    }

    if(analysisMode == "4l2j"){
      
        foundZZCandidate = BuildZZCandidate();
        if(foundZZCandidate == false){
            return false;
        }
        if(jetidx.size()<JetNcut){
            return false;
        }
        bool passDijet = BuildBestDijet();
        if(passDijet == false){
            return false;
        }
        return true;
    }

    std::cerr << "[H4LTools] Unknown analysisMode = " << analysisMode << std::endl;
    return false;
}

float H4LTools::getDg4Constant(float ZZMass){
    return spline_g4->Eval(ZZMass);
}

float H4LTools::getDg2Constant(float ZZMass){
    return spline_g2->Eval(ZZMass);
}

float H4LTools::getDL1Constant(float ZZMass){
    return spline_L1->Eval(ZZMass);
}

float H4LTools::getDL1ZgsConstant(float ZZMass){
    return spline_L1Zgs->Eval(ZZMass);
}
