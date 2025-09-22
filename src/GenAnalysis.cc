#include "../interface/GenAnalysis.h"
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <vector>
int GenAnalysis::motherID(int Genidx){
    int ID=0;
    while(GenPart_pdgId[GenPart_genPartIdxMother[Genidx]]!=2212 || abs(GenPart_pdgId[GenPart_genPartIdxMother[Genidx]])!=21 || abs(GenPart_pdgId[GenPart_genPartIdxMother[Genidx]])>6){
        if(GenPart_pdgId[GenPart_genPartIdxMother[Genidx]]!=GenPart_pdgId[Genidx]){
            ID=GenPart_pdgId[GenPart_genPartIdxMother[Genidx]]; return ID;
        }
        else{
            Genidx=GenPart_genPartIdxMother[Genidx];
        }
    }
    return 2212;
}
int GenAnalysis::mothermotherID(int Genidx){
    int ID=0;
    while(GenPart_pdgId[GenPart_genPartIdxMother[Genidx]]!=2212 || abs(GenPart_pdgId[GenPart_genPartIdxMother[Genidx]])!=21 || abs(GenPart_pdgId[GenPart_genPartIdxMother[Genidx]])>6){
        if(GenPart_pdgId[GenPart_genPartIdxMother[Genidx]]!=GenPart_pdgId[Genidx] && GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[Genidx]]]!=GenPart_pdgId[Genidx] && GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[Genidx]]]!=GenPart_pdgId[GenPart_genPartIdxMother[Genidx]] ){
            ID=GenPart_pdgId[GenPart_genPartIdxMother[GenPart_genPartIdxMother[Genidx]]]; return ID;
        }
        else{
            Genidx=GenPart_genPartIdxMother[Genidx];
        }
    }
    return 2212;
}
void GenAnalysis::SetGenVariables(){
    TLorentzVector GENmom1, GENmom2;
    TLorentzVector LS3_Z1_1, LS3_Z1_2, LS3_Z2_1, LS3_Z2_2, GENgoodj1,GENgoodj2,GEN_H1Vec,GEN_H2Vec;
    int GENmom1_id=-999, GENmom2_id=-999;
    int counter_initParticle=0;
    for(unsigned int genpidx=0; genpidx<nGenPart; genpidx++){
        GENlep_MomId.push_back(-1);
        GENlep_MomMomId.push_back(-1);
        if(GenPart_status[genpidx]==21){
            counter_initParticle++;
            if (counter_initParticle==1){
                 GENmom1.SetPtEtaPhiM(GenPart_pt[genpidx],GenPart_eta[genpidx],GenPart_phi[genpidx],GenPart_mass[genpidx]);
                 GENmom1_id=GenPart_pdgId[genpidx];
             }
             if (counter_initParticle==2){
                 GENmom2.SetPtEtaPhiM(GenPart_pt[genpidx],GenPart_eta[genpidx],GenPart_phi[genpidx],GenPart_mass[genpidx]);
                 GENmom2_id=GenPart_pdgId[genpidx];
             }
             if (counter_initParticle>2) {
                 std::cout<< "Warning: more than 2 initial particles found in the event! "<< std::endl;
                 return; 
             }
        }

        if (abs(GenPart_pdgId[genpidx])==11 || abs(GenPart_pdgId[genpidx])==13 || abs(GenPart_pdgId[genpidx])==15){
            if (!(GenPart_status[genpidx]==1 || abs(GenPart_pdgId[genpidx])==15)) continue;
            if (!(abs(motherID(genpidx))==23 || abs(motherID(genpidx))==443 || abs(motherID(genpidx))==553 || abs(motherID(genpidx))==24 || abs(motherID(genpidx))==111 || abs(motherID(genpidx))==221 || abs(motherID(genpidx))==331 ) ) continue;
            //if ((abs(motherID(genpidx))==111 || abs(motherID(genpidx))==221 || abs(motherID(genpidx))==331)&&(abs(mothermotherID(genpidx))==23)){
            //    std::cout << "Keep lepton: PDG ID = " << GenPart_pdgId[genpidx]
            //      << ", status = " << GenPart_status[genpidx]
            //      << ", mother PDG ID = " << motherID(genpidx) 
            //      << ", mother mother PDG ID = "<< mothermotherID(genpidx) << std::endl;
            //}
            if (abs(GenPart_pdgId[genpidx]) == 15) {
                bool hasStableLepDaughter = false;
                for (unsigned int d = 0; d < GenPart_pdgId.size(); ++d){
                    if (GenPart_genPartIdxMother[d] != (int)genpidx) continue;
                    if ((abs(GenPart_pdgId[d])==11 || abs(GenPart_pdgId[d])==13) && GenPart_status[d]==1){
                        hasStableLepDaughter = true;
                        break;
                    }
                }
                if (!hasStableLepDaughter) continue;
            }
            
            nGENLeptons++;
            // Collect FSR photons
            TLorentzVector lep_dressed;
            lep_dressed.SetPtEtaPhiM(GenPart_pt[genpidx],GenPart_eta[genpidx],GenPart_phi[genpidx],GenPart_mass[genpidx]);
            set<int> gen_fsrset;
            for(size_t k=0; k<GenPart_pt.size();k++){
                if( GenPart_status[k] != 1) continue; // stable particles only
                if( GenPart_pdgId[k] != 22) continue; // only photons
                TLorentzVector thisphoton;
                thisphoton.SetPtEtaPhiM(GenPart_pt[k],GenPart_eta[k],GenPart_phi[k],GenPart_mass[k]);
                double this_dR_lgamma = thisphoton.DeltaR(lep_dressed);
                bool idmatch=false;
                if(GenPart_pdgId[GenPart_genPartIdxMother[k]]==GenPart_pdgId[genpidx]) idmatch=true;
                if (!idmatch) continue;
                if(this_dR_lgamma<((abs(GenPart_pdgId[genpidx])==11)?genIsoConeSizeEl:genIsoConeSizeMu)) {//Check the value
                    gen_fsrset.insert(k);
                    TLorentzVector gamma;
                    gamma.SetPtEtaPhiM(GenPart_pt[k],GenPart_eta[k],GenPart_phi[k],GenPart_mass[k]);
                    lep_dressed = lep_dressed+gamma;
                }
            } // Dressed leptons loop
            GENlep_id.push_back( GenPart_pdgId[genpidx]);
            total_GEN_lepton_count_allEvents++;
            if (abs(GenPart_pdgId[genpidx]) == 15) {
                total_GEN_tau_count_allEvents++;
            }
            //std::cout << "Total GEN leptons (all events): " << total_GEN_lepton_count_allEvents << std::endl;
            //std::cout << "Total GEN taus (all events): " << total_GEN_tau_count_allEvents << std::endl;

            GENlep_status.push_back(GenPart_status[genpidx]);
            GENlep_pt.push_back( lep_dressed.Pt() );
            GENlep_eta.push_back( lep_dressed.Eta() );
            GENlep_phi.push_back( lep_dressed.Phi() );
            GENlep_mass.push_back( lep_dressed.M() );
            GENlep_MomId.pop_back();
            GENlep_MomMomId.pop_back();
            GENlep_MomId.push_back(motherID(genpidx));
            GENlep_MomMomId.push_back(mothermotherID(genpidx));

            TLorentzVector thisLep;
            thisLep.SetPtEtaPhiM(lep_dressed.Pt(),lep_dressed.Eta(),lep_dressed.Phi(),lep_dressed.M());
            // GEN iso calculation
            double this_GENiso=0.0;
            for(size_t j=0; j<GenPart_eta.size();j++){
                if( GenPart_status[j] != 1 ) continue; // stable particles only
                if (abs(GenPart_pdgId[j])==12 || abs(GenPart_pdgId[j])==14 || abs(GenPart_pdgId[j])==16) continue; // exclude neutrinos
                if ((abs(GenPart_pdgId[j])==11 || abs(GenPart_pdgId[j])==13)) continue; // exclude leptons
                if (gen_fsrset.find(j)!=gen_fsrset.end()) continue; // exclude particles which were selected as fsr photons
                TLorentzVector thisiso;
                thisiso.SetPtEtaPhiM(GenPart_pt[j],GenPart_eta[j],GenPart_phi[j],GenPart_mass[j]);
                double this_dRvL =thisLep.DeltaR(thisiso);
                if(this_dRvL<((abs(GenPart_pdgId[genpidx])==11)?genIsoConeSizeEl:genIsoConeSizeMu)) {//check values
                    this_GENiso = this_GENiso + GenPart_pt[j];
                }
            } // GEN iso loop
            this_GENiso = this_GENiso/thisLep.Pt();
            GENlep_RelIso.push_back(this_GENiso);
            // END GEN iso calculation
        }//leptons

        if (GenPart_pdgId[genpidx]==25) {
            GENMH=GenPart_mass[genpidx];
            GENH_pt.push_back(GenPart_pt[genpidx]);
            GENH_eta.push_back(GenPart_eta[genpidx]);
            GENH_phi.push_back(GenPart_phi[genpidx]);
            GENH_mass.push_back(GenPart_mass[genpidx]);
        }

        if ((GenPart_pdgId[genpidx]==23 || GenPart_pdgId[genpidx]==443 || GenPart_pdgId[genpidx]==553) && (GenPart_status[genpidx]>=20 && GenPart_status[genpidx]<30) ){
            GENZ_pt.push_back(GenPart_pt[genpidx]);
            GENZ_eta.push_back(GenPart_eta[genpidx]);
            GENZ_phi.push_back(GenPart_phi[genpidx]);
            GENZ_mass.push_back(GenPart_mass[genpidx]);
            GENZ_MomId.push_back(motherID(genpidx));
            nVECZ++;
        }
    }
    if (GENlep_pt.size()>=4) {

        unsigned int L1_nocuts=99; unsigned int L2_nocuts=99; unsigned int L3_nocuts=99; unsigned int L4_nocuts=99;
        bool passedFiducialSelectionNoCuts = mZ1_mZ2(L1_nocuts, L2_nocuts, L3_nocuts, L4_nocuts, false);//makecuts=false
        if (passedFiducialSelectionNoCuts) {
            TLorentzVector Z1_1, Z1_2, Z2_1, Z2_2;
            Z1_1.SetPtEtaPhiM(GENlep_pt[L1_nocuts],GENlep_eta[L1_nocuts],GENlep_phi[L1_nocuts],GENlep_mass[L1_nocuts]);
            Z1_2.SetPtEtaPhiM(GENlep_pt[L2_nocuts],GENlep_eta[L2_nocuts],GENlep_phi[L2_nocuts],GENlep_mass[L2_nocuts]);
            Z2_1.SetPtEtaPhiM(GENlep_pt[L3_nocuts],GENlep_eta[L3_nocuts],GENlep_phi[L3_nocuts],GENlep_mass[L3_nocuts]);
            Z2_2.SetPtEtaPhiM(GENlep_pt[L4_nocuts],GENlep_eta[L4_nocuts],GENlep_phi[L4_nocuts],GENlep_mass[L4_nocuts]);
            GENmassZZ = (Z1_1+Z1_2+Z2_1+Z2_2).M();
            GENpTZZ = (Z1_1+Z1_2+Z2_1+Z2_2).Pt();
            int genfs;
            if (abs(GENlep_id[L1_nocuts])==abs(GENlep_id[L3_nocuts])) genfs=1;
            else genfs=2;
        }

    }
    /////// DO THE FIDUCIAL VOLUME CALCULATION //////////////
    passedFiducialSelection=false;
    int nFiducialLeptons = 0;
    int nFiducialPtLead=0;
    int nFiducialPtSublead=0;

    int nGENe, nGENmu;
    nGENe=0; nGENmu=0;
    for (unsigned int i=0; i<GENlep_id.size(); ++i) {
        TLorentzVector thisLep;
        thisLep.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);

        if ( ( (abs(GENlep_id[i]) == 13 && thisLep.Pt() > 5.0 && abs(thisLep.Eta()) < 2.4)
            || (abs(GENlep_id[i]) == 11 && thisLep.Pt() > 7.0 && abs(thisLep.Eta()) < 2.5) )
            && GENlep_RelIso[i]<((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu) ) {
            nFiducialLeptons++;
            if (thisLep.Pt()>leadingPtCut) nFiducialPtLead++;
            if (thisLep.Pt()>subleadingPtCut) nFiducialPtSublead++;
        }
        if (abs(GENlep_id[i]) == 13) nGENmu++;
        if (abs(GENlep_id[i]) == 11) nGENe++;
    }
    if(nGENmu>3) {flag4mu++;nGEN4mu++;}
    if(nGENe>3) {flag4e++;nGEN4e++;}
    if((nGENe>1)&&(nGENmu>1)) {flag2e2mu++;nGEN2e2mu++;}
    if (nFiducialLeptons>=4 && nFiducialPtLead>=1 && nFiducialPtSublead>=2 ){
        // START FIDUCIAL EVENT TOPOLOGY CUTS
        unsigned int L1=99; unsigned int L2=99; unsigned int L3=99; unsigned int L4=99;
        unsigned int j1=99; unsigned int j2=99; unsigned int j3=99; unsigned int j4=99;
        unsigned int h1=99; unsigned int h2=99;
        GENmass4l = -1.0;GENmass2j = -1.0; GENmass4e = -1.0; GENmass4mu = -1.0; GENmass2e2mu = -1.0;
        GENmassZ1 = -1.0; GENmassZ2 = -1.0; GENpT4l = -1.0; GENeta4l = 999.; GENrapidity4l = 999.; GENphi4l = 999.;
        GENpT4lj = -1.0; GENpT4ljj=-1.0; GENmass4lj = -1.0; GENmass4ljj=-1.0;

        passedFiducialSelection = mZ1_mZ2(L1, L2, L3, L4, true);//makecuts=true
        if(flag2e2mu){            
            if(flagpassZ1){
                nGEN2e2mupassZ1++;
                if(flagpassFid){
                    nGEN2e2mupassFid++;
                }
            }
        }
        if(flag4e){
            if(flagpassZ1){
                nGEN4epassZ1++;
                if(flagpassFid){
                    nGEN4epassFid++;
                }
            }
        }
        if(flag4mu){
            if(flagpassZ1){
                nGEN4mupassZ1++;
                if(flagpassFid){
                    nGEN4mupassFid++;
                }
            }
        }
        GENlep_Hindex[0] = L1; GENlep_Hindex[1] = L2; GENlep_Hindex[2] = L3; GENlep_Hindex[3] = L4;
        if (passedFiducialSelection) {

            //    TLorentzVector LS3_Z1_1, LS3_Z1_2, LS3_Z2_1, LS3_Z2_2;
            LS3_Z1_1.SetPtEtaPhiM(GENlep_pt[L1],GENlep_eta[L1],GENlep_phi[L1],GENlep_mass[L1]);
            LS3_Z1_2.SetPtEtaPhiM(GENlep_pt[L2],GENlep_eta[L2],GENlep_phi[L2],GENlep_mass[L2]);
            LS3_Z2_1.SetPtEtaPhiM(GENlep_pt[L3],GENlep_eta[L3],GENlep_phi[L3],GENlep_mass[L3]);
            LS3_Z2_2.SetPtEtaPhiM(GENlep_pt[L4],GENlep_eta[L4],GENlep_phi[L4],GENlep_mass[L4]);
            GEN_H1Vec = LS3_Z1_1 + LS3_Z1_2 + LS3_Z2_1 + LS3_Z2_2;

            GENmass4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).M();

            if (abs(GENlep_id[L1])==11 && abs(GENlep_id[L3])==11) {GENmass4e = GENmass4l;}
            if (abs(GENlep_id[L1])==13 && abs(GENlep_id[L3])==13) {GENmass4mu = GENmass4l;}
            if ( (abs(GENlep_id[L1])==11 || abs(GENlep_id[L1])==13) &&
                (abs(GENlep_id[L3])==11 || abs(GENlep_id[L3])==13) &&
                (abs(GENlep_id[L1])!=abs(GENlep_id[L3]) ) ) {GENmass2e2mu = GENmass4l;}
            GENpT4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Pt();
            GENeta4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Eta();
            GENphi4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Phi();
            GENrapidity4l = (LS3_Z1_1+LS3_Z1_2+LS3_Z2_1+LS3_Z2_2).Rapidity();
            GENmassZ1 = (LS3_Z1_1+LS3_Z1_2).M();
            GENmassZ2 = (LS3_Z2_1+LS3_Z2_2).M();

            int tmpIdL1,tmpIdL2,tmpIdL3,tmpIdL4;
            TLorentzVector GENL11P4, GENL12P4, GENL21P4, GENL22P4;
            if(GENlep_id[L1] < 0){ GENL11P4.SetPxPyPzE(LS3_Z1_1.Px(),LS3_Z1_1.Py(),LS3_Z1_1.Pz(),LS3_Z1_1.E()); tmpIdL1 = GENlep_id[L1];}
            else{ GENL11P4.SetPxPyPzE(LS3_Z1_2.Px(),LS3_Z1_2.Py(),LS3_Z1_2.Pz(),LS3_Z1_2.E()); tmpIdL1 = GENlep_id[L2];}
            if(GENlep_id[L2] > 0){ GENL12P4.SetPxPyPzE(LS3_Z1_2.Px(),LS3_Z1_2.Py(),LS3_Z1_2.Pz(),LS3_Z1_2.E()); tmpIdL2 = GENlep_id[L2];}
            else{ GENL12P4.SetPxPyPzE(LS3_Z1_1.Px(),LS3_Z1_1.Py(),LS3_Z1_1.Pz(),LS3_Z1_1.E()); tmpIdL2 = GENlep_id[L1];}
            if(GENlep_id[L3] < 0){ GENL21P4.SetPxPyPzE(LS3_Z2_1.Px(),LS3_Z2_1.Py(),LS3_Z2_1.Pz(),LS3_Z2_1.E()); tmpIdL3 = GENlep_id[L3];}
            else{ GENL21P4.SetPxPyPzE(LS3_Z2_2.Px(),LS3_Z2_2.Py(),LS3_Z2_2.Pz(),LS3_Z2_2.E()); tmpIdL3 = GENlep_id[L4];}
            if(GENlep_id[L4] > 0) { GENL22P4.SetPxPyPzE(LS3_Z2_2.Px(),LS3_Z2_2.Py(),LS3_Z2_2.Pz(),LS3_Z2_2.E()); tmpIdL4 = GENlep_id[L4];}
            else{ GENL22P4.SetPxPyPzE(LS3_Z2_1.Px(),LS3_Z2_1.Py(),LS3_Z2_1.Pz(),LS3_Z2_1.E()); tmpIdL4 = GENlep_id[L3];}

        }
        bool passedMassOS = true; bool passedElMuDeltaR = true; bool passedDeltaR = true;
        unsigned int N=GENlep_pt.size();
        for(unsigned int i = 0; i<N; i++) {
            for(unsigned int j = i+1; j<N; j++) {

                // only consider the leptons from Z1 and Z2
                if (!(i==L1 || i==L2 || i==L3 || i==L4)) continue;
                if (!(j==L1 || j==L2 || j==L3 || j==L4)) continue;

                TLorentzVector li, lj;
                li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
                lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);

                TLorentzVector mll = li+lj;

                if(GENlep_id[i]*GENlep_id[j]<0) {
                    if(mll.M()<=4) { passedMassOS = false; break; }
                }

                if(abs(GENlep_id[i]) != abs(GENlep_id[j])) {
                    double deltaR = li.DeltaR(lj);
                    if(deltaR<=0.02) { passedElMuDeltaR = false; break; }
                }
                double deltaRll = li.DeltaR(lj);
                if(deltaRll<=0.02) { passedDeltaR = false; break; }
            }
        }
        GENZ_DaughtersId[0] = abs(GENlep_id[GENlep_Hindex[0]]);
        GENZ_DaughtersId[1] = abs(GENlep_id[GENlep_Hindex[2]]);//Fix me: maybe not correct!
        if(passedMassOS==false || passedElMuDeltaR==false || passedDeltaR==false) passedFiducialSelection=false;
         if (passedFiducialSelection) {
            // DO GEN JETS

            int GENjet1index=0; int GENjet2index=0; int GENjet1index_2p5=0; int GENjet2index_2p5=0;
            TLorentzVector GENJet1, GENJet2, GENJet1_2p5, GENJet2_2p5;
            vector<int> GEN_goodJetsidx;

            for(unsigned genjetidx=0; genjetidx<GenJet_pt.size(); genjetidx++) {

                double pt = GenJet_pt[genjetidx];  double eta = GenJet_eta[genjetidx];
                //if (pt<30.0 || abs(eta)>2.4) continue;
                if (abs(eta)>4.7) continue; 

                bool inDR_pt30_eta4p7 = false;
                unsigned int N=GENlep_pt.size();
                TLorentzVector thisJ;
                thisJ.SetPtEtaPhiM(GenJet_pt[genjetidx],GenJet_eta[genjetidx],GenJet_phi[genjetidx],GenJet_mass[genjetidx]);
                for(unsigned int i = 0; i<N; i++) {
                    if (!(abs(GENlep_id[i])==11 || abs(GENlep_id[i])==13)) continue;
                    TLorentzVector genlep;
                    genlep.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
                    double dR = genlep.DeltaR(thisJ);
                    if(dR<0.4) {
                        inDR_pt30_eta4p7=true;
                    }
                }

                // count number of gen jets which no gen leptons are inside its cone
                if (!inDR_pt30_eta4p7) {
                    GEN_goodJetsidx.push_back(genjetidx);
                    GENnjets_pt30_eta4p7++;
                    GENjet_pt.push_back(GenJet_pt[genjetidx]);
                    GENjet_eta.push_back(GenJet_eta[genjetidx]);
                    GENjet_phi.push_back(GenJet_phi[genjetidx]);
                    GENjet_mass.push_back(GenJet_mass[genjetidx]);
                    GENjet_hadronFlavour.push_back(GenJet_hadronFlavour[genjetidx]);

                    if (pt>GENpt_leadingjet_pt30_eta4p7) {
                        GENpt_leadingjet_pt30_eta4p7=pt;
                    }
                    if (abs(thisJ.Eta())<2.5) {
                        GENnjets_pt30_eta2p5++;
                        if (pt>GENpt_leadingjet_pt30_eta2p5) {
                            GENpt_leadingjet_pt30_eta2p5=pt;
                        }
                    }
                }
            }// loop over gen jets
            passedFiducialSelection = false;

            int j1 = (GEN_goodJetsidx.size() > 0) ? GEN_goodJetsidx[0] : -1;
            int j2 = (GEN_goodJetsidx.size() > 1) ? GEN_goodJetsidx[1] : -1;
            int j3 = (GEN_goodJetsidx.size() > 2) ? GEN_goodJetsidx[2] : -1;
            int j4 = (GEN_goodJetsidx.size() > 3) ? GEN_goodJetsidx[3] : -1;

            b_jets.clear();
            
            for (int j : {j1, j2, j3, j4}) {
                if (j >= 0 && GENjet_hadronFlavour[j] == 5) {
                        b_jets.push_back(j);
                }
            }

            if (b_jets.size() < 2) return;

            int h1 = b_jets[0];
            int h2 = b_jets[1];
            GENjet_Hindex[0] = h1;
            GENjet_Hindex[1] = h2;
            GENgoodj1.SetPtEtaPhiM(GENjet_pt[h1], GENjet_eta[h1], GENjet_phi[h1], GENjet_mass[h1]);
            GENpTj1 = GENgoodj1.Pt();
            GENetaj1 = GENgoodj1.Eta();
            GENphij1 = GENgoodj1.Phi();
            GENmj1 = GENgoodj1.M();
            GENgoodj2.SetPtEtaPhiM(GENjet_pt[h2], GENjet_eta[h2], GENjet_phi[h2], GENjet_mass[h2]);
            GENpTj2 = GENgoodj2.Pt();
            GENetaj2 = GENgoodj2.Eta();
            GENphij2 = GENgoodj2.Phi();
            GENmj2 = GENgoodj2.M();
            GEN_H2Vec = GENgoodj1 + GENgoodj2 ;
            GENmass2j = GEN_H2Vec.M();
            passedFiducialSelection = true;            
            //std::cout << "[DEBUG] b_jets.size() = " << b_jets.size() << std::endl;
            //std::cout << "[DEBUG] GEN_goodJetsidx = ";
            //for (size_t i = 0; i < GEN_goodJetsidx.size(); ++i) {
            //    std::cout << GEN_goodJetsidx[i] << " ";
            //}
            //std::cout << std::endl;

            //std::cout << "[DEBUG] b_jets = ";
            //for (size_t i = 0; i < b_jets.size(); ++i) {
            //    std::cout << b_jets[i] << " ";
            //}
            //std::cout << std::endl;

            //std::cout << "[DEBUG] h1 = " << h1
            //    << ", h2 = " << h2
            //    << ", total GEN jets = " << GENjet_pt.size()
            //    << std::endl;

            //std::cout << "[DEBUG] GENjet_Hindex[0] = " << GENjet_Hindex[0]
            //    << ", GENjet_Hindex[1] = " << GENjet_Hindex[1]
            //    << std::endl;
        }
    }
    return;
}

bool GenAnalysis::mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts)
{
    // Reference: https://github.com/qyguo/UFHZZAnalysisRun2/blob/e51073652729067c8be7c1334f749aae76931a7b/UFHZZ4LAna/src/UFHZZ4LAna.cc#L9009

    double offshell = 999.0; bool findZ1 = false; bool passZ1 = false;

    L1 = 0; L2 = 0;

    unsigned int N = GENlep_pt.size();

    for(unsigned int i=0; i<N; i++){
        for(unsigned int j=i+1; j<=N; j++){


            if((GENlep_id[i]+GENlep_id[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
            lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);


            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[i]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;

                if ( abs(GENlep_id[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[j]>((abs(GENlep_id[j])==11)?genIsoCutEl:genIsoCutMu)) continue;
            }

            TLorentzVector mll = li+lj;

            if(abs(mll.M()-Zmass)<offshell){
                double mZ1 = mll.M();
                L1 = i; L2 = j; findZ1 = true; offshell = abs(mZ1-Zmass);
            }
        }
    }

    TLorentzVector l1, l2;
    l1.SetPtEtaPhiM(GENlep_pt[L1],GENlep_eta[L1],GENlep_phi[L1],GENlep_mass[L1]);
    l2.SetPtEtaPhiM(GENlep_pt[L2],GENlep_eta[L2],GENlep_phi[L2],GENlep_mass[L2]);
    TLorentzVector ml1l2 = l1+l2;

    if(ml1l2.M()>40 && ml1l2.M()<120 && findZ1) passZ1 = true;
    if (!makeCuts) passZ1 = true;
    if (makeCuts && passZ1) flagpassZ1++;

    double pTL34 = 0.0; bool findZ2 = false;
    //bool m4lwindow=false; double window_lo=70.0; double window_hi=140.0;

    //cout<<"findZ2"<<endl;
    for(unsigned int i=0; i<N; i++){
        if(i==L1 || i==L2) continue; // can not be the lep from Z1
        for(unsigned int j=i+1; j<=N; j++){
            if(j==L1 || j==L2) continue; // can not be the lep from Z1
            if((GENlep_id[i]+GENlep_id[j])!=0) continue;

            TLorentzVector li, lj;
            li.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
            lj.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);
            TLorentzVector Z2 = li+lj;

            if (makeCuts) {
                if ( abs(GENlep_id[i]) == 13 && (li.Pt() < 5.0 || abs(li.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[i]) == 11 && (li.Pt() < 7.0 || abs(li.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[i]>((abs(GENlep_id[i])==11)?genIsoCutEl:genIsoCutMu)) continue;

                if ( abs(GENlep_id[j]) == 13 && (lj.Pt() < 5.0 || abs(lj.Eta()) > 2.4)) continue;
                if ( abs(GENlep_id[j]) == 11 && (lj.Pt() < 7.0 || abs(lj.Eta()) > 2.5)) continue;
                if ( GENlep_RelIso[j]>((abs(GENlep_id[j])==11)?genIsoCutEl:genIsoCutMu)) continue;
            }

            if ( (li.Pt()+lj.Pt())>=pTL34 ) {
                double mZ2 = Z2.M();
                if( (mZ2>12 && mZ2<120) || (!makeCuts) ) {
                    L3 = i; L4 = j; findZ2 = true;
                    pTL34 = li.Pt()+lj.Pt();
                    //if (m4l>window_lo && m4l<window_hi) m4lwindow=true;
                } else {
                    // still assign L3 and L4 to this pair if we don't have a passing Z2 yet
                    if (findZ2 == false) {L3 = i; L4 = j;}
                    //cout<<"is not new GEN cand"<<endl;
                }
            }

        } // lj
    } // li

    unsigned int tmp_;
    if(GENlep_pt[L1]<GENlep_pt[L2])    {tmp_=L1;    L1=L2;    L2=tmp_;}
    if(GENlep_pt[L3]<GENlep_pt[L4])    {tmp_=L3;    L3=L4;    L4=tmp_;}
    if(passZ1 && findZ2 && makeCuts) flagpassFid++;
    if(passZ1 && findZ2) return true;
    else return false;

}
