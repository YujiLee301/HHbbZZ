#ifndef GenAnalysis_h
#define GenAnalysis_h

#include <TLorentzVector.h>
#include <vector>

class GenAnalysis{
    public:
      GenAnalysis(); //Importing Ficducial Space cuts
      std::vector<float> GENlep_pt;
      std::vector<float> GENlep_eta;
      std::vector<float> GENlep_phi;
      std::vector<float> GENlep_mass;
      std::vector<int> GENZ_MomId;
      std::vector<float> GENjet_pt;
      std::vector<float> GENjet_eta;
      std::vector<float> GENjet_phi;
      std::vector<float> GENjet_mass;
      std::vector<int> GENjet_hadronFlavour;
      std::vector<float> GENlep_RelIso;
      std::vector<int> GENlep_id;
      std::vector<int> GENlep_MomId;
      std::vector<int> GENlep_MomMomId;
      float GENmassZZ, GENpTZZ;
      float Zmass=91.1876;
      double genIsoConeSizeEl, genIsoConeSizeMu;
      float genIsoCutEl, genIsoCutMu;
      float GENmass4l,GENmass2j;
      float GENpTj1,GENetaj1,GENphij1,GENmj1;
      float GENpTj2,GENetaj2,GENphij2,GENmj2;
      float GENpT4l,GENrapidity4l;
      float leadingPtCut,subleadingPtCut;
      bool passedFiducialSelection;
      int GENnjets_pt30_eta4p7,nGENLeptons;
      int GENZ_DaughtersId[2];
      int nVECZ;
      int GENlep_Hindex[4];
      int GENjet_Hindex[2];
      int flag4e, flag4mu, flag2e2mu,flagpassZ1,flagpassFid;
      int nGEN4e, nGEN4mu, nGEN2e2mu,nGEN4epassZ1,nGEN4epassFid,nGEN4mupassZ1,nGEN4mupassFid,nGEN2e2mupassZ1,nGEN2e2mupassFid;
      
      //======= K-factor variables =======
      float GEN_sqrt_s_hat = -1;      // sqrt(s_hat) for EWK corrections
      float GEN_t_hat = -1;           // t_hat for EWK corrections  
      int GEN_quark_type = 0;         // Initial state quark type (1=d, 2=u, 3=s, 4=c, 5=b)
      float GEN_dPhiZZ = -1;          // |dPhi(Z1,Z2)| for QCD k-factors
      int GEN_final_state = -1;       // 1=4e/4mu/4tau, 2=2e2mu/2e2tau
      //===========================================
      
      void SetGenParts(float GenPart_pt_, float GenPart_eta_,float GenPart_phi_,float GenPart_mass_,int GenPart_pdgId_,int GenPart_status_,int GenPart_genPartIdxMother_){
        GenPart_pt.push_back(GenPart_pt_);
        GenPart_eta.push_back(GenPart_eta_);
        GenPart_phi.push_back(GenPart_phi_);
        GenPart_mass.push_back(GenPart_mass_);
        GenPart_pdgId.push_back(GenPart_pdgId_);
        GenPart_status.push_back(GenPart_status_);
        GenPart_genPartIdxMother.push_back(GenPart_genPartIdxMother_);
      }

      void SetGenJets(float GenJet_pt_, float GenJet_eta_,float GenJet_phi_,float GenJet_mass_, int GenJet_hadronFlavour_){
        GenJet_pt.push_back(GenJet_pt_);
        GenJet_eta.push_back(GenJet_eta_);
        GenJet_phi.push_back(GenJet_phi_);
        GenJet_mass.push_back(GenJet_mass_);
        GenJet_hadronFlavour.push_back(GenJet_hadronFlavour_);
      }
      //---------DY 0to40----------
      //void SetEventWeights(float lhe_vpt) { 
      //    current_LHE_Vpt = lhe_vpt; 
      //}
      //---------DY 0to40----------

      void Initialize(){
        passedFiducialSelection=false;
        nGENLeptons=0; GENmassZZ= 0; GENpTZZ= 0; GENnjets_pt30_eta4p7=0;
        nVECZ=0;GENmass4l=-99;GENmass2j=-99;
        GENpT4l=0;GENrapidity4l=-99;
        GENpTj1=-99;GENetaj1=-99;GENphij1=-99;GENmj1=-99;
        GENpTj2=-99;GENetaj2=-99;GENphij2=-99;GENmj2=-99;
        GENZ_DaughtersId[0]=0;GENZ_DaughtersId[1]=0;
        for (int i=0; i<4; i++) {GENlep_Hindex[i]=-1;}
        for (int i=0; i<2; i++) {GENjet_Hindex[i]=-1;}

        //======= Reset k-factor variables =======
        GEN_sqrt_s_hat = -1;
        GEN_t_hat = -1;
        GEN_quark_type = 0;
        GEN_dPhiZZ = -1;
        GEN_final_state = -1;
        //=========================================

        GenPart_pt.clear(); GenPart_eta.clear(); GenPart_phi.clear(); GenPart_mass.clear(); GenPart_pdgId.clear();GenPart_status.clear(); GenPart_genPartIdxMother.clear();
        GenJet_pt.clear(); GenJet_eta.clear(); GenJet_phi.clear(); GenJet_mass.clear();GenJet_hadronFlavour.clear();
        GENZ_MomId.clear();
        GENjet_pt.clear();GENjet_eta.clear();GENjet_phi.clear();GENjet_mass.clear();GENjet_hadronFlavour.clear();
        GENlep_eta.clear();GENlep_pt.clear();GENlep_phi.clear();GENlep_mass.clear();GENlep_id.clear();GENlep_MomMomId.clear();GENlep_MomId.clear();GENlep_RelIso.clear();
        flag4e=0; flag4mu=0; flag2e2mu=0;flagpassZ1=0;flagpassFid=0;
      }
      int motherID(int Genidx);
      int mothermotherID(int Genidx);
      void SetGenVariables();
      bool mZ1_mZ2(unsigned int& L1, unsigned int& L2, unsigned int& L3, unsigned int& L4, bool makeCuts);


    private:
      std::vector<float> GenPart_pt;
      std::vector<float> GenPart_eta;
      std::vector<float> GenPart_phi;
      std::vector<float> GenPart_mass;
      std::vector<int> GenPart_pdgId;
      std::vector<int> GenPart_status;
      std::vector<int> GenPart_genPartIdxMother;

      std::vector<float> GenJet_pt;
      std::vector<float> GenJet_eta;
      std::vector<float> GenJet_phi;
      std::vector<float> GenJet_mass;
      std::vector<int> GenJet_hadronFlavour;
};
GenAnalysis::GenAnalysis(){
  // FIXME: Add the values to the yaml file
    genIsoConeSizeEl=0.3; genIsoConeSizeMu=0.3;
    genIsoCutEl=0.35; genIsoCutMu=0.35;
    leadingPtCut=20;subleadingPtCut=10;
    nGEN4e=0; nGEN4mu=0; nGEN2e2mu=0;nGEN4epassZ1=0;nGEN4epassFid=0;nGEN4mupassZ1=0;nGEN4mupassFid=0;nGEN2e2mupassZ1=0;nGEN2e2mupassFid=0;
} 
#endif
