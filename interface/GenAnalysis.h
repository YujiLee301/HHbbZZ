#ifndef GenAnalysis_h
#define GenAnalysis_h

#include <utility>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
#include <TSpline.h>
#include <vector>
//#include "yaml-cpp/yaml.h"
#include "../JHUGenMELA/MELA/interface/Mela.h"

class GenAnalysis{
    public:
      GenAnalysis(); //Importing Ficducial Space cuts
      std::vector<float> GENlep_pt;
      std::vector<float> GENlep_eta;
      std::vector<float> GENlep_phi;
      std::vector<float> GENlep_mass;
      std::vector<float> GENZ_pt;
      std::vector<float> GENZ_eta;
      std::vector<float> GENZ_phi;
      std::vector<float> GENZ_mass;
      std::vector<int> GENZ_MomId;
      std::vector<float> GENH_pt;
      std::vector<float> GENH_eta;
      std::vector<float> GENH_phi;
      std::vector<float> GENH_mass;
      std::vector<float> GENjet_pt;
      std::vector<float> GENjet_eta;
      std::vector<float> GENjet_phi;
      std::vector<float> GENjet_mass;
      std::vector<int> GENjet_hadronFlavour;
      std::vector<float> GENlep_RelIso;
      std::vector<int> GENlep_id;
      std::vector<int> GENlep_status;
      std::vector<int> GENlep_MomId;
      std::vector<int> GENlep_MomMomId;
      float GENMH, GENmassZZ, GENpTZZ;
      float Zmass=91.1876;
      double genIsoConeSizeEl, genIsoConeSizeMu;
      float genIsoCutEl, genIsoCutMu;
      float GENmass4l,GENmass2j,GENmass4e,GENmass4mu, GENmass2e2mu;
      float GENpTj1,GENetaj1,GENphij1,GENmj1;
      float GENpTj2,GENetaj2,GENphij2,GENmj2;
      float GENmassZ1,GENmassZ2,GENpT4l,GENeta4l,GENrapidity4l,GENphi4l;
      float GENpt_leadingjet_pt30_eta4p7,GENpt_leadingjet_pt30_eta2p5;
      float GENpT4lj, GENpT4ljj, GENmass4lj, GENmass4ljj;
      float leadingPtCut,subleadingPtCut;
      bool passedFiducialSelection;
      int GENnjets_pt30_eta4p7,GENnjets_pt30_eta2p5,nGENLeptons;
      int GENZ_DaughtersId[2];
      int nVECZ;
      int GENlep_Hindex[4];
      int GENjet_Hindex[2];
      int flag4e, flag4mu, flag2e2mu,flagpassZ1,flagpassFid;
      int nGEN4e, nGEN4mu, nGEN2e2mu,nGEN4epassZ1,nGEN4epassFid,nGEN4mupassZ1,nGEN4mupassFid, nGEN2e2mupTEtaisocuts,nGEN2e2mupassZ1,nGEN2e2mupassFid;
      void SetGenParts(float GenPart_pt_, float GenPart_eta_,float GenPart_phi_,float GenPart_mass_,int GenPart_pdgId_,int GenPart_status_,int GenPart_statusFlags_,int GenPart_genPartIdxMother_){
        GenPart_pt.push_back(GenPart_pt_);
        GenPart_eta.push_back(GenPart_eta_);
        GenPart_phi.push_back(GenPart_phi_);
        GenPart_mass.push_back(GenPart_mass_);
        GenPart_pdgId.push_back(GenPart_pdgId_);
        GenPart_statusFlags.push_back(GenPart_statusFlags_);
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
      void SetObjectNumGen(unsigned nGenPart_, unsigned nGenJet_){
        nGenPart = nGenPart_;
        nGenJet = nGenJet_;
      }
      void Initialize(){
        passedFiducialSelection=false;
        nGenPart = 0; nGENLeptons=0;nGenJet = 0; GENMH = 0; GENmassZZ= 0; GENpTZZ= 0; GENnjets_pt30_eta4p7=0;GENnjets_pt30_eta2p5=0;
        GENpt_leadingjet_pt30_eta4p7=0; GENpt_leadingjet_pt30_eta2p5=0;nVECZ=0;GENmass4l=-99;GENmass2j=-99;GENmass4e=-99;GENmass4mu=-99; GENmass2e2mu=-99;
        GENmassZ1=0;GENmassZ2=0;GENpT4l=0;GENeta4l=-99;GENrapidity4l=-99;GENphi4l=-99;
        GENpTj1=-99;GENetaj1=-99;GENphij1=-99;GENmj1=-99;
        GENpTj2=-99;GENetaj2=-99;GENphij2=-99;GENmj2=-99;
        GENZ_DaughtersId[0]=0;GENZ_DaughtersId[1]=0;
        for (int i=0; i<4; i++) {GENlep_Hindex[i]=-1;}
        for (int i=0; i<2; i++) {GENjet_Hindex[i]=-1;}
        GenPart_pt.clear(); GenPart_eta.clear(); GenPart_phi.clear(); GenPart_mass.clear(); GenPart_pdgId.clear();GenPart_status.clear();GenPart_statusFlags.clear(); GenPart_genPartIdxMother.clear();
        GenJet_pt.clear(); GenJet_eta.clear(); GenJet_phi.clear(); GenJet_mass.clear();GenJet_hadronFlavour.clear();
        GENZ_phi.clear(); GENZ_pt.clear(); GENZ_eta.clear(); GENZ_mass.clear();GENZ_MomId.clear();
        GENH_phi.clear(); GENH_pt.clear(); GENH_eta.clear(); GENH_mass.clear();
        GENjet_pt.clear();GENjet_eta.clear();GENjet_phi.clear();GENjet_mass.clear();GENjet_hadronFlavour.clear();b_jets.clear();
        GENlep_eta.clear();GENlep_pt.clear();GENlep_phi.clear();GENlep_mass.clear();GENlep_id.clear();GENlep_status.clear();GENlep_MomMomId.clear();GENlep_MomId.clear();GENlep_RelIso.clear();
        flag4e=0; flag4mu=0; flag2e2mu=0;flagpassZ1=0;flagpassFid=0;
      }
      int total_GEN_lepton_count_allEvents = 0;
      int total_GEN_tau_count_allEvents = 0;
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
      std::vector<int> GenPart_statusFlags;
      std::vector<int> GenPart_genPartIdxMother;

      std::vector<float> GenJet_pt;
      std::vector<float> GenJet_eta;
      std::vector<float> GenJet_phi;
      std::vector<float> GenJet_mass;
      std::vector<int> GenJet_hadronFlavour;
      std::vector<int> b_jets;

      unsigned nGenPart, nGenJet;
};
GenAnalysis::GenAnalysis(){
  // FIXME: Add the values to the yaml file
    genIsoConeSizeEl=0.3; genIsoConeSizeMu=0.3;
    genIsoCutEl=0.35; genIsoCutMu=0.35;
    leadingPtCut=20;subleadingPtCut=10;
    nGEN4e=0; nGEN4mu=0; nGEN2e2mu=0;nGEN4epassZ1=0;nGEN4epassFid=0;nGEN4mupassZ1=0;nGEN4mupassFid=0;nGEN2e2mupTEtaisocuts=0;nGEN2e2mupassZ1=0;nGEN2e2mupassFid=0;
} 
#endif
