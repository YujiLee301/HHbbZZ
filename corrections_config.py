#!/usr/bin/env python3
"""
Corrections configuration file for H4L analysis.
Contains all correction module setups organized by year and data tag.
"""

import os
from PhysicsTools.NATModules.modules.electronSF import ElectronSF as electronSF_natlib
from PhysicsTools.NATModules.modules.eleScaleRes import eleScaleRes as eleScaleRes_natlib
from PhysicsTools.NATModules.modules.muonSF import MuonSF as muonSF_natlib
from PhysicsTools.NATModules.modules.muonScaleRes import muonScaleRes as muonScaleRes_natlib
from PhysicsTools.NATModules.modules.puWeightProducer import puWeightProducer as puWeightProducer_natlib
from PhysicsTools.NATModules.modules.jetVetoMap import jetVMAP as jetVMAP_natlib
from PhysicsTools.NATModules.modules.jetCorr import jetJERC as jetJERC_natlib
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import (
    muonScaleRes2016, muonScaleRes2017, muonScaleRes2018
)


def get_electron_sf_2022(data_tag, isMC):
    """Get electron scale factor corrections for 2022."""
    if "pre_EE" in data_tag:
        eleSF = electronSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electron.json.gz")
        era = "2022Re-recoBCD"
    else:
        eleSF = electronSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electron.json.gz")
        era = "2022Re-recoE+PromptFG"
    
    set_name = "Electron-ID-SF"
    
    # Reco SF
    def reco_wp(pt):
        if pt < 20:
            return "RecoBelow20"
        elif pt <= 75:
            return "Reco20to75"
        else:
            return "RecoAbove75"
    
    eleSF.addCorrection(set_name, era, reco_wp, "sf", "Reco_sf")
    eleSF.addCorrection(set_name, era, reco_wp, "sfup", "Reco_sfUp")
    eleSF.addCorrection(set_name, era, reco_wp, "sfdown", "Reco_sfDown")
    
    # ID SF
    for wp in ["Loose", "Medium"]:
        eleSF.addCorrection(set_name, era, wp, "sf", f"ID_{wp}_sf")
        eleSF.addCorrection(set_name, era, wp, "sfup", f"ID_{wp}_sfUp")
        eleSF.addCorrection(set_name, era, wp, "sfdown", f"ID_{wp}_sfDown")
    
    return eleSF


def get_muon_sf_2022(data_tag, isMC):
    """Get muon scale factor corrections for 2022."""
    if "post_EE" not in data_tag:
        muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22/muon_Z.json.gz")
    else:
        muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22EE/muon_Z.json.gz")
    
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
    
    return muSF


def get_muon_scale_res_2022(data_tag, isMC, overwritePt):
    """Get muon scale/resolution corrections for 2022."""
    if "pre_EE" in data_tag:
        muonScale_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22/muon_scalesmearing.json.gz"
    else:
        muonScale_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22EE/muon_scalesmearing.json.gz"
    
    return muonScaleRes_natlib(muonScale_json, is_mc=isMC, overwritePt=overwritePt, minPt=3.)


def get_electron_scale_res_2022(data_tag, isMC, overwritePt):
    """Get electron scale/resolution corrections for 2022."""
    if "pre_EE" in data_tag:
        eleScale_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electronSS_EtDependent.json.gz"
        scaleKey = "EGMScale_Compound_Ele_2022preEE"
        smearKey = "EGMSmearAndSyst_ElePTsplit_2022preEE"
    else:
        eleScale_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electronSS_EtDependent.json.gz"
        scaleKey = "EGMScale_Compound_Ele_2022postEE"
        smearKey = "EGMSmearAndSyst_ElePTsplit_2022postEE"
    
    return eleScaleRes_natlib(eleScale_json, scaleKey, smearKey, overwritePt)


def get_electron_sf_2023(data_tag, isMC):
    """Get electron scale factor corrections for 2023."""
    if "pre_BPix" in data_tag:
        eleSF = electronSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23/electron.json.gz")
        era = "2023PromptC"
    else:
        eleSF = electronSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23BPix/electron.json.gz")
        era = "2023PromptD"
    
    set_name = "Electron-ID-SF"
    
    # Reco SF
    def reco_wp(pt):
        if pt < 20:
            return "RecoBelow20"
        elif pt <= 75:
            return "Reco20to75"
        else:
            return "RecoAbove75"
    
    eleSF.addCorrection(set_name, era, reco_wp, "sf", "Reco_sf")
    eleSF.addCorrection(set_name, era, reco_wp, "sfup", "Reco_sfUp")
    eleSF.addCorrection(set_name, era, reco_wp, "sfdown", "Reco_sfDown")
    
    # ID SF
    for wp in ["Loose", "Medium"]:
        eleSF.addCorrection(set_name, era, wp, "sf", f"ID_{wp}_sf")
        eleSF.addCorrection(set_name, era, wp, "sfup", f"ID_{wp}_sfUp")
        eleSF.addCorrection(set_name, era, wp, "sfdown", f"ID_{wp}_sfDown")
    
    return eleSF


def get_muon_sf_2023(data_tag, isMC):
    """Get muon scale factor corrections for 2023."""
    if "post_BPix" not in data_tag:
        muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23/muon_Z.json.gz")
    else:
        muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23BPix/muon_Z.json.gz")
    
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
    
    return muSF


def get_muon_scale_res_2023(data_tag, isMC, overwritePt):
    """Get muon scale/resolution corrections for 2023."""
    if "pre_BPix" in data_tag:
        muon_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23/muon_scalesmearing.json.gz"
    else:
        muon_json = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2023_Summer23BPix/muon_scalesmearing.json.gz"
    
    return muonScaleRes_natlib(muon_json, is_mc=isMC, overwritePt=overwritePt, minPt=3.)


def get_electron_scale_res_2023(data_tag, isMC, overwritePt):
    """Get electron scale/resolution corrections for 2023."""
    if "post_BPix" not in data_tag:
        json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23/electronSS_EtDependent.json.gz"
        scaleKey = "EGMScale_Compound_Ele_2023preBPIX"
        smearKey = "EGMSmearAndSyst_ElePTsplit_2023preBPIX"
    else:
        json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2023_Summer23BPix/electronSS_EtDependent.json.gz"
        scaleKey = "EGMScale_Compound_Ele_2023postBPIX"
        smearKey = "EGMSmearAndSyst_ElePTsplit_2023postBPIX"
    
    return eleScaleRes_natlib(json_path, scaleKey, smearKey, overwritePt)


def get_electron_sf_2018(isMC):
    """Get electron scale factor corrections for 2018."""
    json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2018_UL/electron.json.gz"
    eleSF = electronSF_natlib(json_path)
    set_name = "UL-Electron-ID-SF"
    eleSF.addCorrection(set_name, "2018", "RecoAbove20", "sf", "RecoAbove20_sf")
    eleSF.addCorrection(set_name, "2018", "RecoBelow20", "sf", "RecoBelow20_sf")
    eleSF.addCorrection(set_name, "2018", "Medium", "sf", "ID_sf")
    return eleSF


def get_muon_sf_2018(isMC):
    """Get muon scale factor corrections for 2018."""
    muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2018_UL/muon_Z.json.gz")
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
    return muSF


def get_muon_scale_res_2018(isMC):
    """Get muon scale/resolution corrections for 2018."""
    return muonScaleRes2018()


def get_electron_scale_res_2018(isMC, overwritePt):
    """Get electron scale/resolution corrections for 2018."""
    json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2018_UL/electronSS.json.gz"
    scaleKey = "EGMScale_Compound_Ele_2018"
    smearKey = "EGMSmearAndSyst_Ele_2018"
    return eleScaleRes_natlib(json_path, scaleKey, smearKey, isMC, overwritePt)


def get_electron_sf_2017(isMC):
    """Get electron scale factor corrections for 2017."""
    json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2017_UL/electron.json.gz"
    eleSF = electronSF_natlib(json_path)
    set_name = "UL-Electron-ID-SF"
    eleSF.addCorrection(set_name, "2017", "RecoAbove20", "sf", "RecoAbove20_sf")
    eleSF.addCorrection(set_name, "2017", "RecoBelow20", "sf", "RecoBelow20_sf")
    eleSF.addCorrection(set_name, "2017", "Medium", "sf", "ID_sf")
    return eleSF


def get_muon_sf_2017(isMC):
    """Get muon scale factor corrections for 2017."""
    muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2017_UL/muon_Z.json.gz")
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
    return muSF


def get_muon_scale_res_2017(isMC):
    """Get muon scale/resolution corrections for 2017."""
    return muonScaleRes2017()


def get_electron_scale_res_2017(isMC, overwritePt):
    """Get electron scale/resolution corrections for 2017."""
    json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2017_UL/electronSS.json.gz"
    scaleKey = "EGMScale_Compound_Ele_2017"
    smearKey = "EGMSmearAndSyst_Ele_2017"
    return eleScaleRes_natlib(json_path, scaleKey, smearKey, isMC, overwritePt)


def get_electron_sf_2016(first_file, isMC):
    """Get electron scale factor corrections for 2016."""
    if "APV" in first_file or "preVFP" in first_file:
        json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2016preVFP_UL/electron.json.gz"
    else:
        json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2016postVFP_UL/electron.json.gz"
    eleSF = electronSF_natlib(json_path)
    set_name = "UL-Electron-ID-SF"
    tag = "2016preVFP" if "preVFP" in first_file or "APV" in first_file else "2016postVFP"
    eleSF.addCorrection(set_name, tag, "RecoAbove20", "sf", "RecoAbove20_sf")
    eleSF.addCorrection(set_name, tag, "RecoBelow20", "sf", "RecoBelow20_sf")
    eleSF.addCorrection(set_name, tag, "Medium", "sf", "ID_sf")
    return eleSF


def get_muon_sf_2016(first_file, isMC):
    """Get muon scale factor corrections for 2016."""
    if "APV" in first_file or "preVFP" in first_file:
        muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2016preVFP_UL/muon_Z.json.gz")
    else:
        muSF = muonSF_natlib("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2016postVFP_UL/muon_Z.json.gz")
    
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "nominal", "MuonSF")
    muSF.addCorrection("NUM_MediumID_DEN_TrackerMuons", "syst", "MuonSFsyst")
    return muSF


def get_muon_scale_res_2016(isMC):
    """Get muon scale/resolution corrections for 2016."""
    return muonScaleRes2016()


def get_electron_scale_res_2016(first_file, isMC, overwritePt):
    """Get electron scale/resolution corrections for 2016."""
    if "APV" in first_file or "preVFP" in first_file:
        json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2016preVFP_UL/electronSS.json.gz"
        scaleKey = "EGMScale_Compound_Ele_2016preVFP"
        smearKey = "EGMSmearAndSyst_Ele_2016preVFP"
    else:
        json_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2016postVFP_UL/electronSS.json.gz"
        scaleKey = "EGMScale_Compound_Ele_2016postVFP"
        smearKey = "EGMSmearAndSyst_Ele_2016postVFP"
    
    return eleScaleRes_natlib(json_path, scaleKey, smearKey, isMC, overwritePt)


def get_pu_weight_2022(data_tag):
    """Get PU weight producer for 2022."""
    if "post_EE" not in data_tag:
        json = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/puWeights/puWeights_2022_Summer22.json.gz" % os.environ['CMSSW_BASE']
        key = "Collisions2022_355100_357900_eraBCD_GoldenJson"
    else:
        json = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/puWeights/puWeights_2022_Summer22EE.json.gz" % os.environ['CMSSW_BASE']
        key = "Collisions2022_359022_362760_eraEFG_GoldenJson"
    
    return puWeightProducer_natlib(json, key)


def get_pu_weight_2023(data_tag):
    """Get PU weight producer for 2023."""
    if "pre_BPix" in data_tag:
        json = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/puWeights/puWeights_2023_Summer23.json.gz" % os.environ['CMSSW_BASE']
        key = "Collisions2023_366403_369802_eraBC_GoldenJson"
    else:
        json = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/nanoAOD_skim/data/puWeights/puWeights_2023_Summer23BPix.json.gz" % os.environ['CMSSW_BASE']
        key = "Collisions2023_369803_370790_eraD_GoldenJson"
    
    return puWeightProducer_natlib(json, key)


def get_jet_veto_map_2022(data_tag):
    """Get jet veto map module for 2022."""
    if "pre_EE" in data_tag:
        json_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/2025-09-23/jetvetomaps.json.gz"
        corrName = "Summer22_23Sep2023_RunCD_V1"
    else:
        json_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jetvetomaps.json.gz"
        corrName = "Summer22EE_23Sep2023_RunEFG_V1"
    
    return jetVMAP_natlib(json_path, corrName, "jetvetomap")


def get_jet_veto_map_2023(data_tag):
    """Get jet veto map module for 2023."""
    if "pre_BPix" in data_tag:
        json_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jetvetomaps.json.gz"
        corrName = "Summer23Prompt23_RunC_V1"
    else:
        json_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/jetvetomaps.json.gz"
        corrName = "Summer23BPixPrompt23_RunD_V1"
    
    return jetVMAP_natlib(json_path, corrName, "jetvetomap")


def get_jet_veto_map_2024(data_tag):
    """Get jet veto map module for 2024."""
    json_path = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/jetvetomaps.json.gz"
    corrName = "Summer24Prompt24_RunBCDEFGHI_V1"
    
    return jetVMAP_natlib(json_path, corrName, "jetvetomap")


def get_jet_correction_2022(data_tag, isMC):
    """Get jet energy correction (JERC) module for 2022.
    Reference: https://github.com/cms-cat/nanoAOD-tools-modules/blob/master/test/example_jetCorr.py
    """
    # JES systematics 11-split scheme
    jes_systematics_11split = [
        "Regrouped_Absolute",
        "Regrouped_Absolute_{year}",
        "Regrouped_BBEC1",
        "Regrouped_BBEC1_{year}",
        "Regrouped_EC2",
        "Regrouped_EC2_{year}",
        "Regrouped_FlavorQCD",
        "Regrouped_HF",
        "Regrouped_HF_{year}",
        "Regrouped_RelativeBal",
        "Regrouped_RelativeSample_{year}",
    ]
    
    if isMC:
        if "pre_EE" in data_tag:
            folderKey = "Run3-22CDSep23-Summer22-NanoAODv12/latest"
            L1Key = "Summer22_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer22_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22_22Sep2023_V3_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer22_22Sep2023_V3_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer22_22Sep2023_V3_MC_{label.format(year=2022)}_AK4PFPuppi" for label in jes_systematics_11split
            ]
            smearKey = "JERSmear"
            JERKey = "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
        else:
            folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/latest"
            L1Key = "Summer22EE_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_V3_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_V3_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer22EE_22Sep2023_V3_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer22EE_22Sep2023_V3_MC_{label.format(year='2022EE')}_AK4PFPuppi" for label in jes_systematics_11split
            ]
            smearKey = "JERSmear"
            JERKey = "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"
    else:
        # Data - JER are not applied to data
        if "pre_EE" in data_tag:
            folderKey = "Run3-22CDSep23-Summer22-NanoAODv12/latest"
            L1Key = "Summer22_22Sep2023_RunCD_V3_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22_22Sep2023_RunCD_V3_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None
            smearKey = None
            JERKey = None
            JERsfKey = None
        else:
            # For post_EE, would need to handle different runs (E, F, G) - using E as default
            folderKey = "Run3-22EFGSep23-Summer22EE-NanoAODv12/latest"
            L1Key = "Summer22EE_22Sep2023_RunE_V3_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer22EE_22Sep2023_RunE_V3_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer22EE_22Sep2023_RunE_V3_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer22EE_22Sep2023_RunE_V3_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None
            smearKey = None
            JERKey = None
            JERsfKey = None
    
    json_JERC = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/%s/jet_jerc.json.gz" % (folderKey)
    json_JERsmear = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz"
    
    # Configuration flags
    overwritePt = True
    usePhiDependentJEC = False  # False for 2022
    useRunDependentJEC = False  # False for 2022
    useJesSplittingScheme11 = False
    scaleKey = scaleKeyRegrouped11 if useJesSplittingScheme11 else scaleTotalKey
    
    # Call jetJERC following reference pattern (positional arguments)
    return jetJERC_natlib(
        json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key,
        scaleKey, smearKey, JERKey, JERsfKey,
        overwritePt, usePhiDependentJEC, useRunDependentJEC
    )


def get_jet_correction_2023(data_tag, isMC):
    """Get jet energy correction (JERC) module for 2023.
    Reference: https://github.com/cms-cat/nanoAOD-tools-modules/blob/master/test/example_jetCorr.py
    """
    # JES systematics 11-split scheme
    jes_systematics_11split = [
        "Regrouped_Absolute",
        "Regrouped_Absolute_{year}",
        "Regrouped_BBEC1",
        "Regrouped_BBEC1_{year}",
        "Regrouped_EC2",
        "Regrouped_EC2_{year}",
        "Regrouped_FlavorQCD",
        "Regrouped_HF",
        "Regrouped_HF_{year}",
        "Regrouped_RelativeBal",
        "Regrouped_RelativeSample_{year}",
    ]
    
    if isMC:
        if "pre_BPix" in data_tag:
            folderKey = "Run3-23CSep23-Summer23-NanoAODv12/latest"
            L1Key = "Summer23Prompt23_V2_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer23Prompt23_V2_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer23Prompt23_V2_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer23Prompt23_V2_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer23Prompt23_V2_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer23Prompt23_V2_MC_{label.format(year='2023')}_AK4PFPuppi" for label in jes_systematics_11split
            ]
            smearKey = "JERSmear"
            JERKey = "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi"
        else:
            folderKey = "Run3-23DSep23-Summer23BPix-NanoAODv12/latest"
            L1Key = "Summer23BPixPrompt23_V3_MC_L1FastJet_AK4PFPuppi"
            L2Key = "Summer23BPixPrompt23_V3_MC_L2Relative_AK4PFPuppi"
            L3Key = "Summer23BPixPrompt23_V3_MC_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer23BPixPrompt23_V3_MC_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = "Summer23BPixPrompt23_V3_MC_Total_AK4PFPuppi"
            scaleKeyRegrouped11 = [
                f"Summer23BPixPrompt23_V3_MC_{label.format(year='2023BPix')}_AK4PFPuppi" for label in jes_systematics_11split
            ]
            smearKey = "JERSmear"
            JERKey = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"
            JERsfKey = "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"
    else:
        # Data - JER are not applied to data
        if "pre_BPix" in data_tag:
            folderKey = "Run3-23CSep23-Summer23-NanoAODv12/latest"
            L1Key = "Summer23Prompt23_V2_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer23Prompt23_V2_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None
            smearKey = None
            JERKey = None
            JERsfKey = None
        else:
            folderKey = "Run3-23DSep23-Summer23BPix-NanoAODv12/latest"
            L1Key = "Summer23BPixPrompt23_V3_DATA_L1FastJet_AK4PFPuppi"
            L2Key = "Summer23BPixPrompt23_V3_DATA_L2Relative_AK4PFPuppi"
            L3Key = "Summer23BPixPrompt23_V3_DATA_L3Absolute_AK4PFPuppi"
            L2L3Key = "Summer23BPixPrompt23_V3_DATA_L2L3Residual_AK4PFPuppi"
            scaleTotalKey = None
            scaleKeyRegrouped11 = None
            smearKey = None
            JERKey = None
            JERsfKey = None
    
    json_JERC = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/%s/jet_jerc.json.gz" % (folderKey)
    json_JERsmear = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz"
    
    # Configuration flags
    overwritePt = True
    usePhiDependentJEC = not ("pre_BPix" in data_tag)  # False for 2023 pre_BPix, True for 2023 post_BPix
    useRunDependentJEC = (not isMC)  # Apply run-dependent JEC only for 2023 data (not MC)
    useJesSplittingScheme11 = False
    scaleKey = scaleKeyRegrouped11 if useJesSplittingScheme11 else scaleTotalKey
    
    # Call jetJERC following reference pattern (positional arguments)
    return jetJERC_natlib(
        json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key,
        scaleKey, smearKey, JERKey, JERsfKey,
        overwritePt, usePhiDependentJEC, useRunDependentJEC
    )


def get_jet_correction_2024(data_tag, isMC):
    """Get jet energy correction (JERC) module for 2024.
    Reference: https://github.com/CJLST/ZZAnalysis/blob/Run3/NanoAnalysis/python/modules/jetJERC.py
    """
    # JES systematics 11-split scheme
    jes_systematics_11split = [
        "Regrouped_Absolute",
        "Regrouped_Absolute_{year}",
        "Regrouped_BBEC1",
        "Regrouped_BBEC1_{year}",
        "Regrouped_EC2",
        "Regrouped_EC2_{year}",
        "Regrouped_FlavorQCD",
        "Regrouped_HF",
        "Regrouped_HF_{year}",
        "Regrouped_RelativeBal",
        "Regrouped_RelativeSample_{year}",
    ]
    
    if isMC:
        folderKey = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest"
        L1Key = "Summer24Prompt24_V1_MC_L1FastJet_AK4PFPuppi"
        L2Key = "Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi"
        L3Key = "Summer24Prompt24_V1_MC_L3Absolute_AK4PFPuppi"
        L2L3Key = "Summer24Prompt24_V1_MC_L2L3Residual_AK4PFPuppi"
        scaleTotalKey = "Summer24Prompt24_V1_MC_Total_AK4PFPuppi"
        scaleKeyRegrouped11 = [
            f"Summer24Prompt24_V1_MC_{label.format(year='2024')}_AK4PFPuppi" for label in jes_systematics_11split
        ]
        smearKey = "JERSmear"
        # It appears the 23BPix keys are used for the following:
        JERKey = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"
        JERsfKey = "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"
    else:
        # Data - JER are not applied to data
        folderKey = "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest"
        L1Key = "Summer24Prompt24_V1_DATA_L1FastJet_AK4PFPuppi"
        L2Key = "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi"
        L3Key = "Summer24Prompt24_V1_DATA_L3Absolute_AK4PFPuppi"
        L2L3Key = "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi"
        scaleTotalKey = None
        scaleKeyRegrouped11 = None
        smearKey = None
        JERKey = None
        JERsfKey = None
    
    json_JERC = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/%s/jet_jerc.json.gz" % (folderKey)
    json_JERsmear = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz"
    
    # Configuration flags
    overwritePt = True
    usePhiDependentJEC = True  # True for 2024 (era >= 2023 and not pre_BPix)
    useRunDependentJEC = (not isMC)  # Apply run-dependent JEC only for 2024 data (not MC)
    useJesSplittingScheme11 = False
    scaleKey = scaleKeyRegrouped11 if useJesSplittingScheme11 else scaleTotalKey
    
    # Call jetJERC following reference pattern (positional arguments)
    return jetJERC_natlib(
        json_JERC, json_JERsmear, L1Key, L2Key, L3Key, L2L3Key,
        scaleKey, smearKey, JERKey, JERsfKey,
        overwritePt, usePhiDependentJEC, useRunDependentJEC
    )


def get_corrections_modules(year, data_tag, first_file, isMC, overwritePt):
    """
    Get all correction modules for a given year and data tag.
    
    Args:
        year: Year (2016, 2017, 2018, 2022, 2023, 2024)
        data_tag: Data tag (e.g., "pre_EE", "post_EE", "pre_BPix", "post_BPix")
        first_file: First file name (used for 2016 VFP determination)
        isMC: Whether processing MC
        overwritePt: Whether to overwrite pT
    
    Returns:
        List of correction modules
    """
    modules = []
    
    if year == 2022:
        modules.append(get_electron_sf_2022(data_tag, isMC))
        modules.append(get_muon_sf_2022(data_tag, isMC))
        modules.append(get_muon_scale_res_2022(data_tag, isMC, overwritePt))
        modules.append(get_electron_scale_res_2022(data_tag, isMC, overwritePt))
        modules.append(get_jet_veto_map_2022(data_tag))
        modules.append(get_jet_correction_2022(data_tag, isMC))
    
    elif year == 2023:
        modules.append(get_electron_sf_2023(data_tag, isMC))
        modules.append(get_muon_sf_2023(data_tag, isMC))
        modules.append(get_muon_scale_res_2023(data_tag, isMC, overwritePt))
        modules.append(get_electron_scale_res_2023(data_tag, isMC, overwritePt))
        modules.append(get_jet_veto_map_2023(data_tag))
        modules.append(get_jet_correction_2023(data_tag, isMC))
    
    elif year == 2024:
        # Note: 2024 corrections may need to be added for electron/muon SF and scale/res
        # For now, only jet corrections are implemented
        modules.append(get_jet_veto_map_2024(data_tag))
        modules.append(get_jet_correction_2024(data_tag, isMC))
    
    elif year == 2018:
        modules.append(get_electron_sf_2018(isMC))
        modules.append(get_muon_sf_2018(isMC))
        modules.append(get_muon_scale_res_2018(isMC))
        modules.append(get_electron_scale_res_2018(isMC, overwritePt))
    
    elif year == 2017:
        modules.append(get_electron_sf_2017(isMC))
        modules.append(get_muon_sf_2017(isMC))
        modules.append(get_muon_scale_res_2017(isMC))
        modules.append(get_electron_scale_res_2017(isMC, overwritePt))
    
    elif year == 2016:
        modules.append(get_electron_sf_2016(first_file, isMC))
        modules.append(get_muon_sf_2016(first_file, isMC))
        modules.append(get_muon_scale_res_2016(isMC))
        modules.append(get_electron_scale_res_2016(first_file, isMC, overwritePt))
    
    return modules


def get_pu_weight_module(year, data_tag):
    """
    Get PU weight module for MC processing.
    
    Args:
        year: Year (2022, 2023)
        data_tag: Data tag
    
    Returns:
        PU weight module or None
    """
    if year == 2022:
        return get_pu_weight_2022(data_tag)
    elif year == 2023:
        return get_pu_weight_2023(data_tag)
    
    return None

